module surface_chion
    ! chion as the yelmox surface model (surface_method = "chion").
    !
    ! chion (fesmc/chion) steps independent snowpack columns one step at a time
    ! and leaves the time loop to the host. This wrapper is that host: it packs
    ! the (nx,ny) fields of grid_surface into chion's column list, drives one
    ! annual cycle of daily steps from monthly climate, and aggregates the result
    ! to the annual fields yelmox lands on Yelmo (dom%schn%ann), in smbpal's
    ! units so couple_surface_to_yelmo treats both models alike.
    !
    ! Inputs chion does not own, supplied here (as smbpal does internally):
    !   - daily top-of-atmosphere insolation (libs/insol), as shortwave_down,
    !     from a per-day latitude table (insol_dlat) interpolated to the columns;
    !   - the annual positive degree days (ITM's critical snow depth), from the
    !     daily temperature with chion's Calov-Greve integral, fixed for the year;
    !   - the snowfall/rainfall split (smbpal's calc_snowfrac).
    ! Monthly -> daily is chion_forcing_monthly (mean-preserving).
    !
    ! Active columns: land or ice, not open ocean (H_ice > 0 or z_bed > z_sl),
    ! refreshed every annual cycle. A column switched on starts from chion's
    ! cold state. Inactive columns: smb = 0, tsrf = annual-mean air temperature.
    !
    ! smb is the surface mass balance, sf + rf - runoff - subl (smbpal's smb),
    ! not chion_get_smb's ice-facing flux. tsrf is computed here once per year
    ! from the annual-mean air temperature and annual net melt (refrz - melt),
    ! as smbpal does: chion's per-step ITM tsrf applies firn_fac to a daily
    ! melt_net, whereas firn_fac is calibrated against the annual total.
    !
    ! The per-column host work (forcing, accumulation, PDDs, insolation lookup)
    ! is OpenMP-parallel over the active columns, like chion's own step.
    !
    ! Only model = "itm" is supported: the host supplies ITM's forcing only.

    use nml,        only : nml_read
    use insolation, only : calc_insol_day
    use chion,      only : wp, dp, wp_acc, MV, chion_class, chion_init, chion_init_state, &
                           chion_update, chion_end, chion_set_active_mask, &
                           chion_get_surface_flux_totals, chion_get_surface, &
                           chion_restart_write, chion_restart_read, &
                           monthly_to_daily_class, monthly_to_daily_init, &
                           monthly_to_daily_controls, interp_monthly_to_day, &
                           pdd_expected_positive_temperature

    implicit none

    private

    integer,  parameter :: nmon       = 12
    integer,  parameter :: nday_mon   = 30       ! 360-day year, as smbpal and the climate
    real(dp), parameter :: insol_dlat = 0.1_dp   ! [deg] latitude spacing of the insolation table

    type surface_chion_param_class
        logical            :: const_insol   ! insolation at const_kabp instead of the model time
        real(wp)           :: const_kabp    ! [kyr BP] insolation time when const_insol
        character(len=512) :: insol_fldr    ! folder of the orbital-parameter tables
        real(wp)           :: sf_a, sf_b    ! snowfall fraction, -0.5*tanh(sf_a*(T-sf_b))+0.5
        real(wp)           :: sigma_pdd     ! [K] temperature standard deviation for the PDDs
        real(wp)           :: firn_fac      ! [K (mm w.e.)-1] firn warming per annual net refreezing
        real(wp)           :: time_equil    ! [yr] cold-start snowpack spin-up
        integer            :: dt_days       ! [d] chion step; divides the 360-day year
    end type surface_chion_param_class

    type surface_chion_ann_class
        ! Annual means on grid_surface. Mass fluxes in [mm w.e./yr].
        real(wp), allocatable :: mask(:,:)      ! [1] chion active (land or ice)
        real(wp), allocatable :: t2m(:,:)       ! [K]
        real(wp), allocatable :: pr(:,:)        ! [mm w.e./yr]
        real(wp), allocatable :: sf(:,:)        ! [mm w.e./yr]
        real(wp), allocatable :: S(:,:)         ! [W m-2] TOA insolation
        real(wp), allocatable :: PDDs(:,:)      ! [K d]
        real(wp), allocatable :: melt(:,:)      ! [mm w.e./yr]
        real(wp), allocatable :: runoff(:,:)    ! [mm w.e./yr]
        real(wp), allocatable :: refrz(:,:)     ! [mm w.e./yr]
        real(wp), allocatable :: smb(:,:)       ! [mm w.e./yr] sf + rf - runoff - subl
        real(wp), allocatable :: tsrf(:,:)      ! [K]
        real(wp), allocatable :: alb_s(:,:)     ! [1] end-of-year surface albedo
        real(wp), allocatable :: H_snow(:,:)    ! [mm w.e.] end-of-year snowpack (itm)
    end type surface_chion_ann_class

    type surface_chion_insol_class
        ! Daily TOA insolation on a regular latitude table, rebuilt only when
        ! the insolation time changes (once per run with const_insol), and the
        ! fixed linear-interpolation stencil of every column into it.
        real(dp)              :: time = huge(1.0_dp)   ! [yr] insolation time of the table
        real(dp), allocatable :: lat(:)                ! (nlat) [deg N]
        real(wp), allocatable :: S(:,:)                ! (nlat,nday) [W m-2]
        integer,  allocatable :: j0(:)                 ! (ncol) table row below the column
        real(wp), allocatable :: w1(:)                 ! (ncol) weight of row j0+1
    end type surface_chion_insol_class

    type surface_chion_class
        type(surface_chion_param_class) :: par
        type(chion_class)               :: chn
        type(monthly_to_daily_class)    :: md
        type(surface_chion_insol_class) :: ins
        type(surface_chion_ann_class)   :: ann
        integer                         :: nx, ny, ncol
    end type surface_chion_class

    public :: surface_chion_class
    public :: surface_chion_init
    public :: surface_chion_spinup
    public :: surface_chion_update
    public :: surface_chion_restart_write
    public :: surface_chion_restart_read
    public :: surface_chion_end

contains

    subroutine surface_chion_init(sc, filename, lats, group, chion_group)
        ! Parameters, chion (cold state) and the annual fields on an (nx,ny) grid.
        type(surface_chion_class), intent(inout) :: sc
        character(len=*),          intent(in)    :: filename
        real(wp),                  intent(in)    :: lats(:,:)      ! [deg N]
        character(len=*),          intent(in)    :: group          ! &surface_chion
        character(len=*),          intent(in)    :: chion_group    ! &chion

        call surface_chion_par_load(sc%par, filename, group)

        sc%nx   = size(lats,1)
        sc%ny   = size(lats,2)
        sc%ncol = sc%nx*sc%ny

        call chion_init(sc%chn, filename, sc%ncol, group=chion_group)

        if (trim(sc%chn%par%model) /= "itm") then
            write(*,*) "surface_chion_init:: error: only chion model = itm is supported &
                       &(the host supplies ITM forcing only); got "//trim(sc%chn%par%model)
            stop 1
        end if

        call monthly_to_daily_init(sc%md, nmon, nday_mon)

        if (sc%par%dt_days < 1 .or. mod(sc%md%nday_year, sc%par%dt_days) /= 0) then
            write(*,*) "surface_chion_init:: error: dt_days must divide the year; dt_days, nday_year = ", &
                       sc%par%dt_days, sc%md%nday_year
            stop 1
        end if

        ! Column icol = i + (j-1)*nx: reshape order.
        sc%chn%forc%latitude_deg = reshape(lats, [sc%ncol])
        call insol_init(sc%ins, sc%chn%forc%latitude_deg, sc%md%nday_year)

        call chion_init_state(sc%chn)

        call ann_alloc(sc%ann, sc%nx, sc%ny)
    end subroutine surface_chion_init

    subroutine surface_chion_spinup(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, time_bp)
        ! Cold-start snowpack: time_equil annual cycles under the initial climate
        ! and geometry (smbpal's 100-yr ITM equilibration).
        type(surface_chion_class), intent(inout) :: sc
        real(wp), intent(in) :: t2m(:,:,:), pr(:,:,:)       ! [K], [mm w.e./d] monthly
        real(wp), intent(in) :: z_srf(:,:), z_bed(:,:), H_ice(:,:), z_sl(:,:)
        real(wp), intent(in) :: time_bp

        integer :: n

        do n = 1, nint(sc%par%time_equil)
            call surface_chion_update(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, time_bp)
        end do
    end subroutine surface_chion_spinup

    subroutine surface_chion_update(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, time_bp)
        ! One annual cycle of daily chion steps from monthly climate, aggregated
        ! to sc%ann.
        type(surface_chion_class), intent(inout) :: sc
        real(wp), intent(in) :: t2m(:,:,:)          ! [K] monthly near-surface air temperature
        real(wp), intent(in) :: pr(:,:,:)           ! [mm w.e./d] monthly precipitation
        real(wp), intent(in) :: z_srf(:,:)          ! [m]
        real(wp), intent(in) :: z_bed(:,:)          ! [m]
        real(wp), intent(in) :: H_ice(:,:)          ! [m]
        real(wp), intent(in) :: z_sl(:,:)           ! [m]
        real(wp), intent(in) :: time_bp             ! [yr] model time, for the insolation

        real(wp), allocatable :: t_ctl(:,:), p_ctl(:,:), mon(:,:), ctl(:,:), col(:)
        real(wp), allocatable :: t_sum(:), pr_sum(:), sf_sum(:), S_sum(:), PDDs(:)
        real(wp), allocatable :: tsrf_c(:), alb_c(:)
        real(wp_acc), allocatable :: melt0(:), runoff0(:), refrz0(:), subl0(:)
        real(wp_acc), allocatable :: melt1(:), runoff1(:), refrz1(:), subl1(:)
        real(wp), allocatable :: melt_a(:), runoff_a(:), refrz_a(:), subl_a(:)
        integer,  allocatable :: idx(:)
        integer  :: na, i, icol, m, day, nday, j0
        real(wp) :: dt, spd, T0, t_d, p_d, sf_d, S_d
        real(dp) :: insol_time

        nday = sc%md%nday_year
        dt   = real(sc%par%dt_days, wp)
        spd  = sc%chn%c%seconds_per_day
        T0   = sc%chn%c%T0

        ! Active columns, from the current geometry: land or ice.
        call chion_set_active_mask(sc%chn, reshape(H_ice > 0.0_wp .or. z_bed > z_sl, [sc%ncol]))
        na  = sc%chn%grd%n_active
        idx = sc%chn%grd%active_idx(1:na)

        insol_time = real(time_bp, dp)
        if (sc%par%const_insol) insol_time = real(sc%par%const_kabp, dp)*1e3_dp
        call insol_update(sc%ins, insol_time, sc%par%dt_days, trim(sc%par%insol_fldr))

        ! Mean-preserving daily control values, (nmon,na) so a column is contiguous.
        allocate(mon(na,nmon), ctl(na,nmon))
        do m = 1, nmon
            col = reshape(t2m(:,:,m), [sc%ncol])
            mon(:,m) = col(idx)
        end do
        call monthly_to_daily_controls(sc%md, mon, ctl)
        t_ctl = transpose(ctl)
        do m = 1, nmon
            col = reshape(pr(:,:,m), [sc%ncol])
            mon(:,m) = col(idx)
        end do
        call monthly_to_daily_controls(sc%md, mon, ctl)
        p_ctl = transpose(ctl)
        deallocate(mon, ctl)

        ! Fixed for the year: geometry and the annual positive degree days.
        col = reshape(z_srf, [sc%ncol])
        sc%chn%forc%surface_height(idx) = col(idx)
        col = reshape(H_ice, [sc%ncol])
        sc%chn%forc%H_ice(idx) = col(idx)

        allocate(PDDs(na))
        !$omp parallel do default(shared) private(i,day,t_d)
        do i = 1, na
            PDDs(i) = 0.0_wp
            do day = 1, nday
                call interp_monthly_to_day(sc%md, t_ctl(:,i), day, t_d)
                PDDs(i) = PDDs(i) + real(pdd_expected_positive_temperature(real(t_d - T0, dp), &
                                                                       real(sc%par%sigma_pdd, dp)), wp)
            end do
            sc%chn%forc%PDDs(idx(i)) = PDDs(i)
        end do
        !$omp end parallel do

        ! Annual flux totals are the difference of chion's running totals.
        allocate(melt0(sc%ncol), runoff0(sc%ncol), refrz0(sc%ncol), subl0(sc%ncol))
        allocate(melt1(sc%ncol), runoff1(sc%ncol), refrz1(sc%ncol), subl1(sc%ncol))
        call chion_get_surface_flux_totals(sc%chn, melt=melt0, runoff=runoff0, refrz=refrz0, subl=subl0)

        allocate(t_sum(na), pr_sum(na), sf_sum(na), S_sum(na))
        t_sum = 0.0_wp;  pr_sum = 0.0_wp;  sf_sum = 0.0_wp;  S_sum = 0.0_wp

        do day = 1, nday, sc%par%dt_days
            sc%chn%forc%day_of_year = real(day, wp)

            !$omp parallel do default(shared) private(i,icol,j0,t_d,p_d,sf_d,S_d)
            do i = 1, na
                icol = idx(i)
                call interp_monthly_to_day(sc%md, t_ctl(:,i), day, t_d)
                call interp_monthly_to_day(sc%md, p_ctl(:,i), day, p_d)
                ! The mean-preserving controls can undershoot zero for a peaked field.
                p_d  = max(p_d, 0.0_wp)
                sf_d = p_d*calc_snowfrac(t_d, sc%par%sf_a, sc%par%sf_b)
                j0   = sc%ins%j0(icol)
                S_d  = (1.0_wp - sc%ins%w1(icol))*sc%ins%S(j0,day) + sc%ins%w1(icol)*sc%ins%S(j0+1,day)

                sc%chn%forc%air_temperature(icol) = t_d
                sc%chn%forc%snowfall_rate(icol)   = sf_d/spd             ! [mm/d] -> [kg m-2 s-1]
                sc%chn%forc%rainfall_rate(icol)   = (p_d - sf_d)/spd
                sc%chn%forc%shortwave_down(icol)  = S_d

                t_sum(i)  = t_sum(i)  + t_d*dt
                pr_sum(i) = pr_sum(i) + p_d*dt
                sf_sum(i) = sf_sum(i) + sf_d*dt
                S_sum(i)  = S_sum(i)  + S_d*dt
            end do
            !$omp end parallel do

            call chion_update(sc%chn, dt)
        end do

        call chion_get_surface_flux_totals(sc%chn, melt=melt1, runoff=runoff1, refrz=refrz1, subl=subl1)

        ! Annual totals [mm w.e./yr] == [kg m-2/yr]; MV (unresolved) -> 0.
        melt_a   = annual_total(melt1(idx),   melt0(idx))
        runoff_a = annual_total(runoff1(idx), runoff0(idx))
        refrz_a  = annual_total(refrz1(idx),  refrz0(idx))
        subl_a   = annual_total(subl1(idx),   subl0(idx))

        allocate(tsrf_c(sc%ncol), alb_c(sc%ncol))
        call chion_get_surface(sc%chn, t_srf=tsrf_c, albedo=alb_c)

        sc%ann%mask = unpack_col(spread(1.0_wp, 1, na), idx, sc, 0.0_wp)
        sc%ann%t2m  = sum(t2m, dim=3)/real(nmon, wp)
        call scatter(sc%ann%t2m,    t_sum/real(nday, wp), idx, sc)
        sc%ann%pr   = sum(pr, dim=3)/real(nmon, wp)*real(nday, wp)
        call scatter(sc%ann%pr,     pr_sum, idx, sc)
        sc%ann%sf     = unpack_col(sf_sum,                       idx, sc, 0.0_wp)
        sc%ann%S      = unpack_col(S_sum/real(nday, wp),         idx, sc, 0.0_wp)
        sc%ann%PDDs   = unpack_col(PDDs,                         idx, sc, 0.0_wp)
        sc%ann%melt   = unpack_col(melt_a,                       idx, sc, 0.0_wp)
        sc%ann%runoff = unpack_col(runoff_a,                     idx, sc, 0.0_wp)
        sc%ann%refrz  = unpack_col(refrz_a,                      idx, sc, 0.0_wp)
        sc%ann%smb    = unpack_col(pr_sum - runoff_a - subl_a,   idx, sc, 0.0_wp)
        sc%ann%alb_s  = unpack_col(alb_c(idx),                   idx, sc, MV)

        ! Surface temperature from the annual means; inactive columns keep t2m.
        sc%ann%tsrf = sc%ann%t2m
        col = reshape(H_ice, [sc%ncol])
        call scatter(sc%ann%tsrf, calc_temp_surf(t_sum/real(nday, wp), col(idx), &
                                                 refrz_a - melt_a, sc%par%firn_fac, T0), idx, sc)

        sc%ann%H_snow = MV
        if (trim(sc%chn%par%model) == "itm") &
            call scatter(sc%ann%H_snow, sc%chn%itm%now%H_snow(idx), idx, sc)
    end subroutine surface_chion_update

    subroutine surface_chion_restart_write(sc, filename, time)
        type(surface_chion_class), intent(in) :: sc
        character(len=*),          intent(in) :: filename
        real(wp),                  intent(in) :: time

        call chion_restart_write(sc%chn, filename, time)
    end subroutine surface_chion_restart_write

    subroutine surface_chion_restart_read(sc, filename)
        ! Snowpack state (and active mask) from a chion restart file; the annual
        ! fields are rebuilt by the next surface_chion_update.
        type(surface_chion_class), intent(inout) :: sc
        character(len=*),          intent(in)    :: filename

        real(wp) :: time_file

        call chion_restart_read(sc%chn, filename, time_file)
    end subroutine surface_chion_restart_read

    subroutine surface_chion_end(sc)
        type(surface_chion_class), intent(inout) :: sc

        call chion_end(sc%chn)
    end subroutine surface_chion_end

    subroutine surface_chion_par_load(par, filename, group)
        type(surface_chion_param_class), intent(out) :: par
        character(len=*),                intent(in)  :: filename
        character(len=*),                intent(in)  :: group

        call nml_read(filename, group, "const_insol", par%const_insol)
        call nml_read(filename, group, "const_kabp",  par%const_kabp)
        call nml_read(filename, group, "insol_fldr",  par%insol_fldr)
        call nml_read(filename, group, "sf_a",        par%sf_a)
        call nml_read(filename, group, "sf_b",        par%sf_b)
        call nml_read(filename, group, "sigma_pdd",   par%sigma_pdd)
        call nml_read(filename, group, "firn_fac",    par%firn_fac)
        call nml_read(filename, group, "time_equil",  par%time_equil)
        call nml_read(filename, group, "dt_days",     par%dt_days)
    end subroutine surface_chion_par_load

    subroutine ann_alloc(ann, nx, ny)
        type(surface_chion_ann_class), intent(inout) :: ann
        integer,                       intent(in)    :: nx, ny

        allocate(ann%mask(nx,ny), ann%t2m(nx,ny), ann%pr(nx,ny), ann%sf(nx,ny))
        allocate(ann%S(nx,ny), ann%PDDs(nx,ny), ann%melt(nx,ny), ann%runoff(nx,ny))
        allocate(ann%refrz(nx,ny), ann%smb(nx,ny), ann%tsrf(nx,ny), ann%alb_s(nx,ny))
        allocate(ann%H_snow(nx,ny))

        ann%mask   = 0.0_wp
        ann%t2m    = 0.0_wp
        ann%pr     = 0.0_wp
        ann%sf     = 0.0_wp
        ann%S      = 0.0_wp
        ann%PDDs   = 0.0_wp
        ann%melt   = 0.0_wp
        ann%runoff = 0.0_wp
        ann%refrz  = 0.0_wp
        ann%smb    = 0.0_wp
        ann%tsrf   = 0.0_wp
        ann%alb_s  = MV
        ann%H_snow = MV
    end subroutine ann_alloc

    function unpack_col(val, idx, sc, fill) result(var)
        ! Active-column values onto the (nx,ny) grid, fill elsewhere.
        real(wp),                  intent(in) :: val(:)
        integer,                   intent(in) :: idx(:)
        type(surface_chion_class), intent(in) :: sc
        real(wp),                  intent(in) :: fill
        real(wp) :: var(sc%nx,sc%ny)

        var = fill
        call scatter(var, val, idx, sc)
    end function unpack_col

    subroutine scatter(var, val, idx, sc)
        ! Overwrite the active columns of an (nx,ny) field.
        real(wp),                  intent(inout) :: var(:,:)
        real(wp),                  intent(in)    :: val(:)
        integer,                   intent(in)    :: idx(:)
        type(surface_chion_class), intent(in)    :: sc

        real(wp), allocatable :: col(:)

        col = reshape(var, [sc%ncol])
        col(idx) = val
        var = reshape(col, [sc%nx, sc%ny])
    end subroutine scatter

    elemental function annual_total(tot1, tot0) result(y)
        ! Difference of two running totals; chion reports MV for a flux the
        ! model does not resolve, which counts as zero here.
        real(wp_acc), intent(in) :: tot1, tot0
        real(wp) :: y

        y = 0.0_wp
        if (tot1 /= real(MV, wp_acc) .and. tot0 /= real(MV, wp_acc)) y = real(tot1 - tot0, wp)
    end function annual_total

    subroutine insol_init(ins, lats, nday)
        ! Regular latitude table (insol_dlat) over [-90,90] and each column's
        ! linear-interpolation stencil into it.
        type(surface_chion_insol_class), intent(inout) :: ins
        real(wp),                        intent(in)    :: lats(:)     ! (ncol) [deg N]
        integer,                         intent(in)    :: nday

        integer :: nlat, j, icol
        real(dp) :: x

        nlat = nint(180.0_dp/insol_dlat) + 1
        allocate(ins%lat(nlat), ins%S(nlat,nday))
        do j = 1, nlat
            ins%lat(j) = -90.0_dp + real(j-1, dp)*insol_dlat
        end do
        ins%S    = 0.0_wp
        ins%time = huge(1.0_dp)

        allocate(ins%j0(size(lats)), ins%w1(size(lats)))
        do icol = 1, size(lats)
            x = (real(lats(icol), dp) + 90.0_dp)/insol_dlat
            ins%j0(icol) = min(max(int(x) + 1, 1), nlat - 1)
            ins%w1(icol) = real(x - real(ins%j0(icol) - 1, dp), wp)
        end do
    end subroutine insol_init

    subroutine insol_update(ins, time, dt_days, fldr)
        ! Daily TOA insolation on the latitude table for the stepped days, only
        ! when the insolation time has changed since the last call.
        type(surface_chion_insol_class), intent(inout) :: ins
        real(dp),                        intent(in)    :: time      ! [yr]
        integer,                         intent(in)    :: dt_days
        character(len=*),                intent(in)    :: fldr

        integer :: day

        if (time == ins%time) return

        do day = 1, size(ins%S,2), dt_days
            ins%S(:,day) = real(calc_insol_day(day, ins%lat, time, fldr=fldr), wp)
        end do
        ins%time = time
    end subroutine insol_update


    elemental function calc_snowfrac(t2m, a, b) result(f)
        ! smbpal calc_snowfrac: snow fraction of total precipitation.
        real(wp), intent(in) :: t2m, a, b
        real(wp) :: f

        f = -0.5_wp*tanh(a*(t2m - b)) + 0.5_wp
    end function calc_snowfrac

    elemental function calc_temp_surf(tann, H_ice, melt_net, fac, T0) result(ts)
        ! smbpal calc_temp_surf, on annual means: net refreezing [mm w.e./yr]
        ! warms the firn, capped at the freezing point on ice.
        real(wp), intent(in) :: tann, H_ice, melt_net, fac, T0
        real(wp) :: ts

        if (H_ice > 0.0_wp) then
            ts = min(T0, tann + fac*max(0.0_wp, melt_net))
        else
            ts = tann
        end if
    end function calc_temp_surf

end module surface_chion
