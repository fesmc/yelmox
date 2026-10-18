module yelmox_climate
    ! Climate backend of a domain, chosen at runtime ([comps] climate):
    !   "snapesm"  -- snapshots blended by indices (snapesm)
    !   "esm"      -- reference climatology plus Earth-system-model anomalies over
    !                 historical/projection periods (esm_forcing)
    !   "rembo"    -- REMBO atmosphere and surface mass balance (climate_rembo,
    !                 built with rembo=1), with the ocean from snapesm (an
    !                 ocean-only reference snapshot)
    !
    ! Presents one interface (climate_init / climate_update) over the backends,
    ! filling a `climate_out_class`; the domain and its coupling read only that.
    ! The driver's indices registry is read by snapesm (its axes) and by the
    ! climate modifiers (climate_mods: anom, ocn_fraction, bias), applied to
    ! `climate_out` after the backend in a fixed order, each checked at init
    ! against what the backend supplies; rembo reads [anom] itself
    ! (rembo_update).

    use precision,     only : wp, dp
    use nml,           only : nml_read
    use ncio
    use coordinates,   only : grid_class
    use timestepping,  only : tstep_class
    use climate_out,   only : climate_out_class, clim_state_class, summer_mean
    use snapesm,       only : snapesm_class, climate_state_class, snapesm_init, snapesm_update, &
                              snapesm_end, snapesm_write_init, snapesm_write_step, &
                              snapesm_lapse_monthly
    use esm_forcing,   only : esm_forcing_class, esm_forcing_init, esm_clim_update, &
                              esm_forcing_update, esm_variability_update, &
                              esm_lapse_monthly
    use marine_shelf,  only : marshelf_class, marshelf_interp_shelf, ocn_variable_extrapolation
    use kryos_forcing, only : indices_class, indices_value
    use climate_mods,  only : climate_caps_class, climate_mods_class, climate_mods_init, &
                              climate_mods_apply, climate_mods_write_init, climate_mods_write_step
    use climate_rembo, only : rembo_clim_init, rembo_clim_update, rembo_clim_restart_write
    use regions,       only : regions_class, regions_basin_ids

    implicit none

    private

    ! Run control of the esm backend: [esm] holds the experiment, the timeline
    ! group of the run phase its periods.
    type esm_ctl_class
        character(len=56)  :: run_type        ! the timeline group, e.g. spinup | transient
        character(len=512) :: par_file        ! esm data configuration
        character(len=56)  :: experiment, esm_name
        logical            :: use_esm, use_smb, use_var, use_hist, use_proj
        real(dp)           :: time_ref(2), time_hist(2), time_proj(2), time_esm_ref(2)
        character(len=56)  :: clim_var
        integer            :: clim_seed
    end type esm_ctl_class

    ! State of the rembo backend beyond the climate products: REMBO's summer
    ! air temperature and the anomalies applied from [anom].
    type rembo_state_class
        real(wp), allocatable :: ta_sum(:,:)     ! [K] near-surface air temperature (summer)
        real(wp) :: dT_sum = 0.0_wp              ! [K] summer air temperature anomaly
        real(wp) :: dT_ann = 0.0_wp              ! [K] annual air temperature anomaly
        real(wp) :: dT_ocn = 0.0_wp              ! [K] ocean temperature anomaly
    end type rembo_state_class

    type yelmox_climate_class
        character(len=16)       :: method = ""   ! snapesm | esm | rembo
        character(len=256)      :: grid_name     ! the climate grid (grid_clim)
        logical                 :: south = .false.   ! southern hemisphere
        type(snapesm_class)     :: snapesm
        type(esm_forcing_class) :: esm
        type(esm_ctl_class)     :: esm_ctl
        type(rembo_state_class) :: rembo
        type(climate_mods_class):: mods          ! modifiers after the backend
        real(wp), allocatable   :: basins(:,:)   ! basin ids on the climate grid (esm: [esm] basins)
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update
    public :: climate_end
    public :: climate_file_base, climate_write_init, climate_write_2D, climate_write_1D
    public :: climate_restart_write

    ! Copy of the fields a state supplies (allocated) into climate_out.
    interface copy_allocated
        module procedure copy_allocated_1D, copy_allocated_2D, copy_allocated_3D
    end interface

contains

    subroutine climate_init(cl, method, filename, domain, grid, time, reg, &
                            south, sfx, timeline_group, smb_direct, idx)
        ! south: the domain lies in the southern hemisphere (seasons, lapse rates).
        ! idx: the driver's indices registry (snapesm and the modifiers check the
        ! indices they name against it).
        ! smb_direct: the surface mass balance is taken from the climate
        ! ([comps] surface_method = "climate"); the esm backend supplies it, and
        ! the rembo backend supplies nothing else (annual fields only).
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: method
        character(len=*), intent(in) :: filename, domain
        type(grid_class), intent(in) :: grid             ! the climate grid
        real(dp),         intent(in) :: time
        type(regions_class), intent(in) :: reg           ! the regions on the climate grid
        logical,          intent(in) :: south
        character(len=*), intent(in) :: sfx              ! namelist group suffix of the domain
        character(len=*), intent(in) :: timeline_group   ! group of the run phase's timeline
        logical,          intent(in) :: smb_direct
        type(indices_class), intent(in) :: idx

        integer :: nx, ny
        type(climate_caps_class) :: caps

        nx        = grid%G%nx
        ny        = grid%G%ny

        cl%method    = trim(method)
        cl%grid_name = trim(grid%name)
        cl%south     = south

        if (smb_direct .and. trim(cl%method) /= "esm" .and. trim(cl%method) /= "rembo") then
            write(*,*) "climate_init:: error: surface_method = climate needs a climate that supplies &
                       &the surface mass balance (esm, rembo); got climate = ", trim(cl%method)
            error stop 1
        end if
        if (.not. smb_direct .and. trim(cl%method) == "rembo") then
            write(*,*) "climate_init:: error: climate = rembo supplies the surface mass balance &
                       &only; set surface_method = climate."
            error stop 1
        end if

        ! The backend, and what it supplies (caps) for the modifiers' contracts.
        ! It supplies the surface mass balance exactly when it is taken from it
        ! (smb_direct: esm with use_smb, rembo; checked above).
        caps%name    = trim(cl%method)
        caps%has_smb = smb_direct
        select case(trim(cl%method))
            case("snapesm")
                call snapesm_init(cl%snapesm, filename, domain, grid, time, &
                                  south, idx, group="snap"//trim(sfx))
                ! The surface reads the atmosphere; an ocean-only reference
                ! serves the rembo ocean only.
                if (.not. cl%snapesm%par%with_atm) then
                    write(*,*) "climate_init:: error: climate = snapesm needs an atmosphere: &
                               &the reference snapshot ([snap"//trim(sfx)//"] ref_name) binds &
                               &no tas/pr."
                    error stop 1
                end if
                caps%has_tsl_ann = cl%snapesm%par%with_atm
                caps%monthly     = cl%snapesm%par%with_atm
            case("esm")
                call esm_init(cl, filename, domain, grid, "esm"//trim(sfx), timeline_group, &
                              smb_direct, reg)
                caps%has_ocn_shelf = .true.
            case("rembo")
                ! REMBO for the atmosphere and smb, snapesm for the ocean: an
                ! ocean-only reference snapshot (to, optionally so).
                call snapesm_init(cl%snapesm, filename, domain, grid, time, &
                                  south, idx, group="snap"//trim(sfx))
                if (cl%snapesm%par%with_atm .or. .not. cl%snapesm%par%with_ocn) then
                    write(*,*) "climate_init:: error: climate = rembo takes the ocean only from &
                               &[snap"//trim(sfx)//"]: the reference snapshot must bind to &
                               &(and optionally so), not tas/pr (REMBO supplies the atmosphere)."
                    error stop 1
                end if
                call rembo_clim_init(time, nx, ny)
                allocate(cl%rembo%ta_sum(nx,ny))
                cl%rembo%ta_sum = 0.0_wp
                caps%monthly    = .false.
                caps%anom_input = .true.
            case default
                write(*,*) "climate_init:: error: climate must be snapesm, esm or rembo; got ", &
                           trim(cl%method)
                error stop 1
        end select

        ! The modifiers, checked against the backend.
        call climate_mods_init(cl%mods, filename, sfx, idx, caps, domain, grid, reg)

    end subroutine climate_init

    subroutine climate_update(cl, out, ts, z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref, &
                              domain, dx, dtt, mshlf, idx, init)
        ! Update the backend on the climate grid and fill `out`, then its
        ! elevation dependence and the modifiers (climate_mods). The geometry
        ! is the domain's, on the climate grid, with its present-day surface
        ! (z_srf_ref); the marine shelf lends its parameters
        ! to the esm ocean, which is interpolated to the shelf base here. idx is
        ! the driver's indices registry, advanced to the forcing time; snapesm
        ! (its axes), the modifiers and rembo ([anom]) read it.
        !
        ! snapesm and rembo are evaluated at the forcing time time_rel
        ! (with tstep_method = "const", time is the elapsed time and time_rel the
        ! fixed forcing time); esm at time. rembo equilibrates at the cold start.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
        real(wp),         intent(in) :: z_srf_ref(:,:)
        character(len=*), intent(in) :: domain
        real(wp),         intent(in) :: dx
        real(dp),         intent(in) :: dtt              ! [yr] time step of the run
        type(marshelf_class),  intent(in) :: mshlf
        type(indices_class),   intent(in) :: idx
        logical,               intent(in), optional :: init

        logical  :: is_init
        real(dp) :: time

        is_init = .false.
        if (present(init)) is_init = init

        time = ts%time_rel

        select case(trim(cl%method))
            case("snapesm")
                call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, idx=idx, dx=dx)
                call snapesm_state_copy(out%now, cl%snapesm%now)
                call snapesm_state_copy(out%ref, cl%snapesm%ref)

            case("esm")
                call esm_update(cl, out, ts%time, dtt, z_srf, H_ice, z_bed, f_grnd, z_sl, &
                                z_srf_ref, mshlf)

            case("rembo")
                call rembo_update(cl, out, ts, time, z_srf, H_ice, z_sl, domain, dx, &
                                  idx, is_init)
        end select

        ! The elevation dependence of the atmosphere, to move it to another surface
        ! (kryos_coupling): its surface z_srf, the monthly lapse rate and the
        ! precipitation sensitivity, as the backend applies them (anom reads
        ! beta_p). rembo: none.
        out%now%z_srf = z_srf
        if (.not. allocated(out%now%beta_p)) allocate(out%now%beta_p(size(z_srf,1),size(z_srf,2)))
        select case(trim(cl%method))
            case("snapesm")
                out%lapse       = snapesm_lapse_monthly(cl%snapesm)
                out%now%beta_p  = cl%snapesm%now%beta_p
            case("esm")
                out%lapse       = esm_lapse_monthly(cl%esm, cl%south)
                out%now%beta_p  = cl%esm%beta_p
            case default
                out%lapse       = 0.0_wp
                out%now%beta_p  = 0.0_wp
        end select

        call climate_mods_apply(cl%mods, out, idx)

    end subroutine climate_update

    subroutine snapesm_state_copy(dst, src)
        ! The fields of a snapesm state that climate_out carries, those it built:
        ! the atmosphere and/or the ocean of the reference snapshot.
        type(clim_state_class),    intent(inout) :: dst
        type(climate_state_class), intent(in)    :: src

        call copy_allocated(dst%tas,     src%tas)
        call copy_allocated(dst%pr,      src%pr)
        call copy_allocated(dst%tsl_ann, src%tsl_ann)
        call copy_allocated(dst%ta_ann,  src%ta_ann)
        call copy_allocated(dst%pr_ann,  src%pr_ann)
        call copy_allocated(dst%to_ann,  src%to_ann)
        call copy_allocated(dst%so_ann,  src%so_ann)
        call copy_allocated(dst%depth,   src%depth)

    end subroutine snapesm_state_copy

    subroutine copy_allocated_1D(dst, src)
        real(wp), allocatable, intent(inout) :: dst(:)
        real(wp), allocatable, intent(in)    :: src(:)
        if (allocated(src)) dst = src
    end subroutine copy_allocated_1D

    subroutine copy_allocated_2D(dst, src)
        real(wp), allocatable, intent(inout) :: dst(:,:)
        real(wp), allocatable, intent(in)    :: src(:,:)
        if (allocated(src)) dst = src
    end subroutine copy_allocated_2D

    subroutine copy_allocated_3D(dst, src)
        real(wp), allocatable, intent(inout) :: dst(:,:,:)
        real(wp), allocatable, intent(in)    :: src(:,:,:)
        if (allocated(src)) dst = src
    end subroutine copy_allocated_3D

    subroutine climate_end(cl)
        ! Teardown of the backend (snapesm, also the ocean of rembo); esm has
        ! none.
        type(yelmox_climate_class), intent(inout) :: cl

        select case(trim(cl%method))
            case("snapesm", "rembo")
                call snapesm_end(cl%snapesm)
        end select

    end subroutine climate_end

    ! ===== output ==========================================================

    function climate_file_base(cl) result(base)
        ! Base name of the climate's output files: <base>.nc (2D, on the climate
        ! grid) and, for esm and rembo, <base>_ts.nc (1D).
        type(yelmox_climate_class), intent(in) :: cl
        character(len=56) :: base

        select case(trim(cl%method))
            case("esm")
                base = "esm"
            case("rembo")
                base = "rembo"
            case default
                base = "snap"
        end select

    end function climate_file_base

    subroutine climate_write_init(cl, out, filename)
        ! The axes of the climate's 2D output file beyond xc, yc and time (created
        ! by the caller): the ocean depth axis, when the climate has ocean
        ! profiles, and the backend's own axes; then the modifiers' static fields.
        type(yelmox_climate_class), intent(in) :: cl
        type(climate_out_class),    intent(in) :: out
        character(len=*),           intent(in) :: filename

        if (allocated(out%now%depth)) &
            call nc_write_dim(filename, "depth", x=out%now%depth, units="m")

        select case(trim(cl%method))
            case("snapesm")
                call snapesm_write_init(cl%snapesm, filename)
        end select

        call climate_mods_write_init(cl%mods, filename)

    end subroutine climate_write_init

    subroutine climate_write_2D(cl, out, filename, ncid, n)
        ! Write one record of the climate (record n of an open 2D file on the
        ! climate grid), as applied (after the modifiers): the annual and,
        ! with monthly fields, summer near-surface air temperature, the annual
        ! precipitation, the annual ocean profiles (when the climate has them),
        ! plus the backend's own fields and the modifiers' diagnostics.
        type(yelmox_climate_class), intent(in) :: cl
        type(climate_out_class),    intent(in) :: out
        character(len=*),           intent(in) :: filename
        integer,                    intent(in) :: ncid, n

        call write2D("t2m_ann", out%now%ta_ann, "K", "Near-surface air temperature (ann)")
        if (allocated(out%now%tas)) &
            call write2D("t2m_sum", summer_mean(out%now%tas, cl%south), "K", &
                         "Near-surface air temperature (sum)")
        if (allocated(out%now%pr_ann)) &
            call write2D("pr_ann", out%now%pr_ann, "mm/a", "Precipitation (ann)")
        if (out%has_smb) then
            ! The climate's own surface fields (with surface_method = "climate",
            ! those landed on Yelmo).
            call write2D("smb_ann",  1e-3_wp*out%now%smb, "m/a water equiv.", "Surface mass balance (ann)")
            call write2D("tsrf_ann", out%now%tsrf, "K", "Surface temperature (ann)")
        end if
        if (allocated(out%now%to_ann)) &
            call write3D("to_ann", out%now%to_ann, "depth", "K", "Ocean temperature (ann)")
        if (allocated(out%now%so_ann)) &
            call write3D("so_ann", out%now%so_ann, "depth", "PSU", "Ocean salinity (ann)")

        select case(trim(cl%method))
            case("snapesm")
                ! The monthly atmosphere, then snapesm's own state (the axes).
                call write3D("tas",     out%now%tas,     "month", "K",    "Near-surface air temperature")
                call write3D("pr",      out%now%pr,      "month", "mm/d", "Precipitation")
                call write2D("tsl_ann", out%now%tsl_ann, "K", "Sea-level air temperature (ann)")
                call snapesm_write_step(cl%snapesm, filename, ncid, n)
            case("esm")
                call write_esm()
            case("rembo")
                ! Annual fields only: REMBO's own summer temperature.
                call write2D("t2m_sum", cl%rembo%ta_sum, "K", "Near-surface air temperature (sum)")
        end select

        call climate_mods_write_step(cl%mods, filename, ncid, n)

    contains

        subroutine write_esm()

            associate(esm => cl%esm)
            call write2D("dts_ann", sum(esm%dts, dim=3)/12.0_wp, "K", "Surface air temperature anomaly")
            if (cl%esm_ctl%use_smb) then
                call write2D("smb_ref_ann", 1e-3_wp*esm%smb_ann, "m/a water equiv.", "Reference SMB (ann)")
                call write2D("dsmb_ann", 1e-3_wp*sum(esm%dsmb, dim=3)/12.0_wp, "m/a water equiv.", &
                             "SMB anomaly (ann)")
                call write2D("dsmbdz",   1e-3_wp*esm%dsmbdz, "m/a m-1 water equiv.", "SMB lapse rate")
            else
                call write2D("dpr_ann", sum(esm%dpr, dim=3)/12.0_wp, "%", "Precipitation anomaly (ann)")
            end if
            call write2D("dto",    esm%dto,    "K",    "Shelf temperature anomaly")
            call write2D("dso",    esm%dso,    "PSU",  "Shelf salinity anomaly")
            call write2D("Qd_ann", esm%Qd_ann, "m3/s", "Subglacial discharge (ann)")
            call write2D("Qd_sum", esm%Qd_sum, "m3/s", "Subglacial discharge (sum)")
            end associate

        end subroutine write_esm

        subroutine write2D(vnm, var, units, long_name)
            character(len=*), intent(in) :: vnm, units, long_name
            real(wp),         intent(in) :: var(:,:)
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                          units=units, long_name=long_name)
        end subroutine write2D

        subroutine write3D(vnm, var, dim3, units, long_name)
            character(len=*), intent(in) :: vnm, dim3, units, long_name
            real(wp),         intent(in) :: var(:,:,:)
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3=dim3, dim4="time", &
                          start=[1,1,1,n], count=[size(var,1),size(var,2),size(var,3),1], &
                          ncid=ncid, units=units, long_name=long_name)
        end subroutine write3D

    end subroutine climate_write_2D

    subroutine climate_write_1D(cl, out, filename, time, H_ice, f_grnd, dx, init)
        ! The backend's 1D diagnostics (esm, rembo), with the geometry on the
        ! climate grid (cell size dx [m]); init creates the file.
        type(yelmox_climate_class), intent(in) :: cl
        type(climate_out_class),    intent(in) :: out
        character(len=*),           intent(in) :: filename
        real(dp),                   intent(in) :: time
        real(wp),                   intent(in) :: H_ice(:,:), f_grnd(:,:)
        real(wp),                   intent(in) :: dx
        logical,                    intent(in) :: init

        integer :: ncid, n

        if (trim(cl%method) /= "esm" .and. trim(cl%method) /= "rembo") return

        if (init) then
            call nc_create(filename)
            call nc_write_dim(filename, "time", x=time, dx=1.0_dp, nx=1, units="year", unlimited=.TRUE.)
        end if

        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        select case(trim(cl%method))
            case("esm")
                call write_esm()
            case("rembo")
                call write_rembo()
        end select
        call nc_close(ncid)

    contains

        subroutine write_esm()
            ! Means of the forcing over the ice (atmosphere) and the floating ice
            ! (ocean, at the shelf base).

            ! Missing value of the ocean means: a field with no valid floating cell.
            real(wp), parameter :: mv_ocn = -9999.0_wp
            logical,  allocatable :: mask_tot(:,:), mask_flt(:,:)
            integer  :: npts_tot
            real(wp) :: t2m_1d, pr_1d, dt_1d, dpr_1d, dt_var_1d, dpr_var_1d
            real(wp) :: to_1d, so_1d, dto_1d, dso_1d, dto_var_1d, dso_var_1d

            associate(esm => cl%esm)

            mask_tot = (H_ice > 0.0_wp)
            mask_flt = (H_ice > 0.0_wp .and. f_grnd == 0.0_wp)
            npts_tot = count(mask_tot)

            if (npts_tot > 0) then
                t2m_1d     = sum(esm%t2m_ann + sum(esm%dts, dim=3)/12.0_wp,  mask=mask_tot) / npts_tot
                pr_1d      = sum(sum(esm%pr*esm%dpr, dim=3)/12.0_wp,          mask=mask_tot) / npts_tot
                dt_1d      = sum(sum(esm%dts, dim=3)/12.0_wp,                 mask=mask_tot) / npts_tot
                dpr_1d     = sum(100.0_wp * sum(esm%dpr, dim=3)/12.0_wp,      mask=mask_tot) / npts_tot
                dt_var_1d  = sum(sum(esm%dts_var, dim=3)/12.0_wp,             mask=mask_tot) / npts_tot
                dpr_var_1d = sum(100.0_wp * sum(esm%dpr_var, dim=3)/12.0_wp,  mask=mask_tot) / npts_tot
            else
                t2m_1d = 0.0_wp; pr_1d = 1.0_wp; dt_1d = 0.0_wp; dpr_1d = 1.0_wp
                dt_var_1d = 0.0_wp; dpr_var_1d = 0.0_wp
            end if

            ! Per field over the floating cells with a value inside a physical range,
            ! so fill values do not enter the mean.
            to_1d      = ocn_mean(out%now%T_shlf, mask_flt,  240.0_wp, 320.0_wp, mv_ocn)
            so_1d      = ocn_mean(out%now%S_shlf, mask_flt,    0.0_wp,  60.0_wp, mv_ocn)
            dto_1d     = ocn_mean(esm%dto,        mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
            dso_1d     = ocn_mean(esm%dso,        mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
            dto_var_1d = ocn_mean(esm%dto_var,    mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
            dso_var_1d = ocn_mean(esm%dso_var,    mask_flt, -100.0_wp, 100.0_wp, mv_ocn)

            end associate

            call write1D("t2m_1d",     t2m_1d,     "K",      "Mean ice surf. Temp.")
            call write1D("pr_1d",      pr_1d,      "mm d-1", "Mean ice surf. Pr.")
            call write1D("dt_1d",      dt_1d,      "K",      "Mean ice surf. Temp. Anomaly")
            call write1D("dpr_1d",     dpr_1d,     "%",      "Mean ice surf. Pr. Anomaly")
            call write1D("dt_var_1d",  dt_var_1d,  "K",      "Mean ice surf. Temp. Anomaly (Variability)")
            call write1D("dpr_var_1d", dpr_var_1d, "%",      "Mean ice surf. Pr. Anomaly (Variability)")
            call write1D("to_1d",      to_1d,      "K",      "Mean ice-shelf draft Temp.")
            call write1D("so_1d",      so_1d,      "PSU",    "Mean ice-shelf draft Sal.")
            call write1D("dto_1d",     dto_1d,     "K",      "Mean ice-shelf Temp. Anomaly")
            call write1D("dso_1d",     dso_1d,     "PSU",    "Mean ice-shelf draft Sal. Anomaly")
            call write1D("dto_var_1d", dto_var_1d, "K",      "Mean ice-shelf Temp. Anomaly (Variability)")
            call write1D("dso_var_1d", dso_var_1d, "PSU",    "Mean ice-shelf Sal. Anomaly (Variability)")

        end subroutine write_esm

        subroutine write_rembo()
            ! The applied anomalies, the surface mass balance integrated over the
            ! ice and the accumulation-area ratio.

            integer  :: ntot
            real(wp) :: smb_tot, aar

            ntot = count(H_ice > 0.0_wp)
            if (ntot > 0) then
                aar     = real(count(H_ice > 0.0_wp .and. out%now%smb > 0.0_wp), wp) / real(ntot, wp)
                ! [m3/a] * [1000 kg/m3] * [1e-12 Gt/kg] = [Gt/a]
                smb_tot = dx**2 * sum(out%now%smb*1e-3_wp, mask=H_ice > 0.0_wp) * 1000.0_wp * 1e-12_wp
            else
                aar     = 0.0_wp
                smb_tot = 0.0_wp
            end if

            call write1D("dT_sum",   cl%rembo%dT_sum, "K",     "Summer air temperature anomaly")
            call write1D("dT_ann",   cl%rembo%dT_ann, "K",     "Annual air temperature anomaly")
            call write1D("dT_ocn",   cl%rembo%dT_ocn, "K",     "Ocean temperature anomaly")
            call write1D("smb_mean", smb_tot,         "Gt/yr", "Surface mass balance over the ice sheet")
            call write1D("aar",      aar,             "1",     "Accumulation area ratio")

        end subroutine write_rembo

        subroutine write1D(vnm, val, units, long_name)
            character(len=*), intent(in) :: vnm, units, long_name
            real(wp),         intent(in) :: val
            call nc_write(filename, vnm, val, dim1="time", start=[n], count=[1], ncid=ncid, &
                          units=units, long_name=long_name)
        end subroutine write1D

    end subroutine climate_write_1D

    function ocn_mean(field, base_mask, lo, hi, mv) result(val)
        ! Mean of `field` over the `base_mask` cells whose value is finite and
        ! inside (lo,hi); `mv` if no cell qualifies (NaN/Inf fail the bounds).
        real(wp), intent(in) :: field(:,:)
        logical,  intent(in) :: base_mask(:,:)
        real(wp), intent(in) :: lo, hi, mv
        real(wp) :: val

        logical, allocatable :: m(:,:)
        integer :: npts

        m    = base_mask .and. field > lo .and. field < hi
        npts = count(m)
        if (npts > 0) then
            val = sum(field, mask=m) / real(npts, wp)
        else
            val = mv
        end if

    end function ocn_mean

    ! ===== restart =========================================================

    subroutine climate_restart_write(cl, fldr, time, z_srf, H_ice, z_sl)
        ! The backend's own restart in the bundle folder `fldr`, with the
        ! geometry on the climate grid: REMBO's (rembo_restart.nc); snapesm
        ! and esm hold no state (snapesm's indices restart with the
        ! driver's registry). REMBO reads its restart itself, as set in its own
        ! parameter file.
        type(yelmox_climate_class), intent(in) :: cl
        character(len=*),           intent(in) :: fldr
        real(dp),                   intent(in) :: time
        real(wp),                   intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)

        select case(trim(cl%method))
            case("rembo")
                call rembo_clim_restart_write(trim(fldr)//"/rembo_restart.nc", time, z_srf, H_ice, z_sl)
        end select

    end subroutine climate_restart_write

    ! ===== esm backend =====================================================

    subroutine esm_init(cl, filename, domain, grid, group, timeline_group, use_smb, reg)
        ! Read [esm] (experiment + physics) and the esm periods of the run phase
        ! (timeline group), seed the random generator for the climate
        ! variability, and initialize esm_forcing on the climate grid.
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: filename, domain, group, timeline_group
        type(grid_class), intent(in) :: grid
        logical,          intent(in) :: use_smb
        type(regions_class), intent(in) :: reg           ! the regions on the climate grid

        integer :: n
        integer, allocatable :: seed(:)
        character(len=56) :: basins

        cl%esm_ctl%run_type = trim(timeline_group)
        cl%esm_ctl%use_smb  = use_smb

        call nml_read(filename, group, "par_file",     cl%esm_ctl%par_file)

        ! Basins of the ocean extrapolation: "<set>", "<set>.group" of the
        ! regions, "None" or "domain" (one basin)
        call nml_read(filename, group, "basins",       basins)
        allocate(cl%basins(grid%G%nx,grid%G%ny))
        cl%basins = real(regions_basin_ids(reg, basins), wp)
        call nml_read(filename, group, "experiment",   cl%esm_ctl%experiment)
        call nml_read(filename, group, "esm_name",     cl%esm_ctl%esm_name)
        call nml_read(filename, group, "use_esm",      cl%esm_ctl%use_esm)
        call nml_read(filename, group, "use_var",      cl%esm_ctl%use_var)
        call nml_read(filename, group, "use_proj",     cl%esm_ctl%use_proj)
        call nml_read(filename, group, "use_hist",     cl%esm_ctl%use_hist)
        call nml_read(filename, group, "lapse",        cl%esm%lapse)
        call nml_read(filename, group, "f_p",          cl%esm%beta_p)
        call nml_read(filename, group, "f_ocn",        cl%esm%f_ocn)
        call nml_read(filename, group, "f_polar",      cl%esm%f_polar)
        call nml_read(filename, group, "dT_threshold", cl%esm%dT_lim)

        call nml_read(filename, timeline_group, "time_ref",     cl%esm_ctl%time_ref)
        call nml_read(filename, timeline_group, "time_hist",    cl%esm_ctl%time_hist)
        call nml_read(filename, timeline_group, "time_proj",    cl%esm_ctl%time_proj)
        call nml_read(filename, timeline_group, "time_esm_ref", cl%esm_ctl%time_esm_ref)
        call nml_read(filename, timeline_group, "clim_var",     cl%esm_ctl%clim_var)
        call nml_read(filename, timeline_group, "clim_seed",    cl%esm_ctl%clim_seed)

        call random_seed(size=n)
        allocate(seed(n))
        seed = cl%esm_ctl%clim_seed
        call random_seed(put=seed)

        call esm_forcing_init(cl%esm, trim(cl%esm_ctl%par_file), domain, grid, &
                              run_type=cl%esm_ctl%run_type, gcm=cl%esm_ctl%esm_name, &
                              experiment=cl%esm_ctl%experiment, use_esm=cl%esm_ctl%use_esm, &
                              use_smb=cl%esm_ctl%use_smb, use_var=cl%esm_ctl%use_var, &
                              use_hist=cl%esm_ctl%use_hist, use_proj=cl%esm_ctl%use_proj)

    end subroutine esm_init

    subroutine esm_update(cl, out, time, dtt, z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref, &
                          mshlf)
        ! The reference climatology at the current surface, the esm anomalies
        ! (historical / projection / homogeneous) and the variability, then the
        ! products: atmosphere, the surface mass balance (surface_method = climate),
        ! the ocean at the shelf base and subglacial discharge.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        real(dp),         intent(in) :: time, dtt
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
        real(wp),         intent(in) :: z_srf_ref(:,:)
        type(marshelf_class), intent(in) :: mshlf

        real(wp), allocatable :: to_ref(:,:,:), so_ref(:,:,:)   ! the reference ocean, extrapolated

        associate(esm => cl%esm, ec => cl%esm_ctl, basins => cl%basins)

        call esm_clim_update(esm, z_srf, time, ec%time_ref, ec%use_smb, cl%south)

        ! Extrapolate the reference ocean into ice-shelf interiors.
        to_ref = esm%to_ref%var(:,:,:,1)
        so_ref = esm%so_ref%var(:,:,:,1)
        if (mshlf%par%extrap_shlf) then
            call ocn_variable_extrapolation(to_ref, H_ice, basins, -esm%to_ref%z, z_bed)
            call ocn_variable_extrapolation(so_ref, H_ice, basins, -esm%so_ref%z, z_bed)
        end if

        call esm_forcing_update(esm, mshlf, time, ec%use_esm, ec%time_ref, ec%time_hist, &
                                ec%time_proj, ec%time_esm_ref, H_ice, basins, z_bed, &
                                f_grnd, z_sl, ec%use_smb, use_ref_atm=.false., use_ref_ocn=.false.)

        call esm_variability_update(esm, mshlf, time, dtt, ec%clim_var, ec%time_ref, H_ice, &
                                    basins, z_bed, f_grnd, z_sl, ec%use_var, &
                                    use_ref_atm=.false., use_ref_ocn=.false.)

        ! Atmosphere.
        out%now%tas    = esm%t2m + esm%dts + esm%dts_var
        out%ref%tas    = esm%t2m
        out%now%ta_ann = sum(out%now%tas, dim=3) / 12.0_wp
        out%ref%ta_ann = esm%t2m_ann
        if (.not. ec%use_smb) then
            out%now%pr     = esm%pr * esm%dpr * esm%dpr_var
            out%ref%pr     = esm%pr
            out%now%pr_ann = sum(out%now%pr, dim=3) / 12.0_wp * 365.0_wp
            out%ref%pr_ann = esm%pr_ann * 365.0_wp
        end if

        ! Surface mass balance at the current surface: the reference smb and its
        ! anomaly (both at the present-day surface), corrected for the elevation
        ! change with the esm smb gradient. Surface temperature from the
        ! near-surface air, at most melting over ice.
        out%has_smb = ec%use_smb
        if (ec%use_smb) then
            out%now%smb = esm%smb_ann + sum(esm%dsmb, dim=3) / 12.0_wp &
                        - esm%dsmbdz*(z_srf_ref - z_srf)
        end if
        out%now%tsrf = sum(esm%t2m + esm%dts + esm%dts_var, dim=3) / 12.0_wp
        where (H_ice > 0.0_wp .and. out%now%tsrf > 273.15_wp) out%now%tsrf = 273.15_wp

        ! Ocean at the shelf base: the reference ocean interpolated to the shelf
        ! base, plus the esm anomalies (themselves at the shelf base).
        out%has_ocn_shelf = .true.
        if (.not. allocated(out%now%T_shlf)) then
            allocate(out%now%T_shlf(size(H_ice,1), size(H_ice,2)))
            allocate(out%now%S_shlf(size(H_ice,1), size(H_ice,2)))
        end if
        call marshelf_interp_shelf(out%now%T_shlf, mshlf, to_ref, &
                                   H_ice, z_bed, f_grnd, z_sl, -esm%to_ref%z)
        call marshelf_interp_shelf(out%now%S_shlf, mshlf, so_ref, &
                                   H_ice, z_bed, f_grnd, z_sl, -esm%so_ref%z)
        out%now%T_shlf  = out%now%T_shlf + esm%dto + esm%dto_var
        out%now%S_shlf  = out%now%S_shlf + esm%dso + esm%dso_var
        out%now%dT_shlf = esm%dto + esm%dto_var
        out%now%dS_shlf = esm%dso + esm%dso_var

        ! Subglacial discharge.
        out%has_Qd = .true.
        out%now%Qd = esm%Qd_ann

        end associate

    end subroutine esm_update

    ! ===== rembo backend ===================================================

    subroutine rembo_update(cl, out, ts, time, z_srf, H_ice, z_sl, domain, dx, &
                            idx, init)
        ! REMBO's atmosphere and surface mass balance, and the ocean from snapesm
        ! (at `time`, as for the snapesm backend). The exception to the
        ! modifiers: REMBO needs the anomaly before it runs, so it reads [anom]
        ! itself (climate_mods skips anom for rembo). The index f maps onto
        ! REMBO's anomalies: the summer air temperature dT_sum = f*f_ta, the
        ! annual dT_ann = 1.3*dT_sum (REMBO's winter factor 1.6), and the ocean
        ! dT_ocn = dT_ann*f_to, added to the snapesm ocean; f_so is unused.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(dp),         intent(in) :: time
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)
        character(len=*), intent(in) :: domain
        real(wp),         intent(in) :: dx
        type(indices_class), intent(in) :: idx
        logical,          intent(in) :: init

        integer :: nx, ny

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        if (cl%mods%anom%active) then
            cl%rembo%dT_sum = indices_value(idx, cl%mods%anom%index) * cl%mods%anom%f_ta
            cl%rembo%dT_ann = 1.3_wp * cl%rembo%dT_sum
            cl%rembo%dT_ocn = cl%rembo%dT_ann * cl%mods%anom%f_to
        else
            cl%rembo%dT_sum = 0.0_wp
            cl%rembo%dT_ann = 0.0_wp
            cl%rembo%dT_ocn = 0.0_wp
        end if

        ! Atmosphere and surface mass balance, at the current surface.
        if (.not. allocated(out%now%ta_ann)) then
            allocate(out%now%ta_ann(nx,ny), out%now%pr_ann(nx,ny))
            allocate(out%now%smb(nx,ny),    out%now%tsrf(nx,ny))
        end if
        call rembo_clim_update(ts%time, ts%time_rel, cl%rembo%dT_sum, z_srf, H_ice, z_sl, init, &
                               out%now%ta_ann, cl%rembo%ta_sum, out%now%pr_ann, &
                               out%now%smb, out%now%tsrf)
        out%has_smb = .true.

        ! Ocean (depth profiles) from snapesm (ocean only), plus dT_ocn.
        call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, idx=idx, dx=dx)
        call snapesm_state_copy(out%now, cl%snapesm%now)
        call snapesm_state_copy(out%ref, cl%snapesm%ref)
        out%now%to_ann = out%now%to_ann + cl%rembo%dT_ocn

    end subroutine rembo_update

end module yelmox_climate
