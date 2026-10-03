module yelmox_climate
    ! Climate backend of a domain, chosen at runtime ([coupling] climate):
    !   "snapclim" -- snapshot/anomaly climate (snapclim)
    !   "snapesm"  -- snapshots blended by indices (snapesm)
    !   "esm"      -- reference climatology plus Earth-system-model anomalies over
    !                 historical/projection periods (esm_forcing)
    !   "rembo"    -- REMBO atmosphere and surface mass balance (climate_rembo,
    !                 built with rembo=1), with the ocean from snapclim
    !
    ! Presents one interface (climate_init / climate_update) over the backends,
    ! filling a `climate_out_class`; the domain and its coupling read only that.
    ! The driver's transient forcing (tsforcing: f_now and the gains f_ta, f_to,
    ! f_so) is handed to the backend, which applies it in its own way.

    use precision,     only : wp
    use nml,           only : nml_read
    use ncio
    use timestepping,  only : tstep_class
    use climate_out,   only : climate_out_class
    use snapclim,      only : snapclim_class, snapclim_init, snapclim_update, snapclim_air_anom
    use snapesm,       only : snapesm_class, snapesm_init, snapesm_update
    use esm_forcing,   only : esm_forcing_class, esm_forcing_init, esm_clim_update, &
                              esm_forcing_update, esm_variability_update, esm_summer_mean
    use marine_shelf,  only : marshelf_class, marshelf_interp_shelf, ocn_variable_extrapolation
    use kryos_forcing, only : tsforcing_class
    use climate_rembo, only : rembo_clim_init, rembo_clim_update, rembo_clim_restart_write

    implicit none

    private

    ! Run control of the esm backend: [esm] holds the experiment, the timeline
    ! group of the run phase its periods.
    type esm_ctl_class
        character(len=56)  :: run_type        ! the timeline group, e.g. spinup | transient
        character(len=512) :: par_file        ! esm data configuration
        character(len=56)  :: experiment, esm_name
        logical            :: use_esm, use_smb, use_var, use_hist, use_proj
        real(wp)           :: time_ref(2), time_hist(2), time_proj(2), time_esm_ref(2)
        character(len=56)  :: clim_var
        integer            :: clim_seed
    end type esm_ctl_class

    ! State of the rembo backend beyond the climate products: REMBO's summer
    ! air temperature and the anomalies applied from the transient forcing.
    type rembo_state_class
        real(wp), allocatable :: ta_sum(:,:)     ! [K] near-surface air temperature (summer)
        real(wp) :: dT_sum = 0.0_wp              ! [K] summer air temperature anomaly
        real(wp) :: dT_ann = 0.0_wp              ! [K] annual air temperature anomaly
        real(wp) :: dT_ocn = 0.0_wp              ! [K] ocean temperature anomaly
    end type rembo_state_class

    type yelmox_climate_class
        character(len=16)       :: method = ""   ! snapclim | snapesm | esm | rembo
        character(len=256)      :: grid_name     ! the climate grid (grid_clim)
        logical                 :: south = .false.   ! southern hemisphere
        type(snapclim_class)    :: snapclim
        type(snapesm_class)     :: snapesm
        type(esm_forcing_class) :: esm
        type(esm_ctl_class)     :: esm_ctl
        type(rembo_state_class) :: rembo
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update
    public :: climate_air_anom
    public :: climate_file_base, climate_write_2D, climate_write_1D
    public :: climate_restart_write

contains

    subroutine climate_init(cl, method, filename, domain, grid_name, nx, ny, time, basins, &
                            south, sfx, timeline_group, smb_direct)
        ! south: the domain lies in the southern hemisphere (seasons, lapse rates).
        ! smb_direct: the surface mass balance is taken from the climate
        ! ([coupling] smb_method = "climate"); the esm backend supplies it, and
        ! the rembo backend supplies nothing else (annual fields only).
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: method
        character(len=*), intent(in) :: filename, domain, grid_name
        integer,          intent(in) :: nx, ny
        real(wp),         intent(in) :: time
        real(wp),         intent(in) :: basins(:,:)
        logical,          intent(in) :: south
        character(len=*), intent(in) :: sfx              ! namelist group suffix of the domain
        character(len=*), intent(in) :: timeline_group   ! group of the run phase's timeline
        logical,          intent(in) :: smb_direct

        cl%method    = trim(method)
        cl%grid_name = trim(grid_name)
        cl%south     = south

        if (smb_direct .and. trim(cl%method) /= "esm" .and. trim(cl%method) /= "rembo") then
            write(*,*) "climate_init:: error: smb_method = climate needs a climate that supplies &
                       &the surface mass balance (esm, rembo); got climate = ", trim(cl%method)
            error stop 1
        end if
        if (.not. smb_direct .and. trim(cl%method) == "rembo") then
            write(*,*) "climate_init:: error: climate = rembo supplies the surface mass balance &
                       &only; set smb_method = climate."
            error stop 1
        end if

        select case(trim(cl%method))
            case("snapclim")
                call snapclim_init(cl%snapclim, filename, domain, grid_name, nx, ny, basins, &
                                   south, group="snap"//trim(sfx))
            case("snapesm")
                call snapesm_init(cl%snapesm, filename, domain, grid_name, nx, ny, time, basins, &
                                  south, group="snap"//trim(sfx))
            case("esm")
                call esm_init(cl, filename, domain, grid_name, "esm"//trim(sfx), timeline_group, &
                              smb_direct)
            case("rembo")
                ! REMBO for the atmosphere and smb, snapclim for the ocean.
                call snapclim_init(cl%snapclim, filename, domain, grid_name, nx, ny, basins, &
                                   south, group="snap"//trim(sfx))
                call rembo_clim_init(time, nx, ny)
                allocate(cl%rembo%ta_sum(nx,ny))
                cl%rembo%ta_sum = 0.0_wp
            case default
                write(*,*) "climate_init:: error: climate must be snapclim, snapesm, esm or rembo; got ", &
                           trim(cl%method)
                error stop 1
        end select

    end subroutine climate_init

    subroutine climate_update(cl, out, ts, z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref, basins, &
                              domain, dx, dtt, mshlf, tsf, init)
        ! Update the backend on the climate grid and fill `out`. The geometry is
        ! the domain's, on the climate grid, with its present-day surface
        ! (z_srf_ref); the marine shelf lends its parameters to the esm ocean,
        ! which is interpolated to the shelf base here.
        !
        ! With an active transient forcing, its spatially homogeneous anomalies
        ! (dTa = f_now*f_ta, dTo = f_now*f_to, dSo = f_now*f_so) go to the backend:
        ! snapclim uses them in place of its own index in its "anom" modes;
        ! snapesm adds them on top in every mode; esm has its own forcing; rembo
        ! maps them its own way (see rembo_update).
        !
        ! The cold start (init) evaluates snapclim/snapesm at time_rel, the time
        ! loop at time; esm always at time. rembo equilibrates at the cold start.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
        real(wp),         intent(in) :: z_srf_ref(:,:)
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in) :: domain
        real(wp),         intent(in) :: dx
        real(wp),         intent(in) :: dtt              ! [yr] time step of the run
        type(marshelf_class),  intent(in) :: mshlf
        type(tsforcing_class), intent(in), optional :: tsf
        logical,               intent(in), optional :: init

        logical  :: forced, is_init
        real(wp) :: time

        forced = .false.
        if (present(tsf)) forced = tsf%active

        is_init = .false.
        if (present(init)) is_init = init

        time = ts%time
        if (is_init) time = ts%time_rel

        select case(trim(cl%method))
            case("snapclim")
                if (forced) then
                    call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, &
                                         dTa=tsf%dTa, dTo=tsf%dTo, dSo=tsf%dSo, dx=dx, basins=basins)
                else
                    call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, &
                                         dx=dx, basins=basins)
                end if

                ! snapclim's reference climate is clim0.
                out%now%tas     = cl%snapclim%now%tas
                out%now%pr      = cl%snapclim%now%pr
                out%now%tsl_ann = cl%snapclim%now%tsl_ann
                out%now%ta_ann  = cl%snapclim%now%ta_ann
                out%now%pr_ann  = cl%snapclim%now%pr_ann
                out%now%to_ann  = cl%snapclim%now%to_ann
                out%now%so_ann  = cl%snapclim%now%so_ann
                out%now%depth   = cl%snapclim%now%depth

                out%ref%tas     = cl%snapclim%clim0%tas
                out%ref%pr      = cl%snapclim%clim0%pr
                out%ref%tsl_ann = cl%snapclim%clim0%tsl_ann
                out%ref%ta_ann  = cl%snapclim%clim0%ta_ann
                out%ref%pr_ann  = cl%snapclim%clim0%pr_ann
                out%ref%to_ann  = cl%snapclim%clim0%to_ann
                out%ref%so_ann  = cl%snapclim%clim0%so_ann
                out%ref%depth   = cl%snapclim%clim0%depth

            case("snapesm")
                if (forced) then
                    call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, &
                                        dTa=tsf%dTa, dTo=tsf%dTo, dSo=tsf%dSo, dx=dx, basins=basins)
                else
                    call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, &
                                        dx=dx, basins=basins)
                end if

                out%now%tas     = cl%snapesm%now%tas
                out%now%pr      = cl%snapesm%now%pr
                out%now%tsl_ann = cl%snapesm%now%tsl_ann
                out%now%ta_ann  = cl%snapesm%now%ta_ann
                out%now%pr_ann  = cl%snapesm%now%pr_ann
                out%now%to_ann  = cl%snapesm%now%to_ann
                out%now%so_ann  = cl%snapesm%now%so_ann
                out%now%depth   = cl%snapesm%now%depth

                out%ref%tas     = cl%snapesm%ref%tas
                out%ref%pr      = cl%snapesm%ref%pr
                out%ref%tsl_ann = cl%snapesm%ref%tsl_ann
                out%ref%ta_ann  = cl%snapesm%ref%ta_ann
                out%ref%pr_ann  = cl%snapesm%ref%pr_ann
                out%ref%to_ann  = cl%snapesm%ref%to_ann
                out%ref%so_ann  = cl%snapesm%ref%so_ann
                out%ref%depth   = cl%snapesm%ref%depth

            case("esm")
                call esm_update(cl, out, ts%time, dtt, z_srf, H_ice, z_bed, f_grnd, z_sl, &
                                z_srf_ref, basins, domain, mshlf)

            case("rembo")
                call rembo_update(cl, out, ts, time, z_srf, H_ice, z_sl, basins, domain, dx, &
                                  forced, is_init, tsf)
        end select

    end subroutine climate_update

    function climate_air_anom(cl, time) result(dT)
        ! [K] The backend's air-temperature anomaly index at `time`, the forcing of
        ! the bipolar ocean box model: snapclim's index series `at` scaled by
        ! dTa_const. snapesm and esm have no such index.
        type(yelmox_climate_class), intent(in) :: cl
        real(wp),                   intent(in) :: time
        real(wp) :: dT

        select case(trim(cl%method))
            case("snapclim")
                dT = snapclim_air_anom(cl%snapclim, time)
            case default
                write(*,*) "climate_air_anom:: error: climate = ", trim(cl%method), &
                           " has no air-temperature anomaly index (needed by the ocean box model)."
                error stop 1
        end select

    end function climate_air_anom

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

    subroutine climate_write_2D(cl, out, filename, ncid, n)
        ! Write one record of the climate (record n of an open 2D file on the
        ! climate grid): the annual near-surface air temperature and
        ! precipitation, plus the backend's own fields.
        type(yelmox_climate_class), intent(in) :: cl
        type(climate_out_class),    intent(in) :: out
        character(len=*),           intent(in) :: filename
        integer,                    intent(in) :: ncid, n

        call write2D("t2m_ann", out%now%ta_ann, "K", "Near-surface air temperature (ann)")
        if (allocated(out%now%pr_ann)) &
            call write2D("pr_ann", out%now%pr_ann, "mm/a", "Precipitation (ann)")

        select case(trim(cl%method))
            case("esm")
                call write_esm()
            case("rembo")
                call write2D("t2m_sum", cl%rembo%ta_sum, "K", "Near-surface air temperature (sum)")
                call write2D("smb_ann", 1e-3_wp*out%now%smb, "m/a water equiv.", &
                             "Surface mass balance (ann)")
        end select

    contains

        subroutine write_esm()

            associate(esm => cl%esm)
            call write2D("t2m_sum", esm%t2m_sum + esm_summer_mean(esm%dts, cl%south), "K", &
                         "Near-surface air temperature (sum)")
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

    end subroutine climate_write_2D

    subroutine climate_write_1D(cl, out, filename, time, H_ice, f_grnd, dx, init)
        ! The backend's 1D diagnostics (esm, rembo), with the geometry on the
        ! climate grid (cell size dx [m]); init creates the file.
        type(yelmox_climate_class), intent(in) :: cl
        type(climate_out_class),    intent(in) :: out
        character(len=*),           intent(in) :: filename
        real(wp),                   intent(in) :: time
        real(wp),                   intent(in) :: H_ice(:,:), f_grnd(:,:)
        real(wp),                   intent(in) :: dx
        logical,                    intent(in) :: init

        integer :: ncid, n

        if (trim(cl%method) /= "esm" .and. trim(cl%method) /= "rembo") return

        if (init) then
            call nc_create(filename)
            call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, units="year", unlimited=.TRUE.)
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
        ! geometry on the climate grid: REMBO's (rembo_restart.nc); the other
        ! backends hold no state.
        type(yelmox_climate_class), intent(in) :: cl
        character(len=*),           intent(in) :: fldr
        real(wp),                   intent(in) :: time
        real(wp),                   intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)

        if (trim(cl%method) == "rembo") &
            call rembo_clim_restart_write(trim(fldr)//"/rembo_restart.nc", time, z_srf, H_ice, z_sl)

    end subroutine climate_restart_write

    ! ===== esm backend =====================================================

    subroutine esm_init(cl, filename, domain, grid_name, group, timeline_group, use_smb)
        ! Read [esm] (experiment + physics) and the esm periods of the run phase
        ! (timeline group), seed the random generator for the climate
        ! variability, and initialize esm_forcing on the climate grid.
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: filename, domain, grid_name, group, timeline_group
        logical,          intent(in) :: use_smb

        integer :: n
        integer, allocatable :: seed(:)

        cl%esm_ctl%run_type = trim(timeline_group)
        cl%esm_ctl%use_smb  = use_smb

        call nml_read(filename, group, "par_file",     cl%esm_ctl%par_file)
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
        call nml_read(filename, group, "grid_src",     cl%esm%grid_src)

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

        call esm_forcing_init(cl%esm, trim(cl%esm_ctl%par_file), domain, grid_name, &
                              run_type=cl%esm_ctl%run_type, gcm=cl%esm_ctl%esm_name, &
                              experiment=cl%esm_ctl%experiment, use_esm=cl%esm_ctl%use_esm, &
                              use_smb=cl%esm_ctl%use_smb, use_var=cl%esm_ctl%use_var, &
                              use_hist=cl%esm_ctl%use_hist, use_proj=cl%esm_ctl%use_proj)

    end subroutine esm_init

    subroutine esm_update(cl, out, time, dtt, z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref, &
                          basins, domain, mshlf)
        ! The reference climatology at the current surface, the esm anomalies
        ! (historical / projection / homogeneous) and the variability, then the
        ! products: atmosphere, the surface mass balance (smb_method = climate),
        ! the ocean at the shelf base and subglacial discharge.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        real(wp),         intent(in) :: time, dtt
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
        real(wp),         intent(in) :: z_srf_ref(:,:)
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in) :: domain
        type(marshelf_class), intent(in) :: mshlf

        associate(esm => cl%esm, ec => cl%esm_ctl)

        call esm_clim_update(esm, z_srf, time, ec%time_ref, ec%use_smb, cl%south)

        ! Extrapolate the reference ocean into ice-shelf interiors.
        if (mshlf%par%extrap_shlf) then
            call ocn_variable_extrapolation(esm%to_ref%var(:,:,:,1), H_ice, basins, -esm%to_ref%z, z_bed)
            call ocn_variable_extrapolation(esm%so_ref%var(:,:,:,1), H_ice, basins, -esm%so_ref%z, z_bed)
        end if

        call esm_forcing_update(esm, mshlf, time, ec%use_esm, ec%time_ref, ec%time_hist, &
                                ec%time_proj, ec%time_esm_ref, domain, H_ice, basins, z_bed, &
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
        call marshelf_interp_shelf(out%now%T_shlf, mshlf, esm%to_ref%var(:,:,:,1), &
                                   H_ice, z_bed, f_grnd, z_sl, -esm%to_ref%z)
        call marshelf_interp_shelf(out%now%S_shlf, mshlf, esm%so_ref%var(:,:,:,1), &
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

    subroutine rembo_update(cl, out, ts, time, z_srf, H_ice, z_sl, basins, domain, dx, &
                            forced, init, tsf)
        ! REMBO's atmosphere and surface mass balance, and the ocean from snapclim
        ! (at `time`, as for the snapclim backend). The transient forcing maps
        ! onto REMBO's anomalies: the summer air temperature dT_sum = f_now*f_ta,
        ! the annual dT_ann = 1.3*dT_sum (REMBO's winter factor 1.6), and the
        ! ocean dT_ocn = dT_ann*f_to, added to a snapclim ocean held at its
        ! reference (ocn_type = "const").
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(wp),         intent(in) :: time
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in) :: domain
        real(wp),         intent(in) :: dx
        logical,          intent(in) :: forced, init
        type(tsforcing_class), intent(in), optional :: tsf

        integer :: nx, ny

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        if (forced) then
            cl%rembo%dT_sum = tsf%dTa
            cl%rembo%dT_ann = 1.3_wp * cl%rembo%dT_sum
            cl%rembo%dT_ocn = cl%rembo%dT_ann * tsf%f_to
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

        ! Ocean (depth profiles) from snapclim; its reference is clim0.
        call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, &
                             dx=dx, basins=basins)
        out%now%to_ann = cl%snapclim%now%to_ann
        out%now%so_ann = cl%snapclim%now%so_ann
        out%now%depth  = cl%snapclim%now%depth
        out%ref%to_ann = cl%snapclim%clim0%to_ann
        out%ref%so_ann = cl%snapclim%clim0%so_ann
        out%ref%depth  = cl%snapclim%clim0%depth
        if (trim(cl%snapclim%par%ocn_type) == "const") &
            out%now%to_ann = out%now%to_ann + cl%rembo%dT_ocn

    end subroutine rembo_update

end module yelmox_climate
