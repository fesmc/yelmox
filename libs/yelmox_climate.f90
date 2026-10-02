module yelmox_climate
    ! Climate backend of a domain, chosen at runtime ([coupling] climate):
    !   "snapclim" -- snapshot/anomaly climate (snapclim)
    !   "snapesm"  -- snapshots blended by indices (snapesm)
    !   "esm"      -- reference climatology plus Earth-system-model anomalies over
    !                 historical/projection periods (esm_forcing)
    !
    ! Presents one interface (climate_init / climate_update) over the backends,
    ! filling a `climate_out_class`; the domain and its coupling read only that.
    ! The driver's transient forcing (tsforcing: f_now and the gains f_ta, f_to,
    ! f_so) is handed to the backend, which applies it in its own way.

    use precision,     only : wp
    use nml,           only : nml_read
    use timestepping,  only : tstep_class
    use climate_out,   only : climate_out_class
    use snapclim,      only : snapclim_class, snapclim_init, snapclim_update, snapclim_air_anom
    use snapesm,       only : snapesm_class, snapesm_init, snapesm_update
    use esm_forcing,   only : esm_forcing_class, esm_forcing_init, esm_clim_update, &
                              esm_forcing_update, esm_variability_update
    use marine_shelf,  only : marshelf_class, marshelf_interp_shelf, ocn_variable_extrapolation
    use kryos_forcing, only : tsforcing_class

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

    type yelmox_climate_class
        character(len=16)       :: method = ""   ! snapclim | snapesm | esm
        character(len=256)      :: grid_name     ! the climate grid (grid_clim)
        type(snapclim_class)    :: snapclim
        type(snapesm_class)     :: snapesm
        type(esm_forcing_class) :: esm
        type(esm_ctl_class)     :: esm_ctl
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update
    public :: climate_air_anom
    public :: climate_ocean_const

contains

    subroutine climate_init(cl, method, filename, domain, grid_name, nx, ny, time, basins, &
                            sfx, timeline_group, smb_direct)
        ! smb_direct: the surface mass balance is taken from the climate
        ! ([coupling] smb_method = "climate"); only the esm backend supplies it.
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: method
        character(len=*), intent(in) :: filename, domain, grid_name
        integer,          intent(in) :: nx, ny
        real(wp),         intent(in) :: time
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in) :: sfx              ! namelist group suffix of the domain
        character(len=*), intent(in) :: timeline_group   ! group of the run phase's timeline
        logical,          intent(in) :: smb_direct

        cl%method    = trim(method)
        cl%grid_name = trim(grid_name)

        if (smb_direct .and. trim(cl%method) /= "esm") then
            write(*,*) "climate_init:: error: smb_method = climate needs a climate that supplies &
                       &the surface mass balance (esm); got climate = ", trim(cl%method)
            error stop 1
        end if

        select case(trim(cl%method))
            case("snapclim")
                call snapclim_init(cl%snapclim, filename, domain, grid_name, nx, ny, basins, &
                                   group="snap"//trim(sfx))
            case("snapesm")
                call snapesm_init(cl%snapesm, filename, domain, grid_name, nx, ny, time, basins, &
                                  group="snap"//trim(sfx))
            case("esm")
                call esm_init(cl, filename, domain, grid_name, "esm"//trim(sfx), timeline_group, &
                              smb_direct)
            case default
                write(*,*) "climate_init:: error: climate must be snapclim, snapesm or esm; got ", &
                           trim(cl%method)
                error stop 1
        end select

    end subroutine climate_init

    subroutine climate_update(cl, out, ts, z_srf, H_ice, z_bed, f_grnd, z_sl, basins, domain, &
                              dx, dtt, mshlf, tsf, init)
        ! Update the backend on the climate grid and fill `out`. The geometry is
        ! the domain's, on the climate grid; the marine shelf lends its parameters
        ! to the esm ocean, which is interpolated to the shelf base here.
        !
        ! With an active transient forcing, its spatially homogeneous anomalies
        ! (dTa = f_now*f_ta, dTo = f_now*f_to, dSo = f_now*f_so) go to the backend:
        ! snapclim uses them in place of its own index in its "anom" modes;
        ! snapesm adds them on top in every mode; esm has its own forcing.
        !
        ! The cold start (init) evaluates snapclim/snapesm at time_rel, the time
        ! loop at time; esm always at time.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
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
                    call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, domain=domain, &
                                         dTa=tsf%dTa, dTo=tsf%dTo, dSo=tsf%dSo, dx=dx, basins=basins)
                else
                    call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, domain=domain, &
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
                    call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, domain=domain, &
                                        dTa=tsf%dTa, dTo=tsf%dTo, dSo=tsf%dSo, dx=dx, basins=basins)
                else
                    call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, domain=domain, &
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
                                basins, domain, mshlf)
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

    logical function climate_ocean_const(cl) result(is_const)
        ! Is the backend's ocean held at its reference state (snapclim ocn_type
        ! = "const")? REMBO then adds its own ocean anomaly.
        type(yelmox_climate_class), intent(in) :: cl

        is_const = .false.
        if (trim(cl%method) == "snapclim") is_const = (trim(cl%snapclim%par%ocn_type) == "const")

    end function climate_ocean_const

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

    subroutine esm_update(cl, out, time, dtt, z_srf, H_ice, z_bed, f_grnd, z_sl, basins, domain, mshlf)
        ! The reference climatology at the current surface, the esm anomalies
        ! (historical / projection / homogeneous) and the variability, then the
        ! products: atmosphere, the surface mass balance (smb_method = climate),
        ! the ocean at the shelf base and subglacial discharge.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        real(wp),         intent(in) :: time, dtt
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in) :: domain
        type(marshelf_class), intent(in) :: mshlf

        associate(esm => cl%esm, ec => cl%esm_ctl)

        call esm_clim_update(esm, z_srf, time, ec%time_ref, ec%use_smb, domain, cl%grid_name)

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

        ! Surface mass balance at the present-day surface, with its anomaly and
        ! elevation gradient; surface temperature from the near-surface air.
        out%has_smb = ec%use_smb
        if (ec%use_smb) then
            out%now%smb     = esm%smb_ann
            out%now%dsmb    = sum(esm%dsmb, dim=3) / 12.0_wp
            out%now%dsmb_dz = esm%dsmbdz
        end if
        out%now%tsrf = sum(esm%t2m + esm%dts + esm%dts_var, dim=3) / 12.0_wp

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

end module yelmox_climate
