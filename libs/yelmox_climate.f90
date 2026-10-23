module yelmox_climate
    ! Climate backend of a domain, chosen at runtime ([comps] climate):
    !   "snapshots" -- snapshots blended by coefficients (snapshots)
    !   "esm"       -- reference climatology plus Earth-system-model anomalies over
    !                  historical/projection periods (esm_forcing)
    !   "rembo"     -- REMBO atmosphere and surface mass balance (climate_rembo,
    !                  built with rembo=1), with the ocean from snapshots (an
    !                  ocean-only reference snapshot)
    !
    ! Presents one interface (climate_init / climate_update) over the backends,
    ! filling a `climate_out_class`; the domain and its coupling read only that.
    ! The forcing of a domain is its terms (climate_terms, [forcing<sfx>]):
    ! the snapshot terms are read by the snapshot backend, the uniform and
    ! replace terms applied to `climate_out` after it, around the derived
    ! coefficients the climate sets (the mean sea-level air anomaly); then the
    ! transforms (climate_mods: taper, bias). Every term and transform is
    ! checked at init against what the backend supplies. rembo takes the
    ! uniform tas terms as its input (rembo_update).

    use precision,     only : wp, dp
    use nml,           only : nml_read
    use ncio
    use coordinates,   only : grid_class
    use timestepping,  only : tstep_class
    use climate_out,   only : climate_out_class, clim_state_class, summer_mean
    use snapshots,     only : snapshots_class, climate_state_class, snapshots_init, snapshots_update, &
                              snapshots_end, snapshots_write_init, snapshots_write_step, &
                              snapshots_lapse_monthly
    use esm_forcing,   only : esm_forcing_class, esm_forcing_init, esm_clim_update, &
                              esm_forcing_update, esm_variability_update, &
                              esm_lapse_monthly
    use marine_shelf,  only : marshelf_class, marshelf_interp_shelf, ocn_variable_extrapolation
    use kryos_forcing, only : coefs_class, coefs_set, coefs_is_derived, &
                              coefs_derived, derived_spec
    use climate_terms, only : climate_caps_class, terms_class, terms_init, terms_check, terms_eval, &
                              terms_save_ocean, terms_apply_atm, terms_apply_ocn, &
                              terms_uniform_dT, terms_write_step
    use climate_mods,  only : climate_mods_class, climate_mods_init, &
                              climate_mods_apply, climate_mods_write_init, climate_mods_write_step
    use climate_rembo, only : rembo_clim_init, rembo_clim_update, rembo_clim_dT_ann, &
                              rembo_clim_restart_write
    use regions,       only : regions_class, regions_basin_ids, regions_select
    use yelmox_nml,    only : VAR_DEFS_MAX, forcing_var_defs, def_file => yelmox_defaults_file

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
    ! air temperature and the anomaly it takes from the uniform tas terms.
    type rembo_state_class
        real(wp), allocatable :: ta_sum(:,:)     ! [K] near-surface air temperature (summer)
        real(wp) :: dT_sum = 0.0_wp              ! [K] summer air temperature anomaly
        real(wp) :: dT_ann = 0.0_wp              ! [K] annual air temperature anomaly (REMBO's, of dT_sum)
    end type rembo_state_class

    ! A derived coefficient this domain computes: the mean sea-level air anomaly
    ! of its climate update over a region (factor * mean(tsl_ann - ref)).
    type derived_air_class
        character(len=64)    :: name
        real(wp)             :: factor
        logical, allocatable :: mask(:,:)
        integer              :: n
    end type derived_air_class

    type yelmox_climate_class
        character(len=16)       :: method = ""   ! snapshots | esm | rembo
        character(len=256)      :: grid_name     ! the climate grid (grid_clim)
        logical                 :: south = .false.   ! southern hemisphere
        type(snapshots_class)   :: snapshots
        type(esm_forcing_class) :: esm
        type(esm_ctl_class)     :: esm_ctl
        type(rembo_state_class) :: rembo
        character(len=64)       :: reference     ! the reference snapshot ([forcing<sfx>] reference)
        type(terms_class)       :: terms         ! the terms ([forcing<sfx>] terms)
        type(derived_air_class), allocatable :: dair(:)   ! the air_anomaly coefficients it sets
        type(climate_mods_class):: mods          ! the transforms after the terms
        real(wp), allocatable   :: basins(:,:)   ! basin ids on the climate grid (esm: [esm] basins)
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update
    public :: climate_apply_ocean
    public :: climate_end
    public :: climate_file_base, climate_write_init, climate_write_2D, climate_write_1D
    public :: climate_restart_write

    ! Copy of the fields a state supplies (allocated) into climate_out.
    interface copy_allocated
        module procedure copy_allocated_1D, copy_allocated_2D, copy_allocated_3D
    end interface

contains

    subroutine climate_init(cl, method, filename, domain, grid, time, reg, &
                            south, sfx, timeline_group, smb_direct, coefs)
        ! south: the domain lies in the southern hemisphere (seasons, lapse rates).
        ! coefs: the driver's coefficient registry (the terms check the
        ! coefficients they name against it).
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
        type(coefs_class),   intent(in) :: coefs

        integer :: nx, ny, k, n
        type(climate_caps_class) :: caps
        type(derived_spec)       :: drv
        character(len=256)       :: var_defs(VAR_DEFS_MAX)

        nx        = grid%G%nx
        ny        = grid%G%ny

        cl%method    = trim(method)
        cl%grid_name = trim(grid%name)
        cl%south     = south

        ! The forcing of the domain: its reference snapshot and its terms.
        call nml_read(filename, "forcing"//trim(sfx), "reference", cl%reference, init=.true., &
                      defaults_file=def_file, defaults_group="forcing")
        var_defs = forcing_var_defs(filename)
        call terms_init(cl%terms, filename, sfx, var_defs, coefs)

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

        ! The backend, and what it supplies (caps) for the contracts of the
        ! terms and transforms. It supplies the surface mass balance exactly
        ! when it is taken from it (smb_direct: esm with use_smb, rembo;
        ! checked above); otherwise its atmosphere drives the surface.
        caps%name = trim(cl%method)
        select case(trim(cl%method))
            case("snapshots")
                call snapshots_init(cl%snapshots, filename, domain, grid, time, &
                                  south, cl%terms, "snapshots"//trim(sfx), cl%reference)
                ! The surface reads the atmosphere; an ocean-only reference
                ! serves the rembo ocean only.
                if (.not. cl%snapshots%par%with_atm) then
                    write(*,*) "climate_init:: error: climate = snapshots needs an atmosphere: &
                               &the reference snapshot ([forcing"//trim(sfx)//"] reference) binds &
                               &no tas/pr."
                    error stop 1
                end if
                caps%atm_drives_srf = .true.
                caps%has_tsl_ann    = .true.
                caps%has_snapshots  = .true.
                caps%has_ocn_prof   = cl%snapshots%par%with_ocn
            case("esm")
                call esm_init(cl, filename, domain, grid, "esm"//trim(sfx), timeline_group, &
                              smb_direct, reg)
                caps%atm_drives_srf = .not. smb_direct
                caps%has_ocn_shelf  = .true.
            case("rembo")
                ! REMBO for the atmosphere and smb, snapshots for the ocean: an
                ! ocean-only reference snapshot (to, optionally so).
                call snapshots_init(cl%snapshots, filename, domain, grid, time, &
                                  south, cl%terms, "snapshots"//trim(sfx), cl%reference)
                if (cl%snapshots%par%with_atm .or. .not. cl%snapshots%par%with_ocn) then
                    write(*,*) "climate_init:: error: climate = rembo takes the ocean only from &
                               &the snapshots: the reference snapshot ([forcing"//trim(sfx)// &
                               "] reference) must bind to (and optionally so), not tas/pr &
                               &(REMBO supplies the atmosphere)."
                    error stop 1
                end if
                call rembo_clim_init(time, nx, ny)
                allocate(cl%rembo%ta_sum(nx,ny))
                cl%rembo%ta_sum = 0.0_wp
                caps%anom_input    = .true.
                caps%has_snapshots = .true.
                caps%has_ocn_prof  = .true.
            case default
                write(*,*) "climate_init:: error: climate must be snapshots, esm or rembo; got ", &
                           trim(cl%method)
                error stop 1
        end select

        ! The terms' contracts against the backend, and the derived air
        ! anomaly coefficients this domain computes (one per term using one), their
        ! region resolved on the climate grid.
        call terms_check(cl%terms, caps)
        n = 0
        do k = 1, cl%terms%n
            if (is_air_anomaly(cl%terms%term(k)%coef)) n = n + 1
        end do
        if (allocated(cl%dair)) deallocate(cl%dair)
        allocate(cl%dair(n))
        n = 0
        do k = 1, cl%terms%n
            associate(t => cl%terms%term(k))
            if (.not. is_air_anomaly(t%coef)) cycle
            if (any(cl%dair(1:n)%name == t%coef)) cycle
            if (.not. caps%has_tsl_ann) then
                write(*,*) "climate_init:: error: ["//trim(t%group)//"] coef "//trim(t%coef)// &
                           " is derived from the mean sea-level air anomaly, which climate = "// &
                           trim(caps%name)//" does not supply (no tsl_ann)."
                error stop 1
            end if
            n = n + 1
            drv = coefs_derived(coefs, t%coef)
            cl%dair(n)%name   = t%coef
            cl%dair(n)%factor = drv%factor
            cl%dair(n)%mask   = regions_select(reg, trim(drv%region))
            cl%dair(n)%n      = count(cl%dair(n)%mask)
            if (cl%dair(n)%n == 0) then
                write(*,*) "climate_init:: error: [coef_"//trim(t%coef)//"] region = '"// &
                           trim(drv%region)//"' selects no cell of the climate grid."
                error stop 1
            end if
            write(*,*) "climate_init:: coefficient "//trim(t%coef)//" = ", drv%factor, &
                       " x the mean sea-level air anomaly over '"//trim(drv%region)//"' (", &
                       cl%dair(n)%n, " of ", size(cl%dair(n)%mask), " cells), set by this climate"
            end associate
        end do
        if (allocated(cl%dair)) cl%dair = cl%dair(1:n)

        ! The transforms, checked against the backend.
        call climate_mods_init(cl%mods, filename, sfx, caps, domain, grid, time)

    contains

        logical function is_air_anomaly(name)
            character(len=*), intent(in) :: name
            is_air_anomaly = .false.
            if (len_trim(name) == 0) return
            if (.not. coefs_is_derived(coefs, name)) return
            drv = coefs_derived(coefs, name)
            is_air_anomaly = trim(drv%from) == "air_anomaly"
        end function is_air_anomaly

    end subroutine climate_init

    subroutine climate_update(cl, out, ts, z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref, &
                              dtt, mshlf, coefs, init)
        ! Evaluate the terms (their coefficients), update the backend on the
        ! climate grid and fill `out`, then its elevation dependence, the
        ! uniform terms on the atmosphere, the derived air-anomaly
        ! coefficients, the ocean terms and the transforms. The
        ! geometry is the domain's, on the climate grid, with its present-day
        ! surface (z_srf_ref); the marine shelf lends its parameters to the esm
        ! ocean, which is interpolated to the shelf base here. coefs is the
        ! driver's coefficient registry, advanced to the forcing time; the
        ! terms read it and the derived coefficients of this domain are set in it.
        !
        ! snapshots and rembo are evaluated at the forcing time time_rel
        ! (with tstep_method = "const", time is the elapsed time and time_rel the
        ! fixed forcing time); esm at time. rembo equilibrates at the cold start.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:)
        real(wp),         intent(in) :: z_srf_ref(:,:)
        real(dp),         intent(in) :: dtt              ! [yr] time step of the run
        type(marshelf_class),  intent(in)    :: mshlf
        type(coefs_class),     intent(inout) :: coefs
        logical,               intent(in), optional :: init

        logical  :: is_init
        real(dp) :: time
        integer  :: k

        is_init = .false.
        if (present(init)) is_init = init

        time = ts%time_rel

        ! The value of every term's coefficient for this update.
        call terms_eval(cl%terms, coefs)

        select case(trim(cl%method))
            case("snapshots")
                call snapshots_update(cl%snapshots, z_srf=z_srf, time=time, terms=cl%terms)
                call snapshots_state_copy(out%now, cl%snapshots%now)
                call snapshots_state_copy(out%ref, cl%snapshots%ref)

            case("esm")
                call esm_update(cl, out, ts%time, dtt, z_srf, H_ice, z_bed, f_grnd, z_sl, &
                                z_srf_ref, mshlf)

            case("rembo")
                call rembo_update(cl, out, ts, time, z_srf, H_ice, z_sl, is_init)
        end select

        ! The elevation dependence of the atmosphere, to move it to another surface
        ! (kryos_coupling): its surface z_srf, the monthly lapse rate and the
        ! precipitation sensitivity, as the backend applies them (the uniform
        ! terms read beta_p). rembo: none.
        out%now%z_srf = z_srf
        if (.not. allocated(out%now%beta_p)) allocate(out%now%beta_p(size(z_srf,1),size(z_srf,2)))
        select case(trim(cl%method))
            case("snapshots")
                out%lapse       = snapshots_lapse_monthly(cl%snapshots)
                out%now%beta_p  = cl%snapshots%now%beta_p
            case("esm")
                out%lapse       = esm_lapse_monthly(cl%esm, cl%south)
                out%now%beta_p  = cl%esm%beta_p
            case default
                out%lapse       = 0.0_wp
                out%now%beta_p  = 0.0_wp
        end select

        ! The terms: the backend's ocean kept as the base of the ocean terms;
        ! the uniform terms on the atmosphere (rembo took them as its input);
        ! the derived air anomaly of this update, and the terms it drives
        ! re-evaluated; the ocean terms; then the transforms.
        call terms_save_ocean(cl%terms, out)
        if (.not. cl%mods%caps%anom_input) call terms_apply_atm(cl%terms, out)
        do k = 1, size(cl%dair)
            call coefs_set(coefs, cl%dair(k)%name, cl%dair(k)%factor * &
                 sum(out%now%tsl_ann - out%ref%tsl_ann, mask=cl%dair(k)%mask) / real(cl%dair(k)%n, wp))
        end do
        call terms_eval(cl%terms, coefs, derived=.true.)
        call terms_apply_ocn(cl%terms, out)
        call climate_mods_apply(cl%mods, out, time)

    end subroutine climate_update

    subroutine climate_apply_ocean(cl, out, coefs)
        ! Re-evaluate the terms with a derived coefficient and re-apply the
        ! ocean terms on the backend's ocean of the last update (idempotent):
        ! for a driver that sets a derived coefficient after the climate update
        ! (the bipolar bridge, the OBM box temperatures). The transforms do not
        ! touch the ocean, so nothing else is redone.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(coefs_class),          intent(in)    :: coefs
        call terms_eval(cl%terms, coefs, derived=.true.)
        call terms_apply_ocn(cl%terms, out)
    end subroutine climate_apply_ocean

    subroutine snapshots_state_copy(dst, src)
        ! The fields of a snapshots state that climate_out carries, those it built:
        ! the atmosphere and/or the ocean of the reference snapshot.
        type(clim_state_class),    intent(inout) :: dst
        type(climate_state_class), intent(in)    :: src

        call copy_allocated(dst%z_srf,   src%z_srf)
        call copy_allocated(dst%tas,     src%tas)
        call copy_allocated(dst%pr,      src%pr)
        call copy_allocated(dst%tsl_ann, src%tsl_ann)
        call copy_allocated(dst%ta_ann,  src%ta_ann)
        call copy_allocated(dst%pr_ann,  src%pr_ann)
        call copy_allocated(dst%to_ann,  src%to_ann)
        call copy_allocated(dst%so_ann,  src%so_ann)
        call copy_allocated(dst%depth,   src%depth)

    end subroutine snapshots_state_copy

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
        ! Teardown of the backend (snapshots, also the ocean of rembo); esm has
        ! none.
        type(yelmox_climate_class), intent(inout) :: cl

        select case(trim(cl%method))
            case("snapshots", "rembo")
                call snapshots_end(cl%snapshots)
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
            case("snapshots")
                call snapshots_write_init(cl%snapshots, filename)
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

        call write2D("t2m_ann", out%now%ta_ann, "K", "Near-surface air temperature (ann)", &
                     prec_abs=1e-3_wp)
        if (allocated(out%now%tas)) &
            call write2D("t2m_sum", summer_mean(out%now%tas, cl%south), "K", &
                         "Near-surface air temperature (sum)", prec_abs=1e-3_wp)
        if (allocated(out%now%pr_ann)) &
            call write2D("pr_ann", out%now%pr_ann, "mm/a", "Precipitation (ann)", &
                         prec_rel=1e-4_wp)
        if (out%has_smb) then
            ! The climate's own surface fields (with surface_method = "climate",
            ! those landed on Yelmo).
            call write2D("smb_ann",  1e-3_wp*out%now%smb, "m/a water equiv.", "Surface mass balance (ann)", &
                         prec_abs=1e-4_wp)
            call write2D("tsrf_ann", out%now%tsrf, "K", "Surface temperature (ann)", &
                         prec_abs=1e-3_wp)
        end if
        if (allocated(out%now%to_ann)) &
            call write3D("to_ann", out%now%to_ann, "depth", "K", "Ocean temperature (ann)", &
                         prec_abs=1e-3_wp)
        if (allocated(out%now%so_ann)) &
            call write3D("so_ann", out%now%so_ann, "depth", "PSU", "Ocean salinity (ann)", &
                         prec_rel=1e-4_wp)

        select case(trim(cl%method))
            case("snapshots")
                ! The monthly atmosphere, then the backend's own state (the axes).
                call write3D("tas",     out%now%tas,     "month", "K",    "Near-surface air temperature", &
                             prec_abs=1e-3_wp)
                call write3D("pr",      out%now%pr,      "month", "mm/d", "Precipitation", &
                             prec_rel=1e-4_wp)
                call write2D("tsl_ann", out%now%tsl_ann, "K", "Sea-level air temperature (ann)", &
                             prec_abs=1e-3_wp)
                call snapshots_write_step(cl%snapshots, filename, ncid, n)
            case("esm")
                call write_esm()
            case("rembo")
                ! Annual fields only: REMBO's own summer temperature.
                call write2D("t2m_sum", cl%rembo%ta_sum, "K", "Near-surface air temperature (sum)", &
                             prec_abs=1e-3_wp)
        end select

        call climate_mods_write_step(cl%mods, filename, ncid, n)
        call terms_write_step(cl%terms, filename, ncid, n)

    contains

        subroutine write_esm()

            associate(esm => cl%esm)
            call write2D("dts_ann", sum(esm%dts, dim=3)/12.0_wp, "K", "Surface air temperature anomaly", &
                         prec_abs=1e-3_wp)
            if (cl%esm_ctl%use_smb) then
                call write2D("smb_ref_ann", 1e-3_wp*esm%smb_ann, "m/a water equiv.", "Reference SMB (ann)", &
                             prec_abs=1e-4_wp)
                call write2D("dsmb_ann", 1e-3_wp*sum(esm%dsmb, dim=3)/12.0_wp, "m/a water equiv.", &
                             "SMB anomaly (ann)", prec_abs=1e-4_wp)
                call write2D("dsmbdz",   1e-3_wp*esm%dsmbdz, "m/a m-1 water equiv.", "SMB lapse rate", &
                             prec_rel=1e-4_wp)
            else
                call write2D("dpr_ann", sum(esm%dpr, dim=3)/12.0_wp, "%", "Precipitation anomaly (ann)", &
                             prec_rel=1e-4_wp)
            end if
            call write2D("dto",    esm%dto,    "K",    "Shelf temperature anomaly", prec_abs=1e-3_wp)
            call write2D("dso",    esm%dso,    "PSU",  "Shelf salinity anomaly",    prec_rel=1e-4_wp)
            call write2D("Qd_ann", esm%Qd_ann, "m3/s", "Subglacial discharge (ann)", prec_rel=1e-4_wp)
            call write2D("Qd_sum", esm%Qd_sum, "m3/s", "Subglacial discharge (sum)", prec_rel=1e-4_wp)
            end associate

        end subroutine write_esm

        ! prec_abs/prec_rel: precision of the written values (ncio bit
        ! rounding), absent for lossless.
        subroutine write2D(vnm, var, units, long_name, prec_abs, prec_rel)
            character(len=*), intent(in) :: vnm, units, long_name
            real(wp),         intent(in) :: var(:,:)
            real(wp),         intent(in), optional :: prec_abs, prec_rel
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                          units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
        end subroutine write2D

        subroutine write3D(vnm, var, dim3, units, long_name, prec_abs, prec_rel)
            character(len=*), intent(in) :: vnm, dim3, units, long_name
            real(wp),         intent(in) :: var(:,:,:)
            real(wp),         intent(in), optional :: prec_abs, prec_rel
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3=dim3, dim4="time", &
                          start=[1,1,1,n], count=[size(var,1),size(var,2),size(var,3),1], &
                          ncid=ncid, units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
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
            ! (ocean, at the shelf base). The air temperature and precipitation
            ! as applied (out%now: variability, terms and transforms included;
            ! the precipitation only when the climate has it, not with the esm
            ! smb); the anomalies those of the esm.

            ! Missing value of the ocean means: a field with no valid floating cell.
            real(wp), parameter :: mv_ocn = -9999.0_wp
            logical,  allocatable :: mask_tot(:,:), mask_flt(:,:)
            logical  :: with_pr
            integer  :: npts_tot
            real(wp) :: t2m_1d, pr_1d, dt_1d, dpr_1d, dt_var_1d, dpr_var_1d
            real(wp) :: to_1d, so_1d, dto_1d, dso_1d, dto_var_1d, dso_var_1d

            associate(esm => cl%esm)

            mask_tot = (H_ice > 0.0_wp)
            mask_flt = (H_ice > 0.0_wp .and. f_grnd == 0.0_wp)
            npts_tot = count(mask_tot)
            with_pr  = allocated(out%now%pr_ann)

            if (npts_tot > 0) then
                t2m_1d     = sum(out%now%ta_ann, mask=mask_tot) / npts_tot
                if (with_pr) pr_1d = sum(out%now%pr_ann/365.0_wp, mask=mask_tot) / npts_tot
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
            if (with_pr) &
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
        ! geometry on the climate grid: REMBO's (rembo_restart.nc); snapshots
        ! and esm hold no state (the coefficients of snapshots restart with the
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
                                f_grnd, z_sl, ec%use_smb)

        call esm_variability_update(esm, mshlf, time, dtt, ec%clim_var, ec%time_ref, H_ice, &
                                    basins, z_bed, f_grnd, z_sl, ec%use_var)

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

    subroutine rembo_update(cl, out, ts, time, z_srf, H_ice, z_sl, init)
        ! REMBO's atmosphere and surface mass balance, and the ocean from snapshots
        ! (at `time`, as for the snapshots backend). REMBO needs the anomaly
        ! before it runs, so the uniform terms on tas are its input (not applied
        ! after it): their sum f_ta*f is the summer air temperature anomaly
        ! dT_sum; dT_ann is the annual mean anomaly REMBO applies with it (its
        ! seasonal cycle, rembo_clim_dT_ann), a diagnostic.
        ! Its ocean comes from the snapshots and the ocean terms (climate_update).
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        type(tstep_class),          intent(in)    :: ts
        real(dp),         intent(in) :: time
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)
        logical,          intent(in) :: init

        integer :: nx, ny

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        ! The summer anomaly REMBO takes as its input: the uniform terms on tas.
        cl%rembo%dT_sum = terms_uniform_dT(cl%terms)
        cl%rembo%dT_ann = rembo_clim_dT_ann(cl%rembo%dT_sum)

        ! Atmosphere and surface mass balance, at the current surface.
        if (.not. allocated(out%now%ta_ann)) then
            allocate(out%now%ta_ann(nx,ny), out%now%pr_ann(nx,ny))
            allocate(out%now%smb(nx,ny),    out%now%tsrf(nx,ny))
        end if
        call rembo_clim_update(ts%time, ts%time_rel, cl%rembo%dT_sum, z_srf, H_ice, z_sl, init, &
                               out%now%ta_ann, cl%rembo%ta_sum, out%now%pr_ann, &
                               out%now%smb, out%now%tsrf)
        out%has_smb = .true.

        ! Ocean (depth profiles) from snapshots (ocean only); its terms are
        ! applied by climate_update.
        call snapshots_update(cl%snapshots, z_srf=z_srf, time=time, terms=cl%terms)
        call snapshots_state_copy(out%now, cl%snapshots%now)
        call snapshots_state_copy(out%ref, cl%snapshots%ref)

    end subroutine rembo_update

end module yelmox_climate
