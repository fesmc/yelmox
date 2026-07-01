module yelmox_domain
    ! Multigrid yelmox: one region's full model state bundled as an ice_domain,
    ! advanced by composable step_* primitives.
    !
    ! Each helper model may live on its own grid; fields that cross a grid
    ! boundary are remapped through the domain's coupler (dom%cpl) at the moment
    ! of coupling. The step_* routines are the reusable pieces -- they are guarded
    ! internally by domain_ctl flags, and yelmox_step composes them. Bipolar runs
    ! are just an array of ice_domain looped through yelmox_step, so north and
    ! south share no state.
    !
    ! Status: bring-up. Types, domain_init (sub-model init on the Yelmo grid +
    ! htopo reference hub + coupler map priming), and the step_* composition are
    ! in place. The per-step remap/update bodies and the initial boundary state
    ! are lifted from yelmox.f90 in subsequent commits.

    use nml,          only : nml_read
    use yelmo,        only : yelmo_class, wp, yelmo_init
    use marine_shelf, only : marshelf_class, marshelf_init
    use fastisostasy, only : isos_class, isos_init, bsl_class, bsl_init
    use snapclim,     only : snapclim_class, snapclim_init
    use smbpal,       only : smbpal_class, smbpal_init
    use sediments,    only : sediments_class, sediments_init
    use geothermal,   only : geothermal_class, geothermal_init
    use htopo,        only : htopo_class, htopo_init
    use coupler,      only : coupler_class, coupler_init, coupler_prime, remap

    implicit none
    private

    type domain_ctl
        ! Run control (from the [ctrl] namelist group).
        character(len=56) :: tstep_method = "const"
        real(wp) :: tstep_const = 0.0_wp
        real(wp) :: time_init   = 0.0_wp
        real(wp) :: time_end    = 0.0_wp
        real(wp) :: time_equil  = 0.0_wp
        real(wp) :: dtt         = 10.0_wp
        real(wp) :: dt_restart  = 0.0_wp
        character(len=56) :: equil_method = "none"
        character(len=56) :: smb_method   = "smbpal"

        ! Which components are active in this domain's coupling sequence.
        logical :: with_ice_sheet    = .true.
        logical :: with_isostasy     = .true.
        logical :: with_marine_shelf = .true.
        logical :: with_climate      = .true.

        ! Domain + grid names (the source of truth for remap keys).
        character(len=256) :: domain     = ""   ! e.g. "Antarctica"
        character(len=256) :: grid_name  = ""   ! hi-res reference (htopo), highest res
        character(len=256) :: grid_yelmo = ""   ! Yelmo grid
        character(len=256) :: grid_mshlf = ""   ! marine-shelf grid (defaults to grid_name)
    end type domain_ctl

    type ice_domain
        type(yelmo_class)      :: yelmo
        type(marshelf_class)   :: mshlf
        type(isos_class)       :: isos
        type(snapclim_class)   :: snp
        type(smbpal_class)     :: smb
        type(bsl_class)        :: bsl
        type(sediments_class)  :: sed
        type(geothermal_class) :: gthrm
        type(htopo_class)      :: topo    ! hi-res geometry reference hub
        type(coupler_class)    :: cpl     ! this region's grid resolution + map cache
        type(domain_ctl)       :: ctl
    end type ice_domain

    public :: domain_ctl, ice_domain
    public :: domain_init, yelmox_step
    public :: step_isostasy, step_icesheet, step_climate, step_marine_shelf

contains

    subroutine domain_init(dom, path_par, time, time_rel)
        ! Initialize all sub-models of one domain (currently all on the Yelmo
        ! grid), load the hi-res reference hub, and prime the Yelmo<->topo maps.
        ! Lifted from yelmox.f90's init phase; the initial boundary state and the
        ! coupling time loop are added in subsequent commits.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: path_par
        real(wp),         intent(in)    :: time       ! model time
        real(wp),         intent(in)    :: time_rel   ! time before present [yr]

        character(len=256) :: domain

        ! --- run control ---
        call domain_ctl_load(dom%ctl, path_par)

        ! --- ice sheet (grid read from file) ---
        call yelmo_init(dom%yelmo, filename=path_par, grid_def="file", time=time)
        domain = trim(dom%yelmo%par%domain)
        dom%ctl%domain     = trim(domain)
        dom%ctl%grid_yelmo = trim(dom%yelmo%par%grid_name)

        ! --- external forcing models (all on the Yelmo grid for now) ---
        call bsl_init(dom%bsl, path_par, time_rel)

        call isos_init(dom%isos, path_par, "isos", dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                       dom%yelmo%grd%dx, dom%yelmo%grd%dy)

        call snapclim_init(dom%snp, path_par, domain, dom%yelmo%par%grid_name, &
                           dom%yelmo%grd%nx, dom%yelmo%grd%ny, dom%yelmo%bnd%basins)

        call smbpal_init(dom%smb, path_par, x=dom%yelmo%grd%xc, y=dom%yelmo%grd%yc, &
                         lats=dom%yelmo%grd%lat)

        call marshelf_init(dom%mshlf, path_par, "marine_shelf", dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                           domain, dom%yelmo%par%grid_name, dom%yelmo%bnd%regions, dom%yelmo%bnd%basins)

        call sediments_init(dom%sed, path_par, dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                            domain, dom%yelmo%par%grid_name)
        dom%yelmo%bnd%H_sed = dom%sed%now%H

        call geothermal_init(dom%gthrm, path_par, dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                             domain, dom%yelmo%par%grid_name)
        dom%yelmo%bnd%Q_geo = dom%gthrm%now%ghf

        ! --- hi-res reference hub + coupler ---
        call htopo_init(dom%topo, path_par, "htopo")
        dom%ctl%grid_name = trim(dom%topo%par%grid_name)
        if (len_trim(dom%ctl%grid_mshlf) == 0) dom%ctl%grid_mshlf = trim(dom%ctl%grid_name)

        ! Grids resolve from maps/grid_<name>.txt; prime the maps we will use.
        call coupler_init(dom%cpl)
        call coupler_prime(dom%cpl, dom%ctl%grid_yelmo, dom%ctl%grid_name, "bilin")  ! downscale
        call coupler_prime(dom%cpl, dom%ctl%grid_name, dom%ctl%grid_yelmo, "con")    ! aggregate

    end subroutine domain_init

    subroutine domain_ctl_load(ctl, path_par)
        type(domain_ctl), intent(inout) :: ctl
        character(len=*), intent(in)    :: path_par

        call nml_read(path_par, "ctrl", "tstep_method",   ctl%tstep_method)
        call nml_read(path_par, "ctrl", "tstep_const",    ctl%tstep_const)
        call nml_read(path_par, "ctrl", "time_init",      ctl%time_init)
        call nml_read(path_par, "ctrl", "time_end",       ctl%time_end)
        call nml_read(path_par, "ctrl", "time_equil",     ctl%time_equil)
        call nml_read(path_par, "ctrl", "dtt",            ctl%dtt)
        call nml_read(path_par, "ctrl", "dt_restart",     ctl%dt_restart)
        call nml_read(path_par, "ctrl", "with_ice_sheet", ctl%with_ice_sheet)
        call nml_read(path_par, "ctrl", "with_isostasy",  ctl%with_isostasy)
        call nml_read(path_par, "ctrl", "equil_method",   ctl%equil_method)
        ctl%smb_method = "smbpal"
        call nml_read(path_par, "ctrl", "smb_method",     ctl%smb_method)
    end subroutine domain_ctl_load

    subroutine yelmox_step(dom, time)
        ! Advance one domain by one coupling step. Component order is fixed here;
        ! each step_* is a no-op when its ctl flag is off.
        type(ice_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time

        call step_isostasy(dom, time)
        call step_icesheet(dom, time)
        call step_climate(dom, time)
        call step_marine_shelf(dom, time)
    end subroutine yelmox_step

    subroutine step_isostasy(dom, time)
        type(ice_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time
        if (.not. dom%ctl%with_isostasy) return
        ! TODO: bsl_update; isos_update on the isostasy grid; aggregate z_bed/z_ss
        !       (con) to the Yelmo grid and set as boundary conditions.
    end subroutine step_isostasy

    subroutine step_icesheet(dom, time)
        type(ice_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time
        if (.not. dom%ctl%with_ice_sheet) return
        ! TODO: yelmo_update on the Yelmo grid.
    end subroutine step_icesheet

    subroutine step_climate(dom, time)
        type(ice_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time
        if (.not. dom%ctl%with_climate) return
        ! TODO: snapclim_update; smbpal/smb_simple update; set Yelmo smb, T_srf.
    end subroutine step_climate

    subroutine step_marine_shelf(dom, time)
        type(ice_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time
        if (.not. dom%ctl%with_marine_shelf) return
        ! TODO: remap inputs (H_ice, z_bed, ...) Yelmo -> mshlf grid (bilin);
        !       marshelf_update_shelf + marshelf_update on the mshlf grid;
        !       aggregate bmb_shlf/T_shlf (con) back to Yelmo as forcing.
    end subroutine step_marine_shelf

end module yelmox_domain
