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
    ! Status: skeleton. Types + empty step_* + composition are in place; the
    ! per-step remap/update bodies and domain_init (grid registration + map
    ! priming) are filled in subsequent commits.

    use yelmo,        only : yelmo_class, wp
    use marine_shelf, only : marshelf_class
    use fastisostasy, only : isos_class, bsl_class
    use snapclim,     only : snapclim_class
    use smbpal,       only : smbpal_class
    use coupler,      only : coupler_class, coupler_init, coupler_add_grid, &
                             coupler_prime, remap

    implicit none
    private

    type domain_ctl
        ! Which components are active in this domain's coupling sequence.
        logical :: with_ice_sheet    = .true.
        logical :: with_isostasy     = .true.
        logical :: with_marine_shelf = .true.
        logical :: with_climate      = .true.

        ! Grid names (registered in dom%cpl); the source of truth for remap keys.
        character(len=256) :: grid_yelmo = ""   ! Yelmo working grid, e.g. "ANT-16KM"
        character(len=256) :: grid_mshlf = ""   ! marine-shelf grid, e.g. "ANT-2KM"
        character(len=256) :: grid_topo  = ""   ! hi-res topography grid (from file)
    end type domain_ctl

    type ice_domain
        type(yelmo_class)    :: yelmo
        type(marshelf_class) :: mshlf
        type(isos_class)     :: isos
        type(snapclim_class) :: snp
        type(smbpal_class)   :: smb
        type(bsl_class)      :: bsl
        type(coupler_class)  :: cpl     ! this region's grid registry + map cache
        type(domain_ctl)     :: ctl
    end type ice_domain

    public :: domain_ctl, ice_domain
    public :: yelmox_step
    public :: step_isostasy, step_icesheet, step_climate, step_marine_shelf

contains

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
