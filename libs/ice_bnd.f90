module ice_bnd
    ! Prescribed ice thickness in the fixed and relaxation cells of Yelmo
    ! ([ice_bnd<sfx>], defaults: [ice_bnd] of input/yelmox_defaults.nml): the
    ! thickness of an ice history (a varslice group, e.g. an IceHistory
    ! product of FesmData) written into Yelmo's reference ice thickness
    ! H_ice_ref, where Yelmo prescribes the ice (bnd%mask_ice =
    ! MASK_ICE_FIXED: H_ice = H_ice_ref) or relaxes it (bnd%tau_relax > 0,
    ! ytopo.topo_rel = -1: towards H_ice_ref). The interior keeps its
    ! H_ice_ref. The cells are those of yelmo_masks (border_width, mask_ice_fixed,
    ! relax, yelmo.mask_border), e.g. the strip at the border of GreenlandPaleo
    ! (docs/_grl-paleo-plan.md, "Prescribed ice at the domain borders").
    !
    ! method:
    !   file    the thickness at the forcing time, interpolated in time
    !           between the slices of the group, the first or last slice
    !           outside its time range (varslice "extrap")
    !   scaled  H = f * H_slice: the slice of the group at slice_time times the
    !           value f of a coefficient of the registry (e.g. before the file
    !           starts: the LGM slice times a glacial coefficient)
    !
    ! The group is read on the Yelmo grid ({grid_name} = the ice grid; remap
    ! and grid_src in the group for a product on another grid).

    use precision,     only : wp, dp
    use nml,           only : nml_read, nml_validate
    use coordinates,   only : grid_class
    use constants,     only : mv
    use varslice,      only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use yelmo_defs,    only : MASK_ICE_FIXED
    use kryos_forcing, only : coefs_class, coefs_value, coefs_has
    use yelmox_nml,    only : VAR_DEFS_MAX, var_defs_find, forcing_var_defs, def_file => yelmox_defaults_file

    implicit none

    private

    type ice_bnd_class
        logical            :: active     = .false.
        character(len=64)  :: group      = ""
        character(len=16)  :: method     = "file"     ! file | scaled
        character(len=256) :: H_ice      = ""         ! var_defs group of the ice thickness [m]
        character(len=64)  :: coef       = ""         ! scaled: the coefficient f
        real(dp)           :: slice_time = 0.0_dp     ! [yr] scaled: time of H_slice
        type(varslice_class)  :: vs                   ! file: the ice history
        real(wp), allocatable :: H_slice(:,:)         ! scaled: the slice at slice_time [m]
    end type ice_bnd_class

    public :: ice_bnd_class
    public :: ice_bnd_init
    public :: ice_bnd_update

contains

    subroutine ice_bnd_init(ibnd, filename, group, domain, grid, coefs)
        ! Read [ice_bnd<sfx>] (group; optional, inactive without it) and open
        ! the ice history on the Yelmo grid (grid). coefs: the driver's
        ! coefficient registry (scaled: the coefficient must be in it).
        type(ice_bnd_class), intent(inout) :: ibnd
        character(len=*),    intent(in)    :: filename, group, domain
        type(grid_class),    intent(in)    :: grid
        type(coefs_class),   intent(in)    :: coefs

        character(len=256) :: var_defs(VAR_DEFS_MAX)
        character(len=1024) :: vdb

        ibnd%group = group
        call nml_validate(filename, def_file, group, defaults_group="ice_bnd")
        call nml_read(filename, group, "active", ibnd%active, defaults_file=def_file, defaults_group="ice_bnd")
        if (.not. ibnd%active) return

        var_defs = forcing_var_defs(filename)
        call nml_read(filename, group, "method",     ibnd%method,     defaults_file=def_file, defaults_group="ice_bnd")
        call nml_read(filename, group, "H_ice",      ibnd%H_ice,      defaults_file=def_file, defaults_group="ice_bnd")
        call nml_read(filename, group, "coef",       ibnd%coef,       defaults_file=def_file, defaults_group="ice_bnd")
        call nml_read(filename, group, "slice_time", ibnd%slice_time, defaults_file=def_file, defaults_group="ice_bnd")

        select case(trim(ibnd%method))
            case("file", "scaled")
                ! ok
            case default
                write(*,*) "ice_bnd_init:: error: ["//trim(group)//"] method must be file or &
                           &scaled; got '"//trim(ibnd%method)//"'."
                error stop 1
        end select
        if (trim(ibnd%method) == "scaled" .and. .not. coefs_has(coefs, ibnd%coef)) then
            write(*,*) "ice_bnd_init:: error: ["//trim(group)//"] method = scaled: coef = '"// &
                       trim(ibnd%coef)//"' is not in [forcing] coefs."
            error stop 1
        end if

        vdb = var_defs_find(var_defs, filename, trim(ibnd%H_ice), "["//trim(group)//"] H_ice")
        call varslice_init_nml(ibnd%vs, trim(vdb), trim(ibnd%H_ice), domain=domain, &
                               grid_name=trim(grid%name), grid=grid)

        if (trim(ibnd%method) == "scaled") then
            call read_slice(ibnd, ibnd%slice_time)
            ibnd%H_slice = ibnd%vs%var(:,:,1,1)
            call varslice_end(ibnd%vs)
            write(*,*) "ice_bnd_init:: "//trim(group)//" active, method scaled: coefficient "// &
                       trim(ibnd%coef)//" times "//trim(ibnd%H_ice)//" at ", ibnd%slice_time, &
                       " yr; range [m]: ", minval(ibnd%H_slice), maxval(ibnd%H_slice)
        else
            write(*,*) "ice_bnd_init:: "//trim(group)//" active, method file: "//trim(ibnd%H_ice)
        end if

    end subroutine ice_bnd_init

    subroutine ice_bnd_update(ibnd, H_ice_ref, mask_ice, tau_relax, time, coefs)
        ! H_ice_ref = the prescribed thickness at the forcing time `time` in the
        ! fixed (mask_ice = MASK_ICE_FIXED) and relaxation (tau_relax > 0)
        ! cells; elsewhere unchanged. coefs: the registry, advanced to `time`.
        type(ice_bnd_class), intent(inout) :: ibnd
        real(wp),            intent(inout) :: H_ice_ref(:,:)
        integer,             intent(in)    :: mask_ice(:,:)
        real(wp),            intent(in)    :: tau_relax(:,:)
        real(dp),            intent(in)    :: time
        type(coefs_class),   intent(in)    :: coefs

        real(wp), allocatable :: H(:,:)
        logical,  allocatable :: strip(:,:)

        if (.not. ibnd%active) return

        select case(trim(ibnd%method))
            case("file")
                call read_slice(ibnd, time)
                H = ibnd%vs%var(:,:,1,1)
            case("scaled")
                H = coefs_value(coefs, ibnd%coef) * ibnd%H_slice
        end select

        ! (scaled: no negative thickness for a coefficient below 0)
        strip = (mask_ice == MASK_ICE_FIXED .or. tau_relax > 0.0_wp)
        where (strip) H_ice_ref = max(H, 0.0_wp)

    end subroutine ice_bnd_update

    subroutine read_slice(ibnd, time)
        ! The thickness of the group at `time` (clamped to its time range) into
        ! ibnd%vs%var; stops on missing values.
        type(ice_bnd_class), intent(inout) :: ibnd
        real(dp),            intent(in)    :: time

        if (ibnd%vs%par%with_time) then
            call varslice_update(ibnd%vs, [time], method="extrap")
        else
            call varslice_update(ibnd%vs)
        end if

        if (any(ibnd%vs%var(:,:,1,1) == mv)) then
            write(*,*) "ice_bnd:: error: ["//trim(ibnd%group)//"] "//trim(ibnd%H_ice)// &
                       " has missing values at time ", time, " (set fill in the group)."
            error stop 1
        end if

    end subroutine read_slice

end module ice_bnd
