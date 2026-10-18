module climate_mods
    ! Climate modifiers: changes applied to a backend's `climate_out` after the
    ! backend update, in a fixed order (docs/_forcing-design.md, section 7):
    !   anom         -- homogeneous anomaly from one index of the driver's
    !                   registry: dTa = f*f_ta, dTo = f*f_to, dSo = f*f_so
    !   ocn_fraction -- ocean temperature anomaly = f_to * the mean sea-level
    !                   air temperature anomaly over a region
    !   bias         -- static bias correction of the atmosphere
    ! (later ocn_obm and taper, between ocn_fraction and bias). Each is an
    ! optional namelist group with `active` (default off), suffixed by the
    ! domain (anom<sfx>, ocn_fraction<sfx>, bias<sfx>); unset keys (an absent
    ! group: all) from the templates anom, ocn_fraction, bias of
    ! input/yelmox_defaults.nml.
    !
    ! Contract: each modifier states what it touches in `now` (`ref` is never
    ! touched) and which backends it refuses; climate_mods_init checks every
    ! active modifier against what the backend supplies (climate_caps_class)
    ! and stops on a mismatch, never a silent partial application.
    !   anom: tas, ta_ann, tsl_ann + dTa; pr, pr_ann * exp(beta_p*dTa); the ocean
    !         profiles to_ann/so_ann + dTo/dSo, or at the shelf base T_shlf/S_shlf
    !         and dT_shlf/dS_shlf. Refused when the backend supplies the smb
    !         (has_smb), which a post-hoc air anomaly cannot reach. Exception:
    !         rembo takes the anomaly as an input (dT_sum) and reads [anom]
    !         itself (yelmox_climate, rembo_update); anom is skipped here for it.
    !         With f_to or f_so /= 0, a climate without any ocean field stops at
    !         the first update (the fields are known only then).
    !   ocn_fraction: to_ann = ref%to_ann + f_to * mean over the region of
    !         (tsl_ann - ref%tsl_ann), at every depth; salinity untouched. Uses
    !         the current update's tsl_ann, which carries anom (no lag). Needs
    !         ocean profiles and tsl_ann: refused for an ocean at the shelf base
    !         (esm) and a backend without tsl_ann (esm, rembo); refused with an
    !         anom f_to /= 0, since it replaces the ocean temperature anomaly.
    !   bias: tas, ta_ann, tsl_ann + dtas; pr * fpr, pr_ann by the ratio of the
    !         annual totals. Refused for a backend with annual fields only (rembo).

    use precision,     only : wp
    use nml,           only : nml_read, nml_validate
    use ncio,          only : nc_write
    use coordinates,   only : grid_class
    use varslice,      only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use regions,       only : regions_class, regions_select
    use climate_out,   only : climate_out_class
    use kryos_forcing, only : indices_class, indices_value, indices_has
    use yelmox_nml,    only : VAR_DEFS_MAX, var_defs_find, def_file => yelmox_defaults_file

    implicit none

    private

    ! What a backend supplies, for the contracts (set by climate_init).
    type climate_caps_class
        character(len=16) :: name = ""           ! the backend, for messages
        logical :: has_smb       = .false.       ! supplies the surface mass balance
        logical :: has_ocn_shelf = .false.       ! ocean at the shelf base (else profiles)
        logical :: has_tsl_ann   = .false.       ! sea-level air temperature (now and ref)
        logical :: monthly       = .true.        ! monthly air temperature and precipitation
        logical :: anom_input    = .false.       ! reads [anom] itself, as an input (rembo)
    end type climate_caps_class

    ! [anom<sfx>]: one index mapped onto homogeneous anomalies by the gains.
    type anom_class
        logical           :: active = .false.
        character(len=64) :: group  = ""
        character(len=64) :: index  = ""        ! the index driving the anomalies
        real(wp) :: f_ta = 1.0_wp               ! [K]   dTa = f * f_ta
        real(wp) :: f_to = 0.0_wp               ! [K]   dTo = f * f_to
        real(wp) :: f_so = 0.0_wp               ! [psu] dSo = f * f_so
    end type anom_class

    ! [ocn_fraction<sfx>]: the ocean temperature anomaly as a fraction of the
    ! mean sea-level air temperature anomaly over a region. region: a selection
    ! expression of the domain's regions on the climate grid ([regions<sfx>]),
    ! resolved once at init: "all" (every cell), or layer:values terms, e.g.
    ! "zone:open_ocean,continental_shelf", "region:Greenland", "~zone:land",
    ! "Zwally2012:11,12", joined by & and | (fesm-utils regions_select).
    type ocn_fraction_class
        logical            :: active = .false.
        character(len=64)  :: group  = ""
        real(wp)           :: f_to   = 0.0_wp   ! [1] dTo = f_to * mean(tsl_ann - ref%tsl_ann)
        character(len=512) :: region = "all"
        logical, allocatable :: mask(:,:)       ! the region's cells
        integer            :: n      = 0        ! number of cells of the region
        real(wp)           :: dTo    = 0.0_wp   ! [K] the anomaly of the last update
    end type ocn_fraction_class

    ! [bias<sfx>]: an additive temperature field and a precipitation factor
    ! field, monthly (12) or annual (1 field, all months), read from var_defs
    ! groups (a list of databases, then the par file) on the climate grid
    ! (online remap allowed), each with a weight:
    !   tas = tas + f_tas*dtas,   pr = pr*(1 + f_pr*(fpr - 1))
    type bias_class
        logical               :: active = .false.
        character(len=64)     :: group  = ""
        real(wp)              :: f_tas = 1.0_wp, f_pr = 1.0_wp
        real(wp), allocatable :: dtas(:,:,:)     ! [K]  (nx,ny,12), 0 without a dtas group
        real(wp), allocatable :: fpr(:,:,:)      ! [1]  (nx,ny,12), 1 without a fpr group
    end type bias_class

    type climate_mods_class
        type(climate_caps_class) :: caps
        type(anom_class)         :: anom
        type(ocn_fraction_class) :: ocn_fraction
        type(bias_class)         :: bias
    end type climate_mods_class

    public :: climate_caps_class
    public :: anom_class, ocn_fraction_class, bias_class
    public :: climate_mods_class
    public :: climate_mods_init
    public :: climate_mods_apply
    public :: climate_mods_write_init, climate_mods_write_step

contains

    subroutine climate_mods_init(mods, filename, sfx, idx, caps, domain, grid, reg)
        ! Read the modifier groups (suffix sfx) and check each active modifier's
        ! contract against what the backend supplies (caps). grid: the climate
        ! grid, with the domain's regions on it (reg). Call after the backend init.
        type(climate_mods_class), intent(inout) :: mods
        character(len=*),         intent(in)    :: filename, sfx, domain
        type(indices_class),      intent(in)    :: idx
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid
        type(regions_class),      intent(in)    :: reg

        mods%caps = caps

        call anom_init(mods%anom, filename, "anom"//trim(sfx), idx, caps)
        call ocn_fraction_init(mods%ocn_fraction, filename, "ocn_fraction"//trim(sfx), caps, &
                               mods%anom, reg)
        call bias_init(mods%bias, filename, "bias"//trim(sfx), caps, domain, grid)

    end subroutine climate_mods_init

    subroutine climate_mods_apply(mods, out, idx)
        ! The active modifiers on `out`, in the fixed order. idx is the driver's
        ! registry, advanced to the forcing time.
        type(climate_mods_class), intent(inout) :: mods
        type(climate_out_class),  intent(inout) :: out
        type(indices_class),      intent(in)    :: idx

        if (mods%anom%active .and. .not. mods%caps%anom_input) &
            call anom_apply(mods%anom, out, indices_value(idx, mods%anom%index))

        if (mods%ocn_fraction%active) call ocn_fraction_apply(mods%ocn_fraction, out)

        if (mods%bias%active) call bias_apply(mods%bias, out)

    end subroutine climate_mods_apply

    subroutine climate_mods_write_init(mods, filename)
        ! The static fields of the active modifiers, into the climate's 2D output
        ! file (xc, yc defined): the bias correction, weighted, as annual means.
        type(climate_mods_class), intent(in) :: mods
        character(len=*),         intent(in) :: filename

        if (mods%bias%active) then
            associate(b => mods%bias)
            call nc_write(filename, "bias_dtas_ann", b%f_tas*sum(b%dtas, dim=3)/12.0_wp, &
                          dim1="xc", dim2="yc", units="K", &
                          long_name="Bias correction of the air temperature (ann)")
            call nc_write(filename, "bias_fpr_ann", &
                          1.0_wp + b%f_pr*(sum(b%fpr, dim=3)/12.0_wp - 1.0_wp), &
                          dim1="xc", dim2="yc", units="1", &
                          long_name="Bias correction factor of the precipitation (ann)")
            end associate
        end if

    end subroutine climate_mods_write_init

    subroutine climate_mods_write_step(mods, filename, ncid, n)
        ! Record n of the active modifiers' diagnostics, into an open 2D output
        ! file: the ocean temperature anomaly of ocn_fraction.
        type(climate_mods_class), intent(in) :: mods
        character(len=*),         intent(in) :: filename
        integer,                  intent(in) :: ncid, n

        if (mods%ocn_fraction%active) &
            call nc_write(filename, "ocn_fraction_dto", mods%ocn_fraction%dTo, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units="K", &
                          long_name="Ocean temperature anomaly of ocn_fraction")

    end subroutine climate_mods_write_step

    ! ===== anom ============================================================

    subroutine anom_init(anom, filename, group, idx, caps)
        ! Read [anom<sfx>] (optional group; inactive without it) and check its
        ! contract: a registry index, and a backend that does not supply the smb
        ! (but rembo, which reads the group itself).
        type(anom_class),         intent(inout) :: anom
        character(len=*),         intent(in)    :: filename, group
        type(indices_class),      intent(in)    :: idx
        type(climate_caps_class), intent(in)    :: caps

        anom%group = group
        call nml_validate(filename, def_file, group, defaults_group="anom")
        call nml_read(filename, group, "active", anom%active, defaults_file=def_file, defaults_group="anom")
        if (.not. anom%active) return

        call nml_read(filename, group, "index",  anom%index,  defaults_file=def_file, defaults_group="anom")
        call nml_read(filename, group, "f_ta",   anom%f_ta,   defaults_file=def_file, defaults_group="anom")
        call nml_read(filename, group, "f_to",   anom%f_to,   defaults_file=def_file, defaults_group="anom")
        call nml_read(filename, group, "f_so",   anom%f_so,   defaults_file=def_file, defaults_group="anom")

        if (.not. indices_has(idx, anom%index)) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] index = '"// &
                       trim(anom%index)//"' is not in [indices] names."
            error stop 1
        end if
        if (caps%has_smb .and. .not. caps%anom_input) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] needs a climate that &
                       &does not supply the surface mass balance (an air anomaly cannot reach &
                       &it); climate = "//trim(caps%name)//" supplies it (esm: use the esm &
                       &smb anomalies or surface_method /= climate)."
            error stop 1
        end if

        write(*,*) "climate_mods_init:: "//trim(group)//" active, index "//trim(anom%index)// &
                   "; f_ta/f_to/f_so =", anom%f_ta, anom%f_to, anom%f_so

    end subroutine anom_init

    subroutine anom_apply(anom, out, f)
        ! Homogeneous anomaly on `now` at the index value f: the air temperature
        ! + dTa, the precipitation * exp(beta_p*dTa) (beta_p: the backend's
        ! sensitivity), the ocean + dTo/dSo, as profiles or, without them, at the
        ! shelf base.
        type(anom_class),        intent(in)    :: anom
        type(climate_out_class), intent(inout) :: out
        real(wp),                intent(in)    :: f

        real(wp) :: dTa, dTo, dSo
        integer  :: m

        dTa = f*anom%f_ta
        dTo = f*anom%f_to
        dSo = f*anom%f_so

        associate(now => out%now)

        now%tas    = now%tas + dTa
        now%ta_ann = now%ta_ann + dTa
        if (allocated(now%tsl_ann)) now%tsl_ann = now%tsl_ann + dTa

        do m = 1, size(now%pr, 3)
            now%pr(:,:,m) = now%pr(:,:,m) * exp(now%beta_p*dTa)
        end do
        if (allocated(now%pr_ann)) now%pr_ann = now%pr_ann * exp(now%beta_p*dTa)

        if (out%has_ocn_shelf) then
            now%T_shlf  = now%T_shlf  + dTo
            now%S_shlf  = now%S_shlf  + dSo
            now%dT_shlf = now%dT_shlf + dTo
            now%dS_shlf = now%dS_shlf + dSo
        else if (allocated(now%to_ann) .or. allocated(now%so_ann)) then
            if (allocated(now%to_ann)) now%to_ann = now%to_ann + dTo
            if (allocated(now%so_ann)) now%so_ann = now%so_ann + dSo
        else if (anom%f_to /= 0.0_wp .or. anom%f_so /= 0.0_wp) then
            write(*,*) "climate_mods_apply:: error: ["//trim(anom%group)//"] f_to or f_so /= 0, &
                       &but the climate has no ocean (neither profiles nor at the shelf base); &
                       &set f_to = f_so = 0."
            error stop 1
        end if

        end associate

    end subroutine anom_apply

    ! ===== ocn_fraction ====================================================

    subroutine ocn_fraction_init(ocf, filename, group, caps, anom, reg)
        ! Read [ocn_fraction<sfx>] (optional group; inactive without it), check its
        ! contract (ocean profiles and tsl_ann; no anom ocean temperature
        ! anomaly, which it would replace) and resolve the region on the
        ! climate grid; an empty region is an error.
        type(ocn_fraction_class), intent(inout) :: ocf
        character(len=*),         intent(in)    :: filename, group
        type(climate_caps_class), intent(in)    :: caps
        type(anom_class),         intent(in)    :: anom
        type(regions_class),      intent(in)    :: reg

        ocf%group = group
        call nml_validate(filename, def_file, group, defaults_group="ocn_fraction")
        call nml_read(filename, group, "active", ocf%active, defaults_file=def_file, defaults_group="ocn_fraction")
        if (.not. ocf%active) return

        call nml_read(filename, group, "f_to",   ocf%f_to,   defaults_file=def_file, defaults_group="ocn_fraction")
        call nml_read(filename, group, "region", ocf%region, defaults_file=def_file, defaults_group="ocn_fraction")

        if (caps%has_ocn_shelf .or. .not. caps%has_tsl_ann) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] needs ocean profiles &
                       &and the sea-level air temperature (snapshot backends); climate = "// &
                       trim(caps%name)//" has the ocean at the shelf base or no tsl_ann."
            error stop 1
        end if
        if (anom%active .and. anom%f_to /= 0.0_wp) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] sets the ocean &
                       &temperature anomaly, ["//trim(anom%group)//"] f_to /= 0 adds one; &
                       &set "//trim(anom%group)//".f_to = 0."
            error stop 1
        end if

        ocf%mask = regions_select(reg, trim(ocf%region))
        ocf%n    = count(ocf%mask)
        if (ocf%n == 0) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] region = '"// &
                       trim(ocf%region)//"' selects no cell of the climate grid."
            error stop 1
        end if

        write(*,*) "climate_mods_init:: "//trim(group)//" active, f_to =", ocf%f_to, &
                   "; region '"//trim(ocf%region)//"':", ocf%n, "of", size(ocf%mask), "cells"

    end subroutine ocn_fraction_init

    subroutine ocn_fraction_apply(ocf, out)
        ! The ocean temperature profiles: the reference plus f_to times the mean
        ! sea-level air temperature anomaly of this update over the region.
        type(ocn_fraction_class), intent(inout) :: ocf
        type(climate_out_class),  intent(inout) :: out

        ocf%dTo = ocf%f_to * sum(out%now%tsl_ann - out%ref%tsl_ann, mask=ocf%mask) &
                  / real(ocf%n, wp)
        out%now%to_ann = out%ref%to_ann + ocf%dTo

    end subroutine ocn_fraction_apply

    ! ===== bias ============================================================

    subroutine bias_init(bias, filename, group, caps, domain, grid)
        ! Read [bias<sfx>] (optional group; inactive without it), check its
        ! contract (monthly fields) and load its fields from the var_defs groups
        ! it names (static, 1 or 12 slices) on the climate grid.
        type(bias_class),         intent(inout) :: bias
        character(len=*),         intent(in)    :: filename, group, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid

        character(len=256) :: var_defs(VAR_DEFS_MAX)
        character(len=256) :: grp_dtas, grp_fpr
        integer :: nx, ny

        bias%group = group
        call nml_validate(filename, def_file, group, defaults_group="bias")
        call nml_read(filename, group, "active", bias%active, defaults_file=def_file, defaults_group="bias")
        if (.not. bias%active) return

        if (.not. caps%monthly) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] needs monthly air &
                       &temperature and precipitation; climate = "//trim(caps%name)// &
                       " has annual fields only."
            error stop 1
        end if

        call nml_read(filename, group, "var_defs", var_defs,   init=.true., &
                      defaults_file=def_file, defaults_group="bias")
        call nml_read(filename, group, "dtas",     grp_dtas,   init=.true., &
                      defaults_file=def_file, defaults_group="bias")
        call nml_read(filename, group, "fpr",      grp_fpr,    init=.true., &
                      defaults_file=def_file, defaults_group="bias")
        call nml_read(filename, group, "f_tas",    bias%f_tas, defaults_file=def_file, defaults_group="bias")
        call nml_read(filename, group, "f_pr",     bias%f_pr,  defaults_file=def_file, defaults_group="bias")

        nx = grid%G%nx
        ny = grid%G%ny
        allocate(bias%dtas(nx,ny,12), bias%fpr(nx,ny,12))
        bias%dtas = 0.0_wp
        bias%fpr  = 1.0_wp
        if (len_trim(grp_dtas) > 0) call load_monthly(bias%dtas, trim(grp_dtas))
        if (len_trim(grp_fpr)  > 0) call load_monthly(bias%fpr,  trim(grp_fpr))

        write(*,*) "climate_mods_init:: "//trim(group)//" active, dtas = ", trim(grp_dtas), &
                   ", fpr = ", trim(grp_fpr), ", f_tas = ", bias%f_tas, ", f_pr = ", bias%f_pr
        write(*,*) "climate_mods_init:: "//trim(group)//" range dtas [K]: ", &
                   minval(bias%dtas), maxval(bias%dtas), ", fpr: ", minval(bias%fpr), maxval(bias%fpr)

    contains

        subroutine load_monthly(var, vgroup)
            ! One field for all months, or 12 monthly fields.
            real(wp),         intent(out) :: var(:,:,:)
            character(len=*), intent(in)  :: vgroup
            type(varslice_class) :: vs
            integer :: m, nm

            call varslice_init_nml(vs, trim(var_defs_find(var_defs, filename, vgroup, &
                                                          "["//trim(group)//"]")), &
                                   vgroup, domain=domain, &
                                   grid_name=trim(grid%name), grid=grid)
            call varslice_update(vs)
            nm = size(vs%var, 3)
            if (nm /= 1 .and. nm /= 12) then
                write(*,*) "climate_mods_init:: error: ", vgroup, " must have 1 or 12 &
                           &slices; got ", nm
                error stop 1
            end if
            do m = 1, 12
                var(:,:,m) = vs%var(:,:,min(m,nm),1)
            end do
            call varslice_end(vs)

        end subroutine load_monthly

    end subroutine bias_init

    subroutine bias_apply(bias, out)
        ! Correct the monthly air temperature and precipitation of `out`, and their
        ! annual means (the temperature additively, the precipitation by the
        ! ratio of the annual totals). The reference climate is unchanged.
        type(bias_class),        intent(in)    :: bias
        type(climate_out_class), intent(inout) :: out

        real(wp), allocatable :: pr_sum0(:,:)
        real(wp) :: dtas_ann(size(bias%dtas,1),size(bias%dtas,2))

        dtas_ann = sum(bias%dtas, dim=3) / 12.0_wp

        out%now%tas    = out%now%tas + bias%f_tas*bias%dtas
        out%now%ta_ann = out%now%ta_ann + bias%f_tas*dtas_ann
        if (allocated(out%now%tsl_ann)) out%now%tsl_ann = out%now%tsl_ann + bias%f_tas*dtas_ann

        pr_sum0     = sum(out%now%pr, dim=3)
        out%now%pr  = out%now%pr * (1.0_wp + bias%f_pr*(bias%fpr - 1.0_wp))
        if (allocated(out%now%pr_ann)) then
            where (pr_sum0 > 0.0_wp) out%now%pr_ann = out%now%pr_ann * sum(out%now%pr, dim=3) / pr_sum0
        end if

    end subroutine bias_apply

end module climate_mods
