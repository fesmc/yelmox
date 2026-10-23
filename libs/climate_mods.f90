module climate_mods
    ! Climate transforms: static changes applied to a backend's `climate_out`
    ! after the backend and the terms (climate_terms), in a fixed order
    ! (docs/_forcing-terms.md, step 5):
    !   taper -- the atmosphere towards the reference with a static weight
    !            field (regional taper)
    !   bias  -- bias correction of the atmosphere (static, or sliced at the
    !            forcing time when its groups have a time axis)
    ! Each is an optional namelist group with `active` (default off), suffixed
    ! by the domain (taper<sfx>, bias<sfx>); unset keys (an absent group: all)
    ! from the templates taper, bias of input/yelmox_defaults.nml. Their
    ! fields are varslice groups, read from the par file, else the first of
    ! the [forcing] var_defs files that has them.
    !
    ! Contract: each transform states what it touches in `now` (`ref` is never
    ! touched); climate_mods_init checks every active transform against what
    ! the backend supplies (climate_caps_class) and stops on a mismatch, never
    ! a silent partial application. Both change the atmosphere after the
    ! backend, so both need an atmosphere that drives the surface (monthly tas
    ! and pr, now and ref: caps_need_atm, the predicate of the uniform terms
    ! too); refused for a backend that supplies the surface mass balance
    ! itself (esm with its smb; rembo).
    !   taper: now = ref + w*(now - ref) for tas, pr, ta_ann, pr_ann and
    !         tsl_ann (those the climate has), w in [0, 1] a static field on
    !         the climate grid, the reference at the current surface
    !         (taper_apply); with outside = "mean", ref + w*(now - ref) +
    !         (1 - w)*the mean anomaly (temperature; precipitation: ratio of the
    !         totals) over the cells with w = 1, per month.
    !   bias: tas, ta_ann, tsl_ann + dtas; pr * fpr, pr_ann by the ratio of the
    !         annual totals.

    use precision,     only : wp, dp
    use nml,           only : nml_read, nml_validate
    use ncio,          only : nc_write
    use coordinates,   only : grid_class
    use varslice,      only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use climate_out,   only : climate_out_class
    use climate_terms, only : climate_caps_class, caps_need_atm
    use yelmox_nml,    only : VAR_DEFS_MAX, var_defs_find, forcing_var_defs, def_file => yelmox_defaults_file

    implicit none

    private

    ! [taper<sfx>]: the atmosphere towards the reference with a weight field w
    ! in [0, 1] (one 2D field), read from a var_defs group on the climate grid
    ! (online remap allowed): now = ref + w*(now - ref); w = 1 keeps the
    ! climate, w = 0 is the reference (outside = "ref") or the reference plus
    ! the mean anomaly of the cells with w = 1 (outside = "mean").
    type taper_class
        logical               :: active  = .false.
        character(len=64)     :: group   = ""
        character(len=256)    :: weight  = ""      ! the var_defs group of w
        character(len=8)      :: outside = "ref"   ! ref | mean: the climate where w = 0
        real(wp), allocatable :: w(:,:)            ! [1] (nx,ny)
        logical,  allocatable :: full(:,:)         ! w = 1: the cells of the mean anomaly
    end type taper_class

    ! [bias<sfx>]: an additive temperature field and a precipitation factor
    ! field, monthly (12) or annual (1 field, all months), read from var_defs
    ! groups on the climate grid (online remap allowed), each with a weight:
    !   tas = tas + f_tas*dtas,   pr = pr*(1 + f_pr*(fpr - 1))
    ! A group with a time axis is kept open and sliced at the forcing time at
    ! every update (bias_update).
    type bias_class
        logical               :: active = .false.
        character(len=64)     :: group  = ""
        real(wp)              :: f_tas = 1.0_wp, f_pr = 1.0_wp
        logical               :: has_dtas = .false., has_fpr = .false.
        logical               :: tdep = .false.   ! a group with a time axis
        type(varslice_class)  :: vs_dtas, vs_fpr
        real(wp), allocatable :: dtas(:,:,:)     ! [K]  (nx,ny,12), 0 without a dtas group
        real(wp), allocatable :: fpr(:,:,:)      ! [1]  (nx,ny,12), 1 without a fpr group
    end type bias_class

    type climate_mods_class
        type(climate_caps_class) :: caps
        type(taper_class)        :: taper
        type(bias_class)         :: bias
    end type climate_mods_class

    public :: taper_class, bias_class
    public :: climate_mods_class
    public :: climate_mods_init
    public :: climate_mods_apply
    public :: climate_mods_write_init, climate_mods_write_step

contains

    subroutine climate_mods_init(mods, filename, sfx, caps, domain, grid, time)
        ! Read the transform groups (suffix sfx) and check each active one's
        ! contract against what the backend supplies (caps). grid: the climate
        ! grid. time: the forcing time, for the bias groups with a time axis.
        ! Call after the backend init.
        type(climate_mods_class), intent(inout) :: mods
        character(len=*),         intent(in)    :: filename, sfx, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid
        real(dp),                 intent(in)    :: time

        mods%caps = caps

        call taper_init(mods%taper, filename, "taper"//trim(sfx), caps, domain, grid)
        call bias_init(mods%bias, filename, "bias"//trim(sfx), caps, domain, grid, time)

    end subroutine climate_mods_init

    subroutine climate_mods_apply(mods, out, time)
        ! The active transforms on `out`, in the fixed order; the bias fields
        ! sliced at `time` (the forcing time) first when they have a time axis.
        type(climate_mods_class), intent(inout) :: mods
        type(climate_out_class),  intent(inout) :: out
        real(dp),                 intent(in)    :: time

        if (mods%taper%active) call taper_apply(mods%taper, out)

        if (mods%bias%active) then
            if (mods%bias%tdep) call bias_update(mods%bias, time)
            call bias_apply(mods%bias, out)
        end if

    end subroutine climate_mods_apply

    subroutine climate_mods_write_init(mods, filename)
        ! The static fields of the active transforms, into the climate's 2D
        ! output file (xc, yc defined): the taper weight; the bias correction,
        ! weighted, as annual means (of the start, for a time-dependent bias).
        type(climate_mods_class), intent(in) :: mods
        character(len=*),         intent(in) :: filename

        if (mods%taper%active) &
            call nc_write(filename, "taper_w", mods%taper%w, dim1="xc", dim2="yc", units="1", &
                          long_name="Weight of the climate towards the reference (taper)")

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
        ! Record n of the active transforms' diagnostics, into an open 2D
        ! output file: the weighted annual temperature correction of a
        ! time-dependent bias, bias_dtas_ann_now (the taper and a static bias
        ! are written once by climate_mods_write_init).
        type(climate_mods_class), intent(in) :: mods
        character(len=*),         intent(in) :: filename
        integer,                  intent(in) :: ncid, n

        if (mods%bias%active .and. mods%bias%tdep) then
            associate(b => mods%bias)
            call nc_write(filename, "bias_dtas_ann_now", b%f_tas*sum(b%dtas, dim=3)/12.0_wp, &
                          dim1="xc", dim2="yc", dim3="time", start=[1,1,n], &
                          count=[size(b%dtas,1),size(b%dtas,2),1], ncid=ncid, units="K", &
                          long_name="Bias correction of the air temperature at the forcing time (ann)", &
                          prec_abs=1e-3_wp)
            end associate
        end if

    end subroutine climate_mods_write_step

    ! ===== taper ===========================================================

    subroutine taper_init(tap, filename, group, caps, domain, grid)
        ! Read [taper<sfx>] (optional group; inactive without it), check its
        ! contract (caps_need_atm) and load the weight from its
        ! var_defs group (static, one field) on the climate grid; a weight
        ! outside [0, 1] (or missing) is an error.
        type(taper_class),        intent(inout) :: tap
        character(len=*),         intent(in)    :: filename, group, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid

        character(len=256) :: var_defs(VAR_DEFS_MAX)
        type(varslice_class) :: vs

        tap%group = group
        call nml_validate(filename, def_file, group, defaults_group="taper")
        call nml_read(filename, group, "active", tap%active, defaults_file=def_file, defaults_group="taper")
        if (.not. tap%active) return

        call caps_need_atm(caps, "["//trim(group)//"]")

        var_defs = forcing_var_defs(filename)
        call nml_read(filename, group, "weight",   tap%weight, defaults_file=def_file, defaults_group="taper")
        call nml_read(filename, group, "outside",  tap%outside, defaults_file=def_file, defaults_group="taper")
        select case(trim(tap%outside))
            case("ref", "mean")
                ! ok
            case default
                write(*,*) "climate_mods_init:: error: ["//trim(group)//"] outside must be ref &
                           &or mean; got '"//trim(tap%outside)//"'."
                error stop 1
        end select

        call varslice_init_nml(vs, trim(var_defs_find(var_defs, filename, trim(tap%weight), &
                                                      "["//trim(group)//"] weight")), &
                               trim(tap%weight), domain=domain, grid_name=trim(grid%name), grid=grid)
        call varslice_update(vs)
        if (size(vs%var, 3) /= 1 .or. size(vs%var, 4) /= 1) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] weight "//trim(tap%weight)// &
                       " must be one 2D field; got ", size(vs%var, 3), " x ", size(vs%var, 4), " slices."
            error stop 1
        end if
        tap%w = vs%var(:,:,1,1)
        call varslice_end(vs)

        if (any(tap%w < 0.0_wp .or. tap%w > 1.0_wp)) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] weight "//trim(tap%weight)// &
                       " must lie in [0, 1] (missing values: set fill in the group); range: ", &
                       minval(tap%w), maxval(tap%w)
            error stop 1
        end if

        tap%full = tap%w >= 1.0_wp
        if (trim(tap%outside) == "mean" .and. .not. any(tap%full)) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] outside = mean needs &
                       &cells with w = 1 (the mean anomaly); weight "//trim(tap%weight)// &
                       " has none."
            error stop 1
        end if

        write(*,*) "climate_mods_init:: "//trim(group)//" active, weight = "//trim(tap%weight)// &
                   ", range: ", minval(tap%w), maxval(tap%w), ", outside = "//trim(tap%outside), &
                   ", cells with w = 1: ", count(tap%full)

    end subroutine taper_init

    subroutine taper_apply(tap, out)
        ! The atmosphere of `out` towards its reference: now = ref + w*(now - ref),
        ! for each field the climate has (now and ref), with the reference at
        ! the surface of `now` (now%z_srf): a reference given at its own surface
        ! (ref%z_srf: snapshots, the elevation of the reference snapshot) is
        ! moved there with the monthly lapse rate and beta_p, as the backend
        ! moves the climate (tas - lapse*dz, pr*exp(-beta_p*lapse*dz)), so that
        ! w = 0 gives the reference climate at the model surface; without
        ! ref%z_srf (esm) the reference is already at the current surface.
        ! The sea-level temperature tsl_ann needs no move.
        ! outside = "mean": where w < 1 the anomaly tends to the mean anomaly
        ! over the cells with w = 1 instead of zero, per month: the temperature
        ! difference dT(m) and the precipitation ratio of the totals f(m);
        ! ta_ann and tsl_ann with the mean of dT, pr_ann with the ratio of the
        ! annual totals.
        type(taper_class),       intent(in)    :: tap
        type(climate_out_class), intent(inout) :: out

        real(wp), allocatable :: dz(:,:), tas_r(:,:,:), pr_r(:,:,:), pr_ann_r(:,:)
        real(wp), allocatable :: dT(:), f(:)
        real(wp) :: dT_ann, f_ann
        integer  :: m, nx, ny, nm

        associate(now => out%now, ref => out%ref, w => tap%w, full => tap%full)

        nx = size(w, 1)
        ny = size(w, 2)
        allocate(dz(nx,ny))
        dz = 0.0_wp
        if (allocated(ref%z_srf)) dz = now%z_srf - ref%z_srf

        dT_ann = 0.0_wp
        if (allocated(now%tas) .and. allocated(ref%tas)) then
            nm = size(now%tas, 3)
            allocate(tas_r(nx,ny,nm), dT(nm))
            dT = 0.0_wp
            do m = 1, nm
                tas_r(:,:,m) = ref%tas(:,:,m) - out%lapse(m)*dz
                if (trim(tap%outside) == "mean") &
                    dT(m) = sum(now%tas(:,:,m) - tas_r(:,:,m), mask=full) / real(count(full), wp)
                now%tas(:,:,m) = tas_r(:,:,m) + w*(now%tas(:,:,m) - tas_r(:,:,m)) + (1.0_wp - w)*dT(m)
            end do
            dT_ann = sum(dT) / real(nm, wp)
            if (allocated(now%ta_ann)) &
                now%ta_ann = sum(tas_r, dim=3)/real(nm, wp) + &
                             w*(now%ta_ann - sum(tas_r, dim=3)/real(nm, wp)) + (1.0_wp - w)*dT_ann
        end if
        if (allocated(now%pr) .and. allocated(ref%pr)) then
            nm = size(now%pr, 3)
            allocate(pr_r(nx,ny,nm), f(nm))
            f = 1.0_wp
            do m = 1, nm
                pr_r(:,:,m) = ref%pr(:,:,m) * exp(-now%beta_p*out%lapse(m)*dz)
            end do
            if (trim(tap%outside) == "mean") then
                do m = 1, nm
                    if (sum(pr_r(:,:,m), mask=full) > 0.0_wp) &
                        f(m) = sum(now%pr(:,:,m), mask=full) / sum(pr_r(:,:,m), mask=full)
                end do
            end if
            f_ann = 1.0_wp
            if (trim(tap%outside) == "mean" .and. sum(sum(pr_r, dim=3), mask=full) > 0.0_wp) &
                f_ann = sum(sum(now%pr, dim=3), mask=full) / sum(sum(pr_r, dim=3), mask=full)
            do m = 1, nm
                now%pr(:,:,m) = w*now%pr(:,:,m) + (1.0_wp - w)*pr_r(:,:,m)*f(m)
            end do
            if (allocated(now%pr_ann)) then
                pr_ann_r   = sum(pr_r, dim=3)/real(nm, wp)*365.0_wp
                now%pr_ann = w*now%pr_ann + (1.0_wp - w)*pr_ann_r*f_ann
            end if
        end if
        if (allocated(now%tsl_ann) .and. allocated(ref%tsl_ann)) &
            now%tsl_ann = ref%tsl_ann + w*(now%tsl_ann - ref%tsl_ann) + (1.0_wp - w)*dT_ann

        end associate

    end subroutine taper_apply

    ! ===== bias ============================================================

    subroutine bias_init(bias, filename, group, caps, domain, grid, time)
        ! Read [bias<sfx>] (optional group; inactive without it), check its
        ! contract (caps_need_atm) and open its fields from the var_defs
        ! groups it names (1 or 12 slices) on the climate grid, sliced at
        ! `time` when they have a time axis.
        type(bias_class),         intent(inout) :: bias
        character(len=*),         intent(in)    :: filename, group, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid
        real(dp),                 intent(in)    :: time

        character(len=256) :: var_defs(VAR_DEFS_MAX)
        character(len=256) :: grp_dtas, grp_fpr
        integer :: nx, ny

        bias%group = group
        call nml_validate(filename, def_file, group, defaults_group="bias")
        call nml_read(filename, group, "active", bias%active, defaults_file=def_file, defaults_group="bias")
        if (.not. bias%active) return

        call caps_need_atm(caps, "["//trim(group)//"]")

        var_defs = forcing_var_defs(filename)
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
        bias%has_dtas = len_trim(grp_dtas) > 0
        bias%has_fpr  = len_trim(grp_fpr)  > 0
        if (bias%has_dtas) call open_field(bias%vs_dtas, trim(grp_dtas))
        if (bias%has_fpr)  call open_field(bias%vs_fpr,  trim(grp_fpr))
        bias%tdep = (bias%has_dtas .and. bias%vs_dtas%par%with_time) .or. &
                    (bias%has_fpr  .and. bias%vs_fpr%par%with_time)
        call bias_update(bias, time)

        write(*,*) "climate_mods_init:: "//trim(group)//" active, dtas = ", trim(grp_dtas), &
                   ", fpr = ", trim(grp_fpr), ", f_tas = ", bias%f_tas, ", f_pr = ", bias%f_pr, &
                   merge(", sliced at the forcing time", "                            ", bias%tdep)
        write(*,*) "climate_mods_init:: "//trim(group)//" range dtas [K]: ", &
                   minval(bias%dtas), maxval(bias%dtas), ", fpr: ", minval(bias%fpr), maxval(bias%fpr)

    contains

        subroutine open_field(vs, vgroup)
            type(varslice_class), intent(inout) :: vs
            character(len=*),     intent(in)    :: vgroup
            call varslice_init_nml(vs, trim(var_defs_find(var_defs, filename, vgroup, &
                                                          "["//trim(group)//"]")), &
                                   vgroup, domain=domain, grid_name=trim(grid%name), grid=grid)
            if (.not. vs%par%with_time) call varslice_update(vs)
        end subroutine open_field

    end subroutine bias_init

    subroutine bias_update(bias, time)
        ! The bias fields at `time`: the groups with a time axis sliced there
        ! (linear between the bracketing records, the first or last outside),
        ! then one field for all months, or 12 monthly fields.
        type(bias_class), intent(inout) :: bias
        real(dp),         intent(in)    :: time

        if (bias%has_dtas) then
            if (bias%vs_dtas%par%with_time) call varslice_update(bias%vs_dtas, [time], method="extrap")
            call fill_monthly(bias%dtas, bias%vs_dtas, "dtas")
        end if
        if (bias%has_fpr) then
            if (bias%vs_fpr%par%with_time) call varslice_update(bias%vs_fpr, [time], method="extrap")
            call fill_monthly(bias%fpr, bias%vs_fpr, "fpr")
        end if

    contains

        subroutine fill_monthly(var, vs, what)
            real(wp),             intent(out) :: var(:,:,:)
            type(varslice_class), intent(in)  :: vs
            character(len=*),     intent(in)  :: what
            integer :: m, nm
            nm = size(vs%var, 3)
            if (nm /= 1 .and. nm /= 12) then
                write(*,*) "climate_mods:: error: ["//trim(bias%group)//"] "//what// &
                           " must have 1 or 12 slices; got ", nm
                error stop 1
            end if
            do m = 1, 12
                var(:,:,m) = vs%var(:,:,min(m,nm),1)
            end do
        end subroutine fill_monthly

    end subroutine bias_update

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
