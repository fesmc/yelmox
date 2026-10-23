module climate_terms
    ! The terms of the climate forcing (docs/_forcing-terms.md): each term is
    ! one coefficient of the driver's registry (kryos_forcing) times one pattern,
    ! added to the reference climate on an explicit list of fields:
    !   snapshots -- the difference of two snapshots, snap1 - snap0 (or snap1
    !                itself without snap0), applied by the snapshot backend
    !                at sea level (snapshots.f90 reads these terms);
    !   uniform   -- a constant pattern: tas + f_ta*f (and pr * exp(beta_p*
    !                f_ta*f)), to + f_to*f, so + f_so*f, applied here (f is
    !                the value of the coefficient, f_ta/f_to/f_so the gains);
    !   replace   -- the field takes the value of the coefficient at every depth
    !                in place of the reference (to, so); the other terms on
    !                the field add to it.
    ! The terms are the list [forcing<sfx>] terms, each a group [term_<name>]
    ! read from the par file, else the first of the [forcing] var_defs files
    ! that has it; unset keys from the template term of
    ! input/yelmox_defaults.nml.
    !
    ! The value f of each term's coefficient is evaluated here once per
    ! update (terms_eval), kept in the term, read by the backends (the
    ! snapshot terms, rembo's input) and the applications below, and written
    ! for every term in one place (terms_write_step).
    !
    ! Order within a climate update (yelmox_climate): terms_eval -> the
    ! backend (with the snapshot terms) -> the uniform terms on the atmosphere
    ! (terms_apply_atm) -> the derived coefficients (coefs_set by the climate,
    ! from the composed atmosphere; terms_eval of the terms they drive) -> the
    ! ocean terms (terms_apply_ocn) -> taper -> bias. The ocean application is
    ! idempotent: it rebuilds the ocean from a copy of the backend's
    ! (terms_save_ocean, after the backend), so a driver that sets a derived
    ! coefficient later in the step (the bipolar bridge, the OBM box
    ! temperatures) re-evaluates and re-applies it.
    !
    ! Contracts, checked at init (terms_check) against what the backend
    ! supplies (climate_caps_class), an error naming the term:
    !   - a snapshots pattern needs the snapshot backend (snapshots, or the
    !     rembo ocean); its fields are bound by snap0, snap1 and the reference
    !     (checked by the backend);
    !   - a uniform term on tas or pr needs an atmosphere that drives the
    !     surface (caps_need_atm, the predicate of the transforms too): not
    !     with a backend that supplies the surface mass balance (esm with its
    !     smb), which a post-hoc anomaly cannot reach; rembo takes the uniform
    !     tas terms as its input (dT_sum, terms_uniform_dT) and ignores pr in
    !     them (it computes its own precipitation; a notice at init, an error
    !     only for a term on pr alone, which would have no effect);
    !   - a derived coefficient is computed after the atmosphere, so its term may
    !     touch to and so only;
    !   - a uniform term on to or so needs an ocean (profiles or at the shelf
    !     base); a replace term: to or so only, one per field, ocean profiles
    !     (not the esm ocean at the shelf base).

    use precision,     only : wp
    use nml,           only : nml_read, nml_validate
    use ncio,          only : nc_write
    use climate_out,   only : climate_out_class
    use kryos_forcing, only : coefs_class, coefs_value, coefs_has, coefs_is_derived
    use yelmox_nml,    only : VAR_DEFS_MAX, var_defs_find, def_file => yelmox_defaults_file

    implicit none

    private

    integer, parameter :: TERMS_MAX = 32   ! max entries of [forcing] terms
    integer, parameter :: FLD_MAX   = 4    ! tas, pr, to, so

    ! What a backend supplies, for the contracts (set by climate_init).
    type climate_caps_class
        character(len=16) :: name = ""           ! the backend, for messages
        ! The surface (smb, tsrf) is computed from the atmosphere of
        ! climate_out: monthly tas and pr, now and ref (snapshots; esm without
        ! its smb). False for a backend that supplies the surface mass balance
        ! itself (esm with its smb; rembo, annual fields only), which a change
        ! of tas or pr after the backend cannot reach (caps_need_atm).
        logical :: atm_drives_srf = .false.
        logical :: has_tsl_ann    = .false.      ! sea-level air temperature (now and ref)
        logical :: anom_input     = .false.      ! takes the uniform tas terms as an input (rembo)
        logical :: has_snapshots  = .false.      ! reads the snapshot terms (snapshots; the rembo ocean)
        logical :: has_ocn_prof   = .false.      ! ocean profiles (to_ann, so_ann)
        logical :: has_ocn_shelf  = .false.      ! ocean at the shelf base (T_shlf, S_shlf)
    end type climate_caps_class

    ! One term: coefficient x pattern on a list of fields.
    type term_class
        character(len=64) :: name    = ""
        character(len=64) :: group   = ""
        character(len=64) :: coef    = ""          ! the coefficient; "" = the constant 1
        character(len=16) :: pattern = "snapshots" ! snapshots | uniform | replace
        character(len=64) :: snap0   = ""          ! snapshots: the snapshot at coefficient 0 ("" = none)
        character(len=64) :: snap1   = ""          ! snapshots: the snapshot at coefficient 1
        character(len=8)  :: fields(FLD_MAX) = ""
        logical  :: on_tas = .false., on_pr = .false., on_to = .false., on_so = .false.
        real(wp) :: f_ta = 1.0_wp                  ! [K]   uniform: tas + f_ta*f (pr * exp(beta_p*f_ta*f))
        real(wp) :: f_to = 1.0_wp                  ! [K]   uniform: to + f_to*f
        real(wp) :: f_so = 1.0_wp                  ! [psu] uniform: so + f_so*f
        logical  :: derived = .false.              ! the coefficient is derived (set by the model)
        real(wp) :: f = 0.0_wp                     ! the coefficient value (1 without one), terms_eval
    end type term_class

    type terms_class
        character(len=64) :: group = ""                ! forcing<sfx>
        type(term_class), allocatable :: term(:)
        integer :: n = 0
        ! The backend's ocean of the last update (terms_apply_ocn is idempotent).
        logical :: ocn_saved = .false.
        real(wp), allocatable :: to_base(:,:,:), so_base(:,:,:)
        real(wp), allocatable :: T_base(:,:), S_base(:,:), dT_base(:,:), dS_base(:,:)
    end type terms_class

    public :: climate_caps_class, caps_need_atm, term_class, terms_class
    public :: terms_init, terms_check, terms_eval, terms_save_ocean, terms_apply_atm, terms_apply_ocn
    public :: terms_uniform_dT, terms_write_step, coef_label

contains

    subroutine terms_init(terms, filename, sfx, var_defs, coefs)
        ! Read [forcing<sfx>] terms and each [term_<name>]; resolve the coefficients
        ! against the registry; the backend-dependent contracts follow in
        ! terms_check once the backend is initialised.
        type(terms_class),   intent(inout) :: terms
        character(len=*),    intent(in)    :: filename, sfx
        character(len=*),    intent(in)    :: var_defs(:)
        type(coefs_class),   intent(in)    :: coefs

        character(len=64)  :: names(TERMS_MAX)
        character(len=512) :: gfile
        integer :: k, j

        terms%group = "forcing"//trim(sfx)
        call nml_validate(filename, def_file, trim(terms%group), defaults_group="forcing")
        call nml_read(filename, trim(terms%group), "terms", names, init=.true., &
                      defaults_file=def_file, defaults_group="forcing")
        terms%n = count(len_trim(names) > 0)
        if (allocated(terms%term)) deallocate(terms%term)
        allocate(terms%term(terms%n))

        do k = 1, terms%n
            associate(t => terms%term(k))
            t%name = trim(names(k))
            if (any(names(1:k-1) == names(k))) call term_error(t, "listed twice in ["// &
                                                                 trim(terms%group)//"] terms.")
            t%group = "term_"//trim(t%name)
            gfile = var_defs_find(var_defs, filename, trim(t%group), "&"//trim(terms%group)//" terms")
            call nml_validate(trim(gfile), def_file, trim(t%group), defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "coef",    t%coef,    init=.true., &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "pattern", t%pattern, init=.true., &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "snap0",   t%snap0,   init=.true., &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "snap1",   t%snap1,   init=.true., &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "fields",  t%fields,  init=.true., &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "f_ta",    t%f_ta, &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "f_to",    t%f_to, &
                          defaults_file=def_file, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "f_so",    t%f_so, &
                          defaults_file=def_file, defaults_group="term")

            select case(trim(t%pattern))
                case("snapshots", "uniform", "replace")
                case default
                    call term_error(t, "pattern = '"//trim(t%pattern)// &
                                       "'; must be snapshots, uniform or replace.")
            end select

            if (count(len_trim(t%fields) > 0) == 0) call term_error(t, "fields is empty.")
            do j = 1, FLD_MAX
                select case(trim(t%fields(j)))
                    case("")
                    case("tas"); t%on_tas = .true.
                    case("pr");  t%on_pr  = .true.
                    case("to");  t%on_to  = .true.
                    case("so");  t%on_so  = .true.
                    case default
                        call term_error(t, "field '"//trim(t%fields(j))//"'; must be tas, pr, to or so.")
                end select
            end do

            if (len_trim(t%coef) > 0) then
                if (.not. coefs_has(coefs, t%coef)) &
                    call term_error(t, "coef = '"//trim(t%coef)//"' is not in [forcing] coefs.")
                t%derived = coefs_is_derived(coefs, t%coef)
            end if

            if (trim(t%pattern) == "snapshots") then
                if (len_trim(t%snap1) == 0) call term_error(t, "pattern snapshots needs snap1.")
            else
                if (len_trim(t%snap0) > 0 .or. len_trim(t%snap1) > 0) &
                    call term_error(t, "snap0/snap1 belong to the snapshots pattern.")
            end if
            if (trim(t%pattern) == "replace" .and. (t%on_tas .or. t%on_pr)) &
                call term_error(t, "pattern replace is for the ocean (to, so) only.")
            if (t%derived .and. (t%on_tas .or. t%on_pr)) &
                call term_error(t, "its coefficient "//trim(t%coef)//" is derived (computed after &
                                   &the atmosphere), so it can only drive to or so.")
            if (len_trim(t%coef) == 0 .and. trim(t%pattern) /= "snapshots") &
                call term_error(t, "a "//trim(t%pattern)//" term needs a coefficient.")

            write(*,*) "terms_init:: "//trim(t%group)//": "//trim(t%pattern)// &
                       ", coefficient "//trim(coef_label(t%coef))//", fields "//trim(fields_list(t))
            end associate
        end do

        do k = 1, terms%n
            associate(t => terms%term(k))
            if (trim(t%pattern) /= "replace") cycle
            do j = 1, k-1
                if (trim(terms%term(j)%pattern) /= "replace") cycle
                if ((t%on_to .and. terms%term(j)%on_to) .or. (t%on_so .and. terms%term(j)%on_so)) &
                    call term_error(t, "a second replace term on the same field as "// &
                                       trim(terms%term(j)%group)//".")
            end do
            end associate
        end do

    end subroutine terms_init

    subroutine terms_check(terms, caps)
        ! The contracts that depend on the backend (see the header).
        type(terms_class),        intent(in) :: terms
        type(climate_caps_class), intent(in) :: caps
        integer :: k

        do k = 1, terms%n
            associate(t => terms%term(k))
            select case(trim(t%pattern))
                case("snapshots")
                    if (.not. caps%has_snapshots) &
                        call term_error(t, "pattern snapshots needs the snapshot backend; climate = "// &
                                           trim(caps%name)//" reads no snapshots.")
                case("uniform")
                    if (caps%anom_input) then
                        ! rembo: the terms on tas are its input (not applied
                        ! after it); it computes its own precipitation.
                        if (t%on_pr) then
                            if (.not. (t%on_tas .or. t%on_to .or. t%on_so)) &
                                call term_error(t, "lists only pr, which climate = "//trim(caps%name)// &
                                                   " ignores (it computes its own precipitation): &
                                                   &the term has no effect.")
                            write(*,*) "terms_check:: ["//trim(t%group)//"] climate = "//trim(caps%name)// &
                                       " ignores pr in this term (it computes its own precipitation)."
                        end if
                    else if (t%on_tas .or. t%on_pr) then
                        call caps_need_atm(caps, "["//trim(t%group)//"] (a uniform term on tas/pr)")
                    end if
                    if ((t%on_to .or. t%on_so) .and. .not. (caps%has_ocn_prof .or. caps%has_ocn_shelf)) &
                        call term_error(t, "lists to or so, but climate = "//trim(caps%name)// &
                                           " has no ocean (neither profiles nor at the shelf base).")
                case("replace")
                    if (.not. caps%has_ocn_prof) &
                        call term_error(t, "pattern replace needs ocean profiles; climate = "// &
                                           trim(caps%name)//" has none (esm: the ocean at the shelf base).")
            end select
            end associate
        end do

    end subroutine terms_check

    subroutine caps_need_atm(caps, who)
        ! The one predicate of the terms and transforms that change the
        ! atmosphere (tas, pr) after the backend: stop unless it drives the
        ! surface (caps%atm_drives_srf). `who` names the term or transform.
        type(climate_caps_class), intent(in) :: caps
        character(len=*),         intent(in) :: who

        if (caps%atm_drives_srf) return

        if (caps%anom_input) then
            write(*,*) "climate_init:: error: "//trim(who)//" changes the atmosphere after the &
                       &backend; climate = "//trim(caps%name)//" computes its own atmosphere and &
                       &surface mass balance (annual fields) and takes only the uniform terms on &
                       &tas, as its input."
        else
            write(*,*) "climate_init:: error: "//trim(who)//" changes the atmosphere after the &
                       &backend, which cannot reach the surface: climate = "//trim(caps%name)// &
                       " supplies the surface mass balance itself (esm: use its smb anomalies, &
                       &or surface_method /= climate)."
        end if
        error stop 1

    end subroutine caps_need_atm

    ! ===== application ======================================================

    subroutine terms_eval(terms, coefs, derived)
        ! The value f of each term's coefficient from the registry (1 for a
        ! term without one): every term, once per update before the backend.
        ! With derived = .true., after the derived coefficients are set (by
        ! the climate, or later in the step by the driver): only the uniform
        ! and replace terms with a derived coefficient, applied after that
        ! (terms_apply_ocn); a snapshot term keeps the value the backend applied.
        type(terms_class), intent(inout) :: terms
        type(coefs_class), intent(in)    :: coefs
        logical,           intent(in), optional :: derived

        logical :: only_derived
        integer :: k

        only_derived = .false.
        if (present(derived)) only_derived = derived

        do k = 1, terms%n
            associate(t => terms%term(k))
            if (only_derived .and. (.not. t%derived .or. trim(t%pattern) == "snapshots")) cycle
            t%f = 1.0_wp
            if (len_trim(t%coef) > 0) t%f = coefs_value(coefs, t%coef)
            end associate
        end do

    end subroutine terms_eval

    subroutine terms_apply_atm(terms, out)
        ! The uniform terms on the atmosphere: tas, ta_ann, tsl_ann + f_ta*f;
        ! pr, pr_ann * exp(beta_p*f_ta*f) (beta_p: the backend's sensitivity).
        type(terms_class),       intent(in)    :: terms
        type(climate_out_class), intent(inout) :: out

        real(wp) :: dTa
        integer  :: k, m

        do k = 1, terms%n
            associate(t => terms%term(k), now => out%now)
            if (trim(t%pattern) /= "uniform") cycle
            if (.not. (t%on_tas .or. t%on_pr)) cycle
            dTa = t%f_ta*t%f
            if (t%on_tas) then
                now%tas    = now%tas + dTa
                now%ta_ann = now%ta_ann + dTa
                if (allocated(now%tsl_ann)) now%tsl_ann = now%tsl_ann + dTa
            end if
            if (t%on_pr) then
                do m = 1, size(now%pr, 3)
                    now%pr(:,:,m) = now%pr(:,:,m) * exp(now%beta_p*dTa)
                end do
                if (allocated(now%pr_ann)) now%pr_ann = now%pr_ann * exp(now%beta_p*dTa)
            end if
            end associate
        end do

    end subroutine terms_apply_atm

    subroutine terms_save_ocean(terms, out)
        ! Keep the backend's ocean of this update (profiles, or at the shelf
        ! base), the base the ocean terms are applied to.
        type(terms_class),       intent(inout) :: terms
        type(climate_out_class), intent(in)    :: out

        terms%ocn_saved = .true.
        if (out%has_ocn_shelf) then
            terms%T_base  = out%now%T_shlf
            terms%S_base  = out%now%S_shlf
            terms%dT_base = out%now%dT_shlf
            terms%dS_base = out%now%dS_shlf
        else
            if (allocated(out%now%to_ann)) terms%to_base = out%now%to_ann
            if (allocated(out%now%so_ann)) terms%so_base = out%now%so_ann
        end if

    end subroutine terms_save_ocean

    subroutine terms_apply_ocn(terms, out)
        ! The ocean terms on the saved base: a replace term sets the field to
        ! its coefficient at every depth, the uniform ones add f_to*f / f_so*f; the
        ! profiles to_ann/so_ann, or at the shelf base T_shlf/S_shlf and
        ! dT_shlf/dS_shlf (the climate has the ocean a term needs: terms_check).
        ! Idempotent (see the header).
        type(terms_class),       intent(in)    :: terms
        type(climate_out_class), intent(inout) :: out

        logical  :: any_ocn
        integer  :: k

        any_ocn = .false.
        do k = 1, terms%n
            if (trim(terms%term(k)%pattern) == "snapshots") cycle
            if (terms%term(k)%on_to .or. terms%term(k)%on_so) any_ocn = .true.
        end do
        if (.not. any_ocn) return

        if (.not. terms%ocn_saved) then
            write(*,*) "terms_apply_ocn:: error: no saved ocean (terms_save_ocean first)."
            error stop 1
        end if

        ! The base.
        if (out%has_ocn_shelf) then
            out%now%T_shlf  = terms%T_base
            out%now%S_shlf  = terms%S_base
            out%now%dT_shlf = terms%dT_base
            out%now%dS_shlf = terms%dS_base
        else
            if (allocated(terms%to_base)) out%now%to_ann = terms%to_base
            if (allocated(terms%so_base)) out%now%so_ann = terms%so_base
        end if

        ! Replace terms first (profiles only, checked at init).
        do k = 1, terms%n
            associate(t => terms%term(k))
            if (trim(t%pattern) /= "replace") cycle
            if (t%on_to) out%now%to_ann = t%f
            if (t%on_so) out%now%so_ann = t%f
            end associate
        end do

        ! Then the uniform ones.
        do k = 1, terms%n
            associate(t => terms%term(k))
            if (trim(t%pattern) /= "uniform") cycle
            if (.not. (t%on_to .or. t%on_so)) cycle
            if (out%has_ocn_shelf) then
                if (t%on_to) then
                    out%now%T_shlf  = out%now%T_shlf  + t%f_to*t%f
                    out%now%dT_shlf = out%now%dT_shlf + t%f_to*t%f
                end if
                if (t%on_so) then
                    out%now%S_shlf  = out%now%S_shlf  + t%f_so*t%f
                    out%now%dS_shlf = out%now%dS_shlf + t%f_so*t%f
                end if
            else
                if (t%on_to) out%now%to_ann = out%now%to_ann + t%f_to*t%f
                if (t%on_so) out%now%so_ann = out%now%so_ann + t%f_so*t%f
            end if
            end associate
        end do

    end subroutine terms_apply_ocn

    real(wp) function terms_uniform_dT(terms) result(dT)
        ! The sum f_ta*f over the uniform terms on tas: the air temperature
        ! anomaly a backend takes as its input (rembo dT_sum).
        type(terms_class),   intent(in) :: terms
        integer :: k
        dT = 0.0_wp
        do k = 1, terms%n
            associate(t => terms%term(k))
            if (trim(t%pattern) == "uniform" .and. t%on_tas) dT = dT + t%f_ta*t%f
            end associate
        end do
    end function terms_uniform_dT

    subroutine terms_write_step(terms, filename, ncid, n)
        ! Record n of every term into an open 2D output file: the value of
        ! its coefficient, term_<name> (1 for a snapshot term without one).
        type(terms_class), intent(in) :: terms
        character(len=*),  intent(in) :: filename
        integer,           intent(in) :: ncid, n
        integer :: k
        do k = 1, terms%n
            associate(t => terms%term(k))
            call nc_write(filename, "term_"//trim(t%name), t%f, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units="1", &
                          long_name="Coefficient "//trim(coef_label(t%coef))//" of the "// &
                                    trim(t%pattern)//" term "//trim(t%name))
            end associate
        end do
    end subroutine terms_write_step

    ! ===== helpers ==========================================================

    function fields_list(t) result(s)
        type(term_class), intent(in) :: t
        character(len=48) :: s
        integer :: j
        s = ""
        do j = 1, FLD_MAX
            if (len_trim(t%fields(j)) > 0) s = trim(s)//" "//trim(t%fields(j))
        end do
        s = adjustl(s)
    end function fields_list

    function coef_label(name) result(s)
        ! The coefficient name, or "1" for none.
        character(len=*), intent(in) :: name
        character(len=64) :: s
        s = "1"
        if (len_trim(name) > 0) s = name
    end function coef_label

    subroutine term_error(t, msg)
        type(term_class), intent(in) :: t
        character(len=*), intent(in) :: msg
        write(*,*) "climate_terms:: error: ["//trim(t%group)//"] "//trim(msg)
        error stop 1
    end subroutine term_error

end module climate_terms
