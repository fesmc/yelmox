module kryos_forcing
    ! Driver-owned forcing: the registry of the coefficients.
    !
    ! A coefficient is a scalar of the forcing time, named by its symbol
    ! ([forcing] coefs = "alpha" "beta"), each defined by a group
    ! [coef_<name>] (from the par file, else the first of the [forcing]
    ! var_defs files that has it), read by name by its consumers (the terms
    ! of the climate, ice_bnd, the OBM forcing). Three kinds:
    !   - tsgen coefficients (method const, ramp-*, sin, series, feedback ...):
    !     advanced once per step on the forcing time with the response
    !     variable (ice volume);
    !   - scaled (method = "scaled"): factor x a tsgen coefficient (coef),
    !     e.g. a regional amplification of a global anomaly, set right after
    !     the tsgen coefficients (a term on the atmosphere may use it);
    !   - derived (method = "derived"): computed each step by the model from
    !     its own state (from = air_anomaly: the mean sea-level air anomaly of
    !     the climate update; obm_north/obm_south: a box temperature of the
    !     bipolar ocean box model), set by its producer with coefs_set after
    !     the atmosphere of the step is composed.
    ! Scaled and derived coefficients have no tsgen state, no restart file and
    ! no kill switch: they are recomputed.
    ! One registry per run, held by the driver like bsl -- deliberately NOT
    ! part of kryos_domain, so every domain of a multi-domain driver reads the
    ! same coefficients. The registry also holds the kill switch (any
    ! coefficient tripped) and the forcing-increment restarts (a bundle each
    ! time a coefficient moves by more than its restart_every_df, folders
    ! restart-<n>); each tsgen coefficient restarts from
    ! coef_<name>_restart.nc (a bundle with the former name
    ! series_<name>_restart.nc is still read).

    use nml,          only : nml_read, nml_validate
    use yelmox_nml,   only : VAR_DEFS_MAX, var_defs_find, def_file => yelmox_defaults_file
    use ncio,         only : nc_write, nc_write_attr
    use precision,    only : wp, dp
    use tsgen,        only : tsgen_class, tsgen_init, tsgen_update, &
                             tsgen_restart_write, tsgen_restart_read

    implicit none
    private

    integer, parameter :: COEF_MAX = 32   ! max entries of [forcing] coefs

    ! The definition of a derived coefficient ([coef_<name>] method = "derived").
    type derived_spec
        character(len=32) :: from    = ""        ! air_anomaly | obm_north | obm_south
        real(wp)          :: factor  = 1.0_wp    ! air_anomaly: value = factor * mean anomaly
        character(len=64) :: region  = "all"     ! air_anomaly: cells of the mean ([regions] selection)
        logical           :: anomaly = .true.    ! obm_*: the box temperature minus its initial value
    end type derived_spec

    ! One coefficient: tsgen, scaled or derived.
    type coef_class
        character(len=64)  :: name
        character(len=8)   :: kind             = "tsgen"  ! tsgen | scaled | derived
        character(len=32)  :: units            = "1"      ! label of <name> in forcing_ts.nc
        real(wp)           :: restart_every_df = 0.0_wp   ! [f] |Δf| between restarts; 0 = off
        real(wp)           :: f_last_restart   = 0.0_wp   ! f at the last restart mark
        real(wp)           :: f                = 0.0_wp   ! the current value (every kind)
        type(tsgen_class)  :: tsg                         ! tsgen only
        character(len=64)  :: source           = ""       ! scaled: the tsgen coefficient it scales
        real(wp)           :: factor           = 1.0_wp   ! scaled: f = factor * f(source)
        integer            :: k_source         = 0        ! scaled: the position of the source
        type(derived_spec) :: drv                         ! derived only
    end type coef_class

    ! The registry.
    type coefs_class
        type(coef_class), allocatable :: coef(:)
        integer :: counter_restart = 0                  ! forcing-increment restarts written
    end type coefs_class

    public :: coefs_class, derived_spec
    public :: coefs_init, coefs_update, coefs_value, coefs_has
    public :: coefs_set, coefs_is_derived, coefs_derived, coefs_varies, coefs_check_used
    public :: coefs_uses_derived, coefs_set_obm
    public :: coefs_kill, coefs_restart_due, coefs_restart_fldr
    public :: coefs_restart_write, coefs_restart_read, coefs_write_step, coefs_write_attrs

contains

    ! ===== coefficient registry ===============================================

    subroutine coefs_init(coefs, path_par, time)
        ! Read [forcing] coefs and initialise each from its group
        ! [coef_<name>] at `time` (the forcing time). The group is read from
        ! the par file if it has it, else from the first of the [forcing]
        ! var_defs files that does (var_defs_find). Per coefficient, besides
        ! the tsgen keys: restart_every_df and units; a scaled coefficient
        ! (method = "scaled") reads coef and factor instead (its source, a
        ! tsgen coefficient of the registry, resolved once all are read); a
        ! derived coefficient (method = "derived") reads from, factor, region,
        ! anomaly and starts at 0. Unset keys from the defaults file (the
        ! template [coef]).
        ! An absent [forcing] group or an empty list give an empty registry.
        ! On a restart run the caller then calls coefs_restart_read.
        type(coefs_class), intent(inout) :: coefs
        character(len=*),  intent(in)    :: path_par
        real(dp),          intent(in)    :: time

        character(len=64)  :: names(COEF_MAX)
        character(len=256) :: var_defs(VAR_DEFS_MAX)
        character(len=128) :: group
        character(len=512) :: gfile      ! the file holding the group
        character(len=32)  :: method
        integer :: k, n

        call nml_validate(path_par, def_file, "forcing")
        call nml_read(path_par, "forcing", "coefs",    names,    init=.true., defaults_file=def_file)
        call nml_read(path_par, "forcing", "var_defs", var_defs, init=.true., defaults_file=def_file)
        n = count(len_trim(names) > 0)

        if (allocated(coefs%coef)) deallocate(coefs%coef)
        allocate(coefs%coef(n))
        coefs%counter_restart = 0

        do k = 1, n
            associate(c => coefs%coef(k))
            c%name = trim(names(k))
            if (any(names(1:k-1) == names(k))) then
                write(*,*) "coefs_init:: error: coefficient '"//trim(c%name)//"' is listed twice in [forcing] coefs."
                error stop 1
            end if
            if (trim(c%name) == "time" .or. trim(c%name) == "bsl") then
                write(*,*) "coefs_init:: error: '"//trim(c%name)//"' is not a coefficient name: "// &
                           "forcing_ts.nc holds it beside the coefficients."
                error stop 1
            end if
            group = "coef_"//trim(c%name)
            gfile = var_defs_find(var_defs, path_par, trim(group), "&forcing coefs")

            call nml_validate(trim(gfile), def_file, trim(group), defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "method", method, init=.true., &
                          defaults_file=def_file, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "units", c%units, init=.true., &
                          defaults_file=def_file, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "restart_every_df", c%restart_every_df, &
                          defaults_file=def_file, defaults_group="coef")

            if (trim(method) == "scaled") then
                c%kind = "scaled"
                call nml_read(trim(gfile), trim(group), "coef",   c%source, init=.true., &
                              defaults_file=def_file, defaults_group="coef")
                call nml_read(trim(gfile), trim(group), "factor", c%factor, &
                              defaults_file=def_file, defaults_group="coef")
                if (c%restart_every_df > 0.0_wp) then
                    write(*,*) "coefs_init:: error: ["//trim(group)//"] restart_every_df > 0 on a "// &
                               "scaled coefficient; the forcing-increment restarts follow tsgen coefficients only."
                    error stop 1
                end if
            else if (trim(method) == "derived") then
                c%kind = "derived"
                call nml_read(trim(gfile), trim(group), "from",    c%drv%from,    init=.true., &
                              defaults_file=def_file, defaults_group="coef")
                call nml_read(trim(gfile), trim(group), "factor",  c%drv%factor, &
                              defaults_file=def_file, defaults_group="coef")
                call nml_read(trim(gfile), trim(group), "region",  c%drv%region,  init=.true., &
                              defaults_file=def_file, defaults_group="coef")
                call nml_read(trim(gfile), trim(group), "anomaly", c%drv%anomaly, &
                              defaults_file=def_file, defaults_group="coef")
                select case(trim(c%drv%from))
                    case("air_anomaly", "obm_north", "obm_south")
                    case default
                        write(*,*) "coefs_init:: error: ["//trim(group)//"] from = '"// &
                                   trim(c%drv%from)//"': derived coefficients come from air_anomaly, "// &
                                   "obm_north or obm_south."
                        error stop 1
                end select
                if (c%restart_every_df > 0.0_wp) then
                    write(*,*) "coefs_init:: error: ["//trim(group)//"] restart_every_df > 0 on a "// &
                               "derived coefficient; the forcing-increment restarts follow tsgen coefficients only."
                    error stop 1
                end if
                c%f = 0.0_wp
                write(*,*) "coefs_init:: "//trim(c%name)//" (derived from "//trim(c%drv%from)// &
                           "), set each step by the model"
            else
                call tsgen_init(c%tsg, trim(gfile), time, group=trim(group), &
                                defaults_file=def_file, defaults_group="coef")
                c%f = c%tsg%f_now
                write(*,*) "coefs_init:: "//trim(c%name)//" ("//trim(c%tsg%par%method)// &
                           "), f =", c%f
            end if
            c%f_last_restart = c%f
            end associate
        end do

        ! The source of each scaled coefficient: a tsgen coefficient of the
        ! registry (a derived one is set after the atmosphere, too late).
        do k = 1, n
            associate(c => coefs%coef(k))
            if (c%kind /= "scaled") cycle
            c%k_source = coef_position(coefs, c%source)
            if (c%k_source == 0) then
                write(*,*) "coefs_init:: error: [coef_"//trim(c%name)//"] coef = '"//trim(c%source)// &
                           "' is not in [forcing] coefs."
                error stop 1
            end if
            if (coefs%coef(c%k_source)%kind /= "tsgen") then
                write(*,*) "coefs_init:: error: [coef_"//trim(c%name)//"] coef = '"//trim(c%source)// &
                           "' is a "//trim(coefs%coef(c%k_source)%kind)//" coefficient; a scaled "// &
                           "coefficient scales a tsgen one."
                error stop 1
            end if
            end associate
        end do
        call scaled_update(coefs)
        do k = 1, n
            associate(c => coefs%coef(k))
            if (c%kind /= "scaled") cycle
            c%f_last_restart = c%f
            write(*,*) "coefs_init:: "//trim(c%name)//" (scaled: ", c%factor, " x "//trim(c%source)// &
                       "), f =", c%f
            end associate
        end do
    end subroutine coefs_init

    subroutine scaled_update(coefs)
        ! Every scaled coefficient from the current value of its source.
        type(coefs_class), intent(inout) :: coefs
        integer :: k
        do k = 1, size(coefs%coef)
            associate(c => coefs%coef(k))
            if (c%kind /= "scaled") cycle
            c%f = c%factor * coefs%coef(c%k_source)%f
            end associate
        end do
    end subroutine scaled_update

    subroutine coefs_update(coefs, time, var)
        ! Advance every tsgen coefficient to `time` (the forcing time). `var` is
        ! the model response variable (ice volume [Gt]) that feedback methods
        ! and the kill switch consume; ignored by open-loop ones. Then the
        ! scaled coefficients from them. Derived coefficients keep their value
        ! until their producer sets them this step.
        type(coefs_class), intent(inout) :: coefs
        real(dp),          intent(in)    :: time
        real(wp),          intent(in)    :: var
        integer :: k
        do k = 1, size(coefs%coef)
            if (coefs%coef(k)%kind /= "tsgen") cycle
            call tsgen_update(coefs%coef(k)%tsg, time, var=var)
            coefs%coef(k)%f = coefs%coef(k)%tsg%f_now
        end do
        call scaled_update(coefs)
    end subroutine coefs_update

    subroutine coefs_set(coefs, name, value)
        ! Set the value of the derived coefficient `name` (its producer, once
        ! per step); stops on an unknown or a tsgen coefficient.
        type(coefs_class), intent(inout) :: coefs
        character(len=*),  intent(in)    :: name
        real(wp),          intent(in)    :: value
        integer :: k
        k = coef_position(coefs, name)
        if (k == 0) then
            write(*,*) "coefs_set:: error: unknown coefficient '"//trim(name)//"' (not in [forcing] coefs)."
            error stop 1
        end if
        if (coefs%coef(k)%kind /= "derived") then
            write(*,*) "coefs_set:: error: coefficient '"//trim(name)//"' is a "// &
                       trim(coefs%coef(k)%kind)//" coefficient, not derived."
            error stop 1
        end if
        coefs%coef(k)%f = value
    end subroutine coefs_set

    logical function coefs_uses_derived(coefs, from) result(uses)
        ! Is any coefficient derived `from` (air_anomaly, obm_north, obm_south)?
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: from
        integer :: k
        uses = .false.
        do k = 1, size(coefs%coef)
            if (coefs%coef(k)%kind == "derived" .and. trim(coefs%coef(k)%drv%from) == trim(from)) uses = .true.
        end do
    end function coefs_uses_derived

    subroutine coefs_set_obm(coefs, from, theta, theta_init)
        ! Set every coefficient derived `from` an OBM box (obm_north,
        ! obm_south) to its box temperature theta [K], or theta - theta_init
        ! with `anomaly`.
        type(coefs_class), intent(inout) :: coefs
        character(len=*),  intent(in)    :: from
        real(wp),          intent(in)    :: theta, theta_init
        integer :: k
        do k = 1, size(coefs%coef)
            associate(c => coefs%coef(k))
            if (.not. (c%kind == "derived" .and. trim(c%drv%from) == trim(from))) cycle
            c%f = theta
            if (c%drv%anomaly) c%f = theta - theta_init
            end associate
        end do
    end subroutine coefs_set_obm

    logical function coefs_has(coefs, name) result(found)
        ! Is `name` a coefficient of the registry?
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name
        found = coef_position(coefs, name) > 0
    end function coefs_has

    logical function coefs_is_derived(coefs, name) result(derived)
        ! Is the coefficient `name` derived (set by the model)? Stops on an unknown name.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name
        derived = coefs%coef(coef_checked(coefs, name, "coefs_is_derived"))%kind == "derived"
    end function coefs_is_derived

    logical function coefs_varies(coefs, name) result(varies)
        ! May the coefficient `name` change in time? A tsgen coefficient unless
        ! it is a const without rate, noise or window; a scaled one as its
        ! source; a derived one always. Stops on an unknown name.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name
        integer :: k
        k = coef_checked(coefs, name, "coefs_varies")
        if (coefs%coef(k)%kind == "scaled") k = coefs%coef(k)%k_source
        if (coefs%coef(k)%kind == "derived") then
            varies = .true.
            return
        end if
        associate(p => coefs%coef(k)%tsg%par)
        varies = trim(p%method) /= "const" .or. p%df_dt_max /= 0.0_wp .or. p%sigma /= 0.0_wp &
                 .or. (p%win_t1 > p%win_t0 .and. p%win_fac /= 1.0_wp)
        end associate
    end function coefs_varies

    subroutine coefs_check_used(coefs, used)
        ! Every coefficient of the registry is used: named by a consumer (the
        ! list `used`, collected by the driver from the terms of its domains,
        ! ice_bnd and the OBM forcing) or the source of a scaled coefficient
        ! that is. Stops on the first unused one: listed in [forcing] coefs
        ! but read by nothing, it is a sign of a setup half changed.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: used(:)
        logical, allocatable :: is_used(:)
        integer :: k
        allocate(is_used(size(coefs%coef)))
        do k = 1, size(coefs%coef)
            is_used(k) = any(used == coefs%coef(k)%name)
        end do
        do k = 1, size(coefs%coef)
            if (coefs%coef(k)%kind == "scaled" .and. is_used(k)) is_used(coefs%coef(k)%k_source) = .true.
        end do
        do k = 1, size(coefs%coef)
            if (is_used(k)) cycle
            write(*,*) "coefs_check_used:: error: coefficient '"//trim(coefs%coef(k)%name)// &
                       "' is in [forcing] coefs, but no term, ice_bnd or the OBM forcing of the run "// &
                       "uses it; remove it from the list (or add the term that uses it)."
            error stop 1
        end do
    end subroutine coefs_check_used

    function coefs_derived(coefs, name) result(drv)
        ! The definition (from, factor, region, anomaly) of the derived
        ! coefficient `name`; stops on an unknown or a tsgen coefficient.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name
        type(derived_spec) :: drv
        integer :: k
        k = coef_checked(coefs, name, "coefs_derived")
        if (coefs%coef(k)%kind /= "derived") then
            write(*,*) "coefs_derived:: error: coefficient '"//trim(name)//"' is a "// &
                       trim(coefs%coef(k)%kind)//" coefficient, not derived."
            error stop 1
        end if
        drv = coefs%coef(k)%drv
    end function coefs_derived

    real(wp) function coefs_value(coefs, name) result(f)
        ! Current value of the coefficient `name`; stops on an unknown name.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name
        f = coefs%coef(coef_checked(coefs, name, "coefs_value"))%f
    end function coefs_value

    integer function coef_checked(coefs, name, caller) result(k)
        ! Position of the coefficient `name`; stops with the caller's name if absent.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name, caller
        k = coef_position(coefs, name)
        if (k == 0) then
            write(*,*) trim(caller)//":: error: unknown coefficient '"//trim(name)// &
                       "' (not in [forcing] coefs)."
            error stop 1
        end if
    end function coef_checked

    integer function coef_position(coefs, name) result(k)
        ! Position of the coefficient `name` in the registry, 0 if absent.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: name
        do k = 1, size(coefs%coef)
            if (trim(coefs%coef(k)%name) == trim(name)) return
        end do
        k = 0
    end function coef_position

    logical function coefs_kill(coefs) result(kill)
        ! The tsgen kill switch of any coefficient (set once the response
        ! equilibrates at a forcing bound).
        type(coefs_class), intent(in) :: coefs
        integer :: k
        kill = .false.
        do k = 1, size(coefs%coef)
            if (coefs%coef(k)%kind /= "tsgen") cycle
            kill = kill .or. coefs%coef(k)%tsg%kill
        end do
    end function coefs_kill

    logical function coefs_restart_due(coefs) result(due)
        ! True when any coefficient has moved by more than its
        ! restart_every_df since its own last mark. Advances the bookkeeping
        ! (the marks of the coefficients that moved, one shared counter) as a
        ! side effect, so the caller writes exactly one bundle per trigger.
        ! Fires on |Δf| in either direction (legacy hyster behaviour); an
        ! oscillating coefficient triggers on each swing.
        type(coefs_class), intent(inout) :: coefs
        integer :: k

        due = .false.
        do k = 1, size(coefs%coef)
            associate(c => coefs%coef(k))
            if (c%restart_every_df <= 0.0_wp) cycle
            if (abs(c%f - c%f_last_restart) > c%restart_every_df) then
                c%f_last_restart = c%f
                due = .true.
            end if
            end associate
        end do
        if (due) coefs%counter_restart = coefs%counter_restart + 1
    end function coefs_restart_due

    function coefs_restart_fldr(coefs) result(fldr)
        ! Sequential forcing-increment restart folder (legacy hyster convention):
        ! restart-<n>.
        type(coefs_class), intent(in) :: coefs
        character(len=512) :: fldr
        character(len=16)  :: cstr
        write(cstr,'(i0)') coefs%counter_restart
        fldr = "restart-"//trim(adjustl(cstr))
    end function coefs_restart_fldr

    subroutine coefs_restart_write(coefs, fldr, time)
        ! The state of every tsgen coefficient in the bundle folder fldr:
        ! coef_<name>_restart.nc. Scaled and derived coefficients are
        ! recomputed, not restarted.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: fldr
        real(dp),          intent(in) :: time
        integer :: k
        do k = 1, size(coefs%coef)
            if (coefs%coef(k)%kind /= "tsgen") cycle
            call tsgen_restart_write(trim(coef_restart_file(fldr, coefs%coef(k)%name)), &
                                     coefs%coef(k)%tsg, time)
        end do
    end subroutine coefs_restart_write

    subroutine coefs_restart_read(coefs, fldr)
        ! Restore every tsgen coefficient from the bundle folder fldr (after
        ! coefs_init) and reset its restart mark. A bundle written before the
        ! rename holds series_<name>_restart.nc instead of coef_<name>_restart.nc;
        ! it is read when the new name is absent. A coefficient without a
        ! restart file (a bundle from a run without it) keeps its cold-start
        ! state, i.e. starts at the restart time. The scaled coefficients then
        ! follow their restored sources.
        type(coefs_class), intent(inout) :: coefs
        character(len=*),  intent(in)    :: fldr
        character(len=1024) :: fnm, fnm_old
        logical :: found, found_old
        integer :: k

        do k = 1, size(coefs%coef)
            associate(c => coefs%coef(k))
            if (c%kind /= "tsgen") cycle
            fnm = coef_restart_file(fldr, c%name)
            inquire(file=trim(fnm), exist=found)
            if (.not. found) then
                fnm_old = trim(fldr)//"/series_"//trim(c%name)//"_restart.nc"
                inquire(file=trim(fnm_old), exist=found_old)
                if (found_old) then
                    write(*,*) "coefs_restart_read:: note: reading the former restart name "// &
                               trim(fnm_old)//" (no "//trim(fnm)//")."
                    fnm   = fnm_old
                    found = .true.
                end if
            end if
            if (found) then
                call tsgen_restart_read(c%tsg, trim(fnm))
                c%f = c%tsg%f_now
                c%f_last_restart = c%f
            else
                write(*,*) "coefs_restart_read:: no restart file "//trim(fnm)// &
                           "; coefficient "//trim(c%name)//" keeps its cold-start state."
            end if
            end associate
        end do
        call scaled_update(coefs)
        do k = 1, size(coefs%coef)
            if (coefs%coef(k)%kind == "scaled") coefs%coef(k)%f_last_restart = coefs%coef(k)%f
        end do
    end subroutine coefs_restart_read

    function coef_restart_file(fldr, name) result(fnm)
        ! Restart file of one coefficient in the bundle folder fldr.
        character(len=*), intent(in) :: fldr, name
        character(len=1024) :: fnm
        fnm = trim(fldr)//"/coef_"//trim(name)//"_restart.nc"
    end function coef_restart_file

    subroutine coefs_write_step(coefs, filename, ncid, n)
        ! Write every coefficient to record n of an open 1D timeseries file
        ! (the run's forcing_ts.nc): its value <name> (in its units) and, for
        ! tsgen coefficients, the rate <name>_df_dt and, with a response
        ! history (feedback methods, kill switch), the windowed response rate
        ! <name>_dv_dt. Their attributes: coefs_write_attrs.
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: filename
        integer,           intent(in) :: ncid, n
        integer :: k

        do k = 1, size(coefs%coef)
            associate(c => coefs%coef(k))
            call nc_write(filename, trim(c%name), c%f, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units=trim(c%units), &
                          long_name="Coefficient "//trim(c%name))
            if (c%kind /= "tsgen") cycle
            call nc_write(filename, trim(c%name)//"_df_dt", c%tsg%df_dt, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units=trim(c%units)//"/yr", &
                          long_name="Rate of change of coefficient "//trim(c%name))
            if (allocated(c%tsg%dv_dt)) &
                call nc_write(filename, trim(c%name)//"_dv_dt", c%tsg%dv_dt_ave, dim1="time", &
                              start=[n], count=[1], ncid=ncid, units="Gt/yr", &
                              long_name="Windowed rate of the response (ice volume) of coefficient "// &
                                        trim(c%name))
            end associate
        end do
    end subroutine coefs_write_step

    subroutine coefs_write_attrs(coefs, filename)
        ! The attributes of each coefficient's value in a closed 1D file
        ! written by coefs_write_step: kind, and method (tsgen), coef + factor
        ! (scaled) or from (derived).
        type(coefs_class), intent(in) :: coefs
        character(len=*),  intent(in) :: filename
        integer :: k

        do k = 1, size(coefs%coef)
            associate(c => coefs%coef(k))
            call nc_write_attr(filename, trim(c%name), "kind", trim(c%kind))
            select case(trim(c%kind))
                case("tsgen")
                    call nc_write_attr(filename, trim(c%name), "method", trim(c%tsg%par%method))
                case("scaled")
                    call nc_write_attr(filename, trim(c%name), "coef", trim(c%source))
                    call nc_write_attr(filename, trim(c%name), "factor", c%factor)
                case("derived")
                    call nc_write_attr(filename, trim(c%name), "from", trim(c%drv%from))
            end select
            end associate
        end do
    end subroutine coefs_write_attrs

end module kryos_forcing
