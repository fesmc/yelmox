module kryos_forcing
    ! Driver-owned forcing: the registry of the scalar series.
    !
    ! Series: named scalars ([forcing] series = "alpha" "beta"), each a group
    ! [series_<name>] (from the par file, else the first of the [forcing]
    ! var_defs files that has it), read by name by their consumers (the terms
    ! of the climate, ice_bnd, the OBM forcing). Two kinds:
    !   - tsgen series (method const, ramp-*, sin, series, feedback ...):
    !     advanced once per step on the forcing time with the response
    !     variable (ice volume);
    !   - derived (method = "derived"): computed each step by the model from
    !     its own state (from = air_anomaly: the mean sea-level air anomaly of
    !     the climate update; obm_north/obm_south: a box temperature of the
    !     bipolar ocean box model), set by its producer with indices_set after
    !     the atmosphere of the step is composed. No tsgen state, no restart
    !     file, no kill switch: it is recomputed.
    ! One registry per run, held by the driver like bsl -- deliberately NOT
    ! part of kryos_domain, so every domain of a multi-domain driver reads the
    ! same series. The registry also holds the kill switch (any series tripped)
    ! and the forcing-increment restarts (a bundle each time a series moves by
    ! more than its restart_every_df, folders restart-<n>); each tsgen series
    ! restarts from series_<name>_restart.nc. (The type and routine names keep
    ! "indices"/"idx" for now.)

    use nml,          only : nml_read, nml_validate
    use yelmox_nml,   only : VAR_DEFS_MAX, var_defs_find, def_file => yelmox_defaults_file
    use ncio,         only : nc_open, nc_close, nc_write, nc_time_index
    use precision,    only : wp, dp
    use tsgen,        only : tsgen_class, tsgen_init, tsgen_update, &
                             tsgen_restart_write, tsgen_restart_read

    implicit none
    private

    integer, parameter :: IDX_MAX = 32   ! max entries of [forcing] series

    ! The definition of a derived series ([series_<name>] method = "derived").
    type derived_spec
        character(len=32) :: from    = ""        ! air_anomaly | obm_north | obm_south
        real(wp)          :: factor  = 1.0_wp    ! air_anomaly: value = factor * mean anomaly
        character(len=64) :: region  = "all"     ! air_anomaly: cells of the mean ([regions] selection)
        logical           :: anomaly = .true.    ! obm_*: the box temperature minus its initial value
    end type derived_spec

    ! One index: a named series, tsgen or derived.
    type index_class
        character(len=64)  :: name
        character(len=32)  :: units            = "1"      ! label of series_<name> in the 1D output
        real(wp)           :: restart_every_df = 0.0_wp   ! [f] |Δf| between restarts; 0 = off
        real(wp)           :: f_last_restart   = 0.0_wp   ! f at the last restart mark
        real(wp)           :: f                = 0.0_wp   ! the current value (both kinds)
        logical            :: derived          = .false.
        type(derived_spec) :: drv
        type(tsgen_class)  :: tsg                         ! tsgen indices only
    end type index_class

    ! The registry.
    type indices_class
        type(index_class), allocatable :: idx(:)
        integer :: counter_restart = 0                  ! forcing-increment restarts written
    end type indices_class

    public :: indices_class, derived_spec
    public :: indices_init, indices_update, indices_value, indices_has
    public :: indices_set, indices_is_derived, indices_derived
    public :: indices_uses_derived, indices_set_obm
    public :: indices_kill, indices_restart_due, indices_restart_fldr
    public :: indices_restart_write, indices_restart_read, indices_write_step

contains

    ! ===== indices registry ===================================================

    subroutine indices_init(idx, path_par, time)
        ! Read [forcing] series and initialise each from its group
        ! [series_<name>] at `time` (the forcing time). The group is read from
        ! the par file if it has it, else from the first of the [forcing]
        ! var_defs files that does (var_defs_find). Per series, besides the
        ! tsgen keys: restart_every_df and units; a derived series (method =
        ! "derived") reads from, factor, region, anomaly instead and starts at
        ! 0. Unset keys from the defaults file (the template series). An absent
        ! [forcing] group or an empty list give an empty registry. On a restart
        ! run the caller then calls indices_restart_read.
        type(indices_class), intent(inout) :: idx
        character(len=*),    intent(in)    :: path_par
        real(dp),            intent(in)    :: time

        character(len=64)  :: names(IDX_MAX)
        character(len=256) :: var_defs(VAR_DEFS_MAX)
        character(len=128) :: group
        character(len=512) :: gfile      ! the file holding the group
        character(len=32)  :: method
        integer :: k, n

        call nml_validate(path_par, def_file, "forcing")
        call nml_read(path_par, "forcing", "series",   names,    init=.true., defaults_file=def_file)
        call nml_read(path_par, "forcing", "var_defs", var_defs, init=.true., defaults_file=def_file)
        n = count(len_trim(names) > 0)

        if (allocated(idx%idx)) deallocate(idx%idx)
        allocate(idx%idx(n))
        idx%counter_restart = 0

        do k = 1, n
            associate(ix => idx%idx(k))
            ix%name = trim(names(k))
            if (any(names(1:k-1) == names(k))) then
                write(*,*) "indices_init:: error: series '"//trim(ix%name)//"' is listed twice in [forcing] series."
                error stop 1
            end if
            group = "series_"//trim(ix%name)
            gfile = var_defs_find(var_defs, path_par, trim(group), "&forcing series")

            call nml_validate(trim(gfile), def_file, trim(group), defaults_group="series")
            call nml_read(trim(gfile), trim(group), "method", method, init=.true., &
                          defaults_file=def_file, defaults_group="series")
            call nml_read(trim(gfile), trim(group), "units", ix%units, init=.true., &
                          defaults_file=def_file, defaults_group="series")
            call nml_read(trim(gfile), trim(group), "restart_every_df", ix%restart_every_df, &
                          defaults_file=def_file, defaults_group="series")

            if (trim(method) == "derived") then
                ix%derived = .true.
                call nml_read(trim(gfile), trim(group), "from",    ix%drv%from,    init=.true., &
                              defaults_file=def_file, defaults_group="series")
                call nml_read(trim(gfile), trim(group), "factor",  ix%drv%factor, &
                              defaults_file=def_file, defaults_group="series")
                call nml_read(trim(gfile), trim(group), "region",  ix%drv%region,  init=.true., &
                              defaults_file=def_file, defaults_group="series")
                call nml_read(trim(gfile), trim(group), "anomaly", ix%drv%anomaly, &
                              defaults_file=def_file, defaults_group="series")
                select case(trim(ix%drv%from))
                    case("air_anomaly", "obm_north", "obm_south")
                    case default
                        write(*,*) "indices_init:: error: ["//trim(group)//"] from = '"// &
                                   trim(ix%drv%from)//"': derived series come from air_anomaly, "// &
                                   "obm_north or obm_south."
                        error stop 1
                end select
                if (ix%restart_every_df > 0.0_wp) then
                    write(*,*) "indices_init:: error: ["//trim(group)//"] restart_every_df > 0 on a "// &
                               "derived series; the forcing-increment restarts follow tsgen series only."
                    error stop 1
                end if
                ix%f = 0.0_wp
                write(*,*) "indices_init:: "//trim(ix%name)//" (derived from "//trim(ix%drv%from)// &
                           "), set each step by the model"
            else
                call tsgen_init(ix%tsg, trim(gfile), time, group=trim(group), &
                                defaults_file=def_file, defaults_group="series")
                ix%f = ix%tsg%f_now
                write(*,*) "indices_init:: "//trim(ix%name)//" ("//trim(ix%tsg%par%method)// &
                           "), f =", ix%f
            end if
            ix%f_last_restart = ix%f
            end associate
        end do
    end subroutine indices_init

    subroutine indices_update(idx, time, var)
        ! Advance every tsgen series to `time` (the forcing time). `var` is the
        ! model response variable (ice volume [Gt]) that feedback methods and
        ! the kill switch consume; ignored by open-loop ones. Derived series
        ! keep their value until their producer sets them this step.
        type(indices_class), intent(inout) :: idx
        real(dp),            intent(in)    :: time
        real(wp),            intent(in)    :: var
        integer :: k
        do k = 1, size(idx%idx)
            if (idx%idx(k)%derived) cycle
            call tsgen_update(idx%idx(k)%tsg, time, var=var)
            idx%idx(k)%f = idx%idx(k)%tsg%f_now
        end do
    end subroutine indices_update

    subroutine indices_set(idx, name, value)
        ! Set the value of the derived series `name` (its producer, once per
        ! step); stops on an unknown or a tsgen series.
        type(indices_class), intent(inout) :: idx
        character(len=*),    intent(in)    :: name
        real(wp),            intent(in)    :: value
        integer :: k
        k = index_position(idx, name)
        if (k == 0) then
            write(*,*) "indices_set:: error: unknown index '"//trim(name)//"' (not in [forcing] series)."
            error stop 1
        end if
        if (.not. idx%idx(k)%derived) then
            write(*,*) "indices_set:: error: index '"//trim(name)//"' is a tsgen series, not derived."
            error stop 1
        end if
        idx%idx(k)%f = value
    end subroutine indices_set

    logical function indices_uses_derived(idx, from) result(uses)
        ! Is any series derived `from` (air_anomaly, obm_north, obm_south)?
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: from
        integer :: k
        uses = .false.
        do k = 1, size(idx%idx)
            if (idx%idx(k)%derived .and. trim(idx%idx(k)%drv%from) == trim(from)) uses = .true.
        end do
    end function indices_uses_derived

    subroutine indices_set_obm(idx, from, theta, theta_init)
        ! Set every series derived `from` an OBM box (obm_north, obm_south) to
        ! its box temperature theta [K], or theta - theta_init with `anomaly`.
        type(indices_class), intent(inout) :: idx
        character(len=*),    intent(in)    :: from
        real(wp),            intent(in)    :: theta, theta_init
        integer :: k
        do k = 1, size(idx%idx)
            associate(ix => idx%idx(k))
            if (.not. (ix%derived .and. trim(ix%drv%from) == trim(from))) cycle
            ix%f = theta
            if (ix%drv%anomaly) ix%f = theta - theta_init
            end associate
        end do
    end subroutine indices_set_obm

    logical function indices_has(idx, name) result(found)
        ! Is `name` a series of the registry?
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        found = index_position(idx, name) > 0
    end function indices_has

    logical function indices_is_derived(idx, name) result(derived)
        ! Is the series `name` derived (set by the model)? Stops on an unknown name.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        derived = idx%idx(index_checked(idx, name, "indices_is_derived"))%derived
    end function indices_is_derived

    function indices_derived(idx, name) result(drv)
        ! The definition (from, factor, region, anomaly) of the derived series
        ! `name`; stops on an unknown or a tsgen series.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        type(derived_spec) :: drv
        integer :: k
        k = index_checked(idx, name, "indices_derived")
        if (.not. idx%idx(k)%derived) then
            write(*,*) "indices_derived:: error: index '"//trim(name)//"' is a tsgen series, not derived."
            error stop 1
        end if
        drv = idx%idx(k)%drv
    end function indices_derived

    real(wp) function indices_value(idx, name) result(f)
        ! Current value of the series `name`; stops on an unknown name.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        f = idx%idx(index_checked(idx, name, "indices_value"))%f
    end function indices_value

    integer function index_checked(idx, name, caller) result(k)
        ! Position of the series `name`; stops with the caller's name if absent.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name, caller
        k = index_position(idx, name)
        if (k == 0) then
            write(*,*) trim(caller)//":: error: unknown index '"//trim(name)// &
                       "' (not in [forcing] series)."
            error stop 1
        end if
    end function index_checked

    integer function index_position(idx, name) result(k)
        ! Position of the series `name` in the registry, 0 if absent.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        do k = 1, size(idx%idx)
            if (trim(idx%idx(k)%name) == trim(name)) return
        end do
        k = 0
    end function index_position

    logical function indices_kill(idx) result(kill)
        ! The tsgen kill switch of any series (set once the response equilibrates
        ! at a forcing bound).
        type(indices_class), intent(in) :: idx
        integer :: k
        kill = .false.
        do k = 1, size(idx%idx)
            if (idx%idx(k)%derived) cycle
            kill = kill .or. idx%idx(k)%tsg%kill
        end do
    end function indices_kill

    logical function indices_restart_due(idx) result(due)
        ! True when any series has moved by more than its restart_every_df since
        ! its own last mark. Advances the bookkeeping (the marks of the series
        ! that moved, one shared counter) as a side effect, so the caller writes
        ! exactly one bundle per trigger. Fires on |Δf| in either direction
        ! (legacy hyster behaviour); an oscillating series triggers on each swing.
        type(indices_class), intent(inout) :: idx
        integer :: k

        due = .false.
        do k = 1, size(idx%idx)
            associate(ix => idx%idx(k))
            if (ix%restart_every_df <= 0.0_wp) cycle
            if (abs(ix%f - ix%f_last_restart) > ix%restart_every_df) then
                ix%f_last_restart = ix%f
                due = .true.
            end if
            end associate
        end do
        if (due) idx%counter_restart = idx%counter_restart + 1
    end function indices_restart_due

    function indices_restart_fldr(idx) result(fldr)
        ! Sequential forcing-increment restart folder (legacy hyster convention):
        ! restart-<n>.
        type(indices_class), intent(in) :: idx
        character(len=512) :: fldr
        character(len=16)  :: cstr
        write(cstr,'(i0)') idx%counter_restart
        fldr = "restart-"//trim(adjustl(cstr))
    end function indices_restart_fldr

    subroutine indices_restart_write(idx, fldr, time)
        ! The state of every tsgen series in the bundle folder fldr:
        ! series_<name>_restart.nc. Derived series are recomputed, not restarted.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: fldr
        real(dp),            intent(in) :: time
        integer :: k
        do k = 1, size(idx%idx)
            if (idx%idx(k)%derived) cycle
            call tsgen_restart_write(trim(index_restart_file(fldr, idx%idx(k)%name)), &
                                     idx%idx(k)%tsg, time)
        end do
    end subroutine indices_restart_write

    subroutine indices_restart_read(idx, fldr)
        ! Restore every tsgen series from the bundle folder fldr (after
        ! indices_init) and reset its restart mark. A series without a restart
        ! file (a bundle from a run without it) keeps its cold-start state, i.e.
        ! starts at the restart time.
        type(indices_class), intent(inout) :: idx
        character(len=*),    intent(in)    :: fldr
        character(len=1024) :: fnm
        logical :: found
        integer :: k

        do k = 1, size(idx%idx)
            associate(ix => idx%idx(k))
            if (ix%derived) cycle
            fnm = index_restart_file(fldr, ix%name)
            inquire(file=trim(fnm), exist=found)
            if (found) then
                call tsgen_restart_read(ix%tsg, trim(fnm))
                ix%f = ix%tsg%f_now
                ix%f_last_restart = ix%f
            else
                write(*,*) "indices_restart_read:: no restart file "//trim(fnm)// &
                           "; series "//trim(ix%name)//" keeps its cold-start state."
            end if
            end associate
        end do
    end subroutine indices_restart_read

    function index_restart_file(fldr, name) result(fnm)
        ! Restart file of one series in the bundle folder fldr.
        character(len=*), intent(in) :: fldr, name
        character(len=1024) :: fnm
        fnm = trim(fldr)//"/series_"//trim(name)//"_restart.nc"
    end function index_restart_file

    subroutine indices_write_step(idx, filename, time)
        ! Append every series to an existing 1D timeseries file at `time`: its
        ! value series_<name> (in its units) and, for tsgen series, the rate
        ! series_<name>_df_dt and, with a response history (feedback methods,
        ! kill switch), the windowed response rate series_<name>_dv_dt. Call
        ! after the file's own 1D write.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: filename
        real(dp),            intent(in) :: time
        integer :: k, ncid, n
        character(len=72) :: vnm

        if (size(idx%idx) == 0) return

        call nc_open(filename, ncid, writable=.true.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)

        do k = 1, size(idx%idx)
            associate(ix => idx%idx(k))
            vnm = "series_"//trim(ix%name)
            call nc_write(filename, trim(vnm), ix%f, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units=trim(ix%units), &
                          long_name="Series "//trim(ix%name))
            if (ix%derived) cycle
            call nc_write(filename, trim(vnm)//"_df_dt", ix%tsg%df_dt, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units=trim(ix%units)//"/yr", &
                          long_name="Rate of change of series "//trim(ix%name))
            if (allocated(ix%tsg%dv_dt)) &
                call nc_write(filename, trim(vnm)//"_dv_dt", ix%tsg%dv_dt_ave, dim1="time", &
                              start=[n], count=[1], ncid=ncid, units="Gt/yr", &
                              long_name="Windowed rate of the response (ice volume) of series "// &
                                        trim(ix%name))
            end associate
        end do

        call nc_close(ncid)
    end subroutine indices_write_step

end module kryos_forcing
