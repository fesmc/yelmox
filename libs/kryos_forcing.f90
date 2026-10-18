module kryos_forcing
    ! Driver-owned forcing: the indices registry.
    !
    ! Indices: named scalar series ([indices] names = "a" "b"), each a tsgen
    ! group [idx_<name>], advanced once per step on the forcing clock with the
    ! response variable (ice volume) and read by name by their consumers (the
    ! snapesm axes, the climate modifiers). One registry per run, held by the
    ! driver like bsl -- deliberately NOT part of kryos_domain, so every domain
    ! of a multi-domain driver reads the same indices. The registry also holds the
    ! kill switch (any index tripped) and the forcing-increment restarts (a
    ! bundle each time an index moves by more than its restart_every_df,
    ! folders restart-<n>); each index restarts from idx_<name>_restart.nc.

    use nml,          only : nml_read, nml_has_param
    use ncio,         only : nc_open, nc_close, nc_write, nc_time_index
    use precision,    only : wp, dp
    use tsgen,        only : tsgen_class, tsgen_init, tsgen_update, &
                             tsgen_restart_write, tsgen_restart_read

    implicit none
    private

    integer, parameter :: IDX_MAX = 32   ! max entries of [indices] names

    ! One index: a named tsgen series.
    type index_class
        character(len=64) :: name
        character(len=32) :: units            = "1"      ! label of idx_<name> in the 1D output
        real(wp)          :: restart_every_df = 0.0_wp   ! [f] |Δf| between restarts; 0 = off
        real(wp)          :: f_last_restart   = 0.0_wp   ! f_now at the last restart mark
        type(tsgen_class) :: tsg
    end type index_class

    ! The registry.
    type indices_class
        type(index_class), allocatable :: idx(:)
        integer :: counter_restart = 0                  ! forcing-increment restarts written
    end type indices_class

    public :: indices_class, indices_init, indices_update, indices_value, indices_has
    public :: indices_kill, indices_restart_due, indices_restart_fldr
    public :: indices_restart_write, indices_restart_read, indices_write_step

contains

    ! ===== indices registry ===================================================

    subroutine indices_init(idx, path_par, time)
        ! Read [indices] names and initialise each index from its tsgen group
        ! [idx_<name>] at `time` (the forcing clock). Optional per index, besides
        ! the tsgen keys: restart_every_df (default 0) and units (default "1").
        ! An absent [indices] group or empty names give an empty registry. On a
        ! restart run the caller then calls indices_restart_read.
        type(indices_class), intent(inout) :: idx
        character(len=*),    intent(in)    :: path_par
        real(dp),            intent(in)    :: time

        character(len=64)  :: names(IDX_MAX)
        character(len=128) :: group
        integer :: k, n

        names = ""
        if (nml_has_param(path_par, "indices", "names")) &
            call nml_read(path_par, "indices", "names", names)
        n = count(len_trim(names) > 0)

        if (allocated(idx%idx)) deallocate(idx%idx)
        allocate(idx%idx(n))
        idx%counter_restart = 0

        do k = 1, n
            associate(ix => idx%idx(k))
            ix%name = trim(names(k))
            if (any(names(1:k-1) == names(k))) then
                write(*,*) "indices_init:: error: index '"//trim(ix%name)//"' is listed twice in [indices] names."
                error stop 1
            end if
            group = "idx_"//trim(ix%name)

            ix%units = "1"
            if (nml_has_param(path_par, trim(group), "units")) &
                call nml_read(path_par, trim(group), "units", ix%units)
            ix%restart_every_df = 0.0_wp
            if (nml_has_param(path_par, trim(group), "restart_every_df")) &
                call nml_read(path_par, trim(group), "restart_every_df", ix%restart_every_df)

            call tsgen_init(ix%tsg, path_par, time, group=trim(group))
            ix%f_last_restart = ix%tsg%f_now

            write(*,*) "indices_init:: "//trim(ix%name)//" ("//trim(ix%tsg%par%method)// &
                       "), f =", ix%tsg%f_now
            end associate
        end do
    end subroutine indices_init

    subroutine indices_update(idx, time, var)
        ! Advance every index to `time` (the forcing clock). `var` is the model
        ! response variable (ice volume [Gt]) that feedback methods and the kill
        ! switch consume; ignored by open-loop ones.
        type(indices_class), intent(inout) :: idx
        real(dp),            intent(in)    :: time
        real(wp),            intent(in)    :: var
        integer :: k
        do k = 1, size(idx%idx)
            call tsgen_update(idx%idx(k)%tsg, time, var=var)
        end do
    end subroutine indices_update

    logical function indices_has(idx, name) result(found)
        ! Is `name` an index of the registry?
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        found = index_position(idx, name) > 0
    end function indices_has

    real(wp) function indices_value(idx, name) result(f)
        ! Current value f_now of the index `name`; stops on an unknown name.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        integer :: k
        k = index_position(idx, name)
        if (k == 0) then
            write(*,*) "indices_value:: error: unknown index '"//trim(name)// &
                       "' (not in [indices] names)."
            error stop 1
        end if
        f = idx%idx(k)%tsg%f_now
    end function indices_value

    integer function index_position(idx, name) result(k)
        ! Position of the index `name` in the registry, 0 if absent.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: name
        do k = 1, size(idx%idx)
            if (trim(idx%idx(k)%name) == trim(name)) return
        end do
        k = 0
    end function index_position

    logical function indices_kill(idx) result(kill)
        ! The tsgen kill switch of any index (set once the response equilibrates
        ! at a forcing bound).
        type(indices_class), intent(in) :: idx
        integer :: k
        kill = .false.
        do k = 1, size(idx%idx)
            kill = kill .or. idx%idx(k)%tsg%kill
        end do
    end function indices_kill

    logical function indices_restart_due(idx) result(due)
        ! True when any index has moved by more than its restart_every_df since
        ! its own last mark. Advances the bookkeeping (the marks of the indices
        ! that moved, one shared counter) as a side effect, so the caller writes
        ! exactly one bundle per trigger. Fires on |Δf| in either direction
        ! (legacy hyster behaviour); an oscillating series triggers on each swing.
        type(indices_class), intent(inout) :: idx
        integer :: k

        due = .false.
        do k = 1, size(idx%idx)
            associate(ix => idx%idx(k))
            if (ix%restart_every_df <= 0.0_wp) cycle
            if (abs(ix%tsg%f_now - ix%f_last_restart) > ix%restart_every_df) then
                ix%f_last_restart = ix%tsg%f_now
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
        ! The state of every index in the bundle folder fldr: idx_<name>_restart.nc.
        type(indices_class), intent(in) :: idx
        character(len=*),    intent(in) :: fldr
        real(dp),            intent(in) :: time
        integer :: k
        do k = 1, size(idx%idx)
            call tsgen_restart_write(trim(index_restart_file(fldr, idx%idx(k)%name)), &
                                     idx%idx(k)%tsg, time)
        end do
    end subroutine indices_restart_write

    subroutine indices_restart_read(idx, fldr)
        ! Restore every index from the bundle folder fldr (after indices_init) and
        ! reset its restart mark. An index without a restart file (a bundle from a
        ! run without it) keeps its cold-start state, i.e. starts at the restart time.
        type(indices_class), intent(inout) :: idx
        character(len=*),    intent(in)    :: fldr
        character(len=1024) :: fnm
        logical :: found
        integer :: k

        do k = 1, size(idx%idx)
            associate(ix => idx%idx(k))
            fnm = index_restart_file(fldr, ix%name)
            inquire(file=trim(fnm), exist=found)
            if (found) then
                call tsgen_restart_read(ix%tsg, trim(fnm))
                ix%f_last_restart = ix%tsg%f_now
            else
                write(*,*) "indices_restart_read:: no restart file "//trim(fnm)// &
                           "; index "//trim(ix%name)//" keeps its cold-start state."
            end if
            end associate
        end do
    end subroutine indices_restart_read

    function index_restart_file(fldr, name) result(fnm)
        ! Restart file of one index in the bundle folder fldr.
        character(len=*), intent(in) :: fldr, name
        character(len=1024) :: fnm
        fnm = trim(fldr)//"/idx_"//trim(name)//"_restart.nc"
    end function index_restart_file

    subroutine indices_write_step(idx, filename, time)
        ! Append every index to an existing 1D timeseries file at `time`: its value
        ! idx_<name> (in its units), the rate idx_<name>_df_dt and, for indices
        ! with a response history (feedback methods, kill switch), the windowed
        ! response rate idx_<name>_dv_dt. Call after the file's own 1D write.
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
            vnm = "idx_"//trim(ix%name)
            call nc_write(filename, trim(vnm), ix%tsg%f_now, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units=trim(ix%units), &
                          long_name="Index "//trim(ix%name))
            call nc_write(filename, trim(vnm)//"_df_dt", ix%tsg%df_dt, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units=trim(ix%units)//"/yr", &
                          long_name="Rate of change of index "//trim(ix%name))
            if (allocated(ix%tsg%dv_dt)) &
                call nc_write(filename, trim(vnm)//"_dv_dt", ix%tsg%dv_dt_ave, dim1="time", &
                              start=[n], count=[1], ncid=ncid, units="Gt/yr", &
                              long_name="Windowed rate of the response (ice volume) of index "// &
                                        trim(ix%name))
            end associate
        end do

        call nc_close(ncid)
    end subroutine indices_write_step

end module kryos_forcing
