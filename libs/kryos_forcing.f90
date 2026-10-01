module kryos_forcing
    ! Driver-owned transient forcing: a tsgen series whose scalar output f_now is
    ! mapped onto snapclim's homogeneous anomalies via the [tsforcing] gains
    ! (dTa=f_now*f_ta, dTo=f_now*f_to, dSo=f_now*f_so; only consumed by snapclim
    ! in "anom" mode). Also holds the forcing-increment restart bookkeeping (write
    ! a bundle each |Δf_now| > restart_every_df, folders named restart-<n>) and
    ! surfaces the tsgen kill flag. One instance per run, held by the driver like
    ! bsl -- deliberately NOT part of kryos_domain, so multi-domain drivers choose
    ! shared vs per-domain forcing.

    use nml,          only : nml_read
    use yelmo,        only : wp
    use tsgen,        only : tsgen_class, tsgen_init, tsgen_update, tsgen_write_step, &
                             tsgen_restart_write, tsgen_restart_read

    implicit none
    private

    type tsforcing_class
        logical  :: active           = .false.
        real(wp) :: f_ta             = 1.0_wp    ! [K]   dTa = f_now * f_ta
        real(wp) :: f_to             = 0.0_wp    ! [K]   dTo = f_now * f_to
        real(wp) :: f_so             = 0.0_wp    ! [psu] dSo = f_now * f_so
        real(wp) :: restart_every_df = 0.0_wp    ! [f] |Δf| between restarts; 0 = off
        type(tsgen_class) :: tsg                 ! the forcing series generator
        real(wp) :: dTa = 0.0_wp                 ! derived anomalies (f_now * gain)
        real(wp) :: dTo = 0.0_wp
        real(wp) :: dSo = 0.0_wp
        real(wp) :: f_last_restart  = 0.0_wp     ! f_now at the last forcing-triggered restart
        integer  :: counter_restart = 0          ! number of forcing-triggered restarts written
    end type tsforcing_class

    public :: tsforcing_class, tsforcing_init, tsforcing_update
    public :: tsforcing_restart_due, tsforcing_restart_fldr, tsforcing_kill
    public :: tsforcing_restart_write, tsforcing_restart_read, tsforcing_write_step

contains

    ! ===== transient forcing (tsgen -> snapclim anomalies) ===================

    subroutine tsforcing_init(tsf, path_par, time)
        ! Read [tsforcing] and initialize the tsgen series + derived anomalies. The
        ! series is always initialized (even when inactive) so tsf%tsg is valid for
        ! any diagnostics/output that reference it regardless of `active`; the
        ! active flag only gates whether the forcing is actually applied. On a
        ! restart run the caller then calls tsforcing_restart_read to resume.
        type(tsforcing_class), intent(inout) :: tsf
        character(len=*),      intent(in)    :: path_par
        real(wp),              intent(in)    :: time

        call nml_read(path_par, "tsforcing", "active",           tsf%active)
        call nml_read(path_par, "tsforcing", "f_ta",             tsf%f_ta)
        call nml_read(path_par, "tsforcing", "f_to",             tsf%f_to)
        call nml_read(path_par, "tsforcing", "f_so",             tsf%f_so)
        call nml_read(path_par, "tsforcing", "restart_every_df", tsf%restart_every_df)

        tsf%dTa = 0.0_wp; tsf%dTo = 0.0_wp; tsf%dSo = 0.0_wp
        tsf%counter_restart = 0

        call tsgen_init(tsf%tsg, path_par, time)
        call tsforcing_set_anom(tsf)
        tsf%f_last_restart = tsf%tsg%f_now

        if (tsf%active) &
            write(*,*) "tsforcing:: active (tsgen); f_ta/f_to/f_so =", &
                        tsf%f_ta, tsf%f_to, tsf%f_so
    end subroutine tsforcing_init

    subroutine tsforcing_update(tsf, time, var)
        ! Advance the forcing series one step and refresh the derived anomalies.
        ! `var` is the model response variable (e.g. ice volume [Gt]) that the
        ! feedback methods and the kill switch consume; ignored by open-loop ones.
        type(tsforcing_class), intent(inout) :: tsf
        real(wp),              intent(in)    :: time
        real(wp),              intent(in)    :: var

        if (.not. tsf%active) return
        call tsgen_update(tsf%tsg, time, var=var)
        call tsforcing_set_anom(tsf)
    end subroutine tsforcing_update

    subroutine tsforcing_set_anom(tsf)
        ! Map the scalar forcing f_now onto the three snapclim anomalies.
        type(tsforcing_class), intent(inout) :: tsf
        tsf%dTa = tsf%tsg%f_now * tsf%f_ta
        tsf%dTo = tsf%tsg%f_now * tsf%f_to
        tsf%dSo = tsf%tsg%f_now * tsf%f_so
    end subroutine tsforcing_set_anom

    logical function tsforcing_kill(tsf) result(kill)
        ! The tsgen kill switch (set once the response equilibrates at a bound);
        ! always .false. when forcing is inactive.
        type(tsforcing_class), intent(in) :: tsf
        kill = tsf%active .and. tsf%tsg%kill
    end function tsforcing_kill

    logical function tsforcing_restart_due(tsf) result(due)
        ! True when the forcing has moved by more than restart_every_df since the
        ! last forcing-triggered restart. Advances the bookkeeping (counter +
        ! last-mark) as a side effect, so the caller writes exactly one bundle per
        ! trigger. Fires on |Δf| in either direction (legacy hyster behaviour); on
        ! an oscillating series it therefore triggers on each swing.
        type(tsforcing_class), intent(inout) :: tsf

        due = .false.
        if (.not. tsf%active) return
        if (tsf%restart_every_df .le. 0.0_wp) return

        if (abs(tsf%tsg%f_now - tsf%f_last_restart) .gt. tsf%restart_every_df) then
            tsf%counter_restart = tsf%counter_restart + 1
            tsf%f_last_restart  = tsf%tsg%f_now
            due = .true.
        end if
    end function tsforcing_restart_due

    function tsforcing_restart_fldr(tsf) result(fldr)
        ! Sequential restart folder name (legacy hyster convention): restart-<n>.
        type(tsforcing_class), intent(in) :: tsf
        character(len=512) :: fldr
        character(len=16)  :: cstr
        write(cstr,'(i0)') tsf%counter_restart
        fldr = "restart-"//trim(adjustl(cstr))
    end function tsforcing_restart_fldr

    subroutine tsforcing_restart_write(tsf, fldr, time)
        ! Persist the tsgen state into fldr/tsgen_restart.nc (no-op if inactive).
        type(tsforcing_class), intent(in) :: tsf
        character(len=*),      intent(in) :: fldr
        real(wp),              intent(in) :: time
        if (.not. tsf%active) return
        call tsgen_restart_write(trim(fldr)//"/tsgen_restart.nc", tsf%tsg, time)
    end subroutine tsforcing_restart_write

    subroutine tsforcing_restart_read(tsf, fldr)
        ! Restore the tsgen state from fldr/tsgen_restart.nc so the series resumes,
        ! then refresh the derived anomalies and the restart mark (no-op if
        ! inactive). Call after tsforcing_init on a restart run.
        type(tsforcing_class), intent(inout) :: tsf
        character(len=*),      intent(in)    :: fldr
        if (.not. tsf%active) return
        call tsgen_restart_read(tsf%tsg, trim(fldr)//"/tsgen_restart.nc")
        call tsforcing_set_anom(tsf)
        tsf%f_last_restart = tsf%tsg%f_now
    end subroutine tsforcing_restart_read

    subroutine tsforcing_write_step(tsf, filename, time)
        ! Append the forcing diagnostics (tsgen_f/df_dt/dv_dt) to an existing 1D
        ! timeseries file at `time` (no-op if inactive). Call after the module's
        ! own 1D write so the time index already exists.
        type(tsforcing_class), intent(in) :: tsf
        character(len=*),      intent(in) :: filename
        real(wp),              intent(in) :: time
        if (.not. tsf%active) return
        call tsgen_write_step(filename, tsf%tsg, time)
    end subroutine tsforcing_write_step

end module kryos_forcing
