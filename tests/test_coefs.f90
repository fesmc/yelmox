program test_coefs
    ! Checks the coefficient registry (libs/kryos_forcing.f90): two tsgen
    ! coefficients (a ramp and a constant), a scaled one and a derived one
    ! from a small par file. The tsgen coefficient advances with the forcing
    ! time and restarts from its file (coef_<name>_restart.nc, or the former
    ! series_<name>_restart.nc of an older bundle); the scaled coefficient
    ! (listed before its source) follows its source, also after a restart,
    ! without a restart file; the derived coefficient starts at 0, is set by
    ! its producer (coefs_set), is not advanced by coefs_update, has no
    ! restart file; all are written to a 1D file. coefs_varies tells the
    ! constant from the others; coefs_check_used counts the source of a used
    ! scaled coefficient as used.
    !
    ! Build and run from the yelmox root (after make yelmox, which builds
    ! kryos_forcing.o and yelmox_nml.o):
    !     make test_coefs; ./libyelmox/bin/test_coefs.x
    ! (reads input/yelmox_defaults.nml; writes to a temporary folder it removes).

    use precision,     only : wp, dp
    use ncio,          only : nc_create, nc_write_dim, nc_write
    use kryos_forcing, only : coefs_class, derived_spec, coefs_init, coefs_update, &
                              coefs_value, coefs_has, coefs_set, coefs_is_derived, &
                              coefs_derived, coefs_kill, coefs_restart_due, &
                              coefs_restart_write, coefs_restart_read, coefs_write_step, &
                              coefs_varies, coefs_check_used

    implicit none

    character(len=*), parameter :: tmp = "libyelmox/test_coefs_tmp"
    character(len=*), parameter :: par = tmp//"/test_coefs.nml"
    type(coefs_class)   :: coefs
    type(derived_spec)  :: drv
    logical :: found
    integer :: nfail, u

    nfail = 0
    call execute_command_line("mkdir -p "//tmp)

    open(newunit=u, file=par, status="replace", action="write")
    write(u,'(a)') "&forcing"
    write(u,'(a)') "    coefs = ""s"" ""c"" ""dto"" ""k"""
    write(u,'(a)') "/"
    write(u,'(a)') "&coef_s"
    write(u,'(a)') "    method = ""scaled"""
    write(u,'(a)') "    coef   = ""c"""
    write(u,'(a)') "    factor = 2.0"
    write(u,'(a)') "/"
    write(u,'(a)') "&coef_k"
    write(u,'(a)') "    method = ""const"""
    write(u,'(a)') "    f_min  = 1.0"
    write(u,'(a)') "/"
    write(u,'(a)') "&coef_c"
    write(u,'(a)') "    method    = ""ramp-slope"""
    write(u,'(a)') "    f_min     = 0.0"
    write(u,'(a)') "    f_max     = 5.0"
    write(u,'(a)') "    df_dt_max = 0.01"
    write(u,'(a)') "    restart_every_df = 0.5"
    write(u,'(a)') "/"
    write(u,'(a)') "&coef_dto"
    write(u,'(a)') "    method = ""derived"""
    write(u,'(a)') "    from   = ""air_anomaly"""
    write(u,'(a)') "    factor = 0.25"
    write(u,'(a)') "    region = ""zone:open_ocean"""
    write(u,'(a)') "    units  = ""K"""
    write(u,'(a)') "/"
    close(u)

    call coefs_init(coefs, par, 0.0_dp)

    call check(size(coefs%coef) == 4, "four coefficients", nfail)
    call check(coefs_has(coefs, "c") .and. coefs_has(coefs, "dto") .and. .not. coefs_has(coefs, "x"), &
               "coefs_has", nfail)
    call check(.not. coefs_is_derived(coefs, "c") .and. coefs_is_derived(coefs, "dto"), &
               "coefs_is_derived", nfail)
    call check(abs(coefs_value(coefs, "c")) < 1e-12_wp, "c starts at f_min = 0", nfail)
    call check(abs(coefs_value(coefs, "dto")) < 1e-12_wp, "dto starts at 0", nfail)
    call check(abs(coefs_value(coefs, "s")) < 1e-12_wp .and. .not. coefs_is_derived(coefs, "s"), &
               "s = 2 c = 0 at the start, not derived", nfail)
    call check(coefs_varies(coefs, "c") .and. coefs_varies(coefs, "s") .and. &
               coefs_varies(coefs, "dto") .and. .not. coefs_varies(coefs, "k"), "coefs_varies", nfail)
    call coefs_check_used(coefs, [character(len=8) :: "s", "dto", "k"])
    call check(.true., "coefs_check_used: c used as the source of s", nfail)

    drv = coefs_derived(coefs, "dto")
    call check(trim(drv%from) == "air_anomaly" .and. abs(drv%factor - 0.25_wp) < 1e-12_wp .and. &
               trim(drv%region) == "zone:open_ocean" .and. drv%anomaly, "derived spec read", nfail)

    ! One step: the ramp advances 0.01 * 100 yr = 1; the derived coefficient
    ! stays until set.
    call coefs_update(coefs, 100.0_dp, var=0.0_wp)
    call check(abs(coefs_value(coefs, "c") - 1.0_wp) < 1e-9_wp, "c = 1 after 100 yr", nfail)
    call check(abs(coefs_value(coefs, "s") - 2.0_wp) < 1e-9_wp, "s = 2 c after the update", nfail)
    call check(abs(coefs_value(coefs, "dto")) < 1e-12_wp, "dto untouched by update", nfail)
    call coefs_set(coefs, "dto", -2.5_wp)
    call check(abs(coefs_value(coefs, "dto") + 2.5_wp) < 1e-12_wp, "dto set", nfail)
    call coefs_update(coefs, 150.0_dp, var=0.0_wp)
    call check(abs(coefs_value(coefs, "dto") + 2.5_wp) < 1e-12_wp, "dto kept across an update", nfail)
    call check(.not. coefs_kill(coefs), "no kill", nfail)

    ! Forcing-increment restart: c moved by 1.5 > 0.5 since the start.
    call check(coefs_restart_due(coefs), "restart due after |df| > 0.5", nfail)
    call check(.not. coefs_restart_due(coefs), "not due again without motion", nfail)

    ! Restart files: c has one, dto has none; the read restores c.
    call coefs_restart_write(coefs, tmp, 150.0_dp)
    inquire(file=tmp//"/coef_c_restart.nc", exist=found)
    call check(found, "coef_c_restart.nc written", nfail)
    inquire(file=tmp//"/coef_dto_restart.nc", exist=found)
    call check(.not. found, "no restart file for the derived coefficient", nfail)
    inquire(file=tmp//"/coef_s_restart.nc", exist=found)
    call check(.not. found, "no restart file for the scaled coefficient", nfail)

    call coefs_init(coefs, par, 0.0_dp)
    call coefs_restart_read(coefs, tmp)
    call check(abs(coefs_value(coefs, "c") - 1.5_wp) < 1e-9_wp, "c restored from its restart file", nfail)
    call check(abs(coefs_value(coefs, "s") - 3.0_wp) < 1e-9_wp, "s follows the restored c", nfail)
    call check(abs(coefs_value(coefs, "dto")) < 1e-12_wp, "dto starts at 0 after a restart", nfail)

    ! An older bundle: only series_c_restart.nc; the read falls back to it.
    call execute_command_line("mv "//tmp//"/coef_c_restart.nc "//tmp//"/series_c_restart.nc")
    call coefs_init(coefs, par, 0.0_dp)
    call check(abs(coefs_value(coefs, "c")) < 1e-12_wp, "c at its cold start before the read", nfail)
    call coefs_restart_read(coefs, tmp)
    call check(abs(coefs_value(coefs, "c") - 1.5_wp) < 1e-9_wp, "c restored from series_c_restart.nc", nfail)

    ! 1D output: every value, df_dt for the tsgen coefficients only.
    call nc_create(tmp//"/ts.nc")
    call nc_write_dim(tmp//"/ts.nc", "time", x=0.0_dp, dx=1.0_dp, nx=1, unlimited=.true.)
    call coefs_write_step(coefs, tmp//"/ts.nc", 150.0_dp)
    call check(nc_has(tmp//"/ts.nc", "coef_c") .and. nc_has(tmp//"/ts.nc", "coef_dto") .and. &
               nc_has(tmp//"/ts.nc", "coef_s") .and. .not. nc_has(tmp//"/ts.nc", "coef_s_df_dt") .and. &
               nc_has(tmp//"/ts.nc", "coef_c_df_dt") .and. .not. nc_has(tmp//"/ts.nc", "coef_dto_df_dt"), &
               "1D output variables", nfail)

    call execute_command_line("rm -rf "//tmp)

    write(*,*)
    if (nfail == 0) then
        write(*,*) "test_coefs: PASS"
    else
        write(*,*) "test_coefs: FAIL", nfail, "checks"
        stop 1
    end if

contains

    subroutine check(ok, what, nfail)
        logical,          intent(in)    :: ok
        character(len=*), intent(in)    :: what
        integer,          intent(inout) :: nfail
        if (ok) then
            write(*,*) "  ok    "//what
        else
            write(*,*) "  FAIL  "//what
            nfail = nfail + 1
        end if
    end subroutine check

    logical function nc_has(filename, vnm) result(has)
        use ncio, only : nc_exists_var
        character(len=*), intent(in) :: filename, vnm
        has = nc_exists_var(filename, vnm)
    end function nc_has

end program test_coefs
