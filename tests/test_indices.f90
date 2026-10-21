program test_indices
    ! Checks the indices registry (libs/kryos_forcing.f90): a tsgen index and
    ! a derived one from a small par file. The tsgen index advances with the
    ! forcing time and restarts from its file; the derived index starts at 0,
    ! is set by its producer (indices_set), is not advanced by indices_update,
    ! has no restart file, and both are written to a 1D file.
    !
    ! Build and run from the yelmox root: make test_indices; ./libyelmox/bin/test_indices.x
    ! (reads input/yelmox_defaults.nml; writes to a temporary folder it removes).

    use precision,     only : wp, dp
    use ncio,          only : nc_create, nc_write_dim, nc_write
    use kryos_forcing, only : indices_class, derived_spec, indices_init, indices_update, &
                              indices_value, indices_has, indices_set, indices_is_derived, &
                              indices_derived, indices_kill, indices_restart_due, &
                              indices_restart_write, indices_restart_read, indices_write_step

    implicit none

    character(len=*), parameter :: tmp = "libyelmox/test_indices_tmp"
    character(len=*), parameter :: par = tmp//"/test_indices.nml"
    type(indices_class) :: idx
    type(derived_spec)  :: drv
    logical :: found
    integer :: nfail, u

    nfail = 0
    call execute_command_line("mkdir -p "//tmp)

    open(newunit=u, file=par, status="replace", action="write")
    write(u,'(a)') "&indices"
    write(u,'(a)') "    names = ""c"" ""dto"""
    write(u,'(a)') "/"
    write(u,'(a)') "&idx_c"
    write(u,'(a)') "    method    = ""ramp-slope"""
    write(u,'(a)') "    f_min     = 0.0"
    write(u,'(a)') "    f_max     = 5.0"
    write(u,'(a)') "    df_dt_max = 0.01"
    write(u,'(a)') "    restart_every_df = 0.5"
    write(u,'(a)') "/"
    write(u,'(a)') "&idx_dto"
    write(u,'(a)') "    method = ""derived"""
    write(u,'(a)') "    from   = ""air_anomaly"""
    write(u,'(a)') "    factor = 0.25"
    write(u,'(a)') "    region = ""zone:open_ocean"""
    write(u,'(a)') "    units  = ""K"""
    write(u,'(a)') "/"
    close(u)

    call indices_init(idx, par, 0.0_dp)

    call check(size(idx%idx) == 2, "two indices", nfail)
    call check(indices_has(idx, "c") .and. indices_has(idx, "dto") .and. .not. indices_has(idx, "x"), &
               "indices_has", nfail)
    call check(.not. indices_is_derived(idx, "c") .and. indices_is_derived(idx, "dto"), &
               "indices_is_derived", nfail)
    call check(abs(indices_value(idx, "c")) < 1e-12_wp, "c starts at f_min = 0", nfail)
    call check(abs(indices_value(idx, "dto")) < 1e-12_wp, "dto starts at 0", nfail)

    drv = indices_derived(idx, "dto")
    call check(trim(drv%from) == "air_anomaly" .and. abs(drv%factor - 0.25_wp) < 1e-12_wp .and. &
               trim(drv%region) == "zone:open_ocean" .and. drv%anomaly, "derived spec read", nfail)

    ! One step: the ramp advances 0.01 * 100 yr = 1; the derived index stays
    ! until set.
    call indices_update(idx, 100.0_dp, var=0.0_wp)
    call check(abs(indices_value(idx, "c") - 1.0_wp) < 1e-9_wp, "c = 1 after 100 yr", nfail)
    call check(abs(indices_value(idx, "dto")) < 1e-12_wp, "dto untouched by update", nfail)
    call indices_set(idx, "dto", -2.5_wp)
    call check(abs(indices_value(idx, "dto") + 2.5_wp) < 1e-12_wp, "dto set", nfail)
    call indices_update(idx, 150.0_dp, var=0.0_wp)
    call check(abs(indices_value(idx, "dto") + 2.5_wp) < 1e-12_wp, "dto kept across an update", nfail)
    call check(.not. indices_kill(idx), "no kill", nfail)

    ! Forcing-increment restart: c moved by 1.5 > 0.5 since the start.
    call check(indices_restart_due(idx), "restart due after |df| > 0.5", nfail)
    call check(.not. indices_restart_due(idx), "not due again without motion", nfail)

    ! Restart files: c has one, dto has none; the read restores c.
    call indices_restart_write(idx, tmp, 150.0_dp)
    inquire(file=tmp//"/idx_c_restart.nc", exist=found)
    call check(found, "idx_c_restart.nc written", nfail)
    inquire(file=tmp//"/idx_dto_restart.nc", exist=found)
    call check(.not. found, "no restart file for the derived index", nfail)

    call indices_init(idx, par, 0.0_dp)
    call indices_restart_read(idx, tmp)
    call check(abs(indices_value(idx, "c") - 1.5_wp) < 1e-9_wp, "c restored from its restart file", nfail)
    call check(abs(indices_value(idx, "dto")) < 1e-12_wp, "dto starts at 0 after a restart", nfail)

    ! 1D output: both values, df_dt for c only.
    call nc_create(tmp//"/ts.nc")
    call nc_write_dim(tmp//"/ts.nc", "time", x=0.0_dp, dx=1.0_dp, nx=1, unlimited=.true.)
    call indices_write_step(idx, tmp//"/ts.nc", 150.0_dp)
    call check(nc_has(tmp//"/ts.nc", "idx_c") .and. nc_has(tmp//"/ts.nc", "idx_dto") .and. &
               nc_has(tmp//"/ts.nc", "idx_c_df_dt") .and. .not. nc_has(tmp//"/ts.nc", "idx_dto_df_dt"), &
               "1D output variables", nfail)

    call execute_command_line("rm -rf "//tmp)

    write(*,*)
    if (nfail == 0) then
        write(*,*) "test_indices: PASS"
    else
        write(*,*) "test_indices: FAIL", nfail, "checks"
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

end program test_indices
