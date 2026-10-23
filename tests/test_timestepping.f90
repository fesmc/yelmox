program test_timestepping
    ! Checks the fesm-utils timestepper as yelmox drives it (tstep_init, then a
    ! loop of output, exit when finished, tstep_update): for the "rel", "cal" and
    ! "const" methods the run ends on time_end after the expected number of steps,
    ! with the calendar and relative times consistent with time_ref and the
    ! constants. The last case steps 0.1 yr from -130 kyr, which ends exactly on
    ! time_end only with the double-precision time.
    !
    ! Build and run from the yelmox root (serial fesm-utils):
    !   gfortran -ffree-line-length-none -Ifesm-utils/include-serial tests/test_timestepping.f90 \
    !       -Lfesm-utils/include-serial -lfesmutils $(nf-config --flibs) -L$(nc-config --libdir) \
    !       -o test_timestepping.x
    !   ./test_timestepping.x

    use precision,    only : dp
    use timestepping, only : tstep_class, tstep_init, tstep_update

    implicit none

    integer :: nfail

    nfail = 0

    !                  method   time_init   time_end   dtt     n      time_rel   time_cal
    call check_run("rel",    -21000.0_dp,  1000.0_dp, 10.0_dp,  2200,  1000.0_dp, 2950.0_dp, nfail)
    call check_run("cal",      1880.0_dp,  2100.0_dp, 10.0_dp,    22,   150.0_dp, 2100.0_dp, nfail)
    call check_run("const",       0.0_dp,   100e3_dp, 10.0_dp, 10000, -21000.0_dp, 1880.0_dp, nfail, &
                   const_rel=-21000.0_dp, const_cal=1880.0_dp)
    call check_run("rel",   -130000.0_dp, -129000.0_dp, 0.1_dp, 10000, -129000.0_dp, -127050.0_dp, nfail)

    write(*,*)
    if (nfail == 0) then
        write(*,*) "test_timestepping: PASS"
    else
        write(*,*) "test_timestepping: FAIL", nfail, "checks"
        stop 1
    end if

contains

    subroutine check_run(method, time_init, time_end, dtt, n_exp, rel_exp, cal_exp, nfail, &
                         const_rel, const_cal)
        ! Run the yelmox time loop for one method and check the final state.
        character(len=*), intent(in)    :: method
        real(dp),         intent(in)    :: time_init, time_end, dtt
        integer,          intent(in)    :: n_exp
        real(dp),         intent(in)    :: rel_exp, cal_exp
        integer,          intent(inout) :: nfail
        real(dp),         intent(in), optional :: const_rel, const_cal

        real(dp), parameter :: tol = 1e-9_dp
        type(tstep_class) :: ts
        character(len=96) :: label

        call tstep_init(ts, time_init, time_end, method=method, units="yr", time_ref=1950.0_dp, &
                        const_rel=const_rel, const_cal=const_cal)

        do
            if (ts%is_finished) exit
            call tstep_update(ts, dtt)
        end do

        write(label, "(a,a,f0.1,a,f0.1,a,f0.2)") trim(method), " ", time_init, " -> ", time_end, " dt ", dtt
        call check(trim(label)//": steps",    real(ts%n,dp), real(n_exp,dp), 0.0_dp, nfail)
        call check(trim(label)//": time",     ts%time,       time_end,        tol,    nfail)
        call check(trim(label)//": time_rel", ts%time_rel,   rel_exp,         tol,    nfail)
        call check(trim(label)//": time_cal", ts%time_cal,   cal_exp,         tol,    nfail)

    end subroutine check_run

    subroutine check(label, got, expected, tol, nfail)
        character(len=*), intent(in)    :: label
        real(dp),         intent(in)    :: got, expected, tol
        integer,          intent(inout) :: nfail

        if (abs(got - expected) <= tol) then
            write(*,"(a,a,a,g0)") "  PASS ", label, " = ", got
        else
            write(*,"(a,a,a,g0,a,g0)") "  FAIL ", label, " = ", got, ", expected ", expected
            nfail = nfail + 1
        end if

    end subroutine check

end program test_timestepping
