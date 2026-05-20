! Aim: Run tests based on Zickfield et al. (2003) to test nautilus.f90
program run_nautilus
    ! use ncio
    ! use nml
    use obm_defs
    use obm
    use nautilus

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Parameters
    real(sp), parameter :: t0=-120.0e3, tend=0.0, dt=1

    ! Local variables
    integer :: n
    real(sp) :: t

    ! Initialize Nautilus object
    type(obm_class) :: nautilus
    call nautilus_init(nautilus, "/home/sergio/entra/proyects/d07_YelmoXBipolar/v1.12.2/par/obm_nautilus.nml","nautilus")

    do n = 0, ceiling((tend-t0)/dt)
        ! Get current time 
        t = t0 + n * dt
        print*, t
    end do

end program run_nautilus
