program test_timestepping

    use nml
    use ncio
    use timestepping
    use yelmo_defs, only : wp
    

    type(tstep_class) :: ts
    real(wp) :: time_init
    real(wp) :: time_end
    real(wp) :: dtt

    select case("cal")
        case("rel","sp","bp")
            time_init = -21.0   ! kyr st
            time_end  =  1.0    ! kyr st
            dtt       =  1.0    ! kyr
            call tstep_init(ts,time_init,time_end,method="sp",time_ref=1950e-3_wp,units="kyr")
        case("cal")
            time_init = 1880.0  ! yr CE
            time_end  = 2100.0  ! yr CE
            dtt       = 10.0    ! yr
            call tstep_init(ts,time_init,time_end,method="cal",time_ref=1950.0_wp,units="yr")
        case("const")
            time_init = 0.0     ! kyr
            time_end  = 100.0   ! kyr
            dtt       = 10.0
            call tstep_init(ts,time_init,time_end,method="const",units="yr",const_rel=-21e3_wp,const_cal=1880._wp)
    end select

    call tstep_print_header(ts)

    ! Advance timesteps
    do while (.not. ts%is_finished)

        ! Update timestep
        call tstep_update(ts,dtt)
        call tstep_print(ts)

        ! Perform model updates here...

    end do

end program
