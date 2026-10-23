module ocean2ice

    use ncio 

    implicit none
    
    private
    public :: obm_temperature_kelvin

    integer,  parameter :: wp = kind(1.0)


contains

function obm_temperature_kelvin(obm_temp, name) result(t_K)
    ! A box-model ocean temperature in K, the unit of the ice-sheet ocean
    ! forcing: nautilus works in degC.
    implicit none
    real(wp),         intent(in) :: obm_temp
    character(len=*), intent(in) :: name
    real(wp) :: t_K
    real(wp), parameter :: deg2kelvin = 273.15

    select case(trim(name))
        case("nautilus")
            t_K = obm_temp + deg2kelvin
        case DEFAULT
            write(*,*) "obm_temperature_kelvin:: error: unknown ocean box model '"//trim(name)//"'."
            error stop 1
    end select

end function obm_temperature_kelvin

end module ocean2ice