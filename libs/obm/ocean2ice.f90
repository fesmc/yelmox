module ocean2ice

    use ncio 

    implicit none
    
    private
    public :: calc_ocean_temperature_field
    public :: calc_ocean_salinity_field

    integer,  parameter :: wp = kind(1.0)


contains

subroutine calc_ocean_temperature_field(orig_otf, obm_temp, name)
    ! This subroutine computes a new ocean temperature field (otf)
    implicit none
    real(wp), allocatable :: orig_otf(:,:,:)
    real(wp) :: obm_temp
    character(len=512):: name
    real(wp), parameter :: deg2kelvin = 273.15

    select case(name)
        case("stommel")
            !orig_otf(:,:,:) = orig_otf(:,:,:) + obm_temp
            orig_otf(:,:,:) = obm_temp
        case("nautilus")
            !orig_otf(:,:,:) = orig_otf(:,:,:) + (obm_temp + deg2kelvin)
            orig_otf(:,:,:) = obm_temp + deg2kelvin
        case DEFAULT
            orig_otf(:,:,:) = -9999
    end select

    return 

end subroutine calc_ocean_temperature_field

subroutine calc_ocean_salinity_field(orig_osf, obm_salt)
    ! This subroutine computes a new ocean salinity field (osf)
    implicit none
    real(wp), allocatable :: orig_osf(:,:,:)
    real(wp) :: obm_salt

    orig_osf(:,:,:) = obm_salt

    return

end subroutine calc_ocean_salinity_field

end module ocean2ice