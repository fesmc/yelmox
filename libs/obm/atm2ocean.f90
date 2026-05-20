module atm2ocean

    use ncio 

    implicit none
    
    private
    public :: update_theta_from_snapclim

    integer,  parameter :: wp = kind(1.0)


contains

subroutine update_theta_from_snapclim(obm_theta, snapclim_theta, mask, hemisphere, name)
    ! Input
    real(wp) :: obm_theta   ! ºC
    real(wp), allocatable :: snapclim_theta(:,:) ! K
    real(wp), allocatable :: mask(:,:)
    character(len=*) :: hemisphere, name

    ! Local variables
    logical, allocatable :: mask_box(:,:)
    integer :: npts_box, nx, ny
    logical, parameter :: use_extended_north=.TRUE.
    real(wp), parameter :: deg2kelvin = 273.15
    
    ! Allocate box mask
    nx = size(mask,1)
    ny = size(mask,2)

    allocate(mask_box(nx,ny))

    select case(hemisphere)
        ! Mask values:
        !   1 = atlantic ocean
        !   2 = arctic ocean
        !   3 = subarctic ocean
        !   4 = indian ocean
        !   5 = pacific ocean
        !   6 = southern ocean

        case("north")

            if (use_extended_north) then
                mask_box = (mask .eq. 1 .or. mask .eq. 3)   ! TThis should be made in the initialization process since it does not change
            else
                mask_box  = (mask .eq. 1)
            end if

            npts_box  = real(count(mask_box),wp)
        
            obm_theta = sum(snapclim_theta, mask=mask_box) / npts_box

        case("south")
            mask_box  = (mask .eq. 1)
            npts_box  = real(count(mask_box),wp)

            obm_theta = sum(snapclim_theta, mask=mask_box) / npts_box

        case DEFAULT
            mask_box  = (mask .gt. 0)
            npts_box  = real(count(mask_box),wp)

            obm_theta = sum(snapclim_theta, mask=mask_box) / npts_box
    end select

    select case(name)
        case("stommel")
            obm_theta = obm_theta
        case("nautilus")
            obm_theta = obm_theta - deg2kelvin
        case DEFAULT
            ! do nothing
    end select

    return

end subroutine update_theta_from_snapclim


end module atm2ocean