module ice_sub_regions

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use ncio 

    implicit none 

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 


    private
    public :: get_ice_sub_region

contains

    subroutine get_ice_sub_region(mask,reg_name,domain,grid_name)

        implicit none

        logical, allocatable, intent(INOUT) :: mask(:,:)
        character(len=*),     intent(IN)    :: reg_name
        character(len=*),     intent(IN)    :: domain
        character(len=*),     intent(IN)    :: grid_name
        
        ! Local variables
        integer :: nx, ny
        character(len=1024) :: regions_mask_fnm
        real(wp), allocatable :: regions_mask(:,:)

        if (.not. allocated(mask)) then
            write(error_unit,*) "get_ice_sub_region:: Error: mask must be allocated before calling this routine."
            stop
        end if

        nx = size(mask,1)
        ny = size(mask,2)

        allocate(regions_mask(nx,ny))

        ! Define specific regions of interest =====================

        select case(trim(domain))

            case("Antarctica")

                ! Load general regions mask from file 
                regions_mask_fnm = "ice_data/Antarctica/"//trim(grid_name)//&
                                    "/"//trim(grid_name)//"_BASINS-nasa.nc"
                call nc_read(regions_mask_fnm,"mask_regions",regions_mask)

                select case(trim(reg_name))

                    case("APIS")    ! (region=3.0 in regions map)
                        mask = .FALSE. 
                        where(abs(regions_mask - 3.0) .lt. 1e-3) mask = .TRUE.
                    
                    case("WAIS")    ! (region=1.0 in regions map)
                        mask = .FALSE. 
                        where(abs(regions_mask - 1.0) .lt. 1e-3) mask = .TRUE.

                    case("EAIS")    ! (region=2.0 in regions map)
                        mask = .FALSE. 
                        where(abs(regions_mask - 2.0) .lt. 1e-3) mask = .TRUE.

                    case DEFAULT

                        write(error_unit,*) "get_ice_sub_region:: Error: region not yet defined."
                        write(error_unit,*) "reg_name = ", trim(reg_name)
                        stop

                end select
            
            case("Greenland")

                ! Load general regions mask from file 
                regions_mask_fnm = "ice_data"//"/"//trim(domain)//"/"//trim(grid_name)//&
                                                        "/"//trim(grid_name)//"_REGIONS.nc"
                call nc_read(regions_mask_fnm,"mask",regions_mask)

                select case(trim(reg_name))

                    case DEFAULT

                        write(error_unit,*) "get_ice_sub_region:: Error: region not yet defined."
                        write(error_unit,*) "reg_name = ", trim(reg_name)
                        stop

                end select

            case("Laurentide","North")

                ! Load general regions mask from file 
                regions_mask_fnm = "ice_data"//"/"//trim(domain)//"/"//trim(grid_name)//&
                                                        "/"//trim(grid_name)//"_REGIONS.nc"
                call nc_read(regions_mask_fnm,"mask",regions_mask)

                select case(trim(reg_name))

                    case("Hudson")  ! (region=1.12 in regions map) 
                        mask = .FALSE. 
                        where(abs(regions_mask - 1.12) .lt. 1e-3) mask = .TRUE.

                    case DEFAULT

                        write(error_unit,*) "get_ice_sub_region:: Error: region not yet defined."
                        write(error_unit,*) "reg_name = ", trim(reg_name)
                        stop

                end select

        end select

        return

    end subroutine get_ice_sub_region





end module ice_sub_regions