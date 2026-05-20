program yelmoxmg

    use nml 
    use ncio 
    use timestepping
    use timer
    use timeout 
    use yelmo 
    use ice_optimization
    use ice_sub_regions

    ! External libraries
    use fastisostasy    ! also reexports barysealevel
    use snapclim
    use marine_shelf
    use smbpal
    use sediments
    use geothermal
    
    implicit none 

    type(tstep_class)      :: ts
    
    type(yelmo_class)      :: yelmo1
    type(bsl_class)        :: bsl
    type(snapclim_class)   :: snp1
    type(marshelf_class)   :: mshlf1
    type(smbpal_class)     :: smbpal1
    type(sediments_class)  :: sed1
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    
contains

    subroutine yelmox_restart_write(bsl,isos,ylmo,mshlf,time,fldr)

        implicit none

        type(bsl_class),      intent(IN) :: bsl
        type(isos_class),     intent(IN) :: isos
        type(yelmo_class),    intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        real(wp),             intent(IN) :: time 
        character(len=*),     intent(IN), optional :: fldr
        
        ! Local variables
        real(wp) :: time_kyr
        character(len=32)   :: time_str
        character(len=1024) :: outfldr

        character(len=56), parameter :: file_bsl   = "bsl_restart.nc"
        character(len=56), parameter :: file_isos  = "isos_restart.nc"
        character(len=56), parameter :: file_yelmo = "yelmo_restart.nc"
        character(len=56), parameter :: file_mshlf = "marine_shelf.nc"

        if (present(fldr)) then
            outfldr = trim(fldr)
        else
            time_kyr = time*1e-3
            write(time_str,"(f20.3)") time_kyr
            outfldr = "./"//"restart-"//trim(adjustl(time_str))//"-kyr"
        end if

        write(*,*) "yelmox_restart_write:: outfldr = ", trim(outfldr)

        ! Make directory (use -p to ignore if directory already exists)
        call execute_command_line('mkdir -p "' // trim(outfldr) // '"')
        
        call bsl_restart_write(bsl,trim(outfldr)//"/"//file_bsl,time)
        call isos_restart_write(isos,trim(outfldr)//"/"//file_isos,time)
        call yelmo_restart_write(ylmo,trim(outfldr)//"/"//file_yelmo,time) 
        call marshelf_restart_write(mshlf,trim(outfldr)//"/"//file_mshlf,time)

        return

    end subroutine yelmox_restart_write

end program yelmoxmg