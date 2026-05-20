

module geothermal

    use nml 
    use ncio 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    real(wp), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(wp), parameter :: pi        = 3.14159265359

    type geothermal_param_class
        logical            :: use_obs
        character(len=512) :: obs_path
        character(len=56)  :: obs_name 
        character(len=56)  :: obs_err_name
        real(wp)           :: f_stdev
        real(wp)           :: ghf_const 
    end type 

    type geothermal_state_class 
        real(wp), allocatable :: ghf(:,:)
        real(wp), allocatable :: ghf_err(:,:)
    end type 

    type geothermal_class
        type(geothermal_param_class) :: par 
        type(geothermal_state_class) :: now 
    end type

    private
    public :: geothermal_class
    public :: geothermal_init 
    public :: geothermal_update
    public :: geothermal_end 

contains 

    subroutine geothermal_init(gthrm,filename,nx,ny,domain,grid_name,group)
        ! Eventually interface should use 'sed' object directly
        ! instead of input/output arrays 

        implicit none 

        type(geothermal_class), intent(OUT) :: gthrm  
        character(len=*),       intent(IN)  :: filename 
        integer,                intent(IN)  :: nx, ny 
        character(len=*),       intent(IN)  :: domain, grid_name
        character(len=*),       intent(IN), optional :: group

        ! Local variables
        character(len=32) :: nml_group
        
        real(wp), parameter :: ghf_min = 0.1  ! [mW/m2] Minimum allowed ghf value

        ! Load geothermal parameters
        call geothermal_par_load(gthrm%par,filename,domain,grid_name,init=.TRUE.,group=group)

        ! Allocate the geothermal object 
        call geothermal_allocate(gthrm%now,nx,ny)

        ! ====================================
        !
        ! Define the geothermal heat flow
        !
        ! ====================================
        
        ! Read in array 
        if (gthrm%par%use_obs) then 
            ! Load geothermal heat flow data from file 
            call nc_read(gthrm%par%obs_path,gthrm%par%obs_name,gthrm%now%ghf)
            write(*,*) "geothermal_init:: geothermal heat flux loaded from: "
            write(*,*) trim(gthrm%par%obs_path)//" : "//trim(gthrm%par%obs_name)

            ! Also read it ghf error (e.g. standard deviation) if available

            if ( (.not.  trim(gthrm%par%obs_err_name) .eq. "None") .and. &
                 (.not.  trim(gthrm%par%obs_err_name) .eq. "none") .and. &
                 (.not.  trim(gthrm%par%obs_err_name) .eq. "")   ) then 
                ! Load stdev field too

                call nc_read(gthrm%par%obs_path,gthrm%par%obs_err_name,gthrm%now%ghf_err)

                ! Offset GHF field by desired sigma level
                gthrm%now%ghf = gthrm%now%ghf + gthrm%par%f_stdev*gthrm%now%ghf_err

                ! Make sure all values stay above set minimum
                where (gthrm%now%ghf .lt. ghf_min) gthrm%now%ghf = ghf_min

            end if
            
        else 
            ! Set geothermal heat flux to constant value
            gthrm%now%ghf = gthrm%par%ghf_const 

        end if 

        ! ====================================
        !
        ! Additional calculations
        !
        ! ====================================

        ! Move this conversion to inside of yelmo_thermodynamics
!         ! Convert units for internal use
!         gthrm%now%ghf = -sec_year/1000.0*gthrm%now%ghf    ! [mW/m2] => [J/m2/a]   

        ! ====================================
        !
        ! Summary
        !
        ! ====================================

        write(*,*) "geothermal_init:: range ghf: ", minval(gthrm%now%ghf), maxval(gthrm%now%ghf)

        return 

    end subroutine geothermal_init


    subroutine geothermal_update(sed,year_bp)

        implicit none 

        type(geothermal_class) :: sed 

        real(wp) :: year_bp 

        ! === To do ====


        return 

    end subroutine geothermal_update


    subroutine geothermal_end(sed)

        implicit none 

        type(geothermal_class) :: sed

        ! Deallocate geothermal state object
        call geothermal_deallocate(sed%now)

        return 

    end subroutine geothermal_end


    subroutine geothermal_par_load(par,filename,domain,grid_name,init,group)

        type(geothermal_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain, grid_name 
        logical, optional :: init 
        logical :: init_pars 
        character(len=*),  intent(IN), optional :: group

        ! Local variables
        character(len=32) :: nml_group
        
        ! Make sure we know the namelist group for the geothermal block
        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "ghf"               ! Default parameter block name
        end if

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,nml_group,"use_obs",         par%use_obs,        init=init_pars)
        call nml_read(filename,nml_group,"obs_path",        par%obs_path,       init=init_pars)
        call nml_read(filename,nml_group,"obs_name",        par%obs_name,       init=init_pars)
        call nml_read(filename,nml_group,"obs_err_name",    par%obs_err_name,   init=init_pars)
        call nml_read(filename,nml_group,"f_stdev",         par%f_stdev,        init=init_pars)
        call nml_read(filename,nml_group,"ghf_const",       par%ghf_const,      init=init_pars)

        ! Replace gridding template values from path
        call parse_path(par%obs_path,domain,grid_name)
        
        return

    end subroutine geothermal_par_load


    ! =======================================================
    !
    ! geothermal memory management
    !
    ! =======================================================

    subroutine geothermal_allocate(now,nx,ny)

        implicit none 

        type(geothermal_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call geothermal_deallocate(now)

        ! Allocate geothermal 
        allocate(now%ghf(nx,ny))
        allocate(now%ghf_err(nx,ny))
        
        ! Initialize to zero
        now%ghf     = 0.0
        now%ghf_err = 0.0
        
        return

    end subroutine geothermal_allocate

    subroutine geothermal_deallocate(now)

        implicit none 

        type(geothermal_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%ghf))     deallocate(now%ghf)
        if (allocated(now%ghf_err)) deallocate(now%ghf_err)

        return

    end subroutine geothermal_deallocate

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    

end module geothermal 


