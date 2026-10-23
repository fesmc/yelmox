

module geothermal

    use nml 
    use ncio 
    use interp2D, only : fill_nearest
    use coordinates, only : grid_class
    use mapping,  only : map_class
    use fesmdata, only : fesmdata_map_init, fesmdata_read_field

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

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

contains 

    subroutine geothermal_init(gthrm,filename,grid,domain,group)
        ! Eventually interface should use 'sed' object directly
        ! instead of input/output arrays 

        implicit none 

        type(geothermal_class), intent(OUT) :: gthrm  
        character(len=*),       intent(IN)  :: filename 
        type(grid_class),       intent(IN)  :: grid         ! Yelmo grid
        character(len=*),       intent(IN)  :: domain
        character(len=*),       intent(IN), optional :: group

        ! Local variables
        character(len=32) :: nml_group
        integer :: nx, ny
        logical :: use_err

        real(wp), parameter :: ghf_min = 0.1  ! [mW/m2] Minimum allowed ghf value

        ! Load geothermal parameters
        call geothermal_par_load(gthrm%par,filename,domain,grid%name,init=.TRUE.,group=group)

        nx = grid%G%nx
        ny = grid%G%ny

        ! Allocate the geothermal object 
        call geothermal_allocate(gthrm%now,nx,ny)

        ! ====================================
        !
        ! Define the geothermal heat flow
        !
        ! ====================================
        
        if (gthrm%par%use_obs) then 
            ! Geothermal heat flux from file
            call read_ghf(gthrm%now%ghf,gthrm%par%obs_path,gthrm%par%obs_name,grid)

            ! Offset by f_stdev times its error field (e.g. standard deviation),
            ! read only when it is used
            use_err = gthrm%par%f_stdev .ne. 0.0_wp .and. &
                      trim(gthrm%par%obs_err_name) .ne. "None" .and. &
                      trim(gthrm%par%obs_err_name) .ne. "none" .and. &
                      trim(gthrm%par%obs_err_name) .ne. ""
            if (use_err) then 
                call read_ghf(gthrm%now%ghf_err,gthrm%par%obs_path,gthrm%par%obs_err_name,grid)
                gthrm%now%ghf = gthrm%now%ghf + gthrm%par%f_stdev*gthrm%now%ghf_err
            end if

            ! Make sure all values stay above set minimum
            where (gthrm%now%ghf .lt. ghf_min) gthrm%now%ghf = ghf_min
            
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

    subroutine read_ghf(var,path,varname,grid)
        ! A heat flux field (ghf or its error) of a file on grid, in mW/m2:
        ! remapped (bilinear) from the grid of the file (FesmData global
        ! attribute grid_name) if it is another grid, else read as it is;
        ! missing values filled; converted to mW/m2.

        implicit none

        real(wp),         intent(OUT) :: var(:,:)
        character(len=*), intent(IN)  :: path, varname
        type(grid_class), intent(IN)  :: grid

        type(map_class) :: map
        logical :: remap
        integer :: nx, ny
        real(wp), allocatable :: tmp(:,:)

        real(wp), parameter :: ghf_mv = -9999.0_wp   ! missing value marker

        call fesmdata_map_init(map,remap,nx,ny,path,varname,"bilin",grid=grid)
        call fesmdata_read_field(path,varname,tmp,nx,ny,remap,map)
        var = tmp

        if (remap) write(*,*) "geothermal_init:: "//trim(varname)//" remapped onto "//trim(grid%name)
        write(*,*) "geothermal_init:: "//trim(varname)//" loaded from: "//trim(path)

        ! A heat flux (and its error) is never negative: negative values are
        ! the file's missing values (e.g. -9999, or -9e33 over the ocean in
        ! the ISMIP7 fields). Fill them from the nearest valid value, so that
        ! cells that are (or become) ice covered get a real heat flux.
        if (any(var .lt. 0.0_wp)) then
            where (var .lt. 0.0_wp) var = ghf_mv
            call fill_nearest(var,ghf_mv)
            write(*,*) "geothermal_init:: "//trim(varname)//": missing values filled from the nearest valid value."
        end if

        ! The model needs mW/m2; the units attribute of the variable says
        ! what the file holds (e.g. the ISMIP7 fields are in W/m2).
        var = var * ghf_to_mW(path,varname)

        return

    end subroutine read_ghf


    function ghf_to_mW(filename,varname) result(fac)
        ! Factor converting the variable varname of filename to mW/m2,
        ! from its units attribute (W/m2 or mW/m2, in any of the usual
        ! spellings). Other units stop the model.

        implicit none

        character(len=*), intent(IN) :: filename, varname
        real(wp) :: fac

        ! Local variables
        character(len=256) :: units, u
        integer :: i

        units = ""
        call nc_read_attr(filename,varname,"units",units)

        ! Normalise the spelling: "mW m**-2", "mW m-2", "mW/m^2" -> "mWm-2", "mW/m2"
        u = ""
        do i = 1, len_trim(units)
            if (index(" *^",units(i:i)) .eq. 0) u = trim(u)//units(i:i)
        end do

        select case(trim(u))
            case("mWm-2","mW/m2")
                fac = 1.0_wp
            case("Wm-2","W/m2")
                fac = 1000.0_wp
                write(*,*) "geothermal_init:: "//trim(varname)//" converted from W/m2 to mW/m2."
            case DEFAULT
                write(*,*) "geothermal_init:: Error: units of the geothermal heat flux not recognized."
                write(*,*) "file:     ", trim(filename)
                write(*,*) "variable: ", trim(varname)
                write(*,*) "units:    ", trim(units)
                write(*,*) "Expected W/m2 or mW/m2 (e.g. 'mW m-2', 'mW m**-2', 'W/m^2')."
                stop
        end select

        return

    end function ghf_to_mW


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
        if (present(init)) init_pars = init

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


