module obm

    use ncio 
    use nml 
    use obm_defs
    use nautilus
    use stommel    

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp
    integer,  parameter :: wp   = sp 

    real(wp), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(wp), parameter :: pi        = 3.14159265359

    private
    public :: obm_init
    public :: obm_update
    public :: write_obm_init
    public :: write_obm_update
    public :: write_obm_restart

contains

    subroutine obm_init(obm_object, path_par, name)
        implicit none 
        type(obm_class) :: obm_object 
        character(len=512):: path_par, name

        select case(name)
            case("stommel")
                call stommel_init(obm_object, path_par, name)
            case("nautilus")
                call nautilus_init(obm_object, path_par, name)
            case DEFAULT
                ! do nothing
        end select

        return

    end subroutine obm_init

    subroutine obm_update(obm_object, dt, name)
        implicit none 
        type(obm_class) :: obm_object 
        character(len=512):: name
        real(wp) :: dt

        select case(name)
            case("stommel")
                call stommel_update(obm_object, dt)
            case("nautilus")
                call nautilus_update(obm_object, dt)
            case DEFAULT
                ! do nothing
        end select

        return

    end subroutine obm_update

    subroutine write_obm_init(filename,time_init, units)

        implicit none 
        character(len=*),  intent(IN)   :: filename, units
        real(wp),          intent(IN)   :: time_init

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"time", x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)
        
        return 

    end subroutine write_obm_init

    subroutine write_obm_update(obm_object, filename, name, t)

        implicit none 
        type(obm_class),    intent(IN) :: obm_object 
        character(len=*),     intent(IN) :: filename, name
        real(wp),             intent(IN) :: t
        real(wp) :: factor

        ! Local variables
        integer    :: ncid, n
        real(wp) :: t_prev 

        select case(name)
                case("stommel")
                        factor = 1
                case("nautilus")
                        factor = 365*24*60*60 * 1e6
                case DEFAULT
                        factor = 1
        end select

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",t_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(t-t_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",t,dim1="time",start=[n],count=[1],ncid=ncid)

        ! ===== obm variables =====

        call nc_write(filename,"tn",obm_object%tn,units="degC",long_name="North Atlantic box temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tt",obm_object%tt,units="degC",long_name="Tropical Atlantic temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ttd",obm_object%ttd,units="degC",long_name="Deep tropical Atlantic temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ts",obm_object%ts,units="degC",long_name="South Atlantic box temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"sn",obm_object%sn,units="psu",long_name="North Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"st",obm_object%st,units="psu",long_name="Tropical Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"std",obm_object%std,units="psu",long_name="Deep tropical Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ss",obm_object%ss,units="psu",long_name="South Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"m",obm_object%m / factor,units="Sv",long_name="Thermo Haline Circulation", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fn",obm_object%fn,units="Sv",long_name="(external) North Atlantic box fresh water flux", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ft",obm_object%ft,units="Sv",long_name="(external) Tropical Atlantic box fresh water flux", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fs",obm_object%fs,units="Sv",long_name="(external) South Atlantic box fresh water flux", &
                dim1="time",start=[n],ncid=ncid)  
        call nc_write(filename,"phin",obm_object%phin, units="Sv",long_name="Freshwater transport from tropical to northern box", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"phit",obm_object%phit, units="Sv",long_name="Freshwater transport from southern to tropical box", &
                dim1="time",start=[n],ncid=ncid) 
        call nc_write(filename,"thetan",obm_object%thetan,units="K",long_name="North Atlantic box air temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"thetat",obm_object%thetat,units="K",long_name="Tropical Atlantic box air temperature", &
                dim1="time",start=[n],ncid=ncid) 
        call nc_write(filename,"thetas",obm_object%thetas,units="K",long_name="South Atlantic box air temperature", &
                dim1="time",start=[n],ncid=ncid) 

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_obm_update

    subroutine write_obm_restart(obm_object, filename, t, units, init)

        implicit none 

        type(obm_class), intent(IN) :: obm_object
        character(len=*),  intent(IN) :: filename, units 
        real(wp),          intent(IN) :: t 
        logical, optional, intent(IN) :: init 

        ! Local variables
        integer  :: ncid, n, nx, ny, nz, nz_ac, nz_r   
        logical  :: initialize_file 
        real(wp) :: t_prev         

        initialize_file = .TRUE. 
        if (present(init)) initialize_file = init

        ! Write all stommel data to file, so that it can be
        ! read later to restart a simulation.
        
        ! == Initialize netcdf file ==============================================

        if (initialize_file) then
            call nc_create(filename)
            call nc_write_dim(filename,"time", x=t,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)
        end if 

        ! == Begin writing data ==============================================
        
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)
        
        if (initialize_file) then 

            ! Current time index to write will be the first and only one 
            n = 1 

        else

            ! Determine current writing time step 
            n = nc_size(filename,"time",ncid)
            call nc_read(filename,"time",t_prev,start=[n],count=[1],ncid=ncid) 
            if (abs(t-t_prev).gt.1e-5) n = n+1 

            ! Update the time step
            call nc_write(filename,"time",t,dim1="time",start=[n],count=[1],ncid=ncid)

        end if 

        call nc_write(filename,"tn",obm_object%tn,units="degC",long_name="North Atlantic box temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tt",obm_object%tt,units="degC",long_name="Tropical Atlantic temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ttd",obm_object%ttd,units="degC",long_name="Deep tropical Atlantic temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ts",obm_object%ts,units="degC",long_name="South Atlantic box temperature", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"sn",obm_object%sn,units="psu",long_name="North Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"st",obm_object%st,units="psu",long_name="Tropical Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"std",obm_object%std,units="psu",long_name="Deep tropical Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ss",obm_object%ss,units="psu",long_name="South Atlantic box salinity", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"m",obm_object%m,units="Sv",long_name="Thermo Haline Circulation", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fn",obm_object%fn,units="Sv",long_name="(external) North Atlantic box fresh water flux", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"ft",obm_object%ft,units="Sv",long_name="(external) Tropical Atlantic box fresh water flux", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fs",obm_object%fs,units="Sv",long_name="(external) South Atlantic box fresh water flux", &
                dim1="time",start=[n],ncid=ncid)  
        call nc_write(filename,"phin",obm_object%phin,units="Sv",long_name="Freshwater transport from tropical to northern box", &
                dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"phit",obm_object%phit,units="Sv",long_name="Freshwater transport from southern to tropical box", &
                dim1="time",start=[n],ncid=ncid) 
                
        ! Close the netcdf file
        call nc_close(ncid)

        ! Write summary 
        write(*,*) 
        write(*,*) "time = ", t, " : saved restart file: ", trim(filename)
        write(*,*) 

        return 

    end subroutine write_obm_restart


end module obm