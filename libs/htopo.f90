module htopo
    ! Hi-resolution geometry hub of a domain (multigrid yelmox).
    !
    ! The hub sits *above* all physics modules (including Yelmo): its grid
    ! (grid_hub) is the finest resolution in the setup, and it is the reference
    ! geometry that the coupler remaps *from* when a coarser module needs
    ! z_bed/H_ice/z_srf/masks.
    !
    ! Two provenance classes of field live here:
    !   * regions, basins, z_bed_sd -- static (code masks and bed roughness),
    !     loaded once from file;
    !   * z_bed, H_ice, z_srf -- dynamic geometry, loaded here as the initial
    !     reference and refreshed each step from the models (refresh_hub).
    !
    ! The file paths and variable names come from the domain definition
    ! (&domain); {domain}/{grid_name} in the paths resolve to the domain name
    ! and the hub grid.

    use nml
    use ncio
    use coords, only : grid_class, grid_cdo_read_desc

    implicit none
    private

    integer, parameter :: wp = kind(1.0)     ! single precision (matches yelmox libs)

    type htopo_par_class
        character(len=256) :: domain
        character(len=256) :: grid_name       ! hub grid (grid_hub), e.g. "ANT-16KM"
        character(len=512) :: topo_path
        character(len=56)  :: topo_names(4)   ! z_bed, H_ice, z_srf, z_bed_sd ("" = none)
        character(len=512) :: regions_path    ! "" = no file (regions = 1)
        character(len=56)  :: regions_var
        character(len=512) :: basins_path     ! "" = no file (basins = 1)
        character(len=56)  :: basins_var
    end type

    type htopo_class
        type(htopo_par_class) :: par
        type(grid_class)      :: grid         ! topo grid, from grid_<name>.txt
        integer               :: nx, ny
        real(wp), allocatable :: z_bed(:,:)   ! [m] bedrock elevation
        real(wp), allocatable :: H_ice(:,:)   ! [m] ice thickness
        real(wp), allocatable :: z_srf(:,:)   ! [m] surface elevation
        real(wp), allocatable :: z_bed_sd(:,:) ! [m] standard deviation of z_bed (static)
        real(wp), allocatable :: regions(:,:) ! region mask
        real(wp), allocatable :: basins(:,:)  ! basin mask
        ! Dynamic geometry refreshed from the models each step (not file-loaded).
        real(wp), allocatable :: f_grnd(:,:)  ! [1] grounded-ice fraction
        real(wp), allocatable :: z_sl(:,:)    ! [m] sea-surface / sea-level height
    end type

    public :: htopo_class, htopo_init
    public :: htopo_write_init, htopo_write_step

contains

    subroutine htopo_init(htopo, filename, group, domain, grid_name, map_fldr)
        ! Load the hub's file paths from the domain definition, resolve its grid
        ! from the disk grid table, and read the reference fields onto that grid.
        type(htopo_class), intent(out) :: htopo
        character(len=*),  intent(in)  :: filename   ! parameter file
        character(len=*),  intent(in)  :: group      ! namelist group, e.g. "domain"
        character(len=*),  intent(in)  :: domain     ! domain name
        character(len=*),  intent(in)  :: grid_name  ! hub grid (grid_hub)
        character(len=*),  intent(in), optional :: map_fldr

        character(len=256) :: mfldr

        mfldr = "maps"
        if (present(map_fldr)) mfldr = trim(map_fldr)

        call htopo_par_load(htopo%par, filename, group, domain, grid_name)

        ! Topo grid definition (nx,ny + coordinates) from grid_<name>.txt.
        call grid_cdo_read_desc(htopo%grid, trim(htopo%par%grid_name), trim(mfldr))
        htopo%nx = htopo%grid%G%nx
        htopo%ny = htopo%grid%G%ny

        allocate(htopo%z_bed(htopo%nx,htopo%ny))
        allocate(htopo%H_ice(htopo%nx,htopo%ny))
        allocate(htopo%z_srf(htopo%nx,htopo%ny))
        allocate(htopo%z_bed_sd(htopo%nx,htopo%ny))
        allocate(htopo%regions(htopo%nx,htopo%ny))
        allocate(htopo%basins(htopo%nx,htopo%ny))
        allocate(htopo%f_grnd(htopo%nx,htopo%ny)); htopo%f_grnd = 0.0_wp
        allocate(htopo%z_sl(htopo%nx,htopo%ny));   htopo%z_sl   = 0.0_wp

        call nc_read(htopo%par%topo_path,    htopo%par%topo_names(1), htopo%z_bed)
        call nc_read(htopo%par%topo_path,    htopo%par%topo_names(2), htopo%H_ice)
        call nc_read(htopo%par%topo_path,    htopo%par%topo_names(3), htopo%z_srf)

        htopo%z_bed_sd = 0.0_wp
        if (len_trim(htopo%par%topo_names(4)) > 0) &
            call nc_read(htopo%par%topo_path, htopo%par%topo_names(4), htopo%z_bed_sd)

        ! Static masks: load from file when a path is given, else default to a
        ! single region/basin (1.0), so paleo domains without mask files run.
        htopo%regions = 1.0_wp
        htopo%basins  = 1.0_wp
        if (len_trim(htopo%par%regions_path) > 0) &
            call nc_read(htopo%par%regions_path, htopo%par%regions_var, htopo%regions)
        if (len_trim(htopo%par%basins_path) > 0) &
            call nc_read(htopo%par%basins_path,  htopo%par%basins_var,  htopo%basins)

    end subroutine htopo_init

    subroutine htopo_par_load(par, filename, group, domain, grid_name)
        type(htopo_par_class), intent(out) :: par
        character(len=*),      intent(in)  :: filename, group
        character(len=*),      intent(in)  :: domain, grid_name

        par%domain    = trim(domain)
        par%grid_name = trim(grid_name)

        ! Blank entries read as "" (nml_read leaves the value untouched)
        par%topo_names   = ""
        par%regions_path = ""
        par%basins_path  = ""

        call nml_read(filename, group, "topo_path",    par%topo_path)
        call nml_read(filename, group, "topo_names",   par%topo_names)
        call nml_read(filename, group, "regions_path", par%regions_path)
        call nml_read(filename, group, "regions_var",  par%regions_var)
        call nml_read(filename, group, "basins_path",  par%basins_path)
        call nml_read(filename, group, "basins_var",   par%basins_var)

        ! Resolve {domain}/{grid_name} against the hub grid.
        call parse_path(par%topo_path,    par%domain, par%grid_name)
        call parse_path(par%basins_path,  par%domain, par%grid_name)
        call parse_path(par%regions_path, par%domain, par%grid_name)

    end subroutine htopo_par_load

    subroutine htopo_write_init(htopo, filename, time_init)
        ! Create a 2D output file on the topo grid, with the static masks.
        type(htopo_class), intent(in) :: htopo
        character(len=*),  intent(in) :: filename
        real(wp),          intent(in) :: time_init

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=htopo%grid%G%x, units="km")
        call nc_write_dim(filename, "yc", x=htopo%grid%G%y, units="km")
        call nc_write_dim(filename, "time", x=time_init, dx=1.0_wp, nx=1, &
                          units="year", unlimited=.TRUE.)

        call nc_write(filename, "regions", htopo%regions, dim1="xc", dim2="yc", &
                      start=[1,1], long_name="Region mask", units="")
        call nc_write(filename, "basins", htopo%basins, dim1="xc", dim2="yc", &
                      start=[1,1], long_name="Basin mask", units="")
    end subroutine htopo_write_init

    subroutine htopo_write_step(htopo, filename, time)
        ! Append the dynamic hi-res geometry at `time`.
        type(htopo_class), intent(in) :: htopo
        character(len=*),  intent(in) :: filename
        real(wp),          intent(in) :: time

        integer :: ncid, n

        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)

        call nc_write(filename, "z_bed", htopo%z_bed, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Bedrock elevation")
        call nc_write(filename, "H_ice", htopo%H_ice, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Ice thickness")
        call nc_write(filename, "z_srf", htopo%z_srf, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Surface elevation")
        call nc_write(filename, "f_grnd", htopo%f_grnd, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="1", &
                      long_name="Grounded-ice fraction")
        call nc_write(filename, "z_sl", htopo%z_sl, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Sea-surface height")
        call nc_close(ncid)
    end subroutine htopo_write_step

    subroutine parse_path(path, domain, grid_name)
        character(len=*), intent(inout) :: path
        character(len=*), intent(in)    :: domain, grid_name
        call nml_replace(path, "{domain}",    trim(domain))
        call nml_replace(path, "{grid_name}", trim(grid_name))
    end subroutine parse_path

end module htopo
