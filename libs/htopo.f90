module htopo
    ! Hi-resolution geometry hub of a domain (multigrid yelmox).
    !
    ! The hub sits *above* all physics modules (including Yelmo): its grid
    ! (grid_hub) is the finest resolution in the setup, and it is the reference
    ! geometry that the coupler remaps *from* when a coarser module needs
    ! z_bed/H_ice/z_srf/masks.
    !
    ! Three kinds of field live here:
    !   * regions, basins, sectors, z_bed_sd -- static (code masks and bed
    !     roughness), loaded once from file;
    !   * z_bed_ref, H_ice_ref, z_srf_ref -- the hi-res reference geometry,
    !     loaded once from file;
    !   * z_bed, H_ice, z_srf, f_grnd, z_sl -- the current geometry, refreshed
    !     each step from the models (couple_yelmo_to_htopo). On a hub finer
    !     than the ice sheet it is the reference plus the models' anomalies
    !     (htopo_update).
    !
    ! The file paths and variable names come from the domain definition
    ! (&domain); {domain}/{grid_name} in the paths resolve to the domain name
    ! and the hub grid. The domain definition also says where ice is allowed
    ! (ice_codes of `regions`), where it relaxes to the reference (relax_codes
    ! of `regions`), and names the regions of interest for 1D output
    ! (region_names/region_codes of one code mask).

    use nml
    use ncio
    use coords,   only : grid_class, grid_cdo_read_desc
    use interp2D, only : fill_nearest
    use phys_constants, only : phys_const_class, phys_const_require, phys_const_get

    implicit none
    private

    integer, parameter :: wp = kind(1.0)     ! single precision (matches yelmox libs)

    real(wp), parameter :: mv      = -9999.0_wp   ! missing value of the topography reads
    real(wp), parameter :: tol_code = 1e-3_wp     ! tolerance to match a mask code

    integer, parameter :: n_codes_max = 20        ! max entries of the code lists

    type htopo_par_class
        character(len=256) :: domain
        character(len=256) :: grid_name       ! hub grid (grid_hub), e.g. "ANT-16KM"
        character(len=512) :: topo_path
        character(len=56)  :: topo_names(4)   ! z_bed, H_ice, z_srf, z_bed_sd ("" = none)
        character(len=512) :: regions_path    ! "" = no file (regions = 1)
        character(len=56)  :: regions_var
        character(len=512) :: basins_path     ! "" = no file (basins = 1)
        character(len=56)  :: basins_var
        character(len=512) :: sectors_path    ! "" = no file (sectors = 1)
        character(len=56)  :: sectors_var
        character(len=16)  :: ice_codes_mode  ! where ice is allowed: all | include | exclude
        real(wp)           :: ice_codes(n_codes_max)      ! codes of `regions`
        character(len=56)  :: region_names(n_codes_max)   ! named regions ("" = none)
        character(len=16)  :: region_mask     ! code mask of the named regions: regions | basins | sectors
        real(wp)           :: region_codes(n_codes_max)   ! one code per named region
        character(len=16)  :: relax_codes_mode ! where ice relaxes to the reference: none | all | include | exclude
        real(wp)           :: relax_codes(n_codes_max)    ! codes of `regions`
        real(wp)           :: relax_tau       ! [yr] relaxation timescale there
        integer            :: n_ice_codes     ! number of ice_codes given
        integer            :: n_relax_codes   ! number of relax_codes given
        integer            :: n_regions       ! number of named regions
        real(wp)           :: rho_ice         ! [kg m-3] ice density (from the domain's constants)
        real(wp)           :: rho_sw          ! [kg m-3] seawater density
    end type

    type htopo_class
        type(htopo_par_class) :: par
        type(grid_class)      :: grid         ! topo grid, from grid_<name>.txt
        integer               :: nx, ny
        ! Reference geometry (static, from file).
        real(wp), allocatable :: z_bed_ref(:,:) ! [m] bedrock elevation
        real(wp), allocatable :: H_ice_ref(:,:) ! [m] ice thickness
        real(wp), allocatable :: z_srf_ref(:,:) ! [m] surface elevation
        real(wp), allocatable :: z_bed_sd(:,:)  ! [m] standard deviation of z_bed (static)
        real(wp), allocatable :: regions(:,:)   ! region mask
        real(wp), allocatable :: basins(:,:)    ! basin mask
        real(wp), allocatable :: sectors(:,:)   ! sector mask (e.g. Antarctic APIS/WAIS/EAIS)
        ! Current geometry, refreshed from the models each step.
        real(wp), allocatable :: z_bed(:,:)     ! [m] bedrock elevation
        real(wp), allocatable :: H_ice(:,:)     ! [m] ice thickness
        real(wp), allocatable :: z_srf(:,:)     ! [m] surface elevation
        real(wp), allocatable :: f_grnd(:,:)    ! [1] grounded-ice fraction
        real(wp), allocatable :: z_sl(:,:)      ! [m] sea-surface / sea-level height
    end type

    public :: htopo_class, htopo_init, htopo_update
    public :: htopo_ice_allowed, htopo_relax_tau, htopo_region_codes
    public :: htopo_write_init, htopo_write_step

contains

    subroutine htopo_init(htopo, filename, group, domain, grid_name, cnst, map_fldr)
        ! Load the hub's file paths from the domain definition, resolve its grid
        ! from the disk grid table, and read the reference fields onto that grid.
        type(htopo_class), intent(out) :: htopo
        character(len=*),  intent(in)  :: filename   ! parameter file
        character(len=*),  intent(in)  :: group      ! namelist group, e.g. "domain"
        character(len=*),  intent(in)  :: domain     ! domain name
        character(len=*),  intent(in)  :: grid_name  ! hub grid (grid_hub)
        type(phys_const_class), intent(in) :: cnst   ! physical constants of the domain
        character(len=*),  intent(in), optional :: map_fldr

        character(len=256) :: mfldr

        mfldr = "maps"
        if (present(map_fldr)) mfldr = trim(map_fldr)

        call htopo_par_load(htopo%par, filename, group, domain, grid_name)

        call phys_const_require(cnst, "htopo_init")
        call phys_const_get(cnst, "rho_ice", htopo%par%rho_ice)
        call phys_const_get(cnst, "rho_sw",  htopo%par%rho_sw)

        ! Topo grid definition (nx,ny + coordinates) from grid_<name>.txt.
        call grid_cdo_read_desc(htopo%grid, trim(htopo%par%grid_name), trim(mfldr))
        htopo%nx = htopo%grid%G%nx
        htopo%ny = htopo%grid%G%ny

        allocate(htopo%z_bed_ref(htopo%nx,htopo%ny))
        allocate(htopo%H_ice_ref(htopo%nx,htopo%ny))
        allocate(htopo%z_srf_ref(htopo%nx,htopo%ny))
        allocate(htopo%z_bed_sd(htopo%nx,htopo%ny))
        allocate(htopo%regions(htopo%nx,htopo%ny))
        allocate(htopo%basins(htopo%nx,htopo%ny))
        allocate(htopo%sectors(htopo%nx,htopo%ny))
        allocate(htopo%f_grnd(htopo%nx,htopo%ny)); htopo%f_grnd = 0.0_wp
        allocate(htopo%z_sl(htopo%nx,htopo%ny));   htopo%z_sl   = 0.0_wp

        call nc_read(htopo%par%topo_path,    htopo%par%topo_names(1), htopo%z_bed_ref, missing_value=mv)
        call nc_read(htopo%par%topo_path,    htopo%par%topo_names(2), htopo%H_ice_ref, missing_value=mv)
        call nc_read(htopo%par%topo_path,    htopo%par%topo_names(3), htopo%z_srf_ref, missing_value=mv)
        call htopo_fill_missing(htopo)

        ! The current geometry starts from the reference.
        htopo%z_bed = htopo%z_bed_ref
        htopo%H_ice = htopo%H_ice_ref
        htopo%z_srf = htopo%z_srf_ref

        htopo%z_bed_sd = 0.0_wp
        if (len_trim(htopo%par%topo_names(4)) > 0) then
            call nc_read(htopo%par%topo_path, htopo%par%topo_names(4), htopo%z_bed_sd, missing_value=mv)
            where (htopo%z_bed_sd == mv) htopo%z_bed_sd = 0.0_wp
        end if

        ! Static masks: load from file when a path is given, else default to a
        ! single region/basin/sector (1.0), so paleo domains without mask files run.
        htopo%regions = 1.0_wp
        htopo%basins  = 1.0_wp
        htopo%sectors = 1.0_wp
        if (len_trim(htopo%par%regions_path) > 0) &
            call nc_read(htopo%par%regions_path, htopo%par%regions_var, htopo%regions)
        if (len_trim(htopo%par%basins_path) > 0) &
            call nc_read(htopo%par%basins_path,  htopo%par%basins_var,  htopo%basins)
        if (len_trim(htopo%par%sectors_path) > 0) &
            call nc_read(htopo%par%sectors_path, htopo%par%sectors_var, htopo%sectors)

    end subroutine htopo_init

    subroutine htopo_update(htopo, dz_bed, dH_ice, z_sl)
        ! Current geometry on a hub finer than the ice sheet: the hi-res reference
        ! plus the models' anomalies (on the hub grid), with ice thickness clipped
        ! at 0. Each hub cell is either fully ice-covered or ice-free, so the
        ! grounded fraction is 0 or 1 from flotation, and the surface follows from
        ! the bed, the ice and sea level.
        type(htopo_class), intent(inout) :: htopo
        real(wp),          intent(in)    :: dz_bed(:,:)   ! [m] bed displacement
        real(wp),          intent(in)    :: dH_ice(:,:)   ! [m] change in ice thickness
        real(wp),          intent(in)    :: z_sl(:,:)     ! [m] sea level

        htopo%z_bed = htopo%z_bed_ref + dz_bed
        htopo%H_ice = max(htopo%H_ice_ref + dH_ice, 0.0_wp)
        htopo%z_sl  = z_sl

        htopo%f_grnd = 0.0_wp
        where (calc_H_grnd(htopo%H_ice, htopo%z_bed, htopo%z_sl, &
                           htopo%par%rho_ice, htopo%par%rho_sw) >= 0.0_wp) htopo%f_grnd = 1.0_wp
        htopo%z_srf = calc_z_srf(htopo%H_ice, htopo%z_bed, htopo%z_sl, &
                                 htopo%par%rho_ice, htopo%par%rho_sw)

    end subroutine htopo_update

    elemental function calc_H_grnd(H_ice, z_bed, z_sl, rho_ice, rho_sw) result(H_grnd)
        ! Ice overburden relative to flotation: >= 0 grounded, < 0 floating. Above
        ! sea level, the bed's height counts too, so ice-free land is grounded.
        real(wp), intent(in) :: H_ice, z_bed, z_sl, rho_ice, rho_sw
        real(wp) :: H_grnd

        if (z_sl > z_bed) then
            H_grnd = H_ice - (rho_sw/rho_ice)*(z_sl - z_bed)
        else
            H_grnd = H_ice + (z_bed - z_sl)
        end if

    end function calc_H_grnd

    elemental function calc_z_srf(H_ice, z_bed, z_sl, rho_ice, rho_sw) result(z_srf)
        ! Surface elevation: the top of grounded ice or of floating ice in
        ! hydrostatic equilibrium, whichever is higher (sea level if ice-free ocean).
        real(wp), intent(in) :: H_ice, z_bed, z_sl, rho_ice, rho_sw
        real(wp) :: z_srf

        z_srf = max(z_bed + H_ice, z_sl + (1.0_wp - rho_ice/rho_sw)*H_ice)

    end function calc_z_srf

    function htopo_ice_allowed(par, regions) result(allowed)
        ! Where ice is allowed, from the region codes on any grid (ice_codes_mode:
        ! "all", "include" = only on ice_codes, "exclude" = everywhere but ice_codes).
        type(htopo_par_class), intent(in) :: par
        real(wp),              intent(in) :: regions(:,:)
        logical :: allowed(size(regions,1),size(regions,2))

        allowed = codes_match(par%ice_codes_mode, par%ice_codes(1:par%n_ice_codes), regions)

    end function htopo_ice_allowed

    function htopo_relax_tau(par, regions) result(tau)
        ! Relaxation timescale of the ice toward the reference (Yelmo tau_relax,
        ! used with ytopo.topo_rel = -1), from the region codes on any grid:
        ! relax_tau where relax_codes_mode selects, -1 (free) elsewhere.
        type(htopo_par_class), intent(in) :: par
        real(wp),              intent(in) :: regions(:,:)
        real(wp) :: tau(size(regions,1),size(regions,2))

        tau = -1.0_wp
        where (codes_match(par%relax_codes_mode, par%relax_codes(1:par%n_relax_codes), regions)) &
            tau = par%relax_tau

    end function htopo_relax_tau

    function codes_match(mode, codes, mask) result(match)
        ! Cells of a code mask selected by mode: "none", "all", "include" (only
        ! on codes) or "exclude" (everywhere but codes).
        character(len=*), intent(in) :: mode
        real(wp),         intent(in) :: codes(:)
        real(wp),         intent(in) :: mask(:,:)
        logical :: match(size(mask,1),size(mask,2))

        integer :: k

        select case(trim(mode))
            case("none")
                match = .false.
            case("all")
                match = .true.
            case("include")
                match = .false.
                do k = 1, size(codes)
                    where (abs(mask - codes(k)) < tol_code) match = .true.
                end do
            case("exclude")
                match = .true.
                do k = 1, size(codes)
                    where (abs(mask - codes(k)) < tol_code) match = .false.
                end do
        end select

    end function codes_match

    function htopo_region_codes(htopo) result(codes)
        ! The code mask the named regions refer to (region_mask), on the hub grid.
        type(htopo_class), intent(in) :: htopo
        real(wp) :: codes(htopo%nx,htopo%ny)

        select case(trim(htopo%par%region_mask))
            case("regions")
                codes = htopo%regions
            case("basins")
                codes = htopo%basins
            case("sectors")
                codes = htopo%sectors
        end select

    end function htopo_region_codes

    subroutine htopo_fill_missing(htopo)
        ! Fill the gaps of the reference geometry (e.g. outside the coverage of
        ! the source dataset): no ice, the bed from the nearest valid cell, and
        ! the surface from the bed and the ice thickness, with sea level at 0.
        type(htopo_class), intent(inout) :: htopo

        integer :: n_bed, n_ice, n_srf

        n_bed = count(htopo%z_bed_ref == mv)
        n_ice = count(htopo%H_ice_ref == mv)
        n_srf = count(htopo%z_srf_ref == mv)
        if (n_bed + n_ice + n_srf == 0) return

        where (htopo%H_ice_ref == mv) htopo%H_ice_ref = 0.0_wp

        if (n_bed > 0) then
            if (n_bed < size(htopo%z_bed_ref)) call fill_nearest(htopo%z_bed_ref, mv)
            if (any(htopo%z_bed_ref == mv)) then
                write(*,*) ""
                write(*,*) "htopo_fill_missing:: error: missing bedrock elevations could not be filled."
                write(*,*) "  topo_path: ", trim(htopo%par%topo_path)
                write(*,*) "  z_bed:     ", trim(htopo%par%topo_names(1))
                write(*,*) "  missing:   ", count(htopo%z_bed_ref == mv), " of ", size(htopo%z_bed_ref)
                stop
            end if
        end if

        where (htopo%z_srf_ref == mv) &
            htopo%z_srf_ref = calc_z_srf(htopo%H_ice_ref, htopo%z_bed_ref, 0.0_wp, &
                                         htopo%par%rho_ice, htopo%par%rho_sw)

        write(*,*) "htopo_init:: filled missing values: z_bed ", n_bed, ", H_ice ", n_ice, &
                   ", z_srf ", n_srf, " of ", size(htopo%z_bed_ref)

    end subroutine htopo_fill_missing

    subroutine htopo_par_load(par, filename, group, domain, grid_name)
        type(htopo_par_class), intent(out) :: par
        character(len=*),      intent(in)  :: filename, group
        character(len=*),      intent(in)  :: domain, grid_name

        par%domain    = trim(domain)
        par%grid_name = trim(grid_name)

        ! Blank entries read as "" (nml_read leaves the value untouched); list
        ! entries not given keep mv / "".
        par%topo_names   = ""
        par%regions_path = ""
        par%basins_path  = ""
        par%sectors_path = ""
        par%ice_codes    = mv
        par%region_names = ""
        par%region_codes = mv
        par%relax_codes  = mv

        call nml_read(filename, group, "topo_path",      par%topo_path)
        call nml_read(filename, group, "topo_names",     par%topo_names)
        call nml_read(filename, group, "regions_path",   par%regions_path)
        call nml_read(filename, group, "regions_var",    par%regions_var)
        call nml_read(filename, group, "basins_path",    par%basins_path)
        call nml_read(filename, group, "basins_var",     par%basins_var)
        call nml_read(filename, group, "sectors_path",   par%sectors_path)
        call nml_read(filename, group, "sectors_var",    par%sectors_var)
        call nml_read(filename, group, "ice_codes_mode", par%ice_codes_mode)
        call nml_read(filename, group, "ice_codes",      par%ice_codes)
        call nml_read(filename, group, "region_names",   par%region_names)
        call nml_read(filename, group, "region_mask",    par%region_mask)
        call nml_read(filename, group, "region_codes",   par%region_codes)
        call nml_read(filename, group, "relax_codes_mode", par%relax_codes_mode)
        call nml_read(filename, group, "relax_codes",    par%relax_codes)
        call nml_read(filename, group, "relax_tau",      par%relax_tau)

        ! Resolve {domain}/{grid_name} against the hub grid.
        call parse_path(par%topo_path,    par%domain, par%grid_name)
        call parse_path(par%basins_path,  par%domain, par%grid_name)
        call parse_path(par%regions_path, par%domain, par%grid_name)
        call parse_path(par%sectors_path, par%domain, par%grid_name)

        par%n_ice_codes   = count(par%ice_codes /= mv)
        par%n_relax_codes = count(par%relax_codes /= mv)
        par%n_regions     = count(len_trim(par%region_names) > 0)

        select case(trim(par%ice_codes_mode))
            case("all")
                ! ice_codes not used
            case("include", "exclude")
                if (par%n_ice_codes == 0) call htopo_par_error(group, &
                    "ice_codes_mode = "//trim(par%ice_codes_mode)//" needs ice_codes.")
            case default
                call htopo_par_error(group, "ice_codes_mode must be all, include or exclude; got "// &
                                     trim(par%ice_codes_mode)//".")
        end select

        select case(trim(par%relax_codes_mode))
            case("none")
                ! relax_codes and relax_tau not used
            case("all", "include", "exclude")
                if (trim(par%relax_codes_mode) /= "all" .and. par%n_relax_codes == 0) &
                    call htopo_par_error(group, "relax_codes_mode = "//trim(par%relax_codes_mode)// &
                                         " needs relax_codes.")
                if (par%relax_tau <= 0.0_wp) call htopo_par_error(group, &
                    "relax_codes_mode = "//trim(par%relax_codes_mode)//" needs relax_tau > 0.")
            case default
                call htopo_par_error(group, "relax_codes_mode must be none, all, include or exclude; got "// &
                                     trim(par%relax_codes_mode)//".")
        end select

        select case(trim(par%region_mask))
            case("regions", "basins", "sectors")
            case default
                call htopo_par_error(group, "region_mask must be regions, basins or sectors; got "// &
                                     trim(par%region_mask)//".")
        end select

        if (any(par%region_codes(1:par%n_regions) == mv)) call htopo_par_error(group, &
            "region_codes needs one code per entry of region_names.")

    end subroutine htopo_par_load

    subroutine htopo_par_error(group, msg)
        character(len=*), intent(in) :: group, msg
        write(*,*) ""
        write(*,*) "htopo_par_load:: error in ["//trim(group)//"]: "//trim(msg)
        stop
    end subroutine htopo_par_error

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
        call nc_write(filename, "sectors", htopo%sectors, dim1="xc", dim2="yc", &
                      start=[1,1], long_name="Sector mask", units="")
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
