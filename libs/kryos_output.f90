module kryos_output
    ! Output of a kryos_domain: one file per module, each on its own grid (the
    ! grid is recorded inside the file, so names carry no grid suffix): 2D as
    ! <module>.nc, 1D timeseries as <module>_ts.nc.

    use ncio
    use coords,       only : grid_class, grid_cdo_read_desc
    use yelmo,        only : wp, yelmo_write_init, yelmo_write_step, yelmo_regions_write
    use marine_shelf, only : marshelf_class
    use fastisostasy, only : isos_class
    use smbpal,       only : smbpal_class
    use htopo,        only : htopo_write_init, htopo_write_step
    use timestepping, only : tstep_class
    use yelmox_climate, only : climate_file_base, climate_write_2D, climate_write_1D
    use cmip_output,  only : cmip_write_init, write_step_2D_cmip, write_step_1D_cmip
    use kryos,        only : kryos_domain, MAP_FLDR, remap, cadence_due

    implicit none
    private

    ! Default Yelmo 2D output variable sets. Public so a driver can reuse or
    ! extend them per context without editing the write routines, e.g.
    !   call domain_write_step(dom, outfldr, time, nms=[YELMO_VARS_2D, "my_var"])
    ! Any name valid for yelmo_write_var may be added. The heavy set adds
    ! bmb_shlf (coupled shelf melt) to the Yelmo default set.
    character(len=56), parameter :: YELMO_VARS_2D(23) = [ character(len=56) :: &
        "H_ice","z_srf","z_bed","mask_bed","uxy_b","uxy_s","uxy_bar", &
        "ux_bar","uy_bar","cb_ref","N_eff","beta","taub","taud","visc_bar", &
        "T_prime_b","hyd_W_til","mb_net","smb","bmb","cmb","z_sl","bmb_shlf" ]

    ! Small set: the minimal fields for frequent monitoring (yelmo_sm.nc).
    character(len=56), parameter :: YELMO_VARS_2D_SM(10) = [ character(len=56) :: &
        "H_ice","z_srf","z_bed","mask_bed","uxy_s","smb","bmb","z_sl", &
        "hyd_W_til","T_prime_b" ]

    public :: YELMO_VARS_2D, YELMO_VARS_2D_SM
    public :: domain_write_init, domain_write_step, domain_write_1D
    public :: domain_write_init_sm, domain_write_step_sm
    public :: domain_write_cmip

contains

    ! ----- output (one file per module, each on its own grid; the grid is
    !        recorded inside the file, so names carry no grid suffix: 2D as
    !        <module>.nc, 1D timeseries as <module>_ts.nc) ---

    function io_fname(outfldr, base) result(fnm)
        character(len=*), intent(in) :: outfldr, base
        character(len=512) :: fnm
        fnm = trim(outfldr)//trim(base)//".nc"
    end function io_fname

    subroutine domain_write_init(dom, outfldr, time)
        ! Create the enabled per-module 2D output files (dims + static fields).
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time

        if (dom%ctl%write_yelmo) &
            call yelmo_write_init(dom%yelmo, trim(io_fname(outfldr,"yelmo")), &
                                  time_init=time, units="years")
        if (dom%ctl%write_htopo) &
            call htopo_write_init(dom%topo, trim(io_fname(outfldr,"htopo")), time_init=time)
        if (dom%ctl%write_isos) &
            call io_dims_init(trim(io_fname(outfldr,"isos")),   dom%ctl%grid_isos,  time)
        if (dom%ctl%write_mshlf) &
            call io_dims_init(trim(io_fname(outfldr,"mshlf")),  dom%ctl%grid_mshlf, time)
        if (dom%ctl%write_smb) &
            call io_dims_init(trim(io_fname(outfldr,"smbpal")), dom%ctl%grid_smb,   time)
        if (dom%ctl%write_clim) &
            call io_dims_init(trim(io_fname(outfldr,climate_file_base(dom%cl))), dom%ctl%grid_clim, time)
    end subroutine domain_write_init

    subroutine domain_write_step(dom, outfldr, time, nms)
        ! Append one time record to each enabled per-module 2D output file. The
        ! Yelmo 2D variable set defaults to YELMO_VARS_2D; pass `nms` to override
        ! or extend it for a specific context (see YELMO_VARS_2D above).
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time
        character(len=*), intent(in), optional :: nms(:)

        character(len=56), allocatable :: yvars(:)

        if (present(nms)) then
            yvars = nms
        else
            yvars = YELMO_VARS_2D
        end if

        if (dom%ctl%write_yelmo) &
            call yelmo_write_step(dom%yelmo, trim(io_fname(outfldr,"yelmo")), &
                                  time, nms=yvars, compare_pd=.FALSE.)
        if (dom%ctl%write_htopo) &
            call htopo_write_step(dom%topo, trim(io_fname(outfldr,"htopo")), time)
        if (dom%ctl%write_isos) &
            call isos_write_step(dom%isos, trim(io_fname(outfldr,"isos")), time)
        if (dom%ctl%write_mshlf) &
            call mshlf_write_step(dom%mshlf, trim(io_fname(outfldr,"mshlf")), time)
        if (dom%ctl%write_smb) &
            call smb_write_step(dom%smb, trim(io_fname(outfldr,"smbpal")), time)
        if (dom%ctl%write_clim) &
            call clim_write_step(dom, trim(io_fname(outfldr,climate_file_base(dom%cl))), time)
    end subroutine domain_write_step

    subroutine domain_write_init_sm(dom, outfldr, time)
        ! Create the small Yelmo 2D output file (yelmo_sm.nc): a reduced field set
        ! for frequent monitoring, written on the sm cadence (tm_2Dsm). Same grid
        ! and file conventions as the heavy yelmo.nc.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time

        if (dom%ctl%write_yelmo) &
            call yelmo_write_init(dom%yelmo, trim(io_fname(outfldr,"yelmo_sm")), &
                                  time_init=time, units="years")
    end subroutine domain_write_init_sm

    subroutine domain_write_step_sm(dom, outfldr, time, nms)
        ! Append one record to yelmo_sm.nc: the core geometry/velocity plus basal
        ! hydrology and temperature -- the minimal set needed to watch a run at a
        ! higher cadence than the heavy yelmo.nc. The variable set defaults to
        ! YELMO_VARS_2D_SM; pass `nms` to override or extend it for a context.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time
        character(len=*), intent(in), optional :: nms(:)

        character(len=56), allocatable :: yvars(:)

        if (present(nms)) then
            yvars = nms
        else
            yvars = YELMO_VARS_2D_SM
        end if

        if (dom%ctl%write_yelmo) &
            call yelmo_write_step(dom%yelmo, trim(io_fname(outfldr,"yelmo_sm")), &
                                  time, nms=yvars, compare_pd=.FALSE.)
    end subroutine domain_write_step_sm

    subroutine domain_write_1D(dom, outfldr, time, init)
        ! Write 1D timeseries: Yelmo regional aggregates + isostasy diagnostics.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time
        logical, intent(in), optional   :: init

        logical :: is_init
        character(len=512) :: fnm_isos
        real(wp), allocatable :: H_ice_c(:,:), f_grnd_c(:,:)

        is_init = .false.
        if (present(init)) is_init = init

        if (dom%ctl%write_yelmo) then
            if (is_init) then
                call yelmo_regions_write(dom%yelmo, time, init=.TRUE., units="years")
            else
                call yelmo_regions_write(dom%yelmo, time)
            end if
        end if

        if (dom%ctl%write_isos) then
            fnm_isos = trim(outfldr)//"isos_ts.nc"
            if (is_init) call isos_write_1D_init(trim(fnm_isos), time)
            call isos_write_1D_step(dom%isos, trim(fnm_isos), time)
        end if

        ! The climate's own 1D diagnostics (esm), over the ice on the climate grid.
        if (dom%ctl%write_clim .and. dom%ctl%with_climate) then
            call remap(dom, dom%topo%H_ice,  dom%ctl%grid_hub, H_ice_c,  dom%ctl%grid_clim, "bilin")
            call remap(dom, dom%topo%f_grnd, dom%ctl%grid_hub, f_grnd_c, dom%ctl%grid_clim, "bilin")
            call climate_write_1D(dom%cl, dom%clim, &
                    trim(outfldr)//trim(climate_file_base(dom%cl))//"_ts.nc", time, &
                    H_ice_c, f_grnd_c, is_init)
        end if
    end subroutine domain_write_1D

    subroutine domain_write_cmip(dom, outfldr, ts, init)
        ! CMIP/ISMIP-formatted output ([output] write_cmip): init creates the
        ! files, otherwise one record every dt_cmip of elapsed time.
        type(kryos_domain), intent(inout) :: dom
        character(len=*),   intent(in)    :: outfldr
        type(tstep_class),  intent(in)    :: ts
        logical, intent(in), optional     :: init

        logical :: is_init

        if (.not. dom%ctl%write_cmip) return

        is_init = .false.
        if (present(init)) is_init = init

        if (is_init) then
            call cmip_write_init(dom%yelmo, trim(io_fname(outfldr,"yelmo_cmip")), &
                                 trim(io_fname(outfldr,"yelmo_ts_cmip")), ts%time)
        else if (cadence_due(ts%time_elapsed, dom%ctl%dt_cmip)) then
            call write_step_2D_cmip(dom%yelmo, dom%mshlf, trim(io_fname(outfldr,"yelmo_cmip")), ts%time)
            call write_step_1D_cmip(dom%yelmo, dom%mshlf, trim(io_fname(outfldr,"yelmo_ts_cmip")), ts%time)
        end if
    end subroutine domain_write_cmip

    ! --- private output helpers ---

    subroutine io_dims_init(filename, grid_name, time)
        ! Create a 2D output file with xc/yc (from the grid table) + time dims.
        character(len=*), intent(in) :: filename, grid_name
        real(wp),         intent(in) :: time
        type(grid_class) :: g
        call grid_cdo_read_desc(g, trim(grid_name), MAP_FLDR)
        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=g%G%x, units="km")
        call nc_write_dim(filename, "yc", x=g%G%y, units="km")
        call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, units="year", unlimited=.TRUE.)
    end subroutine io_dims_init

    subroutine io_var2D(filename, vnm, var, n, ncid, units, long_name)
        character(len=*), intent(in) :: filename, vnm, units, long_name
        real(wp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                      units=units, long_name=long_name)
    end subroutine io_var2D

    subroutine io_ts(filename, vnm, val, n, ncid, units, long_name)
        character(len=*), intent(in) :: filename, vnm, units, long_name
        real(wp),         intent(in) :: val
        integer,          intent(in) :: n, ncid
        call nc_write(filename, vnm, val, dim1="time", start=[n], count=[1], ncid=ncid, &
                      units=units, long_name=long_name)
    end subroutine io_ts

    subroutine isos_write_step(isos, filename, time)
        type(isos_class), intent(in) :: isos
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "z_bed", isos%out%z_bed, n, ncid, "m", "Bedrock elevation")
        call io_var2D(filename, "z_ss",  isos%out%z_ss,  n, ncid, "m", "Sea-surface height")
        call io_var2D(filename, "w",     isos%out%w,     n, ncid, "m", "Viscous displacement")
        call io_var2D(filename, "we",    isos%out%we,    n, ncid, "m", "Elastic displacement")
        call nc_close(ncid)
    end subroutine isos_write_step

    subroutine mshlf_write_step(mshlf, filename, time)
        type(marshelf_class), intent(in) :: mshlf
        character(len=*),     intent(in) :: filename
        real(wp),             intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "bmb_shlf", mshlf%now%bmb_shlf, n, ncid, "m/yr", "Shelf basal mass balance")
        call io_var2D(filename, "T_shlf",   mshlf%now%T_shlf,   n, ncid, "K", "Shelf temperature")
        call io_var2D(filename, "tf_shlf",  mshlf%now%tf_shlf,  n, ncid, "K", "Thermal forcing")
        call nc_close(ncid)
    end subroutine mshlf_write_step

    subroutine smb_write_step(smb, filename, time)
        type(smbpal_class), intent(in) :: smb
        character(len=*),   intent(in) :: filename
        real(wp),           intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "smb",  smb%ann%smb,  n, ncid, "m ie/yr", "Surface mass balance")
        call io_var2D(filename, "tsrf", smb%ann%tsrf, n, ncid, "K", "Surface temperature")
        call nc_close(ncid)
    end subroutine smb_write_step

    subroutine clim_write_step(dom, filename, time)
        ! One record of the climate's 2D file (fields chosen by the backend).
        type(kryos_domain), intent(in) :: dom
        character(len=*),   intent(in) :: filename
        real(wp),           intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call climate_write_2D(dom%cl, dom%clim, filename, ncid, n)
        call nc_close(ncid)
    end subroutine clim_write_step

    subroutine isos_write_1D_init(filename, time)
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        call nc_create(filename)
        call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, units="year", unlimited=.TRUE.)
    end subroutine isos_write_1D_init

    subroutine isos_write_1D_step(isos, filename, time)
        type(isos_class), intent(in) :: isos
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        integer :: ncid, n, np
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        np = size(isos%out%z_bed)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_ts(filename, "bsl",        isos%now%bsl,                       n, ncid, "m", "Barystatic sea level")
        call io_ts(filename, "z_bed_mean", sum(isos%out%z_bed)/real(np,wp),    n, ncid, "m", "Mean bedrock elevation")
        call io_ts(filename, "z_bed_min",  minval(isos%out%z_bed),             n, ncid, "m", "Min bedrock elevation")
        call io_ts(filename, "z_bed_max",  maxval(isos%out%z_bed),             n, ncid, "m", "Max bedrock elevation")
        call io_ts(filename, "w_mean",     sum(isos%out%w)/real(np,wp),        n, ncid, "m", "Mean viscous displacement")
        call io_ts(filename, "we_mean",    sum(isos%out%we)/real(np,wp),       n, ncid, "m", "Mean elastic displacement")
        call nc_close(ncid)
    end subroutine isos_write_1D_step

end module kryos_output
