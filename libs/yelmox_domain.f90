module yelmox_domain
    ! Multigrid yelmox: one region's full model state bundled as an ice_domain,
    ! advanced by composable step_* primitives.
    !
    ! The hi-res topography hub (dom%topo) is the geometry source of truth: it is
    ! refreshed from the prognostic models each step (refresh_htopo), and the
    ! coupling steps remap fields to/from it through the domain's coupler
    ! (dom%cpl). Each helper model runs on its own grid, named in domain_ctl; a
    ! component's grid is just a string, so e.g. marine_shelf can run on the
    ! hi-res hub grid or on the Yelmo grid simply by changing grid_mshlf.
    !
    ! Status: marine_shelf coupled through the hub on a configurable grid. Other
    ! components remain on the Yelmo grid for now.

    use nml,          only : nml_read
    use timestepping, only : tstep_class
    use coords,       only : grid_class, grid_cdo_read_desc
    use yelmo,        only : yelmo_class, wp, yelmo_init, yelmo_update, &
                             yelmo_init_state, yelmo_print_bound, &
                             yelmo_restart_write, yelmo_restart_read
    use marine_shelf, only : marshelf_class, marshelf_init, marshelf_update, &
                             marshelf_update_shelf, marshelf_restart_write, &
                             marshelf_restart_read
    use fastisostasy, only : isos_class, isos_init, isos_update, isos_init_ref, &
                             isos_init_state, isos_restart_write, isos_restart_read, &
                             bsl_class, bsl_init, bsl_update, bsl_restart_write
    use snapclim,     only : snapclim_class, snapclim_init, snapclim_update
    use smbpal,       only : smbpal_class, smbpal_init, smbpal_update_monthly, &
                             smbpal_update_monthly_equil
    use sediments,    only : sediments_class, sediments_init
    use geothermal,   only : geothermal_class, geothermal_init
    use htopo,        only : htopo_class, htopo_init
    use coupler,      only : coupler_class, coupler_init, coupler_prime, remap

    implicit none
    private

    character(len=*), parameter :: MAP_FLDR = "maps"

    type domain_ctl
        ! Run control (from the [ctrl] namelist group).
        character(len=56) :: tstep_method = "const"
        real(wp) :: tstep_const = 0.0_wp
        real(wp) :: time_init   = 0.0_wp
        real(wp) :: time_end    = 0.0_wp
        real(wp) :: time_equil  = 0.0_wp
        real(wp) :: dtt         = 10.0_wp
        real(wp) :: dt_restart  = 0.0_wp
        real(wp) :: dt_clim     = 10.0_wp   ! [yr] snapclim snapshot update frequency
        character(len=56) :: equil_method = "none"
        character(len=56) :: smb_method   = "smbpal"

        ! Which components are active in this domain's coupling sequence.
        logical :: with_ice_sheet    = .true.
        logical :: with_isostasy     = .true.
        logical :: with_marine_shelf = .true.
        logical :: with_climate      = .true.

        ! Domain + grid names (the source of truth for remap keys).
        character(len=256) :: domain     = ""   ! e.g. "Antarctica"
        character(len=256) :: grid_name  = ""   ! hi-res reference (htopo), highest res
        character(len=256) :: grid_yelmo = ""   ! Yelmo grid
        character(len=256) :: grid_mshlf = ""   ! marine-shelf grid ([coupling]; default = grid_name)
        real(wp) :: dx_mshlf = 0.0_wp           ! marine-shelf grid spacing (Yelmo dx units)

        ! Restart bundle folder ([coupling]); "None" = cold start.
        character(len=512) :: restart = "None"
    end type domain_ctl

    type ice_domain
        type(yelmo_class)      :: yelmo
        type(marshelf_class)   :: mshlf
        type(isos_class)       :: isos
        type(snapclim_class)   :: snp
        type(smbpal_class)     :: smb
        type(bsl_class)        :: bsl
        type(sediments_class)  :: sed
        type(geothermal_class) :: gthrm
        type(htopo_class)      :: topo    ! hi-res geometry reference hub
        type(coupler_class)    :: cpl     ! this region's grid resolution + map cache
        type(domain_ctl)       :: ctl
    end type ice_domain

    public :: domain_ctl, ice_domain
    public :: domain_init, domain_init_state, yelmox_step
    public :: domain_restart_write, domain_restart_read
    public :: step_isostasy, step_icesheet, step_climate, refresh_htopo, step_marine_shelf

contains

    subroutine domain_init(dom, path_par, time, time_rel)
        ! Initialize all sub-models of one domain, load the hi-res reference hub,
        ! prime the Yelmo<->hub maps, and place marine_shelf on its configured grid.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: path_par
        real(wp),         intent(in)    :: time       ! model time
        real(wp),         intent(in)    :: time_rel   ! time before present [yr]

        character(len=256)    :: domain
        type(grid_class)      :: grid_m, grid_y
        integer               :: nx_m, ny_m
        real(wp), allocatable :: regions_m(:,:), basins_m(:,:)

        ! --- run control ---
        call domain_ctl_load(dom%ctl, path_par)

        ! --- ice sheet (grid read from file) ---
        call yelmo_init(dom%yelmo, filename=path_par, grid_def="file", time=time)
        domain = trim(dom%yelmo%par%domain)
        dom%ctl%domain     = trim(domain)
        dom%ctl%grid_yelmo = trim(dom%yelmo%par%grid_name)

        ! --- external forcing models (climate/smb/isostasy on the Yelmo grid) ---
        call bsl_init(dom%bsl, path_par, time_rel)

        call isos_init(dom%isos, path_par, "isos", dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                       dom%yelmo%grd%dx, dom%yelmo%grd%dy)

        call snapclim_init(dom%snp, path_par, domain, dom%yelmo%par%grid_name, &
                           dom%yelmo%grd%nx, dom%yelmo%grd%ny, dom%yelmo%bnd%basins)

        call smbpal_init(dom%smb, path_par, x=dom%yelmo%grd%xc, y=dom%yelmo%grd%yc, &
                         lats=dom%yelmo%grd%lat)

        call sediments_init(dom%sed, path_par, dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                            domain, dom%yelmo%par%grid_name)
        dom%yelmo%bnd%H_sed = dom%sed%now%H

        call geothermal_init(dom%gthrm, path_par, dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                             domain, dom%yelmo%par%grid_name)
        dom%yelmo%bnd%Q_geo = dom%gthrm%now%ghf

        ! --- hi-res reference hub + coupler ---
        call htopo_init(dom%topo, path_par, "htopo")
        dom%ctl%grid_name = trim(dom%topo%par%grid_name)
        if (len_trim(dom%ctl%grid_mshlf) == 0) dom%ctl%grid_mshlf = trim(dom%ctl%grid_name)

        ! Grids resolve from maps/grid_<name>.txt; prime the Yelmo<->hub maps.
        call coupler_init(dom%cpl)
        call coupler_prime(dom%cpl, dom%ctl%grid_yelmo, dom%ctl%grid_name, "bilin")  ! Yelmo -> hub
        call coupler_prime(dom%cpl, dom%ctl%grid_name, dom%ctl%grid_yelmo, "con")    ! hub -> Yelmo

        ! --- marine_shelf on its configured grid ---
        call grid_cdo_read_desc(grid_m, trim(dom%ctl%grid_mshlf), MAP_FLDR)
        call grid_cdo_read_desc(grid_y, trim(dom%ctl%grid_yelmo), MAP_FLDR)
        nx_m = grid_m%G%nx
        ny_m = grid_m%G%ny
        ! Grid spacing in Yelmo dx units, scaled by the resolution ratio.
        dom%ctl%dx_mshlf = dom%yelmo%grd%dx * (grid_m%G%dx / grid_y%G%dx)

        ! Region/basin masks on the mshlf grid (from the hub).
        call remap_or_copy_2D(dom, dom%topo%regions, dom%ctl%grid_name, regions_m, dom%ctl%grid_mshlf, "nn")
        call remap_or_copy_2D(dom, dom%topo%basins,  dom%ctl%grid_name, basins_m,  dom%ctl%grid_mshlf, "nn")

        call marshelf_init(dom%mshlf, path_par, "marine_shelf", nx_m, ny_m, &
                           domain, trim(dom%ctl%grid_mshlf), regions_m, basins_m)

    end subroutine domain_init

    subroutine domain_init_state(dom, ts)
        ! Build the initial boundary state and initialize the Yelmo state
        ! variables. Lifted from yelmox.f90's "update initial boundary conditions"
        ! block (minimal core: no smb_simple / domain-special / optimization).
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        ! Sea level + isostasy reference state
        call bsl_update(dom%bsl, ts%time_rel)
        call isos_init_ref(dom%isos, dom%yelmo%bnd%z_bed_ref, dom%yelmo%bnd%H_ice_ref)
        call isos_init_state(dom%isos, dom%yelmo%bnd%z_bed, dom%yelmo%tpo%now%H_ice, ts%time, dom%bsl)
        dom%yelmo%bnd%z_bed = dom%isos%out%z_bed
        dom%yelmo%bnd%z_sl  = dom%isos%out%z_ss

        ! Climate + surface mass balance (note: init uses time_rel for snapclim)
        call snapclim_update(dom%snp, z_srf=dom%yelmo%tpo%now%z_srf, time=ts%time_rel, &
                             domain=dom%ctl%domain, dx=dom%yelmo%grd%dx, basins=dom%yelmo%bnd%basins)

        if (trim(dom%smb%par%abl_method) == "itm") then
            call smbpal_update_monthly_equil(dom%smb, dom%snp%now%tas, dom%snp%now%pr, &
                    dom%yelmo%tpo%now%z_srf, dom%yelmo%tpo%now%H_ice, ts%time_rel, time_equil=100.0_wp)
        end if
        call smbpal_update_monthly(dom%smb, dom%snp%now%tas, dom%snp%now%pr, &
                dom%yelmo%tpo%now%z_srf, dom%yelmo%tpo%now%H_ice, ts%time_rel)
        dom%yelmo%bnd%smb   = dom%smb%ann%smb * dom%yelmo%bnd%c%conv_we_ie * 1e-3
        dom%yelmo%bnd%T_srf = dom%smb%ann%tsrf

        ! Marine shelf: refresh the hub, then run mshlf through it.
        call refresh_htopo(dom)
        call step_marine_shelf(dom, ts)

        ! Initialize state variables (dyn, therm, mat) with a cold base
        call yelmo_print_bound(dom%yelmo%bnd)
        call yelmo_init_state(dom%yelmo, time=ts%time, thrm_method="robin-cold")

    end subroutine domain_init_state

    subroutine domain_restart_write(dom, time, fldr)
        ! Write a restart bundle: a folder (per time, or `fldr`) holding one
        ! restart file per stateful sub-model with fixed names. The hi-res hub is
        ! not written -- it is rebuilt by refresh_htopo from the restored models.
        type(ice_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time
        character(len=*), intent(in), optional :: fldr

        character(len=1024) :: outfldr
        character(len=32)   :: time_str

        if (present(fldr)) then
            outfldr = trim(fldr)
        else
            write(time_str,"(f20.3)") time*1e-3
            outfldr = "restart-"//trim(adjustl(time_str))//"-kyr"
        end if

        call execute_command_line('mkdir -p "'//trim(outfldr)//'"')

        call bsl_restart_write(dom%bsl,      trim(outfldr)//"/bsl_restart.nc",   time)
        call isos_restart_write(dom%isos,    trim(outfldr)//"/isos_restart.nc",  time)
        call yelmo_restart_write(dom%yelmo,  trim(outfldr)//"/yelmo_restart.nc", time)
        call marshelf_restart_write(dom%mshlf, trim(outfldr)//"/marine_shelf.nc", time)

        write(*,*) "domain_restart_write:: wrote bundle "//trim(outfldr)
    end subroutine domain_restart_write

    subroutine domain_restart_read(dom, fldr, time)
        ! Restore all stateful sub-models from a restart bundle folder. One call,
        ! fixed filenames -- no per-file configuration. bsl is not read (it is
        ! reconstructed by bsl_init + bsl_update at the restart time).
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: fldr
        real(wp),         intent(in)    :: time

        call isos_restart_read(dom%isos,       trim(fldr)//"/isos_restart.nc",  time)
        call yelmo_restart_read(dom%yelmo,     trim(fldr)//"/yelmo_restart.nc", time)
        call marshelf_restart_read(dom%mshlf,  trim(fldr)//"/marine_shelf.nc")

        write(*,*) "domain_restart_read:: restored bundle "//trim(fldr)
    end subroutine domain_restart_read

    subroutine domain_ctl_load(ctl, path_par)
        type(domain_ctl), intent(inout) :: ctl
        character(len=*), intent(in)    :: path_par

        call nml_read(path_par, "ctrl", "tstep_method",   ctl%tstep_method)
        call nml_read(path_par, "ctrl", "tstep_const",    ctl%tstep_const)
        call nml_read(path_par, "ctrl", "time_init",      ctl%time_init)
        call nml_read(path_par, "ctrl", "time_end",       ctl%time_end)
        call nml_read(path_par, "ctrl", "time_equil",     ctl%time_equil)
        call nml_read(path_par, "ctrl", "dtt",            ctl%dtt)
        call nml_read(path_par, "ctrl", "dt_restart",     ctl%dt_restart)
        call nml_read(path_par, "ctrl", "with_ice_sheet", ctl%with_ice_sheet)
        call nml_read(path_par, "ctrl", "with_isostasy",  ctl%with_isostasy)
        call nml_read(path_par, "ctrl", "equil_method",   ctl%equil_method)
        ctl%smb_method = "smbpal"
        call nml_read(path_par, "ctrl", "smb_method",     ctl%smb_method)

        ! Component grids + restart bundle ([coupling]).
        ctl%grid_mshlf = ""
        call nml_read(path_par, "coupling", "grid_mshlf", ctl%grid_mshlf)
        ctl%restart = "None"
        call nml_read(path_par, "coupling", "restart",    ctl%restart)
    end subroutine domain_ctl_load

    subroutine yelmox_step(dom, ts)
        ! Advance one domain by one coupling step. Component order is fixed here;
        ! each step_* is a no-op when its ctl flag is off.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        call step_isostasy(dom, ts)
        call step_icesheet(dom, ts)
        call step_climate(dom, ts)
        call refresh_htopo(dom)          ! hi-res geometry mirror, from the models
        call step_marine_shelf(dom, ts)
    end subroutine yelmox_step

    subroutine step_isostasy(dom, ts)
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        if (.not. dom%ctl%with_isostasy) return

        call bsl_update(dom%bsl, ts%time_rel)
        call isos_update(dom%isos, dom%yelmo%tpo%now%H_ice, ts%time, dom%bsl, &
                         dwdt_corr=dom%yelmo%bnd%dzbdt_corr)
        dom%yelmo%bnd%z_bed = dom%isos%out%z_bed
        dom%yelmo%bnd%z_sl  = dom%isos%out%z_ss
    end subroutine step_isostasy

    subroutine step_icesheet(dom, ts)
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        if (.not. dom%ctl%with_ice_sheet) return
        if (ts%n == 0 .and. dom%yelmo%par%use_restart) return

        call yelmo_update(dom%yelmo, ts%time)
    end subroutine step_icesheet

    subroutine step_climate(dom, ts)
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        if (.not. dom%ctl%with_climate) return

        ! snapclim snapshot, updated on the dt_clim cadence
        if (mod(nint(ts%time_elapsed*100), nint(dom%ctl%dt_clim*100)) == 0) then
            call snapclim_update(dom%snp, z_srf=dom%yelmo%tpo%now%z_srf, time=ts%time, &
                                 domain=dom%ctl%domain, dx=dom%yelmo%grd%dx, basins=dom%yelmo%bnd%basins)
        end if

        ! surface mass balance (smbpal)
        call smbpal_update_monthly(dom%smb, dom%snp%now%tas, dom%snp%now%pr, &
                dom%yelmo%tpo%now%z_srf, dom%yelmo%tpo%now%H_ice, ts%time_rel)
        dom%yelmo%bnd%smb   = dom%smb%ann%smb * dom%yelmo%bnd%c%conv_we_ie * 1e-3
        dom%yelmo%bnd%T_srf = dom%smb%ann%tsrf
    end subroutine step_climate

    subroutine refresh_htopo(dom)
        ! Refresh the hi-res geometry hub from the prognostic models (Yelmo grid
        ! -> hub grid, bilinear). The hub is then the geometry source for the
        ! coupling steps. Static masks (regions/basins) are not refreshed.
        type(ice_domain), intent(inout) :: dom

        call remap_or_copy_2D(dom, dom%yelmo%tpo%now%H_ice,  dom%ctl%grid_yelmo, &
                              dom%topo%H_ice,  dom%ctl%grid_name, "bilin")
        call remap_or_copy_2D(dom, dom%yelmo%bnd%z_bed,      dom%ctl%grid_yelmo, &
                              dom%topo%z_bed,  dom%ctl%grid_name, "bilin")
        call remap_or_copy_2D(dom, dom%yelmo%tpo%now%f_grnd, dom%ctl%grid_yelmo, &
                              dom%topo%f_grnd, dom%ctl%grid_name, "bilin")
        call remap_or_copy_2D(dom, dom%yelmo%bnd%z_sl,       dom%ctl%grid_yelmo, &
                              dom%topo%z_sl,   dom%ctl%grid_name, "bilin")
        call remap_or_copy_2D(dom, dom%yelmo%tpo%now%z_srf,  dom%ctl%grid_yelmo, &
                              dom%topo%z_srf,  dom%ctl%grid_name, "bilin")
    end subroutine refresh_htopo

    subroutine step_marine_shelf(dom, ts)
        ! Run marine_shelf on its own grid: geometry/masks from the hub, ocean
        ! forcing from snapclim, outputs aggregated back to the Yelmo grid.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: H_ice_m(:,:), z_bed_m(:,:), f_grnd_m(:,:), z_sl_m(:,:)
        real(wp), allocatable :: regions_m(:,:), basins_m(:,:)
        real(wp), allocatable :: to_m(:,:,:), so_m(:,:,:), dto_m(:,:,:), dto_y(:,:,:)
        real(wp), allocatable :: bmb_y(:,:), Tshlf_y(:,:)
        character(len=256) :: gm, gn, gy

        if (.not. dom%ctl%with_marine_shelf) return

        gm = trim(dom%ctl%grid_mshlf)
        gn = trim(dom%ctl%grid_name)
        gy = trim(dom%ctl%grid_yelmo)

        ! geometry + masks: hub -> mshlf grid
        call remap_or_copy_2D(dom, dom%topo%H_ice,   gn, H_ice_m,   gm, "bilin")
        call remap_or_copy_2D(dom, dom%topo%z_bed,   gn, z_bed_m,   gm, "bilin")
        call remap_or_copy_2D(dom, dom%topo%f_grnd,  gn, f_grnd_m,  gm, "bilin")
        call remap_or_copy_2D(dom, dom%topo%z_sl,    gn, z_sl_m,    gm, "bilin")
        call remap_or_copy_2D(dom, dom%topo%regions, gn, regions_m, gm, "nn")
        call remap_or_copy_2D(dom, dom%topo%basins,  gn, basins_m,  gm, "nn")

        ! ocean forcing (3D): snapclim (Yelmo grid) -> mshlf grid
        call remap_or_copy_3D(dom, dom%snp%now%to_ann, gy, to_m, gm, "bilin")
        call remap_or_copy_3D(dom, dom%snp%now%so_ann, gy, so_m, gm, "bilin")
        dto_y = dom%snp%now%to_ann - dom%snp%clim0%to_ann
        call remap_or_copy_3D(dom, dto_y, gy, dto_m, gm, "bilin")

        ! run marine_shelf on grid_mshlf
        call marshelf_update_shelf(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, basins_m, z_sl_m, &
                dom%ctl%dx_mshlf, dom%snp%now%depth, to_m, so_m, dto_ann=dto_m)
        call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, regions_m, basins_m, &
                z_sl_m, dx=dom%ctl%dx_mshlf)

        ! aggregate outputs -> Yelmo grid (conservative)
        call remap_or_copy_2D(dom, dom%mshlf%now%bmb_shlf, gm, bmb_y,   gy, "con")
        call remap_or_copy_2D(dom, dom%mshlf%now%T_shlf,   gm, Tshlf_y, gy, "con")
        dom%yelmo%bnd%bmb_shlf = bmb_y
        dom%yelmo%bnd%T_shlf   = Tshlf_y
    end subroutine step_marine_shelf

    ! ----- remap helpers: identity-copy when src == dst, else via the coupler ---

    subroutine remap_or_copy_2D(dom, var_src, src, var_dst, dst, method)
        type(ice_domain),      intent(inout) :: dom
        real(wp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst, method
        real(wp), allocatable, intent(inout) :: var_dst(:,:)

        if (trim(src) == trim(dst)) then
            if (allocated(var_dst)) then
                if (size(var_dst,1) /= size(var_src,1) .or. &
                    size(var_dst,2) /= size(var_src,2)) deallocate(var_dst)
            end if
            if (.not. allocated(var_dst)) allocate(var_dst(size(var_src,1), size(var_src,2)))
            var_dst = var_src
        else
            call remap(dom%cpl, var_src, src, var_dst, dst, method=method)
        end if
    end subroutine remap_or_copy_2D

    subroutine remap_or_copy_3D(dom, var_src, src, var_dst, dst, method)
        type(ice_domain),      intent(inout) :: dom
        real(wp),              intent(in)    :: var_src(:,:,:)
        character(len=*),      intent(in)    :: src, dst, method
        real(wp), allocatable, intent(inout) :: var_dst(:,:,:)

        if (trim(src) == trim(dst)) then
            if (allocated(var_dst)) then
                if (size(var_dst,1) /= size(var_src,1) .or. &
                    size(var_dst,2) /= size(var_src,2) .or. &
                    size(var_dst,3) /= size(var_src,3)) deallocate(var_dst)
            end if
            if (.not. allocated(var_dst)) &
                allocate(var_dst(size(var_src,1), size(var_src,2), size(var_src,3)))
            var_dst = var_src
        else
            call remap(dom%cpl, var_src, src, var_dst, dst, method=method)
        end if
    end subroutine remap_or_copy_3D

end module yelmox_domain
