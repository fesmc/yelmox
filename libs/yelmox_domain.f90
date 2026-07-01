module yelmox_domain
    ! Multigrid yelmox: one region's full model state bundled as an ice_domain,
    ! advanced by composable step_* primitives.
    !
    ! Each helper model may live on its own grid; fields that cross a grid
    ! boundary are remapped through the domain's coupler (dom%cpl) at the moment
    ! of coupling. The step_* routines are the reusable pieces -- they are guarded
    ! internally by domain_ctl flags, and yelmox_step composes them. Bipolar runs
    ! are just an array of ice_domain looped through yelmox_step, so north and
    ! south share no state.
    !
    ! Status: minimal core, everything on the Yelmo grid (lift-and-adapt of
    ! yelmox.f90's init + main loop). domain_init sets up the sub-models, the
    ! htopo reference hub and the coupler maps; domain_init_state builds the
    ! initial boundary state; yelmox_step advances one coupling step. Moving
    ! marine_shelf onto the hi-res hub (via the coupler) comes next.

    use nml,          only : nml_read
    use timestepping, only : tstep_class
    use yelmo,        only : yelmo_class, wp, yelmo_init, yelmo_update, &
                             yelmo_init_state, yelmo_print_bound
    use marine_shelf, only : marshelf_class, marshelf_init, marshelf_update, &
                             marshelf_update_shelf
    use fastisostasy, only : isos_class, isos_init, isos_update, isos_init_ref, &
                             isos_init_state, bsl_class, bsl_init, bsl_update
    use snapclim,     only : snapclim_class, snapclim_init, snapclim_update
    use smbpal,       only : smbpal_class, smbpal_init, smbpal_update_monthly, &
                             smbpal_update_monthly_equil
    use sediments,    only : sediments_class, sediments_init
    use geothermal,   only : geothermal_class, geothermal_init
    use htopo,        only : htopo_class, htopo_init
    use coupler,      only : coupler_class, coupler_init, coupler_prime, remap

    implicit none
    private

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
        character(len=256) :: grid_mshlf = ""   ! marine-shelf grid (defaults to grid_name)
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
    public :: step_isostasy, step_icesheet, step_climate, step_marine_shelf

contains

    subroutine domain_init(dom, path_par, time, time_rel)
        ! Initialize all sub-models of one domain (currently all on the Yelmo
        ! grid), load the hi-res reference hub, and prime the Yelmo<->topo maps.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: path_par
        real(wp),         intent(in)    :: time       ! model time
        real(wp),         intent(in)    :: time_rel   ! time before present [yr]

        character(len=256) :: domain

        ! --- run control ---
        call domain_ctl_load(dom%ctl, path_par)

        ! --- ice sheet (grid read from file) ---
        call yelmo_init(dom%yelmo, filename=path_par, grid_def="file", time=time)
        domain = trim(dom%yelmo%par%domain)
        dom%ctl%domain     = trim(domain)
        dom%ctl%grid_yelmo = trim(dom%yelmo%par%grid_name)

        ! --- external forcing models (all on the Yelmo grid for now) ---
        call bsl_init(dom%bsl, path_par, time_rel)

        call isos_init(dom%isos, path_par, "isos", dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                       dom%yelmo%grd%dx, dom%yelmo%grd%dy)

        call snapclim_init(dom%snp, path_par, domain, dom%yelmo%par%grid_name, &
                           dom%yelmo%grd%nx, dom%yelmo%grd%ny, dom%yelmo%bnd%basins)

        call smbpal_init(dom%smb, path_par, x=dom%yelmo%grd%xc, y=dom%yelmo%grd%yc, &
                         lats=dom%yelmo%grd%lat)

        call marshelf_init(dom%mshlf, path_par, "marine_shelf", dom%yelmo%grd%nx, dom%yelmo%grd%ny, &
                           domain, dom%yelmo%par%grid_name, dom%yelmo%bnd%regions, dom%yelmo%bnd%basins)

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

        ! Grids resolve from maps/grid_<name>.txt; prime the maps we will use.
        call coupler_init(dom%cpl)
        call coupler_prime(dom%cpl, dom%ctl%grid_yelmo, dom%ctl%grid_name, "bilin")  ! downscale
        call coupler_prime(dom%cpl, dom%ctl%grid_name, dom%ctl%grid_yelmo, "con")    ! aggregate

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

        ! Marine shelf
        call marshelf_update_shelf(dom%mshlf, dom%yelmo%tpo%now%H_ice, dom%yelmo%bnd%z_bed, &
                dom%yelmo%tpo%now%f_grnd, dom%yelmo%bnd%basins, dom%yelmo%bnd%z_sl, dom%yelmo%grd%dx, &
                dom%snp%now%depth, dom%snp%now%to_ann, dom%snp%now%so_ann, &
                dto_ann=dom%snp%now%to_ann - dom%snp%clim0%to_ann)
        call marshelf_update(dom%mshlf, dom%yelmo%tpo%now%H_ice, dom%yelmo%bnd%z_bed, &
                dom%yelmo%tpo%now%f_grnd, dom%yelmo%bnd%regions, dom%yelmo%bnd%basins, &
                dom%yelmo%bnd%z_sl, dx=dom%yelmo%grd%dx)
        dom%yelmo%bnd%bmb_shlf = dom%mshlf%now%bmb_shlf
        dom%yelmo%bnd%T_shlf   = dom%mshlf%now%T_shlf

        ! Initialize state variables (dyn, therm, mat) with a cold base
        call yelmo_print_bound(dom%yelmo%bnd)
        call yelmo_init_state(dom%yelmo, time=ts%time, thrm_method="robin-cold")

    end subroutine domain_init_state

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
    end subroutine domain_ctl_load

    subroutine yelmox_step(dom, ts)
        ! Advance one domain by one coupling step. Component order is fixed here;
        ! each step_* is a no-op when its ctl flag is off.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        call step_isostasy(dom, ts)
        call step_icesheet(dom, ts)
        call step_climate(dom, ts)
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

    subroutine step_marine_shelf(dom, ts)
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        if (.not. dom%ctl%with_marine_shelf) return

        call marshelf_update_shelf(dom%mshlf, dom%yelmo%tpo%now%H_ice, dom%yelmo%bnd%z_bed, &
                dom%yelmo%tpo%now%f_grnd, dom%yelmo%bnd%basins, dom%yelmo%bnd%z_sl, dom%yelmo%grd%dx, &
                dom%snp%now%depth, dom%snp%now%to_ann, dom%snp%now%so_ann, &
                dto_ann=dom%snp%now%to_ann - dom%snp%clim0%to_ann)
        call marshelf_update(dom%mshlf, dom%yelmo%tpo%now%H_ice, dom%yelmo%bnd%z_bed, &
                dom%yelmo%tpo%now%f_grnd, dom%yelmo%bnd%regions, dom%yelmo%bnd%basins, &
                dom%yelmo%bnd%z_sl, dx=dom%yelmo%grd%dx)
        dom%yelmo%bnd%bmb_shlf = dom%mshlf%now%bmb_shlf
        dom%yelmo%bnd%T_shlf   = dom%mshlf%now%T_shlf
    end subroutine step_marine_shelf

end module yelmox_domain
