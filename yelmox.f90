

program yelmox

    use nml 
    use ncio 
    use timestepping
    use timer
    use timeout 
    use yelmo 
    use yelmo_tools, only : smooth_gauss_2D
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
    
    character(len=256) :: outfldr, file2D, file2D_small, domain
    character(len=256) :: file_isos, file_bsl
    character(len=512) :: path_par
    character(len=512) :: path_lgm
    real(wp) :: dT_now
    real(wp) :: dtt_now
    
    logical,  allocatable :: tmp_mask(:,:)

    type(timeout_class) :: tm_1D, tm_2D, tm_2Dsm

    ! Model timing
    type(timer_class)  :: tmr
    type(timer_class)  :: tmrs
    character(len=512) :: tmr_file 

    type ctrl_params
        character(len=56) :: tstep_method
        real(wp) :: tstep_const
        real(wp) :: time_init
        real(wp) :: time_end
        real(wp) :: time_equil      ! Only for spinup
        real(wp) :: dtt
        real(wp) :: dt_restart
        real(wp) :: dt_clim

        logical  :: with_ice_sheet 
        character(len=56) :: equil_method
        logical  :: write_clim_and_kill

    end type 

    type negis_params
        logical  :: use_negis_par
        real(wp) :: cf_0
        real(wp) :: cf_1
        real(wp) :: cf_centre
        real(wp) :: cf_north
        real(wp) :: cf_south

        real(wp) :: cf_x
        
    end type

    type(ctrl_params)    :: ctl
    type(ice_opt_params) :: opt  
    type(negis_params)   :: ngs 

    ! Internal parameters
    logical  :: running_greenland
    logical  :: greenland_init_marine_H
    logical  :: scale_glacial_smb

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","tstep_method",   ctl%tstep_method)       ! Calendar choice ("const" or "rel")
    call nml_read(path_par,"ctrl","tstep_const",    ctl%tstep_const)        ! Assumed time bp for const method
    call nml_read(path_par,"ctrl","time_init",      ctl%time_init)          ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",       ctl%time_end)           ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",     ctl%time_equil)         ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",            ctl%dtt)                ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt_restart",     ctl%dt_restart)
    call nml_read(path_par,"ctrl","with_ice_sheet", ctl%with_ice_sheet)     ! Include an active ice sheet 
    call nml_read(path_par,"ctrl","equil_method",   ctl%equil_method)       ! What method should be used for spin-up?
    call nml_read(path_par,"ctrl","write_clim_and_kill",   ctl%write_clim_and_kill)       ! Write climate output and kill program?

    ! Get output times
    call timeout_init(tm_1D,  path_par,"tm_1D",  "small",  ctl%time_init,ctl%time_end)
    call timeout_init(tm_2D,  path_par,"tm_2D",  "heavy",  ctl%time_init,ctl%time_end)
    call timeout_init(tm_2Dsm,path_par,"tm_2Dsm","medium", ctl%time_init,ctl%time_end)
    
    ! Hard-coded for now:
    ctl%dt_clim = 10.0      ! [yrs] Frequency to update snapclim snapshot
    
    ! Start timing
    call timer_step(tmr,comp=-1) 
    
    ! === Initialize timestepping ===
    
    call tstep_init(ts,ctl%time_init,ctl%time_end,method=ctl%tstep_method,units="year", &
                                            time_ref=1950.0_wp,const_rel=ctl%tstep_const)

    ! Consistency checks ===

    if (trim(ctl%equil_method) .eq. "opt") then 
        ! Load optimization parameters 

        ! Initially set to zero 
        opt%tf_basins = 0 

        call optimize_par_load(opt,path_par,"opt")

    end if 

    ! Assume program is running from the output folder
    outfldr = "./"
    
    ! Define input and output locations
    file2D              = trim(outfldr)//"yelmo2D.nc"
    file2D_small        = trim(outfldr)//"yelmo2Dsm.nc"

    file_isos           = trim(outfldr)//"fastisostasy.nc"
    file_bsl            = trim(outfldr)//"bsl.nc"

    tmr_file            = trim(outfldr)//"timer_table.txt"

    ! Print summary of run settings 
    write(*,*)
    write(*,*) "timestepping:   ",  trim(ts%method)
    write(*,*) "with_ice_sheet: ",  ctl%with_ice_sheet
    write(*,*) "equil_method:   ",  trim(ctl%equil_method)
    write(*,*)
    write(*,*) "time_init:  ",      ctl%time_init 
    write(*,*) "time_end:   ",      ctl%time_end 
    write(*,*) "dtt:        ",      ctl%dtt 
    write(*,*) "dt_restart: ",      ctl%dt_restart 
    write(*,*) 
    
    if (trim(ts%method) .eq. "const") then 
        write(*,*) "time_equil: ",    ctl%time_equil
    end if 

    write(*,*) "time    = ", ts%time 
    write(*,*) "time_bp = ", ts%time_rel 
    write(*,*) 
    

    ! === Initialize ice sheet model =====

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=ts%time)

    ! Store domain name as a shortcut 
    domain = yelmo1%par%domain 

    ! Ensure optimization fields are allocated and preassigned
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))
    
    opt%cf_min = yelmo1%dyn%par%till_cf_min
    opt%cf_max = yelmo1%dyn%par%till_cf_ref

    ! Define specific regions of interest =====================

    allocate(tmp_mask(yelmo1%grd%nx,yelmo1%grd%ny))
    
    select case(trim(domain))

        case("Antarctica")

            ! Initialize regions
            call yelmo_regions_init(yelmo1,n=3)

            ! APIS
            call get_ice_sub_region(tmp_mask,"APIS",yelmo1%par%domain,yelmo1%par%grid_name)
            call yelmo_region_init(yelmo1%regs(1),"APIS",mask=tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

            ! WAIS
            call get_ice_sub_region(tmp_mask,"WAIS",yelmo1%par%domain,yelmo1%par%grid_name)
            call yelmo_region_init(yelmo1%regs(2),"WAIS",mask=tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

            ! EAIS
            call get_ice_sub_region(tmp_mask,"EAIS",yelmo1%par%domain,yelmo1%par%grid_name)
            call yelmo_region_init(yelmo1%regs(3),"EAIS",mask=tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

        case("Laurentide")

            ! Make sure to set mask_ice to prevent ice from growing in
            ! Greenland (and on grid borders)

            where(abs(yelmo1%bnd%regions - 1.30) .lt. 1e-3) yelmo1%bnd%mask_ice = -1

            yelmo1%bnd%mask_ice(1,:)             = -1
            yelmo1%bnd%mask_ice(yelmo1%grd%nx,:) = -1
            yelmo1%bnd%mask_ice(:,1)             = -1
            yelmo1%bnd%mask_ice(:,yelmo1%grd%ny) = -1
            
            ! Initialize regions
            call yelmo_regions_init(yelmo1,n=1)

            ! Hudson region
            call get_ice_sub_region(tmp_mask,"Hudson",yelmo1%par%domain,yelmo1%par%grid_name)
            call yelmo_region_init(yelmo1%regs(1),"Hudson",mask=tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

        case("Greenland")

            running_greenland = .TRUE.

            ! Should extra ice be imposed over continental shelf to mimic LGM state to start
            greenland_init_marine_H = .TRUE. 
            
            ! Should glacial smb be modified to reduce negative smb values
            scale_glacial_smb = .FALSE. 
            
            ! Should NEGIS parameter modifications be used
            ngs%use_negis_par = .TRUE. 

if (.FALSE.) then
    ! ajr, 2025-01-15, missing these parameters in param files - ask Ilaria!!

            ! Load NEGIS parameters from file, if used
            if (ngs%use_negis_par) then
                
                call nml_read(path_par,"negis","cf_0",       ngs%cf_0)
                call nml_read(path_par,"negis","cf_1",       ngs%cf_1)
                call nml_read(path_par,"negis","cf_centre",  ngs%cf_centre)
                call nml_read(path_par,"negis","cf_north",   ngs%cf_north)
                call nml_read(path_par,"negis","cf_south",   ngs%cf_south)
                
            end if 
end if

            ! Make sure to set mask_ice to prevent ice from growing in
            ! Iceland and Svaalbard (on grid borders)

            where(abs(yelmo1%bnd%regions - 1.20) .lt. 1e-3) yelmo1%bnd%mask_ice = -1
            where(abs(yelmo1%bnd%regions - 1.23) .lt. 1e-3) yelmo1%bnd%mask_ice = -1
            where(abs(yelmo1%bnd%regions - 1.31) .lt. 1e-3) yelmo1%bnd%mask_ice = -1
            
            if (yelmo1%dyn%par%till_method .eq. -1) then 
                ! Initialize cb_ref to constant value to start
                ! (will be updated in the time loop)

                yelmo1%dyn%now%cb_ref = yelmo1%dyn%par%till_cf_ref 

            end if 

        case("Patagonia")

            if (.FALSE.) then
                ! Fix ice thickness to obs everywhere outside of masked region

                    ! Dynamically evolve ice where region mask == 1
                    where(abs(yelmo1%bnd%regions - 1.0) .lt. 1e-3)
                        yelmo1%bnd%mask_ice = 1
                    elsewhere
                        yelmo1%bnd%mask_ice = 0
                    end where

            else
                ! Fix ice thickness to obs on the margins, relax to obs outside of masked region

                ! Evolve ice everywhere except the domain borders
                yelmo1%bnd%mask_ice                  = 1
                yelmo1%bnd%mask_ice(1,:)             = 0
                yelmo1%bnd%mask_ice(yelmo1%grd%nx,:) = 0
                yelmo1%bnd%mask_ice(:,1)             = 0
                yelmo1%bnd%mask_ice(:,yelmo1%grd%ny) = 0

                where(abs(yelmo1%bnd%regions - 1.0) .lt. 1e-3)
                    yelmo1%bnd%tau_relax = -1.0      ! inside icefield: free evolution
                elsewhere
                    yelmo1%bnd%tau_relax = 50.0      ! outside: Relaxation timescale [yrs], relax toward H_ice_ref
                end where
            
            end if

    end select

    ! === Initialize external models (forcing for ice sheet) ======

    ! jablasco
    if (.TRUE.) then
        ! Initialize barysealevel model
        call bsl_init(bsl, path_par, ts%time_rel)

        ! Initialize fastisosaty
        call isos_init(isos1, path_par, "isos", yelmo1%grd%nx, yelmo1%grd%ny, yelmo1%grd%dx, yelmo1%grd%dy)
    end if
    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)
    
    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    yelmo1%bnd%H_sed = sed1%now%H 
    
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    yelmo1%bnd%Q_geo = gthrm1%now%ghf 
    
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

    ! Barystatic sea level
    ! jablasco
    if (.TRUE.) then
        call bsl_update(bsl, year_bp=ts%time_rel)
        call bsl_write_init(bsl, file_bsl, ts%time)

        ! Initialize the isostasy reference state using reference topography fields
        call isos_init_ref(isos1, yelmo1%bnd%z_bed_ref, yelmo1%bnd%H_ice_ref)
        call isos_init_state(isos1, yelmo1%bnd%z_bed, yelmo1%tpo%now%H_ice, ts%time, bsl)
        call isos_write_init_extended(isos1, file_isos, ts%time)

        yelmo1%bnd%z_bed = isos1%out%z_bed
        yelmo1%bnd%z_sl  = isos1%out%z_ss
    end if

    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=ts%time_rel,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,ts%time_rel,time_equil=100.0)
    end if 

    call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,ts%time_rel) 
    yelmo1%bnd%smb   = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

    ! Write out initial boundary conditions if desired
    if (ctl%write_clim_and_kill) then
        call snapclim_write_init(snp1,"yelmox_climate.nc",yelmo1%grd%xc,yelmo1%grd%yc,time=ts%time,units="years")
        call snapclim_write_step(snp1,"yelmox_climate.nc",time=ts%time)
        stop
    end if

    if (trim(yelmo1%par%domain) .eq. "Greenland" .and. scale_glacial_smb) then 
        ! Modify glacial smb
        call calc_glacial_smb(yelmo1%bnd%smb,yelmo1%grd%lat,snp1%now%ta_ann,snp1%clim0%ta_ann)
    end if

    call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%bmb_shlf = -0.5
    yelmo1%bnd%T_shlf   = 270.0

    call yelmo_print_bound(yelmo1%bnd) 
    
    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=ts%time,thrm_method="robin-cold")

    ! ===== basal friction optimization ======
    if (trim(ctl%equil_method) .eq. "opt") then 
        
        ! Ensure that cb_ref will be optimized (till_method == set externally) 
        yelmo1%dyn%par%till_method = -1  

        ! If not using restart, prescribe cb_ref to initial guess 
        if (.not. yelmo1%par%use_restart) then
            yelmo1%dyn%now%cb_ref = opt%cf_init 
        end if 

    end if 
    ! ========================================

    if (.not. yelmo1%par%use_restart) then 
        ! No restart file used, perform various initialization steps

        select case(trim(domain))

        case("Laurentide")
            ! Special start-up steps for Laurentide

            if (trim(ctl%tstep_method) .eq. "const") then
                ! Steady-state simulation, start with lgm state
                call yelmox_init_laurentide_lgm(yelmo1,snp1,smbpal1,ts,path_par,method="ref_lgm", &
                                                        with_ice_sheet=ctl%with_ice_sheet)
            else
                ! Transient simulation - start with no ice thickness
                call yelmox_init_laurentide_lgm(yelmo1,snp1,smbpal1,ts,path_par,method="zero", &
                                                        with_ice_sheet=ctl%with_ice_sheet)
            end if

        case("Greenland")
            ! Special start-up steps for Greenland

            if (ngs%use_negis_par) then 

                ! Ensure till method is correct, since we are updating cb_ref externally
                yelmo1%dyn%par%till_method = -1 

                ! Update cb_ref using negis parameters 
                call negis_update_cb_ref(yelmo1,ngs,ts%time)

            end if 

            if (greenland_init_marine_H) then
                ! Add extra ice-thickness over continental shelf to start with
                ! an LGM-like state

                where(yelmo1%bnd%mask_ice /= -1 .and. yelmo1%tpo%now%H_ice .lt. 600.0 &
                        .and. yelmo1%bnd%z_bed .gt. -500.0)

                        yelmo1%tpo%now%H_ice = 800.0 

                end where

                if (ctl%with_ice_sheet) then
                    ! Run yelmo for a few years with constant boundary conditions
                    ! to synchronize all model fields a bit
                    call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_wp, dt=1.0_wp,topo_fixed=.FALSE.)
                end if 

            end if
                
        case DEFAULT 
            ! Run simple startup equilibration step 
            
            if (ctl%with_ice_sheet) then
                ! Run yelmo for a few years with constant boundary conditions
                ! to synchronize all model fields a bit
                call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_wp, dt=1.0_wp,topo_fixed=.FALSE.)
            end if 

        end select
        
    end if 

    ! ===== Initialize output files ===== 
    
    if (tm_2D%active) then
        call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years")
    end if
    if (tm_2Dsm%active) then
        call yelmo_write_init(yelmo1,file2D_small,time_init=ts%time,units="years")
    end if
    
    call yelmo_regions_write(yelmo1,ts%time,init=.TRUE.,units="years")

    call timer_step(tmr,comp=1,label="initialization") 
    call timer_step(tmrs,comp=-1)
    
    ! ==== Begin main time loop =====

    dtt_now = ctl%dtt
    call tstep_print_header(ts)

    do while (.not. ts%is_finished)

        ! == Update timestep ===

        call tstep_update(ts,dtt_now)
        call tstep_print(ts)
        
        ! Spin-up procedure - only relevant for time_elapsed <= time_equil
        select case(trim(ctl%equil_method))
            
            case("opt")

                if (ts%time_elapsed .le. opt%rel_time2) then 
                    ! Apply relaxation to the model 

                    ! Update model relaxation time scale and error scaling (in [m])
                    call optimize_set_transient_param(opt%rel_tau,ts%time_elapsed,time1=opt%rel_time1,time2=opt%rel_time2, &
                                                    p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)
                    
                    ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                    yelmo1%tpo%par%topo_rel_tau = opt%rel_tau 
                    yelmo1%tpo%par%topo_rel     = 4
                
                else 
                    ! Turn-off relaxation now

                    yelmo1%tpo%par%topo_rel = 0 

                end if 

                ! === Optimization update step =========

                if (opt%opt_cf .and. &
                        (ts%time_elapsed .ge. opt%cf_time_init .and. ts%time_elapsed .le. opt%cf_time_end) ) then
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                                    yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                                    yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                                    opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                                    dt=ctl%dtt,fill_method=opt%fill_method,fill_dist=opt%sigma_err,cb_tgt=yelmo1%dyn%now%cb_tgt)
                    
                end if

                if (opt%opt_tf .and. &
                        (ts%time_elapsed .ge. opt%tf_time_init .and. ts%time_elapsed .le. opt%tf_time_end) ) then
                    ! Perform tf_corr optimization

                    call optimize_tf_corr(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                                                        yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,opt%H_grnd_lim,yelmo1%bnd%basins, &
                                                        opt%basin_fill,opt%tau_m,opt%m_temp,opt%tf_min,opt%tf_max,yelmo1%tpo%par%dx,sigma=opt%tf_sigma,dt=ctl%dtt)
                
                end if 

            case("relax")
                ! ===== relaxation spinup ==================

                if (ts%time_elapsed .lt. ctl%time_equil) then 
                    ! Turn on relaxation for now, to let thermodynamics equilibrate
                    ! without changing the topography too much. Important when 
                    ! effective pressure = f(thermodynamics).

                    yelmo1%tpo%par%topo_rel     = 3
                    yelmo1%tpo%par%topo_rel_tau = 50.0 
                    write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau

                else if (ts%time_elapsed .eq. ctl%time_equil) then 
                    ! Disable relaxation now... 

                    yelmo1%tpo%par%topo_rel     = 0
                    write(*,*) "timelog, relation off..."

                end if 

            case DEFAULT   ! == "none", etc

                ! Pass - do nothing 

        end select 

        call timer_step(tmrs,comp=0) 
        
        ! == ISOSTASY and SEA LEVEL ======================================================
        ! jablasco
        if (.FALSE.) then
                call bsl_update(bsl, ts%time_rel)
                call isos_update(isos1, yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
                yelmo1%bnd%z_bed = isos1%out%z_bed
                yelmo1%bnd%z_sl  = isos1%out%z_ss
        end if

        call timer_step(tmrs,comp=1,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="isostasy") 
        
        ! == ICE SHEET ===================================================

        if (running_greenland .and. ngs%use_negis_par) then 
            ! Update cb_ref using negis parameters 

            call negis_update_cb_ref(yelmo1,ngs,ts%time)

        end if
        
        ! Update Yelmo
        if (ctl%with_ice_sheet .and. (.not. (ts%n .eq. 0 .and. yelmo1%par%use_restart)) ) then
            call yelmo_update(yelmo1,ts%time)
        end if
        
        call timer_step(tmrs,comp=2,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="yelmo")
        
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        
        if (mod(nint(ts%time_elapsed*100),nint(ctl%dt_clim*100))==0) then
                ! Update snapclim
                call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=ts%time,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins) 
        end if 

        ! == SURFACE MASS BALANCE ==============================================

        call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                   yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,ts%time_rel) 
        yelmo1%bnd%smb   = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

        if (trim(yelmo1%par%domain) .eq. "Greenland" .and. scale_glacial_smb) then 
            ! Modify glacial smb
            call calc_glacial_smb(yelmo1%bnd%smb,yelmo1%grd%lat,snp1%now%ta_ann,snp1%clim0%ta_ann)
        end if

        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        if (trim(yelmo1%par%domain) .eq. "Pyrenees") then
                yelmo1%bnd%bmb_shlf = 0.0_wp
                yelmo1%bnd%T_shlf   = 0.0_wp
        else if (trim(yelmo1%par%domain) .eq. "Patagonia") then
                ! For now, set ocean to zero!!!
                yelmo1%bnd%bmb_shlf = 0.0_wp
                yelmo1%bnd%T_shlf   = 0.0_wp
        else
                call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                                snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

                call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                     yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

                yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
                yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf
        end if
        
        call timer_step(tmrs,comp=3,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="climate") 

        ! == MODEL OUTPUT =======================================================

        if (timeout_check(tm_2D,ts%time)) then
            call yelmox_write_step(yelmo1,snp1,mshlf1,smbpal1,file2D,ts%time)
        end if

        if (timeout_check(tm_2Dsm,ts%time)) then 
            call yelmo_write_step(yelmo1,file2D_small,ts%time,compare_pd=.FALSE.)
        end if

        if (timeout_check(tm_1D,ts%time)) then 
            call yelmo_regions_write(yelmo1,ts%time)
        end if 

        if (mod(nint(ts%time*100),nint(ctl%dt_restart*100))==0) then
            call yelmox_restart_write(bsl,isos1,yelmo1,mshlf1,ts%time)
        end if 
 
        call timer_step(tmrs,comp=4,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="io") 
        
        if (mod(ts%time_elapsed,10.0)==0) then
            ! Print timestep timing info and write log table
            call timer_write_table(tmrs,[ts%time,dtt_now]*1e-3,"m",tmr_file,init=ts%time_elapsed .eq. 0.0)
        end if 

        if (mod(ts%time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", ts%time
        end if 
        
    end do 

    ! Stop timing
    call timer_step(tmr,comp=2,time_mod=[ctl%time_init,ts%time]*1e-3,label="timeloop") 
    
    ! Write the restart snapshot for the end of the simulation
    call yelmox_restart_write(bsl,isos1,yelmo1,mshlf1,ts%time)

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=ts%time*1e-3)
    
contains
    
    subroutine yelmox_init_laurentide_lgm(ylmo,snp,smb,ts,path_par,method,with_ice_sheet)

        implicit none

        type(yelmo_class),      intent(INOUT) :: ylmo
        type(snapclim_class),   intent(INOUT) :: snp
        type(smbpal_class),     intent(INOUT) :: smb
        type(tstep_class),      intent(IN)    :: ts
        character(len=*),       intent(IN)    :: path_par
        character(len=*),       intent(IN)    :: method
        logical,                intent(IN)    :: with_ice_sheet

        ! Local variables
        character(len=1024) :: path_lgm
        type(yelmo_class) :: ylmo_eq

        ! Load LGM reconstruction into reference ice thickness
        path_lgm = "ice_data/Laurentide/"//trim(ylmo%par%grid_name)//&
                    "/"//trim(ylmo%par%grid_name)//"_TOPO-ICE-6G_C.nc"
        call nc_read(path_lgm,"dz",ylmo%bnd%H_ice_ref,start=[1,1,1], &
                            count=[ylmo%tpo%par%nx,ylmo%tpo%par%ny,1]) 

        ! Determine initial ice thickness
        select case(trim(method))

        case("const")
            ! Start with some ice cover to speed up initialization

            ylmo%tpo%now%H_ice = 0.0
            where (ylmo%bnd%regions .eq. 1.1 .and. ylmo%bnd%z_bed .gt. 0.0) ylmo%tpo%now%H_ice = 1000.0 
            where (ylmo%bnd%regions .eq. 1.12) ylmo%tpo%now%H_ice = 1000.0 

            ! Apply Gaussian smoothing to keep things stable
            call smooth_gauss_2D(ylmo%tpo%now%H_ice,dx=ylmo%grd%dx,f_sigma=3.0)

            ! Make sure to update topographic info (without loading anything)
            call yelmo_init_topo(ylmo,path_par,ylmo%par%nml_init_topo,ts%time,load_topo=.FALSE.)

        case("ref_lgm")
            ! Set LGM reconstruction as initial ice thickness over North America
            
            where ( ylmo%bnd%z_bed .gt. -500.0 .and. &
                    (   ylmo%bnd%regions .eq. 1.1  .or. &
                        ylmo%bnd%regions .eq. 1.11 .or. &
                        ylmo%bnd%regions .eq. 1.12) )
                ylmo%tpo%now%H_ice = ylmo%bnd%H_ice_ref
            end where 

            ! Apply Gaussian smoothing to keep things stable
            call smooth_gauss_2D(ylmo%tpo%now%H_ice,dx=ylmo%grd%dx,f_sigma=2.0)

            ! Make sure to update topographic info (without loading anything)
            call yelmo_init_topo(ylmo,path_par,ylmo%par%nml_init_topo,ts%time,load_topo=.FALSE.)

        case DEFAULT
            ! Zero ice thickness

            ! Pass - do nothing

        end select 
        
        ! Run Yelmo for briefly to update surface topography fields (but fixed H)
        call yelmo_update_equil(ylmo,ts%time,time_tot=1.0_wp,dt=1.0,topo_fixed=.TRUE.)

        ! Addtional cleanup - remove thin floating ice 
        where( ylmo%tpo%now%mask_bed .eq. 5 .and. ylmo%tpo%now%H_ice .lt. 50.0_wp) ylmo%tpo%now%H_ice = 0.0 
        call yelmo_update_equil(ylmo,ts%time,time_tot=1.0_wp,dt=1.0,topo_fixed=.TRUE.)

        if (trim(method) .eq. "ref_lgm") then
            ! Store new "clean" ice thickness as reference state
            ylmo%bnd%H_ice_ref = ylmo%tpo%now%H_ice
        end if 

        ! Update snapclim to reflect new topography 
        call snapclim_update(snp,z_srf=ylmo%tpo%now%z_srf,time=ts%time,domain=domain,dx=ylmo%grd%dx,basins=ylmo%bnd%basins)

        ! Update smbpal
        call smbpal_update_monthly(smb,snp%now%tas,snp%now%pr, &
                                    ylmo%tpo%now%z_srf,ylmo%tpo%now%H_ice,ts%time_rel) 
        ylmo%bnd%smb   = smb%ann%smb*ylmo%bnd%c%conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
        ylmo%bnd%T_srf = smb%ann%tsrf 

        if (trim(method) .eq. "const") then
            ! Additionally ensure smb is postive for land above 50degN in Laurentide region
            ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
            where (ylmo%bnd%regions .eq. 1.1 .and. ylmo%grd%lat .gt. 50.0 .and. &
                    ylmo%bnd%z_bed .gt. 0.0 .and. ylmo%bnd%smb .lt. 0.0 ) ylmo%bnd%smb = 0.5 
            
            if (with_ice_sheet) then
                ! Run yelmo for several years to ensure stable central ice dome
                call yelmo_update_equil(ylmo,ts%time,time_tot=5e3,dt=5.0,topo_fixed=.FALSE.)
            end if 

        else ! method=="ref_lgm"

            if (with_ice_sheet) then

                ! Set relaxation for equilibrium
                ylmo_eq = ylmo
                !ylmo_eq%tpo%par%calv_flt_method = "zero"
                ! ylmo_eq%tpo%par%topo_rel = 3
                ylmo_eq%dyn%par%beta_min = 100.0

                ! Run yelmo for several years with constant boundary conditions to stabilize fields
                ! ylmo_eq%tpo%par%topo_rel_tau = 10.0
                ! call yelmo_update_equil(ylmo_eq,ts%time,time_tot=1e2,dt=5.0,topo_fixed=.FALSE.)
                !ylmo_eq%tpo%par%topo_rel_tau = 100.0
                !call yelmo_update_equil(ylmo_eq,ts%time,time_tot=1e2,dt=5.0,topo_fixed=.FALSE.)
                call yelmo_update_equil(ylmo_eq,ts%time,time_tot=2e2,dt=5.0,topo_fixed=.FALSE.)

                ! Restore parameters and state back to original yelmo object
                !ylmo_eq%tpo%par%calv_flt_method = ylmo%tpo%par%calv_flt_method
                !ylmo_eq%tpo%par%topo_rel = ylmo%tpo%par%topo_rel 
                !ylmo_eq%tpo%par%topo_rel_tau = ylmo%tpo%par%topo_rel_tau
                ylmo_eq%dyn%par%beta_min = ylmo%dyn%par%beta_min

                ylmo = ylmo_eq

            end if 

        end if 

        return

    end subroutine yelmox_init_laurentide_lgm

    subroutine calc_glacial_smb(smb,lat2D,ta_ann,ta_ann_pd)

        implicit none

        real(wp), intent(INOUT) :: smb(:,:)
        real(wp), intent(IN)    :: lat2D(:,:)
        real(wp), intent(IN)    :: ta_ann(:,:)
        real(wp), intent(IN)    :: ta_ann_pd(:,:)

        ! Local variables 
        integer  :: i, j, nx, ny 
        real(wp) :: t0, tnow
        real(wp) :: at
        real(wp) :: fac

        real(wp), parameter :: dt_lgm  = -8.0 
        real(wp), parameter :: lat_lim = 55.0 
        real(wp), parameter :: fac_lim = 0.90

        nx = size(smb,1)
        ny = size(smb,2) 

        ! Determine a quasi glacial-interglacial index
        ! 0: interglacial
        ! 1: glacial 
        tnow = sum(ta_ann) / real(nx*ny,wp)
        t0   = sum(ta_ann_pd) / real(nx*ny,wp)

        at = (tnow-t0)/dt_lgm
        if (at .lt. 0.0) at = 0.0
        if (at .gt. 1.0) at = 1.0
        
        ! Now determine smb scaling as a function of glacial index and latitude
        ! fac==0: no change to smb
        ! 0<fac<=1: scale smb up to value of fac_lim

        do j = 1, ny 
        do i = 1, nx

            if (smb(i,j) .lt. 0.0 .and. lat2D(i,j) .gt. lat_lim) then
                ! Calculate scaling here

                smb(i,j) = smb(i,j) -  smb(i,j) * at * fac_lim

            end if

        end do
        end do

        return

    end subroutine calc_glacial_smb

    ! ===== NEGIS routines ============================

    subroutine negis_update_cb_ref(ylmo,ngs,time)

        implicit none

        type(yelmo_class),  intent(INOUT) :: ylmo
        type(negis_params), intent(INOUT) :: ngs 
        real(wp), intent(IN) :: time 

        ! Local variables
        integer :: i, j, nx, ny 

        nx = ylmo%grd%nx 
        ny = ylmo%grd%ny 

        if (time .lt. -11e3) then 
            ngs%cf_x = ngs%cf_0
        else
            ! Linear interpolation from cf_0 to cf_1  
            ngs%cf_x = ngs%cf_0 + (time - (-11e3))/ ((0.0) - (-11e3)) * (ngs%cf_1 - ngs%cf_0)
        end if

        if (time .lt. -4e3) then
            ngs%cf_south = 1.0
        else
            ngs%cf_north = 1.0
        end if

        ! Update bed roughness coefficients cb_ref and c_bed (which are independent of velocity)
        ! like normal, using the default function defined in Yelmo:
        call calc_cb_ref(ylmo%dyn%now%cb_ref,ylmo%bnd%z_bed,ylmo%bnd%z_bed_sd,ylmo%bnd%z_sl, &
                ylmo%bnd%H_sed,ylmo%dyn%par%till_f_sed,ylmo%dyn%par%till_sed_min,ylmo%dyn%par%till_sed_max, &
                ylmo%dyn%par%till_cf_ref,ylmo%dyn%par%till_cf_min,ylmo%dyn%par%till_z0,ylmo%dyn%par%till_z1, &
                ylmo%dyn%par%till_n_sd,ylmo%dyn%par%till_scale_zb,ylmo%dyn%par%till_scale_sed)

        ! === Finally, apply NEGIS scaling =============================

        do j = 1, ny 
        do i = 1, nx 

            if(ylmo%bnd%basins(i,j) .eq. 9.1) yelmo1%dyn%now%cb_ref(i,j) = yelmo1%dyn%now%cb_ref(i,j) * ngs%cf_centre
            if(ylmo%bnd%basins(i,j) .eq. 9.2) yelmo1%dyn%now%cb_ref(i,j) = yelmo1%dyn%now%cb_ref(i,j) * ngs%cf_south
            if(ylmo%bnd%basins(i,j) .eq. 9.3) yelmo1%dyn%now%cb_ref(i,j) = yelmo1%dyn%now%cb_ref(i,j) * ngs%cf_north

        end do
        end do

        return

    end subroutine negis_update_cb_ref

    subroutine yelmox_write_step(ylmo,snp,mshlf,srf,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf 
        type(smbpal_class),     intent(IN) :: srf 
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
                
        ! == yelmo_topography ==
        call yelmo_write_var(filename,"H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"mask_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"tau_relax",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_net",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_resid",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"fmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"N_eff",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"dHidt",ylmo,n,ncid)
        call yelmo_write_var(filename,"lsf",ylmo,n,ncid)

        call yelmo_write_var(filename,"cmb_flt",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb_grnd",ylmo,n,ncid)

        !call yelmo_write_var(filename,"eps_eff",ylmo,n,ncid)
        !call yelmo_write_var(filename,"tau_eff",ylmo,n,ncid)

        ! == yelmo_dynamics ==
        call yelmo_write_var(filename,"cb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"c_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_eff_int",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_s",ylmo,n,ncid)
        
        ! == yelmo_material ==
        call yelmo_write_var(filename,"enh_bar",ylmo,n,ncid)
        !call yelmo_write_var(filename,"ATT",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_int",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strn_de",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strn_te",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strn2D_de",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strn2D_div",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strs2D_te",ylmo,n,ncid)

        ! == yelmo_thermodynamics ==
        call yelmo_write_var(filename,"T_prime",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_pmp",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_w",ylmo,n,ncid)
        
        !call yelmo_write_var(filename,"Q_strn",ylmo,n,ncid)
        ! Write Q_strn in K/yr instead of mW/m^3:
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call yelmo_write_var(filename,"Q_ice_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_b",ylmo,n,ncid)
        
        ! == yelmo_boundaries ==
        call yelmo_write_var(filename,"z_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_sl",ylmo,n,ncid)
        !call yelmo_write_var(filename,"smb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"T_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_shlf",ylmo,n,ncid)
        !call yelmo_write_var(filename,"Q_geo",ylmo,n,ncid)     ! Written as static field below

        if (.FALSE.) then
                ! == yelmo_data (comparison with present-day) ==
                call yelmo_write_var(filename,"pd_err_H_ice",ylmo,n,ncid)
                call yelmo_write_var(filename,"pd_err_z_srf",ylmo,n,ncid)
                call yelmo_write_var(filename,"pd_err_uxy_s",ylmo,n,ncid)
                call yelmo_write_var(filename,"pd_err_smb_ref",ylmo,n,ncid)
        end if

        ! == yelmo extra fields ==
if (.FALSE.) then
        write(*,*) "yelmox extra"
        call yelmo_write_var(filename,"ssa_mask_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"ssa_mask_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"dzsdx",ylmo,n,ncid)
        call yelmo_write_var(filename,"dzsdy",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_s",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_s",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta_acy",ylmo,n,ncid)
end if

        ! Static fields
        if (n .le. 1) then 
            call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
            call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW m^-2",long_name="Geothermal heat flow at depth ", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if

        ! == snapclim ==
        call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        !call nc_write(filename,"pr",snp%now%pr*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
        !              dim1="xc",dim2="yc",dim3="month",dim4="time",start=[1,1,1,n],ncid=ncid)
              
        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == smbpal ==
!         call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == marine_shelf ==

        if (.FALSE.) then
                call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                call nc_write(filename,"mask_ocn",mshlf%now%mask_ocn,units="", &
                                long_name="Ocean mask (0: land, 1: grline, 2: fltline, 3: open ocean, 4: deep ocean, 5: lakes)", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
                if (trim(mshlf%par%bmb_method) .eq. "pico") then 
                        call nc_write(filename,"d_shlf",mshlf%pico%now%d_shlf,units="km",long_name="Shelf distance to grounding line", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"d_if",mshlf%pico%now%d_if,units="km",long_name="Shelf distance to ice front", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"boxes",mshlf%pico%now%boxes,units="",long_name="Shelf boxes", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)      
                        call nc_write(filename,"r_shlf",mshlf%pico%now%r_shlf,units="",long_name="Ratio of ice shelf", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"T_box",mshlf%pico%now%T_box,units="K?",long_name="Temperature of boxes", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"S_box",mshlf%pico%now%S_box,units="PSU",long_name="Salinity of boxes", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"A_box",mshlf%pico%now%A_box*1e-6,units="km2",long_name="Box area of ice shelf", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                else 
                        call nc_write(filename,"tf_basin",mshlf%now%tf_basin,units="K",long_name="Mean basin thermal forcing", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing applied correction factor", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

                        call nc_write(filename,"tf_corr_basin",mshlf%now%tf_corr_basin,units="K",long_name="Shelf thermal forcing basin-wide correction factor", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
                                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                end if 
        end if
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step

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

end program yelmox



