

program yelmox

    use nml
    use ncio
    use timer
    use timestepping
    use timeout
    use yelmo
    use ice_optimization

    ! External libraries
    use fastisostasy    ! also reexports barysealevel
    use rembo_sclimate 
    use snapclim
    use marine_shelf
    use sediments
    use geothermal
    
    use hyster
    
    implicit none 

    type(tstep_class)      :: ts
    
    type(yelmo_class)      :: yelmo1
    type(snapclim_class)   :: snp1
    type(bsl_class)        :: bsl
    type(marshelf_class)   :: mshlf1
    type(sediments_class)  :: sed1
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    type(hyster_class)     :: hyst1

    character(len=256) :: outfldr, file1D, file2D, file2D_small, domain
    character(len=256) :: file1D_hyst, file_isos, file_bsl, file_rembo
    character(len=512) :: path_par
    real(wp) :: time_init, time_end, time_equil, dtt, dt_restart
    real(wp) :: dtt_now, deltat_tot
    character(len=56)  :: tstep_method
    real(wp) :: tstep_const
    integer  :: n

    type(timeout_class) :: tm_1D, tm_2D 
    logical :: file_rembo_write_ocn_forcing

    logical :: use_hyster
    logical :: write_restart
    real(4) :: var, dv_dt
    real(4) :: convert_km3_Gt
    real(4) :: dT_summer, dT_ann, dT_ocn

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  
    logical :: lim_pd_ice 
    logical :: with_ice_sheet 
    logical :: optimize 
    logical :: greenland_init_marine_H

    type(ice_opt_params) :: opt 

    ! Model timing
    type(timer_class)  :: tmr
    type(timer_class)  :: tmrs
    character(len=512) :: tmr_file 

    real(wp) :: hyst_f_to 
    real(wp) :: hyst_f_ta

    ! Assume program is running from the output folder
    outfldr = "./"
    
    ! Determine the parameter file from the command line 
    !call yelmo_load_command_line_args(path_par)
    path_par = trim(outfldr)//"yelmo_Greenland_rembo.nml"
    
    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","time_init",      time_init)              ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",       time_end)               ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",     time_equil)             ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",            dtt)                    ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","write_restart",  write_restart)
    call nml_read(path_par,"ctrl","tstep_method",   tstep_method)           ! Calendar choice ("const" or "rel")
    call nml_read(path_par,"ctrl","tstep_const",    tstep_const)            ! Assumed time bp for const method
    call nml_read(path_par,"ctrl","use_hyster",     use_hyster)             ! Use hyster?
    call nml_read(path_par,"ctrl","dT",             dT_summer)              ! Initial summer temperature anomaly
    call nml_read(path_par,"ctrl","lim_pd_ice",     lim_pd_ice)             ! Limit to pd ice extent (apply extra melting outside mask)
    call nml_read(path_par,"ctrl","with_ice_sheet", with_ice_sheet)         ! Active ice sheet? 
    call nml_read(path_par,"ctrl","greenland_init_marine_H", greenland_init_marine_H)   ! Initialize ice thickness with extra marine ice?
    call nml_read(path_par,"ctrl","optimize",       optimize)               ! Optimize basal friction?
    call nml_read(path_par,"ctrl","write_ocn_forcing", file_rembo_write_ocn_forcing)
    
    call nml_read(path_par,"ctrl","f_to",           hyst_f_to)              ! Scale hyster forcing to temp (ocean)
    call nml_read(path_par,"ctrl","f_ta",           hyst_f_ta)              ! Scale hyster forcing to temp (atm)
    
    ! Get output times
    call timeout_init(tm_1D,path_par,"tm_1D","small", time_init,time_end)
    call timeout_init(tm_2D,path_par,"tm_2D","heavy", time_init,time_end)
    !stop 

    if (optimize) then 
        ! Load optimization parameters 

        call optimize_par_load(opt,path_par,"opt")

    end if 

    ! Define input and output locations
    file1D              = trim(outfldr)//"yelmo1D.nc"
    file2D              = trim(outfldr)//"yelmo2D.nc"
    file2D_small        = trim(outfldr)//"yelmo2Dsm.nc"
    file1D_hyst         = trim(outfldr)//"yelmo1D_hyst.nc" 

    file_isos           = trim(outfldr)//"fastisostasy.nc"
    file_bsl            = trim(outfldr)//"bsl.nc"

    file_rembo          = trim(outfldr)//"yelmo-rembo.nc"

    tmr_file            = trim(outfldr)//"timer_table.txt"

    ! How often to write a restart file 
    dt_restart   = 20e3                 ! [yr] 

    ! Start timing
    call timer_step(tmr,comp=-1) 
    
    ! === Initialize timestepping ===
    
    call tstep_init(ts,time_init,time_end,method=tstep_method,units="year", &
                                            time_ref=1950.0_wp,const_rel=tstep_const)

    ! === Initialize ice sheet model =====

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=ts%time)

    ! === Initialize external models (forcing for ice sheet) ======

    ! Store domain name as a shortcut 
    domain = yelmo1%par%domain 

    ! Ensure optimization fields are allocated and preassigned
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))
    
    opt%cf_min = yelmo1%dyn%par%till_cf_min 
    opt%cf_max = yelmo1%dyn%par%till_cf_ref

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! Initialize barysealevel model
    call bsl_init(bsl, path_par, ts%time_rel)

    ! Initialize fastisosaty
    call isos_init(isos1, path_par, "isos", yelmo1%grd%nx, yelmo1%grd%ny, &
        yelmo1%grd%dx, yelmo1%grd%dy)

    ! Initialize the climate model REMBO, including loading parameters from options_rembo 
    call rembo_init(real(ts%time,8))

    ! Initialize hysteresis module for transient forcing experiments 
    call hyster_init(hyst1,path_par,ts%time) 
    convert_km3_Gt = yelmo1%bnd%c%rho_ice *1e-3

    ! Initialize "climate" model (here for ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny, &
                        domain,yelmo1%par%grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Sediments
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    yelmo1%bnd%H_sed = sed1%now%H 
    
    ! Geothermal heat flow
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    yelmo1%bnd%Q_geo = gthrm1%now%ghf 
    
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    ! Barystatic sea level
    call bsl_update(bsl, year_bp=ts%time_rel)
    call bsl_write_init(bsl, file_bsl, ts%time)

    ! Initialize the isostasy reference state using reference topography fields
    call isos_init_ref(isos1, yelmo1%bnd%z_bed_ref, yelmo1%bnd%H_ice_ref)
    call isos_init_state(isos1, yelmo1%bnd%z_bed, yelmo1%tpo%now%H_ice, ts%time, bsl)
    call isos_write_init_extended(isos1, file_isos, ts%time)

    yelmo1%bnd%z_bed = isos1%out%z_bed
    yelmo1%bnd%z_sl  = isos1%out%z_ss

    if (use_hyster) then
        ! Update hysteresis variable 
        var   = yelmo1%reg%V_ice*convert_km3_Gt
        dv_dt = sqrt(sum(yelmo1%tpo%now%dHidt**2)/real(count(yelmo1%tpo%now%f_ice .gt. 0.0),wp))
        call hyster_calc_forcing(hyst1,ts%time,var)
        !call hyster_calc_forcing(hyst1,ts%time,var,dv_dt)
        dT_summer = hyst1%f_now*hyst_f_ta
        dT_ann    = 1.3*dT_summer       ! == 0.5* ( (1.6)*dT_summer + (1.0)*dT_summer), given T_wintfac=1.6
        dT_ocn    = dT_ann*hyst_f_to 
    
    else
        dT_summer = 0.0
        dT_ann    = 0.0
        dT_ocn    = 0.0 
    end if 

    if (yelmo1%par%use_restart) then
        ! Call rembo_update once to get all fields over the year
        ! (restart file only contains information for last day of the year)
        call rembo_update(real(ts%time,8),real(ts%time_rel,8),real(dT_summer,8),real(yelmo1%tpo%now%z_srf,8), &
                                             real(yelmo1%tpo%now%H_ice,8),real(yelmo1%bnd%z_sl,8))
    else
        ! Equilibrate REMBO with correct topography, no temp anomaly, let it equilibrate for several years
        ! This will also calculate the precip correction factor

        call rembo_equilibrate(real(ts%time,8),real(yelmo1%tpo%now%z_srf,8), &
                                real(yelmo1%tpo%now%H_ice,8),real(yelmo1%bnd%z_sl,8),time_tot=10d0)
        
    end if
    
    ! Update surface mass balance and surface temperature from REMBO
    yelmo1%bnd%smb   = rembo_ann%smb    *yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = rembo_ann%T_srf
    
    ! Special treatment for Greenland
    if (lim_pd_ice) then 
        
        ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
        where(mask_noice) yelmo1%bnd%smb = yelmo1%bnd%smb - 4.0 

    end if 
    
    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=ts%time,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

    if (use_hyster .and. trim(snp1%par%ocn_type) .eq. "const") then 
        ! Apply oceanic anomaly from hyster method 
        
        snp1%now%to_ann = snp1%now%to_ann + dT_ocn

    end if 

    call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=ts%time,thrm_method="robin-cold")

    ! ===== basal friction optimization ======
    if (optimize) then 
        
        ! Ensure that cb_ref will be optimized (till_method == set externally) 
        yelmo1%dyn%par%till_method = -1  

        ! If not using restart...
        if (.not. yelmo1%par%use_restart) then

            if (opt%cf_init .gt. 0.0) then 
                ! Prescribe cb_ref to initial guess 
                yelmo1%dyn%now%cb_ref = opt%cf_init 
            else 
                ! Load cb_ref from calculated cb_tgt field
                yelmo1%dyn%now%cb_ref = yelmo1%dyn%now%cb_tgt 
            end if 

        end if 

    end if 
    
    if (greenland_init_marine_H) then
        ! Add extra ice-thickness over continental shelf to start with
        ! an LGM-like state
        ! (Note: this can be done even if running from a restart file...)
        
        ! Increase ice thickness everywhere to start
        yelmo1%tpo%now%H_ice = yelmo1%tpo%now%H_ice*1.2

        ! where(yelmo1%bnd%mask_ice /= -1 .and. yelmo1%tpo%now%H_ice .lt. 600.0 &
        !         .and. yelmo1%bnd%z_bed .gt. -500.0)

        !         yelmo1%tpo%now%H_ice = 800.0 

        ! end where

    end if
        
    if (yelmo1%par%use_restart) then 
        ! If using restart file, set boundary module variables equal to restarted value 

        ! Set boundary module variables equal to restarted value
        isos1%now%z_bed(isos1%domain%icrop1:isos1%domain%icrop2, &
            isos1%domain%jcrop1:isos1%domain%jcrop2)  = yelmo1%bnd%z_bed

    else 
        ! Run yelmo for several years with constant boundary conditions and topo
        ! to equilibrate thermodynamics and dynamics

        if (with_ice_sheet) then 
            call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_wp, dt=1.0_wp,topo_fixed=.FALSE.)
            call yelmo_update_equil(yelmo1,ts%time,time_tot=time_equil,dt=dtt,topo_fixed=.TRUE.)
        end if 

    end if 

    ! Heavy 2D file  
    call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years")

    ! 2D small file 
    ! call yelmo_write_init(yelmo1,file2D_small,time_init=ts%time,units="years")
    
    ! 1D file 
    !call yelmo_regions_write(yelmo1,ts%time,init=.TRUE.,units="years")

    ! Small 1D-2D yelmo-rembo file
    call yelmox_write_init(yelmo1,file_rembo,time_init=ts%time,units="years", &
                    mask=(yelmo1%bnd%mask_ice /= -1),dT_min=hyst1%par%f_min,dT_max=hyst1%par%f_max)

    call timer_step(tmr,comp=1,label="initialization") 
    call timer_step(tmrs,comp=-1)
    
    ! Write model state out to initial set of restart files
    call yelmox_restart_write(isos1,yelmo1,mshlf1,rembo_ann,ts%time)

    ! == Advance timesteps ===

    dtt_now = dtt

    call tstep_print_header(ts)

    do while (.not. ts%is_finished)

        ! Modify dtt and rembo timestepping for transient experiments 
        if (use_hyster) then

            select case(trim(hyst1%par%method))
                case("ramp-time","ramp-time-step")

                    deltat_tot = hyst1%par%dt_init + hyst1%par%dt_ramp + hyst1%par%dt_conv + 100.0

                    if (ts%time_elapsed .lt. deltat_tot) then
                        dtt_now = min(5.0,dtt)
                        rembo_ann%par%dtime_emb = dtt_now
                    else
                        dtt_now = dtt 
                        rembo_ann%par%dtime_emb = 100.0 
                    end if 

                case DEFAULT
                    ! Pass - normally do not change timestepping
            end select
        end if 
        
        ! == Update timestep ===

        call tstep_update(ts,dtt_now)
        call tstep_print(ts)
        
        ! == HYSTER boundary forcing ====================================

        if (use_hyster) then 
        
            ! Update forcing based on hysteresis module
            var   = yelmo1%reg%V_ice*convert_km3_Gt
            dv_dt = sqrt(sum(yelmo1%tpo%now%dHidt**2)/real(count(yelmo1%tpo%now%f_ice .gt. 0.0),wp))
            call hyster_calc_forcing(hyst1,ts%time,var)
            !call hyster_calc_forcing(hyst1,ts%time,var,dv_dt)
            write(*,*) "hyst: ", ts%time, hyst1%dt, hyst1%dv_dt_ave, hyst1%df_dt*1e6, hyst1%f_now 
            
            ! Store in dT_summer for forcing of rembo, etc.
            ! Also get regional dT_ann, given T_wintfac used by rembo
            dT_summer = hyst1%f_now*hyst_f_ta
            dT_ann    = 1.3*dT_summer       ! == 0.5* ( (1.6)*dT_summer + (1.0)*dT_summer), given T_wintfac=1.6
            dT_ocn    = dT_ann*hyst_f_to 

        else
            dT_summer = 0.0
            dT_ann    = 0.0
            dT_ocn    = 0.0 
        end if 

        call timer_step(tmrs,comp=0) 

        ! == ISOSTASY and SEA LEVEL ===========================================
        call bsl_update(bsl, ts%time_rel)
        call isos_update(isos1, yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
        yelmo1%bnd%z_bed = isos1%out%z_bed
        yelmo1%bnd%z_sl  = isos1%out%z_ss

        call timer_step(tmrs,comp=1,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="isostasy") 
        
        ! == Yelmo ice sheet ===================================================
        if (with_ice_sheet .and. (.not. (ts%n .eq. 1 .and. yelmo1%par%use_restart)) ) then

            if (optimize) then 

                ! === Optimization update step =========

                if (opt%opt_cf .and. &
                        (ts%time_elapsed .ge. opt%cf_time_init .and. ts%time_elapsed .le. opt%cf_time_end) ) then  
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                                    yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                                    yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                                    opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                                    dt=dtt_now,fill_method=opt%fill_method,fill_dist=opt%sigma_err,cb_tgt=yelmo1%dyn%now%cb_tgt)
                    
                end if

                if (opt%opt_tf .and. &
                        (ts%time_elapsed .ge. opt%tf_time_init .and. ts%time_elapsed .le. opt%tf_time_end) ) then
                    ! Perform tf_corr optimization

                    call optimize_tf_corr(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                                                        yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,opt%H_grnd_lim,yelmo1%bnd%basins, &
                                                        opt%basin_fill,opt%tau_m,opt%m_temp,opt%tf_min,opt%tf_max,yelmo1%tpo%par%dx,sigma=opt%tf_sigma,dt=dtt_now)
                
                end if 

            
            end if 
            
            ! Update ice sheet to current time 
            call yelmo_update(yelmo1,ts%time)
            
        end if 
        
        call timer_step(tmrs,comp=2,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="yelmo") 
        
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================

        if ( .not. (ts%n .eq. 1 .and. yelmo1%par%use_restart) ) then

            ! call REMBO1
            call rembo_update(real(ts%time,8),real(ts%time_rel,8),real(dT_summer,8),real(yelmo1%tpo%now%z_srf,8), &
                                            real(yelmo1%tpo%now%H_ice,8),real(yelmo1%bnd%z_sl,8))
            
            ! Update surface mass balance and surface temperature from REMBO
            yelmo1%bnd%smb   = rembo_ann%smb    *yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf = rembo_ann%T_srf
            
            ! Special treatment for Greenland
            if (lim_pd_ice) then 
            
                ! Impose additional negative mass balance to no ice points of 4 [m.i.e./a] melting
                where(mask_noice) yelmo1%bnd%smb = yelmo1%bnd%smb - 4.0 

            end if 

            ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
            
            ! Update snapclim
            call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=ts%time,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins) 

            if (use_hyster .and. trim(snp1%par%ocn_type) .eq. "const") then 
                ! Apply oceanic anomaly from hyster method 

                snp1%now%to_ann = snp1%now%to_ann + dT_ocn
                
            end if

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

        if (timeout_check(tm_1D,ts%time)) then  
            ! call yelmo_write_reg_step(yelmo1,file1D,time=ts%time) 
            !call yelmox_write_step_1D(yelmo1,hyst1,file1D_hyst,time=ts%time)
            call yelmox_write_step_small(yelmo1,hyst1,rembo_ann,isos1,mshlf1,file_rembo,ts%time, &
                                            dT_summer,dT_ann,dT_ocn,file_rembo_write_ocn_forcing)
        end if 

        if (timeout_check(tm_2D,ts%time)) then
            call yelmox_write_step(yelmo1,rembo_ann,isos1,mshlf1,file2D,time=ts%time)
        end if 

        if (write_restart .and. mod(ts%time,dt_restart)==0) then 
            call yelmox_restart_write(isos1,yelmo1,mshlf1,rembo_ann,ts%time)
        end if 

        call timer_step(tmrs,comp=4,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="io") 
        
        if (mod(ts%time_elapsed,10.0)==0) then
            ! Print timestep timing info and write log table
            call timer_write_table(tmrs,[ts%time,dtt_now]*1e-3,"m",tmr_file,init=ts%time_elapsed .eq. 0.0)
        end if 

        if (use_hyster .and. hyst1%kill) then 
            write(*,"(a,f12.3,a,f12.3)") "hyster:: kill switch activated. [time, f_now] = ", &
                                                        ts%time, ", ", hyst1%f_now 
            write(*,*) "hyster:: exiting time loop..."
            exit  
        end if 

    end do 

    ! Stop timing
    call timer_step(tmr,comp=2,time_mod=[ts%time_init,ts%time]*1e-3,label="timeloop") 
    
    ! Write the restart files for the end of the simulation
    if (write_restart) then 
        call yelmox_restart_write(isos1,yelmo1,mshlf1,rembo_ann,ts%time)
    end if

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=ts%time*1e-3)
    
contains
    
    subroutine yelmox_write_step_small(ylmo,hyst,rembo,isos,mshlf,filename,time, &
                                                            dT_jja,dT_ann,dT_ocn,write_ocn_forcing)

        implicit none 
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(hyster_class),   intent(IN) :: hyst
        type(rembo_class),    intent(IN) :: rembo
        type(isos_class),     intent(IN) :: isos 
        type(marshelf_class), intent(IN) :: mshlf 
        character(len=*),     intent(IN) :: filename
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: dT_jja
        real(wp), intent(IN) :: dT_ann
        real(wp), intent(IN) :: dT_ocn
        logical, intent(IN), optional :: write_ocn_forcing

        ! Local variables
        integer  :: ncid, n, k 
        real(wp) :: npmb, ntot, aar, smb_tot
        real(wp) :: dHidt_rms, dHidt_rms_1, dHidt_max   
        real(wp) :: dT_axis(1000)
        type(yregions_class) :: reg

        ! Assume region to write is the global region of yelmo 
        reg = ylmo%reg 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! ===== Hyst / forcing variables ===== 

        call nc_write(filename,"hyst_f_now",hyst%f_now,units="K",long_name="hyst: forcing value", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_df_dt",hyst%df_dt*1e6,units="K/(1e6 yr)",long_name="hyst: forcing rate of change", &
                      dim1="time",start=[n],ncid=ncid)
        ! call nc_write(filename,"hyst_dv_dt",hyst%dv_dt_ave,units="Gt/yr",long_name="hyst: volume rate of change", &
        !               dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_dv_dt",hyst%dv_dt_ave,units="m/yr",long_name="hyst: rms thickness rate of change", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! Write volume in volume-dT phase space
        call nc_read(filename,"dT_axis",dT_axis) 
        k = minloc(abs(dT_axis-hyst%f_now),dim=1)
        call nc_write(filename,"V_dT",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[k],ncid=ncid)
        
        ! == yelmo metrics ==

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        ! call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! == 1D Variables ==
        
        call nc_write(filename,"V_ice",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice",reg%A_ice*1e-6,units="1e6 km^2",long_name="Ice area", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_sle",reg%V_sle,units="m sle",long_name="Sea-level equivalent volume", &
                      dim1="time",start=[n],ncid=ncid)
        
        if (count(ylmo%tpo%now%f_ice .gt. 0.0) .gt. 0) then
            dHidt_rms = sqrt(sum(ylmo%tpo%now%dHidt**2)/real(count(ylmo%tpo%now%f_ice .gt. 0.0),wp))
            dHidt_max = maxval(abs(ylmo%tpo%now%dHidt),mask=ylmo%tpo%now%f_ice .gt. 0.0)
        else
            dHidt_rms = 0.0 
            dHidt_max = 0.0 
        end if

        if (count(ylmo%tpo%now%f_ice .gt. 0.0 .and. abs(ylmo%tpo%now%dHidt) .gt. 1e-3) .gt. 0) then
            dHidt_rms_1 = sqrt(sum(ylmo%tpo%now%dHidt**2) / &
                real(count(ylmo%tpo%now%f_ice .gt. 0.0 .and. abs(ylmo%tpo%now%dHidt) .gt. 1e-3),wp))
        else
            dHidt_rms_1 = 0.0
        end if

        ! call nc_write(filename,"rms(dHidt)",dHidt_rms,units="m/yr",long_name="rms ice thickness change", &
        !               dim1="time",start=[n],ncid=ncid)
        ! call nc_write(filename,"rms1(dHidt)",dHidt_rms_1,units="m/yr",long_name="rms ice thickness change", &
        !               dim1="time",start=[n],ncid=ncid)
        ! call nc_write(filename,"max(dHidt)",dHidt_max,units="m/yr",long_name="max. ice thickness change", &
        !               dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"dVidt",ylmo%reg%dVidt,units="km^3/a",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)
        
        if (n .eq. 1) then
            call nc_write(filename,"mask_ice",ylmo%bnd%mask_ice,units="",long_name="Ice mask", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
            call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if 

        ! == yelmo_topography ==
        call yelmo_write_var(filename,"H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"mask_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_net",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb",ylmo,n,ncid)
        
        ! == yelmo_dynamics ==
        call yelmo_write_var(filename,"uxy_s",ylmo,n,ncid)
        
        ! == yelmo_thermodymamics
        call yelmo_write_var(filename,"T_prime_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_w",ylmo,n,ncid)
        
        ! == yelmo_bound ==
        call yelmo_write_var(filename,"z_bed",ylmo,n,ncid)
        
        ! == rembo climate == 
        call nc_write(filename,"ta_ann",rembo%T_ann,units="K",long_name="REMBO Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ta_sum",rembo%T_jja,units="K",long_name="REMBO Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",rembo%pr*1e-3,units="m/a water equiv.",long_name="REMBO Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb_ann",rembo%smb*1e-3,units="m/yr water equiv.",long_name="REMBO Surface mass balance (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == ocean forcing ==
        if (write_ocn_forcing) then 

            call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if 

        ! == ice-sheet wide metrics == 
        ! Note: aar possibily interesting for studies on future melting
        
        ! Get integrated metrics (smb_tot [Gt/yr] and aar [unitless])
        ntot = count(ylmo%tpo%now%H_ice .gt. 0.0)

        if (ntot .gt. 0.0) then 
            npmb = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. rembo%smb .gt. 0.0)
            aar  = real(npmb,prec) / real(ntot,prec) 

            smb_tot = (ylmo%tpo%par%dx**2)*sum(rembo%smb*1e-3,mask=ylmo%tpo%now%H_ice .gt. 0.0)

            ! Convert from m^3/yr => Gt/yr
            ! [m^3/yr] * [1000 kg/m^3] * [1e-12 Gt/kg] == [Gt/yr]
            smb_tot = smb_tot * (1000) *1e-12 
        else
            aar = 0.0
        end if  

        call nc_write(filename,"dT_jja",hyst%f_now,units="K",long_name="Temp. anomaly, regional JJA mean", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dT_ann",dT_ann,units="K",long_name="Temp. anomaly, regional annual mean", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dT_ocn",dT_ocn,units="K",long_name="Temp. anomaly, regional oceanic mean", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"smb_mean",smb_tot,units="Gt/yr",long_name="Mean smb over the ice sheet", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"aar",aar,units="1",long_name="Accumulation area ratio", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step_small

    subroutine yelmox_write_step(ylmo,rembo,isos,mshlf,filename,time)

        implicit none 
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(rembo_class),    intent(IN) :: rembo
        type(isos_class),     intent(IN) :: isos 
        type(marshelf_class), intent(IN) :: mshlf 
        character(len=*),     intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n  
        character(len=12) :: dims3(3)

        ! Define useful dimensions for ncio writing
        dims3(1) = "xc"
        dims3(2) = "yc"
        dims3(3) = "time"

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
        
        if (n .eq. 1) then
            call nc_write(filename,"mask_ice",ylmo%bnd%mask_ice,units="",long_name="Ice mask", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
            call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if
        
        ! == yelmo_topography ==
        call yelmo_write_var(filename,"H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"mask_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_net",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"N_eff",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"dHidt",ylmo,n,ncid)
        
        ! == yelmo_dynamics ==
        call yelmo_write_var(filename,"cb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"c_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_eff_int",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_s",ylmo,n,ncid)
        
        ! == yelmo_material ==
        call yelmo_write_var(filename,"enh_bar",ylmo,n,ncid)
        !call yelmo_write_var(filename,"ATT",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_int",ylmo,n,ncid)
        
        ! == yelmo_thermodynamics ==
        call yelmo_write_var(filename,"T_prime",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_pmp",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_w",ylmo,n,ncid)
        
        ! == yelmo_boundaries ==
        call yelmo_write_var(filename,"z_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_sl",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"T_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_shlf",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_geo",ylmo,n,ncid)
        
        ! == yelmo_data (comparison with present-day) ==
        call yelmo_write_var(filename,"pd_err_H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"pd_err_z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"pd_err_uxy_s",ylmo,n,ncid)
        
        ! == FastIsostasy ==

        call nc_write(filename,"dzbdt",isos%out%dwdt,units="m/a", &
                    long_name="Bedrock uplift rate", dims=dims3,start=[1,1,n],ncid=ncid)
        ! call nc_write(filename, "w_viscous", isos%out%w, units="m", &
        !             long_name="Viscous vertical displacement", dims=dims3,start=[1,1,n], ncid=ncid)
        ! call nc_write(filename, "w_elastic", isos%out%we, units="m", &
        !             long_name="Elastic vertical displacement", dims=dims3,start=[1,1,n], ncid=ncid)
        ! call nc_write(filename, "canom", isos%out%canom_full, units="Pa", &
        !             long_name="Full stress anomaly", dims=dims3,start=[1,1,n], ncid=ncid)
        ! call nc_write(filename, "maskocean", isos%out%maskocean, units="", &
        !             long_name="Ocean mask", dims=dims3,start=[1,1,n], ncid=ncid)
        ! call nc_write(filename, "maskgrounded", isos%out%maskgrounded, units="", &
        !             long_name="Grounded mask", dims=dims3,start=[1,1,n], ncid=ncid)
        ! call nc_write(filename, "maskcontinent", isos%out%maskcontinent, units="", &
        !             long_name="Continent mask", dims=dims3,start=[1,1,n], ncid=ncid)
        
        ! == marine_shelf ==
        
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K", &
                    long_name="Shelf temperature anomaly",dims=dims3,start=[1,1,n],ncid=ncid)

        ! == rembo_annual ==

        call nc_write(filename,"Ta_ann",rembo%T_ann,units="K",long_name="REMBO Near-surface air temperature (ann)", &
                      dims=dims3,start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",rembo%T_jja,units="K",long_name="REMBO Near-surface air temperature (sum)", &
                      dims=dims3,start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",rembo%pr*1e-3,units="m/a water equiv.",long_name="REMBO Precipitation (ann)", &
                      dims=dims3,start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb_ann",rembo%smb*1e-3,units="m/a water equiv.",long_name="REMBO Surface mass balance (ann)", &
                      dims=dims3,start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
        !               dims=dims3,start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
        !               dims=dims3,start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
        !               dims=dims3,start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step

    subroutine yelmox_write_step_1D(ylm,hyst,filename,time)

        implicit none 
        
        type(yelmo_class),  intent(IN) :: ylm
        type(hyster_class), intent(IN) :: hyst 
        character(len=*),   intent(IN) :: filename
        real(wp),           intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, k  
        real(wp) :: dT_axis(1000) 
        type(yregions_class) :: reg

        ! Assume region to write is the global region of yelmo 
        reg = ylm%reg 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)


        ! ===== Hyst / forcing variables ===== 

        call nc_write(filename,"hyst_f_now",hyst%f_now,units="K",long_name="hyst: forcing value", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_df_dt",hyst%df_dt*1e6,units="K/(1e6 yr)",long_name="hyst: forcing rate of change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_dv_dt",hyst%dv_dt,units="Gt/yr",long_name="hyst: volume rate of change", &
                      dim1="time",start=[n],ncid=ncid)

        ! Write volume in volume-dT phase space
        call nc_read(filename,"dT_axis",dT_axis) 
        k = minloc(abs(dT_axis-hyst%f_now),dim=1)
        call nc_write(filename,"V_dT",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[k],ncid=ncid)

        ! ===== Total ice variables =====
        
        call nc_write(filename,"H_ice",reg%H_ice,units="m",long_name="Mean ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"z_srf",reg%z_srf,units="m",long_name="Mean surface elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dHidt",reg%dHidt,units="m/a",long_name="Mean rate ice thickness change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_ice_max",reg%H_ice_max,units="m/a",long_name="Max ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dzsdt",reg%dzsdt,units="m/a",long_name="Mean rate surface elevation change", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"V_ice",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice",reg%A_ice*1e-6,units="1e6 km^2",long_name="Ice area", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dVidt",reg%dVidt,units="km^3/a",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fwf",reg%fwf,units="Sv",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar",reg%uxy_bar,units="m/a",long_name="Mean depth-ave velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s",reg%uxy_s,units="m/a",long_name="Mean surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b",reg%uxy_b,units="m/a",long_name="Mean basal velocity", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"z_bed",reg%z_bed,units="m",long_name="Mean bedrock elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"smb",reg%smb,units="m/a",long_name="Mean surface mass balance", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",reg%T_srf,units="K",long_name="Mean surface temperature", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb",reg%bmb,units="m/a",long_name="Mean total basal mass balance", &
                      dim1="time",start=[n],ncid=ncid)
        


        ! ===== Grounded ice variables =====

        call nc_write(filename,"H_ice_g",reg%H_ice_g,units="m",long_name="Mean ice thickness (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"z_srf_g",reg%z_srf_g,units="m",long_name="Mean surface elevation (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_ice_g",reg%V_ice_g*1e-6,units="1e6 km^3",long_name="Ice volume (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice_g",reg%A_ice_g*1e-6,units="1e6 km^2",long_name="Ice area (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar_g",reg%uxy_bar_g,units="m/a",long_name="Mean depth-ave velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s_g",reg%uxy_s_g,units="m/a",long_name="Mean surface velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b_g",reg%uxy_b_g,units="m/a",long_name="Mean basal velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"f_pmp",reg%f_pmp,units="1",long_name="Temperate fraction (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"H_w",reg%H_w,units="m",long_name="Mean basal water thickness (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"bmb_g",reg%bmb_g,units="m/a",long_name="Mean basal mass balance (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! ===== Floating ice variables =====

        call nc_write(filename,"H_ice_f",reg%H_ice_f,units="m",long_name="Mean ice thickness (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_ice_f",reg%V_ice_f*1e-6,units="1e6 km^3",long_name="Ice volume (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice_f",reg%A_ice_f*1e-6,units="1e6 km^2",long_name="Ice area (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar_f",reg%uxy_bar_f,units="m/a",long_name="Mean depth-ave velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s_f",reg%uxy_s_f,units="m/a",long_name="Mean surface velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b_f",reg%uxy_b_f,units="m/a",long_name="Mean basal velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"z_sl",reg%z_sl,units="m",long_name="Mean sea level (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",reg%bmb_shlf,units="m/a",long_name="Mean basal mass balance (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_shlf",reg%T_shlf,units="K",long_name="Mean marine shelf temperature (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step_1D

    subroutine yelmox_write_init(dom,filename,time_init,units,mask,dT_min,dT_max)

        implicit none 

        type(yelmo_class), intent(IN) :: dom 
        character(len=*),  intent(IN) :: filename, units 
        real(wp),          intent(IN) :: time_init
        logical,           intent(IN) :: mask(:,:) 
        real(wp),          intent(IN) :: dT_min 
        real(wp),          intent(IN) :: dT_max

        ! Local variables
        integer :: n
        real(wp), allocatable :: dT_axis(:) 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",        x=dom%grd%xc*1e-3,      units="kilometers")
        call nc_write_dim(filename,"yc",        x=dom%grd%yc*1e-3,      units="kilometers")
        call nc_write_dim(filename,"zeta",      x=dom%par%zeta_aa,      units="1")
        call nc_write_dim(filename,"time",      x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)
        
        !============================================================
        ! Add temperature axis to 1D hysteresis file 
        allocate(dT_axis(1000))
        do n = 1, 1000 
            dT_axis(n) = dT_min + (dT_max-dT_min)*(n-1)/real(1000-1,wp)
        end do 
        call nc_write_dim(filename,"dT_axis",x=dT_axis,units="degC")
        
        ! Populate variable with missing values for now to initialize it
        dT_axis = missing_value 
        call nc_write(filename,"V_dT",dT_axis,dim1="dT_axis",missing_value=missing_value)
        !============================================================
        
        ! Static information
        call nc_write(filename,"mask", mask,  units="1",long_name="Region mask",dim1="xc",dim2="yc")
        
        return

    end subroutine yelmox_write_init

    subroutine yelmox_restart_write(isos,ylmo,mshlf,rembo_ann,time,fldr)

        implicit none

        type(isos_class),   intent(IN) :: isos
        type(yelmo_class),  intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        type(rembo_class),  intent(IN) :: rembo_ann
        real(wp),           intent(IN) :: time 
        character(len=*),   intent(IN), optional :: fldr
        
        ! Local variables
        real(wp) :: time_kyr
        character(len=32)   :: time_str
        character(len=1024) :: outfldr

        character(len=56), parameter :: file_isos  = "isos_restart.nc"
        character(len=56), parameter :: file_yelmo = "yelmo_restart.nc"
        character(len=56), parameter :: file_mshlf = "marine_shelf.nc"
        character(len=56), parameter :: file_rembo = "rembo_restart.nc"
        
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
        
        call isos_restart_write(isos,trim(outfldr)//"/"//file_isos,time)
        call yelmo_restart_write(ylmo,trim(outfldr)//"/"//file_yelmo,time) 
        call marshelf_restart_write(mshlf,trim(outfldr)//"/"//file_mshlf,time)
        
        call rembo_restart_write(trim(outfldr)//"/"//file_rembo,real(time,dp),real(ylmo%tpo%now%z_srf,dp), &
                    real(ylmo%tpo%now%H_ice,dp),real(ylmo%bnd%z_sl,dp))
        
        return

    end subroutine yelmox_restart_write

end program yelmox



