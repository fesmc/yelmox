

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
    use smbpal
    use sediments
    use geothermal
    
    implicit none 

    type(tstep_class)      :: ts
    
    type(yelmo_class)      :: yelmo1
    type(bsl_class)        :: bsl
    type(snapclim_class)   :: snp1
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

    ! === Initialize external models (forcing for ice sheet) ======

    ! Initialize barysealevel model
    call bsl_init(bsl, path_par, ts%time_rel)

    ! Initialize fastisosaty
    call isos_init(isos1, path_par, "isos", yelmo1%grd%nx, yelmo1%grd%ny, yelmo1%grd%dx, yelmo1%grd%dy)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)
    
    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)

    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    yelmo1%bnd%H_sed = sed1%now%H 
    
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

    ! Barystatic sea level
    call bsl_update(bsl, year_bp=ts%time_rel)
    call bsl_write_init(bsl, file_bsl, ts%time)

    ! Initialize the isostasy reference state using reference topography fields
    call isos_init_ref(isos1, yelmo1%bnd%z_bed_ref, yelmo1%bnd%H_ice_ref)
    call isos_init_state(isos1, yelmo1%bnd%z_bed, yelmo1%tpo%now%H_ice, ts%time, bsl)
    call isos_write_init_extended(isos1, file_isos, ts%time)

    yelmo1%bnd%z_bed = isos1%out%z_bed
    yelmo1%bnd%z_sl  = isos1%out%z_ss

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

    ! Set oceanic boundary conditions to zero for mountain glaciers.
    yelmo1%bnd%bmb_shlf = 0.0_wp
    yelmo1%bnd%T_shlf   = 0.0_wp  

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
        if (ctl%with_ice_sheet) then
            ! Run yelmo for a few years with constant boundary conditions
            ! to synchronize all model fields a bit
            call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)
        end if 
    end if 

    ! ===== Initialize output files ===== 
    
    call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years") 
    call yelmo_write_init(yelmo1,file2D_small,time_init=ts%time,units="years") 
    
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
        
        ! Spin-up procedure - only relevant for time-time_init <= time_equil
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
        call bsl_update(bsl, ts%time_rel)
        call isos_update(isos1, yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
        yelmo1%bnd%z_bed = isos1%out%z_bed
        yelmo1%bnd%z_sl  = isos1%out%z_ss

        call timer_step(tmrs,comp=1,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="isostasy") 
        
        ! == ICE SHEET ===================================================

        ! Update Yelmo
        if (ctl%with_ice_sheet) then
            call yelmo_update(yelmo1,ts%time)
        end if
        
        call timer_step(tmrs,comp=2,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="yelmo")
        
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        
        if (mod(nint(ts%time*100),nint(ctl%dt_clim*100))==0) then
                ! Update snapclim
                call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=ts%time,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins) 
        end if 

        ! == SURFACE MASS BALANCE ==============================================
        call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                   yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,ts%time_rel) 
        yelmo1%bnd%smb   = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        yelmo1%bnd%bmb_shlf = 0.0_wp
        yelmo1%bnd%T_shlf   = 0.0_wp
        
        call timer_step(tmrs,comp=3,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="climate") 

        ! == MODEL OUTPUT =======================================================

        if (timeout_check(tm_2D,ts%time)) then
            call yelmox_write_step(yelmo1,snp1,smbpal1,file2D,ts%time)
        end if

        if (timeout_check(tm_2Dsm,ts%time) .and. .FALSE.) then 
            call yelmo_write_step(yelmo1,file2D_small,ts%time,compare_pd=.FALSE.)
        end if

        if (timeout_check(tm_1D,ts%time)) then 
            call yelmo_regions_write(yelmo1,ts%time)
        end if 

        if (mod(nint(ts%time*100),nint(ctl%dt_restart*100))==0) then
            call yelmox_restart_write(bsl,isos1,yelmo1,ts%time)
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
    call yelmox_restart_write(bsl,isos1,yelmo1,ts%time)

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=ts%time*1e-3)
    
contains
    
    subroutine yelmox_write_step(ylmo,snp,srf,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(snapclim_class),   intent(IN) :: snp 
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
        !call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
                
        ! == yelmo_topography ==
        call yelmo_write_var(filename,"H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"mask_bed",ylmo,n,ncid)
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
        
        ! == yelmo_boundaries ==
        call yelmo_write_var(filename,"z_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_sl",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb",ylmo,n,ncid)
        call yelmo_write_var(filename,"T_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_shlf",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_geo",ylmo,n,ncid)

        if (.FALSE.) then
                ! == yelmo_data (comparison with present-day) ==
                call yelmo_write_var(filename,"pd_err_H_ice",ylmo,n,ncid)
                call yelmo_write_var(filename,"pd_err_z_srf",ylmo,n,ncid)
                call yelmo_write_var(filename,"pd_err_uxy_s",ylmo,n,ncid)
                call yelmo_write_var(filename,"pd_err_smb_ref",ylmo,n,ncid)
        end if

        ! == yelmo extra fields ==
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

        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Static fields
        if (n .le. 1) then 
            call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if

        ! == snapclim ==
        call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_bnd",ylmo%bnd%smb,units="m/a water equiv.",long_name="Boundary SMB", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
              
        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == smbpal ==
!         call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step

    subroutine yelmox_restart_write(bsl,isos,ylmo,time,fldr)

        implicit none

        type(bsl_class),      intent(IN) :: bsl
        type(isos_class),     intent(IN) :: isos
        type(yelmo_class),    intent(IN) :: ylmo
        real(wp),             intent(IN) :: time 
        character(len=*),     intent(IN), optional :: fldr
        
        ! Local variables
        real(wp) :: time_kyr
        character(len=32)   :: time_str
        character(len=1024) :: outfldr

        character(len=56), parameter :: file_bsl   = "bsl_restart.nc"
        character(len=56), parameter :: file_isos  = "isos_restart.nc"
        character(len=56), parameter :: file_yelmo = "yelmo_restart.nc"

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

        return

    end subroutine yelmox_restart_write

end program yelmox
