
program yelmox_esm

    use nml
    use ncio
    use timestepping
    use timer
    use timeout
    use yelmo
    use ice_optimization
    use ice_sub_regions
    
    ! External libraries
    use geothermal
    use esm
    use fastisostasy    ! also reexports barysealevel
    use marine_shelf
    use sediments
    use smbpal
    
    implicit none 

    character(len=256) :: outfldr, file2D, file2Dsm
    character(len=256) :: file1D_esm, file2D_esm
    character(len=256) :: file_isos, file_bsl
    character(len=256) :: domain, grid_name 
    character(len=512) :: path_par  
    character(len=512) :: path_tf_corr 
    character(len=512) :: esm_path_par
    integer  :: n, m
    real(wp) :: time_wt
        
    real(sp) :: convert_km3_Gt
        
    type(tstep_class)           :: ts
            
    type(yelmo_class)           :: yelmo1
    type(bsl_class)             :: bsl
    type(marshelf_class)        :: mshlf1
    type(smbpal_class)          :: smbpal1
    type(sediments_class)       :: sed1
    type(geothermal_class)      :: gthrm1
    type(isos_class)            :: isos1
    type(esm_forcing_class)     :: esm1

    logical,  allocatable :: tmp_mask(:,:)
    type(timeout_class) :: tm_1D, tm_2D, tm_2Dsm

    ! Model timing
    type(timer_class)  :: tmr
    type(timer_class)  :: tmrs
    character(len=512) :: tmr_file 
    
    type ctrl_params
        character(len=56)   :: run_step
        real(wp)            :: time_init
        real(wp)            :: time_end
        real(wp)            :: time_equil   ! Only for spinup 
        real(wp)            :: dtt

        real(wp)            :: time_ref(2)
        real(wp)            :: time_hist(2) 
        real(wp)            :: time_proj(2)  
        real(wp)            :: time_esm_ref(2)  
        character(len=56)   :: clim_var
        integer             :: clim_seed

        character(len=56)   :: tstep_method
        real(wp)            :: tstep_const

        logical             :: kill_shelves
        logical             :: with_ice_sheet 
        character(len=56)   :: equil_method

        character(len=512)  :: esm_par_file
        character(len=56)   :: esm_experiment
        logical             :: esm_use_esm
        character(len=56)   :: esm_name
        logical             :: esm_write_formatted
        real(wp)            :: esm_dt_formatted

        real(wp)            :: isos_tau_1 
        real(wp)            :: isos_tau_2 
        real(wp)            :: isos_sigma 

    end type 

    type(ctrl_params)     :: ctl
    type(ice_opt_params)  :: opt 
    type(esm_experiment_class) :: esmexp 

    ! Set the seed for reproducibility of randomnes (climate variability)
    integer :: seed_size
    integer, allocatable :: seed(:)
    call random_seed(size=seed_size) ! Get the size of the seed array
    allocate(seed(seed_size))        ! Allocate the seed array

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Control parameters 
    call nml_read(path_par,"ctrl","run_step",        ctl%run_step)

    ! esm parameters 
    call nml_read(path_par,"esm","par_file",         ctl%esm_par_file)
    call nml_read(path_par,"esm","experiment",       ctl%esm_experiment)
    call nml_read(path_par,"esm","esm_name",         ctl%esm_name)
    call nml_read(path_par,"esm","use_esm",          ctl%esm_use_esm)

    ! esm cmip conventions
    call nml_read(path_par,"esm","write_formatted",  ctl%esm_write_formatted)
    call nml_read(path_par,"esm","dt_formatted",     ctl%esm_dt_formatted)

    ! esm climate parameters
    call nml_read(path_par,"esm","lapse",            esm1%lapse)
    call nml_read(path_par,"esm","f_p",              esm1%beta_p)
    call nml_read(path_par,"esm","f_ocn",            esm1%f_ocn)
    call nml_read(path_par,"esm","f_polar",          esm1%f_polar)
    call nml_read(path_par,"esm","dT_threshold",     esm1%dT_lim)

    ! esm source grid
    call nml_read(path_par,"esm","grid_src",         esm1%grid_src)

    ! What does this?: needed?
    if (index(ctl%esm_par_file,"ant") .gt. 0) then
        ! Running Antarctica domain, load Antarctica specific parameters
        call esm_experiment_def(esmexp,ctl%esm_experiment,ctl%esm_par_file,"UCM","YELMO",yelmo1%par%domain)
    end if

    ! Read run_step specific control parameters
    call nml_read(path_par,trim(ctl%run_step),"time_init",    ctl%time_init)    ! [yr] Starting time
    call nml_read(path_par,trim(ctl%run_step),"time_end",     ctl%time_end)     ! [yr] Ending time
    call nml_read(path_par,trim(ctl%run_step),"dtt",          ctl%dtt)          ! [yr] Main loop time step 
    call nml_read(path_par,trim(ctl%run_step),"time_equil",   ctl%time_equil)   ! [yr] Years to relax first
    call nml_read(path_par,trim(ctl%run_step),"time_ref",     ctl%time_ref)     ! [yr,yr] Reference period
    call nml_read(path_par,trim(ctl%run_step),"time_hist",    ctl%time_hist)    ! [yr,yr] Historical period
    call nml_read(path_par,trim(ctl%run_step),"time_proj",    ctl%time_proj)    ! [yr,yr] Projection period
    call nml_read(path_par,trim(ctl%run_step),"time_esm_ref", ctl%time_esm_ref) ! [yr,yr] ESM reference period
    call nml_read(path_par,trim(ctl%run_step),"clim_var",     ctl%clim_var)     ! Climate variability
    call nml_read(path_par,trim(ctl%run_step),"clim_seed",    ctl%clim_seed)     ! Climate variability
    call nml_read(path_par,trim(ctl%run_step),"tstep_method", ctl%tstep_method) ! Calendar choice ("const" or "rel")
    call nml_read(path_par,trim(ctl%run_step),"tstep_const",  ctl%tstep_const)  ! Assumed time bp for const method
        
    call nml_read(path_par,trim(ctl%run_step),"kill_shelves",  ctl%kill_shelves)    ! Kill shelves beyond pd?
    call nml_read(path_par,trim(ctl%run_step),"with_ice_sheet",ctl%with_ice_sheet)  ! Active ice sheet? 
    call nml_read(path_par,trim(ctl%run_step),"equil_method",  ctl%equil_method)    ! What method should be used for spin-up?

    if (trim(ctl%equil_method) .eq. "opt") then 
        ! Load optimization parameters 

        call optimize_par_load(opt,path_par,"opt")

    end if 

    ! Get output times
    call timeout_init(tm_1D,  path_par,"tm_1D",  "small",  ctl%time_init,ctl%time_end)
    call timeout_init(tm_2D,  path_par,"tm_2D",  "heavy",  ctl%time_init,ctl%time_end)
    call timeout_init(tm_2Dsm,path_par,"tm_2Dsm","medium", ctl%time_init,ctl%time_end)
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations
    file2D           = trim(outfldr)//"yelmo2D.nc"
    file2Dsm         = trim(outfldr)//"yelmo2Dsm.nc"
    
    file_isos        = trim(outfldr)//"fastisostasy.nc"
    file_bsl         = trim(outfldr)//"bsl.nc"
    
    file1D_esm       = trim(outfldr)//"yelmo1D_esm.nc"
    file2D_esm       = trim(outfldr)//"yelmo2D_esm.nc"

    tmr_file         = trim(outfldr)//"timer_table.txt"

    !  =========================================================
    ! Print summary of run settings 
    write(*,*)
    write(*,*) "run_step:  ", trim(ctl%run_step) 
    write(*,*)
    write(*,*) "time_init: ",   ctl%time_init 
    write(*,*) "time_end:  ",   ctl%time_end 
    write(*,*) "dtt:       ",   ctl%dtt  
    write(*,*) 
    
    write(*,*) "esm_par_file:        ", trim(ctl%esm_par_file)
    write(*,*) "esm_experiment:      ", trim(ctl%esm_experiment)
    write(*,*) "load_esm:            ", ctl%esm_use_esm
    if (ctl%esm_use_esm) write(*,*) "esm_name: ", trim(ctl%esm_name)

    ! set the seed
    seed = ctl%clim_seed                 ! Set the seed value (any integer value as the seed)
    call random_seed(put=seed)       ! Initialize the random number generator with the seed
    
    select case(trim(ctl%run_step))

        case("spinup")

            write(*,*) "time_equil: ",   ctl%time_equil 
            write(*,*) "tstep_const: ",  ctl%tstep_const 
            write(*,*) "time_opt: ",     ctl%time_ref(1),ctl%time_ref(2)

        case("transient")

            write(*,*) "esm_write_formatted: ", ctl%esm_write_formatted
            write(*,*) "esm_file_suffix:     ", trim(ctl%esm_name)
            write(*,*) "time_ref: ",            ctl%time_ref(1),ctl%time_ref(2)
            write(*,*) "time_hist: ",           ctl%time_hist(1),ctl%time_hist(2)
            write(*,*) "time_proj: ",           ctl%time_proj(1),ctl%time_proj(2)
            write(*,*) "time_esm_ref: ",        ctl%time_esm_ref(1),ctl%time_esm_ref(2)
            
    end select

    write(*,*) 
    write(*,*) "with_ice_sheet: ",  ctl%with_ice_sheet
    write(*,*) "equil_method:   ",  trim(ctl%equil_method)
    
    ! Start timing
    call timer_step(tmr,comp=-1) 

    ! === Initialize timestepping ===
    call tstep_init(ts,ctl%time_init,ctl%time_end,method=ctl%tstep_method,units="year", &
                    time_ref=2000.0_wp,const_rel=0.0_wp,const_cal=ctl%tstep_const)
    
    ! === Initialize ice sheet model =====
    
    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=ts%time)

    ! Store domain and grid_name as shortcuts 
    domain    = yelmo1%par%domain 
    grid_name = yelmo1%par%grid_name 

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
    
    end select


    ! === Initialize external models (forcing for ice sheet) ======

    ! Initialize barysealevel model
    call bsl_init(bsl, path_par, ts%time_rel)

    ! Initialize fastisosaty
    call isos_init(isos1, path_par, "isos", yelmo1%grd%nx, yelmo1%grd%ny, yelmo1%grd%dx, yelmo1%grd%dy)
    
    ! Initialize ESM atmospheric and oceanic objects 
    esm_path_par = trim(outfldr)//"/"//trim(ctl%esm_par_file)
    call esm_forcing_init(esm1,esm_path_par,domain,grid_name,gcm=ctl%esm_name,scenario=ctl%esm_experiment,use_esm=ctl%esm_use_esm)

    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)

    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! === Update external modules and pass variables to yelmo boundaries =======

    ! Sediments
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name)
    yelmo1%bnd%H_sed = sed1%now%H 

    ! Geothermal heat flow
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name)
    yelmo1%bnd%Q_geo = gthrm1%now%ghf 

    ! Barystatic sea level
    call bsl_update(bsl, year_bp=ts%time_rel)
    call bsl_write_init(bsl, file_bsl, ts%time)

    ! Initialize the isostasy reference state using reference topography fields
    call isos_init_ref(isos1, yelmo1%bnd%z_bed_ref, yelmo1%bnd%H_ice_ref)
    call isos_init_state(isos1, yelmo1%bnd%z_bed, yelmo1%tpo%now%H_ice, ts%time, bsl)
    call isos_write_init_extended(isos1, file_isos, ts%time)

    yelmo1%bnd%z_bed = isos1%out%z_bed
    yelmo1%bnd%z_sl  = isos1%out%z_ss

    ! Update forcing to present-day reference using esm forcing
    call calc_climate_esm(smbpal1,mshlf1,esm1,yelmo1,ctl, &
                          time=ts%time,time_bp=ts%time_rel, &
                          domain=yelmo1%par%domain,grid_name=yelmo1%par%grid_name)

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,esm1%t2m+esm1%dts,esm1%pr*esm1%dpr, &
            yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,ts%time_rel,time_equil=100.0)
    end if 

    ! Update Yelmo boundary fields
    yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=ts%time,thrm_method="robin-cold")

    ! Set new boundary conditions: if desired kill ice shelves beyond present-day extent
    if (ctl%kill_shelves) then
        where(yelmo1%dta%pd%mask_bed .eq. mask_bed_ocean) yelmo1%bnd%mask_ice = -1
    end if
    
! ================= RUN STEPS ===============================================


    select case(trim(ctl%run_step)) 

    case("spinup")
        ! Model can start from no spinup or equilibration (using restart file), 
        ! here it is run under constant boundary conditions to spinup 
        ! desired state. 

        write(*,*)
        write(*,*) "Performing spinup."
        write(*,*) 

        ! ===== basal friction optimization ======
        if (trim(ctl%equil_method) .eq. "opt") then 
            
            ! Ensure that cb_ref will be optimized (till_method == set externally) 
            yelmo1%dyn%par%till_method = -1  

            ! If not using restart...
            if (.not. yelmo1%par%use_restart) then

                if (opt%cf_init .gt. 0.0) then 
                    ! Prescribe cb_ref to initial guess 
                    yelmo1%dyn%now%cb_ref =  opt%cf_init
                else 
                    ! Load cb_ref from calculated cb_tgt field
                    yelmo1%dyn%now%cb_ref = yelmo1%dyn%now%cb_tgt 

                end if 

            end if 

        end if 

        if (ctl%with_ice_sheet .and. .not. yelmo1%par%use_restart) then 
            ! Run yelmo alone for one or a few years with constant boundary conditions
            ! to sort out inconsistencies from initialization.
            call yelmo_update_equil(yelmo1,ts%time,time_tot=1.0_wp,dt=1.0_wp,topo_fixed=.FALSE.)
        end if 

        if (trim(ctl%equil_method) .eq. "opt") then 
            ! Additional initialization option when running 'opt' spinup...

            if (ctl%with_ice_sheet .and. ctl%time_equil .gt. 0.0) then 
                ! Calculate thermodynamics with fixed ice sheet 
                call yelmo_update_equil(yelmo1,ts%time,time_tot=ctl%time_equil,dt=ctl%dtt,topo_fixed=.TRUE.)
            end if 

        end if 

        write(*,*) "Initialization complete."

        ! Initialize output files for checking progress 
        call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years") 
        call yelmo_write_init(yelmo1,file2Dsm,time_init=ts%time,units="years")  
        call yelmo_regions_write(yelmo1,ts%time,init=.TRUE.,units="years")

        call timer_step(tmr,comp=1,label="initialization") 
        call timer_step(tmrs,comp=-1)
        
        do while (.not. ts%is_finished)

            ! == Update timestep ===

            call tstep_update(ts,ctl%dtt)
            call tstep_print(ts)
            
            !!ajr: only update optimized fields if ice sheet is running
            if (ctl%with_ice_sheet) then
             
            select case(trim(ctl%equil_method))
            
                case("opt")

                    if (ts%time_elapsed .le. opt%rel_time2) then 
                        ! Apply relaxation to the model 

                        ! Update model relaxation time scale and error scaling (in [m])
                        call optimize_set_transient_param(opt%rel_tau,ts%time_elapsed,time1=opt%rel_time1, &
                                                time2=opt%rel_time2,p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)
                        
                        ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                        yelmo1%tpo%par%topo_rel_tau = opt%rel_tau 
                        yelmo1%tpo%par%topo_rel     = 3
                    
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
            
            end if 

            ! ====================================================

            call timer_step(tmrs,comp=0) 
            
            ! == ISOSTASY and SEA LEVEL ===========================================
            ! jablasco: avoid to update bedrock during spinup
            if (.False.) then
                call bsl_update(bsl, ts%time_rel)
                call isos_update(isos1, yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
                yelmo1%bnd%z_bed = isos1%out%z_bed
                yelmo1%bnd%z_sl  = isos1%out%z_ss
            else
                ! Do nothing
            end if
            call timer_step(tmrs,comp=1,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="isostasy") 

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,ts%time)

            call timer_step(tmrs,comp=2,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="yelmo") 

            ! == CLIMATE ===========================================================

            ! Update forcing to present-day reference, but 
            ! adjusting to ice topography
            call calc_climate_esm(smbpal1,mshlf1,esm1,yelmo1,ctl, &
                                  time=ts%time,time_bp=ts%time_rel, &
                                  domain=yelmo1%par%domain,grid_name=yelmo1%par%grid_name) 

            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

            call timer_step(tmrs,comp=3,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="climate") 

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_2Dsm,ts%time)) then 
                call write_step_2D_small(yelmo1,isos1,esm1,mshlf1,smbpal1,file2Dsm,ts%time)
            end if

            if (timeout_check(tm_2D,ts%time)) then
                call write_step_2D_combined(yelmo1,isos1,esm1,mshlf1,smbpal1,file2D,ts%time)
            end if

            if (timeout_check(tm_1D,ts%time)) then
                call yelmo_regions_write(yelmo1,ts%time)
            end if

            ! isostasy checks
            if (.FALSE. .and. (timeout_check(tm_2Dsm, ts%time))) then
                call isos_write_step_extended(isos1, file_isos, ts%time)
            end if 

            call timer_step(tmrs,comp=4,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="io") 
        
            if (mod(ts%time_elapsed,10.0)==0) then
                ! Print timestep timing info and write log table
                call timer_write_table(tmrs,[ts%time,ctl%dtt]*1e-3,"m",tmr_file,init=ts%time_elapsed .eq. 0.0)
            end if 

            if (mod(ts%time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", ts%time
            end if 

        end do

        write(*,*)
        write(*,*) "spinup complete."
        write(*,*)

        ! Write the restart snapshot for the end of the simulation
        call yelmox_restart_write(bsl,isos1,yelmo1,mshlf1,ts%time_rel)

    case("transient")
        ! Here it is assumed that the model has gone through spinup 
        ! and is ready for transient simulations 

        write(*,*)
        write(*,*) "Performing transient."
        write(*,*) 

        ! Initialize output files 
        call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years")
        call yelmo_write_init(yelmo1,file2Dsm,time_init=ts%time,units="years")
        call yelmo_regions_write(yelmo1,ts%time,init=.TRUE.,units="years")
        call yelmo_write_reg_init(yelmo1,file1D_esm,time_init=ts%time,units="years",mask=(yelmo1%bnd%mask_ice /= -1))

        if (ctl%esm_write_formatted) then
            ! Initialize output files for esm
            call yelmo_write_init(yelmo1,file2D_esm,time_init=ts%time,units="years")
            call yelmo_write_reg_init(yelmo1,file1D_esm,time_init=ts%time,units="years",mask=(yelmo1%bnd%mask_ice /= -1))
        end if 

        call timer_step(tmr,comp=1,label="initialization") 
        call timer_step(tmrs,comp=-1)

        ! == Advance timesteps ===

        call tstep_print_header(ts)
        
        do while (.not. ts%is_finished)

            ! == Update timestep ===

            call tstep_update(ts,ctl%dtt)
            call tstep_print(ts)
            call timer_step(tmrs,comp=0) 
            
            ! == ISOSTASY and SEA LEVEL ===========================================
            ! jablasco: avoid to update bedrock to check 16km error
            if (.True.) then
                call bsl_update(bsl, ts%time_rel)
                call isos_update(isos1, yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
                yelmo1%bnd%z_bed = isos1%out%z_bed
                yelmo1%bnd%z_sl  = isos1%out%z_ss
            end if
            call timer_step(tmrs,comp=1,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="isostasy") 

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,ts%time)

            call timer_step(tmrs,comp=2,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="yelmo") 

            ! == CLIMATE and OCEAN ==========================================

            ! Get ESM climate and ocean forcing
            call calc_climate_esm(smbpal1,mshlf1,esm1,yelmo1,ctl,ts%time,ts%time_rel, &
                                  domain=yelmo1%par%domain,grid_name=yelmo1%par%grid_name)
            
            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf   

            call timer_step(tmrs,comp=3,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="climate") 

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_2Dsm,ts%time)) then
                call write_step_2D_small(yelmo1,isos1,esm1,mshlf1,smbpal1,file2Dsm,ts%time)
            end if

            if (timeout_check(tm_2D,ts%time)) then
                call write_step_2D_combined(yelmo1,isos1,esm1,mshlf1,smbpal1,file2D,ts%time)
            end if
           
             
            if (timeout_check(tm_1D,ts%time)) then
                call yelmo_regions_write(yelmo1,ts%time)
                call write_1D_esm(yelmo1,esm1,mshlf1,file1D_esm,ts%time)
            end if 

            ! esm output if desired:
            if (ctl%esm_write_formatted) then
                if (mod(nint(ts%time_elapsed*100),nint(ctl%esm_dt_formatted*100))==0) then
                    call write_step_2D_esm(yelmo1,file2D_esm,ts%time)
                    call write_1D_esm(yelmo1,esm1,mshlf1,file1D_esm,ts%time)
                end if
            end if

            call timer_step(tmrs,comp=4,time_mod=[ts%time-ctl%dtt,ts%time]*1e-3,label="io") 
        
            if (mod(ts%time_elapsed,10.0)==0) then
                ! Print timestep timing info and write log table
                call timer_write_table(tmrs,[ts%time,ctl%dtt]*1e-3,"m",tmr_file,init=ts%time_elapsed .eq. 0.0)
            end if 

            if (mod(ts%time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", ts%time
            end if 
            
        end do 

        write(*,*)
        write(*,*) "Transient complete."
        write(*,*)

        ! Write the restart snapshot for the end of the simulation
        call yelmox_restart_write(bsl,isos1,yelmo1,mshlf1,ts%time)

    end select

    ! Deallocate the seed array
    deallocate(seed) 

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=ts%time*1e-3)
    
contains

    ! === CLIMATE ====
    subroutine calc_climate_esm(smbp,mshlf,esm,ylmo,ctl,time,time_bp,domain,grid_name,use_ref_atm,use_ref_ocn)

        implicit none 
    
        type(smbpal_class),         intent(INOUT) :: smbp
        type(marshelf_class),       intent(INOUT) :: mshlf 
        type(esm_forcing_class),    intent(INOUT) :: esm
        type(yelmo_class),          intent(IN)    :: ylmo
        type(ctrl_params),          intent(IN)    :: ctl
        real(wp),                   intent(IN)    :: time 
        real(wp),                   intent(IN)    :: time_bp
        character(len=*),           intent(IN)    :: domain
        character(len=*),           intent(IN)    :: grid_name 
        logical, optional,          intent(IN)    :: use_ref_atm,use_ref_ocn
    
        ! Step 1: set the reference climatologies (for reference and variability reference)
        call esm_clim_update(esm,ylmo%tpo%now%z_srf,time,ctl%time_ref,domain,grid_name)
        ! extrapolate toward sinterior of ice shelf
        call ocn_variable_extrapolation(esm%to_var_ref%var(:,:,:,1),ylmo%tpo%now%H_ice,ylmo%bnd%basins,-esm%to_var_ref%z,ylmo%bnd%z_bed)
        call ocn_variable_extrapolation(esm%so_var_ref%var(:,:,:,1),ylmo%tpo%now%H_ice,ylmo%bnd%basins,-esm%so_var_ref%z,ylmo%bnd%z_bed)

        ! Step 2: Calculate anomaly fields (forcing)
        call esm_forcing_update(esm,mshlf,time,ctl%esm_use_esm,ctl%time_ref,ctl%time_hist,ctl%time_proj,ctl%time_esm_ref, &
                                ylmo%tpo%now%H_ice,ylmo%bnd%basins,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd,ylmo%bnd%z_sl, &
                                use_ref_atm=.false.,use_ref_ocn=.false.)

        ! Step 3: Calculate the variability anomaly field
        call esm_variability_update(esm,mshlf,time,ctl%dtt,ctl%clim_var,ctl%time_ref, &
                                    ylmo%tpo%now%H_ice,ylmo%bnd%basins,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd,ylmo%bnd%z_sl, &
                                    use_ref_atm=.false.,use_ref_ocn=.false.)
        
        ! Check anomalies
        if (.FALSE.) then
            ! ATM Anomaly check
            write(*,*) "dts = ",     maxval(esm%dts)
            write(*,*) "dts_var = ", maxval(esm%dts_var)
            write(*,*) "dpr = ",     maxval(esm%dpr)
            write(*,*) "dpr_var = ", maxval(esm%dpr_var)
            ! OCN Anomaly check
            write(*,*) "dto = ",     maxval(esm%dto)
            write(*,*) "dto_var = ", maxval(esm%dto_var)
            write(*,*) "dso = ",     maxval(esm%dso)
            write(*,*) "dso_var = ", maxval(esm%dso_var)
        end if
    
        ! === Atmospheric boundary conditions ===
        ! Calculate SMB fields
        call smbpal_update_monthly(smbp,esm%t2m+esm%dts+esm%dts_var,esm%pr*esm%dpr*esm%dpr_var, &
                                    ylmo%tpo%now%z_srf,ylmo%tpo%now%H_ice,time)
    
        ! === Oceanic boundary conditions ===
        ! Interpolate ocean data to the interior of the ice shelf
        call ocn_variable_extrapolation(esm%to_ref%var(:,:,:,1), ylmo%tpo%now%H_ice, ylmo%bnd%basins,-esm%to_ref%z,ylmo%bnd%z_bed)
        call ocn_variable_extrapolation(esm%so_ref%var(:,:,:,1), ylmo%tpo%now%H_ice, ylmo%bnd%basins,-esm%so_ref%z,ylmo%bnd%z_bed)
    
        if (.FALSE.) then
            ! Atm check
            write(*,*) "Ts_ref_max = ",maxval(esm%ts_ref%var(:,:,:,1))
            write(*,*) "Pr_ref_max = ",maxval(esm%pr_ref%var(:,:,:,1))
            write(*,*) "Ts_ref_min = ",minval(esm%ts_ref%var(:,:,:,1))
            write(*,*) "Pr_ref_min = ",minval(esm%pr_ref%var(:,:,:,1))
            ! Ocean check
            write(*,*) "To_max = ",maxval(esm%to_ref%var(:,:,:,1))
            write(*,*) "So_max = ",maxval(esm%so_ref%var(:,:,:,1))
            write(*,*) "To_min = ",minval(esm%to_ref%var(:,:,:,1))
            write(*,*) "So_min = ",minval(esm%so_ref%var(:,:,:,1))
        end if
    
        ! Update oceanic fields at desired levels
        call marshelf_interp_shelf(mshlf%now%T_shlf,mshlf,esm%to_ref%var(:,:,:,1),ylmo%tpo%now%H_ice, &
                                   ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd,ylmo%bnd%z_sl,-esm%to_ref%z)
        call marshelf_interp_shelf(mshlf%now%S_shlf,mshlf,esm%so_ref%var(:,:,:,1),ylmo%tpo%now%H_ice, &
                                        ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd,ylmo%bnd%z_sl,-esm%so_ref%z)
    
        ! Add anomaly to reference ocean
        mshlf%now%T_shlf = mshlf%now%T_shlf + esm%dto + esm%dto_var
        mshlf%now%S_shlf = mshlf%now%S_shlf + esm%dso + esm%dso_var
    
        ! Update bmb_shlf and mask_ocn
        call marshelf_update(mshlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd, &
                             ylmo%bnd%regions,ylmo%bnd%basins,ylmo%bnd%z_sl,dx=ylmo%grd%dx)
    
        return
    
    end subroutine calc_climate_esm
    
    ! Yelmo output routines.

    subroutine write_step_2D_combined(ylmo,isos,esm,mshlf,srf,filename,time)

        implicit none

        type(yelmo_class),       intent(IN) :: ylmo
        type(isos_class),        intent(IN) :: isos
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        type(smbpal_class),      intent(IN) :: srf

        character(len=*),       intent(IN) :: filename
        real(wp),               intent(IN) :: time

        ! Local variables
        integer  :: ncid, n

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
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_grz",ylmo%tpo%now%mask_grz,units="",long_name="Grounding-zone mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lsf",ylmo%tpo%now%lsf,units="",long_name="LSF mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mask_ice",ylmo%bnd%mask_ice,units="",long_name="Ice mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mask_frnt",ylmo%tpo%now%mask_frnt,units="",long_name="Ice-front mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taul_int_acx",ylmo%dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral stress (x)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",ylmo%dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral stress (y)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m/yr^2",long_name="Velocity rate of change", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/a ice equiv.",long_name="Calving mass balance rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/a ice equiv.",long_name="Calving mass balance rate flt", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd_bmb",ylmo%tpo%now%f_grnd_bmb,units="1",long_name="Grounded fraction (bmb)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
                      long_name="Distance to nearest grounding-line point", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal stress", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m/yr^2",long_name="Velocity rate of change", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity (z)", &
                       dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer thickness", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/a ice equiv.",long_name="Net surface mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_ref",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Net basal mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"fmb",ylmo%tpo%now%fmb,units="m/a ice equiv.",long_name="Net margin-front mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        ! External data
        call nc_write(filename,"dzbdt",isos%out%dwdt,units="m/a",long_name="Bedrock uplift rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
    
        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! ESM Atmospheric boundary fields            
        call nc_write(filename,"t2m_ann",esm%t2m_ann+SUM(esm%dts, dim=3)/12.0,units="K",long_name="Near-surface air temperature (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"t2m_sum",esm%t2m_sum+0.333*(esm%dts(:,:,12)+esm%dts(:,:,1)+esm%dts(:,:,2)),units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts_ann",SUM(esm%dts, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",esm%pr_ann*1e-3*SUM(esm%dpr, dim=3)/12.0,units="m/a water equiv.",long_name="Precipitation (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dpr_ann",SUM(esm%dpr, dim=3)/12.0,units="%",long_name="Precipitation anomaly (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Variability fields
        call nc_write(filename,"dts_var_ann",SUM(esm%dts_var, dim=3)/12.0,units="K",long_name="Annual surface air temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dpr_var_ann",SUM(esm%dpr_var, dim=3)/12.0,units="%",long_name="Annual precipitation anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Oceanic boundary conditions
        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto",esm%dto,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso",esm%dso,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto_var",esm%dto_var,units="K",long_name="Shelf temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso_var",esm%dso_var,units="PSU",long_name="Shelf salinity anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dS_shlf",mshlf%now%dS_shlf,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_fp_shlf",mshlf%now%T_fp_shlf,units="K",long_name="Shelf freezing temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mask_ocn",mshlf%now%mask_ocn,units="", &
                        long_name="Ocean mask (0: land, 1: grline, 2: fltline, 3: open ocean, 4: deep ocean, 5: lakes)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"tf_basin",mshlf%now%tf_basin,units="K",long_name="Mean basin thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr_basin",mshlf%now%tf_corr_basin,units="K",long_name="Shelf thermal forcing basin-wide correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
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
        end if

        call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        ! Strain-rate and stress tensors 
        if (.FALSE.) then

            call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Effective strain rate", &
                          dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            call nc_write(filename,"te",ylmo%mat%now%strs%te,units="Pa",long_name="Effective stress", &
                          dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"de2D",ylmo%mat%now%strn2D%de,units="yr^-1",long_name="Effective strain rate", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"div2D",ylmo%mat%now%strn2D%div,units="yr^-1",long_name="Divergence strain rate", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"te2D",ylmo%mat%now%strs2D%te,units="Pa",long_name="Effective stress", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"eps_eig_1",ylmo%mat%now%strn2D%eps_eig_1,units="1/yr",long_name="Eigen strain 1", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eig_2",ylmo%mat%now%strn2D%eps_eig_2,units="1/yr",long_name="Eigen strain 2", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eff",ylmo%tpo%now%eps_eff,units="yr^-1",long_name="Effective calving strain", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"tau_eig_1",ylmo%mat%now%strs2D%tau_eig_1,units="Pa",long_name="Eigen stress 1", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eig_2",ylmo%mat%now%strs2D%tau_eig_2,units="Pa",long_name="Eigen stress 2", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eff",ylmo%tpo%now%tau_eff,units="Pa",long_name="Effective calving stress", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step_2D_combined

    subroutine write_step_2D_small(ylmo,isos,esm,mshlf,srf,filename,time)

        implicit none

        type(yelmo_class),       intent(IN) :: ylmo
        type(isos_class),        intent(IN) :: isos
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        type(smbpal_class),      intent(IN) :: srf

        character(len=*),       intent(IN) :: filename
        real(wp),               intent(IN) :: time

        ! Local variables
        integer  :: ncid, n

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
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! === yelmo forcing ===
        ! ESM Atmospheric boundary fields            
        call nc_write(filename,"t2m_ann",esm%t2m_ann+esm%dts(:,:,1),units="K",long_name="Near-surface air temperature (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"t2m_sum",esm%t2m_sum+esm%dts(:,:,1),units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts",SUM(esm%dts, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts_var",SUM(esm%dts_var, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",esm%pr_ann*1e-3*esm%dpr(:,:,1),units="m/a water equiv.",long_name="Precipitation (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dpr_ann",SUM(esm%dpr, dim=3)/12.0,units="%",long_name="Precipitation anomaly (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dpr_var",SUM(esm%dpr_var, dim=3)/12.0,units="%",long_name="Precipitation anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Oceanic boundary conditions
        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto",esm%dto,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso",esm%dso,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto_var",esm%dto_var,units="K",long_name="Shelf temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso_var",esm%dso_var,units="PSU",long_name="Shelf salinity anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        if (.FALSE.) then
            call nc_write(filename,"so_ref",esm%so_ref%var(:,:,1,1),units="PSU",long_name="Reference oceanic salinity", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"to_ref",esm%to_ref%var(:,:,1,1),units="K",long_name="Reference oceanic temperature", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)    
            call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/a ice equiv.",long_name="Net surface mass balance", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"smb_ref",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)    
        end if

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_small

    ! ===== esm output routines =========

    subroutine write_step_2D_esm(ylmo,filename,time)

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: rho_ice 
        real(wp) :: density_corr
        real(wp) :: m3yr_to_kgs
        real(wp) :: esm_correction
        real(wp) :: yr_to_sec
        
        ! esm variables
        real(wp), allocatable :: bmb_grnd_masked(:,:)
        real(wp), allocatable :: bmb_shlf_masked(:,:)
        real(wp), allocatable :: z_base(:,:)
        real(wp), allocatable :: T_top_ice(:,:)
        real(wp), allocatable :: T_base_grnd(:,:)
        real(wp), allocatable :: T_base_flt(:,:)
        real(wp), allocatable :: flux_grl(:,:)
        
        ! Allocate and initialize local arrays
        allocate(bmb_grnd_masked(ylmo%grd%nx,ylmo%grd%ny))
        allocate(bmb_shlf_masked(ylmo%grd%nx,ylmo%grd%ny))
        allocate(z_base(ylmo%grd%nx,ylmo%grd%ny)) 
        allocate(T_top_ice(ylmo%grd%nx,ylmo%grd%ny))
        allocate(T_base_grnd(ylmo%grd%nx,ylmo%grd%ny))
        allocate(T_base_flt(ylmo%grd%nx,ylmo%grd%ny))
        allocate(flux_grl(ylmo%grd%nx,ylmo%grd%ny))

        bmb_shlf_masked   = 0.0_wp
        bmb_grnd_masked   = 0.0_wp
        z_base            = 0.0_wp
        T_top_ice         = 0.0_wp
        T_base_grnd       = 0.0_wp
        T_base_flt        = 0.0_wp 
        flux_grl          = 0.0_wp

        ! === Data conversion factors ========================================

        rho_ice             = 917.0             ! ice density kg/m3
        m3yr_to_kgs         = 3.2e-5            ! m3/yr of pure water to kg/s
        density_corr        = rho_ice/1000.0    ! ice density correction with pure water
        esm_correction   = m3yr_to_kgs*density_corr
        yr_to_sec           = 31556952.0

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)


        ! Write esm variables 

        ! == yelmo_topography ==
        call nc_write(filename,"lithk",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                        standard_name="land_ice_thickness", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"orog",ylmo%tpo%now%z_srf,units="m",long_name="Surface evelation", &
                        standard_name="surface_altitude", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! ajr, should we just be writing the Yelmo variable ylmo%tpo%now%z_base??
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0) z_base = ylmo%bnd%z_bed
        call nc_write(filename,"base",z_base,units="m",long_name="Base elevation", &
                        standard_name="base_altitude", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"topg",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                        standard_name="bedrock_altitude", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_GHF ==
        call nc_write(filename,"hfgeoubed",ylmo%bnd%Q_geo*1e3,units="W/m^2",long_name="Geothermal heat flux", &
                        standard_name="upward_geothermal_heat_flux_at_ground_level", dim1="xc",dim2="yc",start=[1,1],ncid=ncid)

        ! == yelmo_mass_balances ==
        
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0) bmb_grnd_masked = ylmo%thrm%now%bmb_grnd
        where(ylmo%tpo%now%H_ice .gt. 0.0 .and. ylmo%tpo%now%f_grnd_bmb .eq. 0.0) bmb_shlf_masked = ylmo%bnd%bmb_shlf
        call nc_write(filename,"acabf",ylmo%bnd%smb*esm_correction,units="kg m-2 s-1",long_name="Surface mass balance flux", &
                        standard_name="land_ice_surface_specific_mass_balance_flux", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"libmassbfgr",bmb_grnd_masked*esm_correction,units="kg m-2 s-1",long_name="Basal mass balance flux beneath grounded ice", &
                        standard_name="land_ice_basal_specific_mass_balance_flux", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"libmassbffl",bmb_shlf_masked*esm_correction,units="kg m-2 s-1",long_name="Basal mass balance flux beneath floating ice", &
                        standard_name="float_ice_basal_specific_mass_balance_flux", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_velocities ==
        call nc_write(filename,"dlithkdt",ylmo%tpo%now%dHidt/yr_to_sec,units="m s-1",long_name="Ice thickness imbalance", &
                        standard_name="tendency_of_land_ice_thickness", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! comparison dHdt m/yr
        !call nc_write(filename,"dHdt",ylmo%tpo%now%dHidt,units="m yr-1",long_name="Ice thickness imbalance", &
        !              standard_name="tendency_of_land_ice_thickness", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"xvelsurf",ylmo%dyn%now%ux_s/yr_to_sec,units="m s-1",long_name="Surface velocity in x", &
                        standard_name="land_ice_surface_x_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"yvelsurf",ylmo%dyn%now%uy_s/yr_to_sec,units="m s-1",long_name="Surface velocity in y", &
                        standard_name="land_ice_surface_y_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"zvelsurf",ylmo%dyn%now%uz_s/yr_to_sec,units="m s-1",long_name="Surface velocity in z", &
                        standard_name="land_ice_surface_upward_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"xvelbase",ylmo%dyn%now%ux_b/yr_to_sec,units="m s-1",long_name="Basal velocity in x", &
                        standard_name="land_ice_surface_x_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"yvelbase",ylmo%dyn%now%uy_b/yr_to_sec,units="m s-1",long_name="Basal velocity in y", &
                        standard_name="land_ice_surface_y_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"zvelbase",ylmo%dyn%now%uz_b/yr_to_sec,units="m s-1",long_name="Basal velocity in z", &
                        standard_name="land_ice_basal_upward_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"xvelmean",ylmo%dyn%now%ux_bar/yr_to_sec,units="m s-1",long_name="Mean velocity in x", &
                        standard_name="land_ice_vertical_mean_x_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"yvelmean",ylmo%dyn%now%uy_bar/yr_to_sec,units="m s-1",long_name="Mean velocity in y", &
                        standard_name="land_ice_vertical_mean_y_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! == yelmo_ice_temperature ==
        
        ! Local variables
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0) T_base_grnd = ylmo%thrm%now%T_ice(:,:,1)
        where(ylmo%tpo%now%H_ice .gt. 0.0 .and. ylmo%tpo%now%f_grnd_bmb .eq. 0.0) T_base_flt = ylmo%thrm%now%T_ice(:,:,1)
        where(ylmo%tpo%now%H_ice .gt. 0.0) T_top_ice = ylmo%thrm%now%T_ice(:,:,ylmo%dyn%par%nz_aa)
        call nc_write(filename,"litemptop",T_top_ice,units="K",long_name="Surface temperature", &
                        standard_name="temperature_at_top_of_ice_sheet_model", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"litempbotgr",T_base_grnd,units="K",long_name="Basal temperature beneath grounded ice sheet", &
                        standard_name="land_ice_basal_temperature", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"litempbotfl",T_base_flt,units="K",long_name="Basal temperature beneath floating ice shelf", &
                        standard_name="floating_ice_basal_temperature", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_drag_and_fluxes ==
        ! ajr, flux at gl not calculated correctly?
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0 .and. ylmo%tpo%now%f_grnd_bmb .lt. 1.0) flux_grl = ylmo%tpo%now%dHidt
        call nc_write(filename,"strbasemag",ylmo%dyn%now%taub,units="Pa",long_name="Basal drag", &
                        standard_name="Basal drag", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"licalvf",ylmo%tpo%now%cmb*esm_correction,units="kg m-2 s-1",long_name="Calving flux", &
                        standard_name="land_ice_specific_mass_flux_due_to_calving", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lifmassbf",(ylmo%tpo%now%fmb+ylmo%tpo%now%cmb)*esm_correction,units="kg m-2 s-1",long_name="Ice front melt and calving flux", &
                        standard_name="land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ligroundf",flux_grl*esm_correction,units="kg m-2 s-1",long_name="Grounding line flux", &
                        standard_name="land_ice_specific_mass_flux_at_grounding_line", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_masks ==
        call nc_write(filename,"sftgif",ylmo%tpo%now%f_ice,units="1",long_name="Land ice area fraction", &
                        standard_name="land_ice_area_fraction",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"sftgrf",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice sheet area fraction", &
                        standard_name="grounded_ice_sheet_area_fraction",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"sftflf",MAX(ylmo%tpo%now%f_ice-ylmo%tpo%now%f_grnd,0.0),units="1",long_name="Floating ice sheet area fraction", &
                        standard_name="floating_ice_sheet_area_fraction",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! ATM and OCN test files
        ! call nc_write(filename,"T_atm_snap",esm1%ts%var(:,:,1,1),units="K",long_name="Surface temperature anomaly", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Tf_snap",esm1%tf%var(:,:,1,1),units="K",long_name="TF map", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"T_srf",yelmo1%bnd%T_srf,units="K",long_name="Surface temperature Yelmo", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf temperature", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step_2D_esm

    subroutine write_1D_esm(dom,esm,mshlf,filename,time)

        implicit none

        type(yelmo_class), intent(IN) :: dom
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        ! Local variables
        type(yregions_class) :: reg
        
        integer  :: ncid, n
        real(wp) :: rho_ice
        real(wp) :: density_corr
        real(wp) :: m3yr_to_kgs
        real(wp) :: esm_correction
        real(wp) :: yr_to_sec
        
        real(wp) :: m3_km3
        real(wp) :: m2_km2 

        integer  :: npts_tot
        integer  :: npts_grnd
        integer  :: npts_flt
        integer  :: npts_grl
        integer  :: npts_frnt 

        real(wp) :: dx
        real(wp) :: dy 
        real(wp) :: smb_tot 
        real(wp) :: bmb_tot 
        real(wp) :: bmb_shlf_t 

        real(wp) :: A_ice_grl 
        real(wp) :: flux_grl 
        real(wp) :: A_ice_frnt 
        real(wp) :: calv_flt 
        real(wp) :: flux_frnt 
        
        ! === Climatic variables ===
        ! Atmosphere
        real(wp) :: t2m_1d
        real(wp) :: pr_1d
        real(wp) :: dt_1d,  dt_var_1d
        real(wp) :: dpr_1d, dpr_var_1d

        ! Ocean
        real(wp) :: to_1d
        real(wp) :: so_1d
        real(wp) :: tf_1d
        real(wp) :: dto_1d, dto_var_1d
        real(wp) :: dso_1d, dso_var_1d
        !real(wp) :: dtf_1d

        logical, allocatable :: mask_tot(:,:) 
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:) 
        logical, allocatable :: mask_grl(:,:) 
        logical, allocatable :: mask_frnt(:,:) 
        
        dx = dom%grd%dx 
        dy = dom%grd%dy 

        ! Allocate variables

        allocate(mask_tot(dom%grd%nx,dom%grd%ny))
        allocate(mask_grnd(dom%grd%nx,dom%grd%ny))
        allocate(mask_flt(dom%grd%nx,dom%grd%ny))
        allocate(mask_grl(dom%grd%nx,dom%grd%ny))
        allocate(mask_frnt(dom%grd%nx,dom%grd%ny))

        ! === Data conversion factors ========================================

        rho_ice             = 917.0             ! ice density kg/m3
        m3yr_to_kgs         = 3.2e-5            ! m3/yr of pure water to kg/s
        density_corr        = rho_ice/1000.0    ! ice density correction with pure water
        esm_correction      = m3yr_to_kgs*density_corr
        yr_to_sec           = 31556952.0
        
        m3_km3              = 1e-9 
        m2_km2              = 1e-6 
        
        ! 1. Determine regional values of variables 

        ! Take the global regional data object that 
        ! is calculated over the whole domain at each timestep 
        reg = dom%reg

        ! Assign masks of interest

        mask_tot  = (dom%tpo%now%H_ice .gt. 0.0) 
        mask_grnd = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%f_grnd .gt. 0.0)
        mask_flt  = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%f_grnd .eq. 0.0)
        mask_grl  = (dom%tpo%now%f_grnd_bmb .gt. 0.0 .and. dom%tpo%now%f_grnd_bmb .lt. 1.0)         
        mask_frnt = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%mask_bed .eq. 4)

        ! Determine number of points at grl and frnt
        npts_tot  = count(mask_tot)
        npts_grnd = count(mask_grnd)
        npts_flt  = count(mask_flt)
        npts_grl  = count(mask_grl)      
        npts_frnt = count(mask_frnt) 

        ! Calculate additional variables of interest for esm
        ! Atmosphere
        t2m_1d     = sum(esm%t2m_ann+esm%dts(:,:,1),     mask=mask_tot)/npts_tot
        pr_1d      = sum(esm%pr_ann*1e-3*esm%dpr(:,:,1), mask=mask_tot)/npts_tot
        dt_1d      = sum(esm%dts(:,:,1),                 mask=mask_tot)/npts_tot
        dpr_1d     = sum(100*esm%dpr(:,:,1),             mask=mask_tot)/npts_tot
        dt_var_1d  = sum(esm%dts_var(:,:,1),             mask=mask_tot)/npts_tot
        dpr_var_1d = sum(100*esm%dpr_var(:,:,1),         mask=mask_tot)/npts_tot

        ! Ocean
        to_1d      = sum(mshlf%now%T_shlf,  mask=mask_flt)/npts_flt
        so_1d      = sum(mshlf%now%S_shlf,  mask=mask_flt)/npts_flt
        tf_1d      = sum(mshlf%now%tf_shlf, mask=mask_flt)/npts_flt
        dto_1d     = sum(esm%dto,           mask=mask_flt)/npts_flt
        dso_1d     = sum(esm%dso,           mask=mask_flt)/npts_flt
        dto_var_1d = sum(esm%dto_var,       mask=mask_flt)/npts_flt
        dso_var_1d = sum(esm%dso_var,       mask=mask_flt)/npts_flt

        ! Fluxes
        smb_tot     = sum(dom%bnd%smb,    mask=mask_tot)*(dx*dy)    ! m^3/yr: flux
        bmb_tot     = sum(dom%tpo%now%bmb,mask=mask_tot)*(dx*dy)    ! m^3/yr: flux
        bmb_shlf_t  = sum(dom%tpo%now%bmb,mask=mask_flt)*(dx*dy)    ! m^3/yr: flux

        if (npts_grl .gt. 0) then

            A_ice_grl = count(dom%tpo%now%H_ice .gt. 0.0 .and. mask_grl)*dx*dy*m2_km2 ! [km^2]
            !flux_grl  = sum(dom%tpo%now%H_ice,mask=mask_grl)*(dx*dy)                 ! m^3/yr: flux
            flux_grl  = sum(dom%dyn%now%uxy_bar*dom%tpo%now%H_ice,mask=mask_grl)*dx   ! m^3/yr: flux

            ! ajr, why only *dx above? 

        else

            A_ice_grl = 0.0_wp
            flux_grl  = 0.0_wp

        end if
 
        ! ===== Frontal ice-shelves variables =====

        if (npts_frnt .gt. 0) then

            A_ice_frnt = count(dom%tpo%now%H_ice .gt. 0.0 .and. mask_frnt)*dx*dy*m2_km2         ! [km^2]
            calv_flt   = sum(dom%tpo%now%cmb_flt*dom%tpo%now%H_ice,mask=mask_frnt)*dx          ! m^3/yr: flux [m-1 yr-1]
            flux_frnt  = calv_flt+sum(dom%tpo%now%fmb*dom%tpo%now%H_ice,mask=mask_frnt)*dx      ! m^3/yr: flux [m-1 yr-1]

            ! ajr, why only *dx above? 

        else

            A_ice_frnt = 0.0_wp
            calv_flt   = 0.0_wp 
            flux_frnt  = 0.0_wp

        end if


        ! 2. Begin writing step 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"V_sl",reg%V_sl*1e-6,units="1e6 km^3",long_name="Ice volume above flotation", &
                          dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_sle",reg%V_sle,units="m sle",long_name="Sea-level equivalent volume", &
                        dim1="time",start=[n],ncid=ncid)

        ! Variability climatic fields 
        call nc_write(filename,"dt_var_1d",dt_var_1d,units="K",long_name="Mean ice surf. Temp. Anomaly", &
                standard_name="Mean ice surf. Temp. Anomaly (Variability)",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dpr_var_1d",dpr_var_1d,units="%",long_name="Mean ice surf. Pr. Anomaly", &
                standard_name="Mean ice surf. Pr. Anomaly (Variability)",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dto_var_1d",dto_var_1d,units="K",long_name="Mean ice-shelf Temp. Anomaly", &
                standard_name="Mean ice-shelf Temp. Anomaly (Variability)",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dso_var_1d",dso_var_1d,units="PSU",long_name="Mean ice-shelf draft Sal. Anomaly", &
                standard_name="Mean ice-shelf draft Sal. Anomaly (Variability)",dim1="time",start=[n],ncid=ncid)


        ! ESM climatic fields
        call nc_write(filename,"t2m_1d",t2m_1d,units="K",long_name="Mean ice surf. Temp.", &
                standard_name="Mean ice surf. Temp.",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"pr_1d",pr_1d,units="m yr-1",long_name="Mean ice surf. Pr.", &
                standard_name="Mean ice surf. Pr.",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dt_1d",dt_1d,units="K",long_name="Mean ice surf. Temp. Anomaly", &
                standard_name="Mean ice surf. Temp. Anomaly",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dpr_1d",dpr_1d,units="%",long_name="Mean ice surf. Pr. Anomaly", &
                standard_name="Mean ice surf. Pr. Anomaly",dim1="time",start=[n],ncid=ncid)
                                
        call nc_write(filename,"to_1d",to_1d,units="K",long_name="Mean ice-shelf Temp.", &
                standard_name="Mean ice-shelf draft Temp.",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"so_1d",so_1d,units="PSU",long_name="Mean ice-shelf draft Sal.", &
                standard_name="Mean ice-shelf draft Sal.",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dto_1d",dto_1d,units="K",long_name="Mean ice-shelf Temp. Anomaly", &
                standard_name="Mean ice-shelf Temp. Anomaly",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dso_1d",dso_1d,units="PSU",long_name="Mean ice-shelf draft Sal. Anomaly", &
                standard_name="Mean ice-shelf draft Sal. Anomaly",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tf_1d",tf_1d,units="K",long_name="Mean ice-shelf TF.", &
                standard_name="Mean ice-shelf draft TF",dim1="time",start=[n],ncid=ncid)
                
        call nc_write(filename,"smb_tot",smb_tot,units="m yr-1",long_name="Total SMB flux", &
                standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",bmb_shlf_t,units="m yr-1",long_name="Total BMB flux beneath floating ice", &
                standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
  

        ! ===== Mass ===== 
        call nc_write(filename,"lim",reg%V_ice*rho_ice*1e9,units="kg",long_name="Total ice mass", &
                      standard_name="land_ice_mass",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"limnsw",reg%V_sl*rho_ice*1e9,units="kg",long_name="Mass above flotation", &
                      standard_name="land_ice_mass_not_displacing_sea_water",dim1="time",start=[n],ncid=ncid)

        ! ===== Area ====
        call nc_write(filename,"iareagr",reg%A_ice_g*1e6,units="m^2",long_name="Grounded ice area", &
                      standard_name="grounded_ice_sheet_area",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"iareafl",reg%A_ice_f*1e6,units="m^2",long_name="Floating ice area", &
                      standard_name="floating_ice_shelf_area",dim1="time",start=[n],ncid=ncid)

        ! ==== Fluxes ====
        call nc_write(filename,"tendacabf",smb_tot*esm_correction,units="kg s-1",long_name="Total SMB flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlibmassbf ",bmb_tot*esm_correction,units="kg s-1",long_name="Total BMB flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlibmassbffl",bmb_shlf_t*esm_correction,units="kg s-1",long_name="Total BMB flux beneath floating ice", &
                      standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlicalvf",calv_flt*esm_correction,units="kg s-1",long_name="Total calving flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_calving",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlifmassbf",flux_frnt*esm_correction,units="kg s-1",long_name="Total calving and ice front melting flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendligroundf",flux_grl*esm_correction,units="kg s-1",long_name="Total grounding line flux", &
                      standard_name="tendency_of_grounded_ice_mass",dim1="time",start=[n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_1D_esm

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

end program yelmox_esm
