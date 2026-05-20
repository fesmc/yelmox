program yelmox

    ! == YelmoX-Bipolar, by Sergio Pérez-Montero == !
    ! Some notes:
    !   -> An additional type is built, yelmox. Contins everything that an individual yelmox domain needs to work
    !   -> The coding principle is to reproduce 1:1 the yelmox.f90 in those groups in order to be completely analogue to 1-domain yelmox executables
    !   -> This script is organized in "groups" in order to make the code more clear. 
    !      Each group is a subroutine that contains the specified lines of the original yelmox.f90 code
    !      
    !      Some exceptions are the groups that involve actions related with the coupling between domains and/or the ocean box model
    !      Groups description:
    !       1 - Parameters reading
    !       2 - Output definitions
    !       3 - Initialization:
    !           3.1 - Initialize ice sheet model
    !           3.2 - Initialize external models
    !           3.3 - Initialize boundary conditions  
    !           3.4 - Initialize output files
    !       4 - Loop:
    !           4.1 - Spinup procedure
    !           ... - Update ice sheet and external models (not grouped)
    !           4.2 - Update output

    ! Tools
    use nml 
    use ncio 
    use timestepping
    use timer
    use timeout
    use yelmo 
    use yelmo_tools, only : smooth_gauss_2D
    use ice_optimization
    use ice_sub_regions
    !use basal_dragging, only : calc_lambda_bed_lin, calc_lambda_bed_exp

    ! Yelmo libraries   (this could be removed in the near future! spm) REMOVEEE
    use yelmo_ice
    use yelmo_grid
    use yelmo_timesteps
    use yelmo_topography
    use yelmo_dynamics
    use yelmo_material
    use yelmo_thermodynamics
    use yelmo_regions
    use yelmo_boundaries
    use yelmo_data 
    use mass_conservation

    ! External libraries
    use fastisostasy    ! also reexports barysealevel
    use snapclim
    use marine_shelf
    use smbpal
    use sediments
    use geothermal

    ! Ocean Box Model (OBM) libraries
    use obm_defs
    use stommel
    use nautilus
    use obm
    use ice2ocean
    use ocean2ice
    use atm2ocean

    implicit none

    ! Type definitions

    type negis_params
        logical  :: use_negis_par
        real(wp) :: cf_0
        real(wp) :: cf_1
        real(wp) :: cf_centre
        real(wp) :: cf_north
        real(wp) :: cf_south

        real(wp) :: cf_x
        
    end type

    type(tstep_class)      :: ts

    type(bsl_class)    :: bsl                   ! Barystatic sea level object is common to every domain
    character(len=256) :: file_bsl

    type yelmox_class   ! each yelmox domain contains all these objects

        ! Related with ´yelmo´
        type(yelmo_class)      :: yelmo1
        character(len=256)     :: domain
        logical,  allocatable  :: tmp_mask(:,:)

        ! External models
        type(snapclim_class)   :: snp1
        type(marshelf_class)   :: mshlf1
        type(smbpal_class)     :: smbpal1
        type(sediments_class)  :: sed1
        type(geothermal_class) :: gthrm1
        type(isos_class)       :: isos1

        ! Output files
        character(len=256) :: file2D, file2D_small, file_restart
        character(len=256) :: file_isos
        character(len=512) :: path_lgm

        !character(len=512) :: path_par  ! domain specific nml
        character(len=512) :: path_hydro_mask
        real(wp), allocatable :: hydro_mask(:,:)

        ! Optimization
        type(ice_opt_params) :: opt  

        ! NEGIS
        type(negis_params)   :: ngs 

        ! Internal parameters
        logical  :: running_laurentide
        logical  :: laurentide_init_const_H 
        real(wp) :: laurentide_time_equil 

        logical  :: running_greenland
        logical  :: greenland_init_marine_H
        logical  :: scale_glacial_smb

    end type

    ! Run parameters/variables
    character(len=256) :: outfldr
    real(wp) :: dT_now 
    real(wp) :: dtt_now

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

        ! bipolar-specific control parameters
        logical            :: active_north, active_south
        logical            :: obm2ism, ism2obm, atm2obm, use_yelmox, active_obm
        logical            :: couple_fwf_north, couple_fwf_south
        character(len=512) :: path_par_ocn, obm_name
        character(len=512) :: hydro_mask_north, hydro_mask_south

    end type 

    ! ctrl params
    type(ctrl_params)    :: ctl

    ! OBM files
    character(len=256) :: obm_file, obm_file_restart

    ! Bipolar simulation objects
    type(yelmox_class) :: yelmox_north, yelmox_south  ! Northern and Southern hemisphere objects (each one is a yelmox domain)
    type(obm_class)    :: obm_object                         ! Ocean Box Model object
    character(len=512) :: path_par                    ! Generalistic namelist
    logical            :: hyster_on
    character(len=512) :: hyster_forcing, hyster_forcing_method, fwf_definition
    real(wp)           :: hyster_rate, hyster_positive_branch_time
    real(wp) :: at, ao 
    real(wp) :: dTa_now, dTo_now  
    
    ! #### Start YelmoX-Bipolar simulation ##### !

    ! Group 1: Determine the general parameter file and read parameters
    call read_parameters(path_par, ctl, yelmox_north, yelmox_south, tm_1D, tm_2D, tm_2Dsm)

    ! bipolar ctl parameters
    call nml_read(path_par, "ctrl", "active_north",    ctl%active_north)       ! Is the Northern Hemisphere active?
    call nml_read(path_par, "ctrl", "active_south",    ctl%active_south)       ! Is the Southern Hemisphere active?
    call nml_read(path_par, "ctrl", "active_obm",      ctl%active_obm)
    call nml_read(path_par, "ctrl", "ism2obm",         ctl%ism2obm)
    call nml_read(path_par, "ctrl", "obm2ism",         ctl%obm2ism)
    call nml_read(path_par, "ctrl", "atm2obm",         ctl%atm2obm)
    call nml_read(path_par, "ctrl", "obm_name",        ctl%obm_name)

    if (ctl%ism2obm) then
        call nml_read(path_par, "ctrl", "couple_fwf_north",        ctl%couple_fwf_north)
        call nml_read(path_par, "ctrl", "couple_fwf_south",        ctl%couple_fwf_south)
        call nml_read(path_par, "ctrl", "fwf_definition", fwf_definition)

        if (ctl%couple_fwf_north) then
            call nml_read(path_par,"ctrl","hydro_mask_north", yelmox_north%path_hydro_mask)
        end if
        if (ctl%couple_fwf_south) then
            call nml_read(path_par,"ctrl","hydro_mask_south", yelmox_south%path_hydro_mask)
        end if
    end if

    if (ctl%obm_name .eq. "nautilus") then
        call nml_read(path_par, ctl%obm_name, "hyster_on", hyster_on)
        call nml_read(path_par, ctl%obm_name, "hyster_forcing", hyster_forcing)
        call nml_read(path_par, ctl%obm_name, "hyster_forcing_method", hyster_forcing_method)
        call nml_read(path_par, ctl%obm_name, "hyster_rate", hyster_rate)
        call nml_read(path_par, ctl%obm_name, "hyster_positive_branch_time", hyster_positive_branch_time)
    end if

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
        if (ctl%active_north) then
            yelmox_north%opt%tf_basins = 0 
            call optimize_par_load(yelmox_north%opt,path_par,"opt_north")
        end if
        if (ctl%active_south) then
            yelmox_south%opt%tf_basins = 0 
            call optimize_par_load(yelmox_south%opt,path_par,"opt_south")
        end if

    end if 

    write(*,*) "yelmox_bipolar :: parameters loaded, opt parameters loaded, consistency checked"
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Group 2: Define some input and output locations ! yelmox.f90 lines 196-200
    file_bsl = trim(outfldr)//"bsl.nc"
    tmr_file = trim(outfldr)//"timer_table.txt"
    call define_input_and_output_locations(outfldr, "north", yelmox_north) ! 2D files
    call define_input_and_output_locations(outfldr, "south", yelmox_south)

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

    ! === Initialize Yelmo X domains, Group 3
    ! Initialize barysealevel model (common to all the domains)
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, year_bp=ts%time_rel)
    call bsl_write_init(bsl, file_bsl, ts%time)

    if (ctl%active_north) then
        call initialize_icesheet(path_par, "north", ts, ctl, yelmox_north, "yelmo_north", outfldr)  
        
        ! BIPOLAR: allocate and load hydrographic mask
        if (ctl%couple_fwf_north) then
            allocate(yelmox_north%hydro_mask(yelmox_north%yelmo1%grd%nx,yelmox_north%yelmo1%grd%ny))
            call nc_read(yelmox_north%path_hydro_mask,"mask",yelmox_north%hydro_mask)
        end if

        call initialize_external_models(path_par, yelmox_north,"isos_north","snap_north","smbpal_north","itm_north","marine_shelf_north","sediments_north","geothermal_north")   ! yelmox.f90 lines 396-419
        call initialize_boundary_conditions(ctl, yelmox_north, bsl, ts)   

        call initialize_output_files(outfldr, yelmox_north, ts, "north")

        write(*,*) "Northern Hemisphere domain initialized"
    end if

    if (ctl%active_south) then
        call initialize_icesheet(path_par, "south", ts, ctl, yelmox_south, "yelmo_south", outfldr)  
        
        ! BIPOLAR: allocate and load hydrographic mask
        if (ctl%couple_fwf_south) then
            allocate(yelmox_south%hydro_mask(yelmox_south%yelmo1%grd%nx,yelmox_south%yelmo1%grd%ny))
            call nc_read(yelmox_south%path_hydro_mask,"mask",yelmox_south%hydro_mask)
        end if
        
        call initialize_external_models(path_par, yelmox_south,"isos_south","snap_south","smbpal_south","itm_south","marine_shelf_south","sediments_south","geothermal_south")   ! yelmox.f90 lines 396-419
        call initialize_boundary_conditions(ctl, yelmox_south, bsl, ts)  
        
        call initialize_output_files(outfldr, yelmox_south, ts, "south")

        write(*,*) "Southern Hemisphere domain initialized"
    end if

    ! === Initialize ocean box model =====
    if (ctl%active_obm) then
        obm_file = trim(outfldr)//trim(ctl%obm_name)//".nc"
        obm_file_restart = trim(outfldr)//trim(ctl%obm_name)//"_restart.nc"

        call obm_init(obm_object, path_par, ctl%obm_name)
        call write_obm_init(obm_file, ts%time, "years")
        call write_obm_update(obm_object, obm_file, ctl%obm_name, ts%time)

    end if
    
    call timer_step(tmr,comp=1,label="initialization") 
    call timer_step(tmrs,comp=-1)

    ! ==== Begin main time loop ===== Group 4

    dtt_now = ctl%dtt
    call tstep_print_header(ts)

    do while (.not. ts%is_finished)

        ! == Update timestep ===

        call tstep_update(ts,dtt_now)
        call tstep_print(ts)

        ! Spin-up procedure - only relevant for time-time_init <= time_equil
        if (ctl%active_north) then
            call spinup_procedure(ctl, yelmox_north, ts)
        end if
        if (ctl%active_south) then
            call spinup_procedure(ctl, yelmox_south, ts)
        end if

        call timer_step(tmrs,comp=0) 

        ! == ISOSTASY and SEA LEVEL ======================================================
        call bsl_update(bsl, ts%time_rel)
        if (ctl%active_north) then
            call isos_update(yelmox_north%isos1, yelmox_north%yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmox_north%yelmo1%bnd%dzbdt_corr)
            yelmox_north%yelmo1%bnd%z_bed = yelmox_north%isos1%out%z_bed
            yelmox_north%yelmo1%bnd%z_sl  = yelmox_north%isos1%out%z_ss
        end if
        if (ctl%active_south) then
            call isos_update(yelmox_south%isos1, yelmox_south%yelmo1%tpo%now%H_ice, ts%time, bsl, dwdt_corr=yelmox_south%yelmo1%bnd%dzbdt_corr)
            yelmox_south%yelmo1%bnd%z_bed = yelmox_south%isos1%out%z_bed
            yelmox_south%yelmo1%bnd%z_sl  = yelmox_south%isos1%out%z_ss
        end if

        call timer_step(tmrs,comp=1,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="isostasy") 

        ! == Ocean Box Model ===================================================
        if (ctl%active_obm) then
            write(*,*) "Updating OBM ... time elapsed = ", ts%time_elapsed
            call obm_update(obm_object, dtt_now, ctl%obm_name)
        end if ! use_obm

        ! == ICE SHEET ===================================================
        if (ctl%active_north) then
            write(*,*) "Updating North domain ... time elapsed = ", ts%time_elapsed
            if (yelmox_north%running_greenland .and. yelmox_north%ngs%use_negis_par) then 
                ! Update cb_ref using negis parameters 
                call negis_update_cb_ref(yelmox_north%yelmo1,yelmox_north%ngs,ts%time)
            end if
            
            ! Update Yelmo
            if (ctl%with_ice_sheet .and. (.not. (ts%n .eq. 0 .and. yelmox_north%yelmo1%par%use_restart)) ) then
                call yelmo_update(yelmox_north%yelmo1,ts%time)
            end if
        end if
        if (ctl%active_south) then
            write(*,*) "Updating South domain ... time elapsed = ", ts%time_elapsed
            if (yelmox_south%running_greenland .and. yelmox_south%ngs%use_negis_par) then 
                ! Update cb_ref using negis parameters 
                call negis_update_cb_ref(yelmox_south%yelmo1,yelmox_south%ngs,ts%time)
            end if
            
            ! Update Yelmo
            if (ctl%with_ice_sheet .and. (.not. (ts%n .eq. 0 .and. yelmox_south%yelmo1%par%use_restart)) ) then
                call yelmo_update(yelmox_south%yelmo1,ts%time)
            end if
        end if

        call timer_step(tmrs,comp=2,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="yelmo")

        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        if (mod(nint(ts%time*100),nint(ctl%dt_clim*100))==0) then
                            ! Update snapclim

            if (ctl%active_north) then
                call snapclim_update(yelmox_north%snp1,z_srf=yelmox_north%yelmo1%tpo%now%z_srf,time=ts%time,domain=yelmox_north%domain,dx=yelmox_north%yelmo1%grd%dx,basins=yelmox_north%yelmo1%bnd%basins) 
            end if
            if (ctl%active_south) then
                call snapclim_update(yelmox_south%snp1,z_srf=yelmox_south%yelmo1%tpo%now%z_srf,time=ts%time,domain=yelmox_south%domain,dx=yelmox_south%yelmo1%grd%dx,basins=yelmox_south%yelmo1%bnd%basins) 
            end if

        end if

        ! atm2obm
        if (ctl%active_north .and. ctl%atm2obm) then
            at = interp_linear(yelmox_north%snp1%at%time, yelmox_north%snp1%at%var,xout=ts%time)
            dTa_now = at*yelmox_north%snp1%par%dTa_const
            obm_object%thetan = obm_object%par%thetan_init + dTa_now!*obm_object%par%thermal_ampl_north

            obm_object%thetat = obm_object%par%thetat_init + dTa_now*obm_object%par%thermal_ampl_tropics/obm_object%par%thermal_ampl_north
            
            ! Update vapor flux
            obm_object%phin = obm_object%par%phin_init + obm_object%par%hn*obm_object%par%pnh*dTa_now/obm_object%par%thermal_ampl_north ! hn * pnh * dTG = hn * pnh * dTn/pn
        end if
        if (ctl%active_south .and. ctl%atm2obm) then
            at = series_interp(yelmox_south%snp1%at%time, yelmox_south%snp1%at%var, ts%time)
            dTa_now = at*yelmox_south%snp1%par%dTa_const
            obm_object%thetas = obm_object%par%thetas_init + dTa_now!*obm_object%par%thermal_ampl_south

            obm_object%thetat = obm_object%par%thetat_init + dTa_now*obm_object%par%thermal_ampl_tropics/obm_object%par%thermal_ampl_south ! REDUNDANTEEE
        
            ! Update vapor flux
            obm_object%phit = obm_object%par%phit_init + obm_object%par%hs*obm_object%par%psh*dTa_now/obm_object%par%thermal_ampl_south
        end if

    ! ism2obm
        if (ctl%couple_fwf_north) then
            if (ctl%active_north) then 
                obm_object%fn = calc_fwf(yelmox_north%yelmo1%bnd%c%rho_w,yelmox_north%yelmo1%bnd%c%rho_ice, &
                                    yelmox_north%yelmo1%bnd%c%sec_year, &
                                    yelmox_north%yelmo1%tpo%now%mb_net, &
                                    yelmox_north%yelmo1%tpo%now%smb, &
                                    yelmox_north%yelmo1%tpo%now%bmb, &
                                    yelmox_north%yelmo1%tpo%now%cmb, &
                                    yelmox_north%yelmo1%tpo%now%H_ice, &
                                    yelmox_north%yelmo1%tpo%now%dHidt, &
                                    yelmox_north%yelmo1%tpo%now%f_grnd, &
                                    yelmox_north%yelmo1%tpo%par%dx,yelmox_north%yelmo1%tpo%par%dy,yelmox_north%hydro_mask,"north",fwf_definition)
            end if
        end if

        if (ctl%couple_fwf_south) then
            if (ctl%active_south) then
                obm_object%fs = calc_fwf(yelmox_south%yelmo1%bnd%c%rho_w,yelmox_south%yelmo1%bnd%c%rho_ice, &
                                    yelmox_south%yelmo1%bnd%c%sec_year, &
                                    yelmox_south%yelmo1%tpo%now%mb_net, &
                                    yelmox_south%yelmo1%tpo%now%smb, &
                                    yelmox_south%yelmo1%tpo%now%bmb, &
                                    yelmox_south%yelmo1%tpo%now%cmb, &
                                    yelmox_south%yelmo1%tpo%now%H_ice, &
                                    yelmox_south%yelmo1%tpo%now%dHidt, &
                                    yelmox_south%yelmo1%tpo%now%f_grnd, &
                                    yelmox_south%yelmo1%tpo%par%dx,yelmox_south%yelmo1%tpo%par%dy,yelmox_south%hydro_mask,"south",fwf_definition)
            end if
        end if
            
        if (ctl%obm_name .eq. "nautilus") then
            if (hyster_on) then
                call update_bipolar_hyster_forcing(ts%time, ctl%time_init, obm_object, dtt_now, hyster_positive_branch_time, hyster_rate, hyster_forcing, hyster_forcing_method)
            end if
        end if

        ! obm2ism
        if (ctl%obm2ism) then
            if (ctl%active_north) then
                ! call calc_ocean_temperature_field(yelmox_north%snp1%now%to_ann, obm_object%tn-obm_object%par%tn_init, ctl%obm_name)
                call calc_ocean_temperature_field(yelmox_north%snp1%now%to_ann, obm_object%tn, ctl%obm_name)
                !call calc_ocean_salinity_field(yelmox_north%snp1%now%so_ann, obm_object%sn-obm_object%par%sn_init)
            end if
            if (ctl%active_south) then
                ! call calc_ocean_temperature_field(yelmox_south%snp1%now%to_ann, obm_object%ts-obm_object%par%ts_init, ctl%obm_name)
                call calc_ocean_temperature_field(yelmox_south%snp1%now%to_ann, obm_object%ts, ctl%obm_name)
                !call calc_ocean_salinity_field(yelmox_south%snp1%now%so_ann, obm_object%ss-obm_object%par%ss_init)
            end if
        end if

        ! == SURFACE MASS BALANCE ==============================================

        if (ctl%active_north) then
            call smbpal_update_monthly(yelmox_north%smbpal1,yelmox_north%snp1%now%tas,yelmox_north%snp1%now%pr, &
                                       yelmox_north%yelmo1%tpo%now%z_srf,yelmox_north%yelmo1%tpo%now%H_ice,ts%time_rel) 
            yelmox_north%yelmo1%bnd%smb   = yelmox_north%smbpal1%ann%smb*yelmox_north%yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmox_north%yelmo1%bnd%T_srf = yelmox_north%smbpal1%ann%tsrf 

            if (trim(yelmox_north%yelmo1%par%domain) .eq. "Greenland" .and. yelmox_north%scale_glacial_smb) then 
                ! Modify glacial smb
                call calc_glacial_smb(yelmox_north%yelmo1%bnd%smb,yelmox_north%yelmo1%grd%lat,yelmox_north%snp1%now%ta_ann,yelmox_north%snp1%clim0%ta_ann)
            end if
        
            ! yelmox_north%yelmo1%bnd%smb   = yelmox_north%yelmo1%dta%pd%smb
            ! yelmox_north%yelmo1%bnd%T_srf = yelmox_north%yelmo1%dta%pd%t2m
            
            if (yelmox_north%running_laurentide .and. yelmox_north%laurentide_init_const_H  &
                    .and. ts%time_elapsed .lt. yelmox_north%laurentide_time_equil ) then 
                ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                where (yelmox_north%yelmo1%bnd%regions .eq. 1.1 .and. yelmox_north%yelmo1%grd%lat .gt. 50.0 .and. &
                       yelmox_north%yelmo1%bnd%z_bed .gt. 0.0 .and. yelmox_north%yelmo1%bnd%smb .lt. 0.0 ) yelmox_north%yelmo1%bnd%smb = 0.5 
            end if
        end if
        if (ctl%active_south) then
            call smbpal_update_monthly(yelmox_south%smbpal1,yelmox_south%snp1%now%tas,yelmox_south%snp1%now%pr, &
                                       yelmox_south%yelmo1%tpo%now%z_srf,yelmox_south%yelmo1%tpo%now%H_ice,ts%time_rel) 
            yelmox_south%yelmo1%bnd%smb   = yelmox_south%smbpal1%ann%smb*yelmox_south%yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmox_south%yelmo1%bnd%T_srf = yelmox_south%smbpal1%ann%tsrf 

            if (trim(yelmox_south%yelmo1%par%domain) .eq. "Greenland" .and. yelmox_south%scale_glacial_smb) then 
                ! Modify glacial smb
                call calc_glacial_smb(yelmox_south%yelmo1%bnd%smb,yelmox_south%yelmo1%grd%lat,yelmox_south%snp1%now%ta_ann,yelmox_south%snp1%clim0%ta_ann)
            end if
        
            ! yelmox_south%yelmo1%bnd%smb   = yelmox_south%yelmo1%dta%pd%smb
            ! yelmox_south%yelmo1%bnd%T_srf = yelmox_south%yelmo1%dta%pd%t2m
            
            if (yelmox_south%running_laurentide .and. yelmox_south%laurentide_init_const_H  &
                    .and. ts%time_elapsed .lt. yelmox_south%laurentide_time_equil ) then 
                ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                where (yelmox_south%yelmo1%bnd%regions .eq. 1.1 .and. yelmox_south%yelmo1%grd%lat .gt. 50.0 .and. &
                       yelmox_south%yelmo1%bnd%z_bed .gt. 0.0 .and. yelmox_south%yelmo1%bnd%smb .lt. 0.0 ) yelmox_south%yelmo1%bnd%smb = 0.5 
            end if
        end if

        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        if (ctl%active_north) then
            call marshelf_update_shelf(yelmox_north%mshlf1,yelmox_north%yelmo1%tpo%now%H_ice,yelmox_north%yelmo1%bnd%z_bed,yelmox_north%yelmo1%tpo%now%f_grnd, &
                        yelmox_north%yelmo1%bnd%basins,yelmox_north%yelmo1%bnd%z_sl,yelmox_north%yelmo1%grd%dx,yelmox_north%snp1%now%depth, &
                        yelmox_north%snp1%now%to_ann,yelmox_north%snp1%now%so_ann,dto_ann=yelmox_north%snp1%now%to_ann-yelmox_north%snp1%clim0%to_ann)

            call marshelf_update(yelmox_north%mshlf1,yelmox_north%yelmo1%tpo%now%H_ice,yelmox_north%yelmo1%bnd%z_bed,yelmox_north%yelmo1%tpo%now%f_grnd, &
                                 yelmox_north%yelmo1%bnd%regions,yelmox_north%yelmo1%bnd%basins,yelmox_north%yelmo1%bnd%z_sl,dx=yelmox_north%yelmo1%grd%dx)

            yelmox_north%yelmo1%bnd%bmb_shlf = yelmox_north%mshlf1%now%bmb_shlf
            yelmox_north%yelmo1%bnd%T_shlf   = yelmox_north%mshlf1%now%T_shlf
        end if
        if (ctl%active_south) then
            call marshelf_update_shelf(yelmox_south%mshlf1,yelmox_south%yelmo1%tpo%now%H_ice,yelmox_south%yelmo1%bnd%z_bed,yelmox_south%yelmo1%tpo%now%f_grnd, &
                        yelmox_south%yelmo1%bnd%basins,yelmox_south%yelmo1%bnd%z_sl,yelmox_south%yelmo1%grd%dx,yelmox_south%snp1%now%depth, &
                        yelmox_south%snp1%now%to_ann,yelmox_south%snp1%now%so_ann,dto_ann=yelmox_south%snp1%now%to_ann-yelmox_south%snp1%clim0%to_ann)

            call marshelf_update(yelmox_south%mshlf1,yelmox_south%yelmo1%tpo%now%H_ice,yelmox_south%yelmo1%bnd%z_bed,yelmox_south%yelmo1%tpo%now%f_grnd, &
                                 yelmox_south%yelmo1%bnd%regions,yelmox_south%yelmo1%bnd%basins,yelmox_south%yelmo1%bnd%z_sl,dx=yelmox_south%yelmo1%grd%dx)

            yelmox_south%yelmo1%bnd%bmb_shlf = yelmox_south%mshlf1%now%bmb_shlf
            yelmox_south%yelmo1%bnd%T_shlf   = yelmox_south%mshlf1%now%T_shlf
        end if

        call timer_step(tmrs,comp=3,time_mod=[ts%time-dtt_now,ts%time]*1e-3,label="climate") 

        ! == MODEL OUTPUT =======================================================
        if (ctl%active_north) then
            call update_output(yelmox_north, bsl, ts, tm_1D, tm_2D, tm_2Dsm,"north")
        end if
        if (ctl%active_south) then
            call update_output(yelmox_south, bsl, ts, tm_1D, tm_2D, tm_2Dsm,"south")
        end if

        if (ctl%active_obm) then
            if (timeout_check(tm_1D,ts%time)) then
                call write_obm_update(obm_object, obm_file, ctl%obm_name, ts%time)
            end if

            if (mod(nint(ts%time*100),nint(ctl%dt_restart*100))==0) then
                call write_obm_restart(obm_object, obm_file_restart, ts%time, "years")
            end if
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
    if (ctl%active_north) then
        call yelmox_bipolar_restart_write(bsl,yelmox_north%isos1,yelmox_north%yelmo1,yelmox_north%mshlf1,ts%time,hemisphere="north")
    end if
    if (ctl%active_south) then
        call yelmox_bipolar_restart_write(bsl,yelmox_south%isos1,yelmox_south%yelmo1,yelmox_south%mshlf1,ts%time,hemisphere="south")
    end if
    if (ctl%active_obm) then
        call write_obm_restart(obm_object, obm_file_restart, ts%time, "years")
    end if

    ! Finalize program
    if (ctl%active_north) then
        call yelmo_end(yelmox_north%yelmo1,time=ts%time)
    end if
    if (ctl%active_south) then
        call yelmo_end(yelmox_south%yelmo1,time=ts%time)
    end if

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=ts%time*1e-3)

contains

    ! Yelmo X original subroutines
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

    subroutine negis_update_cb_ref(ylmo,ngs,t)

        implicit none

        type(yelmo_class),  intent(INOUT) :: ylmo
        type(negis_params), intent(INOUT) :: ngs 
        real(wp), intent(IN) :: t 

        ! Local variables
        integer :: i, j, nx, ny 

        nx = ylmo%grd%nx 
        ny = ylmo%grd%ny 

        if (t .lt. -11e3) then 
            ngs%cf_x = ngs%cf_0
        else
            ! Linear interpolation from cf_0 to cf_1  
            ngs%cf_x = ngs%cf_0 + (t - (-11e3))/ ((0.0) - (-11e3)) * (ngs%cf_1 - ngs%cf_0)
        end if

        if (t .lt. -4e3) then
            ngs%cf_south = 1.0
        else
            ngs%cf_north = 1.0
        end if

        ! Update bed roughness coefficients cb_ref and c_bed (which are independent of velocity)
        ! like normal, using the default function defined in Yelmo:
        call calc_cb_ref(ylmo%dyn%now%cb_ref,ylmo%bnd%z_bed,ylmo%bnd%z_bed_sd,ylmo%bnd%z_sl, &
                            ylmo%bnd%H_sed,ylmo%dyn%par%till_f_sed,ylmo%dyn%par%till_sed_min, &
                            ylmo%dyn%par%till_sed_max,ylmo%dyn%par%till_cf_ref,ylmo%dyn%par%till_cf_min, &
                            ylmo%dyn%par%till_z0,ylmo%dyn%par%till_z1,ylmo%dyn%par%till_n_sd, &
                            ylmo%dyn%par%till_scale_zb,ylmo%dyn%par%till_scale_sed)

        ! === Finally, apply NEGIS scaling =============================

        do j = 1, ny 
        do i = 1, nx 

            if(ylmo%bnd%basins(i,j) .eq. 9.1) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_centre
            if(ylmo%bnd%basins(i,j) .eq. 9.2) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_south
            if(ylmo%bnd%basins(i,j) .eq. 9.3) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_north

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
        integer  :: nx, ny, nz, nk, i, j, k
        real(wp), allocatable :: w_ocn_3D(:,:,:), norm_w_ocn_3D(:,:,:), mask_ocn(:,:,:), zetas(:)
        real(wp) :: dx, dy, dz, max_w_ocn_3D, npts_tot 

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
        
        call yelmo_write_var(filename,"cmb_flt",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb_grnd",ylmo,n,ncid)

        ! call yelmo_write_var(filename,"eps_eff",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"tau_eff",ylmo,n,ncid)

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
        ! call yelmo_write_var(filename,"strn2D_de",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"strn2D_div",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"strs2D_te",ylmo,n,ncid)

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
        !call yelmo_write_var(filename,"pd_err_smb",ylmo,n,ncid)
        

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
        ! call yelmo_write_var(filename,"ux_i_bar",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"uy_i_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta_acy",ylmo,n,ncid)

        call nc_write(filename,"err_pd_smb_ref",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_strn_alt_units",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

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
        nx = size(snp%now%ta_ann,1)
        ny = size(snp%now%ta_ann,2)
        npts_tot  = real(nx*ny)

        call nc_write(filename,"meanTa_ann",sum(snp%now%ta_ann)/npts_tot,units="K",long_name="Mean near-surface air temperature (ann)", &
                      dim1="time",start=[n],ncid=ncid)

        nx = size(snp%now%to_ann,1)
        ny = size(snp%now%to_ann,2)
        nk = size(snp%now%to_ann,3)

        allocate(w_ocn_3D(nx,ny,nk))
        allocate(norm_w_ocn_3D(nx,ny,nk))
        allocate(mask_ocn(nx,ny,nk))
        allocate(zetas(nk))

        w_ocn_3D(:,:,:) = 0.0 
        norm_w_ocn_3D(:,:,:) = 0.0 

        mask_ocn = snp%now%to_ann
        zetas = snp%now%depth
        dx = ylmo%tpo%par%dx
        dy = ylmo%tpo%par%dy

        do i = 1,nx
            do j = 1,ny
                do k = 1,nk
                    dz = abs(zetas(k+1) - zetas(k)) 
                    if (mask_ocn(i, j, k) <= 273.15) then
                        w_ocn_3D(i, j, k) = 0.0
                    else
                        w_ocn_3D(i, j, k) = dx*dy*dz
                    endif
                enddo
            enddo
        enddo
        max_w_ocn_3D = maxval(w_ocn_3D)
        
        norm_w_ocn_3D = w_ocn_3D/max(max_w_ocn_3D, 1.0) ! to avoid dividing by zero
        
        ! call nc_write(filename,"w_ocn_3D1",w_ocn_3D(:,:,1),units="--",long_name="", &
        !               dim1="xc", dim2="yc", dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"w_ocn_3Dk",w_ocn_3D(:,:,nk),units="--",long_name="", &
        !               dim1="xc", dim2="yc", dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"norm_w_ocean_3D1",norm_w_ocn_3D(:,:,1),units="--",long_name="", &
        !               dim1="xc", dim2="yc", dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"norm_w_ocean_3Dk",norm_w_ocn_3D(:,:,nk),units="--",long_name="", &
        !               dim1="xc", dim2="yc", dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"meanTo_ann",sum(snp%now%to_ann*norm_w_ocn_3D)/sum(norm_w_ocn_3D),units="K",long_name="Mean ocean temperature (ann)", &
                      dim1="time",start=[n],ncid=ncid)

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

        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
!                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
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

            !call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
            !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if 

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step

    ! Bipolar specific subroutines

    function series_interp(series_time, series_var,time) result(var)
        ! Wrapper for simple `interp_linear` function
        ! for series_types. 
        implicit none 

        real(wp), dimension(:), intent(IN) :: series_time, series_var
        real(wp) :: time 
        real(wp) :: var 
        integer :: nt, i 

        ! Interpolate series object
        var = interp_linear(series_time,series_var,xout=time)

        return 

    end function series_interp

    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), intent(IN) :: xout
        real(wp) :: yout 
        integer :: i, j, n, nout 
        real(wp) :: alph

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))

            end if 
        end if 

        return 

    end function interp_linear

    elemental subroutine calc_temp_anom(temp_now,temp0,dT)
        ! Calculate current climate from clim0 plus anomaly 

        implicit none 

        real(wp), intent(OUT) :: temp_now 
        real(wp), intent(IN)  :: temp0 
        real(wp), intent(IN)  :: dT            ! Current anomaly value to impose
        
        temp_now = temp0 + dT 
        
        return

    end subroutine calc_temp_anom

    function r8_normal_ab ( a, b )

    !*****************************************************************************80
    !
    !! r8_normal_ab() returns a scaled pseudonormal R8.
    !
    !  Discussion:
    !
    !    The normal probability distribution function (PDF) is sampled,
    !    with mean A and standard deviation B.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 August 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = rk ) A, the mean of the PDF.
    !
    !    Input, real ( kind = rk ) B, the standard deviation of the PDF.
    !
    !    Output, real ( kind = rk ) R8_NORMAL_AB, a sample of the normal PDF.
    !
        implicit none
    
        integer, parameter :: rk = kind ( 1.0D+00 )
    
        real ( wp) a
        real ( wp ) b
        real ( kind = rk ) r1
        real ( kind = rk ) r2
        real ( kind = rk ) r8_normal_ab
        real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
        real ( kind = rk ) x
    
        call random_number ( harvest = r1 )
        call random_number ( harvest = r2 )
        x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
    
        r8_normal_ab = a + b * x
    
        return
    end

    subroutine update_bipolar_hyster_forcing(t, t0, ocn_object, dt, branch_time_thr, rate, forcing, forc_method)
        implicit none
        real(wp), intent(IN) :: t, t0
        type(obm_class), intent(INOUT) :: ocn_object
        real(wp) :: dt
        character(len=512) :: forcing, forc_method
        real(wp)           :: rate, branch_time_thr, factor

        select case(trim(forc_method))
            case("triangular")
                if (t <= branch_time_thr) then
                    factor = rate * dt
                else
                    factor = -1 * rate * dt
                end if
            case("sin")
                if (trim(forcing) .eq. "phit") then
                    factor = -ocn_object%phit + rate * sin(2*3.14159265358979*(t-t0)/branch_time_thr)    ! to allow negative values around 0.0
                else if (trim(forcing) .eq. "fn") then
                    factor = -ocn_object%fn + rate * sin(2*3.14159265358979*(t-t0)/branch_time_thr)    ! to allow negative values around 0.0
                end if
            case("noise")
                if (trim(forcing) .eq. "phit") then
                    factor = -ocn_object%phit + r8_normal_ab(0.0, rate)
                else if (trim(forcing) .eq. "fn") then
                    factor = -ocn_object%fn + r8_normal_ab(0.0, rate)   
                end if
            case("linear")
                factor = rate * dt
        end select

        select case(trim(forcing))
            case("phit")
                ocn_object%phit = ocn_object%phit + factor
            case("phin")
                ocn_object%phin = ocn_object%phin + factor
            case("fs")
                ocn_object%fs = ocn_object%fs + factor
            case("fn")
                ocn_object%fn = ocn_object%fn + factor
            case("global_melt")
                ocn_object%fn = ocn_object%fn + factor
                ocn_object%fs = ocn_object%fs + factor
        end select
        return

    end subroutine update_bipolar_hyster_forcing

    subroutine yelmox_bipolar_restart_write(bsl,isos,ylmo,mshlf,time,fldr,hemisphere)

        implicit none

        type(bsl_class),    intent(IN) :: bsl
        type(isos_class),   intent(IN) :: isos
        type(yelmo_class),  intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        real(wp),           intent(IN) :: time 
        character(len=*),   intent(IN), optional :: fldr
        character(len=5),   intent(IN), optional :: hemisphere

        ! Local variables
        real(wp) :: time_kyr
        character(len=32)   :: time_str
        character(len=1024) :: outfldr

        character(len=56) :: file_bsl   = "bsl_restart.nc"
        character(len=56) :: file_isos  = "isos_restart.nc"
        character(len=56) :: file_yelmo = "yelmo_restart.nc"
        character(len=56) :: file_mshlf = "marine_shelf.nc"

        file_bsl   = "bsl_"//trim(hemisphere)//"_restart.nc"
        file_isos  = "isos_"//trim(hemisphere)//"_restart.nc"
        file_mshlf = "marine_shelf_"//trim(hemisphere)//"_restart.nc"
        file_yelmo = "yelmo_"//trim(hemisphere)//"_restart.nc"
        
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

    end subroutine yelmox_bipolar_restart_write


    ! Groups:
    subroutine read_parameters(path_par, ctl, yelmox1, yelmox2, tm_1D, tm_2D, tm_2Dsm)

        implicit none

        character(len=512),  intent(OUT)   :: path_par     ! parameter namelist of the entire simulation
        type(ctrl_params),   intent(INOUT) :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox1, yelmox2
        type(timeout_class), intent(INOUT) :: tm_1D, tm_2D, tm_2Dsm

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

        return

    end subroutine read_parameters

    subroutine define_input_and_output_locations(outf, hemisphere, yelmox)

        implicit none

        character(len=256), intent(IN)    :: outf
        character(len=5), intent(IN)    :: hemisphere
        type(yelmox_class), intent(INOUT) :: yelmox

        yelmox%file2D              = trim(outf)//trim(hemisphere)//"_yelmo2D.nc"
        yelmox%file2D_small        = trim(outf)//trim(hemisphere)//"_yelmo2Dsm.nc"
        yelmox%file_isos           = trim(outf)//trim(hemisphere)//"_fastisostasy.nc"
        
        return

    end subroutine define_input_and_output_locations

    subroutine initialize_icesheet(path_par, hemisphere, ts, ctl, yelmox, group, outfldr)

        implicit none

        ! Arguments
        character(len=512),  intent(IN)    :: path_par     ! parameter namelist of the entire simulation
        character(len=5),    intent(IN)    :: hemisphere
        type(tstep_class),   intent(IN)    :: ts
        type(ctrl_params),   intent(INOUT) :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox
        character(len=*),    intent(IN)    :: group
        character(len=*),    intent(IN)    :: outfldr

        ! === Initialize ice sheet model =====

        ! Initialize data objects and load initial topography
        call yelmo_init(yelmox%yelmo1,filename=path_par,grid_def="file",time=ts%time, group=group)
        !call rewrite_common_parameters(filename=path_par, dom=yelmox%yelmo1)  ! Rewrite those parameters that are common between domains (remove in future version, spm)

        ! Store domain name as a shortcut 
        yelmox%domain = yelmox%yelmo1%par%domain  

        ! Ensure optimization fields are allocated and preassigned
        allocate(yelmox%opt%cf_min(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
        allocate(yelmox%opt%cf_max(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
        
        yelmox%opt%cf_min = yelmox%yelmo1%dyn%par%till_cf_min 
        yelmox%opt%cf_max = yelmox%yelmo1%dyn%par%till_cf_ref

        ! Define specific regions of interest =====================

        ! Shortcut switches for later use
        yelmox%running_laurentide = .FALSE. 
        yelmox%running_greenland  = .FALSE. 

        allocate(yelmox%tmp_mask(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))

        select case(trim(yelmox%domain))

            case("Antarctica")

                ! Initialize regions
                call yelmo_regions_init(yelmox%yelmo1,n=3)

                ! APIS
                call get_ice_sub_region(yelmox%tmp_mask,"APIS",yelmox%yelmo1%par%domain,yelmox%yelmo1%par%grid_name)
                call yelmo_region_init(yelmox%yelmo1%regs(1),"APIS",mask=yelmox%tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

                ! WAIS
                call get_ice_sub_region(yelmox%tmp_mask,"WAIS",yelmox%yelmo1%par%domain,yelmox%yelmo1%par%grid_name)
                call yelmo_region_init(yelmox%yelmo1%regs(2),"WAIS",mask=yelmox%tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

                ! EAIS
                call get_ice_sub_region(yelmox%tmp_mask,"EAIS",yelmox%yelmo1%par%domain,yelmox%yelmo1%par%grid_name)
                call yelmo_region_init(yelmox%yelmo1%regs(3),"EAIS",mask=yelmox%tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

            case("Laurentide")

                yelmox%running_laurentide = .TRUE. 

                yelmox%laurentide_init_const_H = .FALSE.

                if (yelmox%laurentide_init_const_H) then 
                    ! Make sure relaxation spinup is short, but transient spinup
                    ! with modified positive smb over North America is reasonably long.

                    ctl%time_equil        = 10.0 
                    yelmox%laurentide_time_equil = 5e3 

                else 
                    ! When starting from ice-6g, positive smb spinup is not necessary.
                    ! however ctl%time_equil should not be changed, as it may be useful
                    ! for spinning up thermodynamics.

                    yelmox%laurentide_time_equil = 0.0 

                end if 

                ! Make sure to set mask_ice to prevent ice from growing in
                ! Greenland (and on grid borders)

                where(abs(yelmox%yelmo1%bnd%regions - 1.30) .lt. 1e-3) yelmox%yelmo1%bnd%mask_ice = -1

                yelmox%yelmo1%bnd%mask_ice(1,:)                       = -1
                yelmox%yelmo1%bnd%mask_ice(yelmox%yelmo1%grd%nx,:)    = -1
                yelmox%yelmo1%bnd%mask_ice(:,1)                       = -1
                yelmox%yelmo1%bnd%mask_ice(:,yelmox%yelmo1%grd%ny)    = -1
                
                ! Initialize regions
                call yelmo_regions_init(yelmox%yelmo1,n=1)

                ! Hudson region
                call get_ice_sub_region(yelmox%tmp_mask,"Hudson",yelmox%yelmo1%par%domain,yelmox%yelmo1%par%grid_name)
                call yelmo_region_init(yelmox%yelmo1%regs(1),"Hudson",mask=yelmox%tmp_mask,write_to_file=.TRUE.,outfldr=outfldr)

            case("Greenland")

                yelmox%running_greenland = .TRUE.

                ! Should extra ice be imposed over continental shelf to mimic LGM state to start
                yelmox%greenland_init_marine_H = .FALSE. 
                
                ! Should glacial smb be modified to reduce negative smb values
                yelmox%scale_glacial_smb = .FALSE. 
                
                ! Should NEGIS parameter modifications be used
                yelmox%ngs%use_negis_par = .TRUE. 

                ! Load NEGIS parameters from file, if used
                if (yelmox%ngs%use_negis_par) then
                    
                    call nml_read(path_par,"negis","cf_0",       yelmox%ngs%cf_0)
                    call nml_read(path_par,"negis","cf_1",       yelmox%ngs%cf_1)
                    call nml_read(path_par,"negis","cf_centre",  yelmox%ngs%cf_centre)
                    call nml_read(path_par,"negis","cf_north",   yelmox%ngs%cf_north)
                    call nml_read(path_par,"negis","cf_north",   yelmox%ngs%cf_north)
                    
                end if 

                ! Make sure to set mask_ice to prevent ice from growing in
                ! Iceland and Svaalbard (on grid borders)

                where(abs(yelmox%yelmo1%bnd%regions - 1.20) .lt. 1e-3) yelmox%yelmo1%bnd%mask_ice = -1
                where(abs(yelmox%yelmo1%bnd%regions - 1.23) .lt. 1e-3) yelmox%yelmo1%bnd%mask_ice = -1
                where(abs(yelmox%yelmo1%bnd%regions - 1.31) .lt. 1e-3) yelmox%yelmo1%bnd%mask_ice = -1
                
                if (yelmox%yelmo1%dyn%par%till_method .eq. -1) then 
                    ! Initialize cb_ref to constant value to start
                    ! (will be updated in the time loop)

                    yelmox%yelmo1%dyn%now%cb_ref = yelmox%yelmo1%dyn%par%till_cf_ref 

                end if 

        end select

        return

    end subroutine initialize_icesheet

    subroutine initialize_external_models(path_par, yelmox, isos_group, snapclim_group,smbpal_group,itm_group, marshelf_group, sediments_group, geothermal_group)

        implicit none

        ! Arguments
        character(len=512),  intent(IN)    :: path_par     ! parameter namelist of the specific domain
        type(yelmox_class),  intent(INOUT) :: yelmox
        character(len=*),    intent(IN)    :: isos_group, snapclim_group ,smbpal_group, itm_group
        character(len=*),    intent(IN)    :: marshelf_group, sediments_group, geothermal_group

        ! === Initialize external models (forcing for ice sheet) ======

        ! Initialize fastisosaty
        call isos_init(yelmox%isos1, path_par, isos_group, yelmox%yelmo1%grd%nx, yelmox%yelmo1%grd%ny, &
            yelmox%yelmo1%grd%dx, yelmox%yelmo1%grd%dy)

        ! Initialize "climate" model (climate and ocean forcing)
        call snapclim_init(yelmox%snp1,path_par,yelmox%domain,yelmox%yelmo1%par%grid_name,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%yelmo1%bnd%basins,group=snapclim_group)
        ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
        call smbpal_init(yelmox%smbpal1,path_par,x=yelmox%yelmo1%grd%xc,y=yelmox%yelmo1%grd%yc,lats=yelmox%yelmo1%grd%lat,group=smbpal_group,itm_group=itm_group)
        
        ! Initialize marine melt model (bnd%bmb_shlf)
        call marshelf_init(yelmox%mshlf1,path_par,marshelf_group,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%domain,yelmox%yelmo1%par%grid_name,yelmox%yelmo1%bnd%regions,yelmox%yelmo1%bnd%basins)

        ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
        call sediments_init(yelmox%sed1,path_par,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%domain,yelmox%yelmo1%par%grid_name,group=sediments_group)
        yelmox%yelmo1%bnd%H_sed = yelmox%sed1%now%H 
        
        call geothermal_init(yelmox%gthrm1,path_par,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%domain,yelmox%yelmo1%par%grid_name,group=geothermal_group)
        yelmox%yelmo1%bnd%Q_geo    = yelmox%gthrm1%now%ghf 

        return

    end subroutine initialize_external_models

    subroutine initialize_boundary_conditions(ctl, yelmox, bsl, ts)

        implicit none

        ! Arguments
        type(ctrl_params),   intent(IN)    :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox
        type(bsl_class),     intent(INOUT) :: bsl
        type(tstep_class),   intent(IN)    :: ts

        ! === Update initial boundary conditions for current time and yelmo state =====
        ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

        ! Initialize the isostasy reference state using reference topography fields
        call isos_init_ref(yelmox%isos1, yelmox%yelmo1%bnd%z_bed_ref, yelmox%yelmo1%bnd%H_ice_ref)
        call isos_init_state(yelmox%isos1, yelmox%yelmo1%bnd%z_bed, yelmox%yelmo1%tpo%now%H_ice, ts%time, bsl)
        call isos_write_init_extended(yelmox%isos1, yelmox%file_isos, ts%time)

        yelmox%yelmo1%bnd%z_bed = yelmox%isos1%out%z_bed
        yelmox%yelmo1%bnd%z_sl  = yelmox%isos1%out%z_ss

        ! Update snapclim
        call snapclim_update(yelmox%snp1,z_srf=yelmox%yelmo1%tpo%now%z_srf,time=ts%time_rel,domain=yelmox%domain,dx=yelmox%yelmo1%grd%dx,basins=yelmox%yelmo1%bnd%basins)

        ! Equilibrate snowpack for itm
        if (trim(yelmox%smbpal1%par%abl_method) .eq. "itm") then 
            call smbpal_update_monthly_equil(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
                                yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,ts%time_rel,time_equil=100.0)
        end if 

        call smbpal_update_monthly(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
                                yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,ts%time_rel) 
        yelmox%yelmo1%bnd%smb   = yelmox%smbpal1%ann%smb*yelmox%yelmo1%bnd%c%conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
        yelmox%yelmo1%bnd%T_srf = yelmox%smbpal1%ann%tsrf 

        if (trim(yelmox%yelmo1%par%domain) .eq. "Greenland" .and. yelmox%scale_glacial_smb) then 
            ! Modify glacial smb
            call calc_glacial_smb(yelmox%yelmo1%bnd%smb,yelmox%yelmo1%grd%lat,yelmox%snp1%now%ta_ann,yelmox%snp1%clim0%ta_ann)
        end if

        call marshelf_update_shelf(yelmox%mshlf1,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%bnd%z_bed,yelmox%yelmo1%tpo%now%f_grnd, &
                            yelmox%yelmo1%bnd%basins,yelmox%yelmo1%bnd%z_sl,yelmox%yelmo1%grd%dx,yelmox%snp1%now%depth, &
                            yelmox%snp1%now%to_ann,yelmox%snp1%now%so_ann,dto_ann=yelmox%snp1%now%to_ann-yelmox%snp1%clim0%to_ann)

        call marshelf_update(yelmox%mshlf1,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%bnd%z_bed,yelmox%yelmo1%tpo%now%f_grnd, &
                            yelmox%yelmo1%bnd%regions,yelmox%yelmo1%bnd%basins,yelmox%yelmo1%bnd%z_sl,dx=yelmox%yelmo1%grd%dx)

        yelmox%yelmo1%bnd%bmb_shlf = yelmox%mshlf1%now%bmb_shlf
        yelmox%yelmo1%bnd%T_shlf   = yelmox%mshlf1%now%T_shlf  

        call yelmo_print_bound(yelmox%yelmo1%bnd) 
        
        ! Initialize state variables (dyn,therm,mat)
        ! (initialize temps with robin method with a cold base)
        call yelmo_init_state(yelmox%yelmo1,time=ts%time,thrm_method="robin-cold")

        ! ===== basal friction optimization ======
        !if (trim(ctl%equil_method) .eq. "opt") then  !jalv: this should be the case only for opt_cf=true
        if (yelmox%opt%opt_cf) then
            
            ! Ensure that cb_ref will be optimized (till_method == set externally) 
            yelmox%yelmo1%dyn%par%till_method = -1  

            ! If not using restart, prescribe cb_ref to initial guess 
            if (.not. yelmox%yelmo1%par%use_restart) then
            ! jalv: adding cb_target here
                if (yelmox%opt%cf_init .gt. 0.0) then
                    ! Prescribe cb_ref to initial guess
                    yelmox%yelmo1%dyn%now%cb_ref = yelmox%opt%cf_init 
                else
                    ! Load cb_ref from calculated cb_tgt field
                    yelmox%yelmo1%dyn%now%cb_ref = yelmox%yelmo1%dyn%now%cb_tgt
                end if

            end if 

        end if 
        ! ========================================

        if (.not. yelmox%yelmo1%par%use_restart) then 
            ! No restart file used, perform various initialization steps

            if (yelmox%running_laurentide) then 
                ! Start with some ice thickness for testing

                ! Load LGM reconstruction into reference ice thickness
                yelmox%path_lgm = "ice_data/Laurentide/"//trim(yelmox%yelmo1%par%grid_name)//&
                            "/"//trim(yelmox%yelmo1%par%grid_name)//"_TOPO-ICE-6G_C.nc"
                call nc_read(yelmox%path_lgm,"dz",yelmox%yelmo1%bnd%H_ice_ref,start=[1,1,1], &
                                    count=[yelmox%yelmo1%tpo%par%nx,yelmox%yelmo1%tpo%par%ny,1]) 

                ! Start with some ice cover to speed up initialization
                if (yelmox%laurentide_init_const_H) then
                
                    yelmox%yelmo1%tpo%now%H_ice = 0.0
                    where (yelmox%yelmo1%bnd%regions .eq. 1.1 .and. yelmox%yelmo1%bnd%z_bed .gt. 0.0) yelmox%yelmo1%tpo%now%H_ice = 1000.0 
                    where (yelmox%yelmo1%bnd%regions .eq. 1.12) yelmox%yelmo1%tpo%now%H_ice = 1000.0 

                else
                    ! Set LGM reconstruction as initial ice thickness over North America
                    where ( yelmox%yelmo1%bnd%z_bed .gt. -500.0 .and. &
                            (   yelmox%yelmo1%bnd%regions .eq. 1.1  .or. &
                                yelmox%yelmo1%bnd%regions .eq. 1.11 .or. &
                                yelmox%yelmo1%bnd%regions .eq. 1.12) )
                        yelmox%yelmo1%tpo%now%H_ice = yelmox%yelmo1%bnd%H_ice_ref
                    end where 

                    ! Apply Gaussian smoothing to keep things stable
                    call smooth_gauss_2D(yelmox%yelmo1%tpo%now%H_ice,dx=yelmox%yelmo1%grd%dx,f_sigma=2.0)
                
                end if 
                
                ! Load sediment mask 
                yelmox%path_lgm = "ice_data/Laurentide/"//trim(yelmox%yelmo1%par%grid_name)//&
                            "/"//trim(yelmox%yelmo1%par%grid_name)//"_SED-L97.nc"
                call nc_read(yelmox%path_lgm,"z_sed",yelmox%yelmo1%bnd%H_sed) 

                if (ctl%with_ice_sheet) then
                    ! Run Yelmo for briefly to update surface topography
                    call yelmo_update_equil(yelmox%yelmo1,ts%time,time_tot=1.0_prec,dt=1.0,topo_fixed=.TRUE.)

                    ! Addtional cleanup - remove floating ice 
                    where( yelmox%yelmo1%tpo%now%mask_bed .eq. 5) yelmox%yelmo1%tpo%now%H_ice = 0.0 
                    call yelmo_update_equil(yelmox%yelmo1,ts%time,time_tot=1.0_prec,dt=1.0,topo_fixed=.TRUE.)
                end if 

                ! Update snapclim to reflect new topography 
                call snapclim_update(yelmox%snp1,z_srf=yelmox%yelmo1%tpo%now%z_srf,time=ts%time,domain=yelmox%domain,dx=yelmox%yelmo1%grd%dx,basins=yelmox%yelmo1%bnd%basins)

                ! Update smbpal
                call smbpal_update_monthly(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
                                        yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,ts%time_rel) 
                yelmox%yelmo1%bnd%smb   = yelmox%smbpal1%ann%smb*yelmox%yelmo1%bnd%c%conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
                yelmox%yelmo1%bnd%T_srf = yelmox%smbpal1%ann%tsrf 

                if (yelmox%laurentide_init_const_H) then
                    ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                    ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                    where (yelmox%yelmo1%bnd%regions .eq. 1.1 .and. yelmox%yelmo1%grd%lat .gt. 50.0 .and. &
                            yelmox%yelmo1%bnd%z_bed .gt. 0.0 .and. yelmox%yelmo1%bnd%smb .lt. 0.0 ) yelmox%yelmo1%bnd%smb = 0.5 
                    
                    if (ctl%with_ice_sheet) then
                        ! Run yelmo for several years to ensure stable central ice dome
                        call yelmo_update_equil(yelmox%yelmo1,ts%time,time_tot=5e3,dt=5.0,topo_fixed=.FALSE.)
                    end if 

                else 

                    if (ctl%with_ice_sheet) then
                        ! Run yelmo for several years with constant boundary conditions to stabilize fields
                        call yelmo_update_equil(yelmox%yelmo1,ts%time,time_tot=1e2,dt=5.0,topo_fixed=.FALSE.)
                    end if 

                end if 

            else if (yelmox%running_greenland) then
                ! Special start-up steps for Greenland

                if (yelmox%ngs%use_negis_par) then 

                    ! Ensure till method is correct, since we are updating cb_ref externally
                    yelmox%yelmo1%dyn%par%till_method = -1 

                    ! Update cb_ref using negis parameters 
                    call negis_update_cb_ref(yelmox%yelmo1,yelmox%ngs,ts%time)

                end if 

                if (yelmox%greenland_init_marine_H) then
                    ! Add extra ice-thickness over continental shelf to start with
                    ! an LGM-like state

                    where(yelmox%yelmo1%bnd%mask_ice /= -1 .and. yelmox%yelmo1%tpo%now%H_ice .lt. 600.0 &
                            .and. yelmox%yelmo1%bnd%z_bed .gt. -500.0)

                            yelmox%yelmo1%tpo%now%H_ice = 800.0 

                    end where

                    if (ctl%with_ice_sheet) then
                        ! Run yelmo for a few years with constant boundary conditions
                        ! to synchronize all model fields a bit
                        call yelmo_update_equil(yelmox%yelmo1,ts%time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)
                    end if 

                end if
                    
            else 
                ! Run simple startup equilibration step 
                
                if (ctl%with_ice_sheet) then
                    ! Run yelmo for a few years with constant boundary conditions
                    ! to synchronize all model fields a bit
                    call yelmo_update_equil(yelmox%yelmo1,ts%time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)
                end if 

            end if
            
        end if 

        return

    end subroutine initialize_boundary_conditions

    subroutine initialize_output_files(outf, yelmox, ts, hemisphere)

        implicit none

        ! Arguments
        character(len=256), intent(IN)    :: outf
        type(yelmox_class),  intent(INOUT) :: yelmox
        type(tstep_class),   intent(IN)    :: ts
        character(len=5), intent(IN)    :: hemisphere

        call yelmo_write_init(yelmox%yelmo1,yelmox%file2D,time_init=ts%time,units="years") 
        call yelmo_write_init(yelmox%yelmo1,yelmox%file2D_small,time_init=ts%time,units="years") 

        yelmox%yelmo1%reg%fnm      = trim(outf)//trim(hemisphere)//"_yelmo1D.nc"
        call yelmo_regions_write(yelmox%yelmo1,ts%time,init=.TRUE.,units="years")

        return

    end subroutine initialize_output_files

    subroutine spinup_procedure(ctl, yelmox, ts)

        implicit none

        ! Arguments
        type(ctrl_params),   intent(INOUT) :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox
        type(tstep_class),   intent(IN)    :: ts

        select case(trim(ctl%equil_method))
            
            case("opt")

                if (ts%time_elapsed .le. yelmox%opt%rel_time2) then 
                    ! Apply relaxation to the model 

                    ! Update model relaxation time scale and error scaling (in [m])
                    call optimize_set_transient_param(yelmox%opt%rel_tau,ts%time_elapsed,time1=yelmox%opt%rel_time1,time2=yelmox%opt%rel_time2, &
                                                    p1=yelmox%opt%rel_tau1,p2=yelmox%opt%rel_tau2,m=yelmox%opt%rel_m)
                    
                    ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                    yelmox%yelmo1%tpo%par%topo_rel_tau = yelmox%opt%rel_tau 
                    yelmox%yelmo1%tpo%par%topo_rel     = 4
                
                else 
                    ! Turn-off relaxation now

                    yelmox%yelmo1%tpo%par%topo_rel = 0 

                end if 

                ! === Optimization update step =========

                if (yelmox%opt%opt_cf .and. &
                        (ts%time_elapsed .ge. yelmox%opt%cf_time_init .and. ts%time_elapsed .le. yelmox%opt%cf_time_end) ) then
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmox%yelmo1%dyn%now%cb_ref,yelmox%yelmo1%tpo%now%H_ice, &
                                                    yelmox%yelmo1%tpo%now%dHidt,yelmox%yelmo1%bnd%z_bed,yelmox%yelmo1%bnd%z_sl,yelmox%yelmo1%dyn%now%ux_s,yelmox%yelmo1%dyn%now%uy_s, &
                                                    yelmox%yelmo1%dta%pd%H_ice,yelmox%yelmo1%dta%pd%uxy_s,yelmox%yelmo1%dta%pd%H_grnd, &
                                                    yelmox%opt%cf_min,yelmox%opt%cf_max,yelmox%yelmo1%tpo%par%dx,yelmox%opt%sigma_err,yelmox%opt%sigma_vel,yelmox%opt%tau_c,yelmox%opt%H0, &
                                                    dt=ctl%dtt,fill_method=yelmox%opt%fill_method,fill_dist=yelmox%opt%sigma_err,cb_tgt=yelmox%yelmo1%dyn%now%cb_tgt)
                    
                end if

                if (yelmox%opt%opt_tf .and. &
                        (ts%time_elapsed .ge. yelmox%opt%tf_time_init .and. ts%time_elapsed .le. yelmox%opt%tf_time_end) ) then
                    ! Perform tf_corr optimization

                    call optimize_tf_corr(yelmox%mshlf1%now%tf_corr,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%tpo%now%H_grnd,yelmox%yelmo1%tpo%now%dHidt, &
                                          yelmox%yelmo1%dta%pd%H_ice,yelmox%yelmo1%dta%pd%H_grnd,yelmox%opt%H_grnd_lim,yelmox%yelmo1%bnd%basins, &
                                          yelmox%opt%basin_fill,yelmox%opt%tau_m,yelmox%opt%m_temp,yelmox%opt%tf_min,yelmox%opt%tf_max,yelmox%yelmo1%tpo%par%dx,sigma=yelmox%opt%tf_sigma,dt=ctl%dtt)
                
                end if 

            case("relax")
                ! ===== relaxation spinup ==================

                if (ts%time_elapsed .lt. ctl%time_equil) then 
                    ! Turn on relaxation for now, to let thermodynamics equilibrate
                    ! without changing the topography too much. Important when 
                    ! effective pressure = f(thermodynamics).

                    yelmox%yelmo1%tpo%par%topo_rel     = 3
                    yelmox%yelmo1%tpo%par%topo_rel_tau = 50.0 
                    write(*,*) "timelog, tau = ", yelmox%yelmo1%tpo%par%topo_rel_tau

                else if (ts%time_elapsed .eq. ctl%time_equil) then 
                    ! Disable relaxation now... 

                    yelmox%yelmo1%tpo%par%topo_rel     = 0
                    write(*,*) "timelog, relation off..."

                end if 

            case DEFAULT   ! == "none", etc

                ! Pass - do nothing 

        end select 

        return

    end subroutine spinup_procedure

    subroutine update_output(yelmox, bsl, ts, tm_1D, tm_2D, tm_2Dsm,hemisphere)

        implicit none

        type(yelmox_class),  intent(INOUT) :: yelmox
        type(bsl_class),     intent(IN)    :: bsl
        type(tstep_class),   intent(IN)    :: ts
        type(timeout_class), intent(IN)    :: tm_1D, tm_2D, tm_2Dsm
        character(len=5),   intent(IN)     :: hemisphere

        ! == MODEL OUTPUT =======================================================

        if (timeout_check(tm_2D,ts%time)) then
            call yelmox_write_step(yelmox%yelmo1,yelmox%snp1,yelmox%mshlf1,yelmox%smbpal1,yelmox%file2D,time=ts%time)
        end if

        if (timeout_check(tm_2Dsm,ts%time)) then 
            call yelmo_write_step(yelmox%yelmo1,yelmox%file2D_small,ts%time,compare_pd=.FALSE.)
        end if

        if (timeout_check(tm_1D,ts%time)) then 
            call yelmo_regions_write(yelmox%yelmo1,ts%time)
        end if 

        if (mod(nint(ts%time*100),nint(ctl%dt_restart*100))==0) then
            call yelmox_bipolar_restart_write(bsl,yelmox%isos1,yelmox%yelmo1,yelmox%mshlf1,ts%time,hemisphere=hemisphere)
        end if 

        return

    end subroutine update_output

end program yelmox
