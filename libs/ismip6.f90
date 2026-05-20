module ismip6
    ! This module contains routines that help with performing the ISMIP6 suite
    ! of experiments. 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use precision
    use nml  
    use ncio 
    use varslice

    implicit none 
    
    ! Class for holding ice-forcing data from ISMIP6 archives
    type ismip6_forcing_class
        
        ! Experiment information
        character(len=256)     :: gcm 
        character(len=256)     :: scenario
        character(len=256)     :: experiment 
        character(len=256)     :: domain 
        character(len=256)     :: grid_name 
        character(len=256)     :: ctrl_run_type
        logical                :: shlf_collapse

        ! Current state:

        ! Atmospheric fields
        type(varslice_class)   :: ts
        type(varslice_class)   :: pr
        type(varslice_class)   :: smb
        type(varslice_class)   :: dts_dz
        type(varslice_class)   :: dsmb_dz
        
        ! Oceanic fields 
        type(varslice_class)   :: to
        type(varslice_class)   :: so
        type(varslice_class)   :: tf

        ! Collapse mask
        type(varslice_class)   :: mask_shlf

        ! Resources: 

        ! General fields 
        type(varslice_class)   :: basins
        type(varslice_class)   :: z_srf

        ! Atmospheric fields
        type(varslice_class)   :: ts_ref 
        type(varslice_class)   :: pr_ref 
        type(varslice_class)   :: smb_ref

        type(varslice_class)   :: ts_hist 
        type(varslice_class)   :: pr_hist 
        type(varslice_class)   :: smb_hist

        type(varslice_class)   :: ts_proj
        type(varslice_class)   :: pr_proj
        type(varslice_class)   :: smb_proj
        type(varslice_class)   :: dts_dz_proj
        type(varslice_class)   :: dsmb_dz_proj

        ! Oceanic fields 
        type(varslice_class)   :: to_ref
        type(varslice_class)   :: so_ref
        type(varslice_class)   :: tf_ref
        type(varslice_class)   :: tf_cor

        type(varslice_class)   :: to_hist
        type(varslice_class)   :: so_hist
        type(varslice_class)   :: tf_hist

        type(varslice_class)   :: to_proj
        type(varslice_class)   :: so_proj
        type(varslice_class)   :: tf_proj

        ! jablasco
        ! Collapse mask
        type(varslice_class)   :: mask_shlf_proj
        ! iceberg mask
        real(wp), allocatable :: iceberg_mask(:,:)
        
    end type

    type ismip6_ice_var_class
        character(len=56)  :: name 
        character(len=128) :: long_name
        character(len=12)  :: var_type
        character(len=128) :: standard_name 
        character(len=128) :: units_in
        character(len=128) :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
    end type

    type ismip6_experiment_class
        character(len=56)   :: expname
        character(len=56)   :: group
        character(len=56)   :: model
        character(len=256)  :: experiment
        logical             :: shlf_collapse
        character(len=256)  :: file_suffix
    end type
        
    ! Class for holding ice output for writing to standard formats...
    type ismip6_ice_class
        type(ismip6_ice_var_class), allocatable :: vars(:)
    end type 


    private
    public :: ismip6_forcing_class
    public :: ismip6_experiment_class
    public :: ismip6_ice_class

    ! General routines
    public :: ismip6_experiment_def
    public :: ismip6_forcing_init
    public :: ismip6_forcing_update
    
    public :: ismip6_write_init
    public :: ismip6_write_step
    public :: calc_iceberg_island

    public :: nahosmip_ant_forcing_init
    public :: nahosmip_ant_forcing_update 
    
contains
    
    subroutine ismip6_experiment_def(ie,expname,filename,group,model)

        implicit none

        type(ismip6_experiment_class), intent(OUT) :: ie
        character(len=*), intent(IN) :: expname
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group 
        character(len=*), intent(IN) :: model 

        ! Save the experiment name (ctrlAE, expAE01, etc)
        ie%expname = trim(expname)
        ie%group   = trim(group)
        ie%model   = trim(model)

        ! Load parameters associated with this experiment
        call nml_read(filename,ie%expname,"experiment",   ie%experiment)
        call nml_read(filename,ie%expname,"shlf_collapse",ie%shlf_collapse)
        
        ! Define the ouput filename according to protocol
        ie%file_suffix = "AIS_"//trim(group)//"_"//trim(model)//"_"//trim(ie%expname)//".nc"

        return
        
    end subroutine ismip6_experiment_def

    subroutine ismip6_forcing_init(ism,filename,domain,grid_name,gcm,scenario,experiment,shlf_collapse)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN), optional :: gcm
        character(len=*), intent(IN), optional :: scenario
        character(len=*), intent(IN), optional :: experiment
        logical,          intent(IN), optional :: shlf_collapse

        ! Assign domain and grid_name 
        ism%domain    = trim(domain)
        ism%grid_name = trim(grid_name) 

        select case(trim(domain))

            case("Antarctica")
                
                call ismip6_ant_forcing_init(ism,filename,domain,grid_name,gcm,scenario, &
                                                                    experiment,shlf_collapse)

            case("Greenland")

                call ismip6_grl_forcing_init(ism,filename,domain,grid_name,gcm,scenario)

            case DEFAULT

                write(*,*) "ismip6_forcing_init:: Error: domain not recognized."
                write(*,*) "domain = ", trim(domain)
                stop 

        end select

        return 

    end subroutine ismip6_forcing_init

    subroutine ismip6_forcing_update(ism,time,use_ref_atm,use_ref_ocn)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        real(wp), intent(IN) :: time
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 

        select case(trim(ism%domain))

            case("Antarctica")
                
                call ismip6_ant_forcing_update(ism,time,use_ref_atm,use_ref_ocn)

            case("Greenland")

                call ismip6_grl_forcing_update(ism,time,use_ref_atm,use_ref_ocn)
            
            case DEFAULT

                write(*,*) "ismip6_forcing_init:: Error: domain not recognized."
                write(*,*) "domain = ", trim(ism%domain)
                stop 

        end select
        
        return

    end subroutine ismip6_forcing_update

    subroutine ismip6_ant_forcing_init(ism,filename,domain,grid_name,gcm,scenario,experiment,shlf_collapse)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN), optional :: gcm
        character(len=*), intent(IN), optional :: scenario
        character(len=*), intent(IN), optional :: experiment
        logical,          intent(IN), optional :: shlf_collapse

        ! Local variables 
        character(len=256) :: gcm_now
        character(len=256) :: scenario_now
        
        character(len=256) :: group_prefix 

        character(len=256) :: grp_ts_ref 
        character(len=256) :: grp_pr_ref 
        character(len=256) :: grp_smb_ref 
        character(len=256) :: grp_ts_hist 
        character(len=256) :: grp_pr_hist 
        character(len=256) :: grp_smb_hist 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
        character(len=256) :: grp_smb_proj 
        
        character(len=256) :: grp_to_ref 
        character(len=256) :: grp_so_ref 
        character(len=256) :: grp_tf_ref 
        character(len=256) :: grp_tf_cor
        character(len=256) :: grp_to_hist 
        character(len=256) :: grp_so_hist 
        character(len=256) :: grp_tf_hist 
        character(len=256) :: grp_to_proj 
        character(len=256) :: grp_so_proj 
        character(len=256) :: grp_tf_proj

        character(len=256) :: grp_mask_shlf_proj
        
        integer  :: iloc, k 
        real(wp) :: tmp
        real(wp) :: time_par_proj(4) 
        real(wp) :: time_par_proj_msk(4) 

        ! First determine whether gcm+scenario provided or experiment
        ! obtain valid values for gcm and scenario to start.
        if (present(experiment)) then 
            ! Experiment provided, obtain gcm+scenario from it 

            if (trim(experiment) .eq. "ctrl" .or. trim(experiment) .eq. "ctrl0") then 
                gcm_now      = trim(experiment)
                scenario_now = trim(experiment)
            else 
                iloc = index(experiment,"_")
                if (iloc == 0) then
                    write(error_unit,*) "ismip6_ant_forcing_init:: Error: experiment &
                    &argument must be defined as gcm_scenario."
                    write(error_unit,*) "experiment = ", trim(experiment)
                    stop
                end if

                gcm_now = experiment(1:iloc-1)
                scenario_now = experiment(iloc+1:len_trim(experiment))

            end if 

        else if (present(gcm) .and. present(scenario)) then 
            ! gcm and scenario provided go forward as usual 

            gcm_now      = trim(gcm)
            scenario_now = trim(scenario)

        else
            ! Arguments not provided
            write(error_unit,*) "ismip6_ant_forcing_init:: Error: gcm+scenario or experiment must be provided &
            &as arguments."
            stop
        end if

        ! Special case for control runs, use "NorESM1-M_RCP26-repeat"
        if (trim(scenario_now) .eq. "ctrl" .or. trim(scenario_now) .eq. "ctrl0") then
            ! Assign specific gcm+scenario for control runs

            ism%gcm        =  "NorESM1-M"
            ism%scenario   =  "RCP26-repeat"
            ism%experiment = "NorESM1-M_RCP26-repeat"
            
            ! Save control run name ("ctrl0" or "ctrl") here
            ism%ctrl_run_type = trim(scenario_now)

        else
            ! Define the current experiment characteristics

            ism%gcm        = trim(gcm_now)
            ism%scenario   = trim(scenario_now) 
            ism%experiment = trim(ism%gcm)//"_"//trim(ism%scenario) 

            ! Not a control simulation
            ism%ctrl_run_type = "none"

        end if

        ! Finally, whether shlf_collapse should be included 
        ism%shlf_collapse = .FALSE. 
        if (present(shlf_collapse)) ism%shlf_collapse = shlf_collapse

        write(*,*)
        write(*,*) "ismip6_ant_forcing_init:: summary"
        write(*,*) "ctrl_run_type: ", trim(ism%ctrl_run_type)
        write(*,*) "gcm:           ", trim(ism%gcm)
        write(*,*) "scenario:      ", trim(ism%scenario)
        write(*,*) "experiment:    ", trim(ism%experiment)
        write(*,*) "shlf_collapse: ", ism%shlf_collapse
        write(*,*) 

        select case(trim(ism%experiment))

            case("CCSM4_RCP85",                 &
                 "CESM2-WACCM_ssp585",          &
                 "CESM2-WACCM_ssp585-repeat",   &
                 "HadGEM2-ES_RCP85",            &
                 "HadGEM2-ES_RCP85-repeat",     &
                 "NorESM1-M_RCP26-repeat",      &
                 "NorESM1-M_RCP85-repeat",      &
                 "UKESM1-0-LL_ssp126",          &
                 "UKESM1-0-LL_ssp585",          &
                 "UKESM1-0-LL_ssp585-repeat")
                ! Control and scenarios use the same files 
                ! since ctrl specific forcing adapted in update step 

                group_prefix = "gcm_"

                grp_ts_ref   = trim(group_prefix)//"ts_ref"
                grp_pr_ref   = trim(group_prefix)//"pr_ref"
                grp_smb_ref  = trim(group_prefix)//"smb_ref"
                grp_ts_hist  = trim(group_prefix)//"ts_hist"
                grp_pr_hist  = trim(group_prefix)//"pr_hist"
                grp_smb_hist = trim(group_prefix)//"smb_hist"
                grp_ts_proj  = trim(group_prefix)//"ts_proj"
                grp_pr_proj  = trim(group_prefix)//"pr_proj"
                grp_smb_proj = trim(group_prefix)//"smb_proj"
                
                grp_to_ref   = "to_ref"
                grp_so_ref   = "so_ref"
                grp_tf_ref   = "tf_ref"
                grp_tf_cor   = "tf_cor"
                grp_to_hist  = trim(group_prefix)//"to_hist"
                grp_so_hist  = trim(group_prefix)//"so_hist"
                grp_tf_hist  = trim(group_prefix)//"tf_hist"
                grp_to_proj  = trim(group_prefix)//"to_proj"
                grp_so_proj  = trim(group_prefix)//"so_proj"
                grp_tf_proj  = trim(group_prefix)//"tf_proj"

                grp_mask_shlf_proj = trim(group_prefix)//"mask_shlf_proj"
                
            case DEFAULT 

                write(*,*) "ismip6_ant_forcing_init:: Error: exeriment (== gcm_scenario) not recognized."
                write(*,*) "experiment = ", trim(ism%experiment) 
                write(*,*) "gcm        = ", trim(ism%gcm) 
                write(*,*) "scenario   = ", trim(ism%scenario) 

                stop 

        end select

        ! Adjust time_par_proj manually here, since there are inconsistencies between the projection runs
        ! (some runs end in 2299, 2300, or 2301). Hard code choices here to avoid many changes
        ! in parameter file
        select case(trim(ism%experiment))

            case("CESM2-WACCM_ssp585","HadGEM2-ES_RCP85")
                ! Cases that end on year 2299

                time_par_proj     = [1995.0,2299.0,1.0,0.0]
                time_par_proj_msk = [1995.0,2300.0,1.0,0.0]

            case("UKESM1-0-LL_ssp585")
                ! Cases that end on year 2301
            
                time_par_proj     = [-1.0,-1.0,-1.0,0.0]
                time_par_proj_msk = [1995.0,2301.0,1.0,0.0]

            case DEFAULT
                ! Set negative values to time_par so that values are used directly from the file

                time_par_proj     = [-1.0,-1.0,-1.0,0.0]
                time_par_proj_msk = [-1.0,-1.0,-1.0,0.0]

        end select

        ! Note: reference fields are hard-coded in the namelist file to come from the deafult
        ! reference simulation, which is NorESM1-M_rcp26-repeat

        ! Note: for the reference and historical fields, always use the time_par from the namelist file.
        
        ! Initialize all variables from namelist entries 

        ! General fields 
        call varslice_init_nml_ismip6(ism%basins,  filename,"imbie_basins",domain,grid_name,ism%gcm,ism%scenario)
        
        ! Amospheric fields
        call varslice_init_nml_ismip6(ism%ts_ref,  filename,trim(grp_ts_ref), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%pr_ref,  filename,trim(grp_pr_ref), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_ref, filename,trim(grp_smb_ref),domain,grid_name,ism%gcm,ism%scenario)
        
        call varslice_init_nml_ismip6(ism%ts_hist, filename,trim(grp_ts_hist), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%pr_hist, filename,trim(grp_pr_hist), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_hist,filename,trim(grp_smb_hist),domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%ts_proj, filename,trim(grp_ts_proj), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%pr_proj, filename,trim(grp_pr_proj), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_proj,filename,trim(grp_smb_proj),domain,grid_name,ism%gcm,ism%scenario)

        ! Oceanic fields
        call varslice_init_nml_ismip6(ism%to_ref,  filename,trim(grp_to_ref),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%so_ref,  filename,trim(grp_so_ref),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_ref,  filename,trim(grp_tf_ref),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_cor,  filename,trim(grp_tf_cor),domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%to_hist, filename,trim(grp_to_hist),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%so_hist, filename,trim(grp_so_hist),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_hist, filename,trim(grp_tf_hist),domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%to_proj, filename,trim(grp_to_proj),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%so_proj, filename,trim(grp_so_proj),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_proj, filename,trim(grp_tf_proj),domain,grid_name,ism%gcm,ism%scenario)

        if (ism%shlf_collapse) then 
            ! Shelf collapse fields
            call varslice_init_nml_ismip6(ism%mask_shlf_proj, filename,trim(grp_mask_shlf_proj), &
                                                domain,grid_name,ism%gcm,ism%scenario)
        end if

        ! Load time-independent fields

        ! Amospheric fields 
        call varslice_update(ism%ts_ref)
        call varslice_update(ism%pr_ref)
        call varslice_update(ism%smb_ref)

        ! Oceanic fields
        call varslice_update(ism%to_ref)
        call varslice_update(ism%so_ref)
        call varslice_update(ism%tf_ref)
        call varslice_update(ism%tf_cor)

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%to_ref%var,3)
            if (count(ism%to_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(ism%to_ref%var(:,:,k,1),mask=ism%to_ref%var(:,:,k,1) .ne. mv)
                where(ism%to_ref%var(:,:,k,1) .eq. mv) 
                    ism%to_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%so_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%so_ref%var(:,:,k,1),mask=ism%so_ref%var(:,:,k,1) .ne. mv)
                where(ism%so_ref%var(:,:,k,1) .eq. mv) 
                    ism%so_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%tf_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf_ref%var(:,:,k,1),mask=ism%tf_ref%var(:,:,k,1) .ne. mv)
                where(ism%tf_ref%var(:,:,k,1) .eq. mv) 
                    ism%tf_ref%var(:,:,k,1) = tmp
                end where 
            end if
        end do

        ! Initialize iceberg_mask variable in case it is needed
        allocate(ism%iceberg_mask(size(ism%to_ref%var,1),size(ism%to_ref%var,2)))
        ism%iceberg_mask = 0.0

        return 

    end subroutine ismip6_ant_forcing_init

    subroutine ismip6_ant_forcing_update(ism,time,use_ref_atm,use_ref_ocn)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        real(wp), intent(IN) :: time
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 

        ! Local variables 
        integer  :: k 
        real(wp) :: tmp 
        character(len=56) :: slice_method 

        ! Get slices for current time

        slice_method = "extrap" 

        ! === Atmospheric fields ==================================
        
        if (trim(ism%ctrl_run_type) .eq. "ctrl0") then 
            ! For control scenario "ctrl0" with no variability at all, 
            ! override time choices and set control atm and ocn 

            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
        
        else if (time .lt. 1950) then
            ! No data available, impose climatological mean 

            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
        
        else if (time .ge. 1950 .and. time .lt. 2015) then
            ! Get data from historical forcing datasets

            call varslice_update(ism%ts_hist, [time],method=slice_method)
            call varslice_update(ism%pr_hist, [time],method=slice_method)
            call varslice_update(ism%smb_hist,[time],method=slice_method)

            ism%ts  = ism%ts_hist
            ism%pr  = ism%pr_hist
            ism%smb = ism%smb_hist

        else if (time .ge. 2015 .and. time .le. 2300) then
            ! Get data from projection forcing datasets

            call varslice_update(ism%ts_proj, [time],method=slice_method)
            call varslice_update(ism%pr_proj, [time],method=slice_method)
            call varslice_update(ism%smb_proj,[time],method=slice_method) 

            ism%ts  = ism%ts_proj
            ism%pr  = ism%pr_proj
            ism%smb = ism%smb_proj

        else ! time .gt. 2300
            ! Take the average of the last 10 years of projection forcing data

            call varslice_update(ism%ts_proj, [2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%pr_proj, [2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%smb_proj,[2290.0_wp,2300.0_wp],method="range_mean")

            ism%ts  = ism%ts_proj
            ism%pr  = ism%pr_proj
            ism%smb = ism%smb_proj

        end if

        ! === Oceanic fields ==================================
        
        if (trim(ism%ctrl_run_type) .eq. "ctrl0") then 
            ! For control scenario "ctrl0" with no variability at all, 
            ! override time choices and set control atm and ocn 

            ! Set ocn fields to reference values 
            ism%to = ism%to_ref 
            ism%so = ism%so_ref 
            ism%tf = ism%tf_ref 

        else if (time .lt. 1950) then 
            ! No data available, impose climatological mean 

            ism%to = ism%to_ref
            ism%so = ism%so_ref
            ism%tf = ism%tf_ref

        else if (time .ge. 1950 .and. time .lt. 2015) then
            ! Historic 

            ism%to = ism%to_ref
            ism%so = ism%so_ref
            ism%tf = ism%tf_ref

            ! Oceanic fields 
            call varslice_update(ism%to_hist,[time],method=slice_method)
            call varslice_update(ism%so_hist,[time],method=slice_method)
            call varslice_update(ism%tf_hist,[time],method=slice_method)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj

        else if (time .ge. 2015 .and. time .le. 2300) then
            ! Projection period 1 

            call varslice_update(ism%to_proj,[time],method=slice_method)
            call varslice_update(ism%so_proj,[time],method=slice_method)
            call varslice_update(ism%tf_proj,[time],method=slice_method)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj

        else ! time .gt. 2300
            ! Projection period 2 

            call varslice_update(ism%to_proj,[2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%so_proj,[2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%tf_proj,[2290.0_wp,2300.0_wp],method="range_mean")

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj

        end if

        ! === Mask shelf collapse ==========================

        if (ism%shlf_collapse) then 
            if (time .lt. 2015) then
                
                ! Note: the year 2000 is used for this mask, however, any value from 2000 to 2015 
                ! could be used as nothing really changes with the mask in the historical period.
                ! Retreat only begins to be seen some decades later.
                call varslice_update(ism%mask_shlf_proj,[2000.0_wp],method=slice_method)
                ism%mask_shlf = ism%mask_shlf_proj

            else if (time .ge. 2015 .and. time .le. 2300) then

                call varslice_update(ism%mask_shlf_proj,[time],method=slice_method)
                ism%mask_shlf = ism%mask_shlf_proj
            
            else

                call varslice_update(ism%mask_shlf_proj,[2300.0_wp],method=slice_method)
                ism%mask_shlf  = ism%mask_shlf_proj

            end if
        end if 

        ! === Additional calculations ======================

        if (trim(ism%ctrl_run_type) .eq. "ctrl") then 
            ! For control scenario "ctrl", override above choices and 
            ! set control atm and ocn for projection time periods.
            ! This means historical variability is allowed, but 
            ! then the climatological average is imposed for the future.
            ! See protocol at: https://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections2300-Antarctica#Initialization.2C_historical_run.2C_control_run.2C_and_projection_runs

            if (time .ge. 2015) then

                ! Set atm fields to reference values (ie, zero anomaly) 
                ism%ts  = ism%ts_ref 
                ism%pr  = ism%pr_ref 
                ism%smb = ism%smb_ref 

                ! Since atm fields are anomalies, also set actual variable to zero 
                ism%ts%var  = 0.0_wp 
                ism%pr%var  = 0.0_wp 
                ism%smb%var = 0.0_wp 
            
                ! Set ocn fields to reference values 
                ism%to = ism%to_ref 
                ism%so = ism%so_ref 
                ism%tf = ism%tf_ref 

            end if

        end if 
        
        ! If desired, override other choices and use reference atm fields
        if (present(use_ref_atm)) then 
        if (use_ref_atm) then  

            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
        end if 
        end if

        ! If desired, override other choices and use reference ocean fields
        if (present(use_ref_ocn)) then 
        if (use_ref_ocn) then  

            ism%to = ism%to_ref 
            ism%so = ism%so_ref 
            ism%tf = ism%tf_ref 

        end if 
        end if

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%to%var,3)
            if (count(ism%to%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(ism%to%var(:,:,k,1),mask=ism%to%var(:,:,k,1) .ne. mv)
                where(ism%to%var(:,:,k,1) .eq. mv) 
                    ism%to%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%so%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%so%var(:,:,k,1),mask=ism%so%var(:,:,k,1) .ne. mv)
                where(ism%so%var(:,:,k,1) .eq. mv) 
                    ism%so%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%tf%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf%var(:,:,k,1),mask=ism%tf%var(:,:,k,1) .ne. mv)
                where(ism%tf%var(:,:,k,1) .eq. mv) 
                    ism%tf%var(:,:,k,1) = tmp
                end where 
            end if
        end do
        
        return 

    end subroutine ismip6_ant_forcing_update


    subroutine ismip6_grl_forcing_init(ism,filename,domain,grid_name,gcm,scenario)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN) :: gcm
        character(len=*), intent(IN) :: scenario
        
        ! Local variables 
        character(len=256) :: group_prefix 
        character(len=256) :: grp_z_srf

        character(len=256) :: grp_ts_ref 
        character(len=256) :: grp_pr_ref 
        character(len=256) :: grp_smb_ref 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
        character(len=256) :: grp_smb_proj 
        character(len=256) :: grp_dts_dz_proj 
        character(len=256) :: grp_dsmb_dz_proj 

        character(len=256) :: grp_tf_proj 
        
        integer  :: k 
        real(wp) :: tmp 
        
        ! Define the current experiment characteristics
        ism%gcm        = trim(gcm)
        ism%scenario   = trim(scenario) 
        ism%experiment = trim(ism%gcm)//"_"//trim(ism%scenario) 

        ! Define group prefix
        group_prefix = trim(ism%gcm)//"_"//trim(ism%scenario)//"_"

        if (trim(ism%scenario) .eq. "ctrl") then 
            ! Use files from rcp85 scenario as control, just 
            ! to load something. Anomalies are set to zero in any case

            group_prefix = trim(ism%gcm)//"_"//"rcp85"//"_"

        end if 

        grp_z_srf        = trim(group_prefix)//"z_srf"
        
        grp_ts_ref       = trim(group_prefix)//"ts_ref"
        grp_smb_ref      = trim(group_prefix)//"smb_ref"
        grp_ts_proj      = trim(group_prefix)//"ts_proj"
        grp_smb_proj     = trim(group_prefix)//"smb_proj"
        grp_dts_dz_proj  = trim(group_prefix)//"dts_dz_proj"
        grp_dsmb_dz_proj = trim(group_prefix)//"dsmb_dz_proj"
        
        grp_tf_proj      = trim(group_prefix)//"tf_proj"
        
        ! Initialize all variables from namelist entries 

        ! General fields 
        ! call varslice_init_nml_ismip6(ism%basins,       filename,"imbie_basins",        domain,grid_name,ism%gcm,ism%scenario)
        
        call varslice_init_nml_ismip6(ism%z_srf,        filename,trim(grp_z_srf),       domain,grid_name,ism%gcm,ism%scenario)
        
        ! Amospheric fields
        call varslice_init_nml_ismip6(ism%ts_ref,       filename,trim(grp_ts_ref),      domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_ref,      filename,trim(grp_smb_ref),     domain,grid_name,ism%gcm,ism%scenario)
        
        call varslice_init_nml_ismip6(ism%ts_proj,      filename,trim(grp_ts_proj),     domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_proj,     filename,trim(grp_smb_proj),    domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%dts_dz_proj,  filename,trim(grp_dts_dz_proj), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%dsmb_dz_proj, filename,trim(grp_dsmb_dz_proj),domain,grid_name,ism%gcm,ism%scenario)

        ! Oceanic fields
        call varslice_init_nml_ismip6(ism%tf_proj,      filename,trim(grp_tf_proj),     domain,grid_name,ism%gcm,ism%scenario)


        ! Load time-independent fields

        ! General fields 
        call varslice_update(ism%z_srf)

        ! Amospheric fields 
        call varslice_update(ism%ts_ref)
        call varslice_update(ism%smb_ref)

        ! Set remaining variables too, just for consistency, will be set to zeros
        ism%basins      = ism%ts_ref
        ism%basins%var  = 0.0_wp
        
        ism%pr          = ism%ts_ref 
        ism%pr%var      = 0.0_wp
        ism%pr_ref      = ism%ts_ref 
        ism%pr_ref%var  = 0.0_wp 
        ism%pr_proj     = ism%ts_ref
        ism%pr_proj%var = 0.0_wp

        ism%to          = ism%ts_ref 
        ism%to%var      = 0.0_wp
        ism%to_ref      = ism%ts_ref 
        ism%to_ref%var  = 0.0_wp
        ism%to_proj     = ism%ts_ref
        ism%to_proj%var = 0.0_wp
        
        ism%so          = ism%ts_ref 
        ism%so%var      = 0.0_wp
        ism%so_ref      = ism%ts_ref 
        ism%so_ref%var  = 0.0_wp

        ism%so          = ism%ts_ref
        ism%so%var      = 0.0_wp
        ism%so_proj     = ism%ts_ref
        ism%so_proj%var = 0.0_wp

        ism%tf          = ism%tf_proj 
        ism%tf%var      = 0.0_wp
        ism%tf_ref      = ism%tf_proj 
        ism%tf_ref%var  = 0.0_wp
        ism%tf_cor      = ism%tf_proj 
        ism%tf_cor%var  = 0.0_wp

        return 

    end subroutine ismip6_grl_forcing_init


    subroutine ismip6_grl_forcing_update(ism,time,use_ref_atm,use_ref_ocn)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        real(wp), intent(IN) :: time
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 

        ! Local variables 
        integer  :: k 
        real(wp) :: tmp 
        character(len=56) :: slice_method 
        type(varslice_class) :: tf_tmp 

        ! Get slices for current time

        slice_method = "interp" 

        ! === Atmospheric and oceanic fields ==================================
        
        if (trim(ism%scenario) .eq. "ctrl") then 
            ! For control scenario, override time choices and 
            ! set control atm and ocn 

            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
            ! Same for vertical gradient fields 
            ism%dts_dz      = ism%ts_ref 
            ism%dts_dz%var  = 0.0_wp 
            ism%dsmb_dz     = ism%ts_ref 
            ism%dsmb_dz%var = 0.0_wp 
            
            ! Set ocn field anomalies to zero values too
            ism%tf      = ism%tf_ref
            ism%tf%var  = 0.0_wp 

        else if (time .lt. 1950) then 

            ! Prior to 1950, no anomalies 
            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
            ! Same for vertical gradient fields 
            ism%dts_dz      = ism%ts_ref 
            ism%dts_dz%var  = 0.0_wp 
            ism%dsmb_dz     = ism%ts_ref 
            ism%dsmb_dz%var = 0.0_wp 
            
            ! Same for ocean
            ism%tf      = ism%tf_ref
            ism%tf%var  = 0.0_wp 

        else if (time .ge. 1950 .and. time .le. 2100) then 
             
            ! Load combined hist/proj transient variable as normal 
            call varslice_update(ism%ts_proj, [time],method=slice_method)
            call varslice_update(ism%smb_proj,[time],method=slice_method)
            
            call varslice_update(ism%dts_dz_proj, [time],method=slice_method)
            call varslice_update(ism%dsmb_dz_proj,[time],method=slice_method)
            
            ism%ts  = ism%ts_proj 
            ism%smb = ism%smb_proj 
            
            ism%dts_dz  = ism%dts_dz_proj 
            ism%dsmb_dz = ism%dsmb_dz_proj 
            
            ! Calculate anomaly of tf relative to reference period
            
            ! Get historical mean tf values
            tf_tmp = ism%tf_proj 
            call varslice_update(tf_tmp,[1960.0_wp,1989.0_wp],method="range_mean")

            ! Get current tf value
            call varslice_update(ism%tf_proj,[time],method=slice_method)
            
            ! Calculate tf anomaly 
            ism%tf     = ism%tf_proj 
            ism%tf%var = ism%tf_proj%var - tf_tmp%var 
            
        else ! time .gt. 2100

            call varslice_update(ism%ts_proj, [2090.0_wp,2100.0_wp],method="range_mean")
            call varslice_update(ism%smb_proj,[2090.0_wp,2100.0_wp],method="range_mean")

            call varslice_update(ism%dts_dz_proj, [2090.0_wp,2100.0_wp],method="range_mean")
            call varslice_update(ism%dsmb_dz_proj,[2090.0_wp,2100.0_wp],method="range_mean")

            ism%ts  = ism%ts_proj
            ism%smb = ism%smb_proj
            
            ism%dts_dz  = ism%dts_dz_proj 
            ism%dsmb_dz = ism%dsmb_dz_proj 
            
            ! Calculate anomaly of tf relative to reference period
            
            ! Get historical mean tf values
            tf_tmp = ism%tf_proj 
            call varslice_update(tf_tmp,[1960.0_wp,1989.0_wp],method="range_mean")

            ! Get current tf value
            call varslice_update(ism%tf_proj,[2090.0_wp,2100.0_wp],method="range_mean")
            
            ! Calculate tf anomaly 
            ism%tf     = ism%tf_proj 
            ism%tf%var = ism%tf_proj%var - tf_tmp%var 
            
        end if

        ! === Additional calculations ======================

        ! If desired, override other choices and use reference atm fields
        if (present(use_ref_atm)) then 
        if (use_ref_atm) then  

            ism%ts  = ism%ts_ref  
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
            ! Same for vertical gradient fields 
            ism%dts_dz      = ism%ts_ref 
            ism%dts_dz%var  = 0.0_wp 
            ism%dsmb_dz     = ism%ts_ref 
            ism%dsmb_dz%var = 0.0_wp 
            
        end if 
        end if

        ! If desired, override other choices and use reference ocean fields
        if (present(use_ref_ocn)) then 
        if (use_ref_ocn) then  

            ! Set ocn field anomalies to zero values too
            ism%tf     = ism%tf_ref
            ism%tf%var = 0.0_wp 

        end if 
        end if

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%tf%var,3)
            if (count(ism%tf%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf%var(:,:,k,1),mask=ism%tf%var(:,:,k,1) .ne. mv)
                where(ism%tf%var(:,:,k,1) .eq. mv) 
                    ism%tf%var(:,:,k,1) = tmp
                end where 
            end if
        end do

        return 

    end subroutine ismip6_grl_forcing_update

    ! === varslice wrapper routines with ISMIP6 specific options ===================

    subroutine varslice_init_nml_ismip6(vs,filename,group,domain,grid_name,gcm,scenario,verbose)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class),   intent(INOUT) :: vs
        character(len=*),       intent(IN)    :: filename
        character(len=*),       intent(IN)    :: group
        character(len=*),       intent(IN), optional :: domain
        character(len=*),       intent(IN), optional :: grid_name
        character(len=*),       intent(IN)    :: gcm
        character(len=*),       intent(IN)    :: scenario
        logical,                intent(IN), optional :: verbose 
        ! Local variables 
        
        ! First load parameters from nml file 
        call varslice_par_load_ismip6(vs%par,filename,group,domain,grid_name,gcm,scenario,verbose)

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        return 

    end subroutine varslice_init_nml_ismip6

    subroutine varslice_par_load_ismip6(par,filename,group,domain,grid_name,gcm,scenario,verbose)

        type(varslice_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name   
        character(len=*), intent(IN) :: gcm
        character(len=*), intent(IN) :: scenario
        logical, optional :: verbose 
    
        ! Local variables
        logical  :: init_pars
        logical  :: print_summary 
        integer  :: i 
    
        init_pars = .FALSE.
    
        print_summary = .TRUE. 
        if (present(verbose)) print_summary = verbose 
    
        call nml_read(filename,group,"filename",       par%filename,     init=init_pars)
        call nml_read(filename,group,"name",           par%name,         init=init_pars)
        call nml_read(filename,group,"units_in",       par%units_in,     init=init_pars)
        call nml_read(filename,group,"units_out",      par%units_out,    init=init_pars)
        call nml_read(filename,group,"unit_scale",     par%unit_scale,   init=init_pars)   
        call nml_read(filename,group,"unit_offset",    par%unit_offset,  init=init_pars)   
        call nml_read(filename,group,"with_time",      par%with_time,    init=init_pars)   
        call nml_read(filename,group,"time_par",       par%time_par,     init=init_pars)   
            
        ! Parse filename as needed
        call parse_path(par%filename,domain,grid_name)
        call parse_path_ismip6(par%filename,gcm,scenario)

        ! See if multiple files are available
        call get_matching_files(par%filenames, par%filename)
            
        ! Make sure time parameters are consistent time_par=[x0,x1,dx]
        if (par%time_par(3) .eq. 0.0) par%time_par(2) = par%time_par(1) 
    
        if (par%time_par(4) .gt. 1.0) then
            par%with_time_sub = .TRUE.
        else
            par%with_time_sub = .FALSE.
        end if
    
        ! Summary 
        if (print_summary) then  
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "filename      = ", trim(par%filename)
            if (size(par%filenames,1) .gt. 1) then
                write(*,*) "filenames     = "
                do i = 1, size(par%filenames,1)
                    write(*,*) "                  ", trim(par%filenames(i))
                end do
            end if
            write(*,*) "name          = ", trim(par%name)
            write(*,*) "units_in      = ", trim(par%units_in)
            write(*,*) "units_out     = ", trim(par%units_out)
            write(*,*) "unit_scale    = ", par%unit_scale
            write(*,*) "unit_offset   = ", par%unit_offset
            write(*,*) "with_time     = ", par%with_time
            write(*,*) "with_time_sub = ", par%with_time_sub
            if (par%with_time) then
                write(*,*) "time_par    = ", par%time_par
            end if
        end if 
    
        return
    
    end subroutine varslice_par_load_ismip6
    
    subroutine parse_path_ismip6(path,gcm,scenario)

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN)    :: gcm
        character(len=*), intent(IN)    :: scenario

        call nml_replace(path,"{gcm}",        trim(gcm))
        call nml_replace(path,"{scenario}",   trim(scenario))

        return

    end subroutine parse_path_ismip6
    
    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    ! === ISMIP6 OUTPUT ROUTINES ==========

    subroutine ismip6_write_init(filename,xc,yc,time,lon,lat,area,map_name,lambda,phi)

        implicit none 

        character(len=*),   intent(IN) :: filename
        real(wp),           intent(IN) :: xc(:)
        real(wp),           intent(IN) :: yc(:)
        real(wp),           intent(IN) :: time
        real(wp),           intent(IN) :: lon(:,:)
        real(wp),           intent(IN) :: lat(:,:)
        real(wp),           intent(IN) :: area(:,:)
        character(len=*),   intent(IN) :: map_name
        real(wp),           intent(IN) :: lambda
        real(wp),           intent(IN) :: phi 

        ! Local variables 
        character(len=12) :: xnm 
        character(len=12) :: ynm 
        
        xnm = "xc"
        ynm = "yc" 

        ! === Initialize netcdf file and dimensions =========

        ! Create the netcdf file 
        call nc_create(filename)

        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename,xnm,x=xc*1e-3,units="kilometers")
        call nc_write_attr(filename,xnm,"_CoordinateAxisType","GeoX")

        call nc_write_dim(filename,ynm,x=yc*1e-3,units="kilometers")
        call nc_write_attr(filename,ynm,"_CoordinateAxisType","GeoY")
        
        ! Add time axis with current value 
        call nc_write_dim(filename,"time", x=time,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)
        
        ! Projection information 
        call nc_write_map(filename,map_name,dble(lambda),phi=dble(phi))

        ! Lat-lon information
        call nc_write(filename,"lon2D",lon,dim1=xnm,dim2=ynm,grid_mapping=map_name)
        call nc_write_attr(filename,"lon2D","_CoordinateAxisType","Lon")
        call nc_write(filename,"lat2D",lat,dim1=xnm,dim2=ynm,grid_mapping=map_name)
        call nc_write_attr(filename,"lat2D","_CoordinateAxisType","Lat")

        call nc_write(filename,"area",  area*1e-6,  dim1=xnm,dim2=ynm,grid_mapping=map_name,units="km^2")
        call nc_write_attr(filename,"area","coordinates","lat2D lon2D")
        

        return

    end subroutine ismip6_write_init


    subroutine ismip6_write_step(filename,file_nml,time)

        implicit none 

        character(len=*),   intent(IN) :: filename
        character(len=*),   intent(IN) :: file_nml 
        real(wp),           intent(IN) :: time 

        ! Local variables 
        integer    :: ncid, n
        real(wp) :: time_prev 
        type(ismip6_ice_class) :: ismp 

if (.FALSE.) then    
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

end if 
        ! Load up output variable meta information 
        call ismip6_load_ice_var_info(ismp,file_nml,verbose=.TRUE.)

        ! Proceed to calculating and writing each variable 

        ! Variables to output every five years 
        if (time - time_prev .ge. 5.0_wp) then 
            ! Five years have passed, proceed 




        end if 


        ! Next, variables that should be output every year 
        return

    end subroutine ismip6_write_step

    subroutine ismip6_load_ice_var_info(ismp,filename,verbose)

        implicit none 

        type(ismip6_ice_class), intent(OUT) :: ismp 
        character(len=*),       intent(IN)  :: filename 
        logical,                intent(IN)  :: verbose 

        ! Local variables
        integer :: n  
        integer, parameter :: n_variables = 38

        type(ismip6_ice_var_class) :: v 

        ! First initialize ismp object to hold variable meta data 
        if (allocated(ismp%vars)) deallocate(ismp%vars)
        allocate(ismp%vars(n_variables))

        ! Load individual variables by namelist group 
        call ice_var_par_load(ismp%vars(1), filename,var_name="lithk")
        call ice_var_par_load(ismp%vars(2), filename,var_name="orog")
        call ice_var_par_load(ismp%vars(3), filename,var_name="base")
        call ice_var_par_load(ismp%vars(4), filename,var_name="topg")
        call ice_var_par_load(ismp%vars(5), filename,var_name="hfgeoubed")
        call ice_var_par_load(ismp%vars(6), filename,var_name="acabf")
        call ice_var_par_load(ismp%vars(7), filename,var_name="libmassbfgr")
        call ice_var_par_load(ismp%vars(8), filename,var_name="libmassbffl")
        call ice_var_par_load(ismp%vars(9), filename,var_name="dlithkdt")
        call ice_var_par_load(ismp%vars(10),filename,var_name="xvelsurf")
        call ice_var_par_load(ismp%vars(11),filename,var_name="yvelsurf")
        call ice_var_par_load(ismp%vars(12),filename,var_name="zvelsurf")
        call ice_var_par_load(ismp%vars(13),filename,var_name="xvelbase")
        call ice_var_par_load(ismp%vars(14),filename,var_name="yvelbase")
        call ice_var_par_load(ismp%vars(15),filename,var_name="zvelbase")
        call ice_var_par_load(ismp%vars(16),filename,var_name="xvelmean")
        call ice_var_par_load(ismp%vars(17),filename,var_name="yvelmean")
        call ice_var_par_load(ismp%vars(18),filename,var_name="litemptop")
        call ice_var_par_load(ismp%vars(19),filename,var_name="litempbotgr")
        call ice_var_par_load(ismp%vars(20),filename,var_name="litempbotfl")
        call ice_var_par_load(ismp%vars(21),filename,var_name="strbasemag")
        call ice_var_par_load(ismp%vars(22),filename,var_name="licalvf")
        call ice_var_par_load(ismp%vars(23),filename,var_name="lifmassbf")
        call ice_var_par_load(ismp%vars(24),filename,var_name="lifmassbf")
        call ice_var_par_load(ismp%vars(25),filename,var_name="ligroundf")
        call ice_var_par_load(ismp%vars(26),filename,var_name="sftgif")
        call ice_var_par_load(ismp%vars(27),filename,var_name="sftgrf")
        call ice_var_par_load(ismp%vars(28),filename,var_name="sftflf")

        call ice_var_par_load(ismp%vars(29),filename,var_name="lim")
        call ice_var_par_load(ismp%vars(30),filename,var_name="limnsw")
        call ice_var_par_load(ismp%vars(31),filename,var_name="iareagr")
        call ice_var_par_load(ismp%vars(32),filename,var_name="iareafl")
        call ice_var_par_load(ismp%vars(33),filename,var_name="tendacabf")
        call ice_var_par_load(ismp%vars(34),filename,var_name="tendlibmassbf")
        call ice_var_par_load(ismp%vars(35),filename,var_name="tendlibmassbffl")
        call ice_var_par_load(ismp%vars(36),filename,var_name="tendlicalvf")
        call ice_var_par_load(ismp%vars(37),filename,var_name="tendlifmassbf")
        call ice_var_par_load(ismp%vars(38),filename,var_name="tendligroundf")

        if (verbose) then 

            ! === Print summary =========

            write(*,"(a40,a8,a65,a15)") &
                                    "Variable name",    &
                                    "Type",             &
                                    "Standard name",    &
                                    "Unit"
            
            do n = 1, n_variables 
                v = ismp%vars(n)
                write(*,"(a40,a8,a65,a15)") &
                                    trim(v%long_name),      &
                                    trim(v%var_type),       &
                                    trim(v%standard_name),  &
                                    trim(v%units_out) 
            end do 


        end if 

        return 

    end subroutine ismip6_load_ice_var_info

    subroutine ice_var_par_load(ismpv,filename,var_name)
        ! Load parmaeters associated with a given ice variable

        implicit none 

        type(ismip6_ice_var_class), intent(OUT) :: ismpv 
        character(len=*),           intent(IN)  :: filename 
        character(len=*),           intent(IN)  :: var_name

        ! Local variables 
        character(len=56) :: group 

        group = "ismip6_out_"//trim(var_name)

        call nml_read(filename,group,"name",            ismpv%name)
        call nml_read(filename,group,"long_name",       ismpv%long_name)
        call nml_read(filename,group,"var_type",        ismpv%var_type)
        call nml_read(filename,group,"standard_name",   ismpv%standard_name)
        call nml_read(filename,group,"units_in",        ismpv%units_in)
        call nml_read(filename,group,"units_out",       ismpv%units_out)
        call nml_read(filename,group,"unit_scale",      ismpv%unit_scale)
        call nml_read(filename,group,"unit_offset",     ismpv%unit_offset)

        return 

    end subroutine ice_var_par_load

    ! jablasco: subroutine for iceberg islands
    subroutine calc_iceberg_island(iceberg_mask,f_grnd,H_ice)
        ! Calculate distance to the grounding line 

        implicit none

        real(wp), intent(OUT) :: iceberg_mask(:,:) ! Iceberg island (lonely ice shelves)
        real(wp), intent(IN)  :: f_grnd(:,:)       ! [1]  Grounded grid-cell fraction 
        real(wp), intent(IN)  :: H_ice(:,:)        ! [m]  Ice thickness

        ! Local variables 
        integer  :: i, j, nx, ny, q
        integer  :: im1, ip1, jm1, jp1
        integer,  parameter :: iter_max = 1000

        real(wp), allocatable :: mask_ref(:,:)
        real(wp), allocatable :: iceberg_mask_ref(:,:)

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        allocate(mask_ref(nx,ny))
        allocate(iceberg_mask_ref(nx,ny))

        ! 0: Assign mask values. 0 = ocn; 1 = flt; 2 = grd/grl
        ! ======================

        do j = 1, ny
        do i = 1, nx

            ! Grounded point or partially floating point with floating neighbors
            if (f_grnd(i,j) .gt. 0.0) then
                mask_ref(i,j) = 2.0
            else if (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) then
                mask_ref(i,j) = 1.0
            else
                mask_ref(i,j) = 0.0
            end if

        end do
        end do
        
        ! 1. Next, determine the extension of the grounded and connected flt
        ! points ===

        ! Initialize the iceberg_mask_ref
        iceberg_mask_ref = mask_ref

        do q = 1, iter_max

            do j = 1, ny
            do i = 1, nx

            ! Get neighbor indices
            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)

            ! Grounded point or partially floating point with floating neighbors
            if (iceberg_mask_ref(i,j) .eq. 2.0 .and. iceberg_mask_ref(im1,j) .eq. 1.0) then
                iceberg_mask_ref(im1,j) = 2.0
            end if
            if (iceberg_mask_ref(i,j) .eq. 2.0 .and. iceberg_mask_ref(ip1,j) .eq. 1.0) then
                iceberg_mask_ref(ip1,j) = 2.0
            end if
            if (iceberg_mask_ref(i,j) .eq. 2.0 .and. iceberg_mask_ref(i,jm1) .eq. 1.0) then
                iceberg_mask_ref(i,jm1) = 2.0
            end if
            if (iceberg_mask_ref(i,j) .eq. 2.0 .and. iceberg_mask_ref(i,jp1) .eq. 1.0) then
                iceberg_mask_ref(i,jp1) = 2.0
            end if

            end do
            end do
        end do
        
        ! 2. Where numbers are still 1.0 -> icebergs ======================

        do j = 1, ny
        do i = 1, nx

            ! Grounded point or partially floating point with floating neighbors
            if (iceberg_mask_ref(i,j) .eq. 1.0) then
                iceberg_mask(i,j) = 1.0
            else
                iceberg_mask(i,j) = 0.0
            end if

        end do
        end do


        return

    end subroutine calc_iceberg_island










! ==== NAHOSMIP =============

    ! Two routines: nahosmip_ant_forcing_init, nahosmip_ant_forcing_update

    subroutine nahosmip_ant_forcing_init(ism,filename,domain,grid_name,gcm,scenario,experiment)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN), optional :: gcm
        character(len=*), intent(IN), optional :: scenario
        character(len=*), intent(IN), optional :: experiment

        ! Local variables 
        character(len=256) :: gcm_now
        character(len=256) :: scenario_now
        
        character(len=256) :: group_prefix 

        character(len=256) :: grp_ts_ref 
        character(len=256) :: grp_pr_ref 
        character(len=256) :: grp_smb_ref 
        character(len=256) :: grp_ts_hist 
        character(len=256) :: grp_pr_hist 
        character(len=256) :: grp_smb_hist 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
        character(len=256) :: grp_smb_proj 
        
        character(len=256) :: grp_to_ref 
        character(len=256) :: grp_so_ref 
        character(len=256) :: grp_tf_ref 
        character(len=256) :: grp_tf_cor
        character(len=256) :: grp_to_hist 
        character(len=256) :: grp_so_hist 
        character(len=256) :: grp_tf_hist 
        character(len=256) :: grp_to_proj 
        character(len=256) :: grp_so_proj 
        character(len=256) :: grp_tf_proj

        character(len=256) :: grp_mask_shlf_proj
        
        integer  :: iloc, k 
        real(wp) :: tmp
        real(wp) :: time_par_proj(4) 
        real(wp) :: time_par_proj_msk(4) 

        ! First determine whether gcm+scenario provided or experiment
        ! obtain valid values for gcm and scenario to start.
        if (present(experiment)) then 
            ! Experiment provided, obtain gcm+scenario from it 

            if (trim(experiment) .eq. "ctrl" .or. trim(experiment) .eq. "ctrl0") then 
                gcm_now      = trim(experiment)
                scenario_now = trim(experiment)
            else 
                iloc = index(experiment,"_")
                if (iloc == 0) then
                    write(error_unit,*) "nahosmip_ant_forcing_init:: Error: experiment &
                    &argument must be defined as gcm_scenario."
                    write(error_unit,*) "experiment = ", trim(experiment)
                    stop
                end if

                gcm_now = experiment(1:iloc-1)
                scenario_now = experiment(iloc+1:len_trim(experiment))

            end if 

        else if (present(gcm) .and. present(scenario)) then 
            ! gcm and scenario provided go forward as usual 

            gcm_now      = trim(gcm)
            scenario_now = trim(scenario)

        else
            ! Arguments not provided
            write(error_unit,*) "nahosmip_ant_forcing_init:: Error: gcm+scenario or experiment must be provided &
            &as arguments."
            stop
        end if

        ! Special case for control runs, use "NorESM1-M_RCP26-repeat"
        if (trim(scenario_now) .eq. "ctrl" .or. trim(scenario_now) .eq. "ctrl0") then
            ! Assign specific gcm+scenario for control runs

            ism%gcm        =  "NorESM1-M"
            ism%scenario   =  "RCP26-repeat"
            ism%experiment = "NorESM1-M_RCP26-repeat"
            
            ! Save control run name ("ctrl0" or "ctrl") here
            ism%ctrl_run_type = trim(scenario_now)

        else
            ! Define the current experiment characteristics

            ism%gcm        = trim(gcm_now)
            ism%scenario   = trim(scenario_now) 
            ism%experiment = trim(ism%gcm)//"_"//trim(ism%scenario) 

            ! Not a control simulation
            ism%ctrl_run_type = "none"

        end if

        write(*,*)
        write(*,*) "nahosmip_ant_forcing_init:: summary"
        write(*,*) "ctrl_run_type: ", trim(ism%ctrl_run_type)
        write(*,*) "gcm:           ", trim(ism%gcm)
        write(*,*) "scenario:      ", trim(ism%scenario)
        write(*,*) "experiment:    ", trim(ism%experiment)
        write(*,*) 

        select case(trim(ism%experiment))

            case("CCSM4_RCP85",                 &
                 "CESM2-WACCM_ssp585",          &
                 "CESM2-WACCM_ssp585-repeat",   &
                 "HadGEM2-ES_RCP85",            &
                 "HadGEM2-ES_RCP85-repeat",     &
                 "NorESM1-M_RCP26-repeat",      &
                 "NorESM1-M_RCP85-repeat",      &
                 "UKESM1-0-LL_ssp126",          &
                 "UKESM1-0-LL_ssp585",          &
                 "UKESM1-0-LL_ssp585-repeat")
                ! Control and scenarios use the same files 
                ! since ctrl specific forcing adapted in update step 

                group_prefix = "gcm_"

                grp_ts_ref   = trim(group_prefix)//"ts_ref"
                grp_pr_ref   = trim(group_prefix)//"pr_ref"
                grp_smb_ref  = trim(group_prefix)//"smb_ref"
                grp_ts_hist  = trim(group_prefix)//"ts_hist"
                grp_pr_hist  = trim(group_prefix)//"pr_hist"
                grp_smb_hist = trim(group_prefix)//"smb_hist"
                grp_ts_proj  = trim(group_prefix)//"ts_proj"
                grp_pr_proj  = trim(group_prefix)//"pr_proj"
                grp_smb_proj = trim(group_prefix)//"smb_proj"
                
                grp_to_ref   = "to_ref"
                grp_so_ref   = "so_ref"
                grp_tf_ref   = "tf_ref"
                grp_tf_cor   = "tf_cor"
                grp_to_hist  = trim(group_prefix)//"to_hist"
                grp_so_hist  = trim(group_prefix)//"so_hist"
                grp_tf_hist  = trim(group_prefix)//"tf_hist"
                grp_to_proj  = trim(group_prefix)//"to_proj"
                grp_so_proj  = trim(group_prefix)//"so_proj"
                grp_tf_proj  = trim(group_prefix)//"tf_proj"

                grp_mask_shlf_proj = trim(group_prefix)//"mask_shlf_proj"
                
            case DEFAULT 

                write(*,*) "nahosmip_ant_forcing_init:: Error: exeriment (== gcm_scenario) not recognized."
                write(*,*) "experiment = ", trim(ism%experiment) 
                write(*,*) "gcm        = ", trim(ism%gcm) 
                write(*,*) "scenario   = ", trim(ism%scenario) 

                stop 

        end select

        ! Adjust time_par_proj manually here, since there are inconsistencies between the projection runs
        ! (some runs end in 2299, 2300, or 2301). Hard code choices here to avoid many changes
        ! in parameter file
        select case(trim(ism%experiment))

            case("CESM2-WACCM_ssp585","HadGEM2-ES_RCP85")
                ! Cases that end on year 2299

                time_par_proj     = [1995.0,2299.0,1.0,0.0]
                time_par_proj_msk = [1995.0,2300.0,1.0,0.0]

            case("UKESM1-0-LL_ssp585")
                ! Cases that end on year 2301
            
                time_par_proj     = [-1.0,-1.0,-1.0,0.0]
                time_par_proj_msk = [1995.0,2301.0,1.0,0.0]

            case DEFAULT
                ! Set negative values to time_par so that values are used directly from the file

                time_par_proj     = [-1.0,-1.0,-1.0,0.0]
                time_par_proj_msk = [-1.0,-1.0,-1.0,0.0]

        end select

        ! Note: reference fields are hard-coded in the namelist file to come from the default
        ! reference simulation, which is NorESM1-M_rcp26-repeat

        ! Note: for the reference and historical fields, always use the time_par from the namelist file.
        
        ! Initialize all variables from namelist entries 

        ! General fields 
        call varslice_init_nml_ismip6(ism%basins,  filename,"imbie_basins",domain,grid_name,ism%gcm,ism%scenario)
        
        ! Amospheric fields
        call varslice_init_nml_ismip6(ism%ts_ref,  filename,trim(grp_ts_ref), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%pr_ref,  filename,trim(grp_pr_ref), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_ref, filename,trim(grp_smb_ref),domain,grid_name,ism%gcm,ism%scenario)
        
        call varslice_init_nml_ismip6(ism%ts_hist, filename,trim(grp_ts_hist), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%pr_hist, filename,trim(grp_pr_hist), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_hist,filename,trim(grp_smb_hist),domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%ts_proj, filename,trim(grp_ts_proj), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%pr_proj, filename,trim(grp_pr_proj), domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%smb_proj,filename,trim(grp_smb_proj),domain,grid_name,ism%gcm,ism%scenario)

        ! Oceanic fields
        call varslice_init_nml_ismip6(ism%to_ref,  filename,trim(grp_to_ref),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%so_ref,  filename,trim(grp_so_ref),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_ref,  filename,trim(grp_tf_ref),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_cor,  filename,trim(grp_tf_cor),domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%to_hist, filename,trim(grp_to_hist),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%so_hist, filename,trim(grp_so_hist),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_hist, filename,trim(grp_tf_hist),domain,grid_name,ism%gcm,ism%scenario)

        call varslice_init_nml_ismip6(ism%to_proj, filename,trim(grp_to_proj),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%so_proj, filename,trim(grp_so_proj),domain,grid_name,ism%gcm,ism%scenario)
        call varslice_init_nml_ismip6(ism%tf_proj, filename,trim(grp_tf_proj),domain,grid_name,ism%gcm,ism%scenario)

        ! Load time-independent fields

        ! Amospheric fields 
        call varslice_update(ism%ts_ref)
        call varslice_update(ism%pr_ref)
        call varslice_update(ism%smb_ref)

        ! Oceanic fields
        call varslice_update(ism%to_ref)
        call varslice_update(ism%so_ref)
        call varslice_update(ism%tf_ref)
        call varslice_update(ism%tf_cor)

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%to_ref%var,3)
            if (count(ism%to_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(ism%to_ref%var(:,:,k,1),mask=ism%to_ref%var(:,:,k,1) .ne. mv)
                where(ism%to_ref%var(:,:,k,1) .eq. mv) 
                    ism%to_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%so_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%so_ref%var(:,:,k,1),mask=ism%so_ref%var(:,:,k,1) .ne. mv)
                where(ism%so_ref%var(:,:,k,1) .eq. mv) 
                    ism%so_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%tf_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf_ref%var(:,:,k,1),mask=ism%tf_ref%var(:,:,k,1) .ne. mv)
                where(ism%tf_ref%var(:,:,k,1) .eq. mv) 
                    ism%tf_ref%var(:,:,k,1) = tmp
                end where 
            end if
        end do

        ! Initialize iceberg_mask variable in case it is needed
        allocate(ism%iceberg_mask(size(ism%to_ref%var,1),size(ism%to_ref%var,2)))
        ism%iceberg_mask = 0.0

        return 

    end subroutine nahosmip_ant_forcing_init

    subroutine nahosmip_ant_forcing_update(ism,time,use_ref_atm,use_ref_ocn)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        real(wp), intent(IN) :: time
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 

        ! Local variables 
        integer  :: k 
        real(wp) :: tmp 
        character(len=56) :: slice_method 

        ! Get slices for current time

        slice_method = "extrap" 

        ! === Atmospheric fields ==================================
        
        if (trim(ism%ctrl_run_type) .eq. "ctrl0") then 
            ! For control scenario "ctrl0" with no variability at all, 
            ! override time choices and set control atm and ocn 

            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
        
        else if (time .lt. 1950) then
            ! No data available, impose climatological mean 

            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
        
        else if (time .ge. 1950 .and. time .lt. 2015) then
            ! Get data from historical forcing datasets

            call varslice_update(ism%ts_hist, [time],method=slice_method)
            call varslice_update(ism%pr_hist, [time],method=slice_method)
            call varslice_update(ism%smb_hist,[time],method=slice_method)

            ism%ts  = ism%ts_hist
            ism%pr  = ism%pr_hist
            ism%smb = ism%smb_hist

        else if (time .ge. 2015 .and. time .le. 2300) then
            ! Get data from projection forcing datasets

            call varslice_update(ism%ts_proj, [time],method=slice_method)
            call varslice_update(ism%pr_proj, [time],method=slice_method)
            call varslice_update(ism%smb_proj,[time],method=slice_method) 

            ism%ts  = ism%ts_proj
            ism%pr  = ism%pr_proj
            ism%smb = ism%smb_proj

        else ! time .gt. 2300
            ! Take the average of the last 10 years of projection forcing data

            call varslice_update(ism%ts_proj, [2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%pr_proj, [2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%smb_proj,[2290.0_wp,2300.0_wp],method="range_mean")

            ism%ts  = ism%ts_proj
            ism%pr  = ism%pr_proj
            ism%smb = ism%smb_proj

        end if

        ! === Oceanic fields ==================================
        
        if (trim(ism%ctrl_run_type) .eq. "ctrl0") then 
            ! For control scenario "ctrl0" with no variability at all, 
            ! override time choices and set control atm and ocn 

            ! Set ocn fields to reference values 
            ism%to = ism%to_ref 
            ism%so = ism%so_ref 
            ism%tf = ism%tf_ref 

        else if (time .lt. 1950) then 
            ! No data available, impose climatological mean 

            ism%to = ism%to_ref
            ism%so = ism%so_ref
            ism%tf = ism%tf_ref

        else if (time .ge. 1950 .and. time .lt. 2015) then
            ! Historic 

            ism%to = ism%to_ref
            ism%so = ism%so_ref
            ism%tf = ism%tf_ref

            ! Oceanic fields 
            call varslice_update(ism%to_hist,[time],method=slice_method)
            call varslice_update(ism%so_hist,[time],method=slice_method)
            call varslice_update(ism%tf_hist,[time],method=slice_method)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj

        else if (time .ge. 2015 .and. time .le. 2300) then
            ! Projection period 1 

            call varslice_update(ism%to_proj,[time],method=slice_method)
            call varslice_update(ism%so_proj,[time],method=slice_method)
            call varslice_update(ism%tf_proj,[time],method=slice_method)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj

        else ! time .gt. 2300
            ! Projection period 2 

            call varslice_update(ism%to_proj,[2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%so_proj,[2290.0_wp,2300.0_wp],method="range_mean")
            call varslice_update(ism%tf_proj,[2290.0_wp,2300.0_wp],method="range_mean")

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj

        end if

        ! === Mask shelf collapse ==========================

        if (ism%shlf_collapse) then 
            if (time .lt. 2015) then
                
                ! Note: the year 2000 is used for this mask, however, any value from 2000 to 2015 
                ! could be used as nothing really changes with the mask in the historical period.
                ! Retreat only begins to be seen some decades later.
                call varslice_update(ism%mask_shlf_proj,[2000.0_wp],method=slice_method)
                ism%mask_shlf = ism%mask_shlf_proj

            else if (time .ge. 2015 .and. time .le. 2300) then

                call varslice_update(ism%mask_shlf_proj,[time],method=slice_method)
                ism%mask_shlf = ism%mask_shlf_proj
            
            else

                call varslice_update(ism%mask_shlf_proj,[2300.0_wp],method=slice_method)
                ism%mask_shlf  = ism%mask_shlf_proj

            end if
        end if 

        ! === Additional calculations ======================

        if (trim(ism%ctrl_run_type) .eq. "ctrl") then 
            ! For control scenario "ctrl", override above choices and 
            ! set control atm and ocn for projection time periods.
            ! This means historical variability is allowed, but 
            ! then the climatological average is imposed for the future.
            ! See protocol at: https://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections2300-Antarctica#Initialization.2C_historical_run.2C_control_run.2C_and_projection_runs

            if (time .ge. 2015) then

                ! Set atm fields to reference values (ie, zero anomaly) 
                ism%ts  = ism%ts_ref 
                ism%pr  = ism%pr_ref 
                ism%smb = ism%smb_ref 

                ! Since atm fields are anomalies, also set actual variable to zero 
                ism%ts%var  = 0.0_wp 
                ism%pr%var  = 0.0_wp 
                ism%smb%var = 0.0_wp 
            
                ! Set ocn fields to reference values 
                ism%to = ism%to_ref 
                ism%so = ism%so_ref 
                ism%tf = ism%tf_ref 

            end if

        end if 
        
        ! If desired, override other choices and use reference atm fields
        if (present(use_ref_atm)) then 
        if (use_ref_atm) then  

            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
        end if 
        end if

        ! If desired, override other choices and use reference ocean fields
        if (present(use_ref_ocn)) then 
        if (use_ref_ocn) then  

            ism%to = ism%to_ref 
            ism%so = ism%so_ref 
            ism%tf = ism%tf_ref 

        end if 
        end if

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%to%var,3)
            if (count(ism%to%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(ism%to%var(:,:,k,1),mask=ism%to%var(:,:,k,1) .ne. mv)
                where(ism%to%var(:,:,k,1) .eq. mv) 
                    ism%to%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%so%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%so%var(:,:,k,1),mask=ism%so%var(:,:,k,1) .ne. mv)
                where(ism%so%var(:,:,k,1) .eq. mv) 
                    ism%so%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%tf%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf%var(:,:,k,1),mask=ism%tf%var(:,:,k,1) .ne. mv)
                where(ism%tf%var(:,:,k,1) .eq. mv) 
                    ism%tf%var(:,:,k,1) = tmp
                end where 
            end if
        end do
        
        return 

    end subroutine nahosmip_ant_forcing_update

end module ismip6
