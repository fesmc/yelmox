module esm
    ! This module contains routines that help with performing the esm suite
    ! of experiments. 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml  
    use ncio 
    use varslice
    use marine_shelf

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Define default missing value 
    real(wp), parameter :: mv = -9999.0_wp 

    ! Need to call marshelf class for oceanic interpolation into shelves
    type(marshelf_class)        :: mshlf

    ! Class for holding ice-forcing data from esm archives
    type esm_forcing_class
        
        ! Experiment information
        character(len=256)     :: gcm 
        character(len=256)     :: scenario
        character(len=256)     :: experiment 
        character(len=256)     :: domain 
        character(len=256)     :: grid_name 
        character(len=256)     :: ctrl_run_type
        real(wp)               :: lapse(2)
        real(wp)               :: beta_p
        real(wp)               :: f_ocn
        real(wp)               :: f_polar
        real(wp)               :: dT_lim
        character(len=256)     :: grid_src

        ! === Reference climatology ===
        ! Atmosphere
        type(varslice_class)   :: ts_ref
        type(varslice_class)   :: pr_ref
        
        ! Ocean
        type(varslice_class)   :: to_ref, to_ref_src
        type(varslice_class)   :: so_ref

        ! ===  Variability climatologies ===
        ! Atmosphere
        type(varslice_class)   :: ts_var
        type(varslice_class)   :: pr_var
        type(varslice_class)   :: ts_var_ref
        type(varslice_class)   :: pr_var_ref

        ! Ocean
        type(varslice_class)   :: to_var
        type(varslice_class)   :: so_var
        type(varslice_class)   :: to_var_ref
        type(varslice_class)   :: so_var_ref
        
        ! ===      ESM      ===
        ! Atmospheric fields
        type(varslice_class)   :: ts_esm_ref 
        type(varslice_class)   :: pr_esm_ref 

        type(varslice_class)   :: ts_hist 
        type(varslice_class)   :: pr_hist 

        type(varslice_class)   :: ts_proj
        type(varslice_class)   :: pr_proj

        ! Oceanic fields 
        type(varslice_class)   :: to_esm_ref
        type(varslice_class)   :: so_esm_ref

        type(varslice_class)   :: to_hist
        type(varslice_class)   :: so_hist

        type(varslice_class)   :: to_proj
        type(varslice_class)   :: so_proj

        ! General fields 
        type(varslice_class)   :: basins
        type(varslice_class)   :: zs_ref
        type(varslice_class)   :: zs_var    
        type(varslice_class)   :: zs_esm_ref
        type(varslice_class)   :: zs_hist
        type(varslice_class)   :: zs_proj
        
        ! === Diagnostic fields ===

        ! === Atmosphere ===
        ! Monthly fields (but right now constant anomaly between months)
        real(wp), allocatable :: t2m(:,:,:)       ! Monthly surface temperature [K]
        real(wp), allocatable :: pr(:,:,:)        ! Monthly precipitation [mm/yr]
        ! Anomalies
        real(wp), allocatable :: dts(:,:,:)       ! Surface temperature anomaly [K]
        real(wp), allocatable :: dpr(:,:,:)       ! Precipitation relative anomaly [%]
        real(wp), allocatable :: dts_var(:,:,:)   ! Surface temperature anomaly variability[K]
        real(wp), allocatable :: dpr_var(:,:,:)   ! Precipitation relative anomaly variability [%]

        ! ===    Ocean   ===
        real(wp), allocatable :: dto(:,:)         ! Surface temperature anomaly [K]
        real(wp), allocatable :: dso(:,:)         ! Precipitation relative anomaly [%]
        real(wp), allocatable :: dto_var(:,:)     ! Surface temperature anomaly variability [K]
        real(wp), allocatable :: dso_var(:,:)     ! Precipitation relative anomaly variability [%]
        
        ! === Mean fields ===
        real(wp), allocatable :: t2m_sum(:,:)     ! Summer surface temperature [K]
        real(wp), allocatable :: t2m_ann(:,:)     ! Annual surface temperature [K]
        real(wp), allocatable :: pr_ann(:,:)      ! Annual precipitation [mm/yr]

    end type

    type esm_ice_var_class
        character(len=56)  :: name 
        character(len=128) :: long_name
        character(len=12)  :: var_type
        character(len=128) :: standard_name 
        character(len=128) :: units_in
        character(len=128) :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
    end type

    type esm_experiment_class
        character(len=56)   :: expname
        character(len=56)   :: group
        character(len=56)   :: model
        character(len=256)  :: experiment
        character(len=256)  :: file_suffix
    end type
        
    ! Class for holding ice output for writing to standard formats...
    type esm_ice_class
        type(esm_ice_var_class), allocatable :: vars(:)
    end type 

    private
    public :: esm_forcing_class
    public :: esm_experiment_class
    public :: esm_ice_class

    ! General routines
    public :: esm_experiment_def
    public :: esm_forcing_init
    public :: esm_forcing_update
    public :: esm_variability_update
    public :: esm_clim_update
    
    public :: esm_write_init
    public :: esm_write_step
    
contains
    
    subroutine esm_experiment_def(ie,expname,filename,group,model,domain)

        implicit none

        type(esm_experiment_class), intent(OUT) :: ie
        character(len=*), intent(IN) :: expname
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group 
        character(len=*), intent(IN) :: model
        character(len=*), intent(IN) :: domain 

        ! Save the experiment name (ctrlAE, expAE01, etc)
        ie%expname = trim(expname)
        ie%group   = trim(group)
        ie%model   = trim(model)
        ie%model   = trim(domain)

        ! Load parameters associated with this experiment
        call nml_read(filename,ie%expname,"experiment",   ie%experiment)
        
        ! Define the ouput filename according to protocol
        ie%file_suffix = trim(domain)//"_"//trim(group)//"_"//trim(model)//"_"//trim(ie%expname)//".nc"

        return
        
    end subroutine esm_experiment_def

    subroutine esm_forcing_init(esm,filename,domain,grid_name,gcm,scenario,experiment,use_esm,use_hist,use_proj)

        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN), optional :: gcm
        character(len=*), intent(IN), optional :: scenario
        character(len=*), intent(IN), optional :: experiment
        logical,          intent(IN), optional :: use_esm, use_hist, use_proj
    
        ! Local variables 
        character(len=256) :: gcm_now
        character(len=256) :: experiment_now
        character(len=256) :: scenario_now
            
        character(len=256) :: group_prefix 
    
        ! Reference climatology
        character(len=256) :: grp_ts_ref  
        character(len=256) :: grp_pr_ref
        character(len=256) :: grp_zs_ref
        character(len=256) :: grp_to_ref 
        character(len=256) :: grp_so_ref
        ! Variability period
        character(len=256) :: grp_ts_var  
        character(len=256) :: grp_pr_var
        character(len=256) :: grp_zs_var
        character(len=256) :: grp_to_var 
        character(len=256) :: grp_so_var
        ! ESM Reference period
        character(len=256) :: grp_ts_esm_ref  
        character(len=256) :: grp_pr_esm_ref
        character(len=256) :: grp_zs_esm_ref
        character(len=256) :: grp_to_esm_ref 
        character(len=256) :: grp_so_esm_ref
        ! ESM Historical period 
        character(len=256) :: grp_ts_hist 
        character(len=256) :: grp_pr_hist
        character(len=256) :: grp_zs_hist
        character(len=256) :: grp_to_hist 
        character(len=256) :: grp_so_hist
        ! ESM Projection period 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
        character(len=256) :: grp_zs_proj
        character(len=256) :: grp_to_proj 
        character(len=256) :: grp_so_proj

        integer  :: iloc, k 
        real(wp) :: tmp
        real(wp) :: time_par_ref(4),time_par_hist(4),time_par_proj(4),time_par_var(4)
    
        ! First determine whether gcm+scenario provided or experiment
        ! obtain valid values for gcm and scenario to start.
        if (present(experiment)) then 
            ! Experiment provided, obtain gcm+scenario from it 
    
            if (trim(experiment) .eq. "ctrl" .or. trim(experiment) .eq. "ctrl0") then 
                gcm_now        = trim(experiment)
                scenario_now   = trim(experiment)
                experiment_now = "ctrl"
            ! jablasco: change this?
            elseif (trim(experiment) .eq. "1pctCO2") then
                gcm_now      = trim(experiment)
                scenario_now = "transient"
                experiment_now = "transient" 
            elseif (trim(experiment) .eq. "ssp126" .or. trim(experiment) .eq. "ssp585") then
                gcm_now      = trim(experiment)
                scenario_now = "transient"
            else 
                iloc = index(experiment,"_")
                if (iloc == 0) then
                    write(error_unit,*) "esm_ant_forcing_init:: Error: experiment &
                    &argument must be defined as gcm_scenario."
                    write(error_unit,*) "experiment = ", trim(experiment)
                    stop
                end if
    
                gcm_now = experiment(1:iloc-1)
                scenario_now = experiment(iloc+1:len_trim(experiment))
    
            end if 
    
        else if (present(gcm) .and. present(scenario)) then 
            ! gcm and scenario provided go forward as usual 
    
            gcm_now        = trim(gcm)
            scenario_now   = trim(scenario)
            experiment_now = "transient"
    
        else
            ! Arguments not provided 
            write(error_unit,*) "esm_forcing_init:: Error: gcm+scenario or experiment must be provided &
            &as arguments."
            stop
        end if
    
        ! Define the current experiment characteristics
        ! jablasco: que hace esto?
        esm%ctrl_run_type = trim(experiment_now)
        esm%gcm           = trim(gcm_now)
        esm%scenario      = trim(scenario_now) 
        esm%experiment    = trim(esm%gcm)//"_"//trim(esm%scenario) 
    
        write(*,*)
        write(*,*) "esm_forcing_init:: summary"
        write(*,*) "ctrl_run_type: ", trim(esm%ctrl_run_type)
        write(*,*) "gcm:           ", trim(esm%gcm)
        write(*,*) "scenario:      ", trim(esm%scenario)
        write(*,*) "experiment:    ", trim(esm%experiment)
        write(*,*) 
    
        group_prefix = "gcm_"
    
        ! Reference climatology
        grp_ts_ref       = trim(group_prefix)//"ts_ref"
        grp_pr_ref       = trim(group_prefix)//"pr_ref"
        grp_zs_ref       = trim(group_prefix)//"zs_ref"
        grp_to_ref       = trim(group_prefix)//"to_ref"
        grp_so_ref       = trim(group_prefix)//"so_ref"   

        ! Variability climatology
        grp_ts_var       = trim(group_prefix)//"ts_var"
        grp_pr_var       = trim(group_prefix)//"pr_var"
        grp_zs_var       = trim(group_prefix)//"zs_var"
        grp_to_var       = trim(group_prefix)//"to_var"
        grp_so_var       = trim(group_prefix)//"so_var"

        ! ESM Reference climatology (to compute anomalies)
        grp_ts_esm_ref   = trim(group_prefix)//"ts_esm_ref"
        grp_pr_esm_ref   = trim(group_prefix)//"pr_esm_ref"
        grp_zs_esm_ref   = trim(group_prefix)//"zs_esm_ref"
        grp_to_esm_ref   = trim(group_prefix)//"to_esm_ref"
        grp_so_esm_ref   = trim(group_prefix)//"so_esm_ref" 

        ! ESM Historical sims 
        grp_ts_hist  = trim(group_prefix)//"ts_hist"
        grp_pr_hist  = trim(group_prefix)//"pr_hist"
        grp_zs_hist  = trim(group_prefix)//"zs_hist"
        grp_to_hist  = trim(group_prefix)//"to_hist"
        grp_so_hist  = trim(group_prefix)//"so_hist"

        ! ESM projected sims
        grp_ts_proj  = trim(group_prefix)//"ts_proj"
        grp_pr_proj  = trim(group_prefix)//"pr_proj"
        grp_zs_proj  = trim(group_prefix)//"zs_proj"
        grp_to_proj  = trim(group_prefix)//"to_proj"
        grp_so_proj  = trim(group_prefix)//"so_proj"         
 
        ! Initialize all variables from namelist entries 
        ! General fields (needed? switch to reese basins?)
        call varslice_init_nml_esm(esm%basins, filename, "basins", domain, grid_name, esm%gcm, esm%scenario)
            
        ! Climatology
        ! Reference period
        call varslice_init_nml_esm(esm%ts_ref, filename, trim(grp_ts_ref), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%pr_ref, filename, trim(grp_pr_ref), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%zs_ref, filename, trim(grp_zs_ref), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%to_ref, filename, trim(grp_to_ref), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%so_ref, filename, trim(grp_so_ref), domain, grid_name, esm%gcm, esm%scenario)

        ! Initialize variables at other grids if source grid is different
        call varslice_init_nml_esm(esm%to_ref_src, filename, trim(grp_to_ref), domain, grid_name, esm%gcm, esm%scenario)

        ! Variability
        ! Transient dependent field
        call varslice_init_nml_esm(esm%ts_var, filename, trim(grp_ts_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%pr_var, filename, trim(grp_pr_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%to_var, filename, trim(grp_to_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%so_var, filename, trim(grp_so_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%zs_var, filename, trim(grp_zs_var), domain, grid_name, esm%gcm, esm%scenario)
        ! Reference period
        call varslice_init_nml_esm(esm%ts_var_ref, filename, trim(grp_ts_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%pr_var_ref, filename, trim(grp_pr_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%to_var_ref, filename, trim(grp_to_var), domain, grid_name, esm%gcm, esm%scenario)
        call varslice_init_nml_esm(esm%so_var_ref, filename, trim(grp_so_var), domain, grid_name, esm%gcm, esm%scenario)

        ! Transient dependent fields
        if (trim(esm%ctrl_run_type) .eq. "transient") then
            ! ESM reference period
            !if (esm%use_esm) then
            if (use_esm .and. .TRUE.) then
                call varslice_init_nml_esm(esm%ts_esm_ref, filename, trim(grp_ts_esm_ref), domain, grid_name, esm%gcm, esm%scenario)
                call varslice_init_nml_esm(esm%pr_esm_ref, filename, trim(grp_pr_esm_ref), domain, grid_name, esm%gcm, esm%scenario)
                call varslice_init_nml_esm(esm%to_esm_ref, filename, trim(grp_to_esm_ref), domain, grid_name, esm%gcm, esm%scenario)
                call varslice_init_nml_esm(esm%so_esm_ref, filename, trim(grp_so_esm_ref), domain, grid_name, esm%gcm, esm%scenario)
                !call varslice_init_nml_esm(esm%zs_esm_ref, filename, trim(grp_zs_esm_ref), domain, grid_name, esm%gcm, esm%scenario)
            end if
            ! ESM historical period
            !if(trim(esm%use_hist)) then
            if (.FALSE.) then
                call varslice_init_nml_esm(esm%ts_hist, filename,trim(grp_ts_hist), domain, grid_name, esm%gcm, esm%scenario)
                call varslice_init_nml_esm(esm%pr_hist, filename,trim(grp_pr_hist), domain, grid_name, esm%gcm, esm%scenario)
                call varslice_init_nml_esm(esm%to_hist, filename,trim(grp_to_hist), domain, grid_name, esm%gcm, esm%scenario)
                call varslice_init_nml_esm(esm%so_hist, filename,trim(grp_so_hist), domain, grid_name, esm%gcm, esm%scenario)
                !call varslice_init_nml_esm(esm%zs_hist, filename,trim(grp_zs_hist), domain, grid_name, esm%gcm, esm%scenario)
            end if
            ! ESM projection period
            !if(trim(esm%use_proj)) then
            ! jablasco
            if (use_esm .and. .TRUE.) then
            !if (.TRUE.) then    
                call varslice_init_nml_esm(esm%ts_proj, filename,trim(grp_ts_proj), domain,grid_name,esm%gcm,esm%scenario)
                call varslice_init_nml_esm(esm%pr_proj, filename,trim(grp_pr_proj), domain,grid_name,esm%gcm,esm%scenario)
                call varslice_init_nml_esm(esm%to_proj, filename,trim(grp_to_proj), domain,grid_name,esm%gcm,esm%scenario)
                call varslice_init_nml_esm(esm%so_proj, filename,trim(grp_so_proj), domain,grid_name,esm%gcm,esm%scenario)
                !call varslice_init_nml_esm(esm%zs_proj, filename,trim(grp_zs_proj), domain, grid_name, esm%gcm, esm%scenario)
                end if
        end if

        ! Allocate objects (use a time independent field)
        call varslice_update(esm%zs_ref)
        call esm_allocate(esm,size(esm%zs_ref%var,1),size(esm%zs_ref%var,2))
        write(*,*) "dim1, dim2", size(esm%zs_ref%var,1),size(esm%zs_ref%var,2)

        return 
    
    end subroutine esm_forcing_init

    subroutine esm_clim_update(esm,z_srf_ylm,time,time_ref,domain,grid_name)
        ! Routine to update reference climatology to the specific Antarctic elevation and ocean (neccessary?)

        implicit none

        type(esm_forcing_class), intent(INOUT) :: esm
        real(wp),                intent(IN)    :: z_srf_ylm(:,:)
        real(wp),                intent(IN)    :: time
        real(wp),                intent(IN)    :: time_ref(2)
        character(len=*),        intent(IN)    :: domain, grid_name

        ! Local variables 
        integer :: m
        real(wp) :: tmp, lapse
        real(wp), parameter :: pi = 3.14159265359 
        character(len=56)   :: slice_method, ref_grid_name 
        !type(map_scrip_class) :: mps
        logical  :: south 

        ! Get slices for current time
        slice_method = "extrap"

        ! select domain
        south = .FALSE. 
        if (trim(domain).eq."Antarctica") south = .TRUE.

        ! Climatology reference
        ! === Atmospheric fields ===
        call varslice_update(esm%ts_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        call varslice_update(esm%pr_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        ! ===   Oceanic fields   ===
        !call nc_read_attr(esm%to_ref,"grid_name", ref_grid_name)
        !if (trim(file_grid_name) .eq. trim(grid_name) ) then
        !    ! Ref grid and Yelmo grid are the same
        !    call varslice_update(esm%to_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)
        !else
        !    ! Ref's grid is different than Yelmo grid. Load desired time range for the source code.
        !    call varslice_update(esm%to_ref_src, [time_ref(1),time_ref(2)],method="range_mean",rep=1)
        !    ! Load the scrip map from file (should already have been generated via cdo externally)
        !    call map_scrip_init(mps,file_grid_name,grid_name,method="con",fldr="maps",load=.TRUE.)
        !    ! Remap src into the desired target
        !    call varslice_map_to_grid(esm%to_ref,esm%to_ref_src,mps)
        !end if

        call varslice_update(esm%to_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)
        call varslice_update(esm%so_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)  

        ! Variability reference
        ! === Atmospheric fields ===
        call varslice_update(esm%ts_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        call varslice_update(esm%pr_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        ! ===   Oceanic fields   ===
        call varslice_update(esm%to_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)
        call varslice_update(esm%so_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1) 

        ! Convert atmospheric fields to model elevation
        do m = 1,12 
            if(south) then ! Southern Hemisphere
                lapse = (esm%lapse(1)+(esm%lapse(2)-esm%lapse(1))*cos(2*pi*(m*30.4375-30.4375)/365.25))
            else ! Northern Hemisphere
                lapse = (esm%lapse(1)+(esm%lapse(1)-esm%lapse(2))*cos(2*pi*(m*30.4375-30.4375)/365.25))
            end if    
            esm%t2m(:,:,m) = esm%ts_ref%var(:,:,m,1) + lapse*(esm%zs_ref%var(:,:,1,1)-z_srf_ylm)
            esm%pr(:,:,m)  = esm%pr_ref%var(:,:,m,1) * exp(esm%beta_p*lapse*(esm%zs_ref%var(:,:,1,1)-z_srf_ylm))
        end do

        ! Compute diagnostic fields
        esm%t2m_ann = sum(esm%t2m,dim=3) / 12.0
        esm%pr_ann  = sum(esm%pr,dim=3) / 12.0 
        if(south) then
            esm%t2m_sum = (esm%t2m(:,:,1)+esm%t2m(:,:,2)+esm%t2m(:,:,12)) / 3.0
        else
            esm%t2m_sum = (esm%t2m(:,:,6)+esm%t2m(:,:,7)+esm%t2m(:,:,8)) / 3.0
        end if      

        return

    end subroutine esm_clim_update

    subroutine esm_variability_update(esm,mshlf,time,dtt,clim_var,time_ref,H_ice,basins,z_bed,f_grnd,z_sl,use_ref_atm,use_ref_ocn)
        ! Update climatic fields. These will be used as bnd conditions for Yelmo.
        ! Output are anomaly fields with respect to a reference field from the ESM.
    
        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        type(marshelf_class),    intent(IN)    :: mshlf
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: dtt
        character(len=*), intent(IN) :: clim_var
        real(wp), intent(IN) :: time_ref(2)
        real(wp), intent(IN) :: H_ice(:,:),basins(:,:),z_bed(:,:),f_grnd(:,:),z_sl(:,:)
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 
    
        ! Local variables 
        integer  :: m, year, step_idx
        real(wp) :: tmp, lapse
        real(wp) :: rand, year_rand

        ! Initialize anomalies
        esm%dts_var = 0.0_wp
        esm%dpr_var = 1.0_wp
        esm%dto_var = 0.0_wp
        esm%dso_var = 0.0_wp
            
        ! Obtain reference year climatology
        select case(trim(clim_var))
            case("random","rand","white_noise")
                ! Select a random year from the climatology period
                call random_number(rand)
                year_rand = NINT((time_ref(2)-time_ref(1))*rand + time_ref(1))  
                write(*,*) "year_rand = ", year_rand
                    
            case("historic","hist")
                ! Cycle through the selected climatology period
                ! Compute the number of timesteps of size dtt that have passed
                step_idx = INT(time / dtt)      
                year_rand = INT(MOD(step_idx, INT(time_ref(2)-time_ref(1)+1))) + INT(time_ref(1))
                write(*,*) "year_cyclic = ", year_rand  

            case("red_noise")
                ! Select a random year but based on red noise
                ! TO DO
                write(*,*) "Red noise method not available yet."
                STOP

            case DEFAULT
                ! Do nothing for now
        end select

        select case(trim(clim_var))
            case("random","rand","white_noise","historic","hist","red_noise")
                ! === Atmospheric fields === 
                call varslice_update(esm%ts_var,[year_rand],method="interp",rep=12)
                call varslice_update(esm%pr_var,[year_rand],method="interp",rep=12)
                do m = 1, 12 
                    esm%dts_var(:,:,m) = esm%ts_var%var(:,:,m,1)-esm%ts_var_ref%var(:,:,m,1)
                    esm%dpr_var(:,:,m) = esm%pr_var%var(:,:,m,1)/(esm%pr_var_ref%var(:,:,m,1)+1e-8)
                end do
                ! ===   Oceanic fields   ===
                call varslice_update(esm%to_var,[year_rand],method="interp",rep=1)
                call varslice_update(esm%so_var,[year_rand],method="interp",rep=1)

                ! Interpolate ocean data to the interior of ice shelves
                call ocn_variable_extrapolation(esm%to_var%var(:,:,:,1),H_ice,basins,-esm%to_var%z,z_bed)
                call ocn_variable_extrapolation(esm%so_var%var(:,:,:,1),H_ice,basins,-esm%so_var%z,z_bed)

                ! Compute the anomaly at the desired depth level
                call marshelf_interp_shelf(esm%dto_var,mshlf,esm%to_var%var(:,:,:,1)-esm%to_var_ref%var(:,:,:,1),H_ice, &
                                            z_bed,f_grnd,z_sl,-esm%to_var_ref%z)
                call marshelf_interp_shelf(esm%dso_var,mshlf,esm%so_var%var(:,:,:,1)-esm%so_var_ref%var(:,:,:,1),H_ice, &
                                            z_bed,f_grnd,z_sl,-esm%so_var_ref%z) 

            case DEFAULT
                ! Assume no variabiliy anomaly (do nothing)
                
        end select
   
        ! remove ronne variability
        if (.FALSE.) then
                where(basins .eq. 1) esm%dto_var = 0.0_wp
        end if

        if (use_ref_atm) then
            ! set atmosphere to reference values
            esm%dts_var = 0.0_wp
            esm%dpr_var = 1.0_wp
        end if
                    
        if (use_ref_ocn) then
            ! set ocean to reference values
            esm%dto_var = 0.0_wp
            esm%dso_var = 0.0_wp
        end if

        if (.FALSE.) then
            ! variability field
            write(*,*) "maxval ts_var", maxval(esm%ts_esm_ref%var(:,:,:,1))
            write(*,*) "maxval pr_var", maxval(esm%pr_esm_ref%var(:,:,:,1))
            write(*,*) "maxval to_var", maxval(esm%to_esm_ref%var(:,:,:,1))
            write(*,*) "maxval so_var", maxval(esm%so_esm_ref%var(:,:,:,1))
                        
            ! anomaly
            write(*,*) "maxval dts_var",  maxval(esm%dts_var)
            write(*,*) "maxval dpr_var",  maxval(esm%dpr_var)
            write(*,*) "maxval dto_var",  maxval(esm%dto_var)
            write(*,*) "maxval dso_vcar", maxval(esm%dso_var)
                        
        end if
    
        return 
    
    end subroutine esm_variability_update

    subroutine esm_forcing_update(esm,mshlf,time,use_esm,time_ref,time_hist,time_proj,time_esm_ref,H_ice,basins,z_bed,f_grnd,z_sl,use_ref_atm,use_ref_ocn)
        ! Update climatic fields. These will be used as bnd conditions for Yelmo.
        ! Output are anomaly fields with respect to a reference field from the ESM.
    
        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        type(marshelf_class),    intent(IN)    :: mshlf
        real(wp), intent(IN) :: time
        logical,  intent(IN) :: use_esm
        real(wp), intent(IN) :: time_ref(2),time_hist(2),time_proj(2),time_esm_ref(2)
        real(wp), intent(IN) :: H_ice(:,:),basins(:,:),z_bed(:,:),f_grnd(:,:),z_sl(:,:)
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 
    
        ! Local variables 
        integer  :: k, m 
        real(wp) :: tmp, anomaly
        character(len=56) :: slice_method 
    
        ! Get slices for current time
        slice_method = "extrap" 
        anomaly = 0.0_wp
            
        ! Initialize anomalies
        esm%dts = 0.0_wp
        esm%dpr = 1.0_wp
        esm%dto = 0.0_wp
        esm%dso = 0.0_wp 

        select case(trim(esm%ctrl_run_type))
            
            case("ctrl","opt")
                ! If ctrl or opt, run only reference field. (Do nothing)
            
            case("transient")
                if (use_esm) then ! Fields loaded from ESM
                    ! === Reference period   ===
                    ! === Atmospheric fields ===
                    call varslice_update(esm%ts_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=12)
                    call varslice_update(esm%pr_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=12)
                    ! === Oceanic fields ===
                    call varslice_update(esm%to_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=1)
                    call varslice_update(esm%so_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=1)
                        
                    ! Interpolate ocean data to the interior of ice shelves
                    if (.TRUE.) then
                        call ocn_variable_extrapolation(esm%to_esm_ref%var(:,:,:,1),H_ice,basins,-esm%to_esm_ref%z,z_bed)
                        call ocn_variable_extrapolation(esm%so_esm_ref%var(:,:,:,1),H_ice,basins,-esm%so_esm_ref%z,z_bed)
                    end if

                    ! === Historical period ===
                    if (time .le. time_hist(2)) then
                        ! === Atmospheric fields === 
                        call varslice_update(esm%ts_hist,[time],method="extrap",rep=1)
                        call varslice_update(esm%pr_hist,[time],method="extrap",rep=1)
                        do m = 1, 12 
                            esm%dts(:,:,m) = esm%ts_hist%var(:,:,1,1)-esm%ts_esm_ref%var(:,:,1,1)
                            esm%dpr(:,:,m) = esm%pr_hist%var(:,:,1,1)/(esm%pr_esm_ref%var(:,:,1,1)+1e-8)
                        end do
                        ! ===   Oceanic fields   ===
                        call varslice_update(esm%to_hist,[time],method="extrap",rep=1)
                        call varslice_update(esm%so_hist,[time],method="extrap",rep=1)
                            
                        ! Interpolate ocean data to the interior of ice shelves
                        call ocn_variable_extrapolation(esm%to_hist%var(:,:,:,1),H_ice,basins,-esm%to_hist%z,z_bed)
                        call ocn_variable_extrapolation(esm%so_hist%var(:,:,:,1),H_ice,basins,-esm%so_hist%z,z_bed)
                            
                        ! Compute the anomaly at the desired depth level
                        call marshelf_interp_shelf(esm%dto,mshlf,esm%to_hist%var(:,:,:,1)-esm%to_esm_ref%var(:,:,:,1),H_ice, &
                                                    z_bed,f_grnd,z_sl,-esm%to_esm_ref%z)
                        call marshelf_interp_shelf(esm%dso,mshlf,esm%so_hist%var(:,:,:,1)-esm%so_esm_ref%var(:,:,:,1),H_ice, &
                                                    z_bed,f_grnd,z_sl,-esm%so_esm_ref%z) 
                        
                    ! === Projection period ===
                    else if (time .ge. time_proj(1)) then

                        ! === Atmospheric fields ===
                        call varslice_update(esm%ts_proj, [time],method="extrap",rep=12)
                        call varslice_update(esm%pr_proj, [time],method="extrap",rep=12)
                        
                        do m = 1, 12
                            esm%dts(:,:,m) = esm%ts_proj%var(:,:,m,1)-esm%ts_esm_ref%var(:,:,m,1)
                            esm%dpr(:,:,m) = esm%pr_proj%var(:,:,m,1)/(esm%pr_esm_ref%var(:,:,m,1)+1e-8)
                        end do    
                        
                        ! ===   Oceanic fields   ===
                        call varslice_update(esm%to_proj,[time],method="extrap",rep=1)
                        call varslice_update(esm%so_proj,[time],method="extrap",rep=1)
                        
                        ! Interpolate ocean data to the interior
                        call ocn_variable_extrapolation(esm%to_proj%var(:,:,:,1),H_ice,basins,-esm%to_proj%z,z_bed)
                        call ocn_variable_extrapolation(esm%so_proj%var(:,:,:,1),H_ice,basins,-esm%so_proj%z,z_bed)
                        
                        ! Compute the anomaly at the desired depth level
                        esm%dto = 0.0_wp
                        esm%dso = 0.0_wp
                        call marshelf_interp_shelf(esm%dto,mshlf,esm%to_proj%var(:,:,:,1)-esm%to_esm_ref%var(:,:,:,1),H_ice, &
                                                    z_bed,f_grnd,z_sl,-esm%to_esm_ref%z)
                        call marshelf_interp_shelf(esm%dso,mshlf,esm%so_proj%var(:,:,:,1)-esm%so_esm_ref%var(:,:,:,1),H_ice, &
                                                    z_bed,f_grnd,z_sl,-esm%so_esm_ref%z)            
                        
                    ! === Reference period ===
                    ! Only used if there is a gap between the historical and projection period
                    else if (time .gt. time_hist(2) .and. time .lt. time_proj(1)) then
                        ! Do nothing

                    end if

                else ! Compute a spatially homogeneous anomaly
                    select case(trim(esm%experiment))
                        case("1pctCO2_transient")
                            ! jablasco: TO DO change
                            ! sergio test
                            if(time .lt. 1995.0) then !2020.0) then
                                ! Do nothing

                            else
                                !anomaly = 0.02*(time - 2020.0)
                                anomaly = (17.0/(2300.0-1995.0))*(time-1995.0)
                                if (anomaly .ge. esm%dT_lim) anomaly=esm%dT_lim
                                esm%dts = esm%f_polar*anomaly
                                esm%dpr = exp(esm%beta_p*esm%f_polar*anomaly)
                                esm%dto = esm%f_ocn*esm%f_polar*anomaly
                                esm%dso = 0.0_wp
                            end if    
                            
                        case DEFAULT
                            write(*,*) "esm_forcing_update:: Error: transient experiment not recognized: "//trim(esm%experiment)
                            stop
                    end select
                end if
                    
    
            case DEFAULT
                write(*,*) "esm_forcing_update:: Error: ctrl_run_type not recognized: "//trim(esm%ctrl_run_type)
                stop 
     
        end select
            
        if (use_ref_atm) then
            ! set atmosphere to reference values
            esm%dts = 0.0_wp
            esm%dpr = 1.0_wp
        end if
            
        if (use_ref_ocn) then
            ! set ocean to reference values
            esm%dto = 0.0_wp
            esm%dso = 0.0_wp
        end if
    
        if (.FALSE.) then
            ! esm ref
            write(*,*) "minval ts_esm", minval(esm%ts_esm_ref%var(:,:,:,1))
            write(*,*) "minval pr_esm", minval(esm%pr_esm_ref%var(:,:,:,1))
            write(*,*) "minval to_esm", minval(esm%to_esm_ref%var(:,:,:,1))
            write(*,*) "minval so_esm", minval(esm%so_esm_ref%var(:,:,:,1))
                                  
            ! proj
            write(*,*) "minval ts_proj", minval(esm%ts_proj%var(:,:,:,1))
            write(*,*) "minval pr_proj", minval(esm%pr_proj%var(:,:,:,1))
            write(*,*) "minval to_proj", minval(esm%to_proj%var(:,:,:,1))
            write(*,*) "minval so_proj", minval(esm%so_proj%var(:,:,:,1))
                                
            ! anomaly
            write(*,*) "minval dts", minval(esm%dts)
            write(*,*) "minval dpr", minval(esm%dpr)
            write(*,*) "minval dto", minval(esm%dto)
            write(*,*) "minval dso", minval(esm%dso)

            ! esm ref
            write(*,*) "maxval ts_esm", maxval(esm%ts_esm_ref%var(:,:,:,1))
            write(*,*) "maxval pr_esm", maxval(esm%pr_esm_ref%var(:,:,:,1))
            write(*,*) "maxval to_esm", maxval(esm%to_esm_ref%var(:,:,:,1))
            write(*,*) "maxval so_esm", maxval(esm%so_esm_ref%var(:,:,:,1))
                          
            ! proj
            write(*,*) "maxval ts_proj", maxval(esm%ts_proj%var(:,:,:,1))
            write(*,*) "maxval pr_proj", maxval(esm%pr_proj%var(:,:,:,1))
            write(*,*) "maxval to_proj", maxval(esm%to_proj%var(:,:,:,1))
            write(*,*) "maxval so_proj", maxval(esm%so_proj%var(:,:,:,1))
                        
            ! anomaly
            write(*,*) "maxval dts", maxval(esm%dts)
            write(*,*) "maxval dpr", maxval(esm%dpr)
            write(*,*) "maxval dto", maxval(esm%dto)
            write(*,*) "maxval dso", maxval(esm%dso)
                        
        end if
    
        return 
    
    end subroutine esm_forcing_update

    ! === varslice wrapper routines with esm specific options ===================

    subroutine varslice_init_nml_esm(vs,filename,group,domain,grid_name,gcm,scenario,verbose)
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
        call varslice_par_load_esm(vs%par,filename,group,domain,grid_name,gcm,scenario,verbose)

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        return 

    end subroutine varslice_init_nml_esm

    subroutine varslice_par_load_esm(par,filename,group,domain,grid_name,gcm,scenario,verbose)

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
        call parse_path_esm(par%filename,gcm,scenario)

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
    
    end subroutine varslice_par_load_esm
    
    subroutine parse_path_esm(path,gcm,scenario)

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN)    :: gcm
        character(len=*), intent(IN)    :: scenario

        call nml_replace(path,"{gcm}",        trim(gcm))
        call nml_replace(path,"{scenario}",   trim(scenario))

        return

    end subroutine parse_path_esm
    
    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    ! === ESM OUTPUT ROUTINES ==========

    subroutine esm_write_init(filename,xc,yc,time,lon,lat,area,map_name,lambda,phi)

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

    end subroutine esm_write_init

    subroutine esm_write_step(filename,file_nml,time)

        implicit none 

        character(len=*),   intent(IN) :: filename
        character(len=*),   intent(IN) :: file_nml 
        real(wp),           intent(IN) :: time 

        ! Local variables 
        integer    :: ncid, n
        real(wp) :: time_prev 
        type(esm_ice_class) :: esm 

        if (.FALSE.) then    
            ! Open the file for writing
            call nc_open(filename,ncid,writable=.TRUE.)

            ! Determine current writing time step 
            n = nc_size(filename,"time",ncid)
            call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
            if (abs(time-time_prev).gt.1e-5) n = n+1 

        end if 
        ! Load up output variable meta information 
        call esm_load_ice_var_info(esm,file_nml,verbose=.TRUE.)
 
        return

    end subroutine esm_write_step

    subroutine esm_load_ice_var_info(esmi,filename,verbose)

        implicit none 

        type(esm_ice_class), intent(OUT) :: esmi 
        character(len=*),    intent(IN)  :: filename 
        logical,             intent(IN)  :: verbose 

        ! Local variables
        integer :: n  
        integer, parameter :: n_variables = 38

        type(esm_ice_var_class) :: v 

        ! First initialize esm object to hold variable meta data 
        if (allocated(esmi%vars)) deallocate(esmi%vars)
        allocate(esmi%vars(n_variables))

        ! Load individual variables by namelist group 
        call ice_var_par_load(esmi%vars(1), filename,var_name="lithk")
        call ice_var_par_load(esmi%vars(2), filename,var_name="orog")
        call ice_var_par_load(esmi%vars(3), filename,var_name="base")
        call ice_var_par_load(esmi%vars(4), filename,var_name="topg")
        call ice_var_par_load(esmi%vars(5), filename,var_name="hfgeoubed")
        call ice_var_par_load(esmi%vars(6), filename,var_name="acabf")
        call ice_var_par_load(esmi%vars(7), filename,var_name="libmassbfgr")
        call ice_var_par_load(esmi%vars(8), filename,var_name="libmassbffl")
        call ice_var_par_load(esmi%vars(9), filename,var_name="dlithkdt")
        call ice_var_par_load(esmi%vars(10),filename,var_name="xvelsurf")
        call ice_var_par_load(esmi%vars(11),filename,var_name="yvelsurf")
        call ice_var_par_load(esmi%vars(12),filename,var_name="zvelsurf")
        call ice_var_par_load(esmi%vars(13),filename,var_name="xvelbase")
        call ice_var_par_load(esmi%vars(14),filename,var_name="yvelbase")
        call ice_var_par_load(esmi%vars(15),filename,var_name="zvelbase")
        call ice_var_par_load(esmi%vars(16),filename,var_name="xvelmean")
        call ice_var_par_load(esmi%vars(17),filename,var_name="yvelmean")
        call ice_var_par_load(esmi%vars(18),filename,var_name="litemptop")
        call ice_var_par_load(esmi%vars(19),filename,var_name="litempbotgr")
        call ice_var_par_load(esmi%vars(20),filename,var_name="litempbotfl")
        call ice_var_par_load(esmi%vars(21),filename,var_name="strbasemag")
        call ice_var_par_load(esmi%vars(22),filename,var_name="licalvf")
        call ice_var_par_load(esmi%vars(23),filename,var_name="lifmassbf")
        call ice_var_par_load(esmi%vars(24),filename,var_name="lifmassbf")
        call ice_var_par_load(esmi%vars(25),filename,var_name="ligroundf")
        call ice_var_par_load(esmi%vars(26),filename,var_name="sftgif")
        call ice_var_par_load(esmi%vars(27),filename,var_name="sftgrf")
        call ice_var_par_load(esmi%vars(28),filename,var_name="sftflf")

        call ice_var_par_load(esmi%vars(29),filename,var_name="lim")
        call ice_var_par_load(esmi%vars(30),filename,var_name="limnsw")
        call ice_var_par_load(esmi%vars(31),filename,var_name="iareagr")
        call ice_var_par_load(esmi%vars(32),filename,var_name="iareafl")
        call ice_var_par_load(esmi%vars(33),filename,var_name="tendacabf")
        call ice_var_par_load(esmi%vars(34),filename,var_name="tendlibmassbf")
        call ice_var_par_load(esmi%vars(35),filename,var_name="tendlibmassbffl")
        call ice_var_par_load(esmi%vars(36),filename,var_name="tendlicalvf")
        call ice_var_par_load(esmi%vars(37),filename,var_name="tendlifmassbf")
        call ice_var_par_load(esmi%vars(38),filename,var_name="tendligroundf")

        if (verbose) then 

            ! === Print summary =========

            write(*,"(a40,a8,a65,a15)") &
                                    "Variable name",    &
                                    "Type",             &
                                    "Standard name",    &
                                    "Unit"
            
            do n = 1, n_variables 
                v = esmi%vars(n)
                write(*,"(a40,a8,a65,a15)") &
                                    trim(v%long_name),      &
                                    trim(v%var_type),       &
                                    trim(v%standard_name),  &
                                    trim(v%units_out) 
            end do 


        end if 

        return 

    end subroutine esm_load_ice_var_info

    subroutine ice_var_par_load(esmv,filename,var_name)
        ! Load parmaeters associated with a given ice variable

        implicit none 

        type(esm_ice_var_class),    intent(OUT) :: esmv 
        character(len=*),           intent(IN)  :: filename 
        character(len=*),           intent(IN)  :: var_name

        ! Local variables 
        character(len=56) :: group 

        group = "esm_out_"//trim(var_name)

        call nml_read(filename,group,"name",            esmv%name)
        call nml_read(filename,group,"long_name",       esmv%long_name)
        call nml_read(filename,group,"var_type",        esmv%var_type)
        call nml_read(filename,group,"standard_name",   esmv%standard_name)
        call nml_read(filename,group,"units_in",        esmv%units_in)
        call nml_read(filename,group,"units_out",       esmv%units_out)
        call nml_read(filename,group,"unit_scale",      esmv%unit_scale)
        call nml_read(filename,group,"unit_offset",     esmv%unit_offset)

        return 

    end subroutine ice_var_par_load

    ! Initlialize allocatable objects
    subroutine esm_allocate(esm,nx,ny)

        implicit none 
    
        type(esm_forcing_class) :: esm 
        integer :: nx, ny 
    
        ! Make object is deallocated
        call esm_deallocate(esm)
    
        ! Allocate variables
        allocate(esm%t2m(nx,ny,12))
        allocate(esm%pr(nx,ny,12))
        allocate(esm%dts(nx,ny,12))
        allocate(esm%dpr(nx,ny,12))
        allocate(esm%dto(nx,ny))
        allocate(esm%dso(nx,ny))
        allocate(esm%dts_var(nx,ny,12))
        allocate(esm%dpr_var(nx,ny,12))
        allocate(esm%dto_var(nx,ny))
        allocate(esm%dso_var(nx,ny))
        allocate(esm%t2m_ann(nx,ny))
        allocate(esm%t2m_sum(nx,ny))
        allocate(esm%pr_ann(nx,ny))        

        return
    
    end subroutine esm_allocate
    
    subroutine esm_deallocate(esm)
    
        implicit none 
    
        type(esm_forcing_class) :: esm
    
            ! Allocate state objects
            if (allocated(esm%t2m))     deallocate(esm%t2m)
            if (allocated(esm%pr))      deallocate(esm%pr)
            if (allocated(esm%dts))     deallocate(esm%dts)
            if (allocated(esm%dpr))     deallocate(esm%dpr)
            if (allocated(esm%dto))     deallocate(esm%dto)
            if (allocated(esm%dso))     deallocate(esm%dso)
            if (allocated(esm%dts_var)) deallocate(esm%dts_var)
            if (allocated(esm%dpr_var)) deallocate(esm%dpr_var)
            if (allocated(esm%dto_var)) deallocate(esm%dto_var)
            if (allocated(esm%dso_var)) deallocate(esm%dso_var)
            if (allocated(esm%t2m_sum)) deallocate(esm%t2m_sum)
            if (allocated(esm%t2m_ann)) deallocate(esm%t2m_ann)
            if (allocated(esm%pr_ann))  deallocate(esm%pr_ann)

            return
    
        end subroutine esm_deallocate
    
end module esm

