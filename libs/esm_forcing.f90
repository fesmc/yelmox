module esm_forcing
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
        type(varslice_class)   :: smb_ref    
        type(varslice_class)   :: pr_ref
        
        ! Ocean
        type(varslice_class)   :: to_ref, to_ref_src
        type(varslice_class)   :: so_ref

        ! ===  Variability climatologies ===
        ! Atmosphere
        type(varslice_class)   :: ts_var
        type(varslice_class)   :: pr_var
        type(varslice_class)   :: smb_var
        type(varslice_class)   :: ts_var_ref
        type(varslice_class)   :: pr_var_ref
        type(varslice_class)   :: smb_var_ref

        ! Ocean
        type(varslice_class)   :: to_var
        type(varslice_class)   :: so_var
        type(varslice_class)   :: to_var_ref
        type(varslice_class)   :: so_var_ref
        
        ! ===      ESM      ===
        ! Atmospheric fields
        type(varslice_class)   :: ts_esm_ref 
        type(varslice_class)   :: pr_esm_ref 
        type(varslice_class)   :: smb_esm_ref

        type(varslice_class)   :: ts_hist 
        type(varslice_class)   :: pr_hist
        type(varslice_class)   :: smb_hist
        type(varslice_class)   :: dsmbdz_hist 

        type(varslice_class)   :: ts_proj
        type(varslice_class)   :: pr_proj
        type(varslice_class)   :: smb_proj
        type(varslice_class)   :: dsmbdz_proj

        ! Oceanic fields 
        type(varslice_class)   :: to_esm_ref
        type(varslice_class)   :: so_esm_ref
                
        type(varslice_class)   :: to_hist
        type(varslice_class)   :: so_hist
        type(varslice_class)   :: Qd_hist

        type(varslice_class)   :: to_proj
        type(varslice_class)   :: so_proj
        type(varslice_class)   :: Qd_proj

        ! General fields 
        type(varslice_class)   :: zs_ref
        type(varslice_class)   :: zs_esm_ref
        
        ! === Diagnostic fields ===

        ! === Atmosphere ===
        ! Monthly fields (but right now constant anomaly between months)
        real(wp), allocatable :: t2m(:,:,:)       ! Monthly surface temperature [K]
        real(wp), allocatable :: pr(:,:,:)        ! Monthly precipitation [mm/yr]
        real(wp), allocatable :: smb(:,:,:)       ! Monthly SMB [mm/yr]

        ! Anomalies
        real(wp), allocatable :: dts(:,:,:)       ! Surface temperature anomaly [K]
        real(wp), allocatable :: dpr(:,:,:)       ! Precipitation relative anomaly [%]
        real(wp), allocatable :: dsmb(:,:,:)      ! SMB anomaly [mm/yr]
        real(wp), allocatable :: dsmbdz(:,:)      ! SMB height anomaly [mm/yr*m]
        real(wp), allocatable :: dts_var(:,:,:)   ! Surface temperature anomaly variability[K]
        real(wp), allocatable :: dpr_var(:,:,:)   ! Precipitation relative anomaly variability [%]
        real(wp), allocatable :: dsmb_var(:,:,:)  ! SMB anomaly variability [m/yr]

        ! ===    Ocean   ===
        real(wp), allocatable :: dto(:,:)         ! Surface temperature anomaly [K]
        real(wp), allocatable :: dso(:,:)         ! Precipitation relative anomaly [%]
        real(wp), allocatable :: dto_var(:,:)     ! Surface temperature anomaly variability [K]
        real(wp), allocatable :: dso_var(:,:)     ! Precipitation relative anomaly variability [%]
        real(wp), allocatable :: Qd_ann(:,:)      ! Annual mean subglacial discharge [m3/s]
        real(wp), allocatable :: Qd_sum(:,:)      ! Summer mean subglacial discharge [m3/s]
        
        ! === Mean fields ===
        real(wp), allocatable :: t2m_sum(:,:)     ! Summer surface temperature [K]
        real(wp), allocatable :: t2m_ann(:,:)     ! Annual surface temperature [K]
        real(wp), allocatable :: pr_ann(:,:)      ! Annual precipitation [mm/yr]
        real(wp), allocatable :: smb_ann(:,:)     ! Annual SMB [mm/yr]

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
    
    private
    public :: esm_forcing_class
    public :: esm_experiment_class

    ! General routines
    public :: esm_forcing_init
    public :: esm_forcing_update
    public :: esm_variability_update
    public :: esm_clim_update
    
    public :: esm_write_init

contains
    
    subroutine esm_forcing_init(esm,filename,domain,grid_name,run_type,gcm,experiment, &
                                use_esm,use_smb,use_var,&
                                use_hist,time_hist,&
                                use_proj,time_proj)

        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN), optional :: run_type, gcm, experiment
        logical,          intent(IN), optional :: use_esm, use_smb, use_var
        logical,          intent(IN), optional :: use_hist, use_proj
        real(wp),         intent(IN), optional :: time_hist(2), time_proj(2)
    
        ! Local variables 
        character(len=256) :: group_prefix 
    
        ! Reference climatology
        character(len=256) :: grp_ts_ref  
        character(len=256) :: grp_pr_ref
        character(len=256) :: grp_smb_ref
        character(len=256) :: grp_zs_ref
        character(len=256) :: grp_to_ref 
        character(len=256) :: grp_so_ref
        ! Variability period
        character(len=256) :: grp_ts_var  
        character(len=256) :: grp_pr_var
        character(len=256) :: grp_smb_var
        character(len=256) :: grp_to_var 
        character(len=256) :: grp_so_var
        ! ESM Reference period
        character(len=256) :: grp_ts_esm_ref  
        character(len=256) :: grp_pr_esm_ref
        character(len=256) :: grp_smb_esm_ref
        character(len=256) :: grp_zs_esm_ref
        character(len=256) :: grp_to_esm_ref 
        character(len=256) :: grp_so_esm_ref

        ! ESM Historical period 
        character(len=256) :: grp_ts_hist 
        character(len=256) :: grp_pr_hist
        character(len=256) :: grp_smb_hist
        character(len=256) :: grp_dsmbdz_hist
        character(len=256) :: grp_to_hist 
        character(len=256) :: grp_so_hist
        character(len=256) :: grp_Qd_hist

        ! ESM Projection period 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
        character(len=256) :: grp_smb_proj
        character(len=256) :: grp_dsmbdz_proj
        character(len=256) :: grp_to_proj 
        character(len=256) :: grp_so_proj
        character(len=256) :: grp_Qd_proj   ! used for Greenland

        integer  :: iloc, k
        real(wp) :: tmp
        real(wp) :: time_par_ref(4),time_par_hist(4),time_par_proj(4),time_par_var(4)
        character(len=256) :: esm_subs(2,2)   ! {gcm}/{experiment} path substitutions for varslice
    
        ! Define the current experiment characteristics
        esm%ctrl_run_type = trim(run_type)
        esm%gcm           = trim(gcm)
        esm%experiment    = trim(experiment)

        ! {gcm}/{experiment} placeholder substitutions applied to ESM forcing
        ! file paths by the shared varslice reader (via its `subs` argument).
        esm_subs(1,1) = "gcm"
        esm_subs(1,2) = trim(esm%gcm)
        esm_subs(2,1) = "experiment"
        esm_subs(2,2) = trim(esm%experiment)
    
        write(*,*)
        write(*,*) "esm_forcing_init:: summary"
        write(*,*) "ctrl_run_type: ", trim(esm%ctrl_run_type)
        write(*,*) "gcm:           ", trim(esm%gcm)
        write(*,*) "experiment:    ", trim(esm%experiment)
        write(*,*) 
    
        group_prefix = "gcm_"
    
        ! Reference climatology
        grp_ts_ref       = trim(group_prefix)//"ts_ref"
        grp_pr_ref       = trim(group_prefix)//"pr_ref"
        grp_smb_ref      = trim(group_prefix)//"smb_ref"
        grp_zs_ref       = trim(group_prefix)//"zs_ref"
        grp_to_ref       = trim(group_prefix)//"to_ref"
        grp_so_ref       = trim(group_prefix)//"so_ref"   

        ! Variability climatology
        grp_ts_var       = trim(group_prefix)//"ts_var"
        grp_pr_var       = trim(group_prefix)//"pr_var"
        grp_smb_var      = trim(group_prefix)//"smb_var"
        grp_to_var       = trim(group_prefix)//"to_var"
        grp_so_var       = trim(group_prefix)//"so_var"

        ! ESM Reference climatology (to compute anomalies)
        grp_ts_esm_ref   = trim(group_prefix)//"ts_esm_ref"
        grp_pr_esm_ref   = trim(group_prefix)//"pr_esm_ref"
        grp_smb_esm_ref  = trim(group_prefix)//"smb_esm_ref"
        grp_zs_esm_ref   = trim(group_prefix)//"zs_esm_ref"
        grp_to_esm_ref   = trim(group_prefix)//"to_esm_ref"
        grp_so_esm_ref   = trim(group_prefix)//"so_esm_ref" 

        ! ESM Historical sims 
        grp_ts_hist     = trim(group_prefix)//"ts_hist"
        grp_pr_hist     = trim(group_prefix)//"pr_hist"
        grp_smb_hist    = trim(group_prefix)//"smb_hist"
        grp_dsmbdz_hist = trim(group_prefix)//"dsmbdz_hist"
        grp_to_hist     = trim(group_prefix)//"to_hist"
        grp_so_hist     = trim(group_prefix)//"so_hist"
        grp_Qd_hist     = trim(group_prefix)//"sgd_hist"

        ! ESM projected sims
        grp_ts_proj     = trim(group_prefix)//"ts_proj"
        grp_pr_proj     = trim(group_prefix)//"pr_proj"
        grp_smb_proj    = trim(group_prefix)//"smb_proj"
        grp_dsmbdz_proj = trim(group_prefix)//"dsmbdz_proj"
        grp_to_proj     = trim(group_prefix)//"to_proj"
        grp_so_proj     = trim(group_prefix)//"so_proj"
        grp_Qd_proj     = trim(group_prefix)//"sgd_proj"         
     
        ! Climatology
        ! Reference period
        call varslice_init_nml(esm%ts_ref, filename, trim(grp_ts_ref), domain, grid_name, subs=esm_subs)
        if (use_smb) then
            call varslice_init_nml(esm%smb_ref, filename, trim(grp_smb_ref), domain, grid_name, subs=esm_subs)
        else
            call varslice_init_nml(esm%pr_ref, filename, trim(grp_pr_ref), domain, grid_name, subs=esm_subs)
        end if
        call varslice_init_nml(esm%zs_ref, filename, trim(grp_zs_ref), domain, grid_name, subs=esm_subs)
        call varslice_init_nml(esm%to_ref, filename, trim(grp_to_ref), domain, grid_name, subs=esm_subs)
        call varslice_init_nml(esm%so_ref, filename, trim(grp_so_ref), domain, grid_name, subs=esm_subs)

        ! Initialize variables at other grids if source grid is different
        !call varslice_init_nml(esm%to_ref_src, filename, trim(grp_to_ref), domain, grid_name, subs=esm_subs)

        if (use_var) then
            ! Variability
            ! Transient dependent field
            call varslice_init_nml(esm%ts_var, filename, trim(grp_ts_var), domain, grid_name, subs=esm_subs)
            if (use_smb) then
                call varslice_init_nml(esm%smb_var, filename, trim(grp_smb_var), domain, grid_name, subs=esm_subs)
            else
                call varslice_init_nml(esm%pr_var, filename, trim(grp_pr_var), domain, grid_name, subs=esm_subs)
            end if
            call varslice_init_nml(esm%to_var, filename, trim(grp_to_var), domain, grid_name, subs=esm_subs)
            call varslice_init_nml(esm%so_var, filename, trim(grp_so_var), domain, grid_name, subs=esm_subs)
            ! Reference period
            call varslice_init_nml(esm%ts_var_ref, filename, trim(grp_ts_var), domain, grid_name, subs=esm_subs)
            if (use_smb) then
                call varslice_init_nml(esm%smb_var_ref, filename, trim(grp_smb_var), domain, grid_name, subs=esm_subs)
            else
                call varslice_init_nml(esm%pr_var_ref, filename, trim(grp_pr_var), domain, grid_name, subs=esm_subs)
            end if
            call varslice_init_nml(esm%to_var_ref, filename, trim(grp_to_var), domain, grid_name, subs=esm_subs)
            call varslice_init_nml(esm%so_var_ref, filename, trim(grp_so_var), domain, grid_name, subs=esm_subs)
        end if

        ! Transient dependent fields
        if (trim(esm%ctrl_run_type) .eq. "transient" .and. trim(esm%experiment) .ne. "ctrl") then
            ! ESM reference period
            if (use_esm) then
                call varslice_init_nml(esm%ts_esm_ref, filename, trim(grp_ts_esm_ref), domain, grid_name, subs=esm_subs)
                if (use_smb) then
                    call varslice_init_nml(esm%smb_esm_ref, filename, trim(grp_smb_esm_ref), domain, grid_name, subs=esm_subs)
                else
                    call varslice_init_nml(esm%pr_esm_ref, filename, trim(grp_pr_esm_ref), domain, grid_name, subs=esm_subs)
                end if
                call varslice_init_nml(esm%to_esm_ref, filename, trim(grp_to_esm_ref), domain, grid_name, subs=esm_subs)
                call varslice_init_nml(esm%so_esm_ref, filename, trim(grp_so_esm_ref), domain, grid_name, subs=esm_subs)
                call varslice_init_nml(esm%zs_esm_ref, filename, trim(grp_zs_esm_ref), domain, grid_name, subs=esm_subs)

                ! ESM historical period
                if (use_hist) then
                    call varslice_init_nml(esm%ts_hist, filename,trim(grp_ts_hist), domain, grid_name, subs=esm_subs)
                    if (use_smb) then
                        call varslice_init_nml(esm%smb_hist, filename, trim(grp_smb_hist), domain, grid_name, subs=esm_subs)
                        call varslice_init_nml(esm%dsmbdz_hist, filename, trim(grp_dsmbdz_hist), domain, grid_name, subs=esm_subs)
                    else
                        call varslice_init_nml(esm%pr_hist, filename, trim(grp_pr_hist), domain, grid_name, subs=esm_subs)
                    end if
                    call varslice_init_nml(esm%to_hist, filename,trim(grp_to_hist), domain, grid_name, subs=esm_subs)
                    call varslice_init_nml(esm%so_hist, filename,trim(grp_so_hist), domain, grid_name, subs=esm_subs)
                    if (trim(domain).eq."Greenland") then
                        call varslice_init_nml(esm%Qd_hist, filename,trim(grp_Qd_hist), domain,grid_name,subs=esm_subs)
                    end if
                end if
                
                ! ESM projection period
                if (use_proj) then
                    ! atm
                    call varslice_init_nml(esm%ts_proj, filename,trim(grp_ts_proj), domain,grid_name,subs=esm_subs)
                    if (use_smb) then
                        call varslice_init_nml(esm%smb_proj, filename, trim(grp_smb_proj), domain, grid_name, subs=esm_subs)
                        call varslice_init_nml(esm%dsmbdz_proj, filename, trim(grp_dsmbdz_proj), domain, grid_name, subs=esm_subs)
                    else
                        call varslice_init_nml(esm%pr_proj, filename, trim(grp_pr_proj), domain, grid_name, subs=esm_subs)
                    end if
                    ! ocean
                    call varslice_init_nml(esm%to_proj, filename,trim(grp_to_proj), domain,grid_name,subs=esm_subs)
                    call varslice_init_nml(esm%so_proj, filename,trim(grp_so_proj), domain,grid_name,subs=esm_subs)
                    if (trim(domain).eq."Greenland") then
                        call varslice_init_nml(esm%Qd_proj, filename,trim(grp_Qd_proj), domain,grid_name,subs=esm_subs)
                    end if
                end if
            end if
        end if

        ! Allocate objects (use a time independent field)
        call varslice_update(esm%zs_ref)
        call esm_allocate(esm,size(esm%zs_ref%var,1),size(esm%zs_ref%var,2))
        write(*,*) "dim1, dim2", size(esm%zs_ref%var,1),size(esm%zs_ref%var,2)

        return 
    
    end subroutine esm_forcing_init

    subroutine esm_clim_update(esm,z_srf_ylm,time,time_ref,use_smb,domain,grid_name)
        ! Routine to update reference climatology to the specific Antarctic elevation and ocean (neccessary?)

        implicit none

        type(esm_forcing_class), intent(INOUT) :: esm
        real(wp),                intent(IN)    :: z_srf_ylm(:,:)
        real(wp),                intent(IN)    :: time
        real(wp),                intent(IN)    :: time_ref(2)
        logical,                 intent(IN)    :: use_smb
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
        ! Reference fields may be a static 12-month climatology (with_time=False,
        ! e.g. Greenland MAR) or a transient monthly series (with_time=True, e.g.
        ! Antarctica RACMO). range_mean over [time_ref] collapses a transient
        ! series to a 12-month climatology (rep=12) and is a no-op on static data.
        call varslice_update(esm%ts_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        if (use_smb) then
            call varslice_update(esm%smb_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        else
            call varslice_update(esm%pr_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        end if
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

        ! Convert atmospheric fields to model elevation
        do m = 1,12 
            if(south) then ! Southern Hemisphere
                lapse = (esm%lapse(1)+(esm%lapse(2)-esm%lapse(1))*cos(2*pi*(m*30.4375-30.4375)/365.25))
            else ! Northern Hemisphere
                lapse = (esm%lapse(1)+(esm%lapse(1)-esm%lapse(2))*cos(2*pi*(m*30.4375-30.4375)/365.25))
            end if    
            esm%t2m(:,:,m) = esm%ts_ref%var(:,:,m,1) + lapse*(esm%zs_ref%var(:,:,1,1)-z_srf_ylm)
            if (use_smb) then
                esm%smb(:,:,m) = esm%smb_ref%var(:,:,m,1) !No model elevation changes for SMB
            else
                esm%pr(:,:,m)  = esm%pr_ref%var(:,:,m,1) * exp(esm%beta_p*lapse*(esm%zs_ref%var(:,:,1,1)-z_srf_ylm))
            end if
        end do

        ! Compute diagnostic fields
        esm%t2m_ann = sum(esm%t2m,dim=3) / 12.0
        if (use_smb) then
            esm%smb_ann  = sum(esm%smb,dim=3) / 12.0
        else
            esm%pr_ann  = sum(esm%pr,dim=3) / 12.0
        end if

        if(south) then
            esm%t2m_sum = (esm%t2m(:,:,1)+esm%t2m(:,:,2)+esm%t2m(:,:,12)) / 3.0
        else
            esm%t2m_sum = (esm%t2m(:,:,6)+esm%t2m(:,:,7)+esm%t2m(:,:,8)) / 3.0
        end if      

        return

    end subroutine esm_clim_update

    subroutine esm_variability_update(esm,mshlf,time,dtt,clim_var,time_ref,H_ice,basins,z_bed,f_grnd,z_sl,use_var,use_ref_atm,use_ref_ocn)
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
        logical,  intent(IN) :: use_var
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

        ! jablasco: skip for ismip7 (improve)
        if (use_var) then 
        ! Variability reference
        ! === Atmospheric fields ===
        call varslice_update(esm%ts_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        call varslice_update(esm%pr_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=12)
        ! ===   Oceanic fields   ===
        call varslice_update(esm%to_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)
        call varslice_update(esm%so_var_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1) 
        
        call ocn_variable_extrapolation(esm%to_var_ref%var(:,:,:,1),H_ice,basins,-esm%to_var_ref%z,z_bed)
        call ocn_variable_extrapolation(esm%so_var_ref%var(:,:,:,1),H_ice,basins,-esm%so_var_ref%z,z_bed)

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

            case("snapshot","snap")
                ! Cycle through the selected snapshot
                year_rand = time
                write(*,*) "year_snap = ", year_rand    

            case("red_noise")
                ! Select a random year but based on red noise
                ! TO DO
                write(*,*) "Red noise method not available yet."
                STOP

            case DEFAULT
                ! Do nothing for now
        end select

        select case(trim(clim_var))
            case("random","rand","white_noise","historic","hist","red_noise","snapshot","snap")
                ! === Atmospheric fields === 
                call varslice_update(esm%ts_var,[year_rand],method="extrap",rep=12)
                call varslice_update(esm%pr_var,[year_rand],method="extrap",rep=12)
                do m = 1, 12 
                    esm%dts_var(:,:,m) = esm%ts_var%var(:,:,m,1)-esm%ts_var_ref%var(:,:,m,1)
                    esm%dpr_var(:,:,m) = esm%pr_var%var(:,:,m,1)/(esm%pr_var_ref%var(:,:,m,1)+1e-8)
                end do
                ! ===   Oceanic fields   ===
                call varslice_update(esm%to_var,[year_rand],method="extrap",rep=1)
                call varslice_update(esm%so_var,[year_rand],method="extrap",rep=1)

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
   
        ! routine to rome variability in speciic basins. TO DO
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
        end if

        return 
    
    end subroutine esm_variability_update

    subroutine esm_forcing_update(esm,mshlf,time,use_esm,time_ref,time_hist,time_proj,time_esm_ref,&
                                  domain,H_ice,basins,z_bed,f_grnd,z_sl,use_smb,use_ref_atm,use_ref_ocn)
        ! Update climatic fields. These will be used as bnd conditions for Yelmo.
        ! Output are anomaly fields with respect to a reference field from the ESM.
    
        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        type(marshelf_class),    intent(IN)    :: mshlf
        real(wp), intent(IN) :: time
        logical,  intent(IN) :: use_esm
        real(wp), intent(IN) :: time_ref(2),time_hist(2),time_proj(2),time_esm_ref(2)
        character(len=*), intent(IN) :: domain
        real(wp), intent(IN) :: H_ice(:,:),basins(:,:),z_bed(:,:),f_grnd(:,:),z_sl(:,:)
        logical,  intent(IN) :: use_smb  
        logical,  intent(IN), optional :: use_ref_atm, use_ref_ocn
    
        ! Local variables 
        integer  :: k, m 
        real(wp) :: tmp, anomaly
        character(len=56) :: slice_method 
    
        ! Get slices for current time
        slice_method = "extrap" 
        anomaly = 0.0_wp
            
        ! Initialize anomalies
        esm%dts    = 0.0_wp
        esm%dpr    = 1.0_wp
        esm%dsmb   = 0.0_wp
        esm%dsmbdz = 0.0_wp
        esm%dto    = 0.0_wp
        esm%dso    = 0.0_wp 
        esm%Qd_ann = 0.0_wp
        esm%Qd_sum = 0.0_wp

        select case(trim(esm%ctrl_run_type))
            
            case("ctrl","opt","spinup")
                ! ctrl/opt/spinup: run only the reference climatology. (Do nothing)
            
            case("transient")
                if (use_esm) then ! Fields loaded from ESM
                    ! === Reference period   ===
                    ! === Atmospheric fields ===
                    call varslice_update(esm%ts_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=12)
                    if (use_smb) then
                        call varslice_update(esm%smb_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=12)
                    else
                        call varslice_update(esm%pr_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=12)
                    end if
                    ! === Oceanic fields ===
                    call varslice_update(esm%to_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=1)
                    call varslice_update(esm%so_esm_ref,[time_esm_ref(1),time_esm_ref(2)],method="range_mean",rep=1)
                        
                    ! Interpolate ocean data to the interior of ice shelves
                    if (mshlf%par%extrap_shlf) then
                        call ocn_variable_extrapolation(esm%to_esm_ref%var(:,:,:,1),H_ice,basins,-esm%to_esm_ref%z,z_bed)
                        call ocn_variable_extrapolation(esm%so_esm_ref%var(:,:,:,1),H_ice,basins,-esm%so_esm_ref%z,z_bed)
                    end if

                    ! === Historical period ===
                    if (time .le. time_hist(2)) then
                        ! === Atmospheric fields === 
                        call varslice_update(esm%ts_hist,[time],method="extrap",rep=12)
                        if (use_smb) then
                            call varslice_update(esm%smb_hist,[time],method="extrap",rep=12)
                            call varslice_update(esm%dsmbdz_hist,[time],method="extrap",rep=1)
                            esm%dsmbdz(:,:) = esm%dsmbdz_hist%var(:,:,1,1)
                        else
                            call varslice_update(esm%pr_hist,[time],method="extrap",rep=12)
                        end if

                        do m = 1, 12 
                            esm%dts(:,:,m) = esm%ts_hist%var(:,:,m,1)-esm%ts_esm_ref%var(:,:,m,1)
                            if (use_smb) then
                                esm%dsmb(:,:,m) = esm%smb_hist%var(:,:,m,1)-esm%smb_ref%var(:,:,m,1)
                            else
                                esm%dpr(:,:,m) = esm%pr_hist%var(:,:,m,1)/(esm%pr_esm_ref%var(:,:,m,1)+1e-8)
                            end if
                        end do
                        ! ===   Oceanic fields   ===
                        call varslice_update(esm%to_hist,[time],method="extrap",rep=1)
                        call varslice_update(esm%so_hist,[time],method="extrap",rep=1)
                        if (trim(domain).eq."Greenland") then
                            call varslice_update(esm%Qd_hist,[time],method="extrap",rep=12)
                            esm%Qd_ann = sum(esm%Qd_hist%var(:,:,:,1),dim=3) / 12.0
                            esm%Qd_sum = (esm%Qd_hist%var(:,:,6,1)+esm%Qd_hist%var(:,:,7,1)+esm%Qd_hist%var(:,:,8,1)) / 3.0
                        end if
                                         
                        if (mshlf%par%extrap_shlf) then
                            ! Extrapolate ocean data to the interior of ice shelves
                            call ocn_variable_extrapolation(esm%to_hist%var(:,:,:,1),H_ice,basins,-esm%to_hist%z,z_bed)
                            call ocn_variable_extrapolation(esm%so_hist%var(:,:,:,1),H_ice,basins,-esm%so_hist%z,z_bed)
                        end if

                        ! Compute the anomaly at the desired depth level
                        call marshelf_interp_shelf(esm%dto,mshlf,esm%to_hist%var(:,:,:,1)-esm%to_esm_ref%var(:,:,:,1),H_ice, &
                                                    z_bed,f_grnd,z_sl,-esm%to_esm_ref%z)
                        call marshelf_interp_shelf(esm%dso,mshlf,esm%so_hist%var(:,:,:,1)-esm%so_esm_ref%var(:,:,:,1),H_ice, &
                                                    z_bed,f_grnd,z_sl,-esm%so_esm_ref%z) 
                        
                    ! === Projection period ===
                    else if (time .ge. time_proj(1)) then

                        ! === Atmospheric fields ===
                        call varslice_update(esm%ts_proj, [time],method="extrap",rep=12)
                        if (use_smb) then
                            call varslice_update(esm%smb_proj, [time],method="extrap",rep=12)
                            call varslice_update(esm%dsmbdz_proj, [time],method="extrap",rep=1)
                            esm%dsmbdz(:,:) = esm%dsmbdz_proj%var(:,:,1,1)
                        else
                            call varslice_update(esm%pr_proj, [time],method="extrap",rep=12)
                        end if

                        do m = 1, 12
                            esm%dts(:,:,m) = esm%ts_proj%var(:,:,m,1)-esm%ts_esm_ref%var(:,:,m,1)
                            if (use_smb) then
                                esm%dsmb(:,:,m) = esm%smb_proj%var(:,:,m,1)-esm%smb_esm_ref%var(:,:,m,1)
                            else
                                esm%dpr(:,:,m) = esm%pr_proj%var(:,:,m,1)/(esm%pr_esm_ref%var(:,:,m,1)+1e-8)
                            end if
                        end do    
                        
                        ! ===   Oceanic fields   ===
                        call varslice_update(esm%to_proj,[time],method="extrap",rep=1)
                        call varslice_update(esm%so_proj,[time],method="extrap",rep=1)  
                        if (trim(domain).eq."Greenland") then
                            call varslice_update(esm%Qd_proj,[time],method="extrap",rep=12)
                            esm%Qd_ann = sum(esm%Qd_proj%var(:,:,:,1),dim=3) / 12.0
                            esm%Qd_sum = (esm%Qd_proj%var(:,:,6,1)+esm%Qd_proj%var(:,:,7,1)+esm%Qd_proj%var(:,:,8,1)) / 3.0
                        end if
                            
                        if (mshlf%par%extrap_shlf) then
                            ! Interpolate ocean data to the interior
                            call ocn_variable_extrapolation(esm%to_proj%var(:,:,:,1),H_ice,basins,-esm%to_proj%z,z_bed)
                            call ocn_variable_extrapolation(esm%so_proj%var(:,:,:,1),H_ice,basins,-esm%so_proj%z,z_bed)
                        end if

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
                        case("1pctCO2")
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
                        
                        case("ctrl")
                                ! Do nothing

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
            esm%dts  = 0.0_wp
            esm%dpr  = 1.0_wp
            esm%dsmb = 0.0_wp
            esm%dsmbdz = 0.0_wp
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
        allocate(esm%smb(nx,ny,12))
        allocate(esm%dts(nx,ny,12))
        allocate(esm%dpr(nx,ny,12))
        allocate(esm%dsmb(nx,ny,12))
        allocate(esm%dsmbdz(nx,ny))
        allocate(esm%dto(nx,ny))
        allocate(esm%dso(nx,ny))
        allocate(esm%Qd_ann(nx,ny))
        allocate(esm%Qd_sum(nx,ny))
        allocate(esm%dts_var(nx,ny,12))
        allocate(esm%dpr_var(nx,ny,12))
        allocate(esm%dto_var(nx,ny))
        allocate(esm%dso_var(nx,ny))
        allocate(esm%t2m_ann(nx,ny))
        allocate(esm%t2m_sum(nx,ny))
        allocate(esm%pr_ann(nx,ny)) 
        allocate(esm%smb_ann(nx,ny))       

        return
    
    end subroutine esm_allocate
    
    subroutine esm_deallocate(esm)
    
        implicit none 
    
        type(esm_forcing_class) :: esm
    
            ! Allocate state objects
            if (allocated(esm%t2m))     deallocate(esm%t2m)
            if (allocated(esm%pr))      deallocate(esm%pr)
            if (allocated(esm%smb))     deallocate(esm%smb)
            if (allocated(esm%dts))     deallocate(esm%dts)
            if (allocated(esm%dpr))     deallocate(esm%dpr)
            if (allocated(esm%dsmb))    deallocate(esm%dsmb)
            if (allocated(esm%dsmbdz))  deallocate(esm%dsmbdz)
            if (allocated(esm%dto))     deallocate(esm%dto)
            if (allocated(esm%dso))     deallocate(esm%dso)
            if (allocated(esm%dts_var)) deallocate(esm%dts_var)
            if (allocated(esm%dpr_var)) deallocate(esm%dpr_var)
            if (allocated(esm%dto_var)) deallocate(esm%dto_var)
            if (allocated(esm%dso_var)) deallocate(esm%dso_var)
            if (allocated(esm%dso_var)) deallocate(esm%dso_var)
            if (allocated(esm%Qd_ann))  deallocate(esm%Qd_ann)
            if (allocated(esm%Qd_sum))  deallocate(esm%Qd_sum)
            if (allocated(esm%t2m_sum)) deallocate(esm%t2m_sum)
            if (allocated(esm%t2m_ann)) deallocate(esm%t2m_ann)
            if (allocated(esm%pr_ann))  deallocate(esm%pr_ann)
            if (allocated(esm%smb_ann)) deallocate(esm%smb_ann)

            return
    
        end subroutine esm_deallocate

        ! Internal functions
        function int_to_str(i) result(str)
            integer, intent(in) :: i
            character(len=20)   :: str
            write(str, '(i10)') i
            str = trim(str)
        end function
    
end module esm_forcing

