module esm_forcing
    ! This module contains routines that help with performing the esm suite
    ! of experiments. 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml  
    use ncio 
    use varslice
    use coordinates, only : grid_class
    use marine_shelf

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Define default missing value 
    real(wp), parameter :: mv = -9999.0_wp 

    ! Class for holding ice-forcing data from esm archives
    type esm_forcing_class
        
        ! Experiment information
        character(len=256)     :: gcm 
        character(len=256)     :: experiment 
        character(len=256)     :: domain 
        character(len=256)     :: ctrl_run_type
        real(wp)               :: lapse(2)
        real(wp)               :: beta_p
        real(wp)               :: f_ocn
        real(wp)               :: f_polar
        real(wp)               :: dT_lim

        ! === Reference climatology ===
        ! Atmosphere
        type(varslice_class)   :: ts_ref
        type(varslice_class)   :: smb_ref    
        type(varslice_class)   :: pr_ref
        
        ! Ocean
        type(varslice_class)   :: to_ref
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

        ! ===    Ocean   ===
        real(wp), allocatable :: dto(:,:)         ! Surface temperature anomaly [K]
        real(wp), allocatable :: dso(:,:)         ! Precipitation relative anomaly [%]
        real(wp), allocatable :: dto_var(:,:)     ! Surface temperature anomaly variability [K]
        real(wp), allocatable :: dso_var(:,:)     ! Precipitation relative anomaly variability [%]
        real(wp), allocatable :: Qd_ann(:,:)      ! Annual mean subglacial discharge [m3/s]
        real(wp), allocatable :: Qd_sum(:,:)      ! Summer mean subglacial discharge [m3/s]
        logical :: with_sgd_hist = .false.        ! Qd from the sgd_hist group (Greenland; its filename /= "none")
        logical :: with_sgd_proj = .false.        ! Qd from the sgd_proj group (Greenland; its filename /= "none")
        
        ! === Mean fields ===
        real(wp), allocatable :: t2m_ann(:,:)     ! Annual surface temperature [K]
        real(wp), allocatable :: pr_ann(:,:)      ! Annual precipitation [mm/d]
        real(wp), allocatable :: smb_ann(:,:)     ! Annual SMB [mm/yr]

    end type

    private
    public :: esm_forcing_class

    ! General routines
    public :: esm_forcing_init
    public :: esm_forcing_update
    public :: esm_variability_update
    public :: esm_clim_update
    public :: esm_lapse_monthly

contains
    
    subroutine esm_forcing_init(esm,filename,domain,grid,run_type,gcm,experiment, &
                                use_esm,use_smb,use_var,use_hist,use_proj)

        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        type(grid_class), intent(IN) :: grid         ! the climate grid (target of remapped groups)
        character(len=*), intent(IN) :: run_type, gcm, experiment
        logical,          intent(IN) :: use_esm, use_smb, use_var
        logical,          intent(IN) :: use_hist, use_proj
    
        ! Local variables 
        character(len=256) :: grid_name
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

        character(len=256) :: esm_subs(2,2)   ! {gcm}/{experiment} path substitutions for varslice
    
        ! Define the current experiment characteristics
        esm%ctrl_run_type = trim(run_type)
        esm%gcm           = trim(gcm)
        esm%experiment    = trim(experiment)

        ! {gcm}/{experiment} placeholder substitutions applied to ESM forcing
        ! file paths by the shared varslice reader (via its `subs` argument).
        grid_name = trim(grid%name)

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
        call varslice_init_nml(esm%ts_ref, filename, trim(grp_ts_ref), domain, grid_name, subs=esm_subs,grid=grid)
        if (use_smb) then
            call varslice_init_nml(esm%smb_ref, filename, trim(grp_smb_ref), domain, grid_name, subs=esm_subs,grid=grid)
        else
            call varslice_init_nml(esm%pr_ref, filename, trim(grp_pr_ref), domain, grid_name, subs=esm_subs,grid=grid)
        end if
        call varslice_init_nml(esm%zs_ref, filename, trim(grp_zs_ref), domain, grid_name, subs=esm_subs,grid=grid)
        call varslice_init_nml(esm%to_ref, filename, trim(grp_to_ref), domain, grid_name, subs=esm_subs,grid=grid)
        call varslice_init_nml(esm%so_ref, filename, trim(grp_so_ref), domain, grid_name, subs=esm_subs,grid=grid)

        if (use_var) then
            ! Variability
            ! Transient dependent field
            call varslice_init_nml(esm%ts_var, filename, trim(grp_ts_var), domain, grid_name, subs=esm_subs,grid=grid)
            if (use_smb) then
                call varslice_init_nml(esm%smb_var, filename, trim(grp_smb_var), domain, grid_name, &
                        subs=esm_subs,grid=grid)
            else
                call varslice_init_nml(esm%pr_var, filename, trim(grp_pr_var), domain, grid_name, &
                        subs=esm_subs,grid=grid)
            end if
            call varslice_init_nml(esm%to_var, filename, trim(grp_to_var), domain, grid_name, subs=esm_subs,grid=grid)
            call varslice_init_nml(esm%so_var, filename, trim(grp_so_var), domain, grid_name, subs=esm_subs,grid=grid)
            ! Reference period
            call varslice_init_nml(esm%ts_var_ref, filename, trim(grp_ts_var), domain, grid_name, &
                    subs=esm_subs,grid=grid)
            if (use_smb) then
                call varslice_init_nml(esm%smb_var_ref, filename, trim(grp_smb_var), domain, grid_name, &
                        subs=esm_subs,grid=grid)
            else
                call varslice_init_nml(esm%pr_var_ref, filename, trim(grp_pr_var), domain, grid_name, &
                        subs=esm_subs,grid=grid)
            end if
            call varslice_init_nml(esm%to_var_ref, filename, trim(grp_to_var), domain, grid_name, &
                    subs=esm_subs,grid=grid)
            call varslice_init_nml(esm%so_var_ref, filename, trim(grp_so_var), domain, grid_name, &
                    subs=esm_subs,grid=grid)
        end if

        ! Transient dependent fields
        if (trim(esm%ctrl_run_type) .eq. "transient" .and. trim(esm%experiment) .ne. "ctrl") then
            ! ESM reference period
            if (use_esm) then
                call varslice_init_nml(esm%ts_esm_ref, filename, trim(grp_ts_esm_ref), domain, grid_name, &
                        subs=esm_subs,grid=grid)
                if (use_smb) then
                    call varslice_init_nml(esm%smb_esm_ref, filename, trim(grp_smb_esm_ref), domain, grid_name, &
                            subs=esm_subs,grid=grid)
                else
                    call varslice_init_nml(esm%pr_esm_ref, filename, trim(grp_pr_esm_ref), domain, grid_name, &
                            subs=esm_subs,grid=grid)
                end if
                call varslice_init_nml(esm%to_esm_ref, filename, trim(grp_to_esm_ref), domain, grid_name, &
                        subs=esm_subs,grid=grid)
                call varslice_init_nml(esm%so_esm_ref, filename, trim(grp_so_esm_ref), domain, grid_name, &
                        subs=esm_subs,grid=grid)
                call varslice_init_nml(esm%zs_esm_ref, filename, trim(grp_zs_esm_ref), domain, grid_name, &
                        subs=esm_subs,grid=grid)

                ! ESM historical period
                if (use_hist) then
                    call varslice_init_nml(esm%ts_hist, filename,trim(grp_ts_hist), domain, grid_name, &
                            subs=esm_subs,grid=grid)
                    if (use_smb) then
                        call varslice_init_nml(esm%smb_hist, filename, trim(grp_smb_hist), domain, grid_name, &
                                subs=esm_subs,grid=grid)
                        call varslice_init_nml(esm%dsmbdz_hist, filename, trim(grp_dsmbdz_hist), domain, grid_name, &
                                subs=esm_subs,grid=grid)
                    else
                        call varslice_init_nml(esm%pr_hist, filename, trim(grp_pr_hist), domain, grid_name, &
                                subs=esm_subs,grid=grid)
                    end if
                    call varslice_init_nml(esm%to_hist, filename,trim(grp_to_hist), domain, grid_name, &
                            subs=esm_subs,grid=grid)
                    call varslice_init_nml(esm%so_hist, filename,trim(grp_so_hist), domain, grid_name, &
                            subs=esm_subs,grid=grid)
                    call esm_check_ocean_layout(esm%to_hist, esm%to_esm_ref)
                    call esm_check_ocean_layout(esm%so_hist, esm%so_esm_ref)
                    esm%with_sgd_hist = sgd_active(filename, grp_Qd_hist, domain)
                    if (esm%with_sgd_hist) then
                        call varslice_init_nml(esm%Qd_hist, filename,trim(grp_Qd_hist), domain,grid_name, &
                                subs=esm_subs,grid=grid)
                    end if
                end if
                
                ! ESM projection period
                if (use_proj) then
                    ! atm
                    call varslice_init_nml(esm%ts_proj, filename,trim(grp_ts_proj), domain,grid_name, &
                            subs=esm_subs,grid=grid)
                    if (use_smb) then
                        call varslice_init_nml(esm%smb_proj, filename, trim(grp_smb_proj), domain, grid_name, &
                                subs=esm_subs,grid=grid)
                        call varslice_init_nml(esm%dsmbdz_proj, filename, trim(grp_dsmbdz_proj), domain, grid_name, &
                                subs=esm_subs,grid=grid)
                    else
                        call varslice_init_nml(esm%pr_proj, filename, trim(grp_pr_proj), domain, grid_name, &
                                subs=esm_subs,grid=grid)
                    end if
                    ! ocean
                    call varslice_init_nml(esm%to_proj, filename,trim(grp_to_proj), domain,grid_name, &
                            subs=esm_subs,grid=grid)
                    call varslice_init_nml(esm%so_proj, filename,trim(grp_so_proj), domain,grid_name, &
                            subs=esm_subs,grid=grid)
                    call esm_check_ocean_layout(esm%to_proj, esm%to_esm_ref)
                    call esm_check_ocean_layout(esm%so_proj, esm%so_esm_ref)
                    esm%with_sgd_proj = sgd_active(filename, grp_Qd_proj, domain)
                    if (esm%with_sgd_proj) then
                        call varslice_init_nml(esm%Qd_proj, filename,trim(grp_Qd_proj), domain,grid_name, &
                                subs=esm_subs,grid=grid)
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

    function sgd_active(filename, group, domain) result(active)
        ! Subglacial discharge (Qd) is a Greenland forcing; it is read from the
        ! group's file unless its filename is "none" (Qd = 0).
        character(len=*), intent(IN) :: filename, group, domain
        logical :: active
        character(len=512) :: sgd_file

        active = .false.
        if (trim(domain) .ne. "Greenland") return

        call nml_read(filename, trim(group), "filename", sgd_file)
        active = .not. (trim(sgd_file) .eq. "none" .or. trim(sgd_file) .eq. "None")
        if (.not. active) write(*,*) "esm_forcing_init:: "//trim(group)//": filename = none, Qd = 0."
    end function sgd_active

    subroutine esm_clim_update(esm,z_srf_ylm,time,time_ref,use_smb,south)
        ! Routine to update reference climatology to the specific Antarctic elevation and ocean (neccessary?)

        implicit none

        type(esm_forcing_class), intent(INOUT) :: esm
        real(wp),                intent(IN)    :: z_srf_ylm(:,:)
        real(dp),                intent(IN)    :: time
        real(dp),                intent(IN)    :: time_ref(2)
        logical,                 intent(IN)    :: use_smb
        logical,                 intent(IN)    :: south     ! southern hemisphere

        ! Local variables 
        integer :: m
        real(wp) :: lapse, lapse_mon(12)

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
        call varslice_update(esm%to_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)
        call varslice_update(esm%so_ref, [time_ref(1),time_ref(2)],method="range_mean",rep=1)

        ! Convert atmospheric fields to model elevation
        lapse_mon = esm_lapse_monthly(esm, south)
        do m = 1,12 
            lapse = lapse_mon(m)
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

        return

    end subroutine esm_clim_update

    function esm_lapse_monthly(esm, south) result(lapse)
        ! Monthly lapse rate [K/m] of the reference atmosphere: annual lapse(1)
        ! with a seasonal cosine (phase of the month start), summer lapse(2).
        type(esm_forcing_class), intent(IN) :: esm
        logical,                 intent(IN) :: south
        real(wp) :: lapse(12)
        real(wp), parameter :: pi = 3.14159265359
        integer :: m
        do m = 1, 12
            if (south) then
                lapse(m) = (esm%lapse(1)+(esm%lapse(2)-esm%lapse(1))*cos(2*pi*(m*30.4375-30.4375)/365.25))
            else
                lapse(m) = (esm%lapse(1)+(esm%lapse(1)-esm%lapse(2))*cos(2*pi*(m*30.4375-30.4375)/365.25))
            end if
        end do
    end function esm_lapse_monthly

    subroutine esm_variability_update(esm,mshlf,time,dtt,clim_var,time_ref,H_ice,basins,z_bed,f_grnd,z_sl,use_var)
        ! Update climatic fields. These will be used as bnd conditions for Yelmo.
        ! Output are anomaly fields with respect to a reference field from the ESM.
    
        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        type(marshelf_class),    intent(IN)    :: mshlf
        real(dp), intent(IN) :: time
        real(dp), intent(IN) :: dtt
        character(len=*), intent(IN) :: clim_var
        real(dp), intent(IN) :: time_ref(2)
        real(wp), intent(IN) :: H_ice(:,:),basins(:,:),z_bed(:,:),f_grnd(:,:),z_sl(:,:)
        logical,  intent(IN) :: use_var
    
        ! Local variables 
        integer  :: m, step_idx
        real(wp) :: rand
        real(dp) :: year_rand
        real(wp), allocatable :: to_ref(:,:,:), so_ref(:,:,:)   ! the ocean fields, extrapolated
        real(wp), allocatable :: to_now(:,:,:), so_now(:,:,:)

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
        
        to_ref = esm%to_var_ref%var(:,:,:,1)
        so_ref = esm%so_var_ref%var(:,:,:,1)
        call ocn_variable_extrapolation(to_ref,H_ice,basins,-esm%to_var_ref%z,z_bed)
        call ocn_variable_extrapolation(so_ref,H_ice,basins,-esm%so_var_ref%z,z_bed)

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
                to_now = esm%to_var%var(:,:,:,1)
                so_now = esm%so_var%var(:,:,:,1)
                call ocn_variable_extrapolation(to_now,H_ice,basins,-esm%to_var%z,z_bed)
                call ocn_variable_extrapolation(so_now,H_ice,basins,-esm%so_var%z,z_bed)

                ! Compute the anomaly at the desired depth level
                call marshelf_interp_shelf(esm%dto_var,mshlf,to_now-to_ref,H_ice, &
                                            z_bed,f_grnd,z_sl,-esm%to_var_ref%z)
                call marshelf_interp_shelf(esm%dso_var,mshlf,so_now-so_ref,H_ice, &
                                            z_bed,f_grnd,z_sl,-esm%so_var_ref%z)

            case DEFAULT
                ! Assume no variabiliy anomaly (do nothing)
                
        end select
        end if

        return 
    
    end subroutine esm_variability_update

    subroutine esm_forcing_update(esm,mshlf,time,use_esm,time_ref,time_hist,time_proj,time_esm_ref,&
                                  H_ice,basins,z_bed,f_grnd,z_sl,use_smb)
        ! Update climatic fields. These will be used as bnd conditions for Yelmo.
        ! Output are anomaly fields with respect to a reference field from the ESM.
    
        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        type(marshelf_class),    intent(IN)    :: mshlf
        real(dp), intent(IN) :: time
        logical,  intent(IN) :: use_esm
        real(dp), intent(IN) :: time_ref(2),time_hist(2),time_proj(2),time_esm_ref(2)
        real(wp), intent(IN) :: H_ice(:,:),basins(:,:),z_bed(:,:),f_grnd(:,:),z_sl(:,:)
        logical,  intent(IN) :: use_smb  
    
        ! Local variables
        integer  :: m
        real(wp) :: anomaly
        real(wp), allocatable :: to_esm_ref(:,:,:), so_esm_ref(:,:,:)   ! the esm reference ocean, extrapolated

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
                if (trim(esm%experiment) == "ctrl") then
                    ! The ctrl experiment: the reference climatology only, as
                    ! esm_forcing_init loads no ESM fields for it (use_esm or not).

                else if (use_esm) then ! Fields loaded from ESM
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
                    to_esm_ref = esm%to_esm_ref%var(:,:,:,1)
                    so_esm_ref = esm%so_esm_ref%var(:,:,:,1)
                    if (mshlf%par%extrap_shlf) then
                        call ocn_variable_extrapolation(to_esm_ref,H_ice,basins,-esm%to_esm_ref%z,z_bed)
                        call ocn_variable_extrapolation(so_esm_ref,H_ice,basins,-esm%so_esm_ref%z,z_bed)
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
                                esm%dsmb(:,:,m) = esm%smb_hist%var(:,:,m,1)-esm%smb_esm_ref%var(:,:,m,1)
                            else
                                esm%dpr(:,:,m) = esm%pr_hist%var(:,:,m,1)/(esm%pr_esm_ref%var(:,:,m,1)+1e-8)
                            end if
                        end do
                        ! ===   Oceanic fields   ===
                        call esm_ocean_anomaly(esm%dto,mshlf,esm%to_hist,to_esm_ref,time, &
                                               H_ice,basins,z_bed,f_grnd,z_sl)
                        call esm_ocean_anomaly(esm%dso,mshlf,esm%so_hist,so_esm_ref,time, &
                                               H_ice,basins,z_bed,f_grnd,z_sl)
                        if (esm%with_sgd_hist) then
                            call varslice_update(esm%Qd_hist,[time],method="extrap",rep=12)
                            esm%Qd_ann = sum(esm%Qd_hist%var(:,:,:,1),dim=3) / 12.0
                            esm%Qd_sum = (esm%Qd_hist%var(:,:,6,1)+esm%Qd_hist%var(:,:,7,1)+esm%Qd_hist%var(:,:,8,1)) / 3.0
                        end if

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
                        call esm_ocean_anomaly(esm%dto,mshlf,esm%to_proj,to_esm_ref,time, &
                                               H_ice,basins,z_bed,f_grnd,z_sl)
                        call esm_ocean_anomaly(esm%dso,mshlf,esm%so_proj,so_esm_ref,time, &
                                               H_ice,basins,z_bed,f_grnd,z_sl)
                        if (esm%with_sgd_proj) then
                            call varslice_update(esm%Qd_proj,[time],method="extrap",rep=12)
                            esm%Qd_ann = sum(esm%Qd_proj%var(:,:,:,1),dim=3) / 12.0
                            esm%Qd_sum = (esm%Qd_proj%var(:,:,6,1)+esm%Qd_proj%var(:,:,7,1)+esm%Qd_proj%var(:,:,8,1)) / 3.0
                        end if

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

                        case DEFAULT
                            write(*,*) "esm_forcing_update:: Error: transient experiment not recognized: "//trim(esm%experiment)
                            stop

                    end select

                end if        
    
            case DEFAULT
                write(*,*) "esm_forcing_update:: Error: ctrl_run_type not recognized: "//trim(esm%ctrl_run_type)
                stop 
     
        end select
    
        return 
    
    end subroutine esm_forcing_update

    subroutine esm_ocean_anomaly(dvar,mshlf,vs,var_ref,time,H_ice,basins,z_bed,f_grnd,z_sl)
        ! Shelf-depth anomaly of a transient ESM ocean field (historical or
        ! projection) relative to its ESM reference (var_ref, on the same levels:
        ! esm_check_ocean_layout). The field is sliced at the
        ! current time with its declared sub-annual period (varslice_nsub) and
        ! averaged to an annual (x,y,z) field (a single level for depth-less
        ! forcing, e.g. Greenland ISMIP7), optionally extrapolated into the
        ! ice-shelf interiors, and the anomaly is interpolated to the shelf depth.

        implicit none

        real(wp),             intent(INOUT) :: dvar(:,:)
        type(marshelf_class), intent(IN)    :: mshlf
        type(varslice_class), intent(INOUT) :: vs
        real(wp),             intent(IN)    :: var_ref(:,:,:)
        real(dp),             intent(IN)    :: time
        real(wp),             intent(IN)    :: H_ice(:,:),basins(:,:),z_bed(:,:),f_grnd(:,:),z_sl(:,:)

        ! Local variables
        real(wp), allocatable :: var_ann(:,:,:)

        call varslice_update(vs,[time],method="extrap",rep=varslice_nsub(vs))
        var_ann = varslice_sub_mean(vs)

        if (mshlf%par%extrap_shlf) then
            call ocn_variable_extrapolation(var_ann,H_ice,basins,-vs%z,z_bed)
        end if

        call marshelf_interp_shelf(dvar,mshlf,var_ann-var_ref,H_ice, &
                                    z_bed,f_grnd,z_sl,-vs%z)

        return

    end subroutine esm_ocean_anomaly

    subroutine esm_check_ocean_layout(vs,vs_ref)
        ! The transient ocean anomaly is taken level by level against the ESM
        ! reference (esm_ocean_anomaly), so both fields must have the same
        ! spatial layout: rank, horizontal extent and depth levels.

        implicit none

        type(varslice_class), intent(IN) :: vs
        type(varslice_class), intent(IN) :: vs_ref

        ! Local variables
        logical :: same

        same = (vs%par%ndim .eq. vs_ref%par%ndim)              .and. &
               (vs%par%with_time .eqv. vs_ref%par%with_time)   .and. &
               size(vs%x) .eq. size(vs_ref%x) .and. size(vs%y) .eq. size(vs_ref%y) .and. &
               size(vs%z) .eq. size(vs_ref%z)
        if (same) same = all(abs(vs%z-vs_ref%z) .le. 1e-3_wp)

        if (.not. same) then
            write(error_unit,*) "esm_check_ocean_layout:: Error: ocean field layout differs from &
                                &its ESM reference."
            write(error_unit,*) "field:     ", trim(vs%par%name), " (", trim(vs%par%filename), ")"
            write(error_unit,*) "  ndim, nx, ny, nz = ", vs%par%ndim, size(vs%x), size(vs%y), size(vs%z)
            write(error_unit,*) "reference: ", trim(vs_ref%par%name), " (", trim(vs_ref%par%filename), ")"
            write(error_unit,*) "  ndim, nx, ny, nz = ", vs_ref%par%ndim, size(vs_ref%x), size(vs_ref%y), size(vs_ref%z)
            error stop
        end if

        return

    end subroutine esm_check_ocean_layout

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
            if (allocated(esm%t2m_ann)) deallocate(esm%t2m_ann)
            if (allocated(esm%pr_ann))  deallocate(esm%pr_ann)
            if (allocated(esm%smb_ann)) deallocate(esm%smb_ann)

            return
    
        end subroutine esm_deallocate
    
end module esm_forcing

