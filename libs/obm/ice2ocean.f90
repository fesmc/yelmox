module ice2ocean

    use ncio 

    implicit none
    
    private
    public :: calc_fwf

    integer,  parameter :: wp = kind(1.0)


contains

function calc_fwf(rho_water,rho_ice,sec_year,mb,smb,bmb,cmb,Hice,dHidt,f_grnd,dx,dy,mask,hemisphere, fwf_def) result(fwf)
    ! This subroutine is based on yelmo_regions.f90 calc_yregions(...)

    implicit none 

    real(wp) :: fwf

    ! Local variables
    integer  :: nx, ny
    real(wp) :: npts_flt, npts_tot, smb_tot, bmb_tot, cmb_tot, npts_tot_basin

    real(wp) :: m3_km3 = 1e-9 
    real(wp) :: conv_km3a_Sv

    real(wp)   :: rho_water, rho_ice ! bnd%c%rho_w, bnd%c%rho_ice
    real(wp)   :: sec_year ! bnd%c%sec_year
    real(wp), allocatable   :: mb(:,:), smb(:,:), bmb(:,:), cmb(:,:), Hice(:,:), dHidt(:,:) ! tpo%now%mb, tpo%now%dHidt
    real(wp), allocatable   :: f_grnd(:,:)
    real(wp)           :: dx, dy ! tpo%par%dx, tpo%par%dy
    real(wp), allocatable :: mask(:,:) 
    character(len=*) :: hemisphere

    logical, allocatable :: mask_tot(:,:), mask_basin(:,:), mask_flt(:,:), basin_mask(:,:), extended_north_basin_mask(:,:)
    real(wp), allocatable :: H_af(:,:) 
    logical, parameter :: use_extended_north=.TRUE.

    character(len=*) :: fwf_def

    ! Conversion parameter 
    conv_km3a_Sv = 1e-6*(1e9*rho_water/rho_ice)/sec_year   

    ! Grid size 
    nx = size(dHidt,1)
    ny = size(dHidt,2)

    allocate(mask_tot(nx,ny))
    allocate(mask_flt(nx,ny))
    allocate(H_af(nx,ny)) 
    allocate(extended_north_basin_mask(nx,ny))

    ! Define masks
    mask_tot  = (Hice .gt. 0.0) 
    mask_basin = (mask .eq. 1.0) 
    mask_flt  = (dHidt .gt. 0.0 .and. f_grnd .eq. 0.0)
    npts_tot  = real(count(mask_tot),wp)
    npts_tot_basin = real(count(mask_basin),wp)
    npts_flt  = real(count(mask_flt),wp)

    fwf_def = "dVdt" !"dVdt_mask"

    ! ===== Compute fwf =====
    select case(fwf_def)
        case("mb") 
            fwf = -1*sum(mb,mask=mask_tot)*dx*dy/1e6/sec_year 
        case("ablation")
            where (smb .gt. 0.0) smb = 0.0
            where (bmb .gt. 0.0) bmb = 0.0
            where (cmb .gt. 0.0) cmb = 0.0
            smb_tot = sum(smb,mask=mask_tot)  ! only negative values
            bmb_tot = sum(bmb,mask=mask_tot)
            cmb_tot = sum(cmb,mask=mask_tot)
            fwf = -1*(smb_tot+bmb_tot+cmb_tot)*dx*dy/1e6/sec_year
        case("dVdt")
            !fwf = max(-sum(dHidt, mask=mask_tot)*dx*dy*m3_km3*conv_km3a_Sv, 0.0)
            fwf = -sum(dHidt)*dx*dy*m3_km3*conv_km3a_Sv
        case("dVdt_mask")
            ! fwf = max(-sum(dHidt, mask=mask_basin)*dx*dy*m3_km3*conv_km3a_Sv, 0.0)
            if (npts_tot_basin .gt. 0) then 
                fwf = -sum(dHidt, mask=mask_basin)*dx*dy*m3_km3*conv_km3a_Sv
            else
                fwf = 0.0
            end if
        case DEFAULT
            fwf = max(-sum(dHidt)*dx*dy*m3_km3*conv_km3a_Sv, 0.0)
        
        end select
    return

end function calc_fwf

end module ice2ocean
