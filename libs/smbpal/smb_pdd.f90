module smb_pdd
    ! This module contains all subroutines related to the calculation of
    ! surface mass balance using the PDD method.

    use smbpal_precision 

    implicit none

    real(prec), parameter :: sec_day = 86400.0   ! [sec]
    real(prec), parameter :: rho_w   = 1.d3      ! Density of pure water [kg/m3]
    real(prec), parameter :: L_m     = 3.35e5    ! Latent heat of melting [J/kg]

contains 


    elemental subroutine calc_ablation_pdd(melt,runoff,refrz,pdds,acc,mm_snow,mm_ice,f_refrz_max)
        ! Determine total annual ablation based on input pdds
        ! All [mm] quantitites in [mm water equiv.]

        implicit none 

        real(prec), intent(OUT) :: melt                 ! [mm] Total melt rate from energy 
        real(prec), intent(OUT) :: runoff               ! [mm] Actual net melt rate (runoff) based on snowpack 
        real(prec), intent(OUT) :: refrz                ! [mm] Refrozen melt as superimposed ice 
        real(prec), intent(IN)  :: pdds                 ! [K-d] Positive degree days 
        real(prec), intent(IN)  :: acc                  ! [mm] Accumulation 
        real(prec), intent(IN)  :: mm_snow              ! [mm w.e./K-d] Conversion factor PDDs to mm snow melt
        real(prec), intent(IN)  :: mm_ice               ! [mm w.e./K-d] Conversion factor PDDs to mm ice melt 
        real(prec), intent(IN)  :: f_refrz_max          ! [-] Maximum refreezing fraction 

        ! Local variables 
        real(prec) :: refrz_max
        real(prec) :: melt_snow
        real(prec) :: melt_ice 

        ! Define maximum amount of superimposed ice from refreezing that can be formed [mm]
        refrz_max = acc*f_refrz_max
       
        ! Determine the potential snowmelt (mm) from the 
        ! available pdds, then determine the potential ice melt (mm)
        ! from the left over pdds, if there are any.
        melt_snow = pdds*mm_snow 
        melt_ice  = max(0.0, (melt_snow-acc)*mm_ice/mm_snow)
        melt      = melt_snow + melt_ice 
        
        ! Calculate the actual superimposed ice from refreezing 
        refrz = min(melt_snow,refrz_max)

        ! Get the runoff / net melt [mm] - must be positive or zero 
        runoff = melt - refrz 

        return 

    end subroutine calc_ablation_pdd

    subroutine calc_temp_effective(teff, temp, sigma)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Subroutine : e f f e c t i v e T
        ! Author     : Reinhard Calov
        ! Purpose    : Computation of the positive degree days (PDD) with
        !              statistical temperature fluctuations;
        !              based on semi-analytical solution by Reinhard Calov.
        !              This subroutine uses days as time unit, each day
        !              is added individually
        !              (the same sigma as for pdd monthly can be used)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! temp = [degrees Celcius] !!! 

        implicit none

        real(prec), intent(OUT) :: teff(:,:)
        real(prec), intent(IN)  :: temp(:,:)
        real(prec), intent(IN)  :: sigma(:,:)
        
        real(prec), parameter :: inv_sqrt2   = 1.0/sqrt(2.0)
        real(prec), parameter :: inv_sqrt2pi = 1.0/sqrt(2.0*pi)

        real(prec) :: inv_sigma
        real(8) :: val_in
        real(8) :: val_out 
        real(8) :: val_tmp 

        integer :: i, j, nx, ny 

        nx = size(teff,1) 
        ny = size(teff,2) 

        do j = 1, ny 
        do i = 1, nx 

            inv_sigma   = 1.0/sigma(i,j)
            
            ! Perform calculation of erfc at real(8) precision to
            ! avoid potential underflow errors at very low temperatures
            val_in  = -temp(i,j)*inv_sigma*inv_sqrt2
            val_out = erfc(val_in)
            if (dabs(val_out) .lt. 1e-10_dp) val_out = 0.0_dp 

            val_tmp = dexp(-0.5_dp*(temp(i,j)*inv_sigma)**2)
            if (dabs(val_tmp) .lt. 1e-10_dp) val_tmp = 0.0_dp 

            teff(i,j) = sigma(i,j)*inv_sqrt2pi*val_tmp + temp(i,j)*0.5*val_out 

        end do 
        end do 

        ! Result is the assumed/felt/effective positive degrees, 
        ! given the actual temperature (accounting for fluctuations in day/month/etc, 
        ! based on the sigma chosen)
                     
        return

    end subroutine calc_temp_effective

    elemental function calc_temp_effective_old(temp, sigma) result(teff)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Subroutine : e f f e c t i v e T
        ! Author     : Reinhard Calov
        ! Purpose    : Computation of the positive degree days (PDD) with
        !              statistical temperature fluctuations;
        !              based on semi-analytical solution by Reinhard Calov.
        !              This subroutine uses days as time unit, each day
        !              is added individually
        !              (the same sigma as for pdd monthly can be used)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! temp = [degrees Celcius] !!! 

        ! ajr: Note: this version can apparently produce an arithmetic
        ! error due to an underflow when a very negative temperature is
        ! fed into the erfc equation, then -temp*inv_sigma*inv_sqrt2 >> 0,
        ! and the output is a very tiny number, causing a crash. 

        implicit none

        real(prec), intent(IN) :: temp, sigma
        real(prec) :: teff
        
        real(prec) :: temp_c, inv_sigma
        real(prec), parameter :: inv_sqrt2   = 1.0/sqrt(2.0)
        real(prec), parameter :: inv_sqrt2pi = 1.0/sqrt(2.0*pi)

        inv_sigma   = 1.0/sigma

        teff = sigma*inv_sqrt2pi*exp(-0.5*(temp*inv_sigma)**2)  &
                  + temp*0.5*erfc(-temp*inv_sigma*inv_sqrt2)

        ! Result is the assumed/felt/effective positive degrees, 
        ! given the actual temperature (accounting for fluctuations in day/month/etc, 
        ! based on the sigma chosen)
                     
        return

    end function calc_temp_effective_old

end module smb_pdd
