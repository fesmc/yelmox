module yelmox_rembo_output
    ! REMBO-flavor diagnostic output: the combined yelmo+rembo 2D/1D-2D files that
    ! the REMBO Greenland analysis expects. Ported from the legacy driver's
    ! internal write routines; the module takes each sub-model object directly
    ! (yelmo / rembo_ann / isos / mshlf / hyst), so it does not depend on ice_domain.

    use ncio
    use yelmo
    use rembo_sclimate, only : rembo_class
    use fastisostasy,   only : isos_class
    use marine_shelf,   only : marshelf_class
    use hyster,         only : hyster_class

    implicit none
    private

    public :: yelmox_rembo_write_init
    public :: yelmox_rembo_write_step_small
    public :: yelmox_rembo_write_step

contains

    subroutine yelmox_rembo_write_init(dom,filename,time_init,units,mask,dT_min,dT_max)

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
        call nc_write_dim(filename,"xc",   x=real(dom%grd%G%x,wp)*1e-3, units="kilometers")
        call nc_write_dim(filename,"yc",   x=real(dom%grd%G%y,wp)*1e-3, units="kilometers")
        call nc_write_dim(filename,"zeta", x=dom%par%zeta_aa,           units="1")
        call nc_write_dim(filename,"time", x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

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
        call nc_write(filename,"mask", mask, units="1",long_name="Region mask",dim1="xc",dim2="yc")

        return

    end subroutine yelmox_rembo_write_init

    subroutine yelmox_rembo_write_step_small(ylmo,hyst,rembo,isos,mshlf,filename,time, &
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
        logical  :: write_ocn
        type(yregions_class) :: reg

        write_ocn = .FALSE.
        if (present(write_ocn_forcing)) write_ocn = write_ocn_forcing

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
        call yelmo_write_var(filename,"hyd_W_til",ylmo,n,ncid)

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
        if (write_ocn) then

            call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if

        ! == ice-sheet wide metrics ==

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

    end subroutine yelmox_rembo_write_step_small

    subroutine yelmox_rembo_write_step(ylmo,rembo,isos,mshlf,filename,time)

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
        call yelmo_write_var(filename,"visc_int",ylmo,n,ncid)

        ! == yelmo_thermodynamics ==
        call yelmo_write_var(filename,"T_prime",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_pmp",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"hyd_W_til",ylmo,n,ncid)

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

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine yelmox_rembo_write_step

end module yelmox_rembo_output
