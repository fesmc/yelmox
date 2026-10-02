module yelmox_esm_output
    ! Output writers for the multigrid ESM driver (yelmox_esm).
    !
    ! The two 2D writers take the ESM/SMB switch as an explicit use_smb logical.
    !
    !   write_step_2D_combined   standard heavy 2D output
    !   write_step_2D_small      small 2D output
    !   write_1D_esm             ESM 1D timeseries

    use nml
    use ncio
    use yelmo
    use esm_forcing
    use fastisostasy    ! isos_class (reexports barysealevel)
    use marine_shelf
    use smbpal

    implicit none

    private
    public :: write_step_2D_combined
    public :: write_step_2D_small
    public :: write_1D_esm

contains

    subroutine write_step_2D_combined(ylmo,isos,esm,mshlf,srf,use_smb,filename,time)

        implicit none

        type(yelmo_class),       intent(IN) :: ylmo
        type(isos_class),        intent(IN) :: isos
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        type(smbpal_class),      intent(IN) :: srf
        logical,                 intent(IN) :: use_smb

        character(len=*),       intent(IN) :: filename
        real(wp),               intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        logical  :: south

        south = (trim(ylmo%par%domain) .eq. "Antarctica")

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Note: numerics/speed metrics now go to yelmo_metrics.nc (write_metrics);
        ! they are no longer embedded here.

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
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb_flt+ylmo%tpo%now%cmb_grnd,units="m/a ice equiv.",long_name="Calving mass balance rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/a ice equiv.",long_name="Calving mass balance rate flt", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_grnd",ylmo%tpo%now%cmb_grnd,units="m/a ice equiv.",long_name="Calving mass balance rate grnd", &
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
        call nc_write(filename,"hyd_W_til",ylmo%hyd%now%W_til,units="m",long_name="Basal water layer thickness", &
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
        call nc_write(filename,"t2m_sum",esm%t2m_sum+esm_summer_mean(esm%dts,south),units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts_ann",SUM(esm%dts, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        if (use_smb) then
            call nc_write(filename,"dsmb_ann",1e-3*SUM(esm%dsmb, dim=3)/12.0,units="m/a water equiv.",long_name="SMB anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        else
            call nc_write(filename,"pr_ann",SUM(esm%pr*esm%dpr, dim=3)/12.0,units="mm/d water equiv.",long_name="Precipitation (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dpr_ann",SUM(esm%dpr, dim=3)/12.0,units="%",long_name="Precipitation anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            !call nc_write(filename,"dpr_var",SUM(esm%dpr_var, dim=3)/12.0,units="%",long_name="Precipitation anomaly (variability)", &
            !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        end if
        
        ! Oceanic boundary conditions
        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto",esm%dto,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso",esm%dso,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"dto_var",esm%dto_var,units="K",long_name="Shelf temperature anomaly (variability)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"dso_var",esm%dso_var,units="PSU",long_name="Shelf salinity anomaly (variability)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

    subroutine write_step_2D_small(ylmo,isos,esm,mshlf,smbp,use_smb,filename,time)

        implicit none

        type(yelmo_class),       intent(IN) :: ylmo
        type(isos_class),        intent(IN) :: isos
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        type(smbpal_class),      intent(IN) :: smbp
        logical,                 intent(IN) :: use_smb

        character(len=*),        intent(IN) :: filename
        real(wp),                intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        logical  :: south

        south = (trim(ylmo%par%domain) .eq. "Antarctica")

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Note: numerics/speed metrics now go to yelmo_metrics.nc (write_metrics);
        ! they are no longer embedded here.

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
        call nc_write(filename,"cb_tgt",ylmo%dyn%now%cb_tgt,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! === yelmo forcing ===
        ! ESM Atmospheric boundary fields            
        call nc_write(filename,"t2m_ann",esm%t2m_ann+SUM(esm%dts, dim=3)/12.0,units="K",long_name="Near-surface air temperature (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"t2m_sum",esm%t2m_sum+esm_summer_mean(esm%dts,south),units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts",SUM(esm%dts, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts_var",SUM(esm%dts_var, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_ann",ylmo%tpo%now%smb,units="m/a water equiv.",long_name="SMB (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        if (use_smb) then
            call nc_write(filename,"dsmb_ann",1e-3*SUM(esm%dsmb, dim=3)/12.0,units="m/a water equiv.",long_name="SMB anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dsmbdz",1e-3*esm%dsmbdz,units="m/a m-1 water equiv.",long_name="SMB lapse rate", &
                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        else
            call nc_write(filename,"pr_ann",SUM(esm%pr*esm%dpr, dim=3)/12.0,units="mm/d water equiv.",long_name="Precipitation (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dpr_ann",SUM(esm%dpr, dim=3)/12.0,units="%",long_name="Precipitation anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dpr_var",SUM(esm%dpr_var, dim=3)/12.0,units="%",long_name="Precipitation anomaly (variability)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        end if

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
        call nc_write(filename,"Qd_ann",esm%Qd_ann,units="m3/s",long_name="Subglacial discharge (annual)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Qd_sum",esm%Qd_sum,units="m3/s",long_name="Sunglacial discharge (summer)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_small

    ! ===== esm output routines =========
    subroutine write_1D_esm(ylmo, esm, mshlf, filename, time)

        ! Used to plot climatic variable fields
    
        implicit none
    
        type(yelmo_class),       intent(IN) :: ylmo
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        character(len=*),        intent(IN) :: filename
        real(wp),                intent(IN) :: time
    
        ! Local variables
        type(yregions_class) :: reg
    
        integer  :: ncid, n
        real(wp) :: rho_ice, esm_correction
    
        real(wp) :: dx, dy
        integer  :: npts_tot, npts_flt
        real(wp) :: smb_tot, bmb_shlf_t
    
        ! Climatic variables - atmosphere
        real(wp) :: t2m_1d, pr_1d
        real(wp) :: dt_1d,  dt_var_1d, dpr_1d, dpr_var_1d
    
        ! Climatic variables - ocean
        real(wp) :: to_1d
        real(wp) :: so_1d
        real(wp) :: tf_1d
        real(wp) :: dto_1d, dto_var_1d
        real(wp) :: dso_1d, dso_var_1d

        ! Missing-value (N/A) sentinel for the ocean shelf-draft diagnostics.
        ! The marine-shelf module leaves a large out-of-range fill value in cells
        ! where no ocean data is available (e.g. a minority of Greenland floating
        ! cells, whose tf forcing carries a NaN _FillValue). ocn_mean below drops
        ! those PER FIELD (fill is not identical across T_shlf/S_shlf/tf_shlf) via
        ! a physical bound, returning this sentinel only if a field has no valid
        ! floating-ice cell.
        real(wp), parameter :: mv_ocn = -9999.0_wp

        logical, allocatable :: mask_tot(:,:)
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:)
    
        dx = ylmo%grd%G%dx
        dy = ylmo%grd%G%dy
    
        allocate(mask_tot (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_grnd(ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_flt (ylmo%grd%G%nx, ylmo%grd%G%ny))
    
        ! === Unit conversion factors =========================================
        ! Taken from Yelmo's own constants, so the output cannot disagree with
        ! the model that produced it. esm_correction converts a rate in
        ! [<unit> yr-1] of ice to [kg ... s-1]; it collapses to rho_ice/sec_year,
        ! since the old two-step form was (rho_w/sec_year)*(rho_ice/rho_w) with
        ! rho_w/sec_year rounded to 3.2e-5 (1 % high) and sec_year written as
        ! 31556952 (the Gregorian year) rather than Yelmo's own.
        rho_ice        = ylmo%bnd%c%rho_ice
        esm_correction = rho_ice / ylmo%bnd%c%sec_year
    
        ! === Masks ===========================================================
    
        mask_tot  = (ylmo%tpo%now%H_ice .gt. 0.0_wp)
        mask_grnd = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp)
        mask_flt  = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
    
        npts_tot = count(mask_tot)
        npts_flt = count(mask_flt)
    
        ! === Regional object =================================================
    
        reg = ylmo%reg
    
        ! === Integrated fluxes [m yr-1 * m2 -> m3 yr-1] =====================
        ! Total SMB over all ice-covered cells [m3 yr-1]
        smb_tot    = sum(ylmo%bnd%smb,       mask=mask_tot)  * (dx * dy)
    
        ! Total BMB beneath floating ice [m3 yr-1]
        bmb_shlf_t = sum(ylmo%bnd%bmb_shlf, mask=mask_flt)  * (dx * dy)
    
        ! === Spatially averaged climatic fields ==============================
    
        ! Atmosphere (averaged over all ice)
        if (npts_tot .gt. 0.0) then
            t2m_1d     = sum(esm%t2m_ann + sum(esm%dts, dim=3)/12.0_wp,  mask=mask_tot) / npts_tot
            pr_1d      = sum(sum(esm%pr*esm%dpr, dim=3)/12.0_wp,          mask=mask_tot) / npts_tot
            dt_1d      = sum(sum(esm%dts, dim=3)/12.0_wp,                 mask=mask_tot) / npts_tot
            dpr_1d     = sum(100.0_wp * sum(esm%dpr, dim=3)/12.0_wp,      mask=mask_tot) / npts_tot
            dt_var_1d  = sum(sum(esm%dts_var, dim=3)/12.0_wp,             mask=mask_tot) / npts_tot
            dpr_var_1d = sum(100.0_wp * sum(esm%dpr_var, dim=3)/12.0_wp,  mask=mask_tot) / npts_tot
        else
            t2m_1d = 0.0_wp; pr_1d = 1.0_wp; dt_1d = 0.0_wp; dpr_1d = 1.0_wp
            dt_var_1d = 0.0_wp; dpr_var_1d = 0.0_wp
        end if

        ! Ocean (averaged over floating ice, per field, excluding fill cells).
        ! Each shelf field is averaged only over floating cells whose value is
        ! finite and inside a physical bound, so marine-shelf fill cells (large
        ! out-of-range sentinel) no longer contaminate the mean. A field returns
        ! mv_ocn (N/A) if it has no valid floating cell -- e.g. Greenland tf_1d,
        ! whose tf_shlf (~273 K here, an artifact of feeding tf-as-temperature)
        ! falls outside the plausible thermal-forcing range; the Greenland shelf
        ! temperature/salinity (to_1d/so_1d) remain valid.
        to_1d      = ocn_mean(mshlf%now%T_shlf,  mask_flt,  240.0_wp, 320.0_wp, mv_ocn)
        so_1d      = ocn_mean(mshlf%now%S_shlf,  mask_flt,    0.0_wp,  60.0_wp, mv_ocn)
        tf_1d      = ocn_mean(mshlf%now%tf_shlf, mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dto_1d     = ocn_mean(esm%dto,           mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dso_1d     = ocn_mean(esm%dso,           mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dto_var_1d = ocn_mean(esm%dto_var,       mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dso_var_1d = ocn_mean(esm%dso_var,       mask_flt, -100.0_wp, 100.0_wp, mv_ocn)

        ! === Write to file ===================================================
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
            
        ! -- Variability fields -----------------------------------------------
        call nc_write(filename, "dt_var_1d",  dt_var_1d,  units="K",   &
            long_name="Mean ice surf. Temp. Anomaly (Variability)",     &
            standard_name="Mean ice surf. Temp. Anomaly (Variability)", &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dpr_var_1d", dpr_var_1d, units="%",   &
            long_name="Mean ice surf. Pr. Anomaly (Variability)",       &
            standard_name="Mean ice surf. Pr. Anomaly (Variability)",   &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dto_var_1d", dto_var_1d, units="K",   &
            long_name="Mean ice-shelf Temp. Anomaly (Variability)",     &
            standard_name="Mean ice-shelf Temp. Anomaly (Variability)", &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dso_var_1d", dso_var_1d, units="PSU", &
            long_name="Mean ice-shelf Sal. Anomaly (Variability)",      &
            standard_name="Mean ice-shelf Sal. Anomaly (Variability)",  &
            dim1="time", start=[n], ncid=ncid)
    
        ! -- Atmosphere fields ------------------------------------------------
        call nc_write(filename, "t2m_1d",  t2m_1d,  units="K",       &
            long_name="Mean ice surf. Temp.",                          &
            standard_name="Mean ice surf. Temp.",                      &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "pr_1d",   pr_1d,   units="mm d-1",  &
            long_name="Mean ice surf. Pr.",                            &
            standard_name="Mean ice surf. Pr.",                        &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dt_1d",   dt_1d,   units="K",       &
            long_name="Mean ice surf. Temp. Anomaly",                  &
            standard_name="Mean ice surf. Temp. Anomaly",              &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dpr_1d",  dpr_1d,  units="%",       &
            long_name="Mean ice surf. Pr. Anomaly",                    &
            standard_name="Mean ice surf. Pr. Anomaly",                &
            dim1="time", start=[n], ncid=ncid)
    
        ! -- Ocean fields -----------------------------------------------------
        call nc_write(filename, "to_1d",   to_1d,   units="K",       &
            long_name="Mean ice-shelf Temp.",                          &
            standard_name="Mean ice-shelf draft Temp.",                &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "so_1d",   so_1d,   units="PSU",     &
            long_name="Mean ice-shelf draft Sal.",                     &
            standard_name="Mean ice-shelf draft Sal.",                 &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "tf_1d",   tf_1d,   units="K",       &
            long_name="Mean ice-shelf TF.",                            &
            standard_name="Mean ice-shelf draft TF",                   &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dto_1d",  dto_1d,  units="K",       &
            long_name="Mean ice-shelf Temp. Anomaly",                  &
            standard_name="Mean ice-shelf Temp. Anomaly",              &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dso_1d",  dso_1d,  units="PSU",     &
            long_name="Mean ice-shelf draft Sal. Anomaly",             &
            standard_name="Mean ice-shelf draft Sal. Anomaly",         &
            dim1="time", start=[n], ncid=ncid)
    
        ! -- Integrated mass fluxes [kg s-1] ----------------------------------
        call nc_write(filename, "smb_tot",  smb_tot  * esm_correction, units="kg s-1", &
            long_name="Total SMB flux",                                                 &
            standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance",     &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "bmb_shlf", bmb_shlf_t * esm_correction, units="kg s-1", &
            long_name="Total BMB flux beneath floating ice",                              &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",         &
            dim1="time", start=[n], ncid=ncid)
    
        call nc_close(ncid)
    
        return
    
    end subroutine write_1D_esm

    function ocn_mean(field, base_mask, lo, hi, mv) result(val)
        ! Mean of `field` over the `base_mask` cells whose value is finite and
        ! inside the physical range (lo,hi); returns `mv` (missing value / N/A)
        ! if no cell qualifies. NaN/Inf fail the comparisons and so are dropped.
        ! Keeps marine-shelf fill cells (large out-of-range sentinel) out of the
        ! floating-ice ocean-forcing diagnostics, per field.
        implicit none
        real(wp), intent(IN) :: field(:,:)
        logical,  intent(IN) :: base_mask(:,:)
        real(wp), intent(IN) :: lo, hi, mv
        real(wp) :: val
        ! Local variables
        logical, allocatable :: m(:,:)
        integer :: npts

        allocate(m(size(field,1), size(field,2)))
        m = base_mask .and. field .gt. lo .and. field .lt. hi
        npts = count(m)
        if (npts .gt. 0) then
            val = sum(field, mask=m) / real(npts, wp)
        else
            val = mv
        end if

        return
    end function ocn_mean


end module yelmox_esm_output
