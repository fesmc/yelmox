module kryos_regions
    ! Region-specific setup and physics of a kryos_domain: the named regions
    ! for 1D output, the LGM-like marine-ice initial state, the NEGIS
    ! basal-friction modification and the glacial SMB scaling.

    use precision,    only : dp
    use yelmo,        only : yelmo_class, wp
    use yelmo_defs,   only : MASK_ICE_NONE
    use basal_dragging, only : calc_cb_ref
    use kryos,        only : kryos_domain, negis_params, glacial_smb_params, remap

    implicit none
    private

    public :: domain_regions_init, domain_init_marine_ice
    public :: negis_update_cb_ref, calc_glacial_smb

contains

    subroutine domain_regions_init(dom, outfldr)
        ! Output of the domain's regions: Yelmo defines its named regions at init
        ! from the named masks of the regions ([regions] masks, mask_<name>);
        ! here their regional files are placed in outfldr. Without named regions
        ! only the global region is written.
        ! Must be called after domain_init and before the first yelmo_update.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        integer :: i

        ! Hand yelmo its output folder: the single source for all files yelmo
        ! writes internally (regional 1D files below, and yelmo_metrics.nc).
        dom%yelmo%outfldr = trim(outfldr)

        ! Name the regional 1D files (no grid suffix; grid is recorded in-file):
        !   global -> yelmo_ts.nc, sub-region k -> yelmo_ts_<name>.nc
        ! Paths derive from dom%yelmo%outfldr (set above).
        dom%yelmo%reg%fnm = trim(dom%yelmo%outfldr)//"yelmo_ts.nc"
        if (dom%yelmo%par%n_reg > 0) then
            do i = 1, dom%yelmo%par%n_reg
                dom%yelmo%regs(i)%fnm = trim(dom%yelmo%outfldr)//"yelmo_ts_"// &
                                        trim(dom%yelmo%regs(i)%name)//".nc"
            end do
        end if

    end subroutine domain_regions_init

    subroutine domain_init_marine_ice(dom)
        ! LGM-like marine ice at the cold start (init_marine_H): thin
        ! ice (< 600 m) over shallow bed (> -500 m) is thickened to 800 m wherever
        ! ice is allowed.
        type(kryos_domain), intent(inout) :: dom

        where(dom%yelmo%bnd%mask_ice /= MASK_ICE_NONE .and. &
              dom%yelmo%tpo%now%H_ice < 600.0_wp .and. &
              dom%yelmo%bnd%z_bed > -500.0_wp)
            dom%yelmo%tpo%now%H_ice = 800.0_wp
        end where
    end subroutine domain_init_marine_ice

    subroutine negis_update_cb_ref(ylmo, ngs, time)
        ! Northeast Greenland Ice Stream cb_ref modification ([sim] use_negis,
        ! Yelmo's till_method = -1): recompute cb_ref from bed properties
        ! (calc_cb_ref), then scale the NEGIS parts (region_centre/south/north)
        ! by their factors at the forcing time `time`: the centre by cf_centre,
        ! the north part (paleo ice stream) by cf_north before -4 ka, the south
        ! part by cf_south from -4 ka, and all parts by the ramp f_x, cf_0 before
        ! -11 ka, linear to cf_1 at 0 ka, cf_1 after. The factors are functions
        ! of the time only, so a restart continues them. Requires the [negis]
        ! parameters, loaded by domain_negis_init.
        type(yelmo_class),  intent(inout) :: ylmo
        type(negis_params), intent(in)    :: ngs
        real(dp),           intent(in)    :: time

        real(dp), parameter :: t_ramp0 = -11e3_dp   ! [yr] start of the cf_0 -> cf_1 ramp
        real(dp), parameter :: t_ramp1 =   0.0_dp   ! [yr] end of the ramp
        real(dp), parameter :: t_shift =  -4e3_dp   ! [yr] north part -> south part
        real(wp) :: f_x, f_south, f_north
        real(dp) :: w

        w   = min(1.0_dp, max(0.0_dp, (time - t_ramp0) / (t_ramp1 - t_ramp0)))
        f_x = real(ngs%cf_0 + w*(ngs%cf_1 - ngs%cf_0), wp)

        if (time < t_shift) then
            f_north = ngs%cf_north
            f_south = 1.0_wp
        else
            f_north = 1.0_wp
            f_south = ngs%cf_south
        end if

        ! Recompute cb_ref like the standard till function.
        call calc_cb_ref(ylmo%dyn%now%cb_ref, ylmo%bnd%z_bed, ylmo%bnd%z_bed_sd, ylmo%bnd%z_sl, &
                ylmo%bnd%H_sed, ylmo%dyn%par%till_f_sed, ylmo%dyn%par%till_sed_min, ylmo%dyn%par%till_sed_max, &
                ylmo%dyn%par%till_cf_ref, ylmo%dyn%par%till_cf_min, ylmo%dyn%par%till_z0, ylmo%dyn%par%till_z1, &
                ylmo%dyn%par%till_n_sd, ylmo%dyn%par%till_scale_zb, ylmo%dyn%par%till_scale_sed)

        ! Apply the scaling of the NEGIS parts, then the ramp over all of them.
        where (ngs%mask_centre) ylmo%dyn%now%cb_ref = ylmo%dyn%now%cb_ref * ngs%cf_centre
        where (ngs%mask_south)  ylmo%dyn%now%cb_ref = ylmo%dyn%now%cb_ref * f_south
        where (ngs%mask_north)  ylmo%dyn%now%cb_ref = ylmo%dyn%now%cb_ref * f_north
        where (ngs%mask_centre .or. ngs%mask_south .or. ngs%mask_north) &
            ylmo%dyn%now%cb_ref = ylmo%dyn%now%cb_ref * f_x

    end subroutine negis_update_cb_ref

    subroutine calc_glacial_smb(smb, lat2D, ta_ann, ta_ann_pd, gsmb)
        ! Reduce (scale up toward zero) negative surface mass balance during
        ! glacial conditions, above a latitude limit. The
        ! glacial index is derived from the domain-mean cooling.
        real(wp), intent(inout) :: smb(:,:)
        real(wp), intent(in)    :: lat2D(:,:)
        real(wp), intent(in)    :: ta_ann(:,:)
        real(wp), intent(in)    :: ta_ann_pd(:,:)
        type(glacial_smb_params), intent(in) :: gsmb

        integer  :: i, j, nx, ny
        real(wp) :: t0, tnow, at

        nx = size(smb,1)
        ny = size(smb,2)

        ! Quasi glacial-interglacial index (0: interglacial, 1: glacial)
        tnow = sum(ta_ann)    / real(nx*ny,wp)
        t0   = sum(ta_ann_pd) / real(nx*ny,wp)
        at = (tnow-t0)/gsmb%dt_lgm
        if (at .lt. 0.0_wp) at = 0.0_wp
        if (at .gt. 1.0_wp) at = 1.0_wp

        do j = 1, ny
        do i = 1, nx
            if (smb(i,j) .lt. 0.0_wp .and. lat2D(i,j) .gt. gsmb%lat_lim) then
                smb(i,j) = smb(i,j) - smb(i,j) * at * gsmb%fac_lim
            end if
        end do
        end do
    end subroutine calc_glacial_smb

end module kryos_regions
