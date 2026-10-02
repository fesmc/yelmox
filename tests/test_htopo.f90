program test_htopo
    ! Load the hi-res topography reference hub from the ANT-16KM ice_data files
    ! and check the fields came in on the expected grid with sane ranges.

    use htopo
    use ncio
    use phys_constants, only : phys_const_class, phys_const_load

    implicit none

    integer, parameter :: wp = kind(1.0)
    integer, parameter :: n_band = 10           ! width of the missing border band [cells]
    real(wp), parameter :: fill = -9.0e33_wp    ! missing value in the gaps file

    type(htopo_class) :: ht
    type(htopo_class) :: ht_nomask
    type(htopo_class) :: ht_gaps
    type(phys_const_class) :: cnst
    integer :: fails, i, j
    logical, allocatable :: gap(:,:), gap_srf(:,:)
    real(wp), allocatable :: z_srf_exp(:,:), zero(:,:)

    fails = 0

    call phys_const_load(cnst, "input/yelmo_phys_const.nml", group="Earth")

    call htopo_init(ht, "tests/test_htopo.nml", "domain", "Antarctica", "ANT-16KM", cnst, map_fldr="maps")

    write(*,*) "htopo grid   : "//trim(ht%par%grid_name), " nx,ny =", ht%nx, ht%ny
    write(*,*) "z_bed  range :", minval(ht%z_bed_ref),   maxval(ht%z_bed_ref)
    write(*,*) "H_ice  range :", minval(ht%H_ice_ref),   maxval(ht%H_ice_ref)
    write(*,*) "z_srf  range :", minval(ht%z_srf_ref),   maxval(ht%z_srf_ref)
    write(*,*) "regions range:", minval(ht%regions), maxval(ht%regions)
    write(*,*) "basins range :", minval(ht%basins),  maxval(ht%basins)

    if (ht%nx /= 381 .or. ht%ny /= 381) then
        write(*,*) "FAIL: unexpected topo grid size"; fails = fails + 1
    end if
    if (maxval(ht%H_ice_ref) < 1000.0) then
        write(*,*) "FAIL: H_ice looks empty"; fails = fails + 1
    end if
    if (minval(ht%z_bed_ref) > 0.0) then
        write(*,*) "FAIL: z_bed has no ocean floor"; fails = fails + 1
    end if
    if (maxval(ht%basins) < 1.0) then
        write(*,*) "FAIL: basins look empty"; fails = fails + 1
    end if
    if (maxval(ht%z_bed_sd) <= 0.0) then
        write(*,*) "FAIL: z_bed_sd looks empty"; fails = fails + 1
    end if

    ! Sectors and the named regions (APIS/WAIS/EAIS = sectors 3/1/2).
    if (minval(ht%sectors) /= 0.0 .or. maxval(ht%sectors) /= 3.0) then
        write(*,*) "FAIL: sectors range is not 0..3"; fails = fails + 1
    end if
    if (ht%par%n_regions /= 3 .or. trim(ht%par%region_names(2)) /= "WAIS" .or. &
        ht%par%region_codes(2) /= 1.0) then
        write(*,*) "FAIL: named regions not read"; fails = fails + 1
    end if
    if (any(htopo_region_codes(ht) /= ht%sectors)) then
        write(*,*) "FAIL: region_mask = sectors does not give the sectors"; fails = fails + 1
    end if

    ! Ice allowed everywhere except the open ocean (regions code 2.0).
    if (any(htopo_ice_allowed(ht%par, ht%regions) .neqv. (ht%regions /= 2.0))) then
        write(*,*) "FAIL: exclude ice_codes"; fails = fails + 1
    end if
    if (count(.not. htopo_ice_allowed(ht%par, ht%regions)) == 0) then
        write(*,*) "FAIL: no open ocean in the regions"; fails = fails + 1
    end if

    ! Blank mask paths and z_bed_sd name: nothing is read and the masks
    ! default to 1, z_bed_sd to 0.
    call htopo_init(ht_nomask, "tests/test_htopo.nml", "domain_nomask", "Antarctica", "ANT-16KM", cnst, &
                    map_fldr="maps")

    if (minval(ht_nomask%regions) /= 1.0 .or. maxval(ht_nomask%regions) /= 1.0 .or. &
        minval(ht_nomask%basins)  /= 1.0 .or. maxval(ht_nomask%basins)  /= 1.0) then
        write(*,*) "FAIL: blank mask paths did not give masks of 1"; fails = fails + 1
    end if
    if (maxval(abs(ht_nomask%z_bed_sd)) /= 0.0) then
        write(*,*) "FAIL: blank z_bed_sd name did not give z_bed_sd = 0"; fails = fails + 1
    end if
    if (maxval(abs(ht_nomask%z_bed_ref - ht%z_bed_ref)) /= 0.0) then
        write(*,*) "FAIL: blank mask paths changed the topography"; fails = fails + 1
    end if
    if (minval(ht_nomask%sectors) /= 1.0 .or. maxval(ht_nomask%sectors) /= 1.0) then
        write(*,*) "FAIL: blank sectors path did not give sectors of 1"; fails = fails + 1
    end if
    if (ht_nomask%par%n_regions /= 0 .or. ht_nomask%par%n_ice_codes /= 0) then
        write(*,*) "FAIL: blank code lists are not empty"; fails = fails + 1
    end if
    if (.not. all(htopo_ice_allowed(ht_nomask%par, ht%regions))) then
        write(*,*) "FAIL: ice_codes_mode = all does not allow ice everywhere"; fails = fails + 1
    end if

    ! Include: ice only on the given codes.
    ht_nomask%par%ice_codes_mode = "include"
    ht_nomask%par%ice_codes(1:2) = [1.0, 3.0]
    ht_nomask%par%n_ice_codes    = 2
    if (any(htopo_ice_allowed(ht_nomask%par, ht%regions) .neqv. &
            (ht%regions == 1.0 .or. ht%regions == 3.0))) then
        write(*,*) "FAIL: include ice_codes"; fails = fails + 1
    end if

    ! Relaxation: none by default; with exclude, relax_tau everywhere but on
    ! the codes, -1 (free) on them.
    if (any(htopo_relax_tau(ht%par, ht%regions) /= -1.0)) then
        write(*,*) "FAIL: relax_codes_mode = none does not leave tau_relax at -1"; fails = fails + 1
    end if
    ht_nomask%par%relax_codes_mode = "exclude"
    ht_nomask%par%relax_codes(1)   = 1.0
    ht_nomask%par%n_relax_codes    = 1
    ht_nomask%par%relax_tau        = 50.0
    if (any(htopo_relax_tau(ht_nomask%par, ht%regions) /= &
            merge(-1.0, 50.0, ht%regions == 1.0))) then
        write(*,*) "FAIL: exclude relax_codes"; fails = fails + 1
    end if

    ! Data gaps: a copy of the topography with all fields missing in a band
    ! along the x = min border, and z_srf also missing on some ice cells.
    allocate(gap(ht%nx,ht%ny), gap_srf(ht%nx,ht%ny))
    gap = .false.
    gap(1:n_band,:) = .true.
    gap_srf = gap
    do j = 1, ht%ny
    do i = 1, ht%nx
        if (ht%H_ice_ref(i,j) > 0.0 .and. mod(i+j,50) == 0) gap_srf(i,j) = .true.
    end do
    end do

    call write_gaps_file("test_htopo_gaps.nc", ht, gap, gap_srf)
    call htopo_init(ht_gaps, "tests/test_htopo.nml", "domain_gaps", "Antarctica", "ANT-16KM", cnst, &
                    map_fldr="maps")

    ! No ice in the gaps; the bed from the nearest valid cell (the first column
    ! after the band, for rows away from the y borders); the surface from the
    ! bed and the ice thickness at sea level 0.
    allocate(z_srf_exp(ht%nx,ht%ny))
    z_srf_exp = max(ht_gaps%z_bed_ref + ht_gaps%H_ice_ref, (1.0_wp - 910.0_wp/1028.0_wp)*ht_gaps%H_ice_ref)

    if (any(ht_gaps%H_ice_ref /= 0.0 .and. gap)) then
        write(*,*) "FAIL: gap cells have ice"; fails = fails + 1
    end if
    do j = 5, ht%ny-4
        if (any(ht_gaps%z_bed_ref(1:n_band,j) /= ht%z_bed_ref(n_band+1,j))) then
            write(*,*) "FAIL: gap z_bed is not the nearest valid bed, row ", j; fails = fails + 1
            exit
        end if
    end do
    if (any(gap_srf .and. abs(ht_gaps%z_srf_ref - z_srf_exp) > 1e-3)) then
        write(*,*) "FAIL: gap z_srf does not follow z_bed and H_ice"; fails = fails + 1
    end if
    if (any(.not. gap     .and. ht_gaps%z_bed_ref /= ht%z_bed_ref) .or. &
        any(.not. gap     .and. ht_gaps%H_ice_ref /= ht%H_ice_ref) .or. &
        any(.not. gap_srf .and. ht_gaps%z_srf_ref /= ht%z_srf_ref)) then
        write(*,*) "FAIL: valid cells changed"; fails = fails + 1
    end if

    call delete_file("test_htopo_gaps.nc")

    ! Current geometry from the reference plus anomalies (htopo_update). The
    ! current geometry starts from the reference; no anomaly keeps it, with
    ! the grounding and surface from flotation at sea level 0.
    if (any(ht%z_bed /= ht%z_bed_ref) .or. any(ht%H_ice /= ht%H_ice_ref)) then
        write(*,*) "FAIL: current geometry does not start from the reference"; fails = fails + 1
    end if
    allocate(zero(ht%nx,ht%ny)); zero = 0.0_wp
    call htopo_update(ht, zero, zero, zero)
    if (any(ht%z_bed /= ht%z_bed_ref) .or. any(ht%H_ice /= max(ht%H_ice_ref, 0.0_wp)) .or. &
        any(ht%z_sl /= 0.0)) then
        write(*,*) "FAIL: htopo_update without anomalies changed the geometry"; fails = fails + 1
    end if
    if (any(ht%f_grnd == 1.0 .neqv. (ht%z_bed >= 0.0 .or. &
            ht%H_ice - (1028.0_wp/910.0_wp)*(0.0_wp - ht%z_bed) >= 0.0))) then
        write(*,*) "FAIL: grounding does not follow flotation"; fails = fails + 1
    end if
    if (any(ht%z_srf /= max(ht%z_bed + ht%H_ice, (1.0_wp - 910.0_wp/1028.0_wp)*ht%H_ice))) then
        write(*,*) "FAIL: surface does not follow bed, ice and sea level"; fails = fails + 1
    end if
    if (count(ht%f_grnd == 0.0 .and. ht%H_ice > 0.0) == 0) then
        write(*,*) "FAIL: no floating ice in the reference"; fails = fails + 1
    end if

    ! Bed displacement and ice change add to the reference; ice is clipped at 0.
    call htopo_update(ht, zero - 100.0_wp, zero - 500.0_wp, zero)
    if (any(ht%z_bed /= ht%z_bed_ref - 100.0_wp) .or. &
        any(ht%H_ice /= max(ht%H_ice_ref - 500.0_wp, 0.0_wp))) then
        write(*,*) "FAIL: htopo_update anomalies"; fails = fails + 1
    end if

    if (fails > 0) stop 1
    write(*,*) "PASS: test_htopo"

contains

    subroutine write_gaps_file(filename, ht, gap, gap_srf)
        ! Write z_bed, H_ice, z_srf of `ht`, with `fill` where gap / gap_srf.
        character(len=*),  intent(in) :: filename
        type(htopo_class), intent(in) :: ht
        logical,           intent(in) :: gap(:,:), gap_srf(:,:)

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=ht%grid%G%x, units="km")
        call nc_write_dim(filename, "yc", x=ht%grid%G%y, units="km")
        call nc_write(filename, "z_bed", merge(fill, ht%z_bed_ref, gap),     dim1="xc", dim2="yc", &
                      missing_value=fill)
        call nc_write(filename, "H_ice", merge(fill, ht%H_ice_ref, gap),     dim1="xc", dim2="yc", &
                      missing_value=fill)
        call nc_write(filename, "z_srf", merge(fill, ht%z_srf_ref, gap_srf), dim1="xc", dim2="yc", &
                      missing_value=fill)
    end subroutine write_gaps_file

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: u
        open(newunit=u, file=filename, status="old")
        close(u, status="delete")
    end subroutine delete_file

end program test_htopo
