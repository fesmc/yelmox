program test_htopo
    ! Load the hi-res topography reference hub from the ANT-16KM ice_data files
    ! and check the fields came in on the expected grid with sane ranges.

    use htopo
    use ncio

    implicit none

    integer, parameter :: wp = kind(1.0)
    integer, parameter :: n_band = 10           ! width of the missing border band [cells]
    real(wp), parameter :: fill = -9.0e33_wp    ! missing value in the gaps file

    type(htopo_class) :: ht
    type(htopo_class) :: ht_nomask
    type(htopo_class) :: ht_gaps
    integer :: fails, i, j
    logical, allocatable :: gap(:,:), gap_srf(:,:)
    real(wp), allocatable :: z_srf_exp(:,:)

    fails = 0

    call htopo_init(ht, "tests/test_htopo.nml", "domain", "Antarctica", "ANT-16KM", map_fldr="maps")

    write(*,*) "htopo grid   : "//trim(ht%par%grid_name), " nx,ny =", ht%nx, ht%ny
    write(*,*) "z_bed  range :", minval(ht%z_bed),   maxval(ht%z_bed)
    write(*,*) "H_ice  range :", minval(ht%H_ice),   maxval(ht%H_ice)
    write(*,*) "z_srf  range :", minval(ht%z_srf),   maxval(ht%z_srf)
    write(*,*) "regions range:", minval(ht%regions), maxval(ht%regions)
    write(*,*) "basins range :", minval(ht%basins),  maxval(ht%basins)

    if (ht%nx /= 381 .or. ht%ny /= 381) then
        write(*,*) "FAIL: unexpected topo grid size"; fails = fails + 1
    end if
    if (maxval(ht%H_ice) < 1000.0) then
        write(*,*) "FAIL: H_ice looks empty"; fails = fails + 1
    end if
    if (minval(ht%z_bed) > 0.0) then
        write(*,*) "FAIL: z_bed has no ocean floor"; fails = fails + 1
    end if
    if (maxval(ht%basins) < 1.0) then
        write(*,*) "FAIL: basins look empty"; fails = fails + 1
    end if
    if (maxval(ht%z_bed_sd) <= 0.0) then
        write(*,*) "FAIL: z_bed_sd looks empty"; fails = fails + 1
    end if

    ! Blank mask paths and z_bed_sd name: nothing is read and the masks
    ! default to 1, z_bed_sd to 0.
    call htopo_init(ht_nomask, "tests/test_htopo.nml", "domain_nomask", "Antarctica", "ANT-16KM", &
                    map_fldr="maps")

    if (minval(ht_nomask%regions) /= 1.0 .or. maxval(ht_nomask%regions) /= 1.0 .or. &
        minval(ht_nomask%basins)  /= 1.0 .or. maxval(ht_nomask%basins)  /= 1.0) then
        write(*,*) "FAIL: blank mask paths did not give masks of 1"; fails = fails + 1
    end if
    if (maxval(abs(ht_nomask%z_bed_sd)) /= 0.0) then
        write(*,*) "FAIL: blank z_bed_sd name did not give z_bed_sd = 0"; fails = fails + 1
    end if
    if (maxval(abs(ht_nomask%z_bed - ht%z_bed)) /= 0.0) then
        write(*,*) "FAIL: blank mask paths changed the topography"; fails = fails + 1
    end if

    ! Data gaps: a copy of the topography with all fields missing in a band
    ! along the x = min border, and z_srf also missing on some ice cells.
    allocate(gap(ht%nx,ht%ny), gap_srf(ht%nx,ht%ny))
    gap = .false.
    gap(1:n_band,:) = .true.
    gap_srf = gap
    do j = 1, ht%ny
    do i = 1, ht%nx
        if (ht%H_ice(i,j) > 0.0 .and. mod(i+j,50) == 0) gap_srf(i,j) = .true.
    end do
    end do

    call write_gaps_file("test_htopo_gaps.nc", ht, gap, gap_srf)
    call htopo_init(ht_gaps, "tests/test_htopo.nml", "domain_gaps", "Antarctica", "ANT-16KM", &
                    map_fldr="maps")

    ! No ice in the gaps; the bed from the nearest valid cell (the first column
    ! after the band, for rows away from the y borders); the surface from the
    ! bed and the ice thickness at sea level 0.
    allocate(z_srf_exp(ht%nx,ht%ny))
    z_srf_exp = max(ht_gaps%z_bed + ht_gaps%H_ice, (1.0_wp - 910.0_wp/1028.0_wp)*ht_gaps%H_ice)

    if (any(ht_gaps%H_ice /= 0.0 .and. gap)) then
        write(*,*) "FAIL: gap cells have ice"; fails = fails + 1
    end if
    do j = 5, ht%ny-4
        if (any(ht_gaps%z_bed(1:n_band,j) /= ht%z_bed(n_band+1,j))) then
            write(*,*) "FAIL: gap z_bed is not the nearest valid bed, row ", j; fails = fails + 1
            exit
        end if
    end do
    if (any(gap_srf .and. abs(ht_gaps%z_srf - z_srf_exp) > 1e-3)) then
        write(*,*) "FAIL: gap z_srf does not follow z_bed and H_ice"; fails = fails + 1
    end if
    if (any(.not. gap     .and. ht_gaps%z_bed /= ht%z_bed) .or. &
        any(.not. gap     .and. ht_gaps%H_ice /= ht%H_ice) .or. &
        any(.not. gap_srf .and. ht_gaps%z_srf /= ht%z_srf)) then
        write(*,*) "FAIL: valid cells changed"; fails = fails + 1
    end if

    call delete_file("test_htopo_gaps.nc")

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
        call nc_write(filename, "z_bed", merge(fill, ht%z_bed, gap),     dim1="xc", dim2="yc", &
                      missing_value=fill)
        call nc_write(filename, "H_ice", merge(fill, ht%H_ice, gap),     dim1="xc", dim2="yc", &
                      missing_value=fill)
        call nc_write(filename, "z_srf", merge(fill, ht%z_srf, gap_srf), dim1="xc", dim2="yc", &
                      missing_value=fill)
    end subroutine write_gaps_file

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: u
        open(newunit=u, file=filename, status="old")
        close(u, status="delete")
    end subroutine delete_file

end program test_htopo
