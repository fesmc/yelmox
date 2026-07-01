program test_htopo
    ! Load the hi-res topography reference hub from the ANT-16KM ice_data files
    ! and check the fields came in on the expected grid with sane ranges.

    use htopo

    implicit none

    type(htopo_class) :: ht
    integer :: fails

    fails = 0

    call htopo_init(ht, "tests/test_htopo.nml", "htopo", map_fldr="maps")

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

    if (fails > 0) stop 1
    write(*,*) "PASS: test_htopo"

end program test_htopo
