program yelmox_mg
    ! Multigrid yelmox driver.
    !
    ! Bring-up stage: loads the hi-res topography reference hub (htopo) from the
    ! co-located parameter file and reports it. domain_init (sub-model init on
    ! their own grids + coupler map priming) and the coupling time loop are
    ! filled in incrementally, lifting and adapting yelmox.f90 into the
    ! ice_domain / step_* structure (see docs/multigrid.md).

    use nml
    use yelmo, only : yelmo_load_command_line_args
    use htopo, only : htopo_class, htopo_init

    implicit none

    character(len=512) :: path_par
    type(htopo_class)  :: topo

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! Hi-res reference geometry: the top level of the multigrid setup.
    call htopo_init(topo, path_par, "htopo")

    write(*,*) "yelmox_mg: htopo reference loaded"
    write(*,*) "  grid        : "//trim(topo%par%grid_name), topo%nx, topo%ny
    write(*,*) "  H_ice  max  :", maxval(topo%H_ice)
    write(*,*) "  z_bed  range:", minval(topo%z_bed), maxval(topo%z_bed)
    write(*,*) "yelmox_mg: init smoke test complete (domain_init + loop to follow)."

end program yelmox_mg
