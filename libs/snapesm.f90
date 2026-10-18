module snapesm
    ! snapesm — climate boundary forcing built on the fesm-utils primitives.
    !
    !   Loading   -> varslice  (one varslice_class per snapshot x field)
    !   Indices   -> the driver's registry (kryos_forcing), read by name
    !   Wrapper   -> this module (blend + transforms + output state)
    !
    ! Replacement for the monolithic `snapclim`. See docs/climate-snap.md for its
    ! configuration.
    !
    ! Config model: a `var_defs` database nml defines varslice groups (permanent
    ! per-file variable metadata); each climate-state group &<group>_<snapshot> names
    ! the group(s) supplying each field (1 = monthly, 2 = [ann, sum] -> synthesized).
    !
    ! Axes: each axis &<group>_axis_<name> spans two snapshots (lo -> hi) and
    ! names an index of the driver's registry ([indices], kryos_forcing); the
    ! value of the index is the position on the axis (0 at lo, 1 at hi). A field
    ! sums the axes it lists on top of the reference snapshot (Banderas et al.,
    ! 2018, for two axes). The indices restart with the registry, so snapesm
    ! holds no state of its own.
    !
    ! The ocean temperature as a fraction of the air temperature anomaly is the
    ! climate modifier ocn_fraction (climate_mods), on the field `to` held at
    ! the reference (blend = "const").
    !
    ! STATUS: the pipeline (reduce / combine / transform) is ported and
    ! checked against snapclim (tests/test_snapesm_ref.f90); output
    ! (snapesm_write_*) is in place; snapesm_end is a stub.

    use precision, only : wp, sp, dp
    use ncio
    use nml
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use kryos_forcing, only : indices_class, indices_value, indices_has
    use gaussian_filter, only : filter_gaussian
    use coordinates, only : grid_class

    implicit none

    real(wp), parameter :: pi = 3.14159265358979323846_wp

    private

    integer, parameter :: SNP_MAX = 32   ! max entries in a snapshots/fields/axes list
    integer, parameter :: AXS_MAX = 8    ! max axes of one field
    integer, parameter :: NMONTH   = 12

    ! =====================================================================
    ! Types
    ! =====================================================================

    ! Ad-hoc output field, accessed by name (escape hatch for the long tail;
    ! hot fields are named members of climate_state_class instead).
    type named_field_class
        character(len=32)     :: name
        real(wp), allocatable :: var(:,:,:)
    end type

    ! One physical field in the registry (immutable after load).
    type field_spec_class
        character(len=32) :: name           ! "tas","pr","to","so","zs","smb","bmb_shlf",...
        character(len=16) :: kind           ! "atm_temp" | "atm_precip" | "ocn" | "ocn_salt" | "elev"
        character(len=16) :: combine        ! "anomaly" | "absolute"
        character(len=16) :: blend          ! "linear" | "ratio" | "const"
        character(len=64) :: axes(AXS_MAX)  ! names of the axes of this field ("" = none)
        integer           :: n_axis         ! number of axes
        integer           :: axis_k(AXS_MAX) ! their positions in snapesm_class%axis
        logical           :: enabled        ! off by default for exotic fields
        logical           :: apply_lapse
        logical           :: precip_scaling
        logical           :: seasonal_synth
    end type

    ! One climate state: PD, LGM, piControl, ... (config side).
    type snapshot_spec_class
        character(len=64) :: name
        logical           :: is_ref         ! the distinguished reference state
        logical           :: monthly        ! raw atm fields are monthly (else annual+summer -> synthesized)
        real(wp)          :: clim_sigma     ! [km] Gaussian smoothing of the atmosphere at loading (0: none)
        real(wp)          :: ocn_sigma      ! [km] Gaussian smoothing of the ocean at loading (0: none)
    end type

    ! A climate state. Option A: named members for the hot path, unallocated when a
    ! run does not use them (allocated() is the "in use" guard); `extra` for the tail.
    ! Used both as a snapshot's reduced (sea-level) fields and as the blended output.
    type climate_state_class
        ! atmosphere, monthly (nx,ny,12): tas at surface, tsl at sea level
        real(wp), allocatable :: tas(:,:,:), tsl(:,:,:), pr(:,:,:), prcor(:,:,:)
        real(wp), allocatable :: sf(:,:,:), pr_stdev_frac(:,:,:)
        ! atmosphere, annual / summer (nx,ny)
        real(wp), allocatable :: ta_ann(:,:), ta_sum(:,:), pr_ann(:,:)
        real(wp), allocatable :: tsl_ann(:,:), tsl_sum(:,:), prcor_ann(:,:)
        real(wp), allocatable :: beta_p(:,:)          ! precip-temperature sensitivity
        ! ocean
        real(wp), allocatable :: to_ann(:,:,:), so_ann(:,:,:)
        real(wp), allocatable :: depth(:)
        ! geometry / masks
        real(wp), allocatable :: mask(:,:), z_srf(:,:)
        ! exotic / expansion forcing (named for the common ones)
        real(wp), allocatable :: smb(:,:), bmb_shlf(:,:,:)
        ! escape hatch for genuinely ad-hoc future fields (temp-copy access only)
        type(named_field_class), allocatable :: extra(:)
    end type

    ! A field's source binding for one snapshot: one varslice group (already monthly)
    ! or two (annual mean + summer) that `reduce` synthesizes into 12 months. The group
    ! names are references into the var_defs database file.
    type field_binding_class
        type(varslice_class), allocatable :: src(:)   ! size 1 (monthly) or 2 ([ann, sum])
    end type

    ! One loaded snapshot: its spec, a per-field source binding (registry order), and
    ! its reduced (sea-level) climate state used by the blend.
    type snapshot_class
        type(snapshot_spec_class)              :: spec
        type(field_binding_class), allocatable :: bind(:)
        type(climate_state_class)              :: state
    end type

    ! Top-level parameters.
    type snapesm_param_class
        character(len=256) :: domain
        character(len=256) :: grid_name
        character(len=256) :: var_defs      ! path to the varslice variable-database nml
        character(len=64)  :: group         ! base namelist group (already domain-suffixed)
        integer            :: nx, ny
        real(wp)           :: dx_km         ! [km] grid spacing (snapshot smoothing)
        logical            :: south         ! southern hemisphere (seasons, lapse rates)
        character(len=16)  :: combine       ! default combine mode
        character(len=64)  :: ref_name      ! which snapshot is the reference
        real(wp)           :: lapse(2)      ! [annual, summer] lapse rate [K/m]
        ! precip sensitivity parameters (snapclim &snap)
        real(wp)           :: f_p           ! precip-temperature Clausius factor -> beta_p [1/K]
        real(wp)           :: f_p_ne        ! NE-Greenland beta_p multiplier (snapclim hardcodes 1.0)
        real(wp)           :: f_stdev       ! precip-variability multiplier
        integer            :: n_snap, n_field, n_axis
    end type

    ! One axis: two snapshots and the registry index giving the position on it.
    type axis_class
        character(len=64) :: name
        character(len=64) :: index          ! registry index (kryos_forcing)
        integer           :: lo, hi         ! snapshots at index = 0 and 1 (positions in snap)
        real(wp)          :: f              ! position at the last update (output)
    end type

    ! The object callers hold.
    type snapesm_class
        type(snapesm_param_class)          :: par
        type(field_spec_class),  allocatable :: registry(:)  ! field definitions
        type(snapshot_class),    allocatable :: snap(:)      ! was clim0..clim3
        type(axis_class),        allocatable :: axis(:)      ! snapshot axes
        type(climate_state_class)            :: now          ! current output
        type(climate_state_class)            :: ref          ! distinguished reference (was clim0)
    end type

    public :: snapesm_class
    public :: snapesm_init
    public :: snapesm_update
    public :: snapesm_end
    public :: snapesm_write_init
    public :: snapesm_write_step
    public :: snapesm_lapse_monthly

contains

    ! =====================================================================
    ! Lifecycle
    ! =====================================================================

    subroutine snapesm_init(sc, filename, domain, grid, time, south, idx, group)
        ! `grid` is the climate grid: the snapshot fields are read on it, or
        ! remapped onto it online (var_defs groups with remap/grid_src). `time`
        ! primes the time-varying snapshots. `idx` is the driver's registry: each
        ! axis must name one of its indices.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: domain
        type(grid_class),      intent(IN)    :: grid
        real(dp),              intent(IN)    :: time
        logical,               intent(IN)    :: south     ! southern hemisphere
        type(indices_class),   intent(IN)    :: idx
        character(len=*),      intent(IN), optional :: group

        character(len=64) :: base_group

        base_group = "snp"
        if (present(group)) base_group = trim(group)

        ! Read all configuration: top-level, field registry, snapshot specs, axes.
        call snapesm_par_load(sc, filename, trim(base_group), domain, trim(grid%name), idx)
        sc%par%nx = grid%G%nx
        sc%par%ny = grid%G%ny
        sc%par%dx_km = real(grid%G%dx*grid%cs%xy_conv*1e-3_dp, wp)   ! [km] (clim_sigma)
        sc%par%south = south

        ! Load the snapshot fields (varslice), reduce them, and set the reference state.
        call snapesm_load_snapshots(sc, filename, time, grid)

        return
    end subroutine snapesm_init

    subroutine snapesm_update(sc, z_srf, time, idx, dx)
        ! Mirrors snapclim_update. Pipeline: axis positions from the registry
        ! (advanced by the driver to `time`) -> refresh loads -> combine ->
        ! transform. Homogeneous anomalies and the ocean fraction rule are
        ! climate modifiers (climate_mods), applied by the caller afterwards.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        real(dp),              intent(IN)    :: time
        type(indices_class),   intent(IN)    :: idx
        real(wp),              intent(IN), optional :: dx

        integer :: k

        ! 1. The position on each axis: the current value of its index.
        do k = 1, sc%par%n_axis
            sc%axis(k)%f = indices_value(idx, sc%axis(k)%index)
        end do

        ! 2. Refresh time-varying loads and re-reduce those snapshots.
        call snapesm_refresh_loads(sc, time)

        ! 3. Combine the reduced snapshot states into `now`, per field.
        call snapesm_combine(sc)

        ! 4. Transform: inflate sea-level fields to the current model elevation
        !    (tsl -> tas), precip, annual/summer aggregates.
        call snapesm_transform(sc, z_srf)

        return
    end subroutine snapesm_update

    ! ---------------------------------------------------------------------
    ! Update pipeline stages
    ! ---------------------------------------------------------------------

    subroutine snapesm_refresh_loads(sc, time)
        ! Re-slice only the time-varying fields (recon, transient orography, etc.);
        ! static fields were primed once at load. A re-reduce of touched snapshots
        ! follows once the reduction is implemented.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(dp),              intent(IN)    :: time
        integer :: s, f, k
        logical :: touched
        do s = 1, sc%par%n_snap
            if (.not. allocated(sc%snap(s)%bind)) cycle
            touched = .FALSE.
            do f = 1, sc%par%n_field
                if (.not. allocated(sc%snap(s)%bind(f)%src)) cycle
                do k = 1, size(sc%snap(s)%bind(f)%src)
                    if (sc%snap(s)%bind(f)%src(k)%par%with_time) then
                        call varslice_update(sc%snap(s)%bind(f)%src(k), [time], method="extrap")
                        touched = .TRUE.
                    end if
                end do
            end do
            ! Re-reduce only if a time-varying source was refreshed.
            if (touched) call snapesm_reduce_snapshot(sc, s)
        end do
        return
    end subroutine snapesm_refresh_loads

    subroutine snapesm_combine(sc)
        ! Blend the reduced snapshot states into `now`, per field, from the reference
        ! baseline (was clim0). Each axis k of a field (lo_k -> hi_k) adds with
        ! weight a_k = the value of its index:
        !   linear (temp/salinity): now = ref + sum_k a_k*(hi_k - lo_k)
        !   ratio  (precip):        now = ref * prod_k (a_k*(hi_k/lo_k - 1) + 1)
        !   const:                  now = ref
        ! One axis reproduces snapclim calc_temp_1ind/calc_precip_1ind, two axes
        ! calc_temp_2ind/calc_precip_2ind (Banderas et al., 2018).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        integer  :: f, j, k, lo, hi
        real(wp) :: aa
        logical  :: blended

        do f = 1, sc%par%n_field
            if (.not. sc%registry(f)%enabled) cycle

            ! Start from the reference state.
            select case(trim(sc%registry(f)%kind))
                case("atm_temp")
                    if (allocated(sc%now%tsl)    .and. allocated(sc%ref%tsl))    sc%now%tsl    = sc%ref%tsl
                case("atm_precip")
                    if (allocated(sc%now%prcor)  .and. allocated(sc%ref%prcor))  sc%now%prcor  = sc%ref%prcor
                case("ocn")
                    if (allocated(sc%now%to_ann) .and. allocated(sc%ref%to_ann)) sc%now%to_ann = sc%ref%to_ann
                case("ocn_salt")
                    if (allocated(sc%now%so_ann) .and. allocated(sc%ref%so_ann)) sc%now%so_ann = sc%ref%so_ann
                case default
                    cycle   ! "elev" and other non-blended fields: nothing to combine.
            end select

            ! linear/ratio fields add their axes; const fields stay at the
            ! reference.
            blended = trim(sc%registry(f)%blend) == "linear" .or. &
                      trim(sc%registry(f)%blend) == "ratio"
            if (.not. blended) cycle

            do j = 1, sc%registry(f)%n_axis
                k  = sc%registry(f)%axis_k(j)
                lo = sc%axis(k)%lo
                hi = sc%axis(k)%hi
                aa = sc%axis(k)%f

                associate(blend => sc%registry(f)%blend, x0 => sc%snap(lo)%state, x1 => sc%snap(hi)%state)
                select case(trim(sc%registry(f)%kind))
                    case("atm_temp")
                        call blend_axis(blend, sc%now%tsl,    x0%tsl,    x1%tsl,    aa)
                    case("atm_precip")
                        call blend_axis(blend, sc%now%prcor,  x0%prcor,  x1%prcor,  aa)
                    case("ocn")
                        call blend_axis(blend, sc%now%to_ann, x0%to_ann, x1%to_ann, aa)
                    case("ocn_salt")
                        call blend_axis(blend, sc%now%so_ann, x0%so_ann, x1%so_ann, aa)
                end select
                end associate
            end do
        end do

        return
    end subroutine snapesm_combine

    subroutine blend_axis(mode, dst, x0, x1, aa)
        ! Add one axis (x0 -> x1, weight aa) to dst, which holds the reference
        ! plus the axes added so far. No-op when the axis snapshots do not supply
        ! the field (x0 or x1 unallocated).
        implicit none
        character(len=*),      intent(IN)    :: mode
        real(wp), allocatable, intent(INOUT) :: dst(:,:,:)
        real(wp), allocatable, intent(IN)    :: x0(:,:,:), x1(:,:,:)
        real(wp),              intent(IN)    :: aa

        if (.not. (allocated(dst) .and. allocated(x0) .and. allocated(x1))) return

        select case(trim(mode))
            case("ratio")
                where (x0 /= 0.0_wp) dst = dst*(aa*(x1/x0 - 1.0_wp) + 1.0_wp)
            case default   ! "linear"
                dst = dst + aa*(x1 - x0)
        end select

        return
    end subroutine blend_axis

    subroutine snapesm_transform(sc, z_srf)
        ! Inflate the blended sea-level fields to the current model elevation and
        ! compute the annual/summer aggregates (snapclim 779-811). Note the cosine
        ! phase here is m*30-15 (vs m*30-30 in the sea-level reduction).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        integer  :: m
        logical  :: south
        real(wp) :: lapse(NMONTH)

        south = sc%par%south
        lapse = snapesm_lapse_monthly(sc)

        if (allocated(sc%now%z_srf)) sc%now%z_srf = z_srf
        if (allocated(sc%now%mask)) then
            sc%now%mask = 0.0_wp
            where (z_srf > 1.0_wp) sc%now%mask = 1.0_wp
        end if

        ! 3a: monthly tas from tsl via seasonal lapse rate to the current elevation.
        do m = 1, NMONTH
            sc%now%tas(:,:,m) = sc%now%tsl(:,:,m) - z_srf*lapse(m)
        end do

        ! 3b: monthly precip re-sensitized to the current elevation.
        do m = 1, NMONTH
            sc%now%pr(:,:,m) = sc%now%prcor(:,:,m) * &
                exp(sc%now%beta_p*(sc%now%tas(:,:,m) - sc%now%tsl(:,:,m)))
        end do

        ! Step 4: annual and summer aggregates.
        sc%now%tsl_ann = sum(sc%now%tsl, dim=3) / 12.0_wp
        sc%now%ta_ann  = sum(sc%now%tas, dim=3) / 12.0_wp
        if (south) then
            sc%now%tsl_sum = sum(sc%now%tsl(:,:,[12,1,2]), dim=3) / 3.0_wp
            sc%now%ta_sum  = sum(sc%now%tas(:,:,[12,1,2]), dim=3) / 3.0_wp
        else
            sc%now%tsl_sum = sum(sc%now%tsl(:,:,[6,7,8]),  dim=3) / 3.0_wp
            sc%now%ta_sum  = sum(sc%now%tas(:,:,[6,7,8]),  dim=3) / 3.0_wp
        end if
        sc%now%prcor_ann = sum(sc%now%prcor, dim=3) / 12.0_wp * 365.0_wp
        sc%now%pr_ann    = sum(sc%now%pr,    dim=3) / 12.0_wp * 365.0_wp

        return
    end subroutine snapesm_transform

    function snapesm_lapse_monthly(sc) result(lapse)
        ! Monthly lapse rate [K/m] that lifts tsl to the surface (snapclim 788-795):
        ! annual lapse(1) with a seasonal cosine of amplitude lapse(1)-lapse(2),
        ! phase m*30-15, summer at lapse(2).
        implicit none
        type(snapesm_class), intent(IN) :: sc
        real(wp) :: lapse(NMONTH)
        integer  :: m
        real(wp) :: l1, l2
        l1 = sc%par%lapse(1)
        l2 = sc%par%lapse(2)
        do m = 1, NMONTH
            if (sc%par%south) then
                lapse(m) = l1 + (l2-l1)*cos(2*pi*(m*30.0_wp-15.0_wp)/360.0_wp)
            else
                lapse(m) = l1 + (l1-l2)*cos(2*pi*(m*30.0_wp-15.0_wp)/360.0_wp)
            end if
        end do
    end function snapesm_lapse_monthly

    subroutine snapesm_end(sc)
        ! Teardown (snapclim never had one). Free varslice and state arrays.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc

        ! TODO: varslice_end on each snap%bind%src, deallocate snap/axis/registry and states.

        return
    end subroutine snapesm_end

    ! =====================================================================
    ! Diagnostic output
    ! =====================================================================

    subroutine snapesm_write_init(sc, filename)
        ! Add the month axis of the snapesm fields to an existing 2D output file
        ! (xc, yc, time and depth defined by the caller).
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename

        call nc_write_dim(filename, "month", x=1, dx=1, nx=NMONTH, units="month")

        return
    end subroutine snapesm_write_init

    subroutine snapesm_write_step(sc, filename, ncid, n)
        ! Write record n of the snapesm state into an open 2D output file
        ! (axes from snapesm_write_init): the axis positions, the surface
        ! elevation, the monthly sea-level and the summer temperature -- the
        ! fields the enabled registry fills. The backend's state, before the
        ! climate modifiers; the applied atmosphere (tas, pr, tsl_ann) and ocean
        ! are written by the caller (climate_write_2D).
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        integer,               intent(IN) :: ncid, n

        integer :: k

        do k = 1, sc%par%n_axis
            call nc_write(filename, "axis_"//trim(sc%axis(k)%name), sc%axis(k)%f, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units="1", &
                          long_name="Position on axis "//trim(sc%axis(k)%name)//" (index "// &
                                    trim(sc%axis(k)%index)//")")
        end do

        associate(now => sc%now)
        if (allocated(now%z_srf))   call write2D("z_srf",   now%z_srf,   "m",    "Surface elevation")
        if (allocated(now%tsl))     call write3D("tsl",     now%tsl,     "month", "K",    "Sea-level air temperature")
        if (allocated(now%ta_sum))  call write2D("ta_sum",  now%ta_sum,  "K",    "Near-surface air temperature (sum)")
        end associate

        return

    contains

        subroutine write2D(vnm, var, units, long_name)
            character(len=*), intent(IN) :: vnm, units, long_name
            real(wp),         intent(IN) :: var(:,:)
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                          units=units, long_name=long_name)
        end subroutine write2D

        subroutine write3D(vnm, var, dim3, units, long_name)
            character(len=*), intent(IN) :: vnm, dim3, units, long_name
            real(wp),         intent(IN) :: var(:,:,:)
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3=dim3, dim4="time", &
                          start=[1,1,1,n], count=[size(var,1),size(var,2),size(var,3),1], &
                          ncid=ncid, units=units, long_name=long_name)
        end subroutine write3D

    end subroutine snapesm_write_step

    ! =====================================================================
    ! Parameter loading & setup
    ! =====================================================================

    subroutine snapesm_par_load(sc, filename, group, domain, grid_name, idx)
        ! Read the full configuration: top-level knobs, the field registry
        ! (&<group>_field_<name>), the snapshot specs (&<group>_<name>), and the
        ! axes (&<group>_axis_<name>: lo, hi, index), whose indices must be in
        ! the registry idx.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group      ! base group (already domain-suffixed)
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name
        type(indices_class),   intent(IN)    :: idx

        character(len=64) :: names(SNP_MAX)
        character(len=64) :: lo, hi
        character(len=256):: agroup
        integer           :: i, j

        sc%par%group     = group
        sc%par%domain    = domain
        sc%par%grid_name = grid_name

        ! --- top-level --------------------------------------------------------
        sc%par%var_defs  = ""
        sc%par%combine   = "anomaly"
        sc%par%ref_name  = ""
        sc%par%lapse     = 0.0_wp
        sc%par%f_p       = 0.0_wp
        sc%par%f_p_ne    = 1.0_wp
        sc%par%f_stdev   = 0.0_wp
        call nml_read(filename, group, "var_defs",  sc%par%var_defs,  init=.TRUE.)
        call nml_read(filename, group, "combine",   sc%par%combine,   init=.TRUE.)
        call nml_read(filename, group, "ref_name",  sc%par%ref_name,  init=.TRUE.)
        call nml_read(filename, group, "lapse",     sc%par%lapse,     init=.TRUE.)
        call nml_read(filename, group, "f_p",       sc%par%f_p,       init=.TRUE.)
        call nml_read(filename, group, "f_p_ne",    sc%par%f_p_ne,    init=.TRUE.)
        call nml_read(filename, group, "f_stdev",   sc%par%f_stdev,   init=.TRUE.)

        ! --- field registry ---------------------------------------------------
        names = ""
        call nml_read(filename, group, "fields", names, init=.TRUE.)
        sc%par%n_field = count(len_trim(names) > 0)
        if (allocated(sc%registry)) deallocate(sc%registry)
        allocate(sc%registry(sc%par%n_field))
        do i = 1, sc%par%n_field
            call read_field_spec(sc%registry(i), filename, &
                                 trim(group)//"_field_"//trim(names(i)), &
                                 trim(names(i)), sc%par%combine)
        end do

        ! --- snapshots --------------------------------------------------------
        names = ""
        call nml_read(filename, group, "snapshots", names, init=.TRUE.)
        sc%par%n_snap = count(len_trim(names) > 0)
        if (allocated(sc%snap)) deallocate(sc%snap)
        allocate(sc%snap(sc%par%n_snap))
        do i = 1, sc%par%n_snap
            call read_snapshot_spec(sc%snap(i)%spec, filename, &
                                    trim(group)//"_"//trim(names(i)), &
                                    trim(names(i)), sc%par%ref_name)
        end do

        ! --- axes: two snapshots and a registry index -------------------------
        names = ""
        call nml_read(filename, group, "axes", names, init=.TRUE.)
        sc%par%n_axis = count(len_trim(names) > 0)
        if (allocated(sc%axis)) deallocate(sc%axis)
        allocate(sc%axis(sc%par%n_axis))
        do i = 1, sc%par%n_axis
            agroup = trim(group)//"_axis_"//trim(names(i))
            sc%axis(i)%name = trim(names(i))
            call nml_read(filename, trim(agroup), "lo",    lo,                 init=.TRUE.)
            call nml_read(filename, trim(agroup), "hi",    hi,                 init=.TRUE.)
            call nml_read(filename, trim(agroup), "index", sc%axis(i)%index,   init=.TRUE.)
            sc%axis(i)%lo = snapshot_position(sc, lo, "lo of &"//trim(agroup))
            sc%axis(i)%hi = snapshot_position(sc, hi, "hi of &"//trim(agroup))
            if (.not. indices_has(idx, sc%axis(i)%index)) then
                write(*,*) "snapesm:: Error: index '"//trim(sc%axis(i)%index)//"' (index of &"// &
                           trim(agroup)//") is not in [indices] names."
                stop
            end if
            sc%axis(i)%f = 0.0_wp
        end do

        ! --- field axes -> positions in sc%axis -------------------------------
        do i = 1, sc%par%n_field
            do j = 1, sc%registry(i)%n_axis
                sc%registry(i)%axis_k(j) = axis_position(sc, sc%registry(i)%axes(j), &
                                "field "//trim(sc%registry(i)%name))
            end do
        end do

        return
    end subroutine snapesm_par_load

    integer function snapshot_position(sc, name, what) result(s)
        ! Position of the snapshot `name` in sc%snap; stops if it is not listed.
        type(snapesm_class), intent(IN) :: sc
        character(len=*),    intent(IN) :: name, what
        do s = 1, sc%par%n_snap
            if (trim(sc%snap(s)%spec%name) == trim(name)) return
        end do
        write(*,*) "snapesm:: Error: snapshot '"//trim(name)//"' ("//trim(what)// &
                   ") is not in the snapshots list."
        stop
    end function snapshot_position

    integer function axis_position(sc, name, what) result(k)
        ! Position of the axis `name` in sc%axis; stops if it is not listed.
        type(snapesm_class), intent(IN) :: sc
        character(len=*),    intent(IN) :: name, what
        do k = 1, sc%par%n_axis
            if (trim(sc%axis(k)%name) == trim(name)) return
        end do
        write(*,*) "snapesm:: Error: axis '"//trim(name)//"' ("//trim(what)// &
                   ") is not in the axes list."
        stop
    end function axis_position

    subroutine read_field_spec(fs, filename, group, name, combine_default)
        implicit none
        type(field_spec_class), intent(OUT) :: fs
        character(len=*),       intent(IN)  :: filename, group, name, combine_default

        fs%name           = name
        fs%kind           = "atm_temp"
        fs%combine        = combine_default
        fs%blend          = "linear"
        fs%axes           = ""
        fs%n_axis         = 0
        fs%axis_k         = 0
        fs%enabled        = .TRUE.
        fs%apply_lapse    = .FALSE.
        fs%precip_scaling = .FALSE.
        fs%seasonal_synth = .FALSE.
        ! Only the keys the current model uses are read (kind/blend/axes). The
        ! remaining knobs (combine/enabled/*_flags) keep their struct defaults; they
        ! cannot be made optional in the namelist until nml gains an ignore-missing
        ! mode (it errors on any absent key), so they are not read here yet.
        call nml_read(filename, group, "kind",  fs%kind,  init=.TRUE.)
        call nml_read(filename, group, "blend", fs%blend, init=.TRUE.)
        call nml_read(filename, group, "axes",  fs%axes,  init=.TRUE.)
        fs%n_axis = count(len_trim(fs%axes) > 0)

        select case(trim(fs%blend))
            case("linear", "ratio", "const")
                ! ok
            case("fraction")
                write(*,*) "snapesm:: Error: &"//trim(group)//" blend = 'fraction' was removed: &
                           &set blend = 'const' and use the climate modifier [ocn_fraction] &
                           &(f_to, region; see docs/climate-snap.md)."
                stop 1
            case default
                write(*,*) "snapesm:: Error: &"//trim(group)//" blend = '"//trim(fs%blend)// &
                           "'; must be linear, ratio or const."
                stop 1
        end select

        return
    end subroutine read_field_spec

    subroutine read_snapshot_spec(ss, filename, group, name, ref_name)
        implicit none
        type(snapshot_spec_class), intent(OUT) :: ss
        character(len=*),          intent(IN)  :: filename, group, name, ref_name

        ss%name      = name
        ss%is_ref    = .FALSE.
        ss%monthly   = .FALSE.
        call nml_read(filename, group, "monthly",   ss%monthly,   init=.TRUE.)
        ss%clim_sigma = 0.0_wp
        if (nml_has_param(filename, group, "clim_sigma")) &
            call nml_read(filename, group, "clim_sigma", ss%clim_sigma, init=.TRUE.)
        ss%ocn_sigma = 0.0_wp
        if (nml_has_param(filename, group, "ocn_sigma")) &
            call nml_read(filename, group, "ocn_sigma", ss%ocn_sigma, init=.TRUE.)
        ! The reference snapshot is designated solely by the top-level ref_name
        ! (an is_ref namelist key would have to be present in every snapshot group,
        ! given nml's error-on-missing behavior).
        if (trim(name) == trim(ref_name)) ss%is_ref = .TRUE.

        return
    end subroutine read_snapshot_spec

    subroutine snapesm_load_snapshots(sc, filename, time, grid)
        ! For each snapshot and each enabled field, read the varslice group reference(s)
        ! from the state group &<group>_<snapshot> (key = field name; 1 group = already
        ! monthly, 2 groups = [ann, sum] synthesized later), load them from the var_defs
        ! database, prime to `time`, reduce to the sea-level state, and set the reference.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        real(dp),              intent(IN)    :: time
        type(grid_class),      intent(IN)    :: grid     ! target grid (online remap)

        integer            :: s, f, k, nsrc
        character(len=64)  :: refs(2)
        character(len=256) :: sgroup, vdb

        vdb = trim(sc%par%var_defs)

        do s = 1, sc%par%n_snap

            sgroup = trim(sc%par%group)//"_"//trim(sc%snap(s)%spec%name)

            if (allocated(sc%snap(s)%bind)) deallocate(sc%snap(s)%bind)
            allocate(sc%snap(s)%bind(sc%par%n_field))

            do f = 1, sc%par%n_field
                if (.not. sc%registry(f)%enabled) cycle

                ! Group reference(s) for this (snapshot, field): 1 (monthly) or 2 ([ann,sum]).
                refs = ""
                call nml_read(filename, trim(sgroup), trim(sc%registry(f)%name), refs, init=.TRUE.)
                nsrc = count(len_trim(refs) > 0)
                if (nsrc == 0) cycle          ! this snapshot does not supply this field

                allocate(sc%snap(s)%bind(f)%src(nsrc))
                do k = 1, nsrc
                    call varslice_init_nml(sc%snap(s)%bind(f)%src(k), trim(vdb), trim(refs(k)), &
                                           domain=trim(sc%par%domain), grid_name=trim(sc%par%grid_name), &
                                           grid=grid)
                    ! Prime: static read once here, time-varying sliced at `time`.
                    if (sc%snap(s)%bind(f)%src(k)%par%with_time) then
                        call varslice_update(sc%snap(s)%bind(f)%src(k), [time], method="extrap")
                    else
                        call varslice_update(sc%snap(s)%bind(f)%src(k))
                    end if
                end do
            end do

            ! Reduce raw loads -> sea-level state for this snapshot.
            call snapesm_reduce_snapshot(sc, s)

        end do

        ! Set the reference state (anomaly baseline; clim0-equivalent for callers).
        call snapesm_set_ref(sc)

        ! Seed the output state from the reference (allocates now; mirrors snapclim's
        ! `now = clim0`).
        sc%now = sc%ref

        return
    end subroutine snapesm_load_snapshots

    subroutine snapesm_reduce_snapshot(sc, s)
        ! Extract snapshot s's varslice fields into its reduced (sea-level) state.
        ! Faithful port of snapclim read_climate_snapshot (1685-1969) + read_ocean_snapshot.
        ! Monthly snapshots (spec%monthly): raw tas/pr are monthly (pr = sum of sources, sf+rf).
        ! Annual snapshots: tas = cosine-synthesized from [ta_ann, ta_sum]; pr from pr_ann.
        ! Then (both): sea-level lapse reduction, prcor, and ocean vertical interp.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s

        integer  :: nx, ny, m, k, nsrc
        integer  :: f_tas, f_pr, f_zs, f_to, f_so
        logical  :: south
        real(wp) :: zmean, l1, l2

        nx = sc%par%nx
        ny = sc%par%ny
        south = sc%par%south
        l1 = sc%par%lapse(1)
        l2 = sc%par%lapse(2)

        f_zs  = field_bind_index(sc, "zs")
        f_tas = field_bind_index(sc, "tas")
        f_pr  = field_bind_index(sc, "pr")
        f_to  = field_bind_index(sc, "to")
        f_so  = field_bind_index(sc, "so")

        associate(spec => sc%snap(s)%spec, st => sc%snap(s)%state, bind => sc%snap(s)%bind)

        ! --- elevation, mask, beta_p (snapclim 1746-1755) -------------------------
        if (allocated(st%z_srf))  deallocate(st%z_srf)
        if (allocated(st%mask))   deallocate(st%mask)
        if (allocated(st%beta_p)) deallocate(st%beta_p)
        allocate(st%z_srf(nx,ny), st%mask(nx,ny), st%beta_p(nx,ny))
        st%z_srf = bind(f_zs)%src(1)%var(:,:,1,1)
        ! Missing-value fill -> domain mean (GRL inputs are complete; guard on the
        ! -9999 fill sentinel, safe since surface elevation is >= 0).
        if (any(st%z_srf <= -9000.0_wp)) then
            zmean = sum(st%z_srf, mask=st%z_srf > -9000.0_wp) &
                    / real(max(count(st%z_srf > -9000.0_wp),1), wp)
            where (st%z_srf <= -9000.0_wp) st%z_srf = zmean
        end if
        st%mask = 0.0_wp
        where (st%z_srf > 0.0_wp) st%mask = 1.0_wp
        ! beta_p = f_p everywhere (snapclim hardcodes f_p_ne=1.0, so the NE-basin
        ! branch is dormant); f_p_ne is retained in par for future basin dependence.
        st%beta_p = sc%par%f_p * 1.0_wp
        if (sc%par%f_p_ne /= 1.0_wp) st%beta_p = sc%par%f_p  ! (kept uniform; see note)

        ! --- atmosphere: raw monthly tas/pr -------------------------------------
        if (allocated(st%tas))    deallocate(st%tas)
        if (allocated(st%pr))     deallocate(st%pr)
        if (allocated(st%ta_ann)) deallocate(st%ta_ann)
        if (allocated(st%ta_sum)) deallocate(st%ta_sum)
        if (allocated(st%pr_ann)) deallocate(st%pr_ann)
        allocate(st%tas(nx,ny,NMONTH), st%pr(nx,ny,NMONTH))
        allocate(st%ta_ann(nx,ny), st%ta_sum(nx,ny), st%pr_ann(nx,ny))

        if (spec%monthly) then
            ! Monthly climatology (snapclim clim_monthly=True, 1757-1840).
            do m = 1, NMONTH
                st%tas(:,:,m) = bind(f_tas)%src(1)%var(:,:,m,1)
            end do
            ! pr = sum of monthly sources (name3=="sf" -> pr = sf + rf).
            nsrc = size(bind(f_pr)%src)
            st%pr = bind(f_pr)%src(1)%var(:,:,1:NMONTH,1)
            do k = 2, nsrc
                st%pr = st%pr + bind(f_pr)%src(k)%var(:,:,1:NMONTH,1)
            end do
            where (st%pr < 0.0_wp) st%pr = 0.0_wp
            st%ta_ann = sum(st%tas, dim=3) / 12.0_wp
            if (south) then
                st%ta_sum = sum(st%tas(:,:,[12,1,2]), dim=3) / 3.0_wp
            else
                st%ta_sum = sum(st%tas(:,:,[6,7,8]),  dim=3) / 3.0_wp
            end if
            st%pr_ann = sum(st%pr, dim=3) / 12.0_wp * 365.0_wp
        else
            ! Annual+summer inputs, cosine-synthesized to monthly (snapclim 1842-1910).
            st%ta_ann = bind(f_tas)%src(1)%var(:,:,1,1)
            st%ta_sum = bind(f_tas)%src(2)%var(:,:,1,1)
            do m = 1, NMONTH
                if (south) then
                    st%tas(:,:,m) = st%ta_ann + (st%ta_sum-st%ta_ann)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp)
                else
                    st%tas(:,:,m) = st%ta_ann - (st%ta_sum-st%ta_ann)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp)
                end if
            end do
            ! pr_ann raw is [mm/d]; monthly pr = raw (flat), object pr_ann = raw*365 [mm/a].
            st%pr_ann = bind(f_pr)%src(1)%var(:,:,1,1)
            do m = 1, NMONTH
                st%pr(:,:,m) = st%pr_ann
            end do
            st%pr_ann = st%pr_ann * 365.0_wp
        end if

        ! Smooth the snapshot (spec%clim_sigma, e.g. coarse GCM fields), then the
        ! land mask of the smoothed elevation (snapclim read_climate_snapshot).
        if (spec%clim_sigma > 0.0_wp) then
            call filter_gaussian(st%z_srf,  spec%clim_sigma, sc%par%dx_km)
            call filter_gaussian(st%ta_ann, spec%clim_sigma, sc%par%dx_km)
            call filter_gaussian(st%ta_sum, spec%clim_sigma, sc%par%dx_km)
            call filter_gaussian(st%pr_ann, spec%clim_sigma, sc%par%dx_km)
            do m = 1, NMONTH
                call filter_gaussian(st%tas(:,:,m), spec%clim_sigma, sc%par%dx_km)
                call filter_gaussian(st%pr(:,:,m),  spec%clim_sigma, sc%par%dx_km)
            end do
            st%mask = 0.0_wp
            where (st%z_srf > 0.0_wp) st%mask = 1.0_wp
        end if

        ! Precip variability scaling (snapclim 1926-1927); f_stdev=0 -> no-op here.
        ! (pr_stdev_frac == 0 for this config, so the monthly scaling is identity.)

        ! --- sea-level temperatures and elevation-desensitized precip -----------
        if (allocated(st%tsl))       deallocate(st%tsl)
        if (allocated(st%prcor))     deallocate(st%prcor)
        if (allocated(st%tsl_ann))   deallocate(st%tsl_ann)
        if (allocated(st%tsl_sum))   deallocate(st%tsl_sum)
        if (allocated(st%prcor_ann)) deallocate(st%prcor_ann)
        allocate(st%tsl(nx,ny,NMONTH), st%prcor(nx,ny,NMONTH))
        allocate(st%tsl_ann(nx,ny), st%tsl_sum(nx,ny), st%prcor_ann(nx,ny))

        st%tsl_ann   = st%ta_ann + l1*st%z_srf
        st%tsl_sum   = st%ta_sum + l2*st%z_srf
        st%prcor_ann = st%pr_ann / exp(st%beta_p*(st%ta_ann - st%tsl_ann))
        do m = 1, NMONTH
            if (south) then
                st%tsl(:,:,m) = st%tas(:,:,m) + st%z_srf*(l1 + (l2-l1)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp))
            else
                st%tsl(:,:,m) = st%tas(:,:,m) + st%z_srf*(l1 + (l1-l2)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp))
            end if
            st%prcor(:,:,m) = st%pr(:,:,m) / exp(st%beta_p*(st%tas(:,:,m) - st%tsl(:,:,m)))
        end do

        ! --- ocean (only if this snapshot supplies 3-D to/so) -------------------
        if (f_to > 0 .and. f_so > 0) then
            if (allocated(bind(f_to)%src) .and. allocated(bind(f_so)%src)) then
                call reduce_ocean(sc, s, f_to, f_so)
            end if
        end if

        end associate

        return
    end subroutine snapesm_reduce_snapshot

    subroutine reduce_ocean(sc, s, f_to, f_so)
        ! Monthly-average to/so over the 12-month axis, then vertically interpolate
        ! the input depth profiles to the model depth axis (snapclim read_ocean_snapshot
        ! 2078-2150; 42 -> 23 levels, depth0 = abs(file depth)).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s, f_to, f_so

        integer  :: nx, ny, nz0, nzo, i, j, k, nt
        real(dp) :: tlo, thi
        real(wp), allocatable :: depth0(:), depthm(:), to0(:,:,:), so0(:,:,:)

        nx = sc%par%nx
        ny = sc%par%ny

        associate(st => sc%snap(s)%state, vt => sc%snap(s)%bind(f_to)%src(1), &
                                          vs => sc%snap(s)%bind(f_so)%src(1))

        ! Read all 12 months then average to the annual mean (snapclim sum(dim=4)/12).
        ! rep=12 selects the sub-annual axis via whole-year matching (get_indices uses
        ! floor(month)==year), which is robust to the sp/dp rounding of the fractional
        ! month positions that a direct fractional range would trip on.
        tlo = real(floor(minval(vt%time)), dp)   ! the (single) whole year of the axis
        call varslice_update(vt, [tlo,tlo], method="range", rep=12)
        call varslice_update(vs, [tlo,tlo], method="range", rep=12)

        nz0 = size(vt%var,3)
        nt  = size(vt%var,4)
        allocate(depth0(nz0)); depth0 = abs(vt%z(1:nz0))
        allocate(to0(nx,ny,nz0), so0(nx,ny,nz0))
        to0 = sum(vt%var, dim=4) / real(nt,wp)
        so0 = sum(vs%var, dim=4) / real(size(vs%var,4),wp)

        ! Smooth the ocean (spec%ocn_sigma), level by level
        if (sc%snap(s)%spec%ocn_sigma > 0.0_wp) then
            do k = 1, nz0
                call filter_gaussian(to0(:,:,k), sc%snap(s)%spec%ocn_sigma, sc%par%dx_km)
                call filter_gaussian(so0(:,:,k), sc%snap(s)%spec%ocn_sigma, sc%par%dx_km)
            end do
        end if

        call model_depth(depthm)
        nzo = size(depthm)
        if (allocated(st%depth))  deallocate(st%depth)
        if (allocated(st%to_ann)) deallocate(st%to_ann)
        if (allocated(st%so_ann)) deallocate(st%so_ann)
        allocate(st%depth(nzo)); st%depth = depthm
        allocate(st%to_ann(nx,ny,nzo), st%so_ann(nx,ny,nzo))

        do k = 1, nzo
        do j = 1, ny
        do i = 1, nx
            st%to_ann(i,j,k) = interp_linear(depth0, to0(i,j,:), xout=depthm(k))
            st%so_ann(i,j,k) = interp_linear(depth0, so0(i,j,:), xout=depthm(k))
        end do
        end do
        end do

        end associate

        return
    end subroutine reduce_ocean

    subroutine model_depth(depth)
        ! Model ocean depth axis: nzo=23, 0..2000 m in 21 steps + [2500, 3000]
        ! (snapclim snapclim_init 300-307).
        implicit none
        real(wp), allocatable, intent(OUT) :: depth(:)
        integer, parameter :: nzo = 23
        integer :: k
        if (allocated(depth)) deallocate(depth)
        allocate(depth(nzo))
        do k = 1, nzo-2
            depth(k) = (k-1)*2000.0_wp/real(nzo-3,wp)
        end do
        depth(nzo-1:nzo) = [2500.0_wp, 3000.0_wp]
        return
    end subroutine model_depth

    function field_bind_index(sc, name) result(f)
        ! Registry index of the field named `name` (case-sensitive), or 0 if absent.
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: name
        integer :: f, i
        f = 0
        do i = 1, sc%par%n_field
            if (trim(sc%registry(i)%name) == trim(name)) then
                f = i
                return
            end if
        end do
        return
    end function field_bind_index

    subroutine snapesm_set_ref(sc)
        ! Copy the designated reference snapshot's reduced state into sc%ref
        ! (the anomaly baseline). Falls back to the first snapshot if none flagged.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        integer :: s, iref
        iref = 0
        do s = 1, sc%par%n_snap
            if (sc%snap(s)%spec%is_ref) then
                iref = s
                exit
            end if
        end do
        if (iref == 0 .and. sc%par%n_snap >= 1) iref = 1
        if (iref >= 1) sc%ref = sc%snap(iref)%state
        return
    end subroutine snapesm_set_ref

    function interp_linear(x, y, xout) result(yout)
        ! Simple linear point interpolation, clamped to the endpoints. Private copy
        ! (identical to snapclim's) so the ocean vinterp does not depend on which
        ! interp1D module wins the include path in a full yelmox build, and matches
        ! snapclim's ocean interpolation exactly.
        implicit none
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp),               intent(IN) :: xout
        real(wp) :: yout
        integer  :: j, n
        real(wp) :: alph
        n = size(x)
        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n
                if (x(j) .ge. xout) exit
            end do
            if (j .eq. 1) then
                yout = y(1)
            else if (j .eq. n+1) then
                yout = y(n)
            else
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if
        end if
        return
    end function interp_linear

end module snapesm
