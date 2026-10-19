module snapshots
    ! snapshots — climate boundary forcing built on the fesm-utils primitives.
    !
    !   Loading   -> varslice  (one varslice_class per snapshot x field)
    !   Indices   -> the driver's registry (kryos_forcing), read by name
    !   Wrapper   -> this module (blend + transforms + output state)
    !
    ! See docs/climate-snapshots.md for its configuration.
    !
    ! Config model: the `var_defs` databases (a list of nml files, first file
    ! that has a group wins, then the par file itself for run-specific groups;
    ! the shared input/clim_db.nml by default) define varslice groups (permanent per-file variable
    ! metadata); each climate-state group &<group>_<snapshot> names the group(s)
    ! supplying each field (1 = monthly, 2 = [ann, season] -> synthesized). Their
    ! paths take the keywords {grid_global} ([snap] grid_global, the global grid
    ! of the paleoclimate products) and {exp} (the snapshot's exp), besides
    ! {domain} and {grid_name}.
    ! Fields are optional per snapshot (key absent or "" = not bound): the backend
    ! builds the fields the reference snapshot binds, the atmosphere (tas and pr,
    ! with zs) and/or the ocean (to, optionally so); an ocean-only reference serves
    ! the rembo ocean. Every snapshot an axis of a blended field names must bind
    ! that field (checked at init).
    !
    ! Axes: each axis &<group>_axis_<name> spans two snapshots (lo -> hi) and
    ! names an index of the driver's registry ([indices], kryos_forcing); the
    ! value of the index is the position on the axis (0 at lo, 1 at hi). A field
    ! sums the axes it lists on top of the reference snapshot (Banderas et al.,
    ! 2018, for two axes). With lo = "", hi is itself the anomaly (linear) or the
    ! ratio (ratio): an anomaly snapshot binds no zs, its fields are taken as
    ! given, at sea level (no elevation correction of their own). The indices
    ! restart with the registry, so snapshots holds no state of its own.
    !
    ! The ocean temperature as a fraction of the air temperature anomaly is the
    ! climate modifier ocn_fraction (climate_mods), on the field `to` held at
    ! the reference (no axes, or blend = "const").

    use precision, only : wp, sp, dp
    use ncio
    use nml
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use kryos_forcing, only : indices_class, indices_value, indices_has
    use yelmox_nml,    only : VAR_DEFS_MAX, var_defs_find, def_file => yelmox_defaults_file
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
        character(len=16) :: kind           ! from the name: "atm_temp" | "atm_precip" | "ocn" | "ocn_salt" | "elev"
        character(len=16) :: combine        ! "anomaly" | "absolute"
        character(len=16) :: blend          ! "linear" | "ratio" | "const"
        character(len=64) :: axes(AXS_MAX)  ! names of the axes of this field ("" = none)
        integer           :: n_axis         ! number of axes
        integer           :: axis_k(AXS_MAX) ! their positions in snapshots_class%axis
        logical           :: enabled        ! off by default for exotic fields
        logical           :: apply_lapse
        logical           :: precip_scaling
        logical           :: seasonal_synth
    end type

    ! One climate state: PD, LGM, piControl, ... (config side).
    type snapshot_spec_class
        character(len=64) :: name
        logical           :: is_ref         ! the distinguished reference state
        logical           :: monthly        ! raw atm fields are monthly (else annual+season -> synthesized)
        character(len=3)  :: tas_season     ! season of the second tas field: "jja" | "djf"
        real(wp)          :: clim_sigma     ! [km] Gaussian smoothing of the atmosphere at loading (0: none)
        real(wp)          :: ocn_sigma      ! [km] Gaussian smoothing of the ocean at loading (0: none)
        character(len=64) :: exp            ! experiment, the {exp} of its var_defs paths
        character(len=64), allocatable :: refs(:,:)  ! (n_field,2) var_defs groups per field ("" = none)
        integer,           allocatable :: nsrc(:)    ! (n_field) number of groups per field (0 = not bound)
    end type

    ! A climate state. Option A: named members for the hot path, unallocated when a
    ! run does not use them (allocated() is the "in use" guard); `extra` for the tail.
    ! Used both as a snapshot's reduced (sea-level) fields and as the blended output.
    type climate_state_class
        ! atmosphere, monthly (nx,ny,12): tas at surface, tsl at sea level
        real(wp), allocatable :: tas(:,:,:), tsl(:,:,:), pr(:,:,:), prcor(:,:,:)
        ! atmosphere, annual / summer (nx,ny); the summer means of the output only
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
    ! or two (annual mean + season) that `reduce` synthesizes into 12 months. The group
    ! names are references into the var_defs databases.
    type field_binding_class
        type(varslice_class), allocatable :: src(:)   ! size 1 (monthly) or 2 ([ann, season])
    end type

    ! One loaded snapshot: its spec, a per-field source binding (registry order), and
    ! its reduced (sea-level) climate state used by the blend.
    type snapshot_class
        type(snapshot_spec_class)              :: spec
        type(field_binding_class), allocatable :: bind(:)
        type(climate_state_class)              :: state
    end type

    ! Top-level parameters.
    type snapshots_param_class
        character(len=256) :: domain
        character(len=256) :: grid_name
        character(len=256) :: var_defs(VAR_DEFS_MAX) ! varslice databases, searched in order
        character(len=512) :: par_file      ! the par file, searched after var_defs
        character(len=64)  :: grid_global   ! global grid of the paleoclimate products ({grid_global})
        character(len=64)  :: group         ! base namelist group (already domain-suffixed)
        integer            :: nx, ny
        real(wp)           :: dx_km         ! [km] grid spacing (snapshot smoothing)
        logical            :: south         ! southern hemisphere (seasons, lapse rates)
        character(len=16)  :: combine       ! default combine mode
        character(len=64)  :: ref_name      ! which snapshot is the reference
        integer            :: iref          ! its position in snap
        logical            :: with_atm      ! the reference binds the atmosphere (tas, pr)
        logical            :: with_ocn      ! the reference binds the ocean (to)
        real(wp)           :: lapse(2)      ! [annual, summer] lapse rate [K/m]
        ! precip sensitivity parameters
        real(wp)           :: f_p           ! precip-temperature Clausius factor -> beta_p [1/K]
        real(wp)           :: f_p_ne        ! NE-Greenland beta_p multiplier (not active)
        integer            :: n_snap, n_field, n_axis
    end type

    ! One axis: two snapshots and the registry index giving the position on it.
    type axis_class
        character(len=64) :: name
        character(len=64) :: index          ! registry index (kryos_forcing)
        integer           :: lo, hi         ! snapshots at index = 0 and 1 (positions in snap;
                                            ! lo = 0: none, hi is the anomaly/ratio)
        real(wp)          :: f              ! position at the last update (output)
    end type

    ! The object callers hold.
    type snapshots_class
        type(snapshots_param_class)        :: par
        type(field_spec_class),  allocatable :: registry(:)  ! field definitions
        type(snapshot_class),    allocatable :: snap(:)      ! was clim0..clim3
        type(axis_class),        allocatable :: axis(:)      ! snapshot axes
        type(climate_state_class)            :: now          ! current output
        type(climate_state_class)            :: ref          ! distinguished reference (was clim0)
    end type

    public :: snapshots_class
    public :: climate_state_class
    public :: snapshots_init
    public :: snapshots_update
    public :: snapshots_end
    public :: snapshots_write_init
    public :: snapshots_write_step
    public :: snapshots_lapse_monthly

contains

    ! =====================================================================
    ! Lifecycle
    ! =====================================================================

    subroutine snapshots_init(sc, filename, domain, grid, time, south, idx, group)
        ! `grid` is the climate grid: the snapshot fields are read on it, or
        ! remapped onto it online (var_defs groups with remap/grid_src). `time`
        ! primes the time-varying snapshots. `idx` is the driver's registry: each
        ! axis must name one of its indices.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
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

        sc%par%nx = grid%G%nx
        sc%par%ny = grid%G%ny
        sc%par%dx_km = real(grid%G%dx*grid%cs%xy_conv*1e-3_dp, wp)   ! [km] (clim_sigma)
        sc%par%south = south

        ! Read all configuration: top-level, field registry, snapshot specs and
        ! their bindings, axes; then check the bindings against each other.
        call snapshots_par_load(sc, filename, trim(base_group), domain, trim(grid%name), idx)
        call snapshots_check_bindings(sc)

        ! Load the snapshot fields (varslice), reduce them, and set the reference state.
        call snapshots_load_snapshots(sc, time, grid)

        return
    end subroutine snapshots_init

    subroutine snapshots_update(sc, z_srf, time, idx, dx)
        ! Pipeline: axis positions from the registry
        ! (advanced by the driver to `time`) -> refresh loads -> combine ->
        ! transform (the atmosphere, if bound). Homogeneous anomalies and the
        ! ocean fraction rule are climate modifiers (climate_mods), applied by
        ! the caller afterwards.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
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
        call snapshots_refresh_loads(sc, time)

        ! 3. Combine the reduced snapshot states into `now`, per field.
        call snapshots_combine(sc)

        ! 4. Transform: inflate sea-level fields to the current model elevation
        !    (tsl -> tas), precip, annual/summer aggregates.
        if (sc%par%with_atm) call snapshots_transform(sc, z_srf)

        return
    end subroutine snapshots_update

    ! ---------------------------------------------------------------------
    ! Update pipeline stages
    ! ---------------------------------------------------------------------

    subroutine snapshots_refresh_loads(sc, time)
        ! Re-slice only the time-varying fields (recon, transient orography, etc.);
        ! static fields were primed once at load. A snapshot with a refreshed
        ! source is reduced again (and the reference state reset, if it is the
        ! reference).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        real(dp),              intent(IN)    :: time
        integer :: s, f, k
        logical :: touched
        do s = 1, sc%par%n_snap
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
            if (touched) then
                call snapshots_reduce_snapshot(sc, s)
                if (s == sc%par%iref) sc%ref = sc%snap(s)%state
            end if
        end do
        return
    end subroutine snapshots_refresh_loads

    subroutine snapshots_combine(sc)
        ! Blend the reduced snapshot states into `now`, per field bound in the
        ! reference, from the reference baseline (was clim0). Each axis k of a
        ! field (lo_k -> hi_k) adds with weight a_k = the value of its index:
        !   linear (temp/salinity): now = ref + sum_k a_k*(hi_k - lo_k)
        !   ratio  (precip):        now = ref * prod_k (a_k*(hi_k/lo_k - 1) + 1)
        !   const:                  now = ref
        ! An axis with lo = "" adds hi itself: a_k*hi_k (linear) or the factor
        ! a_k*(hi_k - 1) + 1 (ratio). Two axes as in Banderas et al. (2018).
        ! The snapshots of the axes bind the field
        ! (snapshots_check_bindings).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        integer  :: f, j, k, lo, hi
        logical  :: blended
        type(climate_state_class) :: none   ! the lo of a lo = "" axis

        do f = 1, sc%par%n_field
            if (.not. sc%registry(f)%enabled) cycle
            if (sc%snap(sc%par%iref)%spec%nsrc(f) == 0) cycle   ! not bound in the reference: not built

            ! Start from the reference state.
            select case(trim(sc%registry(f)%kind))
                case("atm_temp")
                    sc%now%tsl    = sc%ref%tsl
                case("atm_precip")
                    sc%now%prcor  = sc%ref%prcor
                case("ocn")
                    sc%now%to_ann = sc%ref%to_ann
                case("ocn_salt")
                    sc%now%so_ann = sc%ref%so_ann
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
                if (lo > 0) then
                    call add_axis(sc%registry(f), sc%now, sc%snap(lo)%state, sc%snap(hi)%state, &
                                  sc%axis(k)%f, anomaly=.FALSE.)
                else
                    call add_axis(sc%registry(f), sc%now, none, sc%snap(hi)%state, &
                                  sc%axis(k)%f, anomaly=.TRUE.)
                end if
            end do
        end do

        return
    end subroutine snapshots_combine

    subroutine add_axis(fs, now, x0, x1, aa, anomaly)
        ! Add one axis (x0 -> x1, weight aa; anomaly: x1 alone) to the field fs
        ! of `now`.
        implicit none
        type(field_spec_class),    intent(IN)    :: fs
        type(climate_state_class), intent(INOUT) :: now
        type(climate_state_class), intent(IN)    :: x0, x1
        real(wp),                  intent(IN)    :: aa
        logical,                   intent(IN)    :: anomaly

        select case(trim(fs%kind))
            case("atm_temp")
                call blend_axis(fs%blend, now%tsl,    x0%tsl,    x1%tsl,    aa, anomaly)
            case("atm_precip")
                call blend_axis(fs%blend, now%prcor,  x0%prcor,  x1%prcor,  aa, anomaly)
            case("ocn")
                call blend_axis(fs%blend, now%to_ann, x0%to_ann, x1%to_ann, aa, anomaly)
            case("ocn_salt")
                call blend_axis(fs%blend, now%so_ann, x0%so_ann, x1%so_ann, aa, anomaly)
        end select

        return
    end subroutine add_axis

    subroutine blend_axis(mode, dst, x0, x1, aa, anomaly)
        ! Add one axis (x0 -> x1, weight aa) to dst, which holds the reference
        ! plus the axes added so far. anomaly (lo = ""): x1 is the anomaly
        ! (linear) or the ratio (ratio) itself, x0 is not used.
        implicit none
        character(len=*),      intent(IN)    :: mode
        real(wp), allocatable, intent(INOUT) :: dst(:,:,:)
        real(wp), allocatable, intent(IN)    :: x0(:,:,:), x1(:,:,:)
        real(wp),              intent(IN)    :: aa
        logical,               intent(IN)    :: anomaly

        if (anomaly) then
            select case(trim(mode))
                case("ratio")
                    dst = dst*(aa*(x1 - 1.0_wp) + 1.0_wp)
                case default   ! "linear"
                    dst = dst + aa*x1
            end select
        else
            select case(trim(mode))
                case("ratio")
                    where (x0 /= 0.0_wp) dst = dst*(aa*(x1/x0 - 1.0_wp) + 1.0_wp)
                case default   ! "linear"
                    dst = dst + aa*(x1 - x0)
            end select
        end if

        return
    end subroutine blend_axis

    subroutine snapshots_transform(sc, z_srf)
        ! Inflate the blended sea-level fields to the current model elevation and
        ! compute the annual/summer aggregates. Note the cosine
        ! phase here is m*30-15 (vs m*30-30 in the sea-level reduction).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        integer  :: m
        logical  :: south
        real(wp) :: lapse(NMONTH)

        south = sc%par%south
        lapse = snapshots_lapse_monthly(sc)

        sc%now%z_srf = z_srf
        sc%now%mask  = 0.0_wp
        where (z_srf > 1.0_wp) sc%now%mask = 1.0_wp

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
    end subroutine snapshots_transform

    function snapshots_lapse_monthly(sc) result(lapse)
        ! Monthly lapse rate [K/m] that lifts tsl to the surface:
        ! annual lapse(1) with a seasonal cosine of amplitude lapse(1)-lapse(2),
        ! phase m*30-15, summer at lapse(2).
        implicit none
        type(snapshots_class), intent(IN) :: sc
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
    end function snapshots_lapse_monthly

    subroutine snapshots_end(sc)
        ! Teardown: the varslice sources of every snapshot, the snapshots, axes
        ! and registry, and the output and reference states.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc

        type(climate_state_class) :: empty
        integer :: s, f, k

        do s = 1, sc%par%n_snap
            do f = 1, sc%par%n_field
                if (.not. allocated(sc%snap(s)%bind(f)%src)) cycle   ! not bound
                do k = 1, size(sc%snap(s)%bind(f)%src)
                    call varslice_end(sc%snap(s)%bind(f)%src(k))
                end do
            end do
        end do
        deallocate(sc%snap, sc%axis, sc%registry)

        sc%now = empty
        sc%ref = empty

        return
    end subroutine snapshots_end

    ! =====================================================================
    ! Diagnostic output
    ! =====================================================================

    subroutine snapshots_write_init(sc, filename)
        ! Add the month axis of the snapshots fields to an existing 2D output file
        ! (xc, yc, time and depth defined by the caller).
        implicit none
        type(snapshots_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename

        call nc_write_dim(filename, "month", x=1, dx=1, nx=NMONTH, units="month")

        return
    end subroutine snapshots_write_init

    subroutine snapshots_write_step(sc, filename, ncid, n)
        ! Write record n of the snapshots state into an open 2D output file
        ! (axes from snapshots_write_init): the axis positions and the surface
        ! elevation. The climate as applied (after the climate modifiers) is
        ! written by the caller (climate_write_2D).
        implicit none
        type(snapshots_class), intent(IN) :: sc
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

    end subroutine snapshots_write_step

    ! =====================================================================
    ! Parameter loading & setup
    ! =====================================================================

    subroutine snapshots_par_load(sc, filename, group, domain, grid_name, idx)
        ! Read the full configuration: top-level knobs, the field registry
        ! (&<group>_field_<name>), the snapshot specs and their field bindings
        ! (&<group>_<name>), and the axes (&<group>_axis_<name>: lo, hi, index),
        ! whose indices must be in the registry idx. Keys not set take the
        ! defaults of input/yelmox_defaults.nml, the named groups those of the
        ! templates snap_field, snap_snapshot and snap_axis.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
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
        sc%par%par_file  = filename
        sc%par%domain    = domain
        sc%par%grid_name = grid_name

        ! --- top-level --------------------------------------------------------
        call nml_validate(filename, def_file, group, defaults_group="snap")
        call nml_read(filename, group, "var_defs",    sc%par%var_defs,    init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        call nml_read(filename, group, "grid_global", sc%par%grid_global, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        call nml_read(filename, group, "combine",     sc%par%combine,     init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        call nml_read(filename, group, "ref_name",    sc%par%ref_name,    init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        call nml_read(filename, group, "lapse",       sc%par%lapse,       init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        call nml_read(filename, group, "f_p",         sc%par%f_p,         init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        call nml_read(filename, group, "f_p_ne",      sc%par%f_p_ne,      init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")

        ! --- field registry ---------------------------------------------------
        call nml_read(filename, group, "fields", names, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        sc%par%n_field = count(len_trim(names) > 0)
        if (allocated(sc%registry)) deallocate(sc%registry)
        allocate(sc%registry(sc%par%n_field))
        do i = 1, sc%par%n_field
            call read_field_spec(sc%registry(i), filename, &
                                 trim(group)//"_field_"//trim(names(i)), &
                                 trim(names(i)), sc%par%combine)
        end do

        ! --- snapshots and their bindings ------------------------------------
        call nml_read(filename, group, "snapshots", names, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        sc%par%n_snap = count(len_trim(names) > 0)
        if (allocated(sc%snap)) deallocate(sc%snap)
        allocate(sc%snap(sc%par%n_snap))
        do i = 1, sc%par%n_snap
            call read_snapshot_spec(sc%snap(i)%spec, filename, &
                                    trim(group)//"_"//trim(names(i)), &
                                    trim(names(i)), sc%par%ref_name, sc%registry, sc%par%south)
        end do
        sc%par%iref = snapshot_position(sc, sc%par%ref_name, "ref_name of &"//trim(group))

        ! --- axes: two snapshots (lo may be "") and a registry index ----------
        call nml_read(filename, group, "axes", names, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap")
        sc%par%n_axis = count(len_trim(names) > 0)
        if (allocated(sc%axis)) deallocate(sc%axis)
        allocate(sc%axis(sc%par%n_axis))
        do i = 1, sc%par%n_axis
            agroup = trim(group)//"_axis_"//trim(names(i))
            sc%axis(i)%name = trim(names(i))
            call nml_validate(filename, def_file, trim(agroup), defaults_group="snap_axis")
            call nml_read(filename, trim(agroup), "lo",    lo,               init=.TRUE., &
                          defaults_file=def_file, defaults_group="snap_axis")
            call nml_read(filename, trim(agroup), "hi",    hi,               init=.TRUE., &
                          defaults_file=def_file, defaults_group="snap_axis")
            call nml_read(filename, trim(agroup), "index", sc%axis(i)%index, init=.TRUE., &
                          defaults_file=def_file, defaults_group="snap_axis")
            sc%axis(i)%lo = 0
            if (len_trim(lo) > 0) sc%axis(i)%lo = snapshot_position(sc, lo, "lo of &"//trim(agroup))
            sc%axis(i)%hi = snapshot_position(sc, hi, "hi of &"//trim(agroup))
            if (.not. indices_has(idx, sc%axis(i)%index)) then
                write(*,*) "snapshots:: Error: index '"//trim(sc%axis(i)%index)//"' (index of &"// &
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

        ! --- what the backend builds: the fields of the reference -------------
        sc%par%with_atm = binds(sc, sc%par%iref, "tas")
        sc%par%with_ocn = binds(sc, sc%par%iref, "to")

        return
    end subroutine snapshots_par_load

    subroutine snapshots_check_bindings(sc)
        ! Check the field bindings of the snapshots against each other (stop
        ! with an error on the first mismatch):
        !  - number of var_defs groups: tas 1 (monthly) or 2 ([ann, season]),
        !    pr 1 or 2 (monthly: sum, e.g. snowfall + rainfall) or 1 (annual),
        !    every other field 1;
        !  - a state snapshot (the reference, or on an axis with a lo) binding tas
        !    or pr binds zs (sea-level reduction); a snapshot binding zs and pr
        !    binds tas (the precipitation is reduced with the temperature). Only
        !    the hi of a lo = "" axis may leave out zs: its fields are then taken
        !    as given, at sea level;
        !  - so needs to;
        !  - the reference binds the atmosphere as tas and pr together, or neither;
        !  - no snapshot binds a field the reference does not (it would not be
        !    built);
        !  - every snapshot of an axis of a blended (linear, ratio) field binds it.
        implicit none
        type(snapshots_class), intent(IN) :: sc

        integer :: s, f, j, k, n, iref
        logical :: state
        character(len=256) :: sgroup

        iref = sc%par%iref

        do s = 1, sc%par%n_snap
            associate(spec => sc%snap(s)%spec)
            sgroup = "&"//trim(sc%par%group)//"_"//trim(spec%name)

            do f = 1, sc%par%n_field
                n = spec%nsrc(f)
                if (n == 0) cycle
                select case(trim(sc%registry(f)%name))
                    case("tas")
                        if ((spec%monthly .and. n /= 1) .or. (.not. spec%monthly .and. n /= 2)) &
                            call binding_error(trim(sgroup)//" tas: one group (monthly = True) &
                                               &or two, annual and season (monthly = False).")
                    case("pr")
                        if (.not. spec%monthly .and. n /= 1) &
                            call binding_error(trim(sgroup)//" pr: one group (annual) with &
                                               &monthly = False.")
                    case default
                        if (n /= 1) &
                            call binding_error(trim(sgroup)//" "//trim(sc%registry(f)%name)// &
                                               ": one group.")
                end select
                if (sc%snap(iref)%spec%nsrc(f) == 0) &
                    call binding_error(trim(sgroup)//" binds "//trim(sc%registry(f)%name)// &
                                       ", the reference "//trim(sc%par%ref_name)//" does not; &
                                       &the backend builds only what the reference binds.")
            end do

            ! A state: the reference, or a snapshot of an axis with a lo.
            state = (s == iref)
            do k = 1, sc%par%n_axis
                if (sc%axis(k)%lo > 0 .and. (sc%axis(k)%lo == s .or. sc%axis(k)%hi == s)) &
                    state = .TRUE.
            end do
            if (state .and. (binds(sc, s, "tas") .or. binds(sc, s, "pr")) .and. &
                .not. binds(sc, s, "zs")) &
                call binding_error(trim(sgroup)//" binds tas or pr but not zs (its elevation, &
                                   &for the sea-level reduction); only the hi of an axis with &
                                   &lo = """" may leave it out.")
            if (binds(sc, s, "zs") .and. binds(sc, s, "pr") .and. .not. binds(sc, s, "tas")) &
                call binding_error(trim(sgroup)//" binds pr and zs but not tas (the &
                                   &precipitation is reduced to sea level with the temperature).")
            if (binds(sc, s, "so") .and. .not. binds(sc, s, "to")) &
                call binding_error(trim(sgroup)//" binds so but not to.")

            end associate
        end do

        if (binds(sc, iref, "tas") .neqv. binds(sc, iref, "pr")) &
            call binding_error("the reference &"//trim(sc%par%group)//"_"//trim(sc%par%ref_name)// &
                               " binds one of tas and pr: the atmosphere needs both (or neither).")

        ! The snapshots of the axes of every blended field the reference binds.
        do f = 1, sc%par%n_field
            if (sc%snap(iref)%spec%nsrc(f) == 0) cycle
            if (trim(sc%registry(f)%blend) /= "linear" .and. trim(sc%registry(f)%blend) /= "ratio") cycle
            do j = 1, sc%registry(f)%n_axis
                k = sc%registry(f)%axis_k(j)
                if (sc%axis(k)%lo > 0) call check_axis_snapshot(sc%axis(k)%lo, "lo")
                call check_axis_snapshot(sc%axis(k)%hi, "hi")
            end do
        end do

        return

    contains

        subroutine check_axis_snapshot(s, end_name)
            integer,          intent(IN) :: s
            character(len=*), intent(IN) :: end_name
            if (sc%snap(s)%spec%nsrc(f) == 0) &
                call binding_error("field "//trim(sc%registry(f)%name)//" (blend "// &
                                   trim(sc%registry(f)%blend)//") lists axis "// &
                                   trim(sc%axis(k)%name)//", whose "//end_name//" snapshot "// &
                                   trim(sc%snap(s)%spec%name)//" does not bind it.")
        end subroutine check_axis_snapshot

        subroutine binding_error(msg)
            character(len=*), intent(IN) :: msg
            write(*,*) "snapshots:: Error: "//trim(msg)
            stop 1
        end subroutine binding_error

    end subroutine snapshots_check_bindings

    integer function snapshot_position(sc, name, what) result(s)
        ! Position of the snapshot `name` in sc%snap; stops if it is not listed.
        type(snapshots_class), intent(IN) :: sc
        character(len=*),    intent(IN) :: name, what
        do s = 1, sc%par%n_snap
            if (trim(sc%snap(s)%spec%name) == trim(name)) return
        end do
        write(*,*) "snapshots:: Error: snapshot '"//trim(name)//"' ("//trim(what)// &
                   ") is not in the snapshots list."
        stop
    end function snapshot_position

    integer function axis_position(sc, name, what) result(k)
        ! Position of the axis `name` in sc%axis; stops if it is not listed.
        type(snapshots_class), intent(IN) :: sc
        character(len=*),    intent(IN) :: name, what
        do k = 1, sc%par%n_axis
            if (trim(sc%axis(k)%name) == trim(name)) return
        end do
        write(*,*) "snapshots:: Error: axis '"//trim(name)//"' ("//trim(what)// &
                   ") is not in the axes list."
        stop
    end function axis_position

    subroutine read_field_spec(fs, filename, group, name, combine_default)
        implicit none
        type(field_spec_class), intent(OUT) :: fs
        character(len=*),       intent(IN)  :: filename, group, name, combine_default

        fs%name           = name
        fs%combine        = combine_default
        fs%n_axis         = 0
        fs%axis_k         = 0
        fs%enabled        = .TRUE.
        fs%apply_lapse    = .FALSE.
        fs%precip_scaling = .FALSE.
        fs%seasonal_synth = .FALSE.
        ! Only the keys the current model uses are read (blend/axes, template
        ! snap_field); the remaining knobs (combine/enabled/*_flags) keep their
        ! struct defaults. The field name decides its kind.
        call nml_validate(filename, def_file, group, defaults_group="snap_field")
        call nml_read(filename, group, "blend", fs%blend, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_field")
        call nml_read(filename, group, "axes",  fs%axes,  init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_field")
        fs%n_axis = count(len_trim(fs%axes) > 0)

        select case(trim(name))
            case("tas")
                fs%kind = "atm_temp"
            case("pr")
                fs%kind = "atm_precip"
            case("to")
                fs%kind = "ocn"
            case("so")
                fs%kind = "ocn_salt"
            case("zs")
                fs%kind = "elev"
            case default
                write(*,*) "snapshots:: Error: &"//trim(group)//": unknown field '"// &
                           trim(name)//"' (tas, pr, to, so, zs)."
                stop 1
        end select

        select case(trim(fs%blend))
            case("linear", "ratio", "const")
                ! ok
            case("fraction")
                write(*,*) "snapshots:: Error: &"//trim(group)//" blend = 'fraction' was removed: &
                           &set blend = 'const' and use the climate modifier [ocn_fraction] &
                           &(f_to, region; see docs/climate-snapshots.md)."
                stop 1
            case default
                write(*,*) "snapshots:: Error: &"//trim(group)//" blend = '"//trim(fs%blend)// &
                           "'; must be linear, ratio or const."
                stop 1
        end select

        return
    end subroutine read_field_spec

    subroutine read_snapshot_spec(ss, filename, group, name, ref_name, registry, south)
        ! The spec of one snapshot and its bindings: per registry field, the
        ! var_defs group(s) under the field's name (key absent or "" = not bound).
        ! tas_season: the season of the second tas field of an annual snapshot,
        ! "jja" or "djf" ("", the default: the summer of the hemisphere). Unset
        ! keys from the template snap_snapshot of the defaults file.
        implicit none
        type(snapshot_spec_class), intent(OUT) :: ss
        character(len=*),          intent(IN)  :: filename, group, name, ref_name
        type(field_spec_class),    intent(IN)  :: registry(:)
        logical,                   intent(IN)  :: south

        character(len=64) :: refs(2)
        integer :: f

        ss%name      = name
        ss%is_ref    = .FALSE.
        call nml_validate(filename, def_file, group, defaults_group="snap_snapshot")
        call nml_read(filename, group, "exp",        ss%exp,        init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_snapshot")
        call nml_read(filename, group, "monthly",    ss%monthly,    init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_snapshot")
        call nml_read(filename, group, "clim_sigma", ss%clim_sigma, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_snapshot")
        call nml_read(filename, group, "ocn_sigma",  ss%ocn_sigma,  init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_snapshot")
        call nml_read(filename, group, "tas_season", ss%tas_season, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snap_snapshot")
        if (len_trim(ss%tas_season) == 0) then
            ss%tas_season = "jja"
            if (south) ss%tas_season = "djf"
        end if
        select case(trim(ss%tas_season))
            case("jja", "djf")
                ! ok
            case default
                write(*,*) "snapshots:: Error: &"//trim(group)//" tas_season = '"// &
                           trim(ss%tas_season)//"'; must be jja or djf."
                stop 1
        end select
        ! The reference snapshot is designated solely by the top-level ref_name.
        if (trim(name) == trim(ref_name)) ss%is_ref = .TRUE.

        ! Bindings: 1 group (monthly) or 2 ([ann, season]; pr monthly: summed).
        allocate(ss%refs(size(registry),2), ss%nsrc(size(registry)))
        ss%refs = ""
        do f = 1, size(registry)
            call nml_read(filename, group, trim(registry(f)%name), refs, init=.TRUE., &
                          defaults_file=def_file, defaults_group="snap_snapshot")
            ss%refs(f,:) = refs
            ss%nsrc(f)   = count(len_trim(refs) > 0)
            if (len_trim(refs(1)) == 0 .and. ss%nsrc(f) > 0) then
                write(*,*) "snapshots:: Error: &"//trim(group)//" "//trim(registry(f)%name)// &
                           ": the first group is empty."
                stop 1
            end if
        end do

        return
    end subroutine read_snapshot_spec

    subroutine snapshots_load_snapshots(sc, time, grid)
        ! For each snapshot, load the var_defs groups of every field it binds
        ! (each from the first database of var_defs that has it, else the par
        ! file, with the path
        ! keywords {grid_global} and the snapshot's {exp}), prime them to
        ! `time`, reduce the snapshot to its sea-level state, and set the
        ! reference.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        real(dp),              intent(IN)    :: time
        type(grid_class),      intent(IN)    :: grid     ! target grid (online remap)

        integer            :: s, f, k, nsrc
        character(len=512) :: vdb
        character(len=256) :: subs(2,2)   ! {key} -> value path substitutions

        do s = 1, sc%par%n_snap

            subs(1,:) = [character(len=256) :: "grid_global", sc%par%grid_global]
            subs(2,:) = [character(len=256) :: "exp",         sc%snap(s)%spec%exp]

            if (allocated(sc%snap(s)%bind)) deallocate(sc%snap(s)%bind)
            allocate(sc%snap(s)%bind(sc%par%n_field))

            do f = 1, sc%par%n_field
                if (.not. sc%registry(f)%enabled) cycle

                nsrc = sc%snap(s)%spec%nsrc(f)
                if (nsrc == 0) cycle          ! this snapshot does not bind this field

                allocate(sc%snap(s)%bind(f)%src(nsrc))
                do k = 1, nsrc
                    vdb = var_defs_find(sc%par%var_defs, trim(sc%par%par_file), &
                                        trim(sc%snap(s)%spec%refs(f,k)), &
                                        "&"//trim(sc%par%group)//"_"//trim(sc%snap(s)%spec%name)// &
                                        " "//trim(sc%registry(f)%name))
                    call varslice_init_nml(sc%snap(s)%bind(f)%src(k), trim(vdb), &
                                           trim(sc%snap(s)%spec%refs(f,k)), &
                                           domain=trim(sc%par%domain), grid_name=trim(sc%par%grid_name), &
                                           subs=subs, grid=grid)
                    ! Prime: static read once here, time-varying sliced at `time`.
                    if (sc%snap(s)%bind(f)%src(k)%par%with_time) then
                        call varslice_update(sc%snap(s)%bind(f)%src(k), [time], method="extrap")
                    else
                        call varslice_update(sc%snap(s)%bind(f)%src(k))
                    end if
                end do
            end do

            ! Reduce raw loads -> sea-level state for this snapshot.
            call snapshots_reduce_snapshot(sc, s)

        end do

        ! The reference state (anomaly baseline; clim0-equivalent for callers),
        ! and the output state seeded from it.
        sc%ref = sc%snap(sc%par%iref)%state
        sc%now = sc%ref

        return
    end subroutine snapshots_load_snapshots

    subroutine snapshots_reduce_snapshot(sc, s)
        ! Extract snapshot s's varslice fields into its reduced (sea-level) state:
        ! the atmosphere if it binds tas or pr, the ocean if it binds to.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s

        if (binds(sc, s, "tas") .or. binds(sc, s, "pr")) call reduce_atmosphere(sc, s)
        if (binds(sc, s, "to")) call reduce_ocean(sc, s)

        return
    end subroutine snapshots_reduce_snapshot

    subroutine reduce_atmosphere(sc, s)
        ! The atmosphere of snapshot s, reduced to sea level.
        ! Monthly snapshots (spec%monthly): raw tas/pr are monthly (pr = sum of
        ! sources, sf+rf). Annual snapshots: tas = cosine-synthesized from
        ! [ta_ann, ta_seas] (tas_season); pr from pr_ann.
        ! A snapshot with its elevation zs is then reduced to sea level (lapse
        ! rates, beta_p); one without (the anomaly snapshot of a lo = "" axis) is
        ! taken as given, a sea-level anomaly (tas) or ratio (pr).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s

        integer  :: nx, ny, m, k, nsrc
        integer  :: f_tas, f_pr, f_zs
        logical  :: south, has_tas, has_pr, has_zs
        real(wp) :: zmean, l1, l2
        real(wp), allocatable :: ta_seas(:,:)

        nx = sc%par%nx
        ny = sc%par%ny
        south = sc%par%south
        l1 = sc%par%lapse(1)
        l2 = sc%par%lapse(2)

        f_zs  = field_bind_index(sc, "zs")
        f_tas = field_bind_index(sc, "tas")
        f_pr  = field_bind_index(sc, "pr")
        has_tas = binds(sc, s, "tas")
        has_pr  = binds(sc, s, "pr")
        has_zs  = binds(sc, s, "zs")

        associate(spec => sc%snap(s)%spec, st => sc%snap(s)%state, bind => sc%snap(s)%bind)

        ! --- elevation, mask, beta_p ----------------------------------------------
        if (has_zs) then
            call realloc2(st%z_srf, nx, ny)
            call realloc2(st%mask,  nx, ny)
            call realloc2(st%beta_p, nx, ny)
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
            ! beta_p = f_p everywhere (the NE-basin branch is dormant); f_p_ne is
            ! retained in par for future basin dependence.
            st%beta_p = sc%par%f_p * 1.0_wp
            if (sc%par%f_p_ne /= 1.0_wp) st%beta_p = sc%par%f_p  ! (kept uniform; see note)
        end if

        ! --- air temperature: monthly tas -----------------------------------------
        if (has_tas) then
            call realloc3(st%tas, nx, ny, NMONTH)
            call realloc2(st%ta_ann, nx, ny)
            if (spec%monthly) then
                ! Monthly climatology.
                do m = 1, NMONTH
                    st%tas(:,:,m) = bind(f_tas)%src(1)%var(:,:,m,1)
                end do
                st%ta_ann = sum(st%tas, dim=3) / 12.0_wp
            else
                ! Annual mean + one season, cosine-synthesized to monthly. With c(m) = cos(2*pi*(m-1)/12), +1 in January and -1
                ! in July, the season's field is put at its extreme month:
                !   djf: tas(m) = ta_ann + (ta_seas - ta_ann)*c(m)   (January)
                !   jja: tas(m) = ta_ann - (ta_seas - ta_ann)*c(m)   (July)
                ! The phase follows the calendar, so the summer field of either
                ! hemisphere (jja north, djf south: the default) and a winter field
                ! give the same kind of cycle, the sign of ta_seas - ta_ann setting
                ! which extreme is the maximum.
                st%ta_ann = bind(f_tas)%src(1)%var(:,:,1,1)
                ta_seas   = bind(f_tas)%src(2)%var(:,:,1,1)
                do m = 1, NMONTH
                    if (trim(spec%tas_season) == "djf") then
                        st%tas(:,:,m) = st%ta_ann + (ta_seas-st%ta_ann)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp)
                    else
                        st%tas(:,:,m) = st%ta_ann - (ta_seas-st%ta_ann)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp)
                    end if
                end do
            end if
        end if

        ! --- precipitation: monthly pr --------------------------------------------
        if (has_pr) then
            call realloc3(st%pr, nx, ny, NMONTH)
            call realloc2(st%pr_ann, nx, ny)
            if (spec%monthly) then
                ! pr = sum of monthly sources (name3=="sf" -> pr = sf + rf).
                nsrc = size(bind(f_pr)%src)
                st%pr = bind(f_pr)%src(1)%var(:,:,1:NMONTH,1)
                do k = 2, nsrc
                    st%pr = st%pr + bind(f_pr)%src(k)%var(:,:,1:NMONTH,1)
                end do
                where (st%pr < 0.0_wp) st%pr = 0.0_wp
                st%pr_ann = sum(st%pr, dim=3) / 12.0_wp * 365.0_wp
            else
                ! pr_ann raw is [mm/d]; monthly pr = raw (flat), object pr_ann = raw*365 [mm/a].
                st%pr_ann = bind(f_pr)%src(1)%var(:,:,1,1)
                do m = 1, NMONTH
                    st%pr(:,:,m) = st%pr_ann
                end do
                st%pr_ann = st%pr_ann * 365.0_wp
            end if
        end if

        ! Smooth the snapshot (spec%clim_sigma, e.g. coarse GCM fields), then the
        ! land mask of the smoothed elevation.
        if (spec%clim_sigma > 0.0_wp) then
            if (has_zs) call filter_gaussian(st%z_srf, spec%clim_sigma, sc%par%dx_km)
            if (has_tas) then
                call filter_gaussian(st%ta_ann, spec%clim_sigma, sc%par%dx_km)
                do m = 1, NMONTH
                    call filter_gaussian(st%tas(:,:,m), spec%clim_sigma, sc%par%dx_km)
                end do
            end if
            if (has_pr) then
                call filter_gaussian(st%pr_ann, spec%clim_sigma, sc%par%dx_km)
                do m = 1, NMONTH
                    call filter_gaussian(st%pr(:,:,m), spec%clim_sigma, sc%par%dx_km)
                end do
            end if
            if (has_zs) then
                st%mask = 0.0_wp
                where (st%z_srf > 0.0_wp) st%mask = 1.0_wp
            end if
        end if

        ! --- sea-level temperatures and elevation-desensitized precip -----------
        if (.not. has_zs) then
            ! Taken as given: a sea-level anomaly or ratio.
            if (has_tas) then
                st%tsl     = st%tas
                st%tsl_ann = st%ta_ann
            end if
            if (has_pr) then
                st%prcor     = st%pr
                st%prcor_ann = st%pr_ann
            end if
        else
            if (has_tas) then
                call realloc3(st%tsl, nx, ny, NMONTH)
                call realloc2(st%tsl_ann, nx, ny)
                st%tsl_ann = st%ta_ann + l1*st%z_srf
                do m = 1, NMONTH
                    if (south) then
                        st%tsl(:,:,m) = st%tas(:,:,m) + st%z_srf*(l1 + (l2-l1)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp))
                    else
                        st%tsl(:,:,m) = st%tas(:,:,m) + st%z_srf*(l1 + (l1-l2)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp))
                    end if
                end do
            end if
            if (has_pr) then   ! with tas and zs (snapshots_check_bindings)
                call realloc3(st%prcor, nx, ny, NMONTH)
                call realloc2(st%prcor_ann, nx, ny)
                st%prcor_ann = st%pr_ann / exp(st%beta_p*(st%ta_ann - st%tsl_ann))
                do m = 1, NMONTH
                    st%prcor(:,:,m) = st%pr(:,:,m) / exp(st%beta_p*(st%tas(:,:,m) - st%tsl(:,:,m)))
                end do
            end if
        end if

        end associate

        return
    end subroutine reduce_atmosphere

    subroutine reduce_ocean(sc, s)
        ! The ocean of snapshot s: the annual mean of to (and so, if bound; else
        ! salinity 0), vertically interpolated from the input depth profiles to
        ! the model depth axis (e.g. 75 -> 23 levels, depth0 = abs(file depth)).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s

        integer  :: nx, ny, nz0, nzo, i, j, k, f_to, f_so
        logical  :: has_so
        real(wp), allocatable :: depth0(:), depthm(:), to0(:,:,:), so0(:,:,:)
        character(len=256) :: what

        nx = sc%par%nx
        ny = sc%par%ny
        f_to   = field_bind_index(sc, "to")
        f_so   = field_bind_index(sc, "so")
        has_so = binds(sc, s, "so")
        what   = "&"//trim(sc%par%group)//"_"//trim(sc%snap(s)%spec%name)

        associate(st => sc%snap(s)%state, vt => sc%snap(s)%bind(f_to)%src(1))

        nz0 = size(vt%var,3)
        allocate(depth0(nz0)); depth0 = abs(vt%z(1:nz0))
        to0 = ocean_annual(vt, trim(what)//" to")
        if (has_so) then
            so0 = ocean_annual(sc%snap(s)%bind(f_so)%src(1), trim(what)//" so")
            if (size(so0,3) /= nz0) then
                write(*,*) "snapshots:: Error: "//trim(what)//": so and to have different depth axes."
                stop 1
            end if
        end if

        ! Smooth the ocean (spec%ocn_sigma), level by level
        if (sc%snap(s)%spec%ocn_sigma > 0.0_wp) then
            do k = 1, nz0
                call filter_gaussian(to0(:,:,k), sc%snap(s)%spec%ocn_sigma, sc%par%dx_km)
                if (has_so) call filter_gaussian(so0(:,:,k), sc%snap(s)%spec%ocn_sigma, sc%par%dx_km)
            end do
        end if

        call model_depth(depthm)
        nzo = size(depthm)
        if (allocated(st%depth))  deallocate(st%depth)
        allocate(st%depth(nzo)); st%depth = depthm
        call realloc3(st%to_ann, nx, ny, nzo)
        call realloc3(st%so_ann, nx, ny, nzo)

        do k = 1, nzo
        do j = 1, ny
        do i = 1, nx
            st%to_ann(i,j,k) = interp_linear(depth0, to0(i,j,:), xout=depthm(k))
        end do
        end do
        end do
        if (has_so) then
            do k = 1, nzo
            do j = 1, ny
            do i = 1, nx
                st%so_ann(i,j,k) = interp_linear(depth0, so0(i,j,:), xout=depthm(k))
            end do
            end do
            end do
        else
            st%so_ann = 0.0_wp
        end if

        end associate

        return
    end subroutine reduce_ocean

    function ocean_annual(vs, what) result(var)
        ! The annual mean of an ocean field as loaded (nx,ny,nz,nt): a monthly
        ! climatology (nt = 12, one year) is averaged; an annual field (nt = 1:
        ! no time axis, or the record at the forcing time) is used as is.
        implicit none
        type(varslice_class), intent(IN) :: vs
        character(len=*),     intent(IN) :: what
        real(wp), allocatable :: var(:,:,:)
        integer :: nt

        nt = size(vs%var,4)
        if (nt /= 1 .and. nt /= NMONTH) then
            write(*,*) "snapshots:: Error: "//trim(what)//": an ocean field is a monthly &
                       &climatology (12 slices) or annual (1 slice); got", nt, "slices."
            stop 1
        end if
        var = sum(vs%var, dim=4) / real(nt,wp)

        return
    end function ocean_annual

    subroutine realloc2(var, nx, ny)
        real(wp), allocatable, intent(INOUT) :: var(:,:)
        integer,               intent(IN)    :: nx, ny
        if (allocated(var)) deallocate(var)
        allocate(var(nx,ny))
    end subroutine realloc2

    subroutine realloc3(var, nx, ny, nz)
        real(wp), allocatable, intent(INOUT) :: var(:,:,:)
        integer,               intent(IN)    :: nx, ny, nz
        if (allocated(var)) deallocate(var)
        allocate(var(nx,ny,nz))
    end subroutine realloc3

    subroutine model_depth(depth)
        ! Model ocean depth axis: nzo=23, 0..2000 m in 21 steps + [2500, 3000].
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
        type(snapshots_class), intent(IN) :: sc
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

    logical function binds(sc, s, name)
        ! Snapshot s binds the field `name` (a registry field with a var_defs group).
        implicit none
        type(snapshots_class), intent(IN) :: sc
        integer,             intent(IN) :: s
        character(len=*),    intent(IN) :: name
        integer :: f
        f = field_bind_index(sc, name)
        binds = .FALSE.
        if (f > 0) binds = sc%snap(s)%spec%nsrc(f) > 0
        return
    end function binds

    function interp_linear(x, y, xout) result(yout)
        ! Simple linear point interpolation, clamped to the endpoints. Private copy
        ! so the ocean vinterp does not depend on which interp1D module wins the
        ! include path in a full yelmox build.
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

end module snapshots
