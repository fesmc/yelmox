module snapesm
    ! snapesm — climate boundary forcing built on the fesm-utils primitives.
    !
    !   Loading   -> varslice  (one varslice_class per snapshot x field)
    !   Indices   -> tsgen      (one tsgen_class per driving index)
    !   Wrapper   -> this module (blend + transforms + derive + output state + restart)
    !
    ! Replacement for the monolithic `snapclim`. See docs/snapesm-design.md for the
    ! design rationale, the unified five-knob model, and the migration plan.
    !
    ! Config model: a `var_defs` database nml defines varslice groups (permanent
    ! per-file variable metadata); each climate-state group &<group>_<snapshot> names
    ! the group(s) supplying each field (1 = monthly, 2 = [ann, sum] -> synthesized).
    !
    ! STATUS: config/loading + the pipeline structure compile; the numeric physics
    ! (reduce / weights / transform / derive) is stubbed against snapclim references.

    use precision, only : wp, sp, dp
    use ncio
    use nml
    use interp1D, only : interp_linear
    use series,   only : series_interp1
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use tsgen,    only : tsgen_class, tsgen_init, tsgen_update, &
                         tsgen_restart_write, tsgen_restart_read

    implicit none

    real(wp), parameter :: pi = 3.14159265358979323846_wp

    private

    integer, parameter :: SNP_MAX = 32   ! max entries in a snapshots/fields/indices list
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
        character(len=16) :: blend          ! "linear" | "ratio" | "fraction" | "const"
        character(len=16) :: index          ! name of the driving index (idx) for this field, or ""
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
        real(wp)          :: time           ! reference time [yr] of this snapshot (index normalization)
        real(wp)          :: idx_coord(2)   ! position on the 1-D / 2-D index manifold
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
        character(len=16)  :: combine       ! default combine mode
        integer            :: manifold      ! index-manifold dimension (0/1/2)
        character(len=64)  :: ref_name      ! which snapshot is the reference
        real(wp)           :: lapse(2)      ! [annual, summer] lapse rate [K/m]
        ! precip / ocean sensitivity parameters (snapclim &snap)
        real(wp)           :: f_p           ! precip-temperature Clausius factor -> beta_p [1/K]
        real(wp)           :: f_p_ne        ! NE-Greenland beta_p multiplier (snapclim hardcodes 1.0)
        real(wp)           :: f_stdev       ! precip-variability multiplier
        real(wp)           :: f_to          ! ocean/atm anomaly ratio (fraction ocean rule)
        real(wp)           :: f_hol         ! Holocene index scaling (snap_1ind_new)
        ! bipolar/obm compatibility (obm_coupling reads par%dTa_const today)
        real(wp)           :: dTa_const, dTo_const, dSo_const
        integer            :: n_snap, n_field, n_idx
    end type

    ! The object callers hold.
    type snapesm_class
        type(snapesm_param_class)          :: par
        type(field_spec_class),  allocatable :: registry(:)  ! field definitions
        type(snapshot_class),    allocatable :: snap(:)      ! was clim0..clim3
        type(tsgen_class),       allocatable :: idx(:)       ! was at..bs (each may be nc-channel)
        character(len=64),       allocatable :: idx_name(:)  ! index names (for labels/accessors)
        type(climate_state_class)            :: now          ! current output
        type(climate_state_class)            :: ref          ! distinguished reference (was clim0)
    end type

    public :: snapesm_class
    public :: snapesm_init
    public :: snapesm_update
    public :: snapesm_end
    public :: snapesm_write_init
    public :: snapesm_write_step
    public :: snapesm_restart_write
    public :: snapesm_restart_read

contains

    ! =====================================================================
    ! Lifecycle
    ! =====================================================================

    subroutine snapesm_init(sc, filename, domain, grid_name, nx, ny, time, basins, group)
        ! Mirrors snapclim_init, plus `time` (required to initialize the tsgen indices,
        ! which anchor ramp-type series to their start time).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name
        integer,               intent(IN)    :: nx, ny
        real(wp),              intent(IN)    :: time
        real(wp),              intent(IN)    :: basins(:,:)
        character(len=*),      intent(IN), optional :: group

        character(len=64) :: base_group
        integer           :: k

        base_group = "snp"
        if (present(group)) base_group = trim(group)

        ! Read all configuration: top-level, field registry, snapshot specs, index names.
        call snapesm_par_load(sc, filename, trim(base_group), domain, grid_name)
        sc%par%nx = nx
        sc%par%ny = ny

        ! Initialize the driving indices (tsgen). Group is &<base_group>_idx_<idxname>.
        do k = 1, sc%par%n_idx
            call tsgen_init(sc%idx(k), filename, time, &
                            group=trim(base_group)//"_idx_"//trim(sc%idx_name(k)))
        end do

        ! Load the snapshot fields (varslice), reduce them, and set the reference state.
        call snapesm_load_snapshots(sc, filename, time, basins)

        return
    end subroutine snapesm_init

    subroutine snapesm_update(sc, z_srf, time, domain, dTa, dTo, dSo, dx, basins)
        ! Mirrors snapclim_update. Pipeline: advance indices -> refresh loads ->
        ! combine -> transform -> derive (see docs/snapesm-design.md).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        real(wp),              intent(IN)    :: time
        character(len=*),      intent(IN)    :: domain
        real(wp),              intent(IN), optional :: dTa, dTo, dSo, dx
        real(wp),              intent(IN)    :: basins(:,:)

        real(wp), allocatable :: ta_ann_lag(:,:)

        ! 1. Advance the driving indices (tsgen) to `time`.
        call snapesm_advance_indices(sc, time)

        ! 2. Refresh time-varying loads and re-reduce those snapshots.
        call snapesm_refresh_loads(sc, time)

        ! 3. Combine the reduced snapshot states into `now`, per field.
        call snapesm_combine(sc, time)

        ! Capture ta_ann as it stands BEFORE the transform overwrites it: the
        ! fraction ocean rule (derive) reads the previous step's value (see there).
        if (allocated(sc%now%ta_ann)) ta_ann_lag = sc%now%ta_ann

        ! 4. Transform: inflate sea-level fields to the current model elevation
        !    (tsl -> tas), precip, annual/summer aggregates.
        call snapesm_transform(sc, z_srf)

        ! 5. Derive: coupled rules (e.g. fraction ocean = f_to * atm anomaly).
        call snapesm_derive(sc, ta_ann_lag, dTa, dTo, dSo)

        return
    end subroutine snapesm_update

    ! ---------------------------------------------------------------------
    ! Update pipeline stages
    ! ---------------------------------------------------------------------

    subroutine snapesm_advance_indices(sc, time)
        ! Advance each tsgen index. Our indices are series/analytic (no feedback),
        ! so the model-response argument is a dummy.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: time
        integer  :: k
        real(wp) :: response
        response = 0.0_wp
        do k = 1, sc%par%n_idx
            call tsgen_update(sc%idx(k), time, var=response)
        end do
        return
    end subroutine snapesm_advance_indices

    subroutine snapesm_refresh_loads(sc, time)
        ! Re-slice only the time-varying fields (recon, transient orography, etc.);
        ! static fields were primed once at load. A re-reduce of touched snapshots
        ! follows once the reduction is implemented.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: time
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

    subroutine snapesm_combine(sc, time)
        ! Blend the reduced snapshot states into `now`, per field, from the reference
        ! baseline (was clim0). Each field uses its blend rule and driving index:
        !   linear (temp/salinity): now = ref + aa*(s1 - s0)        [calc_temp_1ind]
        !   ratio  (precip):        now = ref*(aa*(s1/s0 - 1) + 1)  [calc_precip_1ind]
        !   fraction / const:       now = ref  (fraction ocean is finished in derive)
        ! s0,s1 are the manifold-endpoint snapshots (idx_coord min/max among non-ref);
        ! aa is the field's index normalized to those snapshots' times.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: time
        integer  :: f, iref, is0, is1
        real(wp) :: aa

        call find_endpoints(sc, iref, is0, is1)

        do f = 1, sc%par%n_field
            if (.not. sc%registry(f)%enabled) cycle
            aa = 0.0_wp
            if (len_trim(sc%registry(f)%index) > 0) &
                aa = norm_index(sc, trim(sc%registry(f)%index), time, is0, is1)

            select case(trim(sc%registry(f)%kind))
                case("atm_temp")
                    call blend_field(sc%registry(f)%blend, sc%now%tsl, sc%ref%tsl, &
                         sc%snap(is0)%state%tsl, sc%snap(is1)%state%tsl, aa)
                case("atm_precip")
                    call blend_field(sc%registry(f)%blend, sc%now%prcor, sc%ref%prcor, &
                         sc%snap(is0)%state%prcor, sc%snap(is1)%state%prcor, aa)
                case("ocn")
                    call blend_field(sc%registry(f)%blend, sc%now%to_ann, sc%ref%to_ann, &
                         sc%snap(is0)%state%to_ann, sc%snap(is1)%state%to_ann, aa)
                case("ocn_salt")
                    call blend_field(sc%registry(f)%blend, sc%now%so_ann, sc%ref%so_ann, &
                         sc%snap(is0)%state%so_ann, sc%snap(is1)%state%so_ann, aa)
                case default
                    ! "elev" and other non-blended fields: nothing to combine.
            end select
        end do

        return
    end subroutine snapesm_combine

    subroutine blend_field(mode, dst, ref, x0, x1, aa)
        ! Combine a single reduced field. `linear`/`ratio` need the endpoint states
        ! (x0,x1) allocated; `fraction`/`const` only copy the reference. No-op unless
        ! the required arrays are allocated (skips fields a snapshot does not supply).
        implicit none
        character(len=*),      intent(IN)    :: mode
        real(wp), allocatable, intent(INOUT) :: dst(:,:,:)
        real(wp), allocatable, intent(IN)    :: ref(:,:,:), x0(:,:,:), x1(:,:,:)
        real(wp),              intent(IN)    :: aa

        if (.not. (allocated(dst) .and. allocated(ref))) return

        select case(trim(mode))
            case("fraction","const")
                dst = ref
            case("ratio")
                if (.not. (allocated(x0) .and. allocated(x1))) then
                    dst = ref; return
                end if
                where (x0 /= 0.0_wp)
                    dst = ref*(aa*(x1/x0 - 1.0_wp) + 1.0_wp)
                elsewhere
                    dst = ref
                end where
            case default   ! "linear"
                if (.not. (allocated(x0) .and. allocated(x1))) then
                    dst = ref; return
                end if
                dst = ref + aa*(x1 - x0)
        end select

        return
    end subroutine blend_field

    subroutine find_endpoints(sc, iref, is0, is1)
        ! Identify the reference snapshot and the 1-D manifold endpoints (the non-ref
        ! snapshots with the smallest and largest idx_coord). Falls back to the
        ! reference when there are no non-ref snapshots (e.g. the `const` method).
        implicit none
        type(snapesm_class), intent(IN)  :: sc
        integer,              intent(OUT) :: iref, is0, is1
        integer  :: s
        real(wp) :: cmin, cmax
        iref = 0; is0 = 0; is1 = 0
        cmin =  huge(1.0_wp); cmax = -huge(1.0_wp)
        do s = 1, sc%par%n_snap
            if (sc%snap(s)%spec%is_ref) then
                if (iref == 0) iref = s
                cycle
            end if
            if (sc%snap(s)%spec%idx_coord(1) < cmin) then
                cmin = sc%snap(s)%spec%idx_coord(1); is0 = s
            end if
            if (sc%snap(s)%spec%idx_coord(1) > cmax) then
                cmax = sc%snap(s)%spec%idx_coord(1); is1 = s
            end if
        end do
        if (iref == 0) iref = 1
        if (is0  == 0) is0  = iref
        if (is1  == 0) is1  = iref
        return
    end subroutine find_endpoints

    function norm_index(sc, idx_name, time, is0, is1) result(aa)
        ! Normalize a driving index to the manifold endpoints' snapshot times:
        !   aa = (S(time) - S(t0)) / (S(t1) - S(t0))   (snapclim 439-445)
        ! so aa=0 at snapshot is0's time and aa=1 at is1's time. Optional Holocene
        ! scaling (snap_1ind_new): aa *= f_hol when aa>1 within (-12 ka, -1 ka).
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: idx_name
        real(wp),              intent(IN) :: time
        integer,               intent(IN) :: is0, is1
        real(wp) :: aa
        integer  :: k, i
        real(wp) :: t0, t1, snow, s0v, s1v

        k = 0
        do i = 1, sc%par%n_idx
            if (trim(sc%idx_name(i)) == trim(idx_name)) then
                k = i; exit
            end if
        end do

        aa = 0.0_wp
        if (k == 0) return

        t0 = sc%snap(is0)%spec%time
        t1 = sc%snap(is1)%spec%time
        snow = series_interp1(sc%idx(k)%ser, time)
        s0v  = series_interp1(sc%idx(k)%ser, t0)
        s1v  = series_interp1(sc%idx(k)%ser, t1)
        if (s1v /= s0v) aa = (snow - s0v) / (s1v - s0v)

        if (aa > 1.0_wp .and. time > -12.0e3_wp .and. time < -1.0e3_wp) aa = aa*sc%par%f_hol

        return
    end function norm_index

    subroutine snapesm_transform(sc, z_srf)
        ! Inflate the blended sea-level fields to the current model elevation and
        ! compute the annual/summer aggregates (snapclim 779-811). Note the cosine
        ! phase here is m*30-15 (vs m*30-30 in the sea-level reduction).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        integer  :: m
        logical  :: south
        real(wp) :: l1, l2

        south = (trim(sc%par%domain) == "Antarctica")
        l1 = sc%par%lapse(1)
        l2 = sc%par%lapse(2)

        if (allocated(sc%now%z_srf)) sc%now%z_srf = z_srf
        if (allocated(sc%now%mask)) then
            sc%now%mask = 0.0_wp
            where (z_srf > 1.0_wp) sc%now%mask = 1.0_wp
        end if

        ! 3a: monthly tas from tsl via seasonal lapse rate to the current elevation.
        do m = 1, NMONTH
            if (south) then
                sc%now%tas(:,:,m) = sc%now%tsl(:,:,m) - &
                    z_srf*(l1 + (l2-l1)*cos(2*pi*(m*30.0_wp-15.0_wp)/360.0_wp))
            else
                sc%now%tas(:,:,m) = sc%now%tsl(:,:,m) - &
                    z_srf*(l1 + (l1-l2)*cos(2*pi*(m*30.0_wp-15.0_wp)/360.0_wp))
            end if
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

    subroutine snapesm_derive(sc, ta_ann_lag, dTa, dTo, dSo)
        ! Coupled/derived rules applied after the blend. The `fraction` ocean rule
        ! (snapclim 663-675) sets the ocean anomaly to f_to times the domain-mean
        ! atmospheric annual anomaly. snapclim evaluates this BEFORE now%ta_ann is
        ! recomputed, so it reads the previous step's ta_ann -> we pass ta_ann_lag
        ! (now%ta_ann as of the start of this update) to reproduce that lag exactly.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: ta_ann_lag(:,:)
        real(wp),              intent(IN), optional :: dTa, dTo, dSo
        integer  :: f
        real(wp) :: dTo_now

        do f = 1, sc%par%n_field
            if (.not. sc%registry(f)%enabled) cycle
            if (trim(sc%registry(f)%blend) == "fraction" .and. &
                trim(sc%registry(f)%kind)  == "ocn") then
                dTo_now = sc%par%f_to * sum(ta_ann_lag - sc%ref%ta_ann) &
                          / real(sc%par%nx*sc%par%ny, wp)
                if (allocated(sc%now%to_ann) .and. allocated(sc%ref%to_ann)) &
                    sc%now%to_ann = sc%ref%to_ann + dTo_now
            end if
        end do

        return
    end subroutine snapesm_derive

    subroutine snapesm_end(sc)
        ! Teardown (snapclim never had one). Free varslice/tsgen and state arrays.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc

        ! TODO: varslice_end on each snap%bind%src, deallocate snap/idx/registry and states.

        return
    end subroutine snapesm_end

    ! =====================================================================
    ! Diagnostic output
    ! =====================================================================

    subroutine snapesm_write_init(sc, filename, time_init)
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        real(wp),              intent(IN) :: time_init

        ! TODO: create NetCDF file with x/y/month/depth/time dims; write static fields.

        return
    end subroutine snapesm_write_init

    subroutine snapesm_write_step(sc, filename, time)
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        real(wp),              intent(IN) :: time

        ! TODO: append `now` fields, generalized over the enabled registry.

        return
    end subroutine snapesm_write_step

    ! =====================================================================
    ! Restart (+ provenance record)
    ! =====================================================================

    subroutine snapesm_restart_write(sc, filename, time)
        implicit none
        type(snapesm_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        real(wp),              intent(IN) :: time

        ! TODO: tsgen_restart_write per index (real prognostic state) + a provenance
        !       record documenting snapshots/registry/weights (documentation only).

        return
    end subroutine snapesm_restart_write

    subroutine snapesm_restart_read(sc, filename)
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename

        ! TODO: tsgen_restart_read per index; fields are recomputed from time on the
        !       next update, so the provenance record is not reloaded to drive the run.

        return
    end subroutine snapesm_restart_read

    ! =====================================================================
    ! Parameter loading & setup
    ! =====================================================================

    subroutine snapesm_par_load(sc, filename, group, domain, grid_name)
        ! Read the full configuration: top-level knobs, the field registry
        ! (&<group>_field_<name>), the snapshot specs (&<group>_snap_<name>), and
        ! the index names (each an &tsgen_<group>_<name> group, read by tsgen_init).
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group      ! base group (already domain-suffixed)
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name

        character(len=64) :: names(SNP_MAX)
        integer           :: i

        sc%par%group     = group
        sc%par%domain    = domain
        sc%par%grid_name = grid_name

        ! --- top-level --------------------------------------------------------
        sc%par%var_defs  = ""
        sc%par%combine   = "anomaly"
        sc%par%manifold  = 1
        sc%par%ref_name  = ""
        sc%par%lapse     = 0.0_wp
        sc%par%f_p       = 0.0_wp
        sc%par%f_p_ne    = 1.0_wp
        sc%par%f_stdev   = 0.0_wp
        sc%par%f_to      = 0.0_wp
        sc%par%f_hol     = 1.0_wp
        sc%par%dTa_const = 0.0_wp
        sc%par%dTo_const = 0.0_wp
        sc%par%dSo_const = 0.0_wp
        call nml_read(filename, group, "var_defs",  sc%par%var_defs,  init=.TRUE.)
        call nml_read(filename, group, "combine",   sc%par%combine,   init=.TRUE.)
        call nml_read(filename, group, "manifold",  sc%par%manifold,  init=.TRUE.)
        call nml_read(filename, group, "ref_name",  sc%par%ref_name,  init=.TRUE.)
        call nml_read(filename, group, "lapse",     sc%par%lapse,     init=.TRUE.)
        call nml_read(filename, group, "f_p",       sc%par%f_p,       init=.TRUE.)
        call nml_read(filename, group, "f_p_ne",    sc%par%f_p_ne,    init=.TRUE.)
        call nml_read(filename, group, "f_stdev",   sc%par%f_stdev,   init=.TRUE.)
        call nml_read(filename, group, "f_to",      sc%par%f_to,      init=.TRUE.)
        call nml_read(filename, group, "f_hol",     sc%par%f_hol,     init=.TRUE.)
        call nml_read(filename, group, "dTa_const", sc%par%dTa_const, init=.TRUE.)
        call nml_read(filename, group, "dTo_const", sc%par%dTo_const, init=.TRUE.)
        call nml_read(filename, group, "dSo_const", sc%par%dSo_const, init=.TRUE.)

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

        ! --- indices (names only; tsgen_init reads the &tsgen_* groups) --------
        names = ""
        call nml_read(filename, group, "indices", names, init=.TRUE.)
        sc%par%n_idx = count(len_trim(names) > 0)
        if (allocated(sc%idx))      deallocate(sc%idx)
        if (allocated(sc%idx_name)) deallocate(sc%idx_name)
        allocate(sc%idx(sc%par%n_idx))
        allocate(sc%idx_name(sc%par%n_idx))
        do i = 1, sc%par%n_idx
            sc%idx_name(i) = trim(names(i))
        end do

        return
    end subroutine snapesm_par_load

    subroutine read_field_spec(fs, filename, group, name, combine_default)
        implicit none
        type(field_spec_class), intent(OUT) :: fs
        character(len=*),       intent(IN)  :: filename, group, name, combine_default

        fs%name           = name
        fs%kind           = "atm_temp"
        fs%combine        = combine_default
        fs%blend          = "linear"
        fs%index          = ""
        fs%enabled        = .TRUE.
        fs%apply_lapse    = .FALSE.
        fs%precip_scaling = .FALSE.
        fs%seasonal_synth = .FALSE.
        ! Only the keys the current model uses are read (kind/blend/index). The
        ! remaining knobs (combine/enabled/*_flags) keep their struct defaults; they
        ! cannot be made optional in the namelist until nml gains an ignore-missing
        ! mode (it errors on any absent key), so they are not read here yet.
        call nml_read(filename, group, "kind",  fs%kind,  init=.TRUE.)
        call nml_read(filename, group, "blend", fs%blend, init=.TRUE.)
        call nml_read(filename, group, "index", fs%index, init=.TRUE.)

        return
    end subroutine read_field_spec

    subroutine read_snapshot_spec(ss, filename, group, name, ref_name)
        implicit none
        type(snapshot_spec_class), intent(OUT) :: ss
        character(len=*),          intent(IN)  :: filename, group, name, ref_name

        ss%name      = name
        ss%is_ref    = .FALSE.
        ss%monthly   = .FALSE.
        ss%time      = 0.0_wp
        ss%idx_coord = 0.0_wp
        call nml_read(filename, group, "monthly",   ss%monthly,   init=.TRUE.)
        call nml_read(filename, group, "time",      ss%time,      init=.TRUE.)
        call nml_read(filename, group, "idx_coord", ss%idx_coord, init=.TRUE.)
        ! The reference snapshot is designated solely by the top-level ref_name
        ! (an is_ref namelist key would have to be present in every snapshot group,
        ! given nml's error-on-missing behavior).
        if (trim(name) == trim(ref_name)) ss%is_ref = .TRUE.

        return
    end subroutine read_snapshot_spec

    subroutine snapesm_load_snapshots(sc, filename, time, basins)
        ! For each snapshot and each enabled field, read the varslice group reference(s)
        ! from the state group &<group>_<snapshot> (key = field name; 1 group = already
        ! monthly, 2 groups = [ann, sum] synthesized later), load them from the var_defs
        ! database, prime to `time`, reduce to the sea-level state, and set the reference.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        real(wp),              intent(IN)    :: time
        real(wp),              intent(IN)    :: basins(:,:)

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
                                           domain=trim(sc%par%domain), grid_name=trim(sc%par%grid_name))
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
        ! `now = clim0`). This also primes now%ta_ann for the first fraction-ocean lag.
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
        south = (trim(sc%par%domain) == "Antarctica")
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
        real(wp) :: tlo, thi
        real(wp), allocatable :: depth0(:), depthm(:), to0(:,:,:), so0(:,:,:)

        nx = sc%par%nx
        ny = sc%par%ny

        associate(st => sc%snap(s)%state, vt => sc%snap(s)%bind(f_to)%src(1), &
                                          vs => sc%snap(s)%bind(f_so)%src(1))

        ! Read all 12 months then average to the annual mean (snapclim sum(dim=4)/12).
        ! rep=12 selects the sub-annual axis via whole-year matching (get_indices uses
        ! floor(month)==year), which is robust to the sp/dp rounding of the fractional
        ! month positions that a direct fractional range would trip on.
        tlo = real(floor(minval(vt%time)), wp)   ! the (single) whole year of the axis
        call varslice_update(vt, [tlo,tlo], method="range", rep=12)
        call varslice_update(vs, [tlo,tlo], method="range", rep=12)

        nz0 = size(vt%var,3)
        nt  = size(vt%var,4)
        allocate(depth0(nz0)); depth0 = abs(vt%z(1:nz0))
        allocate(to0(nx,ny,nz0), so0(nx,ny,nz0))
        to0 = sum(vt%var, dim=4) / real(nt,wp)
        so0 = sum(vs%var, dim=4) / real(size(vs%var,4),wp)

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

end module snapesm
