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
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use tsgen,    only : tsgen_class, tsgen_init, tsgen_update, &
                         tsgen_restart_write, tsgen_restart_read

    implicit none

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
        character(len=16) :: kind           ! "atm_monthly" | "ocn_3d" | "scalar_2d"
        character(len=16) :: combine        ! "anomaly" | "absolute"
        logical           :: enabled        ! off by default for exotic fields
        logical           :: apply_lapse
        logical           :: precip_scaling
        logical           :: seasonal_synth
    end type

    ! One climate state: PD, LGM, piControl, ... (config side).
    type snapshot_spec_class
        character(len=64) :: name
        logical           :: is_ref         ! the distinguished reference state
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
        real(wp)           :: lapse(2)
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

        real(wp), allocatable :: w(:)

        ! 1. Advance the driving indices (tsgen) to `time`.
        call snapesm_advance_indices(sc, time)

        ! 2. Refresh time-varying loads and re-reduce those snapshots.
        call snapesm_refresh_loads(sc, time)

        ! 3. Combine snapshot states into `now`: now = ref + Sum w_s * snap(s)%state
        !    (anomaly form; w_ref = 0). Weights encode the active method.
        allocate(w(sc%par%n_snap))
        call snapesm_weights(sc, w)
        call snapesm_combine(sc, w)

        ! 4. Transform: inflate sea-level fields to the current model elevation
        !    (tsl -> tas), precip, annual/summer aggregates, ocean vertical interp.
        call snapesm_transform(sc, z_srf)

        ! 5. Derive: optional coupled rules (e.g. ocean anomaly = f_to * atm anomaly).
        call snapesm_derive(sc, dTa, dTo, dSo)

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

    subroutine snapesm_weights(sc, w)
        ! Compute the per-snapshot anomaly weights from the indices and each
        ! snapshot's idx_coord. now = ref + Sum w_s * snap(s)%state, with w_ref = 0.
        !
        ! Each method is a weight pattern, e.g. snap_1ind (snapclim.f90 calc_temp_1ind:
        ! temp_now = temp0 + aa*(temp2 - temp1)) => w = [0, -aa, +aa] on [ref, s1, s2].
        !
        ! TODO(weights): port the per-method weight patterns (1ind/2ind/miocene/abs)
        ! faithfully and reconcile with idx_coord/manifold. For now: `const` (w = 0,
        ! now = ref), which is exact for the const method.
        implicit none
        type(snapesm_class), intent(IN)  :: sc
        real(wp),              intent(OUT) :: w(:)
        w = 0.0_wp
        return
    end subroutine snapesm_weights

    subroutine snapesm_combine(sc, w)
        ! Generic weighted blend of the reduced snapshot states into `now`.
        ! Blends the sea-level / annual fields; surface tas and aggregates are
        ! derived afterwards in snapesm_transform.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: w(:)
        integer :: s
        logical :: first
        do s = 1, sc%par%n_snap
            first = (s == 1)
            call accum3(sc%now%tsl,    sc%ref%tsl,    sc%snap(s)%state%tsl,    w(s), first)
            call accum3(sc%now%prcor,  sc%ref%prcor,  sc%snap(s)%state%prcor,  w(s), first)
            call accum3(sc%now%to_ann, sc%ref%to_ann, sc%snap(s)%state%to_ann, w(s), first)
            call accum3(sc%now%so_ann, sc%ref%so_ann, sc%snap(s)%state%so_ann, w(s), first)
            call accum2(sc%now%beta_p, sc%ref%beta_p, sc%snap(s)%state%beta_p, w(s), first)
        end do
        return
    end subroutine snapesm_combine

    subroutine accum3(dst, ref, src, wgt, first)
        ! dst = ref (on first) then dst += wgt*src  =>  dst = ref + Sum wgt*src.
        ! No-op unless all three arrays are allocated (skips fields not in use).
        implicit none
        real(wp), allocatable, intent(INOUT) :: dst(:,:,:)
        real(wp), allocatable, intent(IN)    :: ref(:,:,:), src(:,:,:)
        real(wp),              intent(IN)    :: wgt
        logical,               intent(IN)    :: first
        if (.not. (allocated(dst) .and. allocated(ref) .and. allocated(src))) return
        if (first) dst = ref
        dst = dst + wgt*src
        return
    end subroutine accum3

    subroutine accum2(dst, ref, src, wgt, first)
        implicit none
        real(wp), allocatable, intent(INOUT) :: dst(:,:)
        real(wp), allocatable, intent(IN)    :: ref(:,:), src(:,:)
        real(wp),              intent(IN)    :: wgt
        logical,               intent(IN)    :: first
        if (.not. (allocated(dst) .and. allocated(ref) .and. allocated(src))) return
        if (first) dst = ref
        dst = dst + wgt*src
        return
    end subroutine accum2

    subroutine snapesm_transform(sc, z_srf)
        ! Inflate sea-level fields to the current model elevation and derive
        ! aggregates. Faithful port target (snapclim.f90:779-812):
        !   south = (domain=="Antarctica")
        !   tas(m)  = tsl(m) - z_srf*(lapse(1) +/- (lapse(2)-lapse(1))*cos(2*pi*(m*30-15)/360))
        !   pr(m)   = prcor(m) * exp(beta_p*(tas(m)-tsl(m)))
        !   tsl_ann = mean(tsl); ta_ann = mean(tas); *_sum = mean over DJF/JJA;
        !   pr_ann  = mean(pr)*365; prcor_ann = mean(prcor)*365
        !
        ! TODO(transform): implement once reduction populates now%tsl/prcor/beta_p.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        if (allocated(sc%now%z_srf)) sc%now%z_srf = z_srf
        return
    end subroutine snapesm_transform

    subroutine snapesm_derive(sc, dTa, dTo, dSo)
        ! Optional coupled/derived rules applied after the blend, e.g. the old
        ! `fraction` ocean rule (to anomaly = f_to * mean atm anomaly), and folding
        ! the driver anomalies dTa/dTo/dSo into the fields.
        ! TODO(derive): port `fraction` and dTa/dTo/dSo application.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        real(wp),              intent(IN), optional :: dTa, dTo, dSo
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
        sc%par%dTa_const = 0.0_wp
        sc%par%dTo_const = 0.0_wp
        sc%par%dSo_const = 0.0_wp
        call nml_read(filename, group, "var_defs",  sc%par%var_defs)
        call nml_read(filename, group, "combine",   sc%par%combine)
        call nml_read(filename, group, "manifold",  sc%par%manifold)
        call nml_read(filename, group, "ref_name",  sc%par%ref_name)
        call nml_read(filename, group, "lapse",     sc%par%lapse)
        call nml_read(filename, group, "dTa_const", sc%par%dTa_const)
        call nml_read(filename, group, "dTo_const", sc%par%dTo_const)
        call nml_read(filename, group, "dSo_const", sc%par%dSo_const)

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
        fs%kind           = "atm_monthly"
        fs%combine        = combine_default
        fs%enabled        = .TRUE.
        fs%apply_lapse    = .FALSE.
        fs%precip_scaling = .FALSE.
        fs%seasonal_synth = .FALSE.
        call nml_read(filename, group, "kind",           fs%kind)
        call nml_read(filename, group, "combine",        fs%combine)
        call nml_read(filename, group, "enabled",        fs%enabled)
        call nml_read(filename, group, "apply_lapse",    fs%apply_lapse)
        call nml_read(filename, group, "precip_scaling", fs%precip_scaling)
        call nml_read(filename, group, "seasonal_synth", fs%seasonal_synth)

        return
    end subroutine read_field_spec

    subroutine read_snapshot_spec(ss, filename, group, name, ref_name)
        implicit none
        type(snapshot_spec_class), intent(OUT) :: ss
        character(len=*),          intent(IN)  :: filename, group, name, ref_name

        ss%name      = name
        ss%is_ref    = .FALSE.
        ss%idx_coord = 0.0_wp
        call nml_read(filename, group, "is_ref",    ss%is_ref)
        call nml_read(filename, group, "idx_coord", ss%idx_coord)
        ! The top-level ref_name also designates the reference snapshot.
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

        return
    end subroutine snapesm_load_snapshots

    subroutine snapesm_reduce_snapshot(sc, s)
        ! Extract raw fields from snapshot s's varslices into its reduced state and
        ! compute the sea-level fields. Faithful port target
        ! (snapclim.f90 read_climate_snapshot, ~1918-1948):
        !   tsl_ann = ta_ann + lapse(1)*z_srf ; tsl_sum = ta_sum + lapse(2)*z_srf
        !   pr      = pr*(1 + f_stdev*pr_stdev_frac)              [variability scaling]
        !   prcor_ann = pr_ann / exp(beta_p*(ta_ann - tsl_ann))
        !   tsl(m)  = tas(m) + z_srf*(lapse(1) +/- (lapse(2)-lapse(1))*cos(2*pi*(m*30-30)/360))
        !   prcor(m)= pr(m) / exp(beta_p*(tas(m) - tsl(m)))
        ! plus: seasonal synthesis for annual-only snapshots, ocean vertical interp to the
        ! model depth axis, and the varslice %var -> state-array extraction (pending
        ! confirmation of the monthly/3D %var layout). Allocates the state arrays it fills.
        implicit none
        type(snapesm_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s
        return
    end subroutine snapesm_reduce_snapshot

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
