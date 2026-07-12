module snapclim2
    ! snapclim2 — climate boundary forcing built on the fesm-utils primitives.
    !
    !   Loading   -> varslice  (one varslice_class per snapshot x field)
    !   Indices   -> tsgen      (one tsgen_class per driving index)
    !   Wrapper   -> this module (blend + transforms + derive + output state + restart)
    !
    ! Replacement for the monolithic `snapclim`. See docs/snapclim2-design.md for the
    ! design rationale, the unified five-knob model, and the migration plan.
    !
    ! STATUS: skeleton. Types and the public API are in place and compile; the
    ! par-load and the load->combine->transform->derive pipeline are stubbed (see TODOs).

    use precision, only : wp, sp, dp
    use ncio
    use nml
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use tsgen,    only : tsgen_class, tsgen_init, tsgen_update, &
                         tsgen_restart_write, tsgen_restart_read

    implicit none

    private

    integer, parameter :: SNP2_MAX = 32   ! max entries in a snapshots/fields/indices list
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

    ! One loaded snapshot: its spec + one varslice per enabled field (registry order).
    type snapshot_class
        type(snapshot_spec_class)         :: spec
        type(varslice_class), allocatable :: fld(:)
    end type

    ! Output state. Option A: named members for the hot path, unallocated when a
    ! run does not use them (allocated() is the "in use" guard); `extra` for the tail.
    type climate_state_class
        ! atmosphere, monthly (nx,ny,12)
        real(wp), allocatable :: tas(:,:,:), pr(:,:,:), sf(:,:,:), pr_stdev_frac(:,:,:)
        ! atmosphere, annual / summer (nx,ny)
        real(wp), allocatable :: ta_ann(:,:), ta_sum(:,:), pr_ann(:,:)
        real(wp), allocatable :: tsl_ann(:,:), tsl_sum(:,:), prcor_ann(:,:)
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

    ! Top-level parameters.
    type snapclim2_param_class
        character(len=256) :: domain
        character(len=256) :: grid_name
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
    type snapclim2_class
        type(snapclim2_param_class)          :: par
        type(field_spec_class),  allocatable :: registry(:)  ! field definitions
        type(snapshot_class),    allocatable :: snap(:)      ! was clim0..clim3
        type(tsgen_class),       allocatable :: idx(:)       ! was at..bs (each may be nc-channel)
        character(len=64),       allocatable :: idx_name(:)  ! index names (for labels/accessors)
        type(climate_state_class)            :: now          ! current output
        type(climate_state_class)            :: ref          ! distinguished reference (was clim0)
    end type

    public :: snapclim2_class
    public :: snapclim2_init
    public :: snapclim2_update
    public :: snapclim2_end
    public :: snapclim2_write_init
    public :: snapclim2_write_step
    public :: snapclim2_restart_write
    public :: snapclim2_restart_read

contains

    ! =====================================================================
    ! Lifecycle
    ! =====================================================================

    subroutine snapclim2_init(sc, filename, domain, grid_name, nx, ny, time, basins, group)
        ! Mirrors snapclim_init, plus `time` (required to initialize the tsgen indices,
        ! which anchor ramp-type series to their start time).
        implicit none
        type(snapclim2_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name
        integer,               intent(IN)    :: nx, ny
        real(wp),              intent(IN)    :: time
        real(wp),              intent(IN)    :: basins(:,:)
        character(len=*),      intent(IN), optional :: group

        character(len=64) :: base_group
        integer           :: k

        base_group = "snp2"
        if (present(group)) base_group = trim(group)

        ! Read all configuration: top-level, field registry, snapshot specs, index names.
        call snapclim2_par_load(sc, filename, trim(base_group), domain, grid_name)
        sc%par%nx = nx
        sc%par%ny = ny

        ! Initialize the driving indices (tsgen). Group is &tsgen_<base_group>_<idxname>.
        do k = 1, sc%par%n_idx
            call tsgen_init(sc%idx(k), filename, time, &
                            label=trim(base_group)//"_"//trim(sc%idx_name(k)))
        end do

        ! Load the snapshot fields (varslice) and set the reference state.
        call snapclim2_load_snapshots(sc, filename, basins)

        return
    end subroutine snapclim2_init

    subroutine snapclim2_update(sc, z_srf, time, domain, dTa, dTo, dSo, dx, basins)
        ! Mirrors snapclim_update. Pipeline: advance indices -> refresh loads ->
        ! combine -> transform -> derive (see docs/snapclim2-design.md).
        implicit none
        type(snapclim2_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        real(wp),              intent(IN)    :: time
        character(len=*),      intent(IN)    :: domain
        real(wp),              intent(IN), optional :: dTa, dTo, dSo, dx
        real(wp),              intent(IN)    :: basins(:,:)

        ! 1. TODO(indices):  tsgen_update(sc%idx(k), time, var) -> weights source
        ! 2. TODO(loads):    varslice_update(sc%snap(s)%fld(f), time, method)
        ! 3. TODO(combine):  now%f = ref%f + Sum w_s*(snap(s)%f - ref%f)   [or absolute]
        ! 4. TODO(transform):lapse (tsl<->tas via source/target z_srf), precip scaling,
        !                    seasonal synthesis, aggregates, ocean vertical interp
        ! 5. TODO(derive):   optional rules, e.g. ocn anom = f_to * mean(atm anom)

        return
    end subroutine snapclim2_update

    subroutine snapclim2_end(sc)
        ! Teardown (snapclim never had one). Free varslice/tsgen and state arrays.
        implicit none
        type(snapclim2_class), intent(INOUT) :: sc

        ! TODO: varslice_end on each snap%fld, deallocate snap/idx/registry and states.

        return
    end subroutine snapclim2_end

    ! =====================================================================
    ! Diagnostic output
    ! =====================================================================

    subroutine snapclim2_write_init(sc, filename, time_init)
        implicit none
        type(snapclim2_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        real(wp),              intent(IN) :: time_init

        ! TODO: create NetCDF file with x/y/month/depth/time dims; write static fields.

        return
    end subroutine snapclim2_write_init

    subroutine snapclim2_write_step(sc, filename, time)
        implicit none
        type(snapclim2_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        real(wp),              intent(IN) :: time

        ! TODO: append `now` fields, generalized over the enabled registry.

        return
    end subroutine snapclim2_write_step

    ! =====================================================================
    ! Restart (+ provenance record)
    ! =====================================================================

    subroutine snapclim2_restart_write(sc, filename, time)
        implicit none
        type(snapclim2_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        real(wp),              intent(IN) :: time

        ! TODO: tsgen_restart_write per index (real prognostic state) + a provenance
        !       record documenting snapshots/registry/weights (documentation only).

        return
    end subroutine snapclim2_restart_write

    subroutine snapclim2_restart_read(sc, filename)
        implicit none
        type(snapclim2_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename

        ! TODO: tsgen_restart_read per index; fields are recomputed from time on the
        !       next update, so the provenance record is not reloaded to drive the run.

        return
    end subroutine snapclim2_restart_read

    ! =====================================================================
    ! Parameter loading & setup
    ! =====================================================================

    subroutine snapclim2_par_load(sc, filename, group, domain, grid_name)
        ! Read the full configuration: top-level knobs, the field registry
        ! (&<group>_field_<name>), the snapshot specs (&<group>_snap_<name>), and
        ! the index names (each an &tsgen_<group>_<name> group, read by tsgen_init).
        implicit none
        type(snapclim2_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group      ! base group (already domain-suffixed)
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name

        character(len=64) :: names(SNP2_MAX)
        integer           :: i

        sc%par%group     = group
        sc%par%domain    = domain
        sc%par%grid_name = grid_name

        ! --- top-level --------------------------------------------------------
        sc%par%combine   = "anomaly"
        sc%par%manifold  = 1
        sc%par%ref_name  = ""
        sc%par%lapse     = 0.0_wp
        sc%par%dTa_const = 0.0_wp
        sc%par%dTo_const = 0.0_wp
        sc%par%dSo_const = 0.0_wp
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
                                    trim(group)//"_snap_"//trim(names(i)), &
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
    end subroutine snapclim2_par_load

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

    subroutine snapclim2_load_snapshots(sc, filename, basins)
        ! Load each snapshot's enabled fields via varslice and set the reference state.
        !
        ! PENDING (varslice {snapshot} binding, option b): each field's varslice config
        ! lives on &<group>_field_<field> with a {snapshot} token in the path, resolved
        ! per snapshot (with an optional &<group>_snap_<snap>_<field> override). This needs
        ! a {snapshot}-style substitution hook in varslice (see docs/snapclim2-design.md);
        ! stubbed until that mechanism is settled.
        implicit none
        type(snapclim2_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        real(wp),              intent(IN)    :: basins(:,:)

        return
    end subroutine snapclim2_load_snapshots

end module snapclim2
