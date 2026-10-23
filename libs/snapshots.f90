module snapshots
    ! snapshots — climate boundary forcing built on the fesm-utils primitives.
    !
    !   Loading   -> varslice  (one varslice_class per snapshot x field)
    !   Indices   -> the driver's registry (kryos_forcing), read by name
    !   Wrapper   -> this module (blend + transforms + output state)
    !
    ! See docs/climate-snapshots.md for its configuration.
    !
    ! Config model: the [forcing] var_defs databases (a list of nml files)
    ! define varslice groups (permanent per-file variable metadata). A group
    ! is taken from the par file if it has it, else from the first database
    ! that does (var_defs_find); the snapshot groups &snapshot_<name> are
    ! looked up the same way. Each snapshot group names the group(s) supplying
    ! each field (1 = monthly, 2 = [ann, season] -> synthesized). Their paths
    ! take the keywords {grid_global} ([snapshots] grid_global, the global grid
    ! of the paleoclimate products) and {exp} (the snapshot's exp), besides
    ! {domain} and {grid_name}.
    ! Fields are optional per snapshot (key absent or "" = not bound): the backend
    ! builds the fields the reference snapshot binds, the atmosphere (tas and pr,
    ! with zs) and/or the ocean (to, optionally so); an ocean-only reference serves
    ! the rembo ocean. Every snapshot a term names must bind the term's fields
    ! (checked at init).
    !
    ! Terms (climate_terms, pattern snapshots): each spans two snapshots
    ! (snap0 -> snap1) and names a coefficient of the driver's registry
    ! ([forcing] coefs, kryos_forcing), or none (the constant 1); the value
    ! of the coefficient is the position between them (0 at snap0, 1 at snap1),
    ! evaluated by climate_terms (terms_eval) and read from the terms here.
    ! The reference is [forcing<sfx>] reference; the snapshots of a run are
    ! the reference and the snap0/snap1 of its terms. A field sums the terms
    ! that list it on top of the reference (Banderas et al., 2018, for two),
    ! additively (tas, to, so) or as a ratio (pr), fixed by field. With
    ! snap0 = "", snap1 is itself the anomaly or the ratio: such a snapshot
    ! binds no zs, its fields are taken as given, at sea level (no elevation
    ! correction of their own). The coefficients restart with the registry, so
    ! snapshots holds no state of its own. The uniform and replace terms are
    ! applied by climate_terms after the backend.

    use precision, only : wp, sp, dp
    use ncio
    use nml
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end
    use climate_terms, only : terms_class
    use yelmox_nml,    only : VAR_DEFS_MAX, var_defs_find, forcing_var_defs, def_file => yelmox_defaults_file
    use gaussian_filter, only : filter_gaussian
    use coordinates, only : grid_class

    implicit none

    real(wp), parameter :: pi = 3.14159265358979323846_wp

    private

    integer, parameter :: SNP_MAX = 32   ! max snapshots / fields
    integer, parameter :: TERM_MAX = 16  ! max terms of one field
    integer, parameter :: NMONTH   = 12

    ! =====================================================================
    ! Types
    ! =====================================================================

    ! One physical field in the registry (immutable after load).
    type field_spec_class
        character(len=32) :: name           ! "tas", "pr", "to", "so", "zs"
        character(len=16) :: kind           ! from the name: "atm_temp" | "atm_precip" | "ocn" | "ocn_salt" | "elev"
        character(len=16) :: blend          ! "linear" | "ratio" | "const", fixed by the field name
        integer           :: n_term          ! number of terms on this field
        integer           :: term_k(TERM_MAX) ! their positions in snapshots_class%term
    end type

    ! One climate state: PD, LGM, piControl, ... (config side).
    type snapshot_spec_class
        character(len=64) :: name
        logical           :: monthly        ! raw atm fields are monthly (else annual+season -> synthesized)
        character(len=3)  :: tas_season     ! season of the second tas field: "jja" | "djf"
        real(wp)          :: clim_sigma     ! [km] Gaussian smoothing of the atmosphere at loading (0: none)
        real(wp)          :: ocn_sigma      ! [km] Gaussian smoothing of the ocean at loading (0: none)
        character(len=64) :: exp            ! experiment, the {exp} of its var_defs paths
        character(len=64), allocatable :: refs(:,:)  ! (n_field,2) var_defs groups per field ("" = none)
        integer,           allocatable :: nsrc(:)    ! (n_field) number of groups per field (0 = not bound)
    end type

    ! A climate state: the fields in use allocated (allocated() is the "in
    ! use" guard). Used both as a snapshot's reduced (sea-level) fields and as
    ! the blended output.
    type climate_state_class
        ! atmosphere, monthly (nx,ny,12): tas at surface, tsl at sea level
        real(wp), allocatable :: tas(:,:,:), tsl(:,:,:), pr(:,:,:), prcor(:,:,:)
        ! atmosphere, annual (nx,ny)
        real(wp), allocatable :: ta_ann(:,:), pr_ann(:,:), tsl_ann(:,:)
        real(wp), allocatable :: beta_p(:,:)          ! precip-temperature sensitivity
        ! ocean
        real(wp), allocatable :: to_ann(:,:,:), so_ann(:,:,:)
        real(wp), allocatable :: depth(:)
        ! geometry
        real(wp), allocatable :: z_srf(:,:)
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
        character(len=256) :: var_defs(VAR_DEFS_MAX) ! [forcing] var_defs: databases, searched in order after the par file
        character(len=512) :: par_file      ! the par file, searched first
        character(len=64)  :: grid_global   ! global grid of the paleoclimate products ({grid_global})
        character(len=64)  :: group         ! base namelist group (already domain-suffixed)
        integer            :: nx, ny
        real(wp)           :: dx_km         ! [km] grid spacing (snapshot smoothing)
        logical            :: south         ! southern hemisphere (seasons, lapse rates)
        character(len=64)  :: ref_name      ! which snapshot is the reference
        integer            :: iref          ! its position in snap
        logical            :: with_atm      ! the reference binds the atmosphere (tas, pr)
        logical            :: with_ocn      ! the reference binds the ocean (to)
        real(wp)           :: lapse(2)      ! [annual, summer] lapse rate [K/m]
        ! precip sensitivity parameters
        real(wp)           :: f_p           ! precip-temperature Clausius factor -> beta_p [1/K]
        integer            :: n_snap, n_field, n_term
    end type

    ! One snapshot term: its position in the domain's terms (name, coefficient
    ! and its value f there) and its two snapshots.
    type snap_term_class
        integer           :: k              ! position in terms%term
        integer           :: snap0, snap1   ! positions in snap (snap0 = 0: none,
                                            ! snap1 is the anomaly/ratio)
    end type

    ! The object callers hold.
    type snapshots_class
        type(snapshots_param_class)        :: par
        type(field_spec_class),  allocatable :: registry(:)  ! field definitions
        type(snapshot_class),    allocatable :: snap(:)      ! was clim0..clim3
        type(snap_term_class),   allocatable :: term(:)      ! the snapshot terms
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

    subroutine snapshots_init(sc, filename, domain, grid, time, south, terms, group, reference)
        ! `grid` is the climate grid: the snapshot fields are read on it, or
        ! remapped onto it online (var_defs groups with remap/grid_src). `time`
        ! primes the time-varying snapshots. `terms` are the domain's terms
        ! (climate_terms), of which those with the snapshots pattern are read
        ! here; `group` the backend group ([snapshots<sfx>]) and `reference`
        ! the reference snapshot ([forcing<sfx>] reference).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: domain
        type(grid_class),      intent(IN)    :: grid
        real(dp),              intent(IN)    :: time
        logical,               intent(IN)    :: south     ! southern hemisphere
        type(terms_class),     intent(IN)    :: terms
        character(len=*),      intent(IN)    :: group, reference

        sc%par%nx = grid%G%nx
        sc%par%ny = grid%G%ny
        sc%par%dx_km = real(grid%G%dx*grid%cs%xy_conv*1e-3_dp, wp)   ! [km] (clim_sigma)
        sc%par%south = south

        ! Read all configuration: top-level, field registry, the terms, the
        ! snapshot specs and their bindings; then check the bindings against
        ! each other.
        call snapshots_par_load(sc, filename, trim(group), domain, trim(grid%name), terms, reference)
        call snapshots_check_bindings(sc, terms)

        ! Load the snapshot fields (varslice), reduce them, and set the reference state.
        call snapshots_load_snapshots(sc, time, grid)

        return
    end subroutine snapshots_init

    subroutine snapshots_update(sc, z_srf, time, terms)
        ! Pipeline: refresh loads -> combine (the position on each term: the
        ! value f of its coefficient in `terms`, evaluated for this update by
        ! terms_eval) -> transform (the atmosphere, if bound). The uniform and
        ! replace terms are applied by the caller afterwards (climate_terms).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        real(dp),              intent(IN)    :: time
        type(terms_class),     intent(IN)    :: terms

        ! 1. Refresh time-varying loads and re-reduce those snapshots.
        call snapshots_refresh_loads(sc, time)

        ! 2. Combine the reduced snapshot states into `now`, per field.
        call snapshots_combine(sc, terms)

        ! 3. Transform: inflate sea-level fields to the current model elevation
        !    (tsl -> tas), precip, annual means.
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

    subroutine snapshots_combine(sc, terms)
        ! Blend the reduced snapshot states into `now`, per field bound in the
        ! reference, from the reference baseline. Each term k on a field
        ! (snap0_k -> snap1_k) adds with weight a_k = the value of its
        ! coefficient (terms, the f of terms_eval):
        !   linear (tas, to, so): now = ref + sum_k a_k*(snap1_k - snap0_k)
        !   ratio  (pr):          now = ref * prod_k (a_k*(snap1_k/snap0_k - 1) + 1)
        !   const  (zs):          now = ref
        ! A term with snap0 = "" adds snap1 itself: a_k*snap1_k (linear) or the
        ! factor a_k*(snap1_k - 1) + 1 (ratio). Two terms as in Banderas et al.
        ! (2018). The snapshots of the terms bind the field
        ! (snapshots_check_bindings).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        type(terms_class),     intent(IN)    :: terms
        integer  :: f, j, k, s0, s1
        logical  :: blended
        real(wp) :: aa
        type(climate_state_class) :: none   ! the snap0 of a term with snap0 = ""

        do f = 1, sc%par%n_field
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

            ! linear/ratio fields add their terms; const fields stay at the
            ! reference.
            blended = trim(sc%registry(f)%blend) == "linear" .or. &
                      trim(sc%registry(f)%blend) == "ratio"
            if (.not. blended) cycle

            do j = 1, sc%registry(f)%n_term
                k  = sc%registry(f)%term_k(j)
                s0 = sc%term(k)%snap0
                s1 = sc%term(k)%snap1
                aa = terms%term(sc%term(k)%k)%f
                if (s0 > 0) then
                    call add_term(sc%registry(f), sc%now, sc%snap(s0)%state, sc%snap(s1)%state, &
                                  aa, anomaly=.FALSE.)
                else
                    call add_term(sc%registry(f), sc%now, none, sc%snap(s1)%state, &
                                  aa, anomaly=.TRUE.)
                end if
            end do
        end do

        return
    end subroutine snapshots_combine

    subroutine add_term(fs, now, x0, x1, aa, anomaly)
        ! Add one term (x0 -> x1, weight aa; anomaly: x1 alone) to the field fs
        ! of `now`.
        implicit none
        type(field_spec_class),    intent(IN)    :: fs
        type(climate_state_class), intent(INOUT) :: now
        type(climate_state_class), intent(IN)    :: x0, x1
        real(wp),                  intent(IN)    :: aa
        logical,                   intent(IN)    :: anomaly

        select case(trim(fs%kind))
            case("atm_temp")
                call blend_term(fs%blend, now%tsl,    x0%tsl,    x1%tsl,    aa, anomaly)
            case("atm_precip")
                call blend_term(fs%blend, now%prcor,  x0%prcor,  x1%prcor,  aa, anomaly)
            case("ocn")
                call blend_term(fs%blend, now%to_ann, x0%to_ann, x1%to_ann, aa, anomaly)
            case("ocn_salt")
                call blend_term(fs%blend, now%so_ann, x0%so_ann, x1%so_ann, aa, anomaly)
        end select

        return
    end subroutine add_term

    subroutine blend_term(mode, dst, x0, x1, aa, anomaly)
        ! Add one term (x0 -> x1, weight aa) to dst, which holds the reference
        ! plus the terms added so far. anomaly (snap0 = ""): x1 is the anomaly
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
    end subroutine blend_term

    subroutine snapshots_transform(sc, z_srf)
        ! Inflate the blended sea-level fields to the current model elevation and
        ! compute the annual means. The monthly lapse rate is the one of the
        ! sea-level reduction (lapse_monthly), so a snapshot at its own
        ! elevation is recovered exactly.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        real(wp),              intent(IN)    :: z_srf(:,:)
        integer  :: m
        real(wp) :: lapse(NMONTH)

        lapse = snapshots_lapse_monthly(sc)

        sc%now%z_srf = z_srf

        ! 3a: monthly tas from tsl via seasonal lapse rate to the current elevation.
        do m = 1, NMONTH
            sc%now%tas(:,:,m) = sc%now%tsl(:,:,m) - z_srf*lapse(m)
        end do

        ! 3b: monthly precip re-sensitized to the current elevation.
        do m = 1, NMONTH
            sc%now%pr(:,:,m) = sc%now%prcor(:,:,m) * &
                exp(sc%now%beta_p*(sc%now%tas(:,:,m) - sc%now%tsl(:,:,m)))
        end do

        ! 3c: annual means.
        sc%now%tsl_ann = sum(sc%now%tsl, dim=3) / 12.0_wp
        sc%now%ta_ann  = sum(sc%now%tas, dim=3) / 12.0_wp
        sc%now%pr_ann  = sum(sc%now%pr,  dim=3) / 12.0_wp * 365.0_wp

        return
    end subroutine snapshots_transform

    function snapshots_lapse_monthly(sc) result(lapse)
        ! Monthly lapse rate [K/m] of the backend (lapse_monthly with its pars).
        implicit none
        type(snapshots_class), intent(IN) :: sc
        real(wp) :: lapse(NMONTH)
        lapse = lapse_monthly(sc%par%lapse(1), sc%par%lapse(2), sc%par%south)
    end function snapshots_lapse_monthly

    pure function lapse_monthly(l1, l2, south) result(lapse)
        ! Monthly lapse rate [K/m] between sea level and the surface: the annual
        ! l1 with a seasonal cosine of amplitude l1-l2 on the calendar of the
        ! season synthesis, c(m) = cos(2*pi*(m-1)/12) (+1 in January, -1 in
        ! July), so the summer lapse rate l2 holds exactly in July (north) or
        ! January (south). Used for the sea-level reduction of a snapshot and
        ! for the inflation of the blend to the model surface: the same rates
        ! both ways, so a snapshot at its own elevation is recovered exactly.
        implicit none
        real(wp), intent(IN) :: l1, l2
        logical,  intent(IN) :: south
        real(wp) :: lapse(NMONTH)
        integer  :: m
        do m = 1, NMONTH
            if (south) then
                lapse(m) = l1 + (l2-l1)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp)
            else
                lapse(m) = l1 + (l1-l2)*cos(2*pi*(m*30.0_wp-30.0_wp)/360.0_wp)
            end if
        end do
    end function lapse_monthly

    subroutine snapshots_end(sc)
        ! Teardown: the varslice sources of every snapshot, the snapshots, terms
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
        deallocate(sc%snap, sc%term, sc%registry)

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
        ! (axes from snapshots_write_init): the surface elevation. The climate
        ! as applied (after the climate modifiers) and the terms are written by
        ! the caller (climate_write_2D).
        implicit none
        type(snapshots_class), intent(IN) :: sc
        character(len=*),      intent(IN) :: filename
        integer,               intent(IN) :: ncid, n

        associate(now => sc%now)
        if (allocated(now%z_srf))   call write2D("z_srf",   now%z_srf,   "m",    "Surface elevation", &
                                                 prec_abs=1e-3_wp)
        end associate

        return

    contains

        subroutine write2D(vnm, var, units, long_name, prec_abs, prec_rel)
            ! prec_abs/prec_rel: precision of the written values (ncio bit
            ! rounding), absent for lossless.
            character(len=*), intent(IN) :: vnm, units, long_name
            real(wp),         intent(IN) :: var(:,:)
            real(wp),         intent(IN), optional :: prec_abs, prec_rel
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                          units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
        end subroutine write2D

    end subroutine snapshots_write_step

    ! =====================================================================
    ! Parameter loading & setup
    ! =====================================================================

    subroutine snapshots_par_load(sc, filename, group, domain, grid_name, terms, reference)
        ! Read the full configuration: the top-level knobs ([snapshots<sfx>],
        ! template snapshots), the field registry (the fields a snapshot may
        ! bind, fixed blend rule by name), the snapshot terms of `terms`
        ! (pattern snapshots: snap0, snap1, coef, fields), and the snapshot
        ! specs and their field bindings (&snapshot_<name>, template snapshot)
        ! of the reference and the terms' snapshots. The snapshot groups are
        ! read from the par file if it has them, else from the first
        ! [forcing] var_defs file that does; [snapshots<sfx>] from the par
        ! file only. Keys not set take the defaults of input/yelmox_defaults.nml.
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group      ! backend group (already domain-suffixed)
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name
        type(terms_class),     intent(IN)    :: terms
        character(len=*),      intent(IN)    :: reference

        character(len=64) :: names(SNP_MAX)
        character(len=256):: sgroup
        character(len=512):: gfile        ! the file holding sgroup
        integer           :: i, j, k, f, n

        sc%par%group     = group
        sc%par%par_file  = filename
        sc%par%domain    = domain
        sc%par%grid_name = grid_name
        sc%par%var_defs  = forcing_var_defs(filename)
        sc%par%ref_name  = reference

        ! --- top-level --------------------------------------------------------
        call nml_validate(filename, def_file, group, defaults_group="snapshots")
        call nml_read(filename, group, "grid_global", sc%par%grid_global, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshots")
        call nml_read(filename, group, "lapse",       sc%par%lapse,       init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshots")
        call nml_read(filename, group, "f_p",         sc%par%f_p,         init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshots")

        ! --- field registry ---------------------------------------------------
        call nml_read(filename, group, "fields", names, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshots")
        sc%par%n_field = count(len_trim(names) > 0)
        if (allocated(sc%registry)) deallocate(sc%registry)
        allocate(sc%registry(sc%par%n_field))
        do i = 1, sc%par%n_field
            call read_field_spec(sc%registry(i), group, trim(names(i)))
        end do

        ! --- the snapshot terms, and the snapshots they and the reference name -
        n = 0
        do k = 1, terms%n
            if (trim(terms%term(k)%pattern) == "snapshots") n = n + 1
        end do
        sc%par%n_term = n
        if (allocated(sc%term)) deallocate(sc%term)
        allocate(sc%term(n))
        names    = ""
        names(1) = reference
        sc%par%n_snap = 1
        i = 0
        do k = 1, terms%n
            associate(t => terms%term(k))
            if (trim(t%pattern) /= "snapshots") cycle
            i = i + 1
            sc%term(i)%k = k
            if (len_trim(t%snap0) > 0) call add_name(t%snap0)
            call add_name(t%snap1)
            end associate
        end do

        ! --- snapshots and their bindings ------------------------------------
        if (allocated(sc%snap)) deallocate(sc%snap)
        allocate(sc%snap(sc%par%n_snap))
        do i = 1, sc%par%n_snap
            sgroup = "snapshot_"//trim(names(i))
            gfile  = var_defs_find(sc%par%var_defs, filename, trim(sgroup), &
                                   "a snapshot of the forcing of "//trim(group))
            call read_snapshot_spec(sc%snap(i)%spec, trim(gfile), trim(sgroup), &
                                    trim(names(i)), sc%registry, sc%par%south)
        end do
        sc%par%iref = snapshot_position(sc, sc%par%ref_name, "reference")

        ! --- the terms' snapshots, and the fields each term adds to -----------
        i = 0
        do k = 1, terms%n
            associate(t => terms%term(k))
            if (trim(t%pattern) /= "snapshots") cycle
            i = i + 1
            sc%term(i)%snap0 = 0
            if (len_trim(t%snap0) > 0) sc%term(i)%snap0 = snapshot_position(sc, t%snap0, "snap0 of &"//trim(t%group))
            sc%term(i)%snap1 = snapshot_position(sc, t%snap1, "snap1 of &"//trim(t%group))
            do j = 1, size(t%fields)
                if (len_trim(t%fields(j)) == 0) cycle
                f = field_position(sc, t%fields(j), "field of &"//trim(t%group))
                if (sc%registry(f)%n_term >= TERM_MAX) then
                    write(*,*) "snapshots:: Error: more than ", TERM_MAX, " terms on field "// &
                               trim(t%fields(j))//"."
                    stop 1
                end if
                sc%registry(f)%n_term = sc%registry(f)%n_term + 1
                sc%registry(f)%term_k(sc%registry(f)%n_term) = i
            end do
            end associate
        end do

        ! --- what the backend builds: the fields of the reference -------------
        sc%par%with_atm = binds(sc, sc%par%iref, "tas")
        sc%par%with_ocn = binds(sc, sc%par%iref, "to")

        return

    contains

        subroutine add_name(name)
            ! Add a snapshot name to `names` once.
            character(len=*), intent(IN) :: name
            integer :: m
            do m = 1, sc%par%n_snap
                if (trim(names(m)) == trim(name)) return
            end do
            if (sc%par%n_snap >= SNP_MAX) then
                write(*,*) "snapshots:: Error: more than ", SNP_MAX, " snapshots."
                stop 1
            end if
            sc%par%n_snap = sc%par%n_snap + 1
            names(sc%par%n_snap) = name
        end subroutine add_name

    end subroutine snapshots_par_load

    subroutine snapshots_check_bindings(sc, terms)
        ! Check the field bindings of the snapshots against each other (stop
        ! with an error on the first mismatch):
        !  - number of var_defs groups: tas 1 (monthly) or 2 ([ann, season]),
        !    pr 1 or 2 (monthly: sum, e.g. snowfall + rainfall) or 1 (annual),
        !    every other field 1;
        !  - a state snapshot (the reference, or of a term with a snap0) binding
        !    tas or pr binds zs (sea-level reduction); a snapshot binding zs and
        !    pr binds tas (the precipitation is reduced with the temperature).
        !    Only the snap1 of a term without snap0 may leave out zs: its fields
        !    are then taken as given, at sea level;
        !  - so needs to;
        !  - the reference binds the atmosphere as tas and pr together, or neither;
        !  - no snapshot binds a field the reference does not (it would not be
        !    built);
        !  - every snapshot of a term on a blended (linear, ratio) field binds it.
        implicit none
        type(snapshots_class), intent(IN) :: sc
        type(terms_class),     intent(IN) :: terms   ! the names of the terms (messages)

        integer :: s, f, j, k, n, iref
        logical :: state
        character(len=256) :: sgroup

        iref = sc%par%iref

        do s = 1, sc%par%n_snap
            associate(spec => sc%snap(s)%spec)
            sgroup = "&snapshot_"//trim(spec%name)

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

            ! A state: the reference, or a snapshot of a term with a snap0.
            state = (s == iref)
            do k = 1, sc%par%n_term
                if (sc%term(k)%snap0 > 0 .and. (sc%term(k)%snap0 == s .or. sc%term(k)%snap1 == s)) &
                    state = .TRUE.
            end do
            if (state .and. (binds(sc, s, "tas") .or. binds(sc, s, "pr")) .and. &
                .not. binds(sc, s, "zs")) &
                call binding_error(trim(sgroup)//" binds tas or pr but not zs (its elevation, &
                                   &for the sea-level reduction); only the snap1 of a term &
                                   &without snap0 may leave it out.")
            if (binds(sc, s, "zs") .and. binds(sc, s, "pr") .and. .not. binds(sc, s, "tas")) &
                call binding_error(trim(sgroup)//" binds pr and zs but not tas (the &
                                   &precipitation is reduced to sea level with the temperature).")
            if (binds(sc, s, "so") .and. .not. binds(sc, s, "to")) &
                call binding_error(trim(sgroup)//" binds so but not to.")

            end associate
        end do

        if (binds(sc, iref, "tas") .neqv. binds(sc, iref, "pr")) &
            call binding_error("the reference &snapshot_"//trim(sc%par%ref_name)// &
                               " binds one of tas and pr: the atmosphere needs both (or neither).")

        ! The snapshots of the terms on every blended field the reference binds.
        do f = 1, sc%par%n_field
            if (sc%snap(iref)%spec%nsrc(f) == 0) cycle
            if (trim(sc%registry(f)%blend) /= "linear" .and. trim(sc%registry(f)%blend) /= "ratio") cycle
            do j = 1, sc%registry(f)%n_term
                k = sc%registry(f)%term_k(j)
                if (sc%term(k)%snap0 > 0) call check_term_snapshot(sc%term(k)%snap0, "snap0")
                call check_term_snapshot(sc%term(k)%snap1, "snap1")
            end do
        end do

        return

    contains

        subroutine check_term_snapshot(s, end_name)
            integer,          intent(IN) :: s
            character(len=*), intent(IN) :: end_name
            if (sc%snap(s)%spec%nsrc(f) == 0) &
                call binding_error("term "//trim(terms%term(sc%term(k)%k)%name)//" lists the field "// &
                                   trim(sc%registry(f)%name)//", but its "//end_name// &
                                   " snapshot "//trim(sc%snap(s)%spec%name)//" does not bind it.")
        end subroutine check_term_snapshot

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
                   ") is not a snapshot of this run."
        stop
    end function snapshot_position

    integer function field_position(sc, name, what) result(f)
        ! Position of the field `name` in the registry; stops if it is not a field.
        type(snapshots_class), intent(IN) :: sc
        character(len=*),    intent(IN) :: name, what
        do f = 1, sc%par%n_field
            if (trim(sc%registry(f)%name) == trim(name)) return
        end do
        write(*,*) "snapshots:: Error: field '"//trim(name)//"' ("//trim(what)// &
                   ") is not one of [snapshots] fields."
        stop 1
    end function field_position

    subroutine read_field_spec(fs, group, name)
        ! One field of the registry: its kind and blend rule follow from its
        ! name (tas, to, so: linear; pr: ratio; zs: const); the terms on it
        ! are added by snapshots_par_load.
        implicit none
        type(field_spec_class), intent(OUT) :: fs
        character(len=*),       intent(IN)  :: group, name

        fs%name   = name
        fs%n_term = 0
        fs%term_k = 0

        select case(trim(name))
            case("tas")
                fs%kind  = "atm_temp"
                fs%blend = "linear"
            case("pr")
                fs%kind  = "atm_precip"
                fs%blend = "ratio"
            case("to")
                fs%kind  = "ocn"
                fs%blend = "linear"
            case("so")
                fs%kind  = "ocn_salt"
                fs%blend = "linear"
            case("zs")
                fs%kind  = "elev"
                fs%blend = "const"
            case default
                write(*,*) "snapshots:: Error: &"//trim(group)//" fields: unknown field '"// &
                           trim(name)//"' (tas, pr, to, so, zs)."
                stop 1
        end select

        return
    end subroutine read_field_spec

    subroutine read_snapshot_spec(ss, filename, group, name, registry, south)
        ! The spec of one snapshot and its bindings: per registry field, the
        ! var_defs group(s) under the field's name (key absent or "" = not bound).
        ! tas_season: the season of the second tas field of an annual snapshot,
        ! "jja" or "djf" ("", the default: the summer of the hemisphere). Unset
        ! keys from the template snap_snapshot of the defaults file.
        implicit none
        type(snapshot_spec_class), intent(OUT) :: ss
        character(len=*),          intent(IN)  :: filename, group, name
        type(field_spec_class),    intent(IN)  :: registry(:)
        logical,                   intent(IN)  :: south

        character(len=64) :: refs(2)
        integer :: f

        ss%name      = name
        call nml_validate(filename, def_file, group, defaults_group="snapshot")
        call nml_read(filename, group, "exp",        ss%exp,        init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshot")
        call nml_read(filename, group, "monthly",    ss%monthly,    init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshot")
        call nml_read(filename, group, "clim_sigma", ss%clim_sigma, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshot")
        call nml_read(filename, group, "ocn_sigma",  ss%ocn_sigma,  init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshot")
        call nml_read(filename, group, "tas_season", ss%tas_season, init=.TRUE., &
                      defaults_file=def_file, defaults_group="snapshot")
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

        ! Bindings: 1 group (monthly) or 2 ([ann, season]; pr monthly: summed).
        allocate(ss%refs(size(registry),2), ss%nsrc(size(registry)))
        ss%refs = ""
        do f = 1, size(registry)
            call nml_read(filename, group, trim(registry(f)%name), refs, init=.TRUE., &
                          defaults_file=def_file, defaults_group="snapshot")
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
        ! (each from the par file if it has it, else the first database of
        ! var_defs that does, with the path
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
                nsrc = sc%snap(s)%spec%nsrc(f)
                if (nsrc == 0) cycle          ! this snapshot does not bind this field

                allocate(sc%snap(s)%bind(f)%src(nsrc))
                do k = 1, nsrc
                    vdb = var_defs_find(sc%par%var_defs, trim(sc%par%par_file), &
                                        trim(sc%snap(s)%spec%refs(f,k)), &
                                        "&snapshot_"//trim(sc%snap(s)%spec%name)// &
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
        ! rates, beta_p); one without (the snap1 of a term with snap0 = "") is
        ! taken as given, a sea-level anomaly (tas) or ratio (pr).
        implicit none
        type(snapshots_class), intent(INOUT) :: sc
        integer,               intent(IN)    :: s

        integer  :: nx, ny, m, k, nsrc
        integer  :: f_tas, f_pr, f_zs
        logical  :: south, has_tas, has_pr, has_zs
        real(wp) :: zmean, l1, l2
        real(wp) :: lapse(NMONTH)
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

        ! --- elevation, beta_p ----------------------------------------------------
        if (has_zs) then
            call realloc2(st%z_srf, nx, ny)
            call realloc2(st%beta_p, nx, ny)
            st%z_srf = bind(f_zs)%src(1)%var(:,:,1,1)
            ! Missing-value fill -> domain mean (GRL inputs are complete; guard on the
            ! -9999 fill sentinel, safe since surface elevation is >= 0).
            if (any(st%z_srf <= -9000.0_wp)) then
                zmean = sum(st%z_srf, mask=st%z_srf > -9000.0_wp) &
                        / real(max(count(st%z_srf > -9000.0_wp),1), wp)
                where (st%z_srf <= -9000.0_wp) st%z_srf = zmean
            end if
            ! beta_p = f_p everywhere
            st%beta_p = sc%par%f_p
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

        ! Smooth the snapshot (spec%clim_sigma, e.g. coarse GCM fields).
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
        end if

        ! --- sea-level temperatures and elevation-desensitized precip -----------
        if (.not. has_zs) then
            ! Taken as given: a sea-level anomaly or ratio.
            if (has_tas) then
                st%tsl     = st%tas
                st%tsl_ann = st%ta_ann
            end if
            if (has_pr) then
                st%prcor = st%pr
            end if
        else
            if (has_tas) then
                call realloc3(st%tsl, nx, ny, NMONTH)
                call realloc2(st%tsl_ann, nx, ny)
                st%tsl_ann = st%ta_ann + l1*st%z_srf
                lapse = lapse_monthly(l1, l2, south)
                do m = 1, NMONTH
                    st%tsl(:,:,m) = st%tas(:,:,m) + st%z_srf*lapse(m)
                end do
            end if
            if (has_pr) then   ! with tas and zs (snapshots_check_bindings)
                call realloc3(st%prcor, nx, ny, NMONTH)
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
