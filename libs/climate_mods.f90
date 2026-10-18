module climate_mods
    ! Climate modifiers: changes applied to a backend's `climate_out` after the
    ! backend update, in a fixed order (docs/_forcing-design.md, section 7):
    !   anom   -- homogeneous anomaly from one index of the driver's registry:
    !             dTa = f*f_ta, dTo = f*f_to, dSo = f*f_so
    ! (later ocn_fraction, ocn_obm, taper, bias). Each is an optional namelist
    ! group with `active` (default off), suffixed by the domain (anom<sfx>).
    !
    ! Contract: each modifier states what it touches in `now` (`ref` is never
    ! touched) and which backends it refuses; climate_mods_init checks every
    ! active modifier against the backend and stops on a mismatch, never a
    ! silent partial application.
    !   anom: tas, ta_ann, tsl_ann + dTa; pr, pr_ann * exp(beta_p*dTa); the ocean
    !         profiles to_ann/so_ann + dTo/dSo, or at the shelf base T_shlf/S_shlf
    !         and dT_shlf/dS_shlf. Refused when the backend supplies the smb
    !         (has_smb), which a post-hoc air anomaly cannot reach. Exception:
    !         rembo takes the anomaly as an input (dT_sum) and reads [anom]
    !         itself (yelmox_climate, rembo_update); anom is skipped here for it.

    use precision,     only : wp
    use nml,           only : nml_read, nml_has_param
    use climate_out,   only : climate_out_class
    use kryos_forcing, only : indices_class, indices_value, indices_has

    implicit none

    private

    ! [anom<sfx>]: one index mapped onto homogeneous anomalies by the gains.
    type anom_class
        logical           :: active = .false.
        character(len=64) :: index  = ""        ! the index driving the anomalies
        real(wp) :: f_ta = 1.0_wp               ! [K]   dTa = f * f_ta
        real(wp) :: f_to = 0.0_wp               ! [K]   dTo = f * f_to
        real(wp) :: f_so = 0.0_wp               ! [psu] dSo = f * f_so
    end type anom_class

    type climate_mods_class
        type(anom_class) :: anom
    end type climate_mods_class

    public :: anom_class
    public :: climate_mods_class
    public :: climate_mods_init
    public :: climate_mods_apply

contains

    subroutine climate_mods_init(mods, filename, sfx, idx, backend, has_smb)
        ! Read the modifier groups (suffix sfx) and check each active modifier's
        ! contract against the backend (its name and whether it supplies the
        ! smb). Call after the backend init.
        type(climate_mods_class), intent(inout) :: mods
        character(len=*),         intent(in)    :: filename, sfx, backend
        type(indices_class),      intent(in)    :: idx
        logical,                  intent(in)    :: has_smb

        call anom_init(mods%anom, filename, "anom"//trim(sfx), idx, backend, has_smb)

    end subroutine climate_mods_init

    subroutine climate_mods_apply(mods, out, idx, backend)
        ! The active modifiers on `out`, in the fixed order. idx is the driver's
        ! registry, advanced to the forcing time.
        type(climate_mods_class), intent(in)    :: mods
        type(climate_out_class),  intent(inout) :: out
        type(indices_class),      intent(in)    :: idx
        character(len=*),         intent(in)    :: backend

        real(wp) :: f

        ! anom (rembo applies it itself)
        if (mods%anom%active .and. trim(backend) /= "rembo") then
            f = indices_value(idx, mods%anom%index)
            call climate_mod_anom(out, f*mods%anom%f_ta, f*mods%anom%f_to, f*mods%anom%f_so)
        end if

    end subroutine climate_mods_apply

    ! ===== anom ============================================================

    subroutine anom_init(anom, filename, group, idx, backend, has_smb)
        ! Read [anom<sfx>] (optional group; inactive without it) and check its
        ! contract: a registry index, and a backend that does not supply the smb
        ! (but rembo, which reads the group itself).
        type(anom_class),    intent(inout) :: anom
        character(len=*),    intent(in)    :: filename, group, backend
        type(indices_class), intent(in)    :: idx
        logical,             intent(in)    :: has_smb

        anom%active = .false.
        if (nml_has_param(filename, group, "active")) &
            call nml_read(filename, group, "active", anom%active)
        if (.not. anom%active) return

        call nml_read(filename, group, "index", anom%index)
        call nml_read(filename, group, "f_ta",  anom%f_ta)
        call nml_read(filename, group, "f_to",  anom%f_to)
        call nml_read(filename, group, "f_so",  anom%f_so)

        if (.not. indices_has(idx, anom%index)) then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] index = '"// &
                       trim(anom%index)//"' is not in [indices] names."
            error stop 1
        end if
        if (has_smb .and. trim(backend) /= "rembo") then
            write(*,*) "climate_mods_init:: error: ["//trim(group)//"] needs a climate that &
                       &does not supply the surface mass balance (an air anomaly cannot reach &
                       &it); climate = "//trim(backend)//" supplies it (esm: use the esm &
                       &smb anomalies or surface_method /= climate)."
            error stop 1
        end if

        write(*,*) "climate_mods_init:: "//trim(group)//" active, index "//trim(anom%index)// &
                   "; f_ta/f_to/f_so =", anom%f_ta, anom%f_to, anom%f_so

    end subroutine anom_init

    subroutine climate_mod_anom(out, dTa, dTo, dSo)
        ! Homogeneous anomaly on `now`: the air temperature + dTa, the
        ! precipitation * exp(beta_p*dTa) (beta_p: the backend's sensitivity),
        ! the ocean + dTo/dSo, as profiles or, without them, at the shelf base.
        type(climate_out_class), intent(inout) :: out
        real(wp),                intent(in)    :: dTa, dTo, dSo

        integer :: m

        associate(now => out%now)

        now%tas    = now%tas + dTa
        now%ta_ann = now%ta_ann + dTa
        if (allocated(now%tsl_ann)) now%tsl_ann = now%tsl_ann + dTa

        do m = 1, size(now%pr, 3)
            now%pr(:,:,m) = now%pr(:,:,m) * exp(now%beta_p*dTa)
        end do
        if (allocated(now%pr_ann)) now%pr_ann = now%pr_ann * exp(now%beta_p*dTa)

        if (out%has_ocn_shelf) then
            now%T_shlf  = now%T_shlf  + dTo
            now%S_shlf  = now%S_shlf  + dSo
            now%dT_shlf = now%dT_shlf + dTo
            now%dS_shlf = now%dS_shlf + dSo
        else
            if (allocated(now%to_ann)) now%to_ann = now%to_ann + dTo
            if (allocated(now%so_ann)) now%so_ann = now%so_ann + dSo
        end if

        end associate

    end subroutine climate_mod_anom

end module climate_mods
