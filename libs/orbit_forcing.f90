module orbit_forcing
    ! The orbit of the insolation ([orbit<sfx>], one per domain): the orbital
    ! elements that apply at a time of the forcing clock. Shared by the surface
    ! models that take insolation (chion's host for itm/bessi, smbpal's itm);
    ! each owns its insol_class (fesmc/insol), initialised once with its
    ! latitudes and the tables in fldr, and recomputes its insolation only when
    ! the elements change. This module only resolves the elements:
    !
    !   method = "clock": the elements at the forcing clock time_rel (Laskar et al. 2004 tables)
    !            "time":  the elements at a fixed time [yr, model time: relative to 1950, past < 0]
    !            "fixed": the given elements ecc, obl [deg], per [deg] (longitude of perihelion)
    !
    ! Read with input/yelmox_defaults.nml: a par file lists only what differs.

    use nml,        only : nml_read, nml_validate
    use insolation, only : dp, insol_class, orbit_pars_class, orbit_pars_at, orbit_pars_fixed
    use yelmox_nml, only : def_file => yelmox_defaults_file

    implicit none

    private

    type orbit_ctl_class
        character(len=16)  :: method    ! clock | time | fixed
        real(dp)           :: time      ! [yr] model time of the orbit (time)
        real(dp)           :: ecc       ! [1]   eccentricity (fixed)
        real(dp)           :: obl       ! [deg] obliquity (fixed)
        real(dp)           :: per       ! [deg] longitude of perihelion (fixed)
        character(len=512) :: fldr      ! folder of the orbital tables
    end type orbit_ctl_class

    public :: orbit_ctl_class
    public :: orbit_pars_class
    public :: orbit_ctl_load
    public :: orbit_pars_now
    public :: orbit_pars_equal

contains

    subroutine orbit_ctl_load(ctl, filename, group)
        ! Read [orbit<sfx>] (optional group; the defaults without it).
        type(orbit_ctl_class), intent(out) :: ctl
        character(len=*),      intent(in)  :: filename, group

        type(orbit_pars_class) :: orbit

        call nml_validate(filename, def_file, group, defaults_group="orbit")
        call nml_read(filename, group, "method", ctl%method, defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "time",   ctl%time,   defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "ecc",    ctl%ecc,    defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "obl",    ctl%obl,    defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "per",    ctl%per,    defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "fldr",   ctl%fldr,   defaults_file=def_file, defaults_group="orbit")

        select case(trim(ctl%method))
            case("clock", "time")
                ! The tables are read by each host's insol_init.
            case("fixed")
                orbit = orbit_pars_fixed(ctl%ecc, ctl%obl, ctl%per)   ! stops on invalid elements
            case default
                write(*,*) "orbit_ctl_load:: error: ["//trim(group)//"] method must be &
                           &clock, time or fixed; got '"//trim(ctl%method)//"'."
                error stop 1
        end select
    end subroutine orbit_ctl_load

    function orbit_pars_now(ctl, is, time_rel) result(orbit)
        ! The orbital elements at the forcing clock time_rel [yr], from the
        ! tables of the host's insol_class (clock, time).
        type(orbit_ctl_class), intent(in) :: ctl
        type(insol_class),     intent(in) :: is
        real(dp),              intent(in) :: time_rel
        type(orbit_pars_class) :: orbit

        select case(trim(ctl%method))
            case("clock")
                orbit = orbit_pars_at(is, time_rel)
            case("time")
                orbit = orbit_pars_at(is, ctl%time)
            case default   ! fixed
                orbit = orbit_pars_fixed(ctl%ecc, ctl%obl, ctl%per)
        end select
    end function orbit_pars_now

    elemental function orbit_pars_equal(a, b) result(same)
        ! The same elements (the insolation need not be recomputed).
        type(orbit_pars_class), intent(in) :: a, b
        logical :: same

        same = (a%ecc == b%ecc .and. a%obl == b%obl .and. a%per == b%per)
    end function orbit_pars_equal

end module orbit_forcing
