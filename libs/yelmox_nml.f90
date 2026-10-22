module yelmox_nml
    ! Namelist helpers shared by the yelmox forcing modules (coefficients,
    ! terms, snapshots, the climate transforms): the defaults file of their groups,
    ! the run's database list ([forcing] var_defs) and the lookup of a named
    ! group (varslice, snapshot, term, coef) in the par file, then in that
    ! list.

    use nml, only : nml_has_group, nml_read

    implicit none

    private

    ! Defaults (schema) of the forcing groups: every key is read with
    ! nml_read(..., defaults_file = yelmox_defaults_file, defaults_group =
    ! its template group), see the file's header.
    character(len=*), parameter, public :: yelmox_defaults_file = "input/yelmox_defaults.nml"

    integer, parameter, public :: VAR_DEFS_MAX = 8   ! max files of a var_defs list

    public :: forcing_var_defs, var_defs_find

contains

    function forcing_var_defs(par_file) result(var_defs)
        ! The run's database list, [forcing] var_defs (the one list every
        ! named group is looked up in after the par file: varslice groups,
        ! snapshots, terms, coefficients); the default of the template when the
        ! group or the key is absent.
        character(len=*), intent(in) :: par_file
        character(len=256) :: var_defs(VAR_DEFS_MAX)
        call nml_read(par_file, "forcing", "var_defs", var_defs, init=.true., &
                      defaults_file=yelmox_defaults_file, defaults_group="forcing")
    end function forcing_var_defs

    function var_defs_find(var_defs, par_file, group, what) result(filename)
        ! The par file itself if it holds the group (searched implicitly first,
        ! so a par file never names itself and its groups override the
        ! databases), else the first file of the list var_defs ("" entries
        ! ignored) that holds it; stops naming the group and the files when
        ! none has it. `what` names the reader in the message (e.g.
        ! "&snapshot_lgm tas").
        character(len=*), intent(in) :: var_defs(:), par_file, group, what
        character(len=max(len(var_defs),len(par_file))) :: filename

        integer :: i

        if (nml_has_group(trim(par_file), trim(group))) then
            filename = par_file
            return
        end if
        do i = 1, size(var_defs)
            if (len_trim(var_defs(i)) == 0) cycle
            if (nml_has_group(trim(var_defs(i)), trim(group))) then
                filename = var_defs(i)
                return
            end if
        end do

        write(*,*) "var_defs_find:: error: group &"//trim(group)//" ("//trim(what)// &
                   ") is in neither the par file nor any of the var_defs files:"
        write(*,*) "    "//trim(par_file)
        do i = 1, size(var_defs)
            if (len_trim(var_defs(i)) > 0) write(*,*) "    "//trim(var_defs(i))
        end do
        error stop 1

    end function var_defs_find

end module yelmox_nml
