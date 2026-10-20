module yelmox_nml
    ! Namelist helpers shared by the yelmox forcing modules (indices, snapshots,
    ! the climate modifiers): the defaults file of their groups, and the
    ! lookup of a named group (varslice, snapshot, axis, index) in the par
    ! file, then in a list of database files (var_defs).

    use nml, only : nml_has_group

    implicit none

    private

    ! Defaults (schema) of the forcing groups: every key is read with
    ! nml_read(..., defaults_file = yelmox_defaults_file, defaults_group =
    ! its template group), see the file's header.
    character(len=*), parameter, public :: yelmox_defaults_file = "input/yelmox_defaults.nml"

    integer, parameter, public :: VAR_DEFS_MAX = 8   ! max files of a var_defs list

    public :: var_defs_find

contains

    function var_defs_find(var_defs, par_file, group, what) result(filename)
        ! The par file itself if it holds the group (searched implicitly first,
        ! so a par file never names itself and its groups override the
        ! databases), else the first file of the list var_defs ("" entries
        ! ignored) that holds it; stops naming the group and the files when
        ! none has it. `what` names the reader in the message (e.g.
        ! "&snap_lgm tas").
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
