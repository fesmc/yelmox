module yelmox_nml
    ! Namelist helpers shared by the yelmox forcing modules (snapesm, the
    ! climate modifiers): the lookup of a varslice group in a list of
    ! database files (var_defs).

    use nml, only : nml_has_group

    implicit none

    private

    integer, parameter, public :: VAR_DEFS_MAX = 8   ! max files of a var_defs list

    public :: var_defs_find

contains

    function var_defs_find(var_defs, group, what) result(filename)
        ! The first file of the list var_defs ("" entries ignored) that holds
        ! the group; stops naming the group and the files when none has it.
        ! `what` names the reader in the message (e.g. "&snap_lgm tas").
        character(len=*), intent(in) :: var_defs(:), group, what
        character(len=len(var_defs)) :: filename

        integer :: i

        do i = 1, size(var_defs)
            if (len_trim(var_defs(i)) == 0) cycle
            if (nml_has_group(trim(var_defs(i)), trim(group))) then
                filename = var_defs(i)
                return
            end if
        end do

        write(*,*) "var_defs_find:: error: group &"//trim(group)//" ("//trim(what)// &
                   ") is in none of the var_defs files:"
        do i = 1, size(var_defs)
            if (len_trim(var_defs(i)) > 0) write(*,*) "    "//trim(var_defs(i))
        end do
        error stop 1

    end function var_defs_find

end module yelmox_nml
