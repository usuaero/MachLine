module myjson_m
#ifdef dnad
    use dnadmod
#ifndef ndv
#define ndv 1
#endif
#endif
    use json_m
    implicit none

    logical :: json_found
#ifdef dnad
    integer, save :: n_design_vars = 1  ! First design variable is always alpha
#endif

    interface myjson_get
        module procedure :: myjson_value_get_real, myjson_file_get_real
#ifdef dnad
        module procedure :: myjson_value_get_dual, myjson_file_get_dual
#endif
        module procedure :: myjson_value_get_integer, myjson_file_get_integer
        module procedure :: myjson_value_get_string, myjson_file_get_string
    end interface myjson_get

#ifdef dnad
    interface json_value_add
        module procedure :: myjson_value_add_dual, myjson_value_add_dual_vec
    end interface
#endif

contains

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_get_real(json, name, value, default_value)
    implicit none
    type(json_value),intent(in),pointer :: json
    character(len=*), intent(in) :: name
    real, intent(out) :: value
    real, intent(in), optional :: default_value

    call json_get(json, name, value, json_found)
    if(json_failed() .or. (.not. json_found)) then
        if (present(default_value)) then
            write(*,*) trim(name),' set to ',default_value
            value = default_value
            call json_clear_exceptions()
        else
            write(*,*) 'Error: Unable to read required value: ',name
            STOP
        end if
    end if
end subroutine myjson_value_get_real

#ifdef dnad
!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_get_dual(json, name, value, default_value)
    implicit none
    type(json_value),intent(in),pointer :: json
    character(len=*), intent(in) :: name
    type(dual), intent(out) :: value
    real, intent(in), optional :: default_value

    real :: rvalue
    real, dimension(:), allocatable :: vec

    call json_get(json, name, rvalue, json_found)
    if((.not.json_failed()) .and. json_found) then
        value = rvalue
    else
        call json_clear_exceptions()
        call json_get(json, name, vec, json_found)
        if(json_found .and. (.not. json_failed())) then
            value = vec(1)  ! This will initialize derivatives to zero
            if(vec(2) /= 0) then
                if(n_design_vars < ndv) then
                    n_design_vars = n_design_vars + 1
                    value%dx(n_design_vars) = vec(2)
                    write(*,*) 'Design Variable ', n_design_vars, ': ', name
                else
                    write(*,*) 'Error: The number of design variables exceeds the compiled limit: ', ndv
                    write(*,*) '       Reduce the number of design variables, or increase the limit by'
                    write(*,*) '       specifying -Dndv=<num> when compiling.'
                    STOP
                end if
            end if
        else
            if (present(default_value)) then
                value = default_value
                write(*,*) trim(name),' set to ', value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ',name
                STOP
            end if
        end if
    end if

end subroutine myjson_value_get_dual

#endif
!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_get_integer(json, name, value, default_value)
    implicit none
    type(json_value),intent(in),pointer :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    integer, intent(in), optional :: default_value

    call json_get(json, name, value, json_found)
    if((.not.json_found) .or. json_failed()) then
        if (present(default_value)) then
            write(*,*) trim(name),' set to ',default_value
            value = default_value
            call json_clear_exceptions()
        else
            write(*,*) 'Error: Unable to read required value: ', name
            STOP
        end if
    end if

end subroutine myjson_value_get_integer

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_get_string(json, name, value, default_value)
    implicit none
    type(json_value), intent(in), pointer :: json
    character(len=*), intent(in) :: name
    character(:), allocatable, intent(out) :: value
    character(len=*), intent(in), optional :: default_value

    call json_get(json, name, value, json_found)
    if((.not.json_found) .or. json_failed()) then
        if (present(default_value)) then
            write(*,*) trim(name), ' set to ', default_value
            value = default_value
            call json_clear_exceptions()
        else
            write(*,*) 'Error: Unable to read required value: ', name
            STOP
        end if
    end if

end subroutine myjson_value_get_string

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_file_get_real(json, name, value, default_value)
    implicit none
    type(json_file) :: json
    character(len=*), intent(in) :: name
    real, intent(out) :: value
    real, intent(in), optional :: default_value

    call json%get(name, value)
    if(json_failed()) then
        if (present(default_value)) then
            write(*,*) trim(name), ' set to ', default_value
            value = default_value
            call json_clear_exceptions()
        else
            write(*,*) 'Error: Unable to read required value: ', trim(name)
            STOP
        end if
    end if

end subroutine myjson_file_get_real

#ifdef dnad
!-----------------------------------------------------------------------------------------------------------
subroutine myjson_file_get_dual(json, name, value, default_value)
    implicit none
    type(json_file) :: json
    character(len=*), intent(in) :: name
    type(dual), intent(out) :: value
    real, intent(in), optional :: default_value

    real :: temp
    real, dimension(:), allocatable :: vec

    call json%get(name, temp)
    if(.not. json_failed()) then
        value = temp  ! Derivatives not specified, initialize to zero
    else
        call json_clear_exceptions()
        call json%get(name, vec)
        if(.not. json_failed()) then
            value = vec(1)  ! This will initialize derivatives to zero
            if(vec(2) /= 0) then
                if(n_design_vars < ndv) then
                    n_design_vars = n_design_vars + 1
                    value%dx(n_design_vars) = vec(2)
                    write(*,*) 'Design Variable ', n_design_vars, ': ', name
                else
                    write(*,*) 'Error: The number of design variables exceeds the compiled limit: ', ndv
                    write(*,*) '       Reduce the number of design variables, or increase the limit by'
                    write(*,*) '       specifying -Dndv=<num> when compiling.'
                    STOP
                end if
            end if
        else
            if (present(default_value)) then
                value = default_value
                write(*,*) trim(name),' set to ', value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ', trim(name)
                STOP
            end if
        end if
    end if

end subroutine myjson_file_get_dual

#endif
!-----------------------------------------------------------------------------------------------------------
subroutine myjson_file_get_integer(json, name, value, default_value)
    implicit none
    type(json_file) :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    integer, intent(in), optional :: default_value

    call json%get(name, value)
    if(json_failed()) then
        if (present(default_value)) then
            write(*,*) trim(name),' set to ',default_value
            value = default_value
            call json_clear_exceptions()
        else
            write(*,*) 'Error: Unable to read required value: ',name
            STOP
        end if
    end if
end subroutine myjson_file_get_integer

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_file_get_string(json, name, value)
    implicit none
    type(json_file) :: json
    character(len=*), intent(in) :: name
    character(:), allocatable, intent(out) :: value

    call json%get(name, value)
    if(json_failed()) then
        write(*,*) 'Error: Unable to read required value: ',name
        STOP
    end if

    value = trim(value)
end subroutine myjson_file_get_string

!-----------------------------------------------------------------------------------------------------------
subroutine json_check()
    if(json_failed()) then
        call print_json_error_message()
        STOP
    end if
end subroutine json_check

!-----------------------------------------------------------------------------------------------------------
subroutine print_json_error_message()
    implicit none
    character(len=:),allocatable :: error_msg
    logical :: status_ok

    !get error message:
    call json_check_for_errors(status_ok, error_msg)

    !print it if there is one:
    if (.not. status_ok) then
        write(*,'(A)') error_msg
        deallocate(error_msg)
        call json_clear_exceptions()
    end if

end subroutine print_json_error_message

#ifdef dnad
!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_add_dual(me, name, val)

    implicit none

    type(json_value), pointer   :: me
    character(len=*),intent(in) :: name
    type(dual),intent(in)       :: val

    real,allocatable,dimension(:) :: vec
    integer :: vec_length

    ! Dual numbers are written to json as a vector containing: [u, du/dx, du/dy, ...]
    vec_length = size(val%dx) + 1
    allocate(vec(vec_length))
    vec(1) = val%x
    vec(2:) = val%dx

    call json_value_add(me, name, vec)

end subroutine myjson_value_add_dual

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_add_dual_vec(me, name, val)

    implicit none

    type(json_value), pointer         :: me
    character(len=*),intent(in)       :: name
    type(dual),dimension(:),intent(in)  :: val

    type(json_value),pointer :: var
    integer :: i

    !create the variable as an array:
    call json_value_create(var)
    call to_array(var,name)

    !populate the array:
    do i=1,size(val)
        call json_value_add(var, '', val(i))
    end do

    !add it:
    call json_value_add(me, var)

    !cleanup:
    nullify(var)

end subroutine myjson_value_add_dual_vec

#endif
end module myjson_m
