module helpers_mod

    implicit none
    
contains

    subroutine check_allocation(stat)

        implicit none

        integer,intent(in) :: stat

        ! Check if stat is nonzero
        if (stat /= 0) then
            write(*,*) "Your computer has insufficient memory for the desired computation. Quitting..."
            stop
        end if
    
    end subroutine check_allocation

    
end module helpers_mod