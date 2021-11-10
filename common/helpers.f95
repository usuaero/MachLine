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


    function mirror_about_plane(vec, plane) result(mirrored_vec)
        ! Returns a version of the given vector mirrored about the given plane
        ! The plane number is the component index which is normal to the plane (i.e. 1: yz plane, etc.)

        implicit none

        real,dimension(3),intent(in) :: vec
        integer,intent(in) :: plane
        real,dimension(3) :: mirrored_vec

        mirrored_vec = vec
        mirrored_vec(plane) = -vec(plane)
        
    end function mirror_about_plane

    
end module helpers_mod