module helpers_mod

    implicit none

    type progress_bar
        ! NOT WORKING

        integer :: N_total, N_complete

        contains

            procedure :: init => progress_bar_init
            procedure :: advance => progress_bar_advance

    end type progress_bar
    
contains

    subroutine check_allocation(stat, what)

        implicit none

        integer,intent(in) :: stat
        character(len=*),intent(in),optional :: what

        ! Check if stat is nonzero
        if (stat /= 0) then
            if (present(what)) then
                write(*,*) "!!! Your computer has insufficient memory to allocate ", what, ". Quitting..."
            else
                write(*,*) "!!! Your computer has insufficient memory. Quitting..."
            end if
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


    subroutine progress_bar_init(this, N_total)
        ! Sets up the progress bar

        implicit none

        class(progress_bar),intent(inout) :: this
        integer,intent(in) :: N_total

        ! Initialize
        this%N_total = N_total
        this%N_complete = 0

        ! Print out initial
        write(*,'(i3, a)',advance='no') 0, "% complete."
    
    end subroutine progress_bar_init


    subroutine progress_bar_advance(this)
        ! Advances the progress bar by one

        implicit none

        class(progress_bar),intent(inout) :: this

        ! Update number complete
        this%N_complete = this%N_complete + 1

        ! Print progress
        write(*,'(tl13, i3, a)',advance='no') 100.0*real(this%N_complete)/real(this%N_total), "% complete."

        ! Write newline
        if (this%N_complete == this%N_total) then
            write(*,*)
        end if
    
    end subroutine progress_bar_advance


    subroutine delete_file(filename)
        ! Deletes the given file if it exists

        implicit none

        character(len=:),allocatable,intent(in) :: filename

        integer :: stat

        open(unit=13, iostat=stat, file=filename)
        if (stat == 0) close(13, status='delete')
    
    end subroutine delete_file

    
end module helpers_mod