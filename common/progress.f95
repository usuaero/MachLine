module progress_mod

    implicit none

    type progress_bar

        integer :: N_total, N_complete

        contains

            procedure :: init => progress_bar_init
            procedure :: advance => progress_bar_advance

    end type progress_bar
    
contains

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
    
end module progress_mod