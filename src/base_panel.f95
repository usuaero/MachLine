! Types for geometric objects
module base_panel_mod

    use linked_list_mod
    use vertex_mod

    implicit none


    type base_panel
        ! A panel with an arbitrary number of sides

        integer :: N ! Number of sides/vertices
        type(vertex_pointer),dimension(:),allocatable :: vertices
        real,dimension(3) :: n_hat ! Normal vector
        real :: A ! Surface area

        contains

            procedure :: calc_area => base_panel_calc_area
            procedure :: calc_normal => base_panel_calc_area

    end type base_panel


    type panel_pointer
        ! A pointer to a base panel, for creating panel arrays (polymorphism!)

        class(base_panel),pointer :: ptr

    end type panel_pointer

    
contains

    subroutine base_panel_calc_area(this)

        implicit none

        class(base_panel),intent(inout) :: this

    end subroutine base_panel_calc_area


    subroutine base_panel_calc_normal(this)

        implicit none

        class(base_panel),intent(inout) :: this

    end subroutine base_panel_calc_normal


    
end module base_panel_mod