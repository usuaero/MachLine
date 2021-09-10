! Types for geometric objects
module tri_panel_mod

    use linked_list_mod
    use base_panel_mod

    implicit none


    type,extends(base_panel) :: tri_panel
        ! A panel with 3 sides

        contains

            procedure :: init => tri_panel_init

    end type tri_panel

    
contains


    subroutine tri_panel_init(this, v1, v2, v3)
        ! Initializes a 3-panel

        implicit none

        class(tri_panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3

        ! Set number of sides
        this%N = 3

        ! Allocate vertex array
        allocate(this%vertices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3

        ! Calculate normal vec

        ! Calculate area

    end subroutine tri_panel_init
    
end module tri_panel_mod