! Types for geometric objects
module quad_panel_mod

    use linked_list_mod
    use vertex_mod
    use base_panel_mod

    implicit none


    type,extends(base_panel) :: quad_panel
        ! A panel with 4 sides

        contains

            procedure :: init => quad_panel_init

    end type quad_panel

    
contains


    subroutine quad_panel_init(this, v1, v2, v3, v4)
        ! Initializes a 4-panel

        implicit none

        class(quad_panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3, v4
        
        ! Set number of sides
        this%N = 4

        ! Allocate vertex array
        allocate(this%vertices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3
        this%vertices(4)%ptr => v4

        ! Calculate normal vector

        ! Calculate area

    end subroutine quad_panel_init
    
end module quad_panel_mod