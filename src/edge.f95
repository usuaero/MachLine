module edge_mod

    implicit none

    type edge
        ! A mesh edge

        integer,dimension(2) :: verts ! Indices of the end vertices in the mesh vertex array
        integer,dimension(2) :: panels ! Indices of the top and bottom panels for this Kutta vertex
        logical :: on_mirror_plane ! Whether this edge lies on the mirror plane
        logical :: sheds_wake ! Whether this edge sheds a wake

        contains

            procedure :: init => edge_init

    end type edge

    
contains


    subroutine edge_init(this, i1, i2, top_panel, bottom_panel)

        implicit none

        class(edge),intent(inout) :: this
        integer,intent(in) :: i1, i2
        integer,intent(in) :: top_panel, bottom_panel

        ! Store indices
        this%verts(1) = i1
        this%verts(2) = i2

        ! Store panels
        this%panels(1) = top_panel
        this%panels(2) = bottom_panel
    
    end subroutine edge_init
    
end module edge_mod