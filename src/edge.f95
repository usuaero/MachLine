module edge_mod

    implicit none

    type edge
        ! A mesh edge

        integer,dimension(2) :: verts ! Indices of the end vertices in the mesh vertex array (ignoring clones)
        integer,dimension(2) :: panels ! Indices of the top and bottom (an odd thing to call these) panels for this edge
        integer :: i_midpoint
        integer,dimension(2) :: edge_index_for_panel ! Index of the edge which this is for each panel; edge should proceed counterclockwise for the top panel
        logical :: on_mirror_plane ! Whether this edge lies on the mirror plane
        logical :: sheds_wake ! Whether this edge sheds a wake
        logical :: discontinuous ! Whether this edge has a jump in doublet strength
        real :: l ! Length

        contains

            procedure :: init => edge_init

    end type edge

    
contains


    subroutine edge_init(this, i1, i2, top_panel, bottom_panel, l)

        implicit none

        class(edge),intent(inout) :: this
        integer,intent(in) :: i1, i2
        integer,intent(in) :: top_panel, bottom_panel
        real,intent(in) :: l

        ! Store indices
        this%verts(1) = i1
        this%verts(2) = i2

        ! Store panels
        this%panels(1) = top_panel
        this%panels(2) = bottom_panel

        ! Store length
        this%l = l

        ! Set defaults
        this%on_mirror_plane = .false.
        this%sheds_wake = .false.
        this%discontinuous = .false.
    
    end subroutine edge_init

    
end module edge_mod