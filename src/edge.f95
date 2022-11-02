module edge_mod

    implicit none

    type edge
        ! A mesh edge

        integer,dimension(2) :: top_verts ! Indices of the end vertices in the mesh vertex array belonging to the top panel
        integer,dimension(2) :: bot_verts ! Indices of the end vertices in the mesh vertex array belonging to the bottom panel
        integer,dimension(2) :: panels ! Indices of the top and bottom (an odd thing to call these) panels for this edge
        integer :: i_midpoint
        integer,dimension(2) :: edge_index_for_panel ! Index of the edge which this is for each panel; edge should proceed counterclockwise for the top panel
        logical :: on_mirror_plane ! Whether this edge lies on the mirror plane
        logical :: sheds_wake ! Whether this edge sheds a wake
        real :: l ! Length

        contains

            procedure :: init => edge_init
            procedure :: get_opposing_panel => edge_get_opposing_panel
            procedure :: touches_vertex => edge_touches_vertex

    end type edge

    
contains


    subroutine edge_init(this, i1, i2, top_panel, bottom_panel, l)

        implicit none

        class(edge),intent(inout) :: this
        integer,intent(in) :: i1, i2
        integer,intent(in) :: top_panel, bottom_panel
        real,intent(in) :: l

        ! Store indices
        this%top_verts(1) = i1
        this%top_verts(2) = i2

        ! Store panels
        this%panels(1) = top_panel
        this%panels(2) = bottom_panel

        ! Store length
        this%l = l

        ! Set defaults
        this%on_mirror_plane = .false.
        this%sheds_wake = .false.
        this%bot_verts = this%top_verts
    
    end subroutine edge_init


    function edge_get_opposing_panel(this, i_panel) result(i_oppose)
        ! Returns the index of the panel opposite this one on the edge

        implicit none
        
        class(edge),intent(in) :: this
        integer,intent(in) :: i_panel

        integer :: i_oppose

        if (i_panel == this%panels(1)) then
            i_oppose = this%panels(2)
        else if (i_panel == this%panels(2)) then
            i_oppose = this%panels(1)
        else
            i_oppose = 0
        end if
        
    end function edge_get_opposing_panel


    function edge_touches_vertex(this, i_vert) result(touches)
        ! Checks whether the edge touches the given vertex

        implicit none
        
        class(edge),intent(in) :: this
        integer,intent(in) :: i_vert

        logical :: touches

        touches = this%top_verts(1) == i_vert .or. this%top_verts(2) == i_vert
        
    end function edge_touches_vertex

    
end module edge_mod