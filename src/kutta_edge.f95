module wake_edge_mod

    use vertex_mod

    implicit none

    type wake_edge
        ! A mesh edge from which a wake is shed

        integer :: i1, i2 ! Indices of the end vertices in the mesh vertex array
        integer :: top_panel, bottom_panel ! Indices of the top and bottom panels for this Kutta vertex
        logical :: on_mirror_plane ! Whether this edge lies on the mirror plane

        contains

            procedure :: init => wake_edge_init

    end type wake_edge

    
contains


    subroutine wake_edge_init(this, i1, i2, top_panel, bottom_panel)

        implicit none

        class(wake_edge),intent(inout) :: this
        integer,intent(in) :: i1, i2
        integer,intent(in) :: top_panel, bottom_panel

        ! Store indices
        this%i1 = i1
        this%i2 = i2

        ! Store panels
        this%top_panel = top_panel
        this%bottom_panel = bottom_panel
    
    end subroutine wake_edge_init
    
end module wake_edge_mod