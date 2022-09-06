module wake_strip_mod

    use vertex_mod
    use panel_mod
    use vertex_mod
    use edge_mod
    use helpers_mod

    implicit none

    type wake_strip

        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        integer :: N_verts, N_panels = 0
        integer :: i_top_parent_1, i_top_parent_2, i_bot_parent_1, i_bot_parent_2

        contains

            procedure :: init => wake_strip_init
            procedure :: init_vertices => wake_strip_init_vertices
            procedure :: init_panels => wake_strip_init_panels

    end type wake_strip
    
contains

subroutine wake_strip_init(this, freestream, starting_edge, mirror_start, mirror_plane, &
                           N_panels_streamwise, trefftz_dist, body_verts)
    ! Initializes this wake strip based on the provided info

    implicit none
    
    class(wake_strip),intent(inout) :: this
    type(flow),intent(in) :: freestream
    type(edge),intent(in) :: starting_edge
    logical,intent(in) :: mirror_start
    integer,intent(in) :: mirror_plane, N_panels_streamwise
    real,intent(in) :: trefftz_dist
    type(vertex),dimension(:),allocatable,intent(in) :: body_verts

    real,dimension(3) :: start_1, start_2
    integer :: i

    ! Determine number of panels and vertices
    N_panels = N_panels_streamwise*2
    N_verts = N_panels + 2

    ! Allocate memory
    allocate(this%panels(N_panels))
    allocate(this%vertices(2,N_panels_streamwise+1)) ! The first dimension tells which side of the strip the vertex is on

    ! Get starting location
    if (mirror_start) then ! Note order is flipped here to maintain convention
        start_1 = mirror_across_plane(body_verts(starting_edge%vertices(2))%loc, mirror_plane)
        start_2 = mirror_across_plane(body_verts(starting_edge%vertices(1))%loc, mirror_plane)
    else
        start_1 = body_verts(starting_edge%vertices(1))%loc
        start_2 = body_verts(starting_edge%vertices(2))%loc
    end if

    ! Initialize starting vertices in wake strip
    call this%vertices(1,1)%init(start_1, 0, 1) ! We don't need to know the indices here
    call this%vertices(2,1)%init(start_2, 0, 1) ! We don't need to know the indices here

    ! Loop through following vertices
    do i=2,N_panels_streamwise
    end do
    
end subroutine wake_strip_init
    
end module wake_strip_mod