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

    real,dimension(3) :: start_1, start_2, loc_1, loc_2
    real :: d1, d2, sep_1, sep_2
    integer :: i, N_body_verts, i1, i2, advance

    ! Determine number of panels and vertices
    N_panels = N_panels_streamwise*2
    N_verts = N_panels + 2
    N_body_verts = size(body_verts)

    ! Allocate memory
    allocate(this%panels(N_panels))
    allocate(this%vertices(2,N_panels_streamwise+1)) ! The first dimension tells which side of the strip the vertex is on

    if (mirror_start) then

        ! Get starting location
        ! Note order is flipped here to maintain convention
        start_1 = mirror_across_plane(body_verts(starting_edge%top_verts(2))%loc, mirror_plane)
        start_2 = mirror_across_plane(body_verts(starting_edge%top_verts(1))%loc, mirror_plane)

        ! Get parent vertices
        this%i_top_parent_1 = starting_edge%top_verts(2) + N_body_verts
        this%i_top_parent_2 = starting_edge%top_verts(1) + N_body_verts
        this%i_bot_parent_1 = starting_edge%bot_verts(2) + N_body_verts
        this%i_bot_parent_2 = starting_edge%bot_verts(1) + N_body_verts

    else

        ! Get starting location
        start_1 = body_verts(starting_edge%top_verts(1))%loc
        start_2 = body_verts(starting_edge%top_verts(2))%loc

        ! Get parent vertices
        this%i_top_parent_1 = starting_edge%top_verts(1)
        this%i_top_parent_2 = starting_edge%top_verts(2)
        this%i_bot_parent_1 = starting_edge%bot_verts(1)
        this%i_bot_parent_2 = starting_edge%bot_verts(2)

    end if

    ! Initialize starting vertices in wake strip
    call this%vertices(1,1)%init(start_1, 0, 1) ! We don't need to know the indices here
    call this%vertices(2,1)%init(start_2, 0, 1) ! We don't need to know the indices here

    ! Calculate distances to Trefftz plane
    d1 = trefftz_dist - inner(start_1, freestream%c_hat_g)
    d2 = trefftz_dist - inner(start_2, freestream%c_hat_g)

    ! Calculate spacing between vertices
    sep_1 = d1/N_panels_streamwise
    sep_2 = d2/N_panels_streamwise

    ! Loop through following vertices
    do i=2,N_panels_streamwise+1

        ! Calculate location of vertices
        loc_1 = start_1 + sep_1*(i-1)*freestream%c_hat_g
        loc_2 = start_2 + sep_2*(i-1)*freestream%c_hat_g

        ! Initialize vertices
        call this%vertices(1,i)%init(loc_1, 0, 1)
        call this%vertices(2,i)%init(loc_2, 0, 1)

        ! Set parents
        if (mirror_start) then
            this%vertices(1,i)%top_parent = i_top_parent_1 + N_body_verts
            this%vertices(1,i)%bot_parent = i_bot_parent_1 + N_body_verts
            this%vertices(2,i)%top_parent = i_top_parent_2 + N_body_verts
            this%vertices(2,i)%bot_parent = i_bot_parent_2 + N_body_verts
        else
            this%vertices(1,i)%top_parent = i_top_parent_1
            this%vertices(1,i)%bot_parent = i_bot_parent_1
            this%vertices(2,i)%top_parent = i_top_parent_2
            this%vertices(2,i)%bot_parent = i_bot_parent_2
        end if

    end do

    ! Create panels
    i1 = 1
    i2 = 1
    do i=1,this%N_panels
        
        ! Initialize
        advance = 0

        ! See if one is at the end
        if (i1 == N_panels_streamwise+1) then
            advance = 2
        else if (i2 == N_panels_streamwise+1) then
            advance = 1
        else

            ! Check lengths of hypotenuses
            d1 = dist(this%vertices(1,i1+1)%loc, this%vertices(2,i2)%loc)
            d2 = dist(this%vertices(1,i1)%loc, this%vertices(2,i2+1)%loc)

            ! Pick which one
            if (d1 < d2) then
                advance = 1
            else
                advance = 2
            end if
        end if

        ! Advance
        if (advance == 1) then
            call this%panels(i)%init(this%vertices(1,i1), this%vertices(1,i1+1), this%vertices(2,i2), i)
            i1 = i1 + 1
        else
            call this%panels(i)%init(this%vertices(1,i1), this%vertices(2,i2+1), this%vertices(2,i2), i)
            i2 = i2 + 1
        end if

    end do
    
end subroutine wake_strip_init
    
end module wake_strip_mod