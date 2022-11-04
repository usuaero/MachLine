module wake_strip_mod

    use panel_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod

    implicit none

    type, extends(mesh) :: wake_strip

        integer :: i_top_parent_1, i_top_parent_2, i_bot_parent_1, i_bot_parent_2
        integer :: i_top_parent, i_bot_parent
        integer :: i_top_parent_mid, i_bot_parent_mid
        logical :: on_mirror_plane

        contains

            procedure :: init => wake_strip_init
            procedure :: init_vertices => wake_strip_init_vertices
            procedure :: init_panels => wake_strip_init_panels
            procedure :: init_panel => wake_strip_init_panel
            procedure :: init_midpoints => wake_strip_init_midpoints

    end type wake_strip
    
contains


    subroutine wake_strip_init(this, freestream, starting_edge, mirror_start, mirror_plane, &
                               N_panels_streamwise, trefftz_dist, body_verts, wake_mirrored, initial_panel_order, N_body_panels)
        ! Initializes this wake strip based on the provided info

        implicit none

        class(wake_strip),intent(inout) :: this
        type(flow),intent(in) :: freestream
        type(edge),intent(in) :: starting_edge
        logical,intent(in) :: mirror_start
        integer,intent(in) :: mirror_plane, N_panels_streamwise, initial_panel_order, N_body_panels
        real,intent(in) :: trefftz_dist
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        logical,intent(in) :: wake_mirrored

        real,dimension(3) :: start_1, start_2
        integer :: N_body_verts, i

        ! Get number of vertices on the body
        N_body_verts = size(body_verts)

        ! Set mirroring
        this%on_mirror_plane = starting_edge%on_mirror_plane
        this%mirror_plane = mirror_plane
        this%mirrored = wake_mirrored .and. .not. this%on_mirror_plane

        ! Get starting locations and parent vertices
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

            ! Get parent panels
            this%i_top_parent = starting_edge%panels(1) + N_body_panels
            this%i_top_parent = starting_edge%panels(2) + N_body_panels

        else

            ! Get starting location
            start_1 = body_verts(starting_edge%top_verts(1))%loc
            start_2 = body_verts(starting_edge%top_verts(2))%loc

            ! Get parent vertices
            this%i_top_parent_1 = starting_edge%top_verts(1)
            this%i_top_parent_2 = starting_edge%top_verts(2)
            this%i_bot_parent_1 = starting_edge%bot_verts(1)
            this%i_bot_parent_2 = starting_edge%bot_verts(2)

            ! Get parent panels
            this%i_top_parent = starting_edge%panels(1)
            this%i_top_parent = starting_edge%panels(2)

        end if

        ! Initialize vertices
        call this%init_vertices(freestream, mirror_start, N_panels_streamwise, trefftz_dist, start_1, start_2, body_verts)

        ! Intialize panels
        call this%init_panels(N_panels_streamwise)

        ! Initialize other panel properties
        do i=1,this%N_panels
            call this%panels(i)%init_with_flow(freestream, this%mirrored, mirror_plane)
            call this%panels(i)%set_distribution(initial_panel_order, this%panels, mirrored=this%mirrored)
        end do

    end subroutine wake_strip_init


    subroutine wake_strip_init_vertices(this, freestream, mirror_start, N_panels_streamwise, trefftz_dist, &
                                        start_1, start_2, body_verts)
        ! Initializes this wake strip's vertices based on the provided info

        implicit none

        class(wake_strip),intent(inout) :: this
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_start
        integer,intent(in) :: N_panels_streamwise
        real,intent(in) :: trefftz_dist
        real,dimension(3),intent(in) :: start_1, start_2
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts

        real,dimension(3) :: loc_1, loc_2, loc
        real :: d1, d2, sep_1, sep_2
        integer :: i, N_body_verts

        ! Allocate memory
        this%N_verts = N_panels_streamwise*2 + 2
        N_body_verts = size(body_verts)
        allocate(this%vertices(this%N_verts))

        ! Initialize starting vertices in wake strip
        call this%vertices(1)%init(start_1, 1, 1)
        call this%vertices(2)%init(start_2, 2, 1)

        ! Set parents
        this%vertices(1)%top_parent = this%i_top_parent_1
        this%vertices(1)%bot_parent = this%i_bot_parent_1
        this%vertices(2)%top_parent = this%i_top_parent_2
        this%vertices(2)%bot_parent = this%i_bot_parent_2

        ! Calculate distances to Trefftz plane
        d1 = trefftz_dist - inner(start_1, freestream%c_hat_g)
        d2 = trefftz_dist - inner(start_2, freestream%c_hat_g)

        ! Calculate spacing between vertices
        sep_1 = d1 / N_panels_streamwise
        sep_2 = d2 / N_panels_streamwise

        ! Loop through following vertices
        do i=3,this%N_verts

            ! Calculate location of vertices
            if (modulo(i, 2) == 0) then
                loc = start_2 + sep_2*(i-1)/2*freestream%c_hat_g
            else
                loc = start_1 + sep_1*i/2*freestream%c_hat_g
            end if

            ! Initialize vertices
            call this%vertices(i)%init(loc, i, 1)

            ! Set parents
            if (modulo(i, 2) == 0) then
                this%vertices(i)%top_parent = this%i_top_parent_2
                this%vertices(i)%bot_parent = this%i_bot_parent_2
            else
                this%vertices(i)%top_parent = this%i_top_parent_1
                this%vertices(i)%bot_parent = this%i_bot_parent_1
            end if

        end do

    end subroutine wake_strip_init_vertices


    subroutine wake_strip_init_panels(this, N_panels_streamwise)
        ! Initializes this wake strip's panels based on the provided info

        implicit none

        class(wake_strip),intent(inout) :: this
        integer,intent(in) :: N_panels_streamwise

        real :: d1, d2
        integer :: i, j, i1, i2, advance, N_skipped
        logical,dimension(:),allocatable :: skipped_panels
        type(panel),dimension(:),allocatable :: temp_panels

        ! Determine number of panels
        this%N_panels = N_panels_streamwise*2
        allocate(this%panels(this%N_panels))
        allocate(skipped_panels(this%N_panels), source=.false.)

        ! Create panels
        i1 = 1
        i2 = 2
        do i=1,this%N_panels

            ! Initialize
            advance = 0

            ! See if one is at the end
            if (i1 == this%N_verts - 1) then
                advance = 2
            else if (i2 == this%N_verts) then
                advance = 1
            else

                ! Check lengths of hypotenuses
                d1 = dist(this%vertices(i1+2)%loc, this%vertices(i2)%loc)
                d2 = dist(this%vertices(i1)%loc, this%vertices(i2+2)%loc)

                ! Pick which one
                if (d1 < d2) then
                    advance = 1
                else
                    advance = 2
                end if
            end if

            ! Advance
            if (advance == 1) then

                ! Initialize
                call this%init_panel(i, i1, i1+2, i2, skipped_panels)
                
                ! Increment index
                i1 = i1 + 2
            else
                
                ! Initialize
                call this%init_panel(i, i1, i2+2, i2, skipped_panels)

                ! Increment index
                i2 = i2 + 2
            end if

        end do

        ! Get rid of skipped panels
        N_skipped = count(skipped_panels)
        allocate(temp_panels(this%N_panels - N_skipped))

        ! Move non-skipped panels over
        j = 0
        do i=1,this%N_panels

            ! Check that this panel was not skipped
            if (.not. skipped_panels(i)) then

                ! Update index
                j = j + 1

                ! Copy over
                temp_panels(j) = this%panels(i)

            end if

        end do

        ! Move allocation
        call move_alloc(temp_panels, this%panels)

        ! Update number of panels
        this%N_panels = this%N_panels - N_skipped

    end subroutine wake_strip_init_panels


    subroutine wake_strip_init_panel(this, i_panel, i1, i2, i3, skipped_panels)
        ! Initializes the specified panel

        implicit none
        
        class(wake_strip),intent(inout) :: this
        integer,intent(in) :: i_panel, i1, i2, i3
        logical,dimension(:),allocatable,intent(inout) :: skipped_panels

        ! Check for zero area
        if (this%has_zero_area(i1, i2, i3)) then
            skipped_panels(i_panel) = .true.
        else
            call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel)
        end if

        ! Set that it's in the wake
        this%panels(i_panel)%in_wake = .true.
    
    end subroutine wake_strip_init_panel


    subroutine wake_strip_init_midpoints(this)
        ! Initializes the mesh midpoints

        implicit none
        
        class(wake_strip),target,intent(inout) :: this

        integer :: N_mids, i_mid, i, edge_prev
        real,dimension(3) :: loc

        ! Determine how many midpoints we need
        N_mids = 2*this%N_panels + 1

        ! Allocate more space
        call this%allocate_new_vertices(N_mids)

        ! Loop through panels (since we know how the mesh is basically structured)
        i_mid = this%N_verts - N_mids
        do i=1,this%N_panels

            ! Create midpoint at the base of this panel
            i_mid = i_mid + 1

            ! Get location
            loc = 0.5*(this%panels(i)%get_vertex_loc(1) + this%panels(i)%get_vertex_loc(3))

            ! Initialize
            call this%vertices(i_mid)%init(loc, i_mid, 2)

            ! Store adjacent panels
            call this%vertices(i_mid)%panels%append(i)
            if (i > 1) call this%vertices(i_mid)%panels%append(i-1)

            ! Store adjacent vertices
            call this%vertices(i_mid)%adjacent_vertices%append(this%panels(i)%get_vertex_index(1))
            call this%vertices(i_mid)%adjacent_vertices%append(this%panels(i)%get_vertex_index(3))

            ! Point panel to it
            this%panels(i)%midpoints(3)%ptr => this%vertices(i_mid)

            ! Point previous panel to it
            if (i > 1) then
                if (edge_prev == 1) then
                    this%panels(i-1)%midpoints(2)%ptr => this%vertices(i_mid)
                else
                    this%panels(i-1)%midpoints(1)%ptr => this%vertices(i_mid)
                end if
            end if

            ! Store parents
            this%vertices(i_mid)%top_parent = this%i_top_parent_mid
            this%vertices(i_mid)%bot_parent = this%i_bot_parent_mid

            ! Initialize side midpoint
            i_mid = i_mid + 1
            
            ! Check which side we're on
            ! Side 1
            if (this%panels(i)%get_vertex_index(2) - this%panels(i)%get_vertex_index(1) == 2) then

                ! Get location
                loc = 0.5*(this%panels(i)%get_vertex_loc(1) + this%panels(i)%get_vertex_loc(2))

                ! Initialize
                call this%vertices(i_mid)%init(loc, i_mid, 2)

                ! Store adjacent panels
                call this%vertices(i_mid)%panels%append(i)

                ! Store adjacent vertices
                call this%vertices(i_mid)%adjacent_vertices%append(this%panels(i)%get_vertex_index(1))
                call this%vertices(i_mid)%adjacent_vertices%append(this%panels(i)%get_vertex_index(2))

                ! Point panel to it
                this%panels(i)%midpoints(1)%ptr => this%vertices(i_mid)
                edge_prev = 1

                ! Store parents
                this%vertices(i_mid)%top_parent = this%i_top_parent_1
                this%vertices(i_mid)%bot_parent = this%i_bot_parent_1

            ! Side 2
            else

                ! Get location
                loc = 0.5*(this%panels(i)%get_vertex_loc(3) + this%panels(i)%get_vertex_loc(2))

                ! Initialize
                call this%vertices(i_mid)%init(loc, i_mid, 2)

                ! Store adjacent panels
                call this%vertices(i_mid)%panels%append(i)

                ! Store adjacent vertices
                call this%vertices(i_mid)%adjacent_vertices%append(this%panels(i)%get_vertex_index(3))
                call this%vertices(i_mid)%adjacent_vertices%append(this%panels(i)%get_vertex_index(2))

                ! Point panel to it
                this%panels(i)%midpoints(2)%ptr => this%vertices(i_mid)
                edge_prev = 2

                ! Store parents
                this%vertices(i_mid)%top_parent = this%i_top_parent_2
                this%vertices(i_mid)%bot_parent = this%i_bot_parent_2

            end if

        end do

        ! Initialize final midpoint
        i_mid = i_mid + 1

        ! Get location
        loc = 0.5*(this%vertices(this%N_verts-N_mids)%loc + this%vertices(this%N_verts-1-N_mids)%loc)

        ! Initialize
        call this%vertices(i_mid)%init(loc, i_mid, 2)

        ! Store adjacent panels
        call this%vertices(i_mid)%panels%append(this%N_panels)

        ! Store adjacent vertices
        call this%vertices(i_mid)%adjacent_vertices%append(this%N_verts-1-N_mids)
        call this%vertices(i_mid)%adjacent_vertices%append(this%N_verts-N_mids)

        ! Point previous panel to it
        if (edge_prev == 1) then
            this%panels(this%N_panels)%midpoints(2)%ptr => this%vertices(i_mid)
        else
            this%panels(this%N_panels)%midpoints(1)%ptr => this%vertices(i_mid)
        end if

        ! Store parents
        this%vertices(i_mid)%top_parent = this%i_top_parent_mid
        this%vertices(i_mid)%bot_parent = this%i_bot_parent_mid

        this%midpoints_created = .true.
        
    end subroutine wake_strip_init_midpoints
    
end module wake_strip_mod