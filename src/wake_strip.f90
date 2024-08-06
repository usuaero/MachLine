module wake_strip_mod

    use panel_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod

    implicit none

    type, extends(mesh) :: wake_strip

        integer :: i_top_parent_1, i_top_parent_2, i_bot_parent_1, i_bot_parent_2
        integer :: i_top_parent, i_bot_parent
        logical :: on_mirror_plane

        logical :: calc_adjoint

        contains

            procedure :: init => wake_strip_init
            procedure :: init_vertices => wake_strip_init_vertices
            procedure :: init_panels => wake_strip_init_panels
            procedure :: init_panel => wake_strip_init_panel

    end type wake_strip
    
contains


    subroutine wake_strip_init(this, freestream, starting_edge, mirror_start, mirror_plane, &
                               N_panels_streamwise, trefftz_dist, body_verts, wake_mirrored, &
                               initial_panel_order, N_body_panels, calc_adjoint)
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
        logical,intent(inout),optional :: calc_adjoint

        real,dimension(3) :: start_1, start_2
        integer :: N_body_verts, i

        logical :: adjoint
        type(sparse_matrix) :: d_start_1, d_start_2

        ! see if calc_adjoint
        this%calc_adjoint = calc_adjoint

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
            this%i_bot_parent = starting_edge%panels(2) + N_body_panels

        else

            ! Get starting location
            start_1 = body_verts(starting_edge%top_verts(1))%loc
            start_2 = body_verts(starting_edge%top_verts(2))%loc

            if (this%calc_adjoint) then
                d_start_1 = body_verts(starting_edge%top_verts(1))%d_loc
                d_start_2 = body_verts(starting_edge%top_verts(2))%d_loc
            end if

            ! Get parent vertices
            this%i_top_parent_1 = starting_edge%top_verts(1)
            this%i_top_parent_2 = starting_edge%top_verts(2)
            this%i_bot_parent_1 = starting_edge%bot_verts(1)
            this%i_bot_parent_2 = starting_edge%bot_verts(2)

            ! Get parent panels
            this%i_top_parent = starting_edge%panels(1)
            this%i_bot_parent = starting_edge%panels(2)

        end if

        ! Initialize vertices
        if (this%calc_adjoint) then
            call this%init_vertices(freestream, N_panels_streamwise, trefftz_dist, start_1, start_2,&
            body_verts, d_start_1, d_start_2)
        else
            call this%init_vertices(freestream, N_panels_streamwise, trefftz_dist, start_1, start_2,&
            body_verts)
        end if

        
        ! Intialize panels
        call this%init_panels(N_panels_streamwise)
        
        ! Initialize other panel properties
        do i=1,this%N_panels
            call this%panels(i)%init_with_flow(freestream, this%mirrored, mirror_plane)
            
            call this%panels(i)%set_distribution(1, this%panels, this%vertices, this%mirrored) ! With the current formulation, wake panels are always linear
            
            if (this%calc_adjoint) call this%panels(i)%init_with_flow_adjoint(freestream)
        end do
        
    end subroutine wake_strip_init


    subroutine wake_strip_init_vertices(this, freestream, N_panels_streamwise, trefftz_dist, &
                                start_1, start_2, body_verts, d_start_1, d_start_2)
        ! Initializes this wake strip's vertices based on the provided info

        implicit none

        class(wake_strip),intent(inout) :: this
        type(flow),intent(in) :: freestream
        integer,intent(in) :: N_panels_streamwise
        real,intent(in) :: trefftz_dist
        real,dimension(3),intent(in) :: start_1, start_2
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        
        type(sparse_matrix),intent(inout),optional :: d_start_1, d_start_2

        real,dimension(3) :: loc
        real :: d1, d2, sep_1, sep_2
        integer :: i, N_original_verts

        type(sparse_vector) :: d_d1, d_d2, d_sep_1, d_sep_2
        type(sparse_matrix) :: d_loc


        ! Allocate memory
        this%N_verts = N_panels_streamwise*2 + 2
        
        ! Get number of original vertices of mesh that was read in (before vertex cloning)
        ! this is used to initialize the adjoint size
        N_original_verts = body_verts(1)%d_loc%full_num_cols/3

        allocate(this%vertices(this%N_verts))

        ! Initialize starting vertices in wake strip
        call this%vertices(1)%init(start_1, 1)
        call this%vertices(2)%init(start_2, 2)

        ! Set parents
        this%vertices(1)%top_parent = this%i_top_parent_1
        this%vertices(1)%bot_parent = this%i_bot_parent_1
        this%vertices(2)%top_parent = this%i_top_parent_2
        this%vertices(2)%bot_parent = this%i_bot_parent_2

        if (this%calc_adjoint) then
            ! I think this should be set to N_original verts
            call this%vertices(1)%init_adjoint(N_original_verts, wake_vertex = .true.)
            call this%vertices(2)%init_adjoint(N_original_verts, wake_vertex = .true.)
        end if

        ! Calculate distances to Trefftz plane
        d1 = trefftz_dist - inner(start_1, freestream%c_hat_g)
        d2 = trefftz_dist - inner(start_2, freestream%c_hat_g)

        
        ! Calculate spacing between vertices
        sep_1 = d1 / N_panels_streamwise
        sep_2 = d2 / N_panels_streamwise
        
        if (this%calc_adjoint) then
            d_d1 = d_start_1%broadcast_vector_dot_element(-freestream%c_hat_g)
            d_d2 = d_start_2%broadcast_vector_dot_element(-freestream%c_hat_g)

            call d_sep_1%init_from_sparse_vector(d_d1)
            call d_sep_1%broadcast_element_times_scalar(1./N_panels_streamwise)

            call d_sep_2%init_from_sparse_vector(d_d2)
            call d_sep_2%broadcast_element_times_scalar(1./N_panels_streamwise)

        end if

        ! Loop through following vertices
        do i=3,this%N_verts

            ! Calculate location of vertices
            if (modulo(i, 2) == 0) then

                ! put the location of the vertex evenly spaced along the panel edge (i think)
                loc = start_2 + sep_2*(i-2)/2*freestream%c_hat_g

                ! if calc_adjoint, then calc d_loc
                if (this%calc_adjoint) then

                    d_loc = d_sep_2%broadcast_element_times_vector((i-2)/2*freestream%c_hat_g)
                    
                    call d_loc%sparse_add(d_start_2)

                end if

            else

                ! put the location of the vertex evenly spaced along the panel edge (i think)
                loc = start_1 + sep_1*(i-1)/2*freestream%c_hat_g

                ! if calc_adjoint, then calc d_loc
                if (this%calc_adjoint) then
                    d_loc = d_sep_1%broadcast_element_times_vector((i-1)/2*freestream%c_hat_g)
                    call d_loc%sparse_add(d_start_1)
                end if

            end if

            ! Initialize vertices
            call this%vertices(i)%init(loc, i)

            ! Set parents
            if (modulo(i, 2) == 0) then
                this%vertices(i)%top_parent = this%i_top_parent_2
                this%vertices(i)%bot_parent = this%i_bot_parent_2
            else
                this%vertices(i)%top_parent = this%i_top_parent_1
                this%vertices(i)%bot_parent = this%i_bot_parent_1
            end if

            ! if calc_adjoint, set vertex attribute d_loc as the d_loc calculated in this do loop
            if (this%calc_adjoint) then
                call this%vertices(i)%d_loc%init_from_sparse_matrix(d_loc)
                
                ! deallocate d_loc columns for next loop iteration
                deallocate(d_loc%columns)
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

                ! Check lengths of hypotenuses (dist = norm2(a-b))
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
            call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel, in_wake=.true.)
        
            ! if calc adjoint, init panel adjoint
            if (this%calc_adjoint) then
                call this%panels(i_panel)%init_adjoint()
            end if

        end if

    end subroutine wake_strip_init_panel

    
end module wake_strip_mod
