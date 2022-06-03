! Class for modeling wake meshes
module wake_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use linked_list_mod
    use helpers_mod
    use vertex_mod
    use panel_mod
    use math_mod
    use flow_mod
    use edge_mod

    implicit none


    type wake_mesh

        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        integer :: N_verts, N_panels

        contains

            procedure :: init => wake_mesh_init
            procedure :: init_vertices => wake_mesh_init_vertices
            procedure :: init_panels => wake_mesh_init_panels
            procedure :: init_midpoints => wake_mesh_init_midpoints

    end type wake_mesh


contains


    subroutine wake_mesh_init(this, body_edges, body_verts, N_body_panels, freestream, asym_flow, mirror_plane, &
                              N_panels_streamwise, trefftz_dist)
        ! Creates the vertices and panels. Handles vertex association.

        implicit none

        class(wake_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        integer,intent(in) :: N_body_panels
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_panels_streamwise
        real,intent(in) :: trefftz_dist

        integer :: i, j
        integer :: N_wake_edge_verts, N_wake_edges, N_mids
        integer,dimension(:),allocatable :: wake_edge_indices, wake_edge_verts
        logical,dimension(:),allocatable :: is_wake_edge_vertex

        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Initializing wake with a Trefftz distance of ", trefftz_dist, "..."

        ! Count up wake-shedding edges
        N_wake_edges = 0
        do i=1,size(body_edges)
            if (body_edges(i)%sheds_wake) N_wake_edges = N_wake_edges + 1
        end do

        ! Get indices of wake-shedding edges and vertices
        allocate(wake_edge_indices(N_wake_edges))
        allocate(is_wake_edge_vertex(size(body_verts)))
        j = 0
        do i=1,size(body_edges)
            if (body_edges(i)%sheds_wake) then

                ! Store edge index
                j = j + 1
                wake_edge_indices(j) = i

                ! Store that the vertices are wake-shedding
                is_wake_edge_vertex(body_edges(i)%verts(1)) = .true.
                is_wake_edge_vertex(body_edges(i)%verts(2)) = .true.

            end if
        end do

        ! Determine how many vertices are along the wake-shedding edges
        N_wake_edge_verts = count(is_wake_edge_vertex)

        ! Get wake-shedding edge vertex indices
        allocate(wake_edge_verts(N_wake_edge_verts))
        j = 0
        do i=1,size(body_verts)
            if (is_wake_edge_vertex(i)) then
                j = j + 1
                wake_edge_verts(j) = i
                body_verts(i)%index_in_wake_vertices = j
            end if
        end do

        ! Determine necessary number of vertices
        if (asym_flow) then
            this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)*2
        else
            this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)
        end if

        ! Determine necessary number of panels
        if (asym_flow) then
            this%N_panels = N_wake_edges*N_panels_streamwise*4
        else
            this%N_panels = N_wake_edges*N_panels_streamwise*2
        end if

        ! Determine necessary number of midpoints
        if (doublet_order == 2) then
            if (asym_flow) then
                N_mids = (N_wake_edge_verts*N_panels_streamwise + N_wake_edges*(2*N_panels_streamwise + 1) )*2
            else
                N_mids = N_wake_edge_verts*N_panels_streamwise + N_wake_edges*(2*N_panels_streamwise + 1)
            end if
        else
            N_mids = 0
        end if

        ! Initialize vertices
        allocate(this%vertices(this%N_verts + N_mids))
        call this%init_vertices(body_verts, freestream, wake_edge_verts, asym_flow, trefftz_dist, N_panels_streamwise, mirror_plane)

        ! Initialize panels
        allocate(this%panels(this%N_panels))
        call this%init_panels(body_edges, body_verts, N_panels_streamwise, wake_edge_indices, asym_flow, N_body_panels)

        ! Initialize midpoints (if needed)
        if (doublet_order == 2) then

            call this%init_midpoints(body_edges, body_verts, wake_edge_indices)
        
            ! Include midpoints in the total number of vertices
            this%N_verts = this%N_verts + N_mids

        end if

        ! Initialize freestream-dependent properties of panels once the midpoints have been created
        ! The mirror of wake panels will never need to be initialized
        do i=1,this%N_panels
            call this%panels(i)%init_with_flow(freestream, .false., mirror_plane)
            call this%panels(i)%set_influencing_verts()
        end do

        if (verbose) write(*,'(a, i7, a, i7, a)') "Done. Created ", this%N_verts, " wake vertices and ", &
                                                  this%N_panels, " wake panels."

    end subroutine wake_mesh_init


    subroutine wake_mesh_init_vertices(this, body_verts, freestream, wake_edge_verts, asym_flow, &
                                       trefftz_dist, N_panels_streamwise, mirror_plane)
        ! Initializes the vertices for the mesh

        class(wake_mesh),intent(inout) :: this
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(in) :: wake_edge_verts
        logical,intent(in) :: asym_flow
        real,intent(in) :: trefftz_dist
        integer,intent(in) :: N_panels_streamwise, mirror_plane

        integer :: i_vert, i, j, i_top_parent, i_bot_parent, i_mirrored_vert, N_wake_edge_verts
        real :: distance, vertex_separation, mirrored_distance, mirrored_vertex_separation
        real,dimension(3) :: start, loc, mirrored_start

        ! Determine vertex placement
        i_vert = 0
        N_wake_edge_verts = size(wake_edge_verts)
        do i=1,N_wake_edge_verts

            ! Check this is not a midpoint vertex
            i_top_parent = wake_edge_verts(i)
            if (body_verts(i_top_parent)%vert_type == 1) then

                ! Get indices
                i_bot_parent = body_verts(i_top_parent)%i_wake_partner

                ! Determine distance from origin to wake-shedding vertex in the direction of the freestream flow
                start = body_verts(i_top_parent)%loc
                distance = trefftz_dist - inner(start, freestream%c_hat_g)

                ! Double-check
                if (dist(start, body_verts(i_bot_parent)%loc) > 1.e-12) then
                    write(*,*) "!!! Wake edge vertices are not identical. Quitting..."
                    stop
                end if

                ! Determine vertex separation
                vertex_separation = distance/N_panels_streamwise

                ! Same for mirror
                if (asym_flow) then

                    ! Determine start location
                    mirrored_start = mirror_across_plane(start, mirror_plane)
                    mirrored_distance = trefftz_dist-inner(mirrored_start, freestream%c_hat_g)

                    ! Determine vertex separation
                    mirrored_vertex_separation = mirrored_distance/N_panels_streamwise

                end if

                ! Loop down the streamwise direction to place vertices
                do j=1,N_panels_streamwise+1

                    ! Determine location
                    i_vert = i_vert + 1
                    loc = start + vertex_separation*(j-1)*freestream%c_hat_g

                    ! Initialize vertex
                    call this%vertices(i_vert)%init(loc, i_vert, 1)

                    ! Set parent index
                    this%vertices(i_vert)%top_parent = i_top_parent
                    this%vertices(i_vert)%bot_parent = i_bot_parent

                    ! Initialize mirror
                    if (asym_flow) then

                        ! Determine location
                        i_mirrored_vert = i_vert + this%N_verts/2
                        mirrored_start = mirror_across_plane(start, mirror_plane)
                        loc = mirrored_start + mirrored_vertex_separation*(j-1)*freestream%c_hat_g

                        ! Initialize vertex
                        call this%vertices(i_mirrored_vert)%init(loc, i_mirrored_vert, 1)

                        ! Set parent index
                        this%vertices(i_mirrored_vert)%top_parent = i_top_parent + size(body_verts)
                        this%vertices(i_mirrored_vert)%bot_parent = i_bot_parent + size(body_verts)

                    end if

                end do
            end if
        end do
    
        
    end subroutine wake_mesh_init_vertices


    subroutine wake_mesh_init_panels(this, body_edges, body_verts, N_panels_streamwise, wake_edge_indices, asym_flow, N_body_panels)
        ! Initializes the wake panels

        class(wake_mesh),intent(inout) :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        integer,intent(in) :: N_panels_streamwise, N_body_panels
        integer,dimension(:),allocatable,intent(in) :: wake_edge_indices
        logical,intent(in) :: asym_flow

        integer :: i, j, k, i_start, i_stop, i_panel, i1, i2, i3, N_wake_edges

        ! Loop through wake-shedding edges to intialize panels
        N_wake_edges = size(wake_edge_indices)
        do k=1,N_wake_edges

            ! Get index of edge
            i = wake_edge_indices(k)

            ! Determine which wake-shedding vertices this panel lies between
            i_start = body_verts(body_edges(i)%verts(1))%index_in_wake_vertices
            i_stop = body_verts(body_edges(i)%verts(2))%index_in_wake_vertices

            ! Create panels heading downstream
            do j=1,N_panels_streamwise

                ! Determine index of first triangular panel
                i_panel = (k-1)*N_panels_streamwise*2+2*j-1

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_start-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1

                ! Initialize
                call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel)

                ! Specify this panel is in the wake
                this%panels(i_panel)%in_wake = .true.

                ! Create mirror
                if (asym_flow) then

                    ! Determine index
                    i_panel = i_panel + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_start-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2

                    ! Initialize (order of vertices is reversed to maintain panel orientation through mirror)
                    call this%panels(i_panel)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_panel)

                    ! Specify this panel is in the wake
                    this%panels(i_panel)%in_wake = .true.

                end if

                ! Determine index of second triangular panel
                i_panel = (k-1)*N_panels_streamwise*2+2*j

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j

                ! Initialize
                call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel)

                ! Specify this panel is in the wake
                this%panels(i_panel)%in_wake = .true.

                ! Create mirror
                if (asym_flow) then

                    ! Determine index
                    i_panel = i_panel + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+this%N_verts/2

                    ! Initialize (again, order is reversed)
                    call this%panels(i_panel)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_panel)

                    ! Specify this panel is in the wake
                    this%panels(i_panel)%in_wake = .true.

                end if

            end do
        end do
    
        
    end subroutine wake_mesh_init_panels


    subroutine wake_mesh_init_midpoints(this, body_edges, body_verts, wake_edge_indices)
        ! Initializes the midpoints for the wake mesh

        class(wake_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        integer,dimension(:),allocatable,intent(in) :: wake_edge_indices

        integer :: i_mid, i, j, k, m, n, m1, n1, temp, i_edge, i_midpoint_parent, i_pot_edge
        real :: distance
        integer,dimension(2) :: shared_verts
        logical :: already_found_shared
        real,dimension(3) :: loc

        i_mid = this%N_verts
        do i=1,this%N_panels

            ! Loop through each potential neighbor
            neighbor_loop: do j=i+1,this%N_panels

                ! Check if we've found all neighbors for this panel
                if (all(this%panels(i)%abutting_panels /= 0)) then
                    exit neighbor_loop
                end if

                ! Initialize for this panel pair
                already_found_shared = .false.

                ! Check if the panels are abutting
                abutting_loop: do m=1,this%panels(i)%N
                    do n=1,this%panels(j)%N

                        ! Get distance between vertices
                        distance = dist(this%panels(i)%get_vertex_loc(m), this%panels(j)%get_vertex_loc(n))

                        ! Check distance
                        if (distance < 1.e-12) then

                            ! Previously found a shared vertex, so we have abutting panels
                            if (already_found_shared) then

                                ! Store second shared vertex
                                shared_verts(2) = this%panels(i)%get_vertex_index(m)

                                ! Check order; edge should proceed counterclockwise about panel i
                                if (m1 == 1 .and. m == this%panels(i)%N) then
                                    temp = shared_verts(1)
                                    shared_verts(1) = shared_verts(2)
                                    shared_verts(2) = temp
                                end if

                                ! Initialize midpoint
                                i_mid = i_mid + 1
                                loc = 0.5*(this%vertices(shared_verts(1))%loc + this%vertices(shared_verts(2))%loc)
                                call this%vertices(i_mid)%init(loc, i_mid, 2)

                                ! Store adjacent panels

                                ! Store that i is adjacent to j
                                if ( (n1 == 1 .and. n == this%panels(j)%N) .or. (n == 1 .and. n1 == this%panels(j)%N) ) then
                                    i_edge = this%panels(j)%N
                                else
                                    n1 = min(n, n1)
                                    i_edge = n1
                                end if
                                this%panels(j)%abutting_panels(i_edge) = i
                                this%panels(j)%midpoints(i_edge)%ptr => this%vertices(i_mid)

                                ! Store that j is adjacent to i
                                if (m1 == 1 .and. m == this%panels(i)%N) then ! Nth edge
                                    i_edge = m
                                else ! 1st or 2nd edge
                                    i_edge = m1
                                end if
                                this%panels(i)%abutting_panels(i_edge) = j
                                this%panels(i)%midpoints(i_edge)%ptr => this%vertices(i_mid)

                                ! Store parent indices for midpoint
                                ! If the endpoints have the same parents, then the midpoint will too
                                if (this%vertices(shared_verts(1))%top_parent == this%vertices(shared_verts(2))%top_parent) then
                                    this%vertices(i_mid)%top_parent = this%vertices(shared_verts(1))%top_parent
                                    this%vertices(i_mid)%bot_parent = this%vertices(shared_verts(1))%bot_parent
                                
                                ! Otherwise, this midpoint's parents are also midpoints
                                else

                                    ! Loop through wake-shedding edges to find this one's parent
                                    do k=1,size(wake_edge_indices)

                                        ! Get edge index
                                        i_pot_edge = wake_edge_indices(k)

                                        ! Check the vertices match
                                        if ((this%vertices(shared_verts(1))%top_parent == body_edges(i_pot_edge)%verts(1).and.&
                                             this%vertices(shared_verts(2))%top_parent == body_edges(i_pot_edge)%verts(2)).or.&
                                            (this%vertices(shared_verts(2))%top_parent == body_edges(i_pot_edge)%verts(1).and.&
                                             this%vertices(shared_verts(1))%top_parent == body_edges(i_pot_edge)%verts(2)))then

                                            ! Get midpoint index
                                            i_midpoint_parent = body_edges(i_pot_edge)%i_midpoint

                                            ! Set parents
                                            this%vertices(i_mid)%top_parent = i_midpoint_parent
                                            this%vertices(i_mid)%bot_parent = body_verts(i_midpoint_parent)%i_wake_partner

                                        end if
                                    end do

                                end if

                                ! Break out of loop
                                exit abutting_loop

                            ! First shared vertex
                            else

                                already_found_shared = .true.
                                shared_verts(1) = this%panels(i)%get_vertex_index(m)
                                m1 = m
                                n1 = n

                            end if
                        end if

                    end do

                end do abutting_loop

            end do neighbor_loop

            ! Check for edges on empty space
            do j=1,this%panels(i)%N

                if (this%panels(i)%abutting_panels(j) == 0) then

                    ! Initialize midpoint
                    i_mid = i_mid + 1
                    loc = 0.5*(this%panels(i)%get_vertex_loc(j) + this%panels(i)%get_vertex_loc(modulo(j, this%panels(i)%N)+1))
                    call this%vertices(i_mid)%init(loc, i_mid, 2)

                    ! Point panel to it
                    this%panels(i)%midpoints(j)%ptr => this%vertices(i_mid)

                    ! Get endpoints
                    shared_verts(1) = this%panels(i)%get_vertex_index(j)
                    shared_verts(2) = this%panels(i)%get_vertex_index(modulo(j, this%panels(i)%N)+1)

                    ! Store parent indices for midpoint
                    ! If the endpoints have the same parents, then the midpoint will too
                    if (this%vertices(shared_verts(1))%top_parent == this%vertices(shared_verts(2))%top_parent) then
                        this%vertices(i_mid)%top_parent = this%vertices(shared_verts(1))%top_parent
                        this%vertices(i_mid)%bot_parent = this%vertices(shared_verts(1))%bot_parent
                    
                    ! Otherwise, this midpoint's parents are also midpoints
                    else

                        ! Loop through wake-shedding edges to find this one's parent
                        do k=1,size(wake_edge_indices)

                            ! Get edge index
                            i_pot_edge = wake_edge_indices(k)

                            ! Check the vertices match
                            if ((this%vertices(shared_verts(1))%top_parent == body_edges(i_pot_edge)%verts(1).and.&
                                 this%vertices(shared_verts(2))%top_parent == body_edges(i_pot_edge)%verts(2)).or.&
                                (this%vertices(shared_verts(2))%top_parent == body_edges(i_pot_edge)%verts(1).and.&
                                 this%vertices(shared_verts(1))%top_parent == body_edges(i_pot_edge)%verts(2)))then

                                ! Get midpoint index
                                i_midpoint_parent = body_edges(i_pot_edge)%i_midpoint

                                ! Set parents
                                this%vertices(i_mid)%top_parent = i_midpoint_parent
                                this%vertices(i_mid)%bot_parent = body_verts(i_midpoint_parent)%i_wake_partner

                            end if
                        end do

                    end if

                end if

            end do

        end do
        
    end subroutine wake_mesh_init_midpoints


end module wake_mesh_mod