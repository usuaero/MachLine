! Class for modeling wake meshes
module wake_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use linked_list_mod
    use helpers_mod
    use base_geom_mod
    use panel_mod
    use math_mod
    use flow_mod
    use vtk_mod
    use wake_strip_mod
    use mesh_mod

    implicit none


    type, extends(mesh) ::  wake_mesh

        type(wake_strip),allocatable,dimension(:) :: strips
        integer :: N_strips = 0
        integer :: N_max_strip_verts = 0
        integer :: N_max_strip_panels = 0

        contains

            procedure :: init => wake_mesh_init
            procedure :: init_strips => wake_mesh_init_strips
            procedure :: init_vertices => wake_mesh_init_vertices
            procedure :: init_panels => wake_mesh_init_panels
            procedure :: init_panel => wake_mesh_init_panel
            procedure :: init_midpoints => wake_mesh_init_midpoints
            procedure :: determine_midpoint_parents => wake_mesh_determine_midpoint_parents
            procedure :: write_wake => wake_mesh_write_wake
            procedure :: write_strips => wake_mesh_write_strips

    end type wake_mesh


contains


    subroutine wake_mesh_init(this, body_edges, body_verts, freestream, asym_flow, mirror_plane, N_panels_streamwise, &
                              trefftz_dist, body_mirrored)
        ! Initializes the wake mesh

        implicit none

        class(wake_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_panels_streamwise
        real,intent(in) :: trefftz_dist
        logical,intent(in) :: body_mirrored

        integer :: i, j
        integer :: N_wake_edge_verts, N_wake_edges
        integer,dimension(:),allocatable :: wake_edge_indices, wake_edge_verts
        logical,dimension(:),allocatable :: is_wake_edge_vertex

        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Initializing wake with a Trefftz distance of ", trefftz_dist, "..."

        ! Set whether the wake will be mirrored
        this%mirrored = body_mirrored .and. .not. asym_flow
        this%mirror_plane = mirror_plane

        ! Initialize strips
        call this%init_strips(body_edges, body_verts, freestream, asym_flow, mirror_plane, N_panels_streamwise, trefftz_dist)

        !! Count up wake-shedding edges
        !N_wake_edges = 0
        !do i=1,size(body_edges)
        !    if (body_edges(i)%sheds_wake) N_wake_edges = N_wake_edges + 1
        !end do

        !! Get indices of wake-shedding edges and vertices
        !allocate(wake_edge_indices(N_wake_edges))
        !allocate(is_wake_edge_vertex(size(body_verts)), source=.false.)
        !j = 0
        !do i=1,size(body_edges)
        !    if (body_edges(i)%sheds_wake) then

        !        ! Store edge index
        !        j = j + 1
        !        wake_edge_indices(j) = i

        !        ! Store that the vertices are wake-shedding
        !        is_wake_edge_vertex(body_edges(i)%top_verts(1)) = .true.
        !        is_wake_edge_vertex(body_edges(i)%top_verts(2)) = .true.

        !    end if
        !end do

        !! Determine how many vertices are along the wake-shedding edges
        !N_wake_edge_verts = count(is_wake_edge_vertex)

        !! Get wake-shedding edge vertex indices
        !allocate(wake_edge_verts(N_wake_edge_verts))
        !j = 0
        !do i=1,size(body_verts)
        !    if (is_wake_edge_vertex(i)) then
        !        j = j + 1
        !        wake_edge_verts(j) = i
        !        body_verts(i)%index_in_wake_vertices = j
        !    end if
        !end do
        !deallocate(is_wake_edge_vertex)

        !! Determine necessary number of vertices
        !if (asym_flow) then
        !    this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)*2
        !else
        !    this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)
        !end if

        !! Determine necessary number of panels
        !if (asym_flow) then
        !    this%N_panels = N_wake_edges*N_panels_streamwise*4
        !else
        !    this%N_panels = N_wake_edges*N_panels_streamwise*2
        !end if

        !! Initialize vertices
        !call this%init_vertices(body_verts, freestream, wake_edge_verts, asym_flow, trefftz_dist, N_panels_streamwise, mirror_plane)

        !! Initialize panels
        !call this%init_panels(body_edges, body_verts, N_panels_streamwise, wake_edge_indices, asym_flow, N_body_panels)

        !! Initialize midpoints (if needed)
        !if (doublet_order == 2) then
        !    call this%init_midpoints(body_edges, body_verts, wake_edge_indices, asym_flow)
        !end if

        !! Initialize freestream-dependent properties of panels once the midpoints have been created
        !! The mirror of wake panels will never need to be initialized
        !do i=1,this%N_panels
        !    call this%panels(i)%init_with_flow(freestream, .false., mirror_plane)
        !    call this%panels(i)%set_influencing_verts()
        !end do

        if (verbose) write(*,'(a i7 a i7 a i7 a)') "Done. Created ", this%N_verts, " wake vertices and ", &
                                              this%N_panels, " wake panels distributed between ", &
                                              this%N_strips, " strips."

    end subroutine wake_mesh_init


    subroutine wake_mesh_init_strips(this, body_edges, body_verts, freestream, asym_flow, mirror_plane, &
                              N_panels_streamwise, trefftz_dist)
        ! Creates the strips for this wake

        implicit none

        class(wake_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_panels_streamwise
        real,intent(in) :: trefftz_dist

        integer :: i, i_strip, i_start_edge
        type(list) :: wake_shedding_edges

        ! Loop through edges to find which ones shed a wake and how many there are
        this%N_strips = 0
        do i=1,size(body_edges)

            ! Check if it sheds a wake
            if (body_edges(i)%sheds_wake) then
                this%N_strips = this%N_strips + 1
                call wake_shedding_edges%append(i)

                ! If the flow is asymmetric, then wake-shedding edges not on the mirror plane will also be used to generate a wake strip
                if (asym_flow .and. .not. body_edges(i)%on_mirror_plane) then
                    this%N_strips = this%N_strips + 1
                end if
            end if

        end do

        ! Allocate strip storage
        allocate(this%strips(this%N_strips))
    
        ! Initialize strips
        i = 0
        i_strip = 0
        do while (i_strip < this%N_strips)

            ! Get starting edge index
            i = i + 1
            call wake_shedding_edges%get(i, i_start_edge)

            ! Initialize strip
            i_strip = i_strip + 1
            call this%strips(i_strip)%init(freestream, body_edges(i_start_edge), .false., mirror_plane, &
                                           N_panels_streamwise, trefftz_dist, body_verts, this%mirrored)

            ! Check if we need to create a mirror strip
            if (asym_flow .and. .not. body_edges(i_start_edge)%on_mirror_plane) then
                i_strip = i_strip + 1
                call this%strips(i_strip)%init(freestream, body_edges(i_start_edge), .true., mirror_plane, &
                                               N_panels_streamwise, trefftz_dist, body_verts, this%mirrored)
            end if
        end do

        ! Find out the maximum number of vertices and panels between the strips and totals
        this%N_verts = 0
        this%N_panels = 0
        do i=1,this%N_strips

            ! Find maxima
            this%N_max_strip_panels = max(this%N_max_strip_panels, this%strips(i)%N_panels)
            this%N_max_strip_verts = max(this%N_max_strip_verts, this%strips(i)%N_verts)

            ! Sum totals
            this%N_verts = this%N_verts + this%strips(i)%N_verts
            this%N_panels = this%N_panels + this%strips(i)%N_panels
        end do
        
    end subroutine wake_mesh_init_strips


    subroutine wake_mesh_init_vertices(this, body_verts, freestream, wake_edge_verts, asym_flow, &
                                       trefftz_dist, N_panels_streamwise, mirror_plane)
        ! Initializes the vertices for the mesh
        ! Obsolete

        class(wake_mesh),intent(inout) :: this
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(in) :: wake_edge_verts
        logical,intent(in) :: asym_flow
        real,intent(in) :: trefftz_dist
        integer,intent(in) :: N_panels_streamwise, mirror_plane

        integer :: i_vert, i, j, i_top_parent, i_bot_parent, i_mirrored_vert, N_wake_edge_verts, N_body_verts
        real :: distance, vertex_separation, mirrored_distance, mirrored_vertex_separation
        real,dimension(3) :: start, loc, mirrored_start

        ! Allocate memory
        allocate(this%vertices(this%N_verts))

        ! Determine vertex placement
        i_vert = 0
        N_body_verts = size(body_verts)
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

                ! Double-check starting location for starting points not on the mirror plane
                if (i_bot_parent < N_body_verts) then
                    if (dist(start, body_verts(i_bot_parent)%loc) > 1.e-12) then
                        write(*,*) "!!! Wake edge vertices are not identical. Quitting..."
                        stop
                    end if
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
                        this%vertices(i_mirrored_vert)%top_parent = i_top_parent + N_body_verts
                        this%vertices(i_mirrored_vert)%bot_parent = i_bot_parent + N_body_verts

                    end if

                end do
            end if
        end do
    
        
    end subroutine wake_mesh_init_vertices


    subroutine wake_mesh_init_panels(this, body_edges, body_verts, N_panels_streamwise, wake_edge_indices, asym_flow, N_body_panels)
        ! Initializes the wake panels
        ! Obsolete

        class(wake_mesh),intent(inout) :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        integer,intent(in) :: N_panels_streamwise, N_body_panels
        integer,dimension(:),allocatable,intent(in) :: wake_edge_indices
        logical,intent(in) :: asym_flow

        integer :: i, j, k, i_start, i_stop, i_panel, i1, i2, i3, N_wake_edges, N_skipped
        logical,dimension(:),allocatable :: skipped_panels
        type(panel),dimension(:),allocatable :: temp_panels

        ! Allocate memory
        allocate(this%panels(this%N_panels))

        ! Initialize storage for skipping zero-area panels
        allocate(skipped_panels(this%N_panels), source=.false.)

        ! Loop through wake-shedding edges to intialize panels
        N_wake_edges = size(wake_edge_indices)
        do k=1,N_wake_edges

            ! Get index of edge
            i = wake_edge_indices(k)

            ! Determine which wake-shedding vertices this panel lies between
            i_start = body_verts(body_edges(i)%top_verts(1))%index_in_wake_vertices
            i_stop = body_verts(body_edges(i)%top_verts(2))%index_in_wake_vertices

            ! Create panels heading downstream
            do j=1,N_panels_streamwise

                ! Determine index of first triangular panel
                i_panel = (k-1)*N_panels_streamwise*2+2*j-1

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_start-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1

                ! Create panel
                call this%init_panel(i_panel, i1, i2, i3, skipped_panels, .false.)

                ! Create mirror
                if (asym_flow) then

                    ! Determine index
                    i_panel = i_panel + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_start-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2

                    ! Create panel
                    call this%init_panel(i_panel, i1, i2, i3, skipped_panels, .true.)

                end if

                ! Determine index of second triangular panel
                i_panel = (k-1)*N_panels_streamwise*2+2*j

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j

                ! Create panel
                call this%init_panel(i_panel, i1, i2, i3, skipped_panels, .false.)

                ! Create mirror
                if (asym_flow) then

                    ! Determine index
                    i_panel = i_panel + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+this%N_verts/2

                    ! Create panel
                    call this%init_panel(i_panel, i1, i2, i3, skipped_panels, .true.)

                end if

            end do
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
        
    end subroutine wake_mesh_init_panels


    subroutine wake_mesh_init_panel(this, i_panel, i1, i2, i3, skipped_panels, mirror)
        ! Creates the specified panel
        ! Obsolete

        implicit none
        
        class(wake_mesh),intent(inout) :: this
        integer,intent(in) :: i_panel, i1, i2, i3
        logical,dimension(:),allocatable,intent(inout) :: skipped_panels
        logical,intent(in) :: mirror

        ! Check for zero area
        if (this%has_zero_area(i1, i2, i3)) then
            skipped_panels(i_panel) = .true.
        else

            ! Initialize object
            if (mirror) then
                call this%panels(i_panel)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_panel)
            else
                call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel)
            end if

            ! Specify this panel is in the wake
            this%panels(i_panel)%in_wake = .true.

        end if
        
    end subroutine wake_mesh_init_panel


    subroutine wake_mesh_init_midpoints(this, body_edges, body_verts, wake_edge_indices, asym_flow)
        ! Initializes the midpoints for the wake mesh
        ! Obsolete

        class(wake_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        integer,dimension(:),allocatable,intent(in) :: wake_edge_indices
        logical,intent(in) :: asym_flow

        integer :: i_mid, i, j, k, m, n, m1, n1, temp, i_edge, i_pot_edge, j_next, i_start, N_mids
        real :: distance
        integer,dimension(2) :: shared
        integer,dimension(:,:),allocatable :: shared_verts
        integer,dimension(:),allocatable :: edge_for_panel_i, edge_for_panel_j, panel_i, panel_j
        logical :: already_found_shared
        real,dimension(3) :: loc

        ! Initialize
        N_mids = 0
        i_mid = this%N_verts
        allocate(shared_verts(2,3*this%N_panels))
        allocate(edge_for_panel_i(3*this%N_panels))
        allocate(edge_for_panel_j(3*this%N_panels))
        allocate(panel_i(3*this%N_panels))
        allocate(panel_j(3*this%N_panels))

        ! Loop through panels
        panel_loop: do i=1,this%N_panels

            ! Loop through each potential neighbor
            neighbor_loop: do j=i+1,this%N_panels

                ! Check if we've found all neighbors for this panel
                if (all(this%panels(i)%abutting_panels /= 0)) cycle panel_loop

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

                                ! Update number of midpoints
                                N_mids = N_mids + 1

                                ! Store second shared vertex
                                shared(2) = this%panels(i)%get_vertex_index(m)

                                ! Check order; edge should proceed counterclockwise about panel i
                                if (m1 == 1 .and. m == this%panels(i)%N) then
                                    temp = shared(1)
                                    shared(1) = shared(2)
                                    shared(2) = temp
                                end if

                                ! Store in arrays
                                shared_verts(:,N_mids) = shared
                                panel_i(N_mids) = i
                                panel_j(N_mids) = j

                                ! Determine which edge this is for each panel
                                if ( (n1 == 1 .and. n == this%panels(j)%N) .or. (n == 1 .and. n1 == this%panels(j)%N) ) then
                                    i_edge = this%panels(j)%N
                                else
                                    n1 = min(n, n1)
                                    i_edge = n1
                                end if
                                edge_for_panel_j(N_mids) = i_edge
                                this%panels(j)%abutting_panels(i_edge) = i

                                if (m1 == 1 .and. m == this%panels(i)%N) then ! Nth edge
                                    i_edge = m
                                else ! 1st or 2nd edge
                                    i_edge = m1
                                end if
                                edge_for_panel_i(N_mids) = i_edge
                                this%panels(i)%abutting_panels(i_edge) = j

                                ! Break out of loop
                                exit abutting_loop

                            ! First shared vertex
                            else

                                already_found_shared = .true.
                                shared(1) = this%panels(i)%get_vertex_index(m)
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

                    j_next = modulo(j, this%panels(i)%N)+1

                    ! Update number of midpoints
                    N_mids = N_mids + 1

                    ! Store info
                    shared_verts(1,N_mids) = this%panels(i)%get_vertex_index(j)
                    shared_verts(2,N_mids) = this%panels(i)%get_vertex_index(j_next)
                    panel_i(N_mids) = i
                    panel_j(N_mids) = 0
                    edge_for_panel_i(N_mids) = j

                end if

            end do

        end do panel_loop

        ! Add more vertices
        call this%allocate_new_vertices(N_mids)

        ! Initialize midpoints
        i_start = this%N_verts - N_mids
        do i=1,N_mids

            i_mid = i_start + i

            ! Get location
            loc = 0.5*(this%vertices(shared_verts(1,i))%loc + this%vertices(shared_verts(2,i))%loc)
            call this%vertices(i_mid)%init(loc, i_mid, 2)

            ! Store parent indices for midpoint
            call this%determine_midpoint_parents(i_mid, shared_verts(:,i), body_edges, body_verts, wake_edge_indices, asym_flow)

            ! Store for panels
            this%panels(panel_i(i))%midpoints(edge_for_panel_i(i))%ptr => this%vertices(i_mid)
            if (panel_j(i) > 0) then
                this%panels(panel_j(i))%midpoints(edge_for_panel_j(i))%ptr => this%vertices(i_mid)
            end if

        end do

        this%midpoints_created = .true.

    end subroutine wake_mesh_init_midpoints


    subroutine wake_mesh_determine_midpoint_parents(this, i_mid, shared, body_edges, body_verts, wake_edge_indices, asym_flow)
        ! Determines which mesh vertices are the parents of the given midpoint
        ! Obsolete

        implicit none
        
        class(wake_mesh),intent(inout) :: this
        integer,intent(in) :: i_mid
        integer,dimension(2),intent(in) :: shared
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        integer,dimension(:),allocatable,intent(in) :: wake_edge_indices
        logical,intent(in) :: asym_flow

        integer :: k, i_pot_edge, N_body_verts

        N_body_verts = size(body_verts)

        ! If this midpoints parents are a vertex, then both neighboring vertices will have the same parents
        if (this%vertices(shared(1))%top_parent == this%vertices(shared(2))%top_parent) then

            this%vertices(i_mid)%top_parent = this%vertices(shared(1))%top_parent
            this%vertices(i_mid)%bot_parent = this%vertices(shared(1))%bot_parent
        
        ! Otherwise, this midpoint's parents are also midpoints
        else

            ! Loop through wake-shedding edges to find this one's parent
            potential_edge_loop: do k=1,size(wake_edge_indices)

                ! Get edge index
                i_pot_edge = wake_edge_indices(k)

                ! Check if the vertices match on original mesh
                if ((this%vertices(shared(1))%top_parent == body_edges(i_pot_edge)%top_verts(1).and.&
                     this%vertices(shared(2))%top_parent == body_edges(i_pot_edge)%top_verts(2)).or.&
                    (this%vertices(shared(2))%top_parent == body_edges(i_pot_edge)%top_verts(1).and.&
                     this%vertices(shared(1))%top_parent == body_edges(i_pot_edge)%top_verts(2)))then

                    ! Set parents
                    this%vertices(i_mid)%top_parent = body_edges(i_pot_edge)%top_midpoint
                    this%vertices(i_mid)%bot_parent = body_edges(i_pot_edge)%bot_midpoint

                    exit potential_edge_loop

                ! Check if the vertices match the mirror
                else if (asym_flow .and. &
                ((this%vertices(shared(1))%top_parent == body_edges(i_pot_edge)%top_verts(1)+N_body_verts.and.&
                  this%vertices(shared(2))%top_parent == body_edges(i_pot_edge)%top_verts(2)+N_body_verts).or.&
                 (this%vertices(shared(2))%top_parent == body_edges(i_pot_edge)%top_verts(1)+N_body_verts.and.&
                  this%vertices(shared(1))%top_parent == body_edges(i_pot_edge)%top_verts(2)+N_body_verts)))then

                    ! Get midpoint index

                    ! Set parents
                    this%vertices(i_mid)%top_parent = body_edges(i_pot_edge)%top_midpoint + N_body_verts
                    this%vertices(i_mid)%bot_parent = body_edges(i_pot_edge)%bot_midpoint + N_body_verts

                end if
            end do potential_edge_loop

        end if
    
    end subroutine wake_mesh_determine_midpoint_parents


    subroutine wake_mesh_write_wake(this, wake_file, exported, mu)
        ! Writes the wake out to file
        ! Obsolete

        implicit none
        
        class(wake_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: wake_file
        logical,intent(out) :: exported
        real,dimension(:),allocatable,intent(in),optional :: mu

        type(vtk_out) :: wake_vtk
        integer :: i
        real,dimension(:),allocatable :: mu_on_wake

        ! Clear old file
        call delete_file(wake_file)

        if (this%N_panels > 0) then

            ! Write out geometry
            call wake_vtk%begin(wake_file)
            call wake_vtk%write_points(this%vertices)
            call wake_vtk%write_panels(this%panels, subdivide=doublet_order==2, mirror=.false.)
            call wake_vtk%write_cell_normals(this%panels)

            if (present(mu)) then

                ! Calculate doublet strengths
                allocate(mu_on_wake(this%N_verts))
                do i=1,this%N_verts
                    mu_on_wake(i) = mu(this%vertices(i)%top_parent) - mu(this%vertices(i)%bot_parent)
                end do

                ! Write doublet strengths
                call wake_vtk%write_point_scalars(mu_on_wake, "mu")
            end if

            ! Finish up
            call wake_vtk%finish()
            exported = .true.

        else
            exported = .false.
        end if
        
    end subroutine wake_mesh_write_wake


    subroutine wake_mesh_write_strips(this, wake_file, exported, mu)
        ! Writes the wake strips out to file

        implicit none
        
        class(wake_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: wake_file
        logical,intent(out) :: exported
        real,dimension(:),allocatable,intent(in),optional :: mu

        type(vtk_out) :: wake_vtk
        integer :: i, j, k, N_verts, N_panels, shift
        real,dimension(:),allocatable :: mu_on_wake
        real,dimension(:,:),allocatable :: verts

        ! Clear old file
        call delete_file(wake_file)

        if (this%N_strips > 0) then

            ! Get total number of vertices and panels
            N_verts = 0
            N_panels = 0
            do i=1,this%N_strips
                N_verts = N_verts + this%strips(i)%N_verts
                N_panels = N_panels + this%strips(i)%N_panels
            end do

            ! Get all vertices
            allocate(verts(3,N_verts))
            i = 0
            do k=1,this%N_strips
                do j=1,this%strips(k)%N_verts
                    i = i + 1
                    verts(:,i) = this%strips(k)%vertices(j)%loc
                end do
            end do

            ! Initialize and write out vertices
            call wake_vtk%begin(wake_file)
            call wake_vtk%write_points(verts)

            ! Write out panels
            shift = 0
            do k=1,this%N_strips
                call wake_vtk%write_panels(this%strips(k)%panels, subdivide=doublet_order==2, mirror=.false., &
                                           vertex_index_shift=shift, N_total_panels=N_panels)
                shift = shift + this%strips(k)%N_verts
            end do

            if (present(mu)) then

                ! Calculate doublet strengths
                allocate(mu_on_wake(N_verts))
                i = 0
                do k=1,this%N_strips
                    do j=1,this%strips(k)%N_verts
                        i = i + 1
                        mu_on_wake(i) = mu(this%strips(k)%vertices(j)%top_parent) - mu(this%strips(k)%vertices(j)%bot_parent)
                    end do
                end do

                ! Write doublet strengths
                call wake_vtk%write_point_scalars(mu_on_wake, "mu")
            end if

            ! Finish up
            call wake_vtk%finish()
            exported = .true.

        else
            exported = .false.
        end if
        
    end subroutine wake_mesh_write_strips


end module wake_mesh_mod