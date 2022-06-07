! A surface mesh type encompassing a body, wakes, and shocks
module surface_mesh_mod

    use omp_lib
    use json_mod
    use json_xtnsn_mod
    use vtk_mod
    use stl_mod
    use tri_mod
    use vertex_mod
    use panel_mod
    use flow_mod
    use math_mod
    use edge_mod
    use wake_mesh_mod
    use sort_mod
    use helpers_mod

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels, N_cp, N_edges
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        type(edge),allocatable,dimension(:) :: edges
        type(wake_mesh) :: wake
        real :: C_wake_shedding_angle, trefftz_distance, C_min_panel_angle
        integer :: N_wake_panels_streamwise
        logical :: wake_present, append_wake
        real,dimension(:,:),allocatable :: cp, cp_mirrored
        real,dimension(:),allocatable :: phi_cp, phi_cp_sigma, phi_cp_mu ! Induced potentials at control points
        real,dimension(:),allocatable :: Phi_u ! Total potential on outer surface
        real,dimension(:),allocatable :: C_p_pg, C_p_lai, C_p_kt ! Corrected surface pressure coefficients
        real,dimension(:),allocatable :: C_p_inc, C_p_ise, C_p_2nd, C_p_sln, C_p_lin ! Surface pressure coefficients
        real,dimension(:,:),allocatable :: V, dC_f ! Surface velocities and pressure forces
        real :: control_point_offset
        logical :: mirrored ! Whether the mesh is to be mirrored about any planes
        integer :: mirror_plane ! Index of the plane across which the mesh is mirrored (1: yz, 2: xz, 3: xy); this is the index of the normal to that plane
        logical :: asym_flow ! Whether the flow is asymmetric about the mirror plane
        logical :: found_discontinuous_edges
        real,dimension(:),allocatable :: mu, sigma ! Singularity strengths
        real :: S_ref ! Reference parameters

        contains

            procedure :: init => surface_mesh_init
            procedure :: init_with_flow => surface_mesh_init_with_flow
            procedure :: output_results => surface_mesh_output_results
            procedure :: locate_adjacent_panels => surface_mesh_locate_adjacent_panels
            procedure :: store_adjacent_vertices => surface_mesh_store_adjacent_vertices
            procedure :: characterize_edges => surface_mesh_characterize_edges
            procedure :: find_vertices_on_mirror => surface_mesh_find_vertices_on_mirror
            procedure :: get_indices_to_panel_vertices => surface_mesh_get_indices_to_panel_vertices
            procedure :: allocate_new_vertices => surface_mesh_allocate_new_vertices
            procedure :: clone_vertices => surface_mesh_clone_vertices
            procedure :: set_up_mirroring => surface_mesh_set_up_mirroring
            procedure :: get_vertex_sorting_indices => surface_mesh_get_vertex_sorting_indices
            procedure :: rearrange_vertices => surface_mesh_rearrange_vertices
            procedure :: calc_vertex_normals => surface_mesh_calc_vertex_normals
            procedure :: init_wake => surface_mesh_init_wake
            procedure :: update_supersonic_trefftz_distance => surface_mesh_update_supersonic_trefftz_distance
            procedure :: update_subsonic_trefftz_distance => surface_mesh_update_subsonic_trefftz_distance
            procedure :: place_interior_control_points => surface_mesh_place_interior_control_points

    end type surface_mesh
    

contains


    subroutine surface_mesh_init(this, settings)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(json_value),pointer,intent(inout) :: settings

        character(len=:),allocatable :: extension
        integer :: loc, i
        character(len=:),allocatable :: mesh_file, mirror_plane
        logical :: file_exists
        real :: wake_shedding_angle

        ! Set singularity orders
        call json_xtnsn_get(settings, 'singularity_order.doublet', doublet_order, 1)
        call json_xtnsn_get(settings, 'singularity_order.source', source_order, 0)
        if (verbose) write(*,'(a, i1, a, i1, a)') "     User has selected ", doublet_order, &
                                                  "-order doublet and ", source_order, "-order source panels."

        ! Check source distribution
        if (source_order < 0 .or. source_order > 1) then
            write(*,*) "!!! Only constant or linear source distributions are allowed. Quitting..."
            stop
        end if

        ! Check doublet distribution
        if (doublet_order < 1 .or. doublet_order > 2) then
            write(*,*) "!!! Only linear or quadratic doublet distributions are allowed. Quitting..."
            stop
        end if

        ! Get mesh file
        call json_get(settings, 'file', mesh_file)
        mesh_file = trim(mesh_file)

        ! Check mesh file exists
        inquire(file=mesh_file, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) "!!! Mesh file", mesh_file, "does not exist. Quitting..."
            stop
        end if

        ! Determine the type of mesh file
        loc = index(mesh_file, '.')
        extension = mesh_file(loc:len(mesh_file))

        ! Load mesh file
        if (verbose) write(*,'(a, a, a)',advance='no') "     Reading surface mesh in from file ", mesh_file, "..."
        if (extension == '.vtk') then
            call load_surface_vtk(mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        else if (extension == '.stl') then
            call load_surface_stl(mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        else if (extension == '.tri') then
            call load_surface_tri(mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        else
            write(*,*) "MachLine cannot read ", extension, " type mesh files. Quitting..."
            stop
        end if
        if (verbose) write(*,*) "Done."

        ! Display mesh info
        if (verbose) write(*,'(a, i7, a, i7, a)') "     Surface mesh has ", this%N_verts, " vertices and ", &
                                                  this%N_panels, " panels."

        ! Get mirroring
        call json_xtnsn_get(settings, 'mirror_about', mirror_plane, "none")
        this%mirrored = .true.
        select case (mirror_plane)
        case ("xy")
            this%mirror_plane = 3
            if (verbose) write(*,*) "    Mesh set to mirror about xy plane."
        case ("xz")
            this%mirror_plane = 2
            if (verbose) write(*,*) "    Mesh set to mirror about xz plane."
        case ("yz")
            this%mirror_plane = 1
            if (verbose) write(*,*) "    Mesh set to mirror about yz plane."
        case default
            this%mirror_plane = 0
            this%mirrored = .false.
        end select

        ! Check if the user wants a wake
        call json_xtnsn_get(settings, 'wake_model.wake_present', this%wake_present, .true.)
        call json_xtnsn_get(settings, 'wake_model.append_wake', this%append_wake, this%wake_present)

        ! Check that we're not trying to append a wake which is not present
        if (.not. this%wake_present .and. this%append_wake) then
            this%append_wake = .false.
        end if

        ! Store settings for wake models
        if (this%wake_present) then
            call json_xtnsn_get(settings, 'wake_model.wake_shedding_angle', wake_shedding_angle, 90.0) ! Maximum allowable angle between panel normals without having separation
            this%C_wake_shedding_angle = cos(wake_shedding_angle*pi/180.0)

            if (this%append_wake) then
                call json_xtnsn_get(settings, 'wake_model.trefftz_distance', this%trefftz_distance, -1.0) ! Distance from origin to wake termination
                call json_xtnsn_get(settings, 'wake_model.N_panels', this%N_wake_panels_streamwise, 1)
            end if
        end if

        ! Store references
        call json_xtnsn_get(settings, 'reference.area', this%S_ref, 1.0)

        ! Locate which vertices are on the mirror plane
        if (this%mirrored) then
            call this%find_vertices_on_mirror()
        end if

        ! Locate adjacent panels
        call this%locate_adjacent_panels()

    end subroutine surface_mesh_init


    subroutine surface_mesh_find_vertices_on_mirror(this)
        ! Locates which vertices are on the mirror plane

        implicit none

        class(surface_mesh),intent(inout) :: this

        integer :: i

        ! Search for vertices lying on mirror plane
        do i=1,this%N_verts

            ! Check coordinate normal to mirror plane
            if (abs(this%vertices(i)%loc(this%mirror_plane))<1e-12) then

                ! The vertex is on the mirror plane
                this%vertices(i)%on_mirror_plane = .true.

            end if
            
        end do
    
    end subroutine surface_mesh_find_vertices_on_mirror


    subroutine surface_mesh_locate_adjacent_panels(this)
        ! Loops through panels to determine which are adjacent

        implicit none

        class(surface_mesh),intent(inout),target :: this

        integer :: i, j, m, n, m1, n1, temp, i_edge, N_edges, i_mid, N_orig_verts
        logical :: already_found_shared, watertight
        real :: distance
        integer,dimension(2) :: i_endpoints
        integer,dimension(this%N_panels*3) :: panel1, panel2, vertex1, vertex2, edge_index1, edge_index2 ! By definition, we will have no more than 3*N_panels edges; we should have much less, but some people...
        real,dimension(this%N_panels*3) :: edge_length
        logical,dimension(this%N_panels*3) :: on_mirror_plane
        real,dimension(3) :: loc
        

        if (verbose) write(*,'(a)',advance='no') "     Locating adjacent panels..."

        ! Initialize a few things
        on_mirror_plane = .false.

        ! Loop through each panel
        this%N_edges = 0
        !$OMP parallel private(j, already_found_shared, distance, i_endpoints, m, m1, n, n1, temp, i_edge)

        !$OMP do schedule(dynamic)
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

                        ! Check if vertices have the same index
                        if (this%panels(i)%get_vertex_index(m) == this%panels(j)%get_vertex_index(n)) then

                            ! Previously found a shared vertex, so we have abutting panels
                            if (already_found_shared) then

                                ! Store second shared vertex
                                i_endpoints(2) = this%panels(i)%get_vertex_index(m)

                                ! Check order; edge should proceed counterclockwise about panel i
                                if (m1 == 1 .and. m == this%panels(i)%N) then
                                    temp = i_endpoints(1)
                                    i_endpoints(1) = i_endpoints(2)
                                    i_endpoints(2) = temp
                                end if

                                !$OMP critical
                                
                                ! Update number of edges
                                this%N_edges = this%N_edges + 1
                                i_edge = this%N_edges

                                ! Store vertices being adjacent to one another
                                call this%store_adjacent_vertices(i_endpoints, i_edge)

                                ! Store adjacent panels and panel edges
                                ! This stores the adjacent panels and edges according to the index of that edge
                                ! for the current panel

                                ! Store that i is adjacent to j
                                ! This one is more complicated because we don't know that n1 will be less than n; just the nature of the nested loop.
                                ! Basically, if one is 1 and the other is N, then we're dealing with edge N for panel j.
                                ! Otherwise, we're dealing with abs(n1-n) being 1, meaning edge min(n1, n).
                                if ( (n1 == 1 .and. n == this%panels(j)%N) .or. (n == 1 .and. n1 == this%panels(j)%N) ) then
                                    this%panels(j)%abutting_panels(this%panels(j)%N) = i
                                    edge_index2(i_edge) = this%panels(j)%N
                                else
                                    n1 = min(n, n1)
                                    this%panels(j)%abutting_panels(n1) = i
                                    edge_index2(i_edge) = n1
                                end if

                                ! Store that j is adjacent to i
                                if (m1 == 1 .and. m == this%panels(i)%N) then ! Nth edge
                                    this%panels(i)%abutting_panels(m) = j
                                    edge_index1(i_edge) = m
                                else ! 1st or 2nd edge
                                    this%panels(i)%abutting_panels(m1) = j
                                    edge_index1(i_edge) = m1
                                end if

                                ! Store information in arrays for later storage in edge objects
                                panel1(i_edge) = i
                                panel2(i_edge) = j
                                vertex1(i_edge) = i_endpoints(1)
                                vertex2(i_edge) = i_endpoints(2)
                                edge_length(i_edge) = dist(this%vertices(i_endpoints(1))%loc, this%vertices(i_endpoints(2))%loc)

                                !$OMP end critical
                                
                                ! Break out of loop
                                exit abutting_loop

                            ! First shared vertex
                            else

                                already_found_shared = .true.
                                i_endpoints(1) = this%panels(i)%get_vertex_index(m)
                                m1 = m
                                n1 = n

                            end if
                        end if

                    end do

                end do abutting_loop

            end do neighbor_loop

            ! For each panel, check if it abuts the mirror plane
            if (this%mirrored) then

                ! Initialize checks
                already_found_shared = .false.

                ! Loop through vertices
                mirror_loop: do m=1,this%panels(i)%N

                    ! Check if we've found all neighbors for this panel
                    if (all(this%panels(i)%abutting_panels /= 0)) then
                        exit mirror_loop
                    end if

                    ! Check if vertex is on the mirror plane
                    n = this%panels(i)%get_vertex_index(m)
                    if (this%vertices(n)%on_mirror_plane) then

                        ! Previously found a vertex on mirror plane, so the panels are abutting
                        if (already_found_shared) then

                            ! Store the second shared vertex
                            i_endpoints(2) = n

                            ! Check order
                            if (m1 == 1 .and. m == 3) then
                                temp = i_endpoints(1)
                                i_endpoints(1) = i_endpoints(2)
                                i_endpoints(2) = temp
                            end if

                            !$OMP critical
                            
                            ! Update number of edges
                            this%N_edges = this%N_edges + 1
                            i_edge = this%N_edges

                            ! Store adjacent vertices
                            call this%store_adjacent_vertices(i_endpoints, i_edge)

                            ! Store adjacent panel
                            if (m-m1 == 1) then
                                this%panels(i)%abutting_panels(m1) = i+this%N_panels
                                edge_index1(i_edge) = m1
                            else
                                this%panels(i)%abutting_panels(m) = i+this%N_panels
                                edge_index1(i_edge) = m
                            end if

                            ! Store in arrays for later storage in edge objects
                            panel1(i_edge) = i
                            panel2(i_edge) = i+this%N_panels
                            vertex1(i_edge) = i_endpoints(1)
                            vertex2(i_edge) = i_endpoints(2)
                            on_mirror_plane(i_edge) = .true.
                            edge_index2(i_edge) = 0 ! Just a placeholder since the second panel doesn't technically exist
                            edge_length(i_edge) = dist(this%vertices(i_endpoints(1))%loc, this%vertices(i_endpoints(2))%loc)

                            !$OMP end critical

                            ! Break out of loop
                            exit mirror_loop

                        ! First vertex on the mirror plane
                        else

                            already_found_shared = .true.
                            i_endpoints(1) = n
                            m1 = m

                        end if
                    end if

                end do mirror_loop

            end if

        end do

        ! Check for panels abutting empty space and add those edges.
        !$OMP single
        watertight = .true.
        do i=1,this%N_panels

            ! Check for an edge with no abutting panel
            do j=1,this%panels(i)%N
                if (this%panels(i)%abutting_panels(j) == 0) then

                    ! Mark that the mesh is not watertight
                    watertight = .false.

                    ! Get endpoint indices
                    i_endpoints(1) = this%panels(i)%get_vertex_index(j)
                    i_endpoints(2) = this%panels(i)%get_vertex_index(mod(j, this%panels(i)%N) + 1)

                    ! Set up an edge
                    this%N_edges = this%N_edges + 1
                    i_edge = this%N_edges
                    panel1(i_edge) = i
                    panel2(i_edge) = 0 ! Placeholder
                    vertex1(i_edge) = i_endpoints(1)
                    vertex2(i_edge) = i_endpoints(2)
                    on_mirror_plane(i_edge) = .false.
                    edge_index1(i_edge) = j
                    edge_index2(i_edge) = 0

                    ! Store adjacent vertices
                    call this%store_adjacent_vertices(i_endpoints, i_edge)

                end if
            end do
        end do
        if (run_checks .and. .not. watertight) write(*,*) "!!! The supplied mesh is not watertight. &
                                                               Solution quality may be adversely affected."

        ! Allocate edge storage
        allocate(this%edges(this%N_edges))
        !$OMP end single

        ! Initialize edges
        !$OMP do schedule(static)
        do i=1,this%N_edges

            ! Initialize
            call this%edges(i)%init(vertex1(i), vertex2(i), panel1(i), panel2(i), edge_length(i))

            ! Store more information
            this%edges(i)%on_mirror_plane = on_mirror_plane(i)
            this%edges(i)%edge_index_for_panel(1) = edge_index1(i)
            this%edges(i)%edge_index_for_panel(2) = edge_index2(i)

        end do

        !$OMP end parallel

        ! Initialize midpoint vertex objects
        if (doublet_order == 2) then

            ! Allocate more space
            call this%allocate_new_vertices(this%N_edges)
            N_orig_verts = this%N_verts - this%N_edges

            do i=1,this%N_edges

                ! Determine location
                loc = 0.5*(this%vertices(vertex1(i))%loc + this%vertices(vertex2(i))%loc)

                ! Initialize vertex object
                i_mid = N_orig_verts + i
                call this%vertices(i_mid)%init(loc, i_mid, 2)
                this%vertices(i_mid)%l_avg = edge_length(i)

                ! Add adjacent vertices
                call this%vertices(i_mid)%adjacent_vertices%append(vertex1(i))
                call this%vertices(i_mid)%adjacent_vertices%append(vertex2(i))

                ! Add adjacent panels
                call this%vertices(i_mid)%panels%append(panel1(i))
                if (panel2(i) > 0 .and. panel2(i) <= this%N_panels) then
                    call this%vertices(i_mid)%panels%append(panel2(i))
                end if

                ! Add edge
                call this%vertices(i_mid)%adjacent_edges%append(i)

                ! Point edge to it
                this%edges(i)%i_midpoint = i_mid

                ! Point panels to it
                this%panels(panel1(i))%midpoints(edge_index1(i))%ptr => this%vertices(i_mid)

                if (panel2(i) > 0 .and. panel2(i) <= this%N_panels) then
                    this%panels(panel2(i))%midpoints(edge_index2(i))%ptr => this%vertices(i_mid)
                end if

            end do
        end if

        if (verbose) write(*,"(a, i7, a)") "Done. Found ", this%N_edges, " edges."

        ! Print info for added midpoint vertices
        if (verbose .and. doublet_order == 2) write(*,"(a, i7, a, i7, a)") "     For a quadratic doublet distribution, ", &
                                                                           this%N_edges, " vertices were added. Mesh now has ", &
                                                                           this%N_verts, " vertices."
    
    end subroutine surface_mesh_locate_adjacent_panels


    subroutine surface_mesh_store_adjacent_vertices(this, i_verts, i_edge)
        ! Stores that the given two vertices are adjacent (i.e. they share an edge)

        implicit none

        class(surface_mesh), intent(inout) :: this
        integer,dimension(2), intent(in) :: i_verts
        integer,intent(in) :: i_edge

        ! Store that the vertices are adjacent
        if (.not. this%vertices(i_verts(1))%adjacent_vertices%is_in(i_verts(2))) then
            call this%vertices(i_verts(1))%adjacent_vertices%append(i_verts(2))
        end if
        if (.not. this%vertices(i_verts(2))%adjacent_vertices%is_in(i_verts(1))) then
            call this%vertices(i_verts(2))%adjacent_vertices%append(i_verts(1))
        end if

        ! Store that they touch this edge
        call this%vertices(i_verts(1))%adjacent_edges%append(i_edge)
        call this%vertices(i_verts(2))%adjacent_edges%append(i_edge)
        
    end subroutine surface_mesh_store_adjacent_vertices


    subroutine surface_mesh_init_with_flow(this, freestream, wake_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream
        character(len=:),allocatable,intent(in) :: wake_file

        integer :: i

        ! Check flow symmetry condition
        this%asym_flow = .false.
        if (this%mirrored) then
            if (.not. freestream%sym_about(this%mirror_plane)) then
                this%asym_flow = .true.
            end if
        end if

        ! Calculate panel coordinate transformations
        !$OMP parallel do schedule(static)
        do i=1,this%N_panels
            call this%panels(i)%init_with_flow(freestream, this%asym_flow, this%mirror_plane)
        end do

        ! Figure out wake-shedding edges, discontinuous edges, etc.
        ! Edge-characterization is only necessary for flows with wakes
        ! According to Davis, sharp, subsonic, leading edges in supersonic flow must have discontinuous doublet strength.
        ! I don't know why this would be, except in the case of leading-edge vortex separation. But Davis doesn't
        ! model leading-edge vortices. Wake-shedding trailing edges are still discontinuous in supersonic flow. Supersonic
        ! leading edges should have continuous doublet strength.
        if (this%wake_present) then
            call this%characterize_edges(freestream)
        end if

        ! Set up mirroring
        if (this%mirrored) then
            call this%set_up_mirroring()
        end if

        ! Clone necessary vertices
        if (this%wake_present .or. freestream%supersonic) then
            call this%clone_vertices()
        end if

        ! For supersonic flows, rearrange the vertices to proceed in the freestream direction
        if (freestream%supersonic) then
            call this%rearrange_vertices(freestream)
        end if

        ! Calculate normals (for placing control points)
        call this%calc_vertex_normals()

        ! Initialize wake
        call this%init_wake(freestream, wake_file)

        ! Set up influencing vertex arrays
        !$OMP parallel do schedule(static)
        do i=1,this%N_panels
            call this%panels(i)%set_influencing_verts()
        end do
    
    end subroutine surface_mesh_init_with_flow


    subroutine surface_mesh_characterize_edges(this, freestream)
        ! Locates wake-shedding edges and supersonic/subsonic edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        integer :: i, j, k, m, n, temp, top_panel, bottom_panel, i_vert_1, i_vert_2, mid, N_wake_edges
        real :: C_angle, C_min_angle
        real,dimension(3) :: second_normal

        if (verbose) write(*,'(a)',advance='no') "     Characterizing edges..."

        ! Initialize
        N_wake_edges = 0
        this%found_discontinuous_edges = .false.

        !$OMP parallel private(i, j, second_normal, C_angle, i_vert_1, i_vert_2, mid) &
        !$OMP & default(none) shared(this, freestream, doublet_order, N_wake_edges) reduction(min : C_min_angle)

        ! Loop through each edge
        !$OMP do schedule(dynamic)
        do k=1,this%N_edges

            ! Get info
            i = this%edges(k)%panels(1)
            j = this%edges(k)%panels(2)

            ! Skip edge on empty space
            if (j == 0) cycle

            ! Get normal for panel j (dependent on mirroring)
            if (this%edges(k)%on_mirror_plane) then
                second_normal = mirror_across_plane(this%panels(i)%n_g, this%mirror_plane)
            else
                second_normal = this%panels(j)%n_g
            end if

            ! Calculate angle between panels (this is the flow-turning angle; it is the most straightforward to compute)
            C_angle = inner(this%panels(i)%n_g, second_normal)

            ! Update minimum angle
            C_min_angle = min(C_angle, C_min_angle)

            ! Determine if this edge is wake-shedding

            ! Check the angle between the panels
            if (C_angle < this%C_wake_shedding_angle) then

                ! Check angle of panel normal with freestream
                if (inner(this%panels(i)%n_g, freestream%V_inf) > 0.0 .or. inner(second_normal, freestream%V_inf) > 0.0) then

                    ! Having passed the previous two checks, we've found a wake-shedding edge

                    ! Set that we've found discontinuous edges
                    this%found_discontinuous_edges = .true.

                    ! Get vertex indices (simplifies later code)
                    i_vert_1 = this%edges(k)%verts(1)
                    i_vert_2 = this%edges(k)%verts(2)

                    ! Set the character of the edge
                    this%edges(k)%sheds_wake = .true.
                    this%edges(k)%discontinuous = .true.

                    ! Update information for midpoint vertex (unique for the edge, so this doesn't need to be inside the critical block)
                    if (doublet_order == 2) then
                        mid = this%edges(k)%i_midpoint
                        this%vertices(mid)%N_wake_edges = 1
                        this%vertices(mid)%N_discont_edges = 1
                        this%vertices(mid)%clone = .true.
                    end if

                    !$OMP critical

                    ! Update number of wake-shedding edges
                    N_wake_edges = N_wake_edges + 1

                    ! Update number of wake edges touching the first vertex
                    this%vertices(i_vert_1)%N_wake_edges = this%vertices(i_vert_1)%N_wake_edges + 1
                    this%vertices(i_vert_1)%N_discont_edges = this%vertices(i_vert_1)%N_discont_edges + 1

                    ! Update number of wake edges touching the second vertex
                    this%vertices(i_vert_2)%N_wake_edges = this%vertices(i_vert_2)%N_wake_edges + 1
                    this%vertices(i_vert_2)%N_discont_edges = this%vertices(i_vert_2)%N_discont_edges + 1

                    !$OMP end critical

                end if
            end if

        end do

        ! If a given vertex is touching at least two wake-shedding edges, it will need to be cloned
        !$OMP do schedule(static)
        do i=1,this%N_verts
            if (this%vertices(i)%N_wake_edges >= 2) this%vertices(i)%clone = .true.
        end do

        !$OMP end parallel

        ! Store minimum angle
        this%C_min_panel_angle = C_min_angle

        if (verbose) write(*,'(a, i3, a, i3, a)') "Done. Found ", N_wake_edges, " wake-shedding edges."

    end subroutine surface_mesh_characterize_edges


    subroutine surface_mesh_set_up_mirroring(this)
        ! Sets up information necessary for mirroring the mesh

        implicit none

        class(surface_mesh),intent(inout) :: this

        integer :: i, j, m, n, mid

        if (verbose) write(*,'(a)',advance='no') "     Setting up mesh mirroring..."

        ! If a vertex on the mirror plane doesn't belong to a discontinuous edge, then its mirror will not be unique
        do i=1,this%N_verts

            if (this%vertices(i)%on_mirror_plane .and. this%vertices(i)%N_discont_edges == 0) then
                this%vertices(i)%mirrored_is_unique = .false.
            end if
            
        end do

        ! Check if any discontinuous edges touch the mirror plane
        do i=1,this%N_edges

            ! Check if it is discontinuous
            if (this%edges(i)%discontinuous) then

                ! Get vertex indices
                m = this%edges(i)%verts(1)
                n = this%edges(i)%verts(2)
                if (doublet_order == 2) mid = this%edges(i)%i_midpoint

                ! If a given discontinuous edge has only one of its endpoints lying on the mirror plane, then that endpoint has another
                ! adjacent edge, since the edge will be mirrored across that plane. This endpoint will need a clone, but it's mirrored
                ! vertex will be the same

                ! Check for edge touching mirror plane at only one end
                if (this%vertices(m)%on_mirror_plane .neqv. this%vertices(n)%on_mirror_plane) then

                    ! Only add discontinuous edge if the endpoint is on the mirror plane
                    if (this%vertices(m)%on_mirror_plane) then

                        this%vertices(m)%N_wake_edges = this%vertices(m)%N_wake_edges + 1
                        this%vertices(m)%N_discont_edges = this%vertices(m)%N_discont_edges + 1
                        this%vertices(m)%clone = .true.
                        this%vertices(m)%mirrored_is_unique = .false.

                    else

                        this%vertices(n)%N_wake_edges = this%vertices(n)%N_wake_edges + 1
                        this%vertices(n)%N_discont_edges = this%vertices(n)%N_discont_edges + 1
                        this%vertices(n)%clone = .true.
                        this%vertices(n)%mirrored_is_unique = .false.

                    end if

                    ! The midpoint will need to be cloned and its mirror will be unique
                    if (doublet_order == 2) then
                        this%vertices(mid)%clone = .true.
                        this%vertices(mid)%mirrored_is_unique = .true.
                    end if

                ! If the discontinuous edge has both endpoints lying on the mirror plane, then these vertices need no clones, but the mirrored
                ! vertices will still be unique
                else if (this%edges(i)%on_mirror_plane) then

                    ! The mirrored vertices will function as clones
                    this%vertices(m)%clone = .false.
                    this%vertices(n)%clone = .false.
                    this%vertices(m)%mirrored_is_unique = .true.
                    this%vertices(n)%mirrored_is_unique = .true.

                    ! Same with the midpoint
                    if (doublet_order == 2) then
                        this%vertices(mid)%clone = .false.
                        this%vertices(mid)%mirrored_is_unique = .true.
                    end if

                end if
            end if
        end do
        if (verbose) write(*,*) "Done."

    end subroutine surface_mesh_set_up_mirroring


    subroutine surface_mesh_get_indices_to_panel_vertices(this, i_vertices)
        ! Returns of the list of indices which point to the vertices (and midpoints) of each panel

        implicit none

        class(surface_mesh),intent(in) :: this
        integer,dimension(:,:),allocatable,intent(out) :: i_vertices

        integer :: i, j

        ! Allocate space
        allocate(i_vertices(8,this%N_panels))

        ! Get vertex indices for each panel since we will lose this information as soon as this%vertices is reallocated
        do i=1,this%N_panels
            do j=1,this%panels(i)%N

                ! Get vertex indices
                i_vertices(j,i) = this%panels(i)%get_vertex_index(j)
                
                ! Get midpoint indices
                if (doublet_order == 2) then
                    i_vertices(j+this%panels(i)%N,i) = this%panels(i)%get_midpoint_index(j)
                end if

            end do
        end do
        
    end subroutine surface_mesh_get_indices_to_panel_vertices


    subroutine surface_mesh_allocate_new_vertices(this, N_new_verts, i_rearrange)
        ! Adds the specified number of vertex objects to the end of the surface mesh's vertex array.
        ! Handles moving panel pointers to the new allocation of previously-existing vertices.
        ! i_rearrange may be used to rearrange the vertices for a desired behavior.

        implicit none

        class(surface_mesh),intent(inout),target :: this
        integer,intent(in) :: N_new_verts
        integer,dimension(:),allocatable,intent(in),optional :: i_rearrange
        
        type(vertex),dimension(:),allocatable :: temp_vertices
        integer :: i, j
        integer,dimension(:,:),allocatable :: i_vertices
        integer,dimension(:),allocatable :: i_rearrange_inv

        ! Get panel vertex indices
        call this%get_indices_to_panel_vertices(i_vertices)

        ! Allocate more space
        allocate(temp_vertices(this%N_verts + N_new_verts))
        allocate(i_rearrange_inv(this%N_verts))

        ! Copy over vertices
        if (present(i_rearrange)) then

            do i=1,this%N_verts

                ! Copy vertices
                temp_vertices(i) = this%vertices(i_rearrange(i))
                temp_vertices(i)%index = i

                ! Get inverse mapping
                i_rearrange_inv(i_rearrange(i)) = i
            end do

        else

            ! Copy vertices
            temp_vertices(1:this%N_verts) = this%vertices

            ! Get inverse mapping (identity)
            do i=1,this%N_verts
                i_rearrange_inv(i) = i
            end do

        end if

        ! Move allocation
        call move_alloc(temp_vertices, this%vertices)

        ! Fix vertex pointers in panel objects (necessary because this%vertices got reallocated)
        do i=1,this%N_panels
            do j=1,this%panels(i)%N

                ! Fix vertex pointers
                this%panels(i)%vertices(j)%ptr => this%vertices(i_rearrange_inv(i_vertices(j,i)))

                ! Fix midpoint pointers
                if (doublet_order == 2) then
                    this%panels(i)%midpoints(j)%ptr => this%vertices(i_rearrange_inv(i_vertices(j+this%panels(i)%N,i)))
                end if

            end do
        end do

        deallocate(i_vertices)

        ! Fix wake dependencies
        if (present(i_rearrange) .and. this%wake%N_panels > 0) then
            do i=1,this%wake%N_verts

                ! Fix top parent
                if (this%wake%vertices(i)%top_parent > this%N_verts) then
                    this%wake%vertices(i)%top_parent = i_rearrange_inv(this%wake%vertices(i)%top_parent - this%N_verts) &
                                                       + this%N_verts + N_new_verts
                else
                    this%wake%vertices(i)%top_parent = i_rearrange_inv(this%wake%vertices(i)%top_parent)
                end if

                ! Fix bottom parent
                if (this%wake%vertices(i)%bot_parent > this%N_verts) then
                    this%wake%vertices(i)%bot_parent = i_rearrange_inv(this%wake%vertices(i)%bot_parent - this%N_verts) &
                                                       + this%N_verts + N_new_verts
                else
                    this%wake%vertices(i)%bot_parent = i_rearrange_inv(this%wake%vertices(i)%bot_parent)
                end if
            end do
        end if

        ! Fix edge pointers

        ! Update number of vertices
        this%N_verts = this%N_verts + N_new_verts
        
    end subroutine surface_mesh_allocate_new_vertices


    subroutine surface_mesh_clone_vertices(this)
        ! Takes vertices which lie within discontinuous edges and splits them into two vertices.
        ! Handles rearranging of necessary dependencies.

        implicit none

        class(surface_mesh),intent(inout),target :: this

        integer :: i, j, k, m, n, N_clones, i_jango, i_boba, N_discont_verts, i_bot_panel, i_abutting_panel, i_adj_vert, i_top_panel
        integer :: i_edge
        type(vertex),dimension(:),allocatable :: temp_vertices
        integer,dimension(:),allocatable :: i_clones

        ! Check whether any discontinuities exist
        if (this%found_discontinuous_edges) then

            if (verbose) write(*,'(a)',advance='no') "     Cloning vertices at discontinuous edges..."

            ! Determine number of vertices which need to be cloned
            N_clones = 0
            do i=1,this%N_verts
                if (this%vertices(i)%clone) N_clones = N_clones + 1
            end do

            ! Add space for new vertices
            call this%allocate_new_vertices(N_clones)
            allocate(i_clones(N_clones))

            ! Initialize clones
            j = 1
            do i_jango=1,this%N_verts

                ! Check if this vertex needs to be cloned
                if (this%vertices(i_jango)%clone) then

                    ! Get index for the clone
                    i_boba = this%N_verts - N_clones + j ! Will be at position N_verts-N_clones+j in the new vertex array

                    ! Initialize clone
                    call this%vertices(i_boba)%init(this%vertices(i_jango)%loc, i_boba, 1)
                    i_clones(j) = i_boba

                    ! Specify wake partners
                    this%vertices(i_jango)%i_wake_partner = i_boba
                    this%vertices(i_boba)%i_wake_partner = i_jango

                    ! Store number of adjacent wake-shedding and discontinuous edges (probably unecessary at this point, but let's be consistent)
                    this%vertices(i_boba)%N_wake_edges = this%vertices(i_jango)%N_wake_edges
                    this%vertices(i_boba)%N_discont_edges = this%vertices(i_jango)%N_discont_edges

                    ! Copy over mirroring properties
                    this%vertices(i_boba)%mirrored_is_unique = this%vertices(i_jango)%mirrored_is_unique
                    this%vertices(i_boba)%on_mirror_plane = this%vertices(i_jango)%on_mirror_plane

                    ! Copy over adjacent panels
                    do k=1,this%vertices(i_jango)%panels%len()

                        ! Get adjacent panel index from original vertex
                        call this%vertices(i_jango)%panels%get(k, i_abutting_panel)

                        ! Copy to clone
                        call this%vertices(i_boba)%panels%append(i_abutting_panel)

                        ! Copy to original vertex's panels_not_across_wake_edge list (bottom panels will be removed)
                        call this%vertices(i_jango)%panels_not_across_wake_edge%append(i_abutting_panel)

                    end do

                    ! Copy over adjacent vertices
                    do k=1,this%vertices(i_jango)%adjacent_vertices%len()

                        ! Get adjacent panel index from original vertex
                        call this%vertices(i_jango)%adjacent_vertices%get(k, i_adj_vert)

                        ! Copy to new vertex
                        call this%vertices(i_boba)%adjacent_vertices%append(i_adj_vert)

                    end do

                    ! Remove bottom panels from top vertex and give them to the bottom vertex
                    ! Loop through edges adjacent to this vertex
                    do n=1,this%vertices(i_jango)%adjacent_edges%len()

                        ! Get edge index
                        call this%vertices(i_jango)%adjacent_edges%get(n, i_edge)

                        ! Copy to new vertex
                        call this%vertices(i_boba)%adjacent_edges%append(i_edge)

                        ! Check if this is a wake-shedding edge
                        if (this%edges(i_edge)%sheds_wake) then

                            ! Get bottom panel index
                            i_top_panel = this%edges(i_edge)%panels(1)
                            i_bot_panel = this%edges(i_edge)%panels(2)

                            ! Make sure this bottom panel is not a mirrored panel
                            if (i_bot_panel <= this%N_panels) then

                                ! Remove bottom panel index from original vertex
                                call this%vertices(i_jango)%panels_not_across_wake_edge%delete(i_bot_panel)

                                ! Add to clone
                                if (.not. this%vertices(i_boba)%panels_not_across_wake_edge%is_in(i_bot_panel)) then
                                    call this%vertices(i_boba)%panels_not_across_wake_edge%append(i_bot_panel)
                                end if

                                ! If there are any panels attached to this vertex and abutting the bottom panel, shift them over as well
                                do m=1,this%panels(i_bot_panel)%N

                                    ! Get the index of the panel abutting this bottom panel
                                    i_abutting_panel = this%panels(i_bot_panel)%abutting_panels(m)

                                    ! Check if it is not the top panel
                                    if (i_abutting_panel /= i_top_panel) then

                                        ! Make sure the abutting panel is not a mirrored panel
                                        if (i_abutting_panel <= this%N_panels) then

                                            ! See if this panel touches the vertex
                                            if (this%panels(i_abutting_panel)%touches_vertex(i_jango)) then

                                                ! Remove from original vertex
                                                call this%vertices(i_jango)%panels_not_across_wake_edge%delete(i_abutting_panel)

                                                ! Add to cloned vertex
                                                if (.not.this%vertices(i_boba)%panels_not_across_wake_edge%is_in(i_abutting_panel))&
                                                    then
                                                    call this%vertices(i_boba)%panels_not_across_wake_edge%append(i_abutting_panel)
                                                end if

                                            end if
                                        end if
                                    end if
                                end do

                            end if

                        end if

                    end do

                    ! Update bottom panels to point to cloned vertex
                    do k=1,this%vertices(i_boba)%panels_not_across_wake_edge%len()

                        ! Get panel index
                        call this%vertices(i_boba)%panels_not_across_wake_edge%get(k, i_bot_panel)

                        ! Update (doesn't need to be done for mirrored panels)
                        if (i_bot_panel <= this%N_panels) then
                            call this%panels(i_bot_panel)%point_to_new_vertex(this%vertices(i_boba))
                        end if

                    end do

                    ! Update clone index
                    j = j + 1

                else

                    ! If this vertex did not need to be cloned, but it is on the mirror plane and its mirror is unique
                    ! then the wake strength will be determined by its mirror as well in the case of an asymmetric flow.
                    if (this%mirrored .and. this%asym_flow .and.  this%vertices(i_jango)%on_mirror_plane .and. &
                        this%vertices(i_jango)%mirrored_is_unique) then

                        this%vertices(i_jango)%i_wake_partner = i_jango + this%N_verts

                    end if

                end if

            end do

            ! Store that the cloned vertices are clones
            do i=1,N_clones
                this%vertices(i_clones(i))%clone = .true.
            end do

            if (verbose) write(*,'(a, i4, a, i7, a)') "Done. Cloned ", N_clones, " vertices. Mesh now has ", &
                                                      this%N_verts, " vertices."

        end if

    end subroutine surface_mesh_clone_vertices


    subroutine surface_mesh_get_vertex_sorting_indices(this, freestream, i_sorted)
        ! Determines the vertex indices which will sort them in the conpressibility direction

        implicit none

        class(surface_mesh),intent(in) :: this
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(out) :: i_sorted

        real,dimension(:),allocatable :: x
        integer :: i, j, k

        ! Allocate the compressibility distance array
        allocate(x(this%N_verts))

        ! Add compressibility distance of each vertex
        ! What this really should is the compressibility distance of the most-upstream vertex
        ! of the panels touching this vertex, since influences may be transmitted upstream
        ! via a panel.
        do i=1,this%N_verts

            ! Initialize with the compressibility distance of this vertex
            x(i) = inner(freestream%c_hat_g, this%vertices(i)%loc)

            !! Loop through neighboring vertices
            !do j=1,this%vertices(i)%adjacent_vertices%len()
            !    call this%vertices(i)%adjacent_vertices%get(j, k)
            !    x(i) = max(inner(freestream%c_hat_g, this%vertices(k)%loc), x(i))
            !end do

        end do

        ! Get sorted indices
        call insertion_arg_sort(x, i_sorted)
        
    end subroutine surface_mesh_get_vertex_sorting_indices


    subroutine surface_mesh_rearrange_vertices(this, freestream)
        ! Rearranges the vertices to proceed in the freestream direction

        implicit none

        class(surface_mesh),intent(inout),target :: this
        type(flow),intent(in) :: freestream

        integer,dimension(:),allocatable :: i_sorted

        if (verbose) write(*,'(a)',advance='no') "     Sorting vertices for efficient matrix structure..."

        ! Get sorted indices
        call this%get_vertex_sorting_indices(freestream, i_sorted)

        ! Move the vertices around
        call this%allocate_new_vertices(0, i_sorted)

        if (verbose) write(*,*) "Done."
        
    end subroutine surface_mesh_rearrange_vertices


    subroutine surface_mesh_calc_vertex_normals(this)
        ! Initializes the normal vectors associated with each vertex.
        ! Must be called only once wake-shedding edge vertices have been cloned.

        implicit none

        class(surface_mesh),intent(inout) :: this

        real,dimension(3) :: n_avg
        integer :: i, j, i_panel

        if (verbose) write(*,'(a)',advance='no') "     Calculating vertex normals..."

        ! Loop through vertices
        !$OMP parallel do private(n_avg, i, i_panel) schedule(dynamic)
        do j=1,this%N_verts

            ! Loop through neighboring panels and compute the average of their normal vectors
            n_avg = 0
            do i=1,this%vertices(j)%panels%len()
                call this%vertices(j)%panels%get(i, i_panel)
                n_avg = n_avg + this%panels(i_panel)%n_g
            end do

            ! For vertices on the mirror plane, the component normal to the plane should be zeroed
            if (this%mirrored .and. this%vertices(j)%on_mirror_plane) then
                n_avg(this%mirror_plane) = 0.
            end if

            ! Normalize and store
            this%vertices(j)%n_g = n_avg/norm2(n_avg)

            ! Calculate mirrored normal for mirrored vertex
            if (this%mirrored) then
                this%vertices(j)%n_g_mir = mirror_across_plane(this%vertices(j)%n_g, this%mirror_plane)
            end if

            ! Calculate average edge lengths for each vertex
            call this%vertices(j)%calc_average_edge_length(this%vertices)

        end do

        if (verbose) write(*,*) "Done."

    end subroutine surface_mesh_calc_vertex_normals


    subroutine surface_mesh_init_wake(this, freestream, wake_file)
        ! Handles wake initialization

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream
        character(len=:),allocatable,intent(in) :: wake_file

        type(vtk_out) :: wake_vtk

        if (this%append_wake) then

            ! Update default Trefftz distance
            if (this%trefftz_distance < 0.) then

                ! Supersonic
                if (freestream%supersonic) then
                    call this%update_supersonic_trefftz_distance(freestream)

                ! Subsonic
                else
                    call this%update_subsonic_trefftz_distance(freestream)

                end if
            end if

            ! Initialize wake
            call this%wake%init(this%edges, this%vertices, this%N_panels, freestream, &
                                this%asym_flow, this%mirror_plane, this%N_wake_panels_streamwise, this%trefftz_distance)
        
            ! Export wake geometry
            if (wake_file /= 'none') then

                ! Clear old file
                call delete_file(wake_file)

                ! Write new geometry
                if (this%wake%N_panels > 0) then

                    call wake_vtk%begin(wake_file)
                    call wake_vtk%write_points(this%wake%vertices)
                    call wake_vtk%write_panels(this%wake%panels)
                    call wake_vtk%write_cell_normals(this%wake%panels)
                    call wake_vtk%finish()

                end if
            end if

        else
            
            ! Set parameters to let later code know the wake is not being modeled
            this%wake%N_panels = 0
            this%wake%N_verts = 0

        end if
    
    end subroutine surface_mesh_init_wake


    subroutine surface_mesh_update_supersonic_trefftz_distance(this, freestream)
        ! Determines the appropriate Trefftz distance based on the mesh geometry

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real :: dist, max_dist
        integer :: i

        ! Loop through mesh vertices, looking for the most downstream
        max_dist = 0.
        do i=1,this%N_verts

            ! Calculate distance
            dist = inner(this%vertices(i)%loc, freestream%c_hat_g)

            ! Check maximum
            if (dist > max_dist) then
                max_dist = dist
            end if

        end do

        ! Set Trefftz distance
        this%trefftz_distance = max_dist
    
    end subroutine surface_mesh_update_supersonic_trefftz_distance


    subroutine surface_mesh_update_subsonic_trefftz_distance(this, freestream)
        ! Determines the appropriate Trefftz distance based on the mesh geometry

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real :: back, front, x
        integer :: i

        ! Loop through vertices to calculate most downstream and upstream distances
        front = inner(freestream%c_hat_g, this%vertices(1)%loc)
        back = front
        do i=2,this%N_verts
            x = inner(freestream%c_hat_g, this%vertices(i)%loc)
            front = min(front, x)
            back = max(back, x)
        end do

        ! Calculate Trefftz distance
        this%trefftz_distance = 20.*abs(front-back)
    
    end subroutine surface_mesh_update_subsonic_trefftz_distance


    subroutine surface_mesh_place_interior_control_points(this, offset)

        implicit none

        class(surface_mesh),intent(inout) :: this
        real,intent(in) :: offset

        integer :: i, j, i_panel
        real,dimension(3) :: n_avg
        real :: offset_ratio

        ! Specify number of control points
        this%N_cp = this%N_verts

        ! Allocate memory
        allocate(this%cp(3,this%N_verts))

        ! Calculate offset ratio such that the control point will remain within the body based on the minimum detected angle between panels
        if (this%wake_present) then
            offset_ratio = 0.5*sqrt(0.5*(1. + this%C_min_panel_angle))
        end if

        ! Loop through vertices
        !$OMP parallel do private(j, n_avg, i_panel) schedule(dynamic) shared(this, offset, offset_ratio) default(none)
        do i=1,this%N_verts

            ! If the vertex has been cloned, it needs to be shifted off the normal slightly so that it is unique from its counterpart
            if (this%vertices(i)%clone) then

                ! Loop through panels associated with this clone to get their average normal vector
                n_avg = 0.
                do j=1,this%vertices(i)%panels_not_across_wake_edge%len()

                    ! Get panel index
                    call this%vertices(i)%panels_not_across_wake_edge%get(j, i_panel)

                    ! Add normal vector
                    n_avg = n_avg + this%panels(i_panel)%n_g

                end do

                ! Add effect of mirrored panels
                if (this%vertices(i)%on_mirror_plane) then
                    n_avg(this%mirror_plane) = 0.
                end if

                ! Normalize
                n_avg = n_avg/norm2(n_avg)

                ! Place control point
                this%cp(:,i) = this%vertices(i)%loc - offset * (this%vertices(i)%n_g - offset_ratio*n_avg)*this%vertices(i)%l_avg

            ! If it has no clone, then placement simply follows the normal vector
            else

                this%cp(:,i) = this%vertices(i)%loc - offset*this%vertices(i)%n_g*this%vertices(i)%l_avg

            end if

        end do

        ! Calculate mirrored points, if necessary
        if (this%mirrored) then

            ! Allocate memory
            allocate(this%cp_mirrored(3,this%N_cp))

            ! Calculate mirrors
            !$OMP parallel do schedule(static)
            do i=1,this%N_cp
                this%cp_mirrored(:,i) = mirror_across_plane(this%cp(:,i), this%mirror_plane)
            end do

        end if

    end subroutine surface_mesh_place_interior_control_points


    subroutine surface_mesh_output_results(this, body_file, wake_file, control_point_file, mirrored_body_file, &
                                           mirrored_control_point_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: body_file, wake_file, control_point_file
        character(len=:),allocatable,intent(in) :: mirrored_body_file, mirrored_control_point_file

        real,dimension(:),allocatable :: mu_on_wake, panel_inclinations
        type(vtk_out) :: body_vtk, wake_vtk, cp_vtk
        integer :: i

        ! Write out data for body
        if (body_file /= 'none') then

            ! Clear old file
            call delete_file(body_file)

            ! Get panel inclinations
            allocate(panel_inclinations(this%N_panels))
            do i=1,this%N_panels
                panel_inclinations(i) = this%panels(i)%r
            end do

            ! Write geometry
            call body_vtk%begin(body_file)
            call body_vtk%write_points(this%vertices)
            call body_vtk%write_panels(this%panels)

            ! Pressures
            if (allocated(this%C_p_inc)) then
                call body_vtk%write_cell_scalars(this%C_p_inc(1:this%N_panels), "C_p_inc")
            end if
            if (allocated(this%C_p_ise)) then
                call body_vtk%write_cell_scalars(this%C_p_ise(1:this%N_panels), "C_p_ise")
            end if
            if (allocated(this%C_p_2nd)) then
                call body_vtk%write_cell_scalars(this%C_p_2nd(1:this%N_panels), "C_p_2nd")
            end if
            if (allocated(this%C_p_lin)) then
                call body_vtk%write_cell_scalars(this%C_p_lin(1:this%N_panels), "C_p_lin")
            end if
            if (allocated(this%C_p_sln)) then
                call body_vtk%write_cell_scalars(this%C_p_sln(1:this%N_panels), "C_p_sln")
            end if

            ! Corrected pressures
            if (allocated(this%C_p_pg)) then
                call body_vtk%write_cell_scalars(this%C_p_pg(1:this%N_panels), "C_p_PG")
            end if
            if (allocated(this%C_p_kt)) then
                call body_vtk%write_cell_scalars(this%C_p_kt(1:this%N_panels), "C_p_KT")
            end if
            if (allocated(this%C_p_lai)) then
                call body_vtk%write_cell_scalars(this%C_p_lai(1:this%N_panels), "C_p_L")
            end if

            ! Constant sources
            if (source_order == 0) then
                call body_vtk%write_cell_scalars(this%sigma(1:this%N_panels), "sigma")
            end if

            ! Other
            call body_vtk%write_cell_scalars(panel_inclinations, "inclination")
            call body_vtk%write_cell_vectors(this%v(:,1:this%N_panels), "v")
            call body_vtk%write_cell_vectors(this%dC_f(:,1:this%N_panels), "dC_f")

            ! Linear sources
            if (source_order == 1) then
                call body_vtk%write_point_scalars(this%sigma(1:this%N_verts), "sigma")
            end if

            call body_vtk%write_point_scalars(this%mu(1:this%N_verts), "mu")
            call body_vtk%write_point_scalars(this%Phi_u(1:this%N_verts), "Phi_u")
            call body_vtk%finish()

            if (verbose) write(*,*) "    Surface results written to: ", body_file
        end if

        ! Write out data for mirrored body
        if (mirrored_body_file /= 'none' .and. this%asym_flow) then

            ! Clear old file
            call delete_file(mirrored_body_file)

            ! Get panel inclinations
            if (.not. allocated(panel_inclinations)) then
                allocate(panel_inclinations(this%N_panels))
            end if
            do i=1,this%N_panels
                panel_inclinations(i) = this%panels(i)%r_mir
            end do

            ! Write geometry
            call body_vtk%begin(mirrored_body_file)
            call body_vtk%write_points(this%vertices, this%mirror_plane)
            call body_vtk%write_panels(this%panels)

            ! Pressures
            if (allocated(this%C_p_inc)) then
                call body_vtk%write_cell_scalars(this%C_p_inc(this%N_panels+1:this%N_panels*2), "C_p_inc")
            end if
            if (allocated(this%C_p_ise)) then
                call body_vtk%write_cell_scalars(this%C_p_ise(this%N_panels+1:this%N_panels*2), "C_p_ise")
            end if
            if (allocated(this%C_p_2nd)) then
                call body_vtk%write_cell_scalars(this%C_p_2nd(this%N_panels+1:this%N_panels*2), "C_p_2nd")
            end if
            if (allocated(this%C_p_lin)) then
                call body_vtk%write_cell_scalars(this%C_p_lin(this%N_panels+1:this%N_panels*2), "C_p_lin")
            end if
            if (allocated(this%C_p_sln)) then
                call body_vtk%write_cell_scalars(this%C_p_sln(this%N_panels+1:this%N_panels*2), "C_p_sln")
            end if

            ! Corrected pressures
            if (allocated(this%C_p_pg)) then
                call body_vtk%write_cell_scalars(this%C_p_pg(this%N_panels+1:this%N_panels*2), "C_p_PG")
            end if
            if (allocated(this%C_p_kt)) then
                call body_vtk%write_cell_scalars(this%C_p_kt(this%N_panels+1:this%N_panels*2), "C_p_KT")
            end if
            if (allocated(this%C_p_lai)) then
                call body_vtk%write_cell_scalars(this%C_p_lai(this%N_panels+1:this%N_panels*2), "C_p_L")
            end if

            ! Constant sources
            if (source_order == 0) then
                call body_vtk%write_cell_scalars(this%sigma(this%N_panels+1:this%N_panels*2), "sigma")
            end if

            ! Other
            call body_vtk%write_cell_vectors(this%v(:,this%N_panels+1:this%N_panels*2), "v")
            call body_vtk%write_cell_vectors(this%dC_f(:,this%N_panels+1:this%N_panels*2), "dC_f")

            ! Linear sources
            if (source_order == 1) then
                call body_vtk%write_point_scalars(this%sigma(this%N_verts+1:this%N_verts*2), "sigma")
            end if

            call body_vtk%write_point_scalars(this%mu(this%N_cp+1:this%N_cp*2), "mu")
            call body_vtk%write_point_scalars(this%Phi_u(this%N_verts+1:this%N_verts*2), "Phi_u")
            call body_vtk%finish()

            if (verbose) write(*,*) "    Mirrored surface results written to: ", mirrored_body_file

        end if
        
        ! Write out data for wake
        if (wake_file /= 'none') then

            ! Clear old file
            call delete_file(wake_file)

            if (this%wake%N_panels > 0) then

                ! Write out geometry
                call wake_vtk%begin(wake_file)
                call wake_vtk%write_points(this%wake%vertices)
                call wake_vtk%write_panels(this%wake%panels)
                call wake_vtk%write_cell_normals(this%wake%panels)

                ! Calculate doublet strengths
                allocate(mu_on_wake(this%wake%N_verts))
                do i=1,this%wake%N_verts
                    mu_on_wake(i) = this%mu(this%wake%vertices(i)%top_parent)-this%mu(this%wake%vertices(i)%bot_parent)
                end do

                ! Write doublet strengths
                call wake_vtk%write_point_scalars(mu_on_wake, "mu")

                ! Finish up
                call wake_vtk%finish()
                if (verbose) write(*,*) "    Wake results written to: ", wake_file

            else
                if (verbose) write(*,*) "    No wake to export."

            end if
        end if
        
        ! Write out data for control points
        if (control_point_file /= 'none') then

            ! Clear old file
            call delete_file(control_point_file)

            ! Write out data
            call cp_vtk%begin(control_point_file)
            call cp_vtk%write_points(this%cp)
            call cp_vtk%write_vertices(this%N_cp)
            call cp_vtk%write_point_scalars(this%phi_cp(1:this%N_cp), "phi")
            call cp_vtk%write_point_scalars(this%phi_cp_mu(1:this%N_cp), "phi_mu")
            call cp_vtk%write_point_scalars(this%phi_cp_sigma(1:this%N_cp), "phi_sigma")
            call cp_vtk%finish()

            if (verbose) write(*,*) "    Control point results written to: ", control_point_file
        end if
        
        ! Write out data for mirrored control points
        if (mirrored_control_point_file /= 'none' .and. this%asym_flow) then

            ! Clear old file
            call delete_file(mirrored_control_point_file)

            ! Write out data
            call cp_vtk%begin(mirrored_control_point_file)
            call cp_vtk%write_points(this%cp_mirrored)
            call cp_vtk%write_vertices(this%N_cp)
            call cp_vtk%write_point_scalars(this%phi_cp(this%N_cp+1:this%N_cp*2), "phi")
            call cp_vtk%write_point_scalars(this%phi_cp_mu(this%N_cp+1:this%N_cp*2), "phi_mu")
            call cp_vtk%write_point_scalars(this%phi_cp_sigma(this%N_cp+1:this%N_cp*2), "phi_sigma")
            call cp_vtk%finish()

            if (verbose) write(*,*) "    Mirrored control point results written to: ", mirrored_control_point_file
        end if
    
    end subroutine surface_mesh_output_results


end module surface_mesh_mod