! A surface mesh type encompassing a body
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

        integer :: N_verts, N_panels, N_cp, N_edges, N_true_verts ! as in, not midpoints
        integer :: N_subinc, N_supinc
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
        real,dimension(:,:),allocatable :: V_cells, V_verts_avg, V_verts_std, dC_f ! Surface velocities and pressure forces
        real :: control_point_offset
        logical :: mirrored ! Whether the mesh is to be mirrored about any planes
        integer :: mirror_plane ! Index of the plane across which the mesh is mirrored (1: yz, 2: xz, 3: xy); this is the index of the normal to that plane
        logical :: asym_flow ! Whether the flow is asymmetric about the mirror plane
        logical :: found_discontinuous_edges, midpoints_created
        real,dimension(:),allocatable :: mu, sigma ! Singularity strengths
        real :: S_ref ! Reference parameters

        contains

            ! Basic initialization
            procedure :: init => surface_mesh_init
            procedure :: locate_adjacent_panels => surface_mesh_locate_adjacent_panels
            procedure :: store_adjacent_vertices => surface_mesh_store_adjacent_vertices
            procedure :: check_panels_adjacent => surface_mesh_check_panels_adjacent

            ! Initialization based on flow properties
            procedure :: init_with_flow => surface_mesh_init_with_flow
            procedure :: count_panel_inclinations => surface_mesh_count_panel_inclinations
            procedure :: characterize_edges => surface_mesh_characterize_edges
            procedure :: clone_vertices => surface_mesh_clone_vertices
            procedure :: init_vertex_clone => surface_mesh_init_vertex_clone
            procedure :: set_up_mirroring => surface_mesh_set_up_mirroring
            procedure :: rearrange_vertices_streamwise => surface_mesh_rearrange_vertices_streamwise
            procedure :: shuffle_midpoints_in => surface_mesh_shuffle_midpoints_in
            procedure :: calc_vertex_normals => surface_mesh_calc_vertex_normals
            procedure :: calc_point_dod => surface_mesh_calc_point_dod

            ! Helpers
            procedure :: find_vertices_on_mirror => surface_mesh_find_vertices_on_mirror
            procedure :: get_indices_to_panel_vertices => surface_mesh_get_indices_to_panel_vertices
            procedure :: allocate_new_vertices => surface_mesh_allocate_new_vertices
            procedure :: get_vertex_sorting_indices => surface_mesh_get_vertex_sorting_indices

            ! Wake stuff
            procedure :: init_wake => surface_mesh_init_wake
            procedure :: update_supersonic_trefftz_distance => surface_mesh_update_supersonic_trefftz_distance
            procedure :: update_subsonic_trefftz_distance => surface_mesh_update_subsonic_trefftz_distance

            ! Control points
            procedure :: place_interior_control_points => surface_mesh_place_interior_control_points

            ! Post-processing
            procedure :: get_induced_potentials_at_point => surface_mesh_get_induced_potentials_at_point
            procedure :: output_results => surface_mesh_output_results
            procedure :: write_body => surface_mesh_write_body
            procedure :: write_body_mirror => surface_mesh_write_body_mirror
            procedure :: write_control_points => surface_mesh_write_control_points
            procedure :: write_mirrored_control_points => surface_mesh_write_mirrored_control_points

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

        ! Initialize a few things
        this%midpoints_created = .false.

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

        ! Store some mesh parameters
        this%N_true_verts = this%N_verts

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
            call this%vertices(i)%set_whether_on_mirror_plane(this%mirror_plane)
        end do
    
    end subroutine surface_mesh_find_vertices_on_mirror


    subroutine surface_mesh_locate_adjacent_panels(this)
        ! Loops through panels to determine which are adjacent

        implicit none

        class(surface_mesh),intent(inout),target :: this

        integer :: i, j, m, n, m1, n1, temp, i_edge, N_edges, i_mid, N_orig_verts, edge_index_i, edge_index_j
        logical :: already_found_shared
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
        !$OMP parallel private(j, already_found_shared, distance, i_endpoints, m, m1, n, n1, temp, i_edge) &
        !$OMP & private(edge_index_i, edge_index_j)

        !$OMP do schedule(dynamic)
        do i=1,this%N_panels

            ! For each panel, check if it abuts the mirror plane
            if (this%mirrored) then

                ! Perform check
                if (this%panels(i)%check_abutting_mirror_plane(this%N_panels, i_endpoints, edge_index_i)) then

                    !$OMP critical
                    
                    ! Update number of edges
                    this%N_edges = this%N_edges + 1
                    i_edge = this%N_edges

                    ! Store adjacent vertices
                    call this%store_adjacent_vertices(i_endpoints, i_edge)

                    ! Store edge index for this panel
                    edge_index1(i_edge) = edge_index_i

                    ! Store in arrays for later storage in edge objects
                    panel1(i_edge) = i
                    panel2(i_edge) = i + this%N_panels
                    vertex1(i_edge) = i_endpoints(1)
                    vertex2(i_edge) = i_endpoints(2)
                    on_mirror_plane(i_edge) = .true.
                    edge_index2(i_edge) = 0 ! Just a placeholder since the second panel doesn't technically exist
                    edge_length(i_edge) = dist(this%vertices(i_endpoints(1))%loc, this%vertices(i_endpoints(2))%loc)

                    !$OMP end critical

                end if

            end if

            ! Loop through each potential neighbor
            neighbor_loop: do j=i+1,this%N_panels

                ! Check if we've found all neighbors for this panel
                if (all(this%panels(i)%abutting_panels /= 0)) exit neighbor_loop

                ! Check if these are abutting
                if (this%check_panels_adjacent(i, j, i_endpoints, edge_index_i, edge_index_j)) then

                    !$OMP critical
                    
                    ! Update number of edges
                    this%N_edges = this%N_edges + 1
                    i_edge = this%N_edges

                    ! Store vertices being adjacent to one another
                    call this%store_adjacent_vertices(i_endpoints, i_edge)

                    ! Store adjacent panels and panel edges
                    edge_index1(i_edge) = edge_index_i
                    edge_index2(i_edge) = edge_index_j

                    ! Store information in arrays for later storage in edge objects
                    panel1(i_edge) = i
                    panel2(i_edge) = j
                    vertex1(i_edge) = i_endpoints(1)
                    vertex2(i_edge) = i_endpoints(2)
                    edge_length(i_edge) = dist(this%vertices(i_endpoints(1))%loc, this%vertices(i_endpoints(2))%loc)

                    !$OMP end critical

                end if

            end do neighbor_loop

        end do

        ! Check for panels abutting empty space and add those edges.
        !$OMP single
        do i=1,this%N_panels

            ! Check for an edge with no abutting panel
            do j=1,this%panels(i)%N
                if (this%panels(i)%abutting_panels(j) == 0) then

                    ! Warn that the mesh is not watertight
                    if (run_checks) then
                        write(*,*) "!!! Panel ", i, " is missing a neighbor, meaning the supplied mesh may not be watertight."
                        write(*,*) "!!! Solution quality may be adversely affected."
                    end if

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
                call this%vertices(i_mid)%panels_not_across_wake_edge%append(panel1(i))
                if (panel2(i) > 0 .and. panel2(i) <= this%N_panels) then
                    call this%vertices(i_mid)%panels%append(panel2(i))
                    call this%vertices(i_mid)%panels_not_across_wake_edge%append(panel2(i))
                end if

                ! Add edge
                call this%vertices(i_mid)%adjacent_edges%append(i)

                ! Check if midpoint is on mirror plane
                if (this%mirrored) then
                    call this%vertices(i_mid)%set_whether_on_mirror_plane(this%mirror_plane)
                end if

                ! Point edge to it
                this%edges(i)%i_midpoint = i_mid

                ! Point panels to it
                this%panels(panel1(i))%midpoints(edge_index1(i))%ptr => this%vertices(i_mid)

                if (panel2(i) > 0 .and. panel2(i) <= this%N_panels) then
                    this%panels(panel2(i))%midpoints(edge_index2(i))%ptr => this%vertices(i_mid)
                end if

            end do

            ! Set flag
            this%midpoints_created = .true.

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


    function surface_mesh_check_panels_adjacent(this, i, j, i_endpoints, edge_index_i, edge_index_j) result(adjacent)
        ! Checks whether panels i and j are adjacent

        class(surface_mesh),intent(inout) :: this
        integer,intent(in) :: i, j
        integer,dimension(2),intent(out) :: i_endpoints
        integer,intent(out) :: edge_index_i, edge_index_j

        logical :: adjacent

        logical :: already_found_shared
        integer :: m, n, m1, n1, temp

        ! Initialize
        adjacent = .false.

        ! Initialize for this panel pair
        already_found_shared = .false.

        ! Check if the panels are abutting
        abutting_loop: do m=1,this%panels(i)%N
            do n=1,this%panels(j)%N

                ! Check if vertices have the same index
                if (this%panels(i)%get_vertex_index(m) == this%panels(j)%get_vertex_index(n)) then

                    ! Previously found a shared vertex, so we have abutting panels
                    if (already_found_shared) then

                        adjacent = .true.

                        ! Store second shared vertex
                        i_endpoints(2) = this%panels(i)%get_vertex_index(m)

                        ! Check order; edge should proceed counterclockwise about panel i
                        if (m1 == 1 .and. m == this%panels(i)%N) then
                            temp = i_endpoints(1)
                            i_endpoints(1) = i_endpoints(2)
                            i_endpoints(2) = temp
                        end if

                        ! Store adjacent panels and panel edges
                        ! This stores the adjacent panels and edges according to the index of that edge
                        ! for the current panel

                        ! Store that i is adjacent to j
                        ! This one is more complicated because we don't know that n1 will be less than n; just the nature of the nested loop.
                        ! Basically, if one is 1 and the other is N, then we're dealing with edge N for panel j.
                        ! Otherwise, we're dealing with abs(n1-n) being 1, meaning edge min(n1, n).
                        if ( (n1 == 1 .and. n == this%panels(j)%N) .or. (n == 1 .and. n1 == this%panels(j)%N) ) then
                            this%panels(j)%abutting_panels(this%panels(j)%N) = i
                            edge_index_j = this%panels(j)%N
                        else
                            n1 = min(n, n1)
                            this%panels(j)%abutting_panels(n1) = i
                            edge_index_j = n1
                        end if

                        ! Store that j is adjacent to i
                        if (m1 == 1 .and. m == this%panels(i)%N) then ! Nth edge
                            this%panels(i)%abutting_panels(m) = j
                            edge_index_i = m
                        else ! 1st or 2nd edge
                            this%panels(i)%abutting_panels(m1) = j
                            edge_index_i = m1
                        end if

                        return

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
    
    end function surface_mesh_check_panels_adjacent


    subroutine surface_mesh_init_with_flow(this, freestream, body_file, wake_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream
        character(len=:),allocatable,intent(in) :: body_file, wake_file

        integer :: i
        real,dimension(:),allocatable :: panel_inclinations
        type(vtk_out) :: body_vtk

        ! Check flow symmetry condition
        this%asym_flow = .false.
        if (this%mirrored) then
            if (.not. freestream%sym_about(this%mirror_plane)) then
                this%asym_flow = .true.
            end if
        end if

        ! Initialize properties of panels dependent upon the flow
        if (verbose) write(*,'(a)',advance='no') "     Calculating panel properties..."
        !$OMP parallel do schedule(static)
        do i=1,this%N_panels
            call this%panels(i)%init_with_flow(freestream, this%asym_flow, this%mirror_plane)
        end do
        if (verbose) write(*,*) "Done."

        ! Determine number of sub- and superinclined panels
        call this%count_panel_inclinations()

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

        ! Calculate normals (for placing control points)
        call this%calc_vertex_normals()

        ! Initialize wake
        call this%init_wake(freestream, wake_file)

        ! Set up influencing vertex arrays
        if (verbose) write(*,"(a)",advance='no') "     Setting influencing vertices for panels..."
        !$OMP parallel do schedule(static)
        do i=1,this%N_panels
            call this%panels(i)%set_influencing_verts()
        end do

        if (verbose) write(*,*) "Done."

        ! Write out body file to ensure it's been parsed correctly
        if (body_file /= 'none') then
            call this%write_body(body_file, .false.)
        end if
    
    end subroutine surface_mesh_init_with_flow


    subroutine surface_mesh_count_panel_inclinations(this)
        ! Counts the number of sub- and superinclined panels

        implicit none
        
        class(surface_mesh), intent(inout) :: this

        integer :: i
    
        this%N_subinc = 0
        this%N_supinc = 0
        do i=1,this%N_panels

            ! Original panel
            if (this%panels(i)%r > 0.) then
                this%N_subinc = this%N_subinc + 1
            else
                this%N_supinc = this%N_supinc + 1
            end if

            ! Mirrored panel
            if (this%asym_flow) then
                if (this%panels(i)%r_mir > 0.) then
                    this%N_subinc = this%N_subinc + 1
                else
                    this%N_supinc = this%N_supinc + 1
                end if
            end if

        end do
        
    end subroutine surface_mesh_count_panel_inclinations


    subroutine surface_mesh_characterize_edges(this, freestream)
        ! Locates wake-shedding edges and supersonic/subsonic edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        integer :: i, j, k, m, n, temp, top_panel, bottom_panel, i_vert_1, i_vert_2, mid, N_wake_edges
        real :: C_angle, C_min_angle
        real,dimension(3) :: second_normal, cross_result
        real,dimension(3) :: t_hat_g, d

        if (verbose) write(*,'(a)',advance='no') "     Characterizing edges..."

        ! Initialize
        N_wake_edges = 0
        this%found_discontinuous_edges = .false.

        !$OMP parallel private(i, j, second_normal, C_angle, i_vert_1, i_vert_2, mid) &
        !$OMP & private(cross_result, d, t_hat_g) &
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
                    
                    ! Get vertex indices (simplifies later code)
                    i_vert_1 = this%edges(k)%verts(1)
                    i_vert_2 = this%edges(k)%verts(2)
                    
                    ! Calculate tangent in global coords
                    d = this%vertices(i_vert_2)%loc - this%vertices(i_vert_1)%loc
                    t_hat_g = d / norm2(d)
                    
                    cross_result = cross(this%panels(i)%n_g, second_normal)
                    
                    ! Check sign between normal vectors cross product and edge tangent
                    if (inner(cross_result, t_hat_g) > 0.) then
                        
                        ! Having passed the previous three checks, we've found a wake-shedding edge
                        
                        ! Set that we've found discontinuous edges
                        this%found_discontinuous_edges = .true.
                        
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
                if (this%midpoints_created) then
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

            ! Fix wake partners
            do i=1,this%N_verts
                temp_vertices(i)%i_wake_partner = i_rearrange_inv(temp_vertices(i)%i_wake_partner)
            end do

            ! Fix edge endpoints (and midpoints if necessary)
            do i=1,this%N_edges

                ! Endpoints
                this%edges(i)%verts(1) = i_rearrange_inv(this%edges(i)%verts(1))
                this%edges(i)%verts(2) = i_rearrange_inv(this%edges(i)%verts(2))

                ! Midpoints
                if (doublet_order == 2) then
                    this%edges(i)%i_midpoint = i_rearrange_inv(this%edges(i)%i_midpoint)
                end if

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
                if (this%midpoints_created) then
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

        class(surface_mesh),intent(inout) :: this

        integer :: i, j, N_clones, i_jango, i_boba
        integer,dimension(:),allocatable :: i_rearrange, i_rearrange_inv

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

            ! Allocate rearranged indices array
            allocate(i_rearrange(this%N_verts), source=0)
            allocate(i_rearrange_inv(this%N_verts), source=0)

            ! Initialize clones
            j = 0
            do i_jango=1,this%N_verts-N_clones ! Only need to loop through original vertices here

                ! Check if this vertex needs to be cloned
                if (this%vertices(i_jango)%clone) then

                    ! Get index for the clone
                    j = j + 1
                    i_boba = this%N_verts - N_clones + j ! Will be at position N_verts-N_clones+j in the new vertex array

                    ! Get rearranged indices
                    i_rearrange_inv(i_jango) = i_jango + j - 1
                    i_rearrange_inv(i_boba) = i_jango + j

                    ! Initialize clone
                    call this%init_vertex_clone(i_jango, i_boba)

                else

                    ! Get rearranged indices
                    i_rearrange_inv(i_jango) = i_jango + j

                    ! If this vertex did not need to be cloned, but it is on the mirror plane and its mirror is unique
                    ! then the wake strength will be determined by its mirror as well in the case of an asymmetric flow.
                    if (this%mirrored .and. this%asym_flow .and.  this%vertices(i_jango)%on_mirror_plane .and. &
                        this%vertices(i_jango)%mirrored_is_unique) then

                        this%vertices(i_jango)%i_wake_partner = i_jango + this%N_verts

                    end if

                end if

            end do

            ! Get inverse mapping
            do i=1,this%N_verts
                i_rearrange(i_rearrange_inv(i)) = i
            end do

            ! Rearrange vertices so clones are next to clones
            call this%allocate_new_vertices(0, i_rearrange)

            if (verbose) write(*,'(a, i4, a, i7, a)') "Done. Cloned ", N_clones, " vertices. Mesh now has ", &
                                                      this%N_verts, " vertices."

        end if

    end subroutine surface_mesh_clone_vertices


    subroutine surface_mesh_init_vertex_clone(this, i_jango, i_boba)
        ! Clones the vertex at i_jango into i_boba

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        integer,intent(in) :: i_jango, i_boba

        integer :: k, m, n, i_bot_panel, i_abutting_panel, i_adj_vert, i_top_panel, i_edge

        ! Basic initialization
        call this%vertices(i_boba)%init(this%vertices(i_jango)%loc, i_boba, this%vertices(i_jango)%vert_type)
        this%vertices(i_boba)%clone = .true.

        ! If this is a true vertex (not a midpoint), add to the count
        if (this%vertices(i_jango)%vert_type==1) this%N_true_verts = this%N_true_verts + 1

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

                            ! Make sure the abutting panel is not a mirrored/nonexistant panel
                            if (i_abutting_panel > 0 .and. i_abutting_panel <= this%N_panels) then

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
        
    end subroutine surface_mesh_init_vertex_clone


    subroutine surface_mesh_get_vertex_sorting_indices(this, freestream, i_sorted)
        ! Determines the vertex indices which will sort them in the conpressibility direction

        implicit none

        class(surface_mesh),intent(in) :: this
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(out) :: i_sorted

        real,dimension(:),allocatable :: x
        integer :: i

        ! Allocate the compressibility distance array
        allocate(x(this%N_verts))

        ! Add compressibility distance of each vertex
        do i=1,this%N_verts

            ! Initialize with the compressibility distance of this vertex
            x(i) = -inner(freestream%c_hat_g, this%vertices(i)%loc)

        end do

        ! Get sorted indices
        call insertion_arg_sort(x, i_sorted)
        
    end subroutine surface_mesh_get_vertex_sorting_indices


    subroutine surface_mesh_rearrange_vertices_streamwise(this, freestream)
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
        
    end subroutine surface_mesh_rearrange_vertices_streamwise


    subroutine surface_mesh_shuffle_midpoints_in(this)
        ! Rearranges the vertices (originals and midpoints) so that midpoints are in a better location in the array

        class(surface_mesh), intent(inout) :: this

        integer,dimension(:),allocatable :: i_sorted
        real,dimension(:),allocatable :: i_near_vert
        integer :: N_inserted, i, j

        if (verbose) write(*,'(a)',advance='no') "     Sorting midpoints into the vertex array..."

        ! Allocate index array
        allocate(i_near_vert(this%N_verts))

        ! Loop through midpoints
        N_inserted = 0
        do i=1,this%N_verts

            ! Check this is a midpoint
            if (this%vertices(i)%vert_type==2) then

                ! Get first neighboring vertex
                call this%vertices(i)%adjacent_vertices%get(1, j)
                i_near_vert(i) = j

            ! If it's not a midpoint, it's its own near neighbor
            else
                i_near_vert(i) = i
            end if
        end do

        ! Sort
        call insertion_arg_sort(i_near_vert, i_sorted)

        ! Move the vertices around
        call this%allocate_new_vertices(0, i_sorted)

        if (verbose) write(*,*) "Done."
    
    end subroutine surface_mesh_shuffle_midpoints_in


    subroutine surface_mesh_calc_vertex_normals(this)
        ! Initializes the normal vectors associated with each vertex.
        ! Must be called only once wake-shedding edge vertices have been cloned.

        implicit none

        class(surface_mesh),intent(inout) :: this

        real,dimension(3) :: n_avg, n_avg_plane
        integer :: i, j, i_panel
        real,dimension(:),allocatable :: constraint_array
        real :: x

        if (verbose) write(*,'(a)',advance='no') "     Calculating vertex normals..."

        ! Loop through vertices
        !$OMP parallel do private(n_avg, i, i_panel, constraint_array, n_avg_plane) schedule(dynamic)
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

            ! Get the inner products with the panel normals and the computed vertex normal
            allocate(constraint_array(this%vertices(j)%panels%len()))
            do i=1,this%vertices(j)%panels%len()
                call this%vertices(j)%panels%get(i, i_panel)
                constraint_array(i) = inner(this%panels(i_panel)%n_g, n_avg)
            end do

            ! Check whether the vertex normal points wholly outside the mesh
            do while (any(constraint_array < 0.))
                ! Get the normal vector describing the average plane for the panels for which the constraint is not satisfied
                n_avg_plane = 0.
                do i=1,this%vertices(j)%panels%len()

                    ! Loop through neighboring panels
                    if (constraint_array(i) < 0.) then
                        call this%vertices(j)%panels%get(i, i_panel)
                        n_avg_plane = n_avg_plane + this%panels(i_panel)%n_g
                    end if
                    
                end do

                ! Normalize
                n_avg_plane = n_avg_plane/norm2(n_avg_plane)

                ! Calculate new vertex normal
                n_avg = n_avg - 2.*inner(n_avg, n_avg_plane)*n_avg_plane

                ! Recalculate constraint array
                do i=1,this%vertices(j)%panels%len()
                    call this%vertices(j)%panels%get(i, i_panel)
                    constraint_array(i) = inner(this%panels(i_panel)%n_g, n_avg)
                end do

            end do

            deallocate(constraint_array)

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


    function surface_mesh_calc_point_dod(this, point, freestream) result(dod_info)
        ! Calculates the domain of dependence for the point

        implicit none
        
        class(surface_mesh),intent(in) :: this
        real,dimension(3),intent(in) :: point
        type(flow),intent(in) :: freestream

        type(dod),dimension(:),allocatable :: dod_info

        logical,dimension(:),allocatable :: verts_in_dod
        integer :: k, stat
        real,dimension(3) :: vert_loc, mirrored_vert_loc

        ! Allocate DoD storage
        if (this%mirrored) then
            allocate(dod_info(2*this%N_panels), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(verts_in_dod(2*this%N_verts), source=.false., stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        else
            allocate(dod_info(this%N_panels), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(verts_in_dod(this%N_verts), source=.false., stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        end if

        ! If the freestream is supersonic, calculate domain of dependence info
        if (freestream%supersonic) then

            ! Loop through body vertices
            do k=1,this%N_verts

                ! Get this vertex location
                vert_loc = this%vertices(k)%loc

                ! Check if this vertex is in the DoD
                verts_in_dod(k) = freestream%point_in_dod(vert_loc, point)

                if (this%mirrored) then

                    ! Get the mirrored vertex location
                    mirrored_vert_loc = mirror_across_plane(vert_loc, this%mirror_plane)

                    ! Check if this vertex is in the DoD
                    verts_in_dod(k+this%N_verts) = freestream%point_in_dod(mirrored_vert_loc, point)

                end if
            end do

            ! Loop through body panels
            do k=1,this%N_panels

                ! Check if the panel is in the DoD
                dod_info(k) = this%panels(k)%check_dod(point, freestream, verts_in_dod)

                if (this%mirrored) then

                    ! Check DoD for mirrored panel
                    dod_info(k+this%N_panels) = this%panels(k)%check_dod(point, freestream, verts_in_dod, &
                                                                         .true., this%mirror_plane)
                end if
            end do

        end if
        
    end function surface_mesh_calc_point_dod


    subroutine surface_mesh_init_wake(this, freestream, wake_file)
        ! Handles wake initialization

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream
        character(len=:),allocatable,intent(in) :: wake_file

        logical :: dummy

        if (this%append_wake .and. this%found_discontinuous_edges) then

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
                call this%wake%write_wake(wake_file, dummy)
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
            max_dist = max(dist, max_dist)

            ! Check for mirror
            if (this%asym_flow) then

                ! Calculate distance
                dist = inner(mirror_across_plane(this%vertices(i)%loc, this%mirror_plane), freestream%c_hat_g)

                ! Check maximum
                max_dist = max(dist, max_dist)

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
                !this%cp(:,i) = this%vertices(i)%loc - offset*this%vertices(i)%n_g*this%vertices(i)%l_avg*3.

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


    subroutine surface_mesh_get_induced_potentials_at_point(this, point, freestream, phi_d, phi_s)
        ! Calculates the source- and doublet-induced potentials at the given point

        implicit none
        
        class(surface_mesh),intent(in) :: this
        real,dimension(3),intent(in) :: point
        type(flow),intent(in) :: freestream
        real,intent(out) :: phi_d, phi_s

        integer :: i, j, k
        real :: phi_d_panel, phi_s_panel
        logical,dimension(:),allocatable :: verts_in_dod
        type(dod),dimension(:),allocatable :: dod_info

        ! Calculate domain of dependence
        dod_info = this%calc_point_dod(point, freestream)

        ! Loop through panels
        phi_s = 0.
        phi_d = 0.
        do k=1,this%N_panels

            ! Check DoD
            if (dod_info(k)%in_dod) then
            
                ! Calculate influence
                call this%panels(k)%calc_potentials(point, freestream, dod_info(k), .false., &
                                                    this%sigma, this%mu, phi_s_panel, phi_d_panel)
                phi_d = phi_d + phi_d_panel
                phi_s = phi_s + phi_s_panel

            end if

            ! Calculate mirrored influences
            if (this%mirrored) then

                if (dod_info(k+this%N_panels)%in_dod) then

                    if (this%asym_flow) then

                        ! Get influence of mirrored panel
                        call this%panels(k)%calc_potentials(point, freestream, dod_info(k+this%N_panels), .true., &
                                                            this%sigma, this%mu, phi_s_panel, phi_d_panel)
                        phi_d = phi_d + phi_d_panel
                        phi_s = phi_s + phi_s_panel

                    else

                        ! Get influence of panel on mirrored point
                        call this%panels(k)%calc_potentials(mirror_across_plane(point, this%mirror_plane), freestream, &
                                                            dod_info(k+this%N_panels), .false., &
                                                            this%sigma, this%mu, phi_s_panel, phi_d_panel)
                        phi_d = phi_d + phi_d_panel
                        phi_s = phi_s + phi_s_panel

                    end if

                end if

            end if

        end do
        
    end subroutine surface_mesh_get_induced_potentials_at_point


    subroutine surface_mesh_output_results(this, body_file, wake_file, control_point_file, mirrored_body_file, &
                                           mirrored_control_point_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: body_file, wake_file, control_point_file
        character(len=:),allocatable,intent(in) :: mirrored_body_file, mirrored_control_point_file

        logical :: wake_exported

        ! Write out data for body
        if (body_file /= 'none') then
            call this%write_body(body_file, solved=.true.)
            if (verbose) write(*,'(a30 a)') "    Surface: ", body_file
        end if

        ! Write out data for mirrored body
        if (mirrored_body_file /= 'none' .and. this%asym_flow) then
            call this%write_body_mirror(mirrored_body_file, solved=.true.)
            if (verbose) write(*,'(a30 a)') "    Mirrored surface: ", mirrored_body_file
        end if
        
        ! Write out data for wake
        if (wake_file /= 'none') then
            call this%wake%write_wake(wake_file, wake_exported, this%mu)

            if (wake_exported) then
                if (verbose) write(*,'(a30 a)') "    Wake: ", wake_file
            else
                if (verbose) write(*,'(a30 a)') "    Wake: ", "no wake to export"
            end if
        end if
        
        ! Write out data for control points
        if (control_point_file /= 'none') then
            call this%write_control_points(control_point_file, solved=.true.)

            if (verbose) write(*,'(a30 a)') "    Control points: ", control_point_file
        end if
        
        ! Write out data for mirrored control points
        if (mirrored_control_point_file /= 'none' .and. this%asym_flow) then
            call this%write_mirrored_control_points(mirrored_control_point_file, solved=.true.)

            if (verbose) write(*,'(a30 a)') "    Mirrored control points: ", mirrored_control_point_file
        end if
    
    end subroutine surface_mesh_output_results


    subroutine surface_mesh_write_body(this, body_file, solved)
        ! Writes the body and results (if solved) out to file

        implicit none
        
        class(surface_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: body_file
        logical,intent(in) :: solved

        type(vtk_out) :: body_vtk
        integer :: i, N_cells
        real,dimension(:),allocatable :: panel_inclinations
        real,dimension(:,:),allocatable :: cents

        ! Clear old file
        call delete_file(body_file)

        ! Determine number of cells to export
        if (doublet_order == 2) then
            N_cells = this%N_panels*4
        else
            N_cells = this%N_panels
        end if

        ! Get panel inclinations and centroids
        allocate(panel_inclinations(this%N_panels))
        allocate(cents(3,this%N_panels))
        do i=1,this%N_panels
            panel_inclinations(i) = this%panels(i)%r
            cents(:,i) = this%panels(i)%centr
        end do

        ! Write geometry
        call body_vtk%begin(body_file)
        call body_vtk%write_points(this%vertices)
        call body_vtk%write_panels(this%panels, subdivide=doublet_order==2, mirror=.false.)
        call body_vtk%write_cell_normals(this%panels)
        call body_vtk%write_cell_scalars(panel_inclinations, "inclination", .true.)
        call body_vtk%write_cell_vectors(cents, "centroid", .true.)

        if (solved) then

            ! Pressures
            if (allocated(this%C_p_inc)) then
                call body_vtk%write_cell_scalars(this%C_p_inc(1:N_cells), "C_p_inc", .false.)
            end if
            if (allocated(this%C_p_ise)) then
                call body_vtk%write_cell_scalars(this%C_p_ise(1:N_cells), "C_p_ise", .false.)
            end if
            if (allocated(this%C_p_2nd)) then
                call body_vtk%write_cell_scalars(this%C_p_2nd(1:N_cells), "C_p_2nd", .false.)
            end if
            if (allocated(this%C_p_lin)) then
                call body_vtk%write_cell_scalars(this%C_p_lin(1:N_cells), "C_p_lin", .false.)
            end if
            if (allocated(this%C_p_sln)) then
                call body_vtk%write_cell_scalars(this%C_p_sln(1:N_cells), "C_p_sln", .false.)
            end if

            ! Corrected pressures
            if (allocated(this%C_p_pg)) then
                call body_vtk%write_cell_scalars(this%C_p_pg(1:N_cells), "C_p_PG", .false.)
            end if
            if (allocated(this%C_p_kt)) then
                call body_vtk%write_cell_scalars(this%C_p_kt(1:N_cells), "C_p_KT", .false.)
            end if
            if (allocated(this%C_p_lai)) then
                call body_vtk%write_cell_scalars(this%C_p_lai(1:N_cells), "C_p_L", .false.)
            end if

            ! Constant sources
            if (source_order == 0) then
                call body_vtk%write_cell_scalars(this%sigma(1:this%N_panels), "sigma", .true.)
            end if

        end if

        ! Other
        if (solved) then
            call body_vtk%write_cell_vectors(this%V_cells(:,1:N_cells), "v", .false.)
            call body_vtk%write_cell_vectors(this%dC_f(:,1:N_cells), "dC_f", .false.)

            ! Linear sources
            if (source_order == 1) then
                call body_vtk%write_point_scalars(this%sigma(1:this%N_verts), "sigma")
            end if

            call body_vtk%write_point_scalars(this%mu(1:this%N_verts), "mu")
            call body_vtk%write_point_scalars(this%Phi_u(1:this%N_verts), "Phi_u")

            ! Quadratic doublets
            if (doublet_order == 2) then
                call body_vtk%write_point_vectors(this%V_verts_avg(:,1:this%N_verts), "v_avg")
                call body_vtk%write_point_vectors(this%V_verts_std(:,1:this%N_verts), "v_std_dev")
            end if
        end if

        ! Finalize
        call body_vtk%finish()
    
    end subroutine surface_mesh_write_body


    subroutine surface_mesh_write_body_mirror(this, mirrored_body_file, solved)
        ! Writes the body mirror and results (if solved) out to file

        implicit none
        
        class(surface_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: mirrored_body_file
        logical,intent(in) :: solved

        type(vtk_out) :: body_vtk
        integer :: i, N_cells
        real,dimension(:),allocatable :: panel_inclinations
        real,dimension(:,:),allocatable :: cents

        ! Clear old file
        call delete_file(mirrored_body_file)

        ! Determine number of cells to export
        if (doublet_order == 2) then
            N_cells = this%N_panels*4
        else
            N_cells = this%N_panels
        end if

        ! Get panel inclinations
        if (.not. allocated(panel_inclinations)) then
            allocate(panel_inclinations(this%N_panels))
            allocate(cents(3,this%N_panels))
        end if
        do i=1,this%N_panels
            panel_inclinations(i) = this%panels(i)%r_mir
            cents(:,i) = this%panels(i)%centr_mir
        end do

        ! Write geometry
        call body_vtk%begin(mirrored_body_file)
        call body_vtk%write_points(this%vertices, this%mirror_plane)
        call body_vtk%write_panels(this%panels, subdivide=doublet_order==2, mirror=.true.)
        call body_vtk%write_cell_normals(this%panels, this%mirror_plane)
        call body_vtk%write_cell_scalars(panel_inclinations, "inclination", .true.)
        call body_vtk%write_cell_vectors(cents, "centroid", .true.)

        ! Pressures
        if (allocated(this%C_p_inc)) then
            call body_vtk%write_cell_scalars(this%C_p_inc(N_cells+1:N_cells*2), "C_p_inc", .false.)
        end if
        if (allocated(this%C_p_ise)) then
            call body_vtk%write_cell_scalars(this%C_p_ise(N_cells+1:N_cells*2), "C_p_ise", .false.)
        end if
        if (allocated(this%C_p_2nd)) then
            call body_vtk%write_cell_scalars(this%C_p_2nd(N_cells+1:N_cells*2), "C_p_2nd", .false.)
        end if
        if (allocated(this%C_p_lin)) then
            call body_vtk%write_cell_scalars(this%C_p_lin(N_cells+1:N_cells*2), "C_p_lin", .false.)
        end if
        if (allocated(this%C_p_sln)) then
            call body_vtk%write_cell_scalars(this%C_p_sln(N_cells+1:N_cells*2), "C_p_sln", .false.)
        end if

        ! Corrected pressures
        if (allocated(this%C_p_pg)) then
            call body_vtk%write_cell_scalars(this%C_p_pg(N_cells+1:N_cells*2), "C_p_PG", .false.)
        end if
        if (allocated(this%C_p_kt)) then
            call body_vtk%write_cell_scalars(this%C_p_kt(N_cells+1:N_cells*2), "C_p_KT", .false.)
        end if
        if (allocated(this%C_p_lai)) then
            call body_vtk%write_cell_scalars(this%C_p_lai(N_cells+1:N_cells*2), "C_p_L", .false.)
        end if

        ! Constant sources
        if (source_order == 0) then
            call body_vtk%write_cell_scalars(this%sigma(this%N_panels+1:this%N_panels*2), "sigma", .true.)
        end if

        ! Other
        call body_vtk%write_cell_vectors(this%V_cells(:,N_cells+1:N_cells*2), "v", .false.)
        call body_vtk%write_cell_vectors(this%dC_f(:,N_cells+1:N_cells*2), "dC_f", .false.)

        ! Linear sources
        if (source_order == 1) then
            call body_vtk%write_point_scalars(this%sigma(this%N_verts+1:this%N_verts*2), "sigma")
        end if

        call body_vtk%write_point_scalars(this%mu(this%N_cp+1:this%N_cp*2), "mu")
        call body_vtk%write_point_scalars(this%Phi_u(this%N_verts+1:this%N_verts*2), "Phi_u")

        ! Quadratic doublets
        if (doublet_order == 2) then
            call body_vtk%write_point_vectors(this%V_verts_avg(:,this%N_verts+1:this%N_verts*2), "v_avg")
            call body_vtk%write_point_vectors(this%V_verts_std(:,this%N_verts+1:this%N_verts*2), "v_std_dev")
        end if

        call body_vtk%finish()

    end subroutine surface_mesh_write_body_mirror


    subroutine surface_mesh_write_control_points(this, control_point_file, solved)
        ! Writes the control points (and results if solved) out to file

        implicit none
        
        class(surface_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: control_point_file
        logical,intent(in) :: solved

        type(vtk_out) :: cp_vtk

        ! Clear old file
        call delete_file(control_point_file)

        ! Write out data
        call cp_vtk%begin(control_point_file)
        call cp_vtk%write_points(this%cp)
        call cp_vtk%write_vertices(this%N_cp)
        
        ! Results
        if (solved) then
            call cp_vtk%write_point_scalars(this%phi_cp(1:this%N_cp), "phi")
            call cp_vtk%write_point_scalars(this%phi_cp_mu(1:this%N_cp), "phi_mu")
            call cp_vtk%write_point_scalars(this%phi_cp_sigma(1:this%N_cp), "phi_sigma")
        end if

        call cp_vtk%finish()
    
        
    end subroutine surface_mesh_write_control_points


    subroutine surface_mesh_write_mirrored_control_points(this, control_point_file, solved)
        ! Writes the control points (and results if solved) out to file

        implicit none
        
        class(surface_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: control_point_file
        logical,intent(in) :: solved

        type(vtk_out) :: cp_vtk

        ! Clear old file
        call delete_file(control_point_file)

        ! Write out data
        call cp_vtk%begin(control_point_file)
        call cp_vtk%write_points(this%cp_mirrored)
        call cp_vtk%write_vertices(this%N_cp)
        
        ! Results
        if (solved) then
            call cp_vtk%write_point_scalars(this%phi_cp(this%N_cp+1:this%N_cp*2), "phi")
            call cp_vtk%write_point_scalars(this%phi_cp_mu(this%N_cp+1:this%N_cp*2), "phi_mu")
            call cp_vtk%write_point_scalars(this%phi_cp_sigma(this%N_cp+1:this%N_cp*2), "phi_sigma")
        end if

        call cp_vtk%finish()
    
        
    end subroutine surface_mesh_write_mirrored_control_points

end module surface_mesh_mod