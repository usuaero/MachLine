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

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels, N_cp, N_edges, N_wake_edges
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        type(edge),allocatable,dimension(:) :: edges
        type(wake_mesh) :: wake
        integer,allocatable,dimension(:) :: wake_edge_verts
        real :: C_wake_shedding_angle, trefftz_distance, C_min_panel_angle
        integer :: N_wake_panels_streamwise
        logical :: wake_present, append_wake
        real,dimension(:,:),allocatable :: cp, cp_mirrored
        real,dimension(:),allocatable :: phi_cp, phi_cp_sigma, phi_cp_mu ! Induced potentials at control points
        real,dimension(:),allocatable :: C_p_inc, C_p_ise, C_p_2nd ! Surface pressure coefficients
        real,dimension(:,:),allocatable :: V, dC_f ! Surface velocities and pressure forces
        real :: control_point_offset
        logical :: mirrored ! Whether the mesh is to be mirrored about any planes
        integer :: mirror_plane ! Index of the plane across which the mesh is mirrored (1: yz, 2: xz, 3: xy); this is the index of the normal to that plane
        logical :: asym_flow ! Whether the flow is asymmetric about the mirror plane
        real,dimension(:),allocatable :: mu, sigma ! Singularity strengths
        real :: S_ref ! Reference parameters

        contains

            procedure :: init => surface_mesh_init
            procedure :: init_with_flow => surface_mesh_init_with_flow
            procedure :: output_results => surface_mesh_output_results
            procedure :: locate_adjacent_panels => surface_mesh_locate_adjacent_panels
            procedure :: characterize_edges => surface_mesh_characterize_edges
            procedure :: find_vertices_on_mirror => surface_mesh_find_vertices_on_mirror
            procedure :: clone_vertices => surface_mesh_clone_vertices
            procedure :: set_up_mirroring => surface_mesh_set_up_mirroring
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
        write(*,'(a, i1, a, i1, a)') "     User has selected ", doublet_order, &
                                     "-order doublet panels and ", source_order, "-order source panels."

        ! Check
        if (doublet_order /= 1 .or. source_order /= 0) then
            write(*,*) "!!! Such distributions are not currently available."
            write(*,*) "!!! Defaulting to a linear doublet distribution and a constant source distribution."
            doublet_order = 1
            source_order = 0
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
        write(*,'(a, a, a)',advance='no') "     Reading surface mesh in from file ", mesh_file, "..."
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
        write(*,*) "Done."

        ! Display mesh info
        write(*,'(a, i7, a, i7, a)') "     Surface mesh has ", this%N_verts, " vertices and ", this%N_panels, " panels."

        ! Get mirroring
        call json_xtnsn_get(settings, 'mirror_about', mirror_plane, "none")
        this%mirrored = .true.
        select case (mirror_plane)
        case ("xy")
            this%mirror_plane = 3
            write(*,*) "    Mesh set to mirror about xy plane."
        case ("xz")
            this%mirror_plane = 2
            write(*,*) "    Mesh set to mirror about xz plane."
        case ("yz")
            this%mirror_plane = 1
            write(*,*) "    Mesh set to mirror about yz plane."
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

        class(surface_mesh),intent(inout) :: this

        integer :: i, j, m, n, m1, n1, temp, i_edge, i_panel1, i_panel2, i_vert1, i_vert2, edge_on_mirror, i_edge1, i_edge2, N_edges
        logical :: already_found_shared, dummy
        real :: distance
        integer,dimension(2) :: shared_verts
        integer,dimension(this%N_panels*4) :: panel1, panel2, vertex1, vertex2, edge_index1, edge_index2
        logical,dimension(this%N_panels*4) :: on_mirror_plane

        write(*,'(a)',advance='no') "     Locating adjacent panels..."

        ! Loop through each panel
        N_edges = 0
        !$OMP parallel private(j, already_found_shared, distance, shared_verts, m, m1, n, n1, temp, i_edge) &
        !$OMP private(i_panel1, i_panel2, i_vert1, i_vert2, edge_on_mirror, i_edge1, i_edge2)

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

                        ! Get distance between vertices. This is more robust than checking vertex indices; mesh may not be ideal.
                        ! TODO : refine mesh import to ensure we can do this purely off of indices
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

                                !$OMP critical
                                
                                ! Update number of edges
                                N_edges = N_edges + 1
                                i_edge = N_edges

                                ! Store vertices being adjacent to one another
                                if (.not. this%vertices(shared_verts(1))%adjacent_vertices%is_in(shared_verts(2))) then
                                    call this%vertices(shared_verts(1))%adjacent_vertices%append(shared_verts(2))
                                end if
                                if (.not. this%vertices(shared_verts(2))%adjacent_vertices%is_in(shared_verts(1))) then
                                    call this%vertices(shared_verts(2))%adjacent_vertices%append(shared_verts(1))
                                end if

                                ! Store that this edge touches the two end vertices
                                call this%vertices(shared_verts(1))%adjacent_edges%append(i_edge)
                                call this%vertices(shared_verts(2))%adjacent_edges%append(i_edge)

                                ! Store adjacent panels and panel edges
                                ! This stores the adjacent panels and edges according to the index of that edge
                                ! for the current panel

                                ! Store that i is adjacent to j
                                ! This one is more complicated because we don't know that n1 will be less than n; just the nature of the nested loop.
                                ! Basically, if one is 1 and the other is N, then we're dealing with edge N for panel j.
                                ! Otherwise, we're dealing with abs(n1-n) being 1, meaning edge min(n1, n).
                                if ( (n1 == 1 .and. n == this%panels(j)%N) .or. (n == 1 .and. n1 == this%panels(j)%N) ) then
                                    this%panels(j)%abutting_panels(this%panels(j)%N) = i
                                    this%panels(j)%edges(this%panels(j)%N) = i_edge
                                    edge_index2(i_edge) = this%panels(j)%N
                                else
                                    n1 = min(n, n1)
                                    this%panels(j)%abutting_panels(n1) = i
                                    this%panels(j)%edges(n1) = i_edge
                                    edge_index2(i_edge) = n1
                                end if

                                ! Store that j is adjacent to i
                                if (m1 == 1 .and. m == this%panels(i)%N) then ! Nth edge
                                    this%panels(i)%abutting_panels(m) = j
                                    this%panels(i)%edges(m) = i_edge
                                    edge_index1(i_edge) = m
                                else ! 1st or 2nd edge
                                    this%panels(i)%abutting_panels(m1) = j
                                    this%panels(i)%edges(m1) = i_edge
                                    edge_index1(i_edge) = m1
                                end if

                                ! Store information in arrays for later storage in edge objects
                                panel1(i_edge) = i
                                panel2(i_edge) = j
                                vertex1(i_edge) = shared_verts(1)
                                vertex2(i_edge) = shared_verts(2)
                                on_mirror_plane(i_edge) = .false.

                                !$OMP end critical
                                
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

            ! For each panel, check if it abuts the mirror plane
            if (this%mirrored) then

                ! Initialize checks
                already_found_shared = .false.

                ! Loop through vertices
                mirror_loop: do m=1,this%panels(i)%N

                    ! Check if vertex is on the mirror plane
                    n = this%panels(i)%vertex_indices(m)
                    if (this%vertices(n)%on_mirror_plane) then

                        ! Previously found a vertex on mirror plane, so the panels are abutting
                        if (already_found_shared) then

                            ! Store the second shared vertex
                            shared_verts(2) = n

                            ! Check order
                            if (m1 == 1 .and. m == 3) then
                                temp = shared_verts(1)
                                shared_verts(1) = shared_verts(2)
                                shared_verts(2) = temp
                            end if

                            !$OMP critical
                            
                            ! Update number of edges
                            N_edges = N_edges + 1
                            i_edge = N_edges

                            ! Store adjacent vertices
                            if (.not. this%vertices(shared_verts(1))%adjacent_vertices%is_in(shared_verts(2))) then
                                call this%vertices(shared_verts(1))%adjacent_vertices%append(shared_verts(2))
                            end if
                            if (.not. this%vertices(shared_verts(2))%adjacent_vertices%is_in(shared_verts(1))) then
                                call this%vertices(shared_verts(2))%adjacent_vertices%append(shared_verts(1))
                            end if

                            ! Store adjacent panel
                            if (m-m1 == 1) then
                                this%panels(i)%abutting_panels(m1) = i+this%N_panels
                                this%panels(i)%edges(m1) = i_edge
                                edge_index1(i_edge) = m1
                            else
                                this%panels(i)%abutting_panels(m) = i+this%N_panels
                                this%panels(i)%edges(m) = i_edge
                                edge_index1(i_edge) = m
                            end if

                            ! Store in arrays for later storage in edge objects
                            panel1(i_edge) = i
                            panel2(i_edge) = i+this%N_panels
                            vertex1(i_edge) = shared_verts(1)
                            vertex2(i_edge) = shared_verts(2)
                            on_mirror_plane(i_edge) = .true.
                            edge_index2(i_edge) = 0 ! Just a placeholder since the second panel doesn't technically exist

                            !$OMP end critical

                            ! Break out of loop
                            exit mirror_loop

                        ! First vertex on the mirror plane
                        else

                            already_found_shared = .true.
                            shared_verts(1) = n
                            m1 = m

                        end if
                    end if

                end do mirror_loop

            end if

        end do

        ! Check that no panel abuts empty space (i.e. non-watertight mesh)
        !$OMP do schedule(static)
        do i=1,this%N_panels
            if (any(this%panels(i)%abutting_panels == 0)) then
                write(*,*)
                write(*,*) "!!! The supplied mesh is not watertight. Panel", i, "is missing at least one neighbor. Quitting..."
                write(*,*) this%panels(i)%abutting_panels
                stop
            end if
        end do

        ! Allocate edge storage
        !$OMP single
        this%N_edges = N_edges
        allocate(this%edges(this%N_edges))
        !$OMP end single

        ! Initialize edges
        !$OMP do schedule(static)
        do i=1,this%N_edges

            ! Initialize
            call this%edges(i)%init(vertex1(i), vertex2(i), panel1(i), panel2(i))

            ! Store more information
            this%edges(i)%on_mirror_plane = on_mirror_plane(i)
            this%edges(i)%edge_index_for_panel(1) = edge_index1(i)
            this%edges(i)%edge_index_for_panel(2) = edge_index2(i)

        end do

        !$OMP end parallel

        write(*,"(a, i7, a)") "Done. Found ", this%N_edges, " edges."
    
    end subroutine surface_mesh_locate_adjacent_panels


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
            call this%panels(i)%calc_transforms(freestream)
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

        ! Clone necessary vertices and calculate normals (for placing control points)
        if (this%wake_present .or. freestream%supersonic) then
            call this%clone_vertices()
        end if
        call this%calc_vertex_normals()

        ! INitialize wake
        call this%init_wake(freestream, wake_file)
    
    end subroutine surface_mesh_init_with_flow


    subroutine surface_mesh_characterize_edges(this, freestream)
        ! Locates wake-shedding edges and supersonic/subsonic edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        integer :: i, j, k, m, n, temp, top_panel, bottom_panel, i_vert_1, i_vert_2, N_wake_edge_verts
        integer,dimension(:),allocatable :: wake_edge_verts
        real :: C_angle, C_min_angle
        real,dimension(3) :: second_normal

        write(*,'(a)',advance='no') "     Characterizing edges..."

        ! Initialize
        allocate(wake_edge_verts(this%N_verts/4)) ! I sure hope you're not trying to run a mesh where every fourth vertex has a wake emanating from it...
        this%N_wake_edges = 0
        N_wake_edge_verts = 0
        C_min_angle = 100.

        ! Loop through each edge
        !$OMP parallel do private(i, j, second_normal, C_angle, i_vert_1, i_vert_2) reduction(min : C_min_angle) &
        !$OMP & default(none) shared(this, freestream, wake_edge_verts, N_wake_edge_verts)
        do k=1,this%N_edges

            ! Get info
            i = this%edges(k)%panels(1)
            j = this%edges(k)%panels(2)

            ! Get normal for panel j (dependent on mirroring)
            if (this%edges(k)%on_mirror_plane) then
                second_normal = mirror_about_plane(this%panels(i)%n_g, this%mirror_plane)
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

                    ! Set the character of the edge
                    this%edges(k)%sheds_wake = .true.
                    this%edges(k)%discontinuous = .true.

                    !$OMP critical

                    ! Update number of wake-shedding edges
                    this%N_wake_edges = this%N_wake_edges + 1

                    ! If this vertex does not already belong to a wake-shedding edge, add it to the list of wake edge vertices
                    if (this%vertices(i_vert_1)%N_wake_edges == 0) then 
                        N_wake_edge_verts = N_wake_edge_verts + 1
                        wake_edge_verts(N_wake_edge_verts) = i_vert_1
                        this%vertices(i_vert_1)%index_in_wake_vertices = N_wake_edge_verts

                    ! If it does already belong to a wake-shedding edge, then we may now conclude it is 'in' an edge
                    ! Because of this, it will likely need to be cloned unless it's on a mirror plane
                    else if (.not. this%edges(k)%on_mirror_plane) then
                        this%vertices(i_vert_1)%needs_clone = .true.
                    end if

                    ! Update number of wake edges touching this vertex
                    this%vertices(i_vert_1)%N_wake_edges = this%vertices(i_vert_1)%N_wake_edges + 1
                    this%vertices(i_vert_1)%N_discont_edges = this%vertices(i_vert_1)%N_discont_edges + 1

                    ! Do the same for the other vertex
                    if (this%vertices(i_vert_2)%N_wake_edges == 0) then 
                        N_wake_edge_verts = N_wake_edge_verts + 1
                        wake_edge_verts(N_wake_edge_verts) = i_vert_2
                        this%vertices(i_vert_2)%index_in_wake_vertices = N_wake_edge_verts

                    ! If it does already belong to a wake-shedding edge, then we may now conclude it is 'in' an edge
                    ! Because of this, it will likely need to be cloned unless it's on a mirror plane
                    else if (.not. this%edges(k)%on_mirror_plane) then
                        this%vertices(i_vert_2)%needs_clone = .true.
                    end if

                    ! Update number of wake edges touching this vertex
                    this%vertices(i_vert_2)%N_wake_edges = this%vertices(i_vert_2)%N_wake_edges + 1
                    this%vertices(i_vert_2)%N_discont_edges = this%vertices(i_vert_2)%N_discont_edges + 1

                    !$OMP end critical

                end if
            end if

        end do

        ! Store minimum angle
        this%C_min_panel_angle = C_min_angle

        ! Allocate wake vertices array
        allocate(this%wake_edge_verts, source=wake_edge_verts(1:N_wake_edge_verts))

        write(*,'(a, i3, a, i3, a)') "Done. Found ", this%N_wake_edges, " wake-shedding edges."

    end subroutine surface_mesh_characterize_edges


    subroutine surface_mesh_set_up_mirroring(this)
        ! Sets up information necessary for mirroring the mesh

        implicit none

        class(surface_mesh),intent(inout) :: this

        integer :: i, j, m, n

        write(*,'(a)',advance='no') "     Setting up mesh mirror..."

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

                ! Get endpoint indices
                m = this%edges(i)%verts(1)
                n = this%edges(i)%verts(2)

                ! If a given discontinuous edge has only one of its endpoints lying on the mirror plane, then that endpoint has another
                ! adjacent edge, since the edge will be mirrored across that plane. This endpoint will need a clone, but it's mirrored
                ! vertex will be the same

                ! Check for edge touching mirror plane at only one end
                if (this%vertices(m)%on_mirror_plane .neqv. this%vertices(n)%on_mirror_plane) then

                    ! Only add discontinuous edge if the endpoint is on the mirror plane
                    if (this%vertices(m)%on_mirror_plane) then

                        this%vertices(m)%N_wake_edges = this%vertices(m)%N_wake_edges + 1
                        this%vertices(m)%N_discont_edges = this%vertices(m)%N_discont_edges + 1
                        this%vertices(m)%needs_clone = .true.
                        this%vertices(m)%mirrored_is_unique = .false.

                    else

                        this%vertices(n)%N_wake_edges = this%vertices(n)%N_wake_edges + 1
                        this%vertices(n)%N_discont_edges = this%vertices(n)%N_discont_edges + 1
                        this%vertices(n)%needs_clone = .true.
                        this%vertices(n)%mirrored_is_unique = .false.

                    end if

                ! If the discontinuous edge has both endpoints lying on the mirror plane, then these vertices need no clones, but the mirrored
                ! vertices will still be unique
                else if (this%edges(i)%on_mirror_plane) then

                    ! The mirrored vertices will function as clones
                    this%vertices(m)%needs_clone = .false.
                    this%vertices(n)%needs_clone = .false.
                    this%vertices(m)%mirrored_is_unique = .true.
                    this%vertices(n)%mirrored_is_unique = .true.

                end if
            end if
        end do
        write(*,*) "Done."

    end subroutine surface_mesh_set_up_mirroring


    subroutine surface_mesh_clone_vertices(this)
        ! Takes vertices which lie within discontinuous edges and splits them into two vertices.
        ! Handles rearranging of necessary dependencies.

        implicit none

        class(surface_mesh),intent(inout),target :: this

        integer :: i, j, k, m, n, N_clones, i_jango, i_boba, N_discont_verts, i_bot_panel, i_abutting_panel, i_adj_vert, i_top_panel
        integer :: i_edge
        type(vertex),dimension(:),allocatable :: temp_vertices

        write(*,'(a)',advance='no') "     Cloning vertices at discontinuous edges..."

        ! Allocate array which will store which discontinuous vertices need to be cloned
        N_discont_verts = size(this%wake_edge_verts)

        ! Determine number of vertices which need to be cloned
        N_clones = 0
        do i=1,N_discont_verts

            if (this%vertices(this%wake_edge_verts(i))%needs_clone) then

                ! Update the number of needed clones
                N_clones = N_clones + 1

            end if
        end do

        ! Extend allocation of mesh vertex array
        allocate(temp_vertices, source=this%vertices)
        deallocate(this%vertices)
        allocate(this%vertices(this%N_verts + N_clones))

        ! Place existing vertices in new array
        this%vertices(1:this%N_verts) = temp_vertices
        deallocate(temp_vertices)

        ! Update number of vertices
        this%N_verts = this%N_verts + N_clones

        ! Fix vertex pointers in panel objects (necessary because this%vertices got reallocated)
        do i=1,this%N_panels
            do j=1,this%panels(i)%N
                this%panels(i)%vertices(j)%ptr => this%vertices(this%panels(i)%vertex_indices(j))
            end do
        end do

        ! Initialize clones
        j = 1
        do i=1,N_discont_verts

            ! Get index of vertex to be cloned
            i_jango = this%wake_edge_verts(i)

            ! Check if this vertex needs to be cloned
            if (this%vertices(i_jango)%needs_clone) then

                ! Get index for the clone
                i_boba = this%N_verts - N_clones + j ! Will be at position N_verts-N_clones+j in the new vertex array

                ! Initialize clone
                call this%vertices(i_boba)%init(this%vertices(i_jango)%loc, i_boba)

                ! Specify wake partners
                this%vertices(i_jango)%i_wake_partner = i_boba
                this%vertices(i_boba)%i_wake_partner = i_jango

                ! Store number of adjacent wake-shedding edges (probably unecessary at this point, but let's be consistent)
                this%vertices(i_boba)%N_wake_edges = this%vertices(i_jango)%N_wake_edges

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
                                            if (.not. this%vertices(i_boba)%panels_not_across_wake_edge%is_in(i_abutting_panel)) &
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
                        call this%panels(i_bot_panel)%point_to_vertex_clone(this%vertices(i_boba))
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

        write(*,'(a, i4, a, i7, a)') "Done. Cloned ", N_clones, " vertices. Mesh now has ", this%N_verts, " vertices."

    end subroutine surface_mesh_clone_vertices


    subroutine surface_mesh_calc_vertex_normals(this)
        ! Initializes the normal vectors associated with each vertex.
        ! Must be called only once wake-shedding edge vertices have been cloned.

        implicit none

        class(surface_mesh),intent(inout) :: this

        real,dimension(3) :: vec_sum
        integer :: i, j, i_panel

        write(*,'(a)',advance='no') "     Calculating vertex normals..."

        ! Loop through vertices
        !$OMP parallel do private(vec_sum, i, i_panel) schedule(dynamic)
        do j=1,this%N_verts

            ! Loop through neighboring panels and compute the average of their normal vectors
            vec_sum = 0
            do i=1,this%vertices(j)%panels%len()
                call this%vertices(j)%panels%get(i, i_panel)
                vec_sum = vec_sum + this%panels(i_panel)%n_g
            end do

            ! For vertices on the mirror plane, the component normal to the plane should be zeroed
            if (this%mirrored .and. this%vertices(j)%on_mirror_plane) then
                vec_sum(this%mirror_plane) = 0.
            end if

            ! Normalize and store
            this%vertices(j)%n_g = vec_sum/norm2(vec_sum)

            ! Calculate average edge lengths for each vertex
            call this%vertices(j)%calc_average_edge_length(this%vertices)

        end do

        write(*,*) "Done."

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
            call this%wake%init(freestream, this%wake_edge_verts, &
                                this%edges, this%N_wake_edges, &
                                this%N_wake_panels_streamwise, this%vertices, &
                                this%trefftz_distance, this%mirrored .and. this%asym_flow, &
                                this%mirror_plane, this%N_panels)

            ! Clean up
            deallocate(this%wake_edge_verts)
        
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
        real,dimension(3) :: normal
        real :: offset_ratio

        if (doublet_order == 1) then

            ! Specify number of control points
            this%N_cp = this%N_verts

            ! Allocate memory
            allocate(this%cp(3,this%N_verts))

            ! Calculate offset ratio such that the control point will remain within the body based on the minimum detected angle between panels
            if (this%wake_present) then
                offset_ratio = 0.5*sqrt(0.5*(1. + this%C_min_panel_angle))
            end if

            ! Loop through vertices
            !$OMP parallel do private(j, normal, i_panel) schedule(dynamic) shared(this, offset, offset_ratio) default(none)
            do i=1,this%N_verts

                ! If the vertex is in a wake edge, it needs to be shifted off the normal slightly so that it is unique from its counterpart
                if (this%vertices(i)%N_wake_edges > 1) then

                    ! Loop through panels associated with this clone to get their average normal vector
                    normal = 0.
                    do j=1,this%vertices(i)%panels_not_across_wake_edge%len()

                        ! Get panel index
                        call this%vertices(i)%panels_not_across_wake_edge%get(j, i_panel)

                        ! Add normal vector
                        normal = normal + this%panels(i_panel)%n_g

                    end do

                    ! Add effect of mirrored panels
                    if (this%vertices(i)%on_mirror_plane) then
                        normal(this%mirror_plane) = 0.
                    end if

                    ! Normalize
                    normal = normal/norm2(normal)

                    ! Place control point
                    this%cp(:,i) = this%vertices(i)%loc &
                                               - offset * (this%vertices(i)%n_g - offset_ratio * normal)*this%vertices(i)%l_avg

                ! If it's not in a wake-shedding edge (i.e. has no clone), then placement simply follows the normal vector
                else

                    this%cp(:,i) = this%vertices(i)%loc - offset*this%vertices(i)%n_g*this%vertices(i)%l_avg

                end if

            end do

        end if

        ! Calculate mirrored points, if necessary
        if (this%mirrored) then

            ! Allocate memory
            allocate(this%cp_mirrored(3,this%N_cp))

            ! Calculate mirrors
            !$OMP parallel do schedule(static)
            do i=1,this%N_cp
                this%cp_mirrored(:,i) = mirror_about_plane(this%cp(:,i), this%mirror_plane)
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

            ! Write data
            call body_vtk%begin(body_file)
            call body_vtk%write_points(this%vertices)
            call body_vtk%write_panels(this%panels)
            call body_vtk%write_cell_scalars(this%sigma(1:this%N_panels), "sigma")
            if (allocated(this%C_p_inc)) then
                call body_vtk%write_cell_scalars(this%C_p_inc(1:this%N_panels), "C_p_inc")
            end if
            if (allocated(this%C_p_ise)) then
                call body_vtk%write_cell_scalars(this%C_p_ise(1:this%N_panels), "C_p_ise")
            end if
            if (allocated(this%C_p_2nd)) then
                call body_vtk%write_cell_scalars(this%C_p_2nd(1:this%N_panels), "C_p_2nd")
            end if
            call body_vtk%write_cell_scalars(panel_inclinations, "inclination")
            call body_vtk%write_cell_vectors(this%v(:,1:this%N_panels), "v")
            call body_vtk%write_cell_vectors(this%dC_f(:,1:this%N_panels), "dC_f")
            call body_vtk%write_point_scalars(this%mu(1:this%N_cp), "mu")
            call body_vtk%finish()

            write(*,*) "    Surface results written to: ", body_file
        end if

        ! Write out data for mirrored body
        if (mirrored_body_file /= 'none' .and. this%asym_flow) then

            ! Clear old file
            call delete_file(mirrored_body_file)

            ! Write geometry
            call body_vtk%begin(mirrored_body_file)
            call body_vtk%write_points(this%vertices, this%mirror_plane)
            call body_vtk%write_panels(this%panels)

            ! Write source strengths
            call body_vtk%write_cell_scalars(this%sigma(this%N_panels+1:this%N_panels*2), "sigma")

            ! Write pressures
            if (allocated(this%C_p_inc)) then
                call body_vtk%write_cell_scalars(this%C_p_inc(this%N_panels+1:this%N_panels*2), "C_p_inc")
            end if
            if (allocated(this%C_p_ise)) then
                call body_vtk%write_cell_scalars(this%C_p_ise(this%N_panels+1:this%N_panels*2), "C_p_ise")
            end if
            if (allocated(this%C_p_2nd)) then
                call body_vtk%write_cell_scalars(this%C_p_2nd(this%N_panels+1:this%N_panels*2), "C_p_2nd")
            end if

            ! Write flow properties
            call body_vtk%write_cell_vectors(this%v(:,this%N_panels+1:this%N_panels*2), "v")
            call body_vtk%write_cell_vectors(this%dC_f(:,this%N_panels+1:this%N_panels*2), "dC_f")
            call body_vtk%write_point_scalars(this%mu(this%N_cp+1:this%N_cp*2), "mu")
            call body_vtk%finish()

            write(*,*) "    Mirrored surface results written to: ", mirrored_body_file

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
                write(*,*) "    Wake results written to: ", wake_file

            else
                write(*,*) "    No wake to export."

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

            write(*,*) "    Control point results written to: ", control_point_file
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

            write(*,*) "    Mirrored control point results written to: ", mirrored_control_point_file
        end if
    
    end subroutine surface_mesh_output_results


end module surface_mesh_mod