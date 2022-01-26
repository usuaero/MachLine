! A surface mesh type encompassing a body, wakes, and shocks
module surface_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use vtk_mod
    use vertex_mod
    use panel_mod
    use flow_mod
    use math_mod
    use edge_mod
    use wake_mesh_mod

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels, N_cp, N_edges, N_wake_edges
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        type(edge),allocatable,dimension(:) :: edges
        type(wake_mesh) :: wake
        integer,allocatable,dimension(:) :: wake_edge_top_verts, wake_edge_bot_verts, wake_edge_indices
        real :: wake_shedding_angle, C_wake_shedding_angle, trefftz_distance, C_min_panel_angle
        integer :: N_wake_panels_streamwise
        logical :: wake_present, append_wake
        real,dimension(:,:),allocatable :: control_points, cp_mirrored
        real,dimension(:),allocatable :: phi_cp, phi_cp_sigma, phi_cp_mu ! Induced potentials at control points
        real,dimension(:),allocatable :: C_p_inc, C_p_ise, C_p_2nd ! Surface pressure coefficients
        real,dimension(:,:),allocatable :: V ! Surface velocities
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
            procedure :: update_supersonic_trefftz_distance => surface_mesh_update_supersonic_trefftz_distance
            procedure :: place_interior_control_points => surface_mesh_place_interior_control_points

    end type surface_mesh
    

contains


    subroutine surface_mesh_init(this, settings)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(json_value),pointer,intent(inout) :: settings
        character(len=:),allocatable :: extension
        integer :: loc
        character(len=:),allocatable :: mesh_file, mirror_plane

        ! Set singularity orders
        call json_xtnsn_get(settings, 'singularity_order.doublet', doublet_order, 1)
        call json_xtnsn_get(settings, 'singularity_order.source', source_order, 0)
        write(*,'(a, i1, a, i1, a)') "     User has selected: ", doublet_order, &
                                     "-order doublet panels and ", source_order, "-order source panels."

        ! Check
        if (doublet_order /= 1 .or. source_order /= 0) then
            write(*,*) "    !!! Such distributions are not currently available."
            write(*,*) "    !!! Defaulting a linear doublet distribution and a constant source distribution."
            doublet_order = 1
            source_order = 0
        end if

        ! Get mesh file
        call json_get(settings, 'file', mesh_file)
        mesh_file = trim(mesh_file)
        write(*,*) "    Reading surface mesh in from file: ", mesh_file

        ! Determine the type of mesh file
        loc = index(mesh_file, '.')
        extension = mesh_file(loc:len(mesh_file))

        ! Load vtk
        if (extension == '.vtk') then
            call load_surface_vtk(mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        else
            write(*,*) "TriPan cannot read ", extension, " type mesh files. Quitting..."
            stop
        end if

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
            call json_xtnsn_get(settings, 'wake_model.wake_shedding_angle', this%wake_shedding_angle, 90.0) ! Maximum allowable angle between panel normals without having separation
            this%C_wake_shedding_angle = cos(this%wake_shedding_angle*pi/180.0)

            if (this%append_wake) then
                call json_xtnsn_get(settings, 'wake_model.trefftz_distance', this%trefftz_distance, 100.0) ! Distance from origin to wake termination
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


        integer :: i, j, m, n, m1, n1
        logical :: already_found_shared, dummy
        real :: distance
        integer,dimension(2) :: shared_verts
        type(list) :: panel1, panel2, vertex1, vertex2, on_mirror_plane, edge_index1, edge_index2

        write(*,'(a)',advance='no') "     Locating adjacent panels..."

        ! Loop through each panel
        do i=1,this%N_panels

            ! Loop through each potential neighbor
            do j=i+1,this%N_panels

                ! Initialize for this panel pair
                already_found_shared = .false.

                ! Check if the panels are abutting
                abutting_loop: do m=1,this%panels(i)%N
                    do n=1,this%panels(j)%N

                        ! Get distance between vertices. This is more robust than checking vertex indices; mesh may not be ideal.
                        distance = dist(this%panels(i)%get_vertex_loc(m), this%panels(j)%get_vertex_loc(n))

                        ! Check distance
                        if (distance < 1e-12) then

                            ! Previously found a shared vertex
                            if (already_found_shared) then

                                ! Store second shared vertex
                                shared_verts(2) = this%panels(i)%get_vertex_index(m)

                                ! Store in lists for later storage in edge objects
                                call panel1%append(i)
                                call panel2%append(j)
                                call vertex1%append(shared_verts(1))
                                call vertex2%append(shared_verts(2))
                                call on_mirror_plane%append(.false.)

                                ! Store adjacent vertices
                                if (.not. this%vertices(shared_verts(1))%adjacent_vertices%is_in(shared_verts(2))) then
                                    call this%vertices(shared_verts(1))%adjacent_vertices%append(shared_verts(2))
                                end if
                                if (.not. this%vertices(shared_verts(2))%adjacent_vertices%is_in(shared_verts(1))) then
                                    call this%vertices(shared_verts(2))%adjacent_vertices%append(shared_verts(1))
                                end if

                                ! Store adjacent panels and panel edges
                                ! This stores the adjacent panels and edges according to the index of that edge
                                ! for the current panel
                                if (m-m1 == 1) then
                                    this%panels(i)%abutting_panels(m1) = j
                                    this%panels(i)%edges(m1) = panel1%len()
                                    call edge_index1%append(m1)
                                else
                                    this%panels(i)%abutting_panels(m) = j
                                    this%panels(i)%edges(m) = panel1%len()
                                    call edge_index1%append(m)
                                end if
                                if (n-n1 == 1) then
                                    this%panels(j)%abutting_panels(n1) = i
                                    this%panels(i)%edges(n1) = panel1%len()
                                    call edge_index2%append(n1)
                                else
                                    this%panels(j)%abutting_panels(n) = i
                                    this%panels(i)%edges(n) = panel1%len()
                                    call edge_index2%append(n)
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

            end do

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

                            ! Store in lists for later storage in arrays
                            call panel1%append(i)
                            call panel2%append(i+this%N_panels)
                            call vertex1%append(shared_verts(1))
                            call vertex2%append(shared_verts(2))
                            call on_mirror_plane%append(.true.)

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
                                call edge_index1%append(m1)
                            else
                                this%panels(i)%abutting_panels(m) = i+this%N_panels
                                call edge_index1%append(m)
                            end if
                            call edge_index2%append(0) ! Really meaningless since the second panel doesn't technically exist

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

        ! Allocate edge storage
        this%N_edges = panel1%len()
        allocate(this%edges(this%N_edges))

        ! Initialize edges
        do i=1,this%N_edges

            ! Get information
            call vertex1%get(i, shared_verts(1))
            call vertex2%get(i, shared_verts(2))
            call panel1%get(i, m)
            call panel2%get(i, n)
            call on_mirror_plane%get(i, dummy)

            ! Initialize
            call this%edges(i)%init(shared_verts(1), shared_verts(2), m, n)

            ! Store more information
            call edge_index1%get(i, m)
            call edge_index2%get(i, n)
            this%edges(i)%on_mirror_plane = dummy
            this%edges(i)%edge_index_for_panel(1) = m
            this%edges(i)%edge_index_for_panel(2) = n

        end do

        write(*,"(a, i7, a)") "Done. Found ", this%N_edges, " edges."
    
    end subroutine surface_mesh_locate_adjacent_panels


    subroutine surface_mesh_init_with_flow(this, freestream, wake_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream
        character(len=:),allocatable,intent(in) :: wake_file

        integer :: i
        type(vtk_out) :: wake_vtk

        ! Check flow symmetry condition
        this%asym_flow = .false.
        if (this%mirrored) then
            if (.not. freestream%sym_about(this%mirror_plane)) then
                this%asym_flow = .true.
            end if
        end if

        ! Calculate panel coordinate transformations
        do i=1,this%N_panels
            call this%panels(i)%calc_transforms(freestream)
        end do

        ! Figure out wake-shedding edges, discontinuous edges, etc.
        ! Edge-characterization is only necessary for flows with wakes or supersonic flows, as discontinuities only appear in these flows
        if (this%wake_present .or. freestream%supersonic) then
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

        ! Handle wake creation
        if (this%append_wake) then

            ! For supersonic flows, calculate the Trefftz distance based on the mesh
            if (freestream%supersonic) then
                call this%update_supersonic_trefftz_distance(freestream)
            end if

            ! Initialize wake
            call this%wake%init(freestream, this%wake_edge_top_verts, &
                                this%wake_edge_bot_verts, this%edges, this%wake_edge_indices, &
                                this%N_wake_panels_streamwise, this%vertices, &
                                this%trefftz_distance, this%mirrored .and. this%asym_flow, &
                                this%mirror_plane)

            ! Clean up
            deallocate(this%wake_edge_top_verts)
            deallocate(this%wake_edge_bot_verts)
        
            ! Export wake geometry
            if (wake_file /= 'none') then

                ! Clear old file
                call delete_file(wake_file)

                ! Write new geometry
                if (this%wake%N_panels > 0) then

                    call wake_vtk%begin(wake_file)
                    call wake_vtk%write_points(this%wake%vertices)
                    call wake_vtk%write_panels(this%wake%panels)
                    call wake_vtk%finish()

                end if
            end if

        else
            
            ! Set parameters to let later code know there is no actual wake
            this%wake%N_panels = 0
            this%wake%N_verts = 0

        end if
    
    end subroutine surface_mesh_init_with_flow


    subroutine surface_mesh_characterize_edges(this, freestream)
        ! Locates wake-shedding edges and supersonic/subsonic edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        integer :: i, j, k, m, n, temp, top_panel, bottom_panel, vert1, vert2
        type(list) :: wake_edge_verts, wake_edges
        real :: C_angle
        real,dimension(3) :: second_normal

        write(*,'(a)',advance='no') "     Characterizing edges..."

        ! We need to store the minimum angle between two panels in order to place control points within the body at edges having discontinuities
        this%C_min_panel_angle = 1.

        ! Loop through each edge
        this%N_wake_edges = 0
        do k=1,this%N_edges

            ! Get info
            i = this%edges(k)%panels(1)
            j = this%edges(k)%panels(2)

            ! Get normal for panel j (dependent on mirroring)
            if (this%edges(k)%on_mirror_plane) then
                second_normal = mirror_about_plane(this%panels(i)%normal, this%mirror_plane)
            else
                second_normal = this%panels(j)%normal
            end if

            ! Calculate angle between panels (this is the flow-turning angle; it is the most straightforward to compute)
            C_angle = inner(this%panels(i)%normal, second_normal)

            ! Update minimum angle
            this%C_min_panel_angle = min(C_angle, this%C_min_panel_angle)

            ! Determine if this edge is wake-shedding; this depends on the angle between the panels
            ! and the angles made by the panel normals with the freestream
            if (this%wake_present) then

                ! Check angle between panels
                if (C_angle < this%C_wake_shedding_angle) then

                    ! Check angle of panel normal with freestream
                    if (inner(this%panels(i)%normal, freestream%V_inf) > 0.0 .or. &
                        inner(second_normal, freestream%V_inf) > 0.0) then

                        ! Set the character of the edge
                        this%edges(k)%sheds_wake = .true.
                        this%edges(k)%discontinuous = .true.

                        ! Update number of wake-shedding edges
                        this%N_wake_edges = this%N_wake_edges + 1

                        ! Store the index of this edge as being wake-shedding
                        call wake_edges%append(k)

                        ! Get vertex indices (simplifies later code)
                        vert1 = this%edges(k)%verts(1)
                        vert2 = this%edges(k)%verts(2)

                        ! If this vertex does not already belong to a wake-shedding edge, add it to the list of wake edge vertices
                        if (this%vertices(vert1)%N_wake_edges == 0) then 

                            ! Add the first time
                            call wake_edge_verts%append(vert1)
                            this%vertices(vert1)%index_in_wake_vertices = wake_edge_verts%len()

                        else if (.not. this%edges(k)%on_mirror_plane) then

                            ! It is in an edge, so it will likely need to be cloned
                            ! Unless it's on a mirror plane
                            this%vertices(vert1)%needs_clone = .true.

                        end if

                        ! Update number of wake edges touching this vertex
                        this%vertices(vert1)%N_wake_edges = this%vertices(vert1)%N_wake_edges + 1
                        this%vertices(vert1)%N_discont_edges = this%vertices(vert1)%N_discont_edges + 1

                        ! Do the same for the other vertex
                        if (this%vertices(vert2)%N_wake_edges == 0) then 

                            ! Add the first time
                            call wake_edge_verts%append(vert2)
                            this%vertices(vert2)%index_in_wake_vertices = wake_edge_verts%len()

                        else if (.not. this%edges(k)%on_mirror_plane) then

                            ! It is in an edge, so it will likely need to be cloned
                            ! Unless it's on a mirror plane
                            this%vertices(vert2)%needs_clone = .true.

                        end if

                        ! Update number of wake edges touching this vertex
                        this%vertices(vert2)%N_wake_edges = this%vertices(vert2)%N_wake_edges + 1
                        this%vertices(vert2)%N_discont_edges = this%vertices(vert2)%N_discont_edges + 1

                    end if
                end if
            end if

            ! Determine if this edge is discontinuous in supersonic flow
            if (freestream%supersonic) then

                ! Update edge inclination
                this%edges(k)%inclination = this%panels(this%edges(i)%panels(1))%q(this%edges(i)%edge_index_for_panel(1))

                ! According to Davis, sharp, subsonic, leading edges in supersonic flow must have discontinuous doublet strength.
                ! I don't know why this would be, except in the case of leading-edge vortex separation. But Davis doesn't
                ! model leading-edge vortices. Wake-shedding trailing edges are still discontinuous in supersonic flow. Supersonic
                ! leading edges should have continuous doublet strength.

            end if

        end do

        ! Allocate wake top vertices array
        allocate(this%wake_edge_top_verts(wake_edge_verts%len()))

        do i=1,wake_edge_verts%len()

            ! Store vertices into array
            call wake_edge_verts%get(i, m)
            this%wake_edge_top_verts(i) = m

        end do

        ! Allocate wake bottom vertices array
        allocate(this%wake_edge_bot_verts, source=this%wake_edge_top_verts)

        ! Store wake edge indices in array
        allocate(this%wake_edge_indices(this%N_wake_edges))
        do i=1,this%N_wake_edges

            ! Store index
            call wake_edges%get(i, m)
            this%wake_edge_indices(i) = m

        end do

        write(*,'(a, i3, a, i3, a)') "Done. Found ", this%N_wake_edges, " wake-shedding edges and ", &
                                     0, " other discontinuous edges."

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
        do j=1,this%N_wake_edges

            ! Get index of edge
            i = this%wake_edge_indices(j)

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

        end do
        write(*,*) "Done."

    end subroutine surface_mesh_set_up_mirroring


    subroutine surface_mesh_clone_vertices(this)
        ! Takes vertices which lie within discontinuous edges and splits them into two vertices.
        ! Handles rearranging of necessary dependencies.

        implicit none

        class(surface_mesh),intent(inout),target :: this

        integer :: i, j, k, m, n, N_clones, ind, new_ind, N_discont_verts, bottom_panel_ind, abutting_panel_ind, adj_vert_ind
        type(vertex),dimension(:),allocatable :: cloned_vertices, temp_vertices

        write(*,'(a)',advance='no') "     Cloning vertices at discontinuous edges..."

        ! Allocate array which will store which discontinuous vertices need to be cloned
        N_discont_verts = size(this%wake_edge_top_verts)

        ! Determine number of vertices which need to be cloned
        N_clones = 0
        do i=1,N_discont_verts

            if (this%vertices(this%wake_edge_top_verts(i))%needs_clone) then

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
            ind = this%wake_edge_top_verts(i)

            ! Check if this vertex needs to be cloned
            if (this%vertices(ind)%needs_clone) then

                ! Get index for the vertex clone
                new_ind = this%N_verts-N_clones+j ! Will be at position N_verts-N_clones+j in the new vertex array

                ! Initialize new vertex
                call this%vertices(new_ind)%init(this%vertices(ind)%loc, new_ind)

                ! Store clone's index in list of bottom vertices
                this%wake_edge_bot_verts(i) = new_ind

                ! Store number of adjacent wake-shedding edges (probably unecessary at this point, but let's be consistent)
                this%vertices(new_ind)%N_wake_edges = this%vertices(ind)%N_wake_edges

                ! Copy over mirroring properties
                this%vertices(new_ind)%mirrored_is_unique = this%vertices(ind)%mirrored_is_unique
                this%vertices(new_ind)%on_mirror_plane = this%vertices(ind)%on_mirror_plane

                ! Copy over adjacent panels
                do k=1,this%vertices(ind)%panels%len()

                    ! Get adjacent panel index from original vertex
                    call this%vertices(ind)%panels%get(k, abutting_panel_ind)

                    ! Copy to new vertex
                    call this%vertices(new_ind)%panels%append(abutting_panel_ind)

                    ! Copy to original vertex's panels_not_across_wake_edge list (bottom panels will be removed)
                    call this%vertices(ind)%panels_not_across_wake_edge%append(abutting_panel_ind)

                end do

                ! Copy over adjacent vertices
                do k=1,this%vertices(ind)%adjacent_vertices%len()

                    ! Get adjacent panel index from original vertex
                    call this%vertices(ind)%adjacent_vertices%get(k, adj_vert_ind)

                    ! Copy to new vertex
                    call this%vertices(new_ind)%adjacent_vertices%append(adj_vert_ind)

                end do

                ! Remove bottom panels from top vertex and give them to the bottom vertex
                do n=1,this%N_wake_edges

                    ! Get edge index
                    k = this%wake_edge_indices(n)

                    ! Check if this vertex belongs to this wake-shedding edge
                    if (this%edges(k)%verts(1) == ind .or. this%edges(k)%verts(2) == ind) then

                        ! Get bottom panel index
                        bottom_panel_ind = this%edges(k)%panels(2)

                        ! Remove bottom panel index from original vertex
                        call this%vertices(ind)%panels_not_across_wake_edge%delete(bottom_panel_ind)

                        ! Add to cloned vertex
                        if (.not. this%vertices(new_ind)%panels_not_across_wake_edge%is_in(bottom_panel_ind)) then
                            call this%vertices(new_ind)%panels_not_across_wake_edge%append(bottom_panel_ind)
                        end if

                        ! If there are any panels attached to this vertex and abutting the bottom panel, shift them over as well
                        do m=1,size(this%panels(bottom_panel_ind)%abutting_panels)

                            ! Get the index of the panel abutting this bottom panel
                            abutting_panel_ind = this%panels(bottom_panel_ind)%abutting_panels(m)

                            ! Check if it touches the current index
                            if (this%panels(abutting_panel_ind)%touches_vertex(ind)) then

                                ! Remove from original vertex
                                call this%vertices(ind)%panels_not_across_wake_edge%delete(abutting_panel_ind)

                                ! Add to cloned vertex
                                if (.not. this%vertices(new_ind)%panels_not_across_wake_edge%is_in(abutting_panel_ind)) then
                                    call this%vertices(new_ind)%panels_not_across_wake_edge%append(abutting_panel_ind)
                                end if

                            end if
                        end do

                    end if

                end do

                ! Update bottom panels to point to cloned vertex
                do k=1,this%vertices(new_ind)%panels_not_across_wake_edge%len()

                    ! Get panel index
                    call this%vertices(new_ind)%panels_not_across_wake_edge%get(k, bottom_panel_ind)

                    ! Update
                    call this%panels(bottom_panel_ind)%point_to_vertex_clone(this%vertices(new_ind))

                end do

                ! Update clone index
                j = j + 1

            else

                ! If this vertex did not need to be cloned, but it is on the mirror plane and its mirror is unique
                ! then the wake strength is partially determined by its mirror. But this is only in the case of an asymmetric flow.
                if (this%mirrored .and. this%asym_flow .and.  this%vertices(ind)%on_mirror_plane .and. &
                    this%vertices(ind)%mirrored_is_unique) then

                    this%wake_edge_bot_verts(i) = ind + this%N_verts

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
        real,dimension(3) :: vec_sum, normal
        integer :: i, j, N, ind

        write(*,'(a)',advance='no') "     Calculating vertex normals..."

        ! Loop through vertices
        do j=1,this%N_verts

            ! Loop through neighboring panels and compute the average of their normal vectors
            N = this%vertices(j)%panels%len()
            vec_sum = 0
            do i=1,N
                call this%vertices(j)%panels%get(i, ind)
                vec_sum = vec_sum + this%panels(ind)%normal
            end do

            ! For vertices on the mirror plane, the component normal to the plane should be zeroed
            if (this%mirrored .and. this%vertices(j)%on_mirror_plane) then

                vec_sum(this%mirror_plane) = 0.

            end if

            ! Normalize and store
            this%vertices(j)%normal = vec_sum/norm(vec_sum)

            ! Calculate average edge lengths for each vertex
            call this%vertices(j)%calc_average_edge_length(this%vertices)

        end do

        write(*,*) "Done."

    end subroutine surface_mesh_calc_vertex_normals


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


    subroutine surface_mesh_place_interior_control_points(this, offset)

        implicit none

        class(surface_mesh),intent(inout) :: this
        real,intent(in) :: offset

        integer :: i, j, N, ind
        real,dimension(3) :: sum
        real :: C_theta_2, offset_ratio

        if (doublet_order == 1) then

            ! Specify number of control points
            this%N_cp = this%N_verts

            ! Allocate memory
            allocate(this%control_points(this%N_verts,3))

            ! Calculate offset ratio such that the control point will remain within the body based on the minimum detected wake-shedding angle
            offset_ratio = 0.5*sqrt(0.5*(1.0+this%C_min_panel_angle))

            ! Loop through vertices
            do i=1,this%N_verts

                ! If the vertex is in a wake edge, it needs to be shifted off the normal slightly so that it is unique from its counterpart
                if (this%vertices(i)%N_wake_edges > 1) then

                    ! Loop through panels associated with this clone to get their average normal vector
                    N = this%vertices(i)%panels_not_across_wake_edge%len()
                    sum = 0
                    do j=1,N

                        ! Get panel index
                        call this%vertices(i)%panels_not_across_wake_edge%get(j, ind)

                        ! Add normal vector
                        sum = sum + this%panels(ind)%normal

                    end do

                    ! Normalize
                    sum = sum/norm(sum)

                    ! Place control point
                    this%control_points(i,:) = this%vertices(i)%loc &
                                               - offset * (this%vertices(i)%normal - offset_ratio * sum)*this%vertices(i)%l_avg

                ! If it's not in a wake-shedding edge (i.e. has no clone), then placement simply follows the normal vector
                else

                    this%control_points(i,:) = this%vertices(i)%loc-offset*this%vertices(i)%normal*this%vertices(i)%l_avg

                end if

            end do

        end if

        ! Calculate mirrored points, if necessary
        if (this%mirrored) then

            ! Allocate memory
            allocate(this%cp_mirrored(this%N_cp, 3))

            ! Calculate mirrors
            do i=1,this%N_cp
                this%cp_mirrored(i,:) = mirror_about_plane(this%control_points(i,:), this%mirror_plane)
            end do

        end if

    end subroutine surface_mesh_place_interior_control_points


    subroutine surface_mesh_output_results(this, body_file, wake_file, control_point_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: body_file, wake_file, control_point_file

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
            call body_vtk%write_cell_vectors(this%v(1:this%N_panels,:), "v")
            call body_vtk%write_point_scalars(this%mu(1:this%N_cp), "mu")
            call body_vtk%finish()

            write(*,*) "    Surface results written to: ", body_file
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
            call cp_vtk%write_points(this%control_points)
            call cp_vtk%write_vertices(this%control_points)
            call cp_vtk%write_point_scalars(this%phi_cp(1:this%N_cp), "phi")
            call cp_vtk%write_point_scalars(this%phi_cp_mu(1:this%N_cp), "phi_mu")
            call cp_vtk%write_point_scalars(this%phi_cp_sigma(1:this%N_cp), "phi_sigma")
            call cp_vtk%finish()

            write(*,*) "    Control point results written to: ", control_point_file
        end if
    
    end subroutine surface_mesh_output_results


end module surface_mesh_mod