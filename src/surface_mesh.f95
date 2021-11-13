! A surface mesh type encompassing a body, wakes, and shocks
module surface_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use vtk_mod
    use vertex_mod
    use panel_mod
    use flow_mod
    use math_mod
    use wake_edge_mod
    use wake_mesh_mod

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels, N_cp, N_unique_mirrored_verts
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        type(wake_mesh) :: wake
        integer,allocatable,dimension(:) :: wake_edge_top_verts, wake_edge_bot_verts
        type(wake_edge),allocatable,dimension(:) :: wake_edges
        real :: wake_shedding_angle, C_wake_shedding_angle, trefftz_distance, C_min_wake_shedding_angle
        integer :: N_wake_panels_streamwise
        real,dimension(:,:),allocatable :: control_points
        real,dimension(:),allocatable :: phi_cp, phi_cp_sigma, phi_cp_mu ! Induced potentials at control points
        real,dimension(:),allocatable :: C_p ! Surface pressure coefficients
        real,dimension(:,:),allocatable :: V ! Surface velocities
        real :: control_point_offset
        logical :: mirrored ! Whether the mesh is to be mirrored about any planes
        integer :: mirror_plane ! Index of the plane across which the mesh is mirrored (1: yz, 2: xz, 3: xy); this is the index of the normal to that plane
        logical :: mirrored_and_asym ! Whether the mesh is mirrored and the flow is asymmetric about that plane
        real,dimension(:),allocatable :: mu, sigma ! Singularity strengths

        contains

            procedure :: init => surface_mesh_init
            procedure :: init_with_flow => surface_mesh_init_with_flow
            procedure :: output_results => surface_mesh_output_results
            procedure :: locate_wake_shedding_edges => surface_mesh_locate_wake_shedding_edges
            procedure :: clone_wake_shedding_vertices => surface_mesh_clone_wake_shedding_vertices
            procedure :: set_up_mirroring => surface_mesh_set_up_mirroring
            procedure :: calc_vertex_normals => surface_mesh_calc_vertex_normals
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
        write(*,*)
        write(*,*) "    User has selected:", doublet_order, "-order doublet panels and", source_order, "-order source panels."

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
        write(*,*)
        write(*,*) "    Reading surface mesh in from file: ", mesh_file

        ! Determine the type of mesh file
        loc = index(mesh_file, '.')
        extension = mesh_file(loc:len(mesh_file))

        ! Load vtk
        if (extension == '.vtk') then
            call load_surface_vtk(mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        else
            write(*,*) "MFTran cannot read ", extension, " type mesh files. Quitting..."
            stop
        end if

        ! Display mesh info
        write(*,*)
        write(*,*) "    Surface mesh has", this%N_verts, "vertices and", this%N_panels, "panels."

        ! Get mirroring
        call json_xtnsn_get(settings, 'mirror_about', mirror_plane, "none")
        this%mirrored = .true.
        select case (mirror_plane)
        case ("xy")
            this%mirror_plane = 3
        case ("xz")
            this%mirror_plane = 2
        case ("yz")
            this%mirror_plane = 1
        case default
            this%mirror_plane = 0
            this%mirrored = .false.
        end select

        ! Store settings for wake models
        call json_xtnsn_get(settings, 'wake_model.wake_shedding_angle', this%wake_shedding_angle, 90.0) ! Maximum allowable angle between panel normals without having separation
        call json_xtnsn_get(settings, 'wake_model.trefftz_distance', this%trefftz_distance, 100.0) ! Distance from origin to wake termination
        call json_xtnsn_get(settings, 'wake_model.N_panels', this%N_wake_panels_streamwise, 20)
        this%C_wake_shedding_angle = cos(this%wake_shedding_angle*pi/180.0)

    end subroutine surface_mesh_init


    subroutine surface_mesh_init_with_flow(this, freestream)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        ! Figure out wake-shedding edges
        call this%locate_wake_shedding_edges(freestream)

        ! Check mirroring and flow symmetry condition
        this%mirrored_and_asym = .false.
        if (this%mirrored) then
            if (.not. freestream%sym_about(this%mirror_plane)) then
                this%mirrored_and_asym = .true.
            end if
        end if

        ! Set up mirroring
        if (this%mirrored) then
            call this%set_up_mirroring()
        end if

        ! Clone necessary vertices and calculate normals (for placing control points)
        call this%clone_wake_shedding_vertices()
        call this%calc_vertex_normals()

        ! Initialize wake
        call this%wake%init(freestream, this%wake_edge_top_verts, &
                            this%wake_edge_bot_verts, this%wake_edges, &
                            this%N_wake_panels_streamwise, this%vertices, &
                            this%trefftz_distance, this%mirrored_and_asym, &
                            this%mirror_plane)

        ! Clean up
        deallocate(this%wake_edge_top_verts)
        deallocate(this%wake_edge_bot_verts)
    
    end subroutine surface_mesh_init_with_flow


    subroutine surface_mesh_locate_wake_shedding_edges(this, freestream)
        ! Locates wake-shedding edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream

        integer :: i, j, m, n, mm, temp, top_panel, bottom_panel
        integer :: N_wake_edges = 0
        integer,dimension(2) :: shared_verts
        type(list) :: wake_edge_starts, wake_edge_stops, top_panels, bottom_panels, wake_edge_verts
        logical :: abutting, already_found_shared, is_wake_edge
        real :: distance, C_angle

        write(*,*)
        write(*,'(a)',advance='no') "     Locating wake-shedding edges..."

        ! We need to store the minimum angle between two panels in order to place control points within the body
        this%C_min_wake_shedding_angle = this%C_wake_shedding_angle

        ! Loop through each pair of panels
        do i=1,this%N_panels
            do j=i+1,this%N_panels

                ! Initialize for this panel pair
                already_found_shared = .false.
                abutting = .false.

                ! Check if the panels are abutting
                abutting_loop: do m=1,this%panels(i)%N
                    do n=1,this%panels(j)%N

                        ! Get distance between vertices
                        distance = dist(this%panels(i)%get_vertex_loc(m), this%panels(j)%get_vertex_loc(n)) ! More robust than checking vertex indices; mesh may not be ideal

                        ! Check distance
                        if (distance < 1e-10) then

                            ! Previously found a shared vertex
                            if (already_found_shared) then
                                abutting = .true.
                                shared_verts(2) = this%panels(i)%get_vertex_index(m)
                                exit abutting_loop

                            ! First shared vertex
                            else
                                already_found_shared = .true.
                                shared_verts(1) = this%panels(i)%get_vertex_index(m)
                                mm = m

                            end if
                        end if

                    end do
                end do abutting_loop

                ! Perform checks on panels we know are abutting
                if (abutting) then

                    ! Store adjacent vertices
                    if (.not. this%vertices(shared_verts(1))%adjacent_vertices%is_in(shared_verts(2))) then
                        call this%vertices(shared_verts(1))%adjacent_vertices%append(shared_verts(2))
                    end if
                    if (.not. this%vertices(shared_verts(2))%adjacent_vertices%is_in(shared_verts(1))) then
                        call this%vertices(shared_verts(2))%adjacent_vertices%append(shared_verts(1))
                    end if

                    ! Reinitialize check for wake edge
                    is_wake_edge = .false.

                    ! Check angle between panels
                    C_angle = inner(this%panels(i)%normal, this%panels(j)%normal)
                    if (C_angle < this%C_wake_shedding_angle) then

                        ! Check angle of panel normal with freestream
                        if (inner(this%panels(i)%normal, freestream%V_inf) > 0.0 .or. &
                            inner(this%panels(j)%normal, freestream%V_inf) > 0.0) then

                            ! At this point, the panels are abutting, have a small enough angle,
                            ! and at least one of them is pointing downstream. Thus, this is a
                            ! wake-shedding edge.

                            ! Update minimum angle
                            this%C_min_wake_shedding_angle = min(C_angle, this%C_min_wake_shedding_angle)

                            ! Check order the vertices were stored in
                            if (mm == 1 .and. m == this%panels(i)%N) then

                                ! Rearrange so the wake-shedding edge proceeds in the counterclockwise direction around the "top" panel
                                ! I'll use "top" and "bottom" to refer to panels neighboring a wake-shedding edge. These terms are arbitrary but consistent.
                                temp = shared_verts(2)
                                shared_verts(2) = shared_verts(1)
                                shared_verts(1) = temp

                            end if

                            ! Update number of wake-shedding edges
                            N_wake_edges = N_wake_edges + 1
                            is_wake_edge = .true.

                            ! Store in starts and stops list
                            call wake_edge_starts%append(shared_verts(1))
                            call wake_edge_stops%append(shared_verts(2))

                            ! Store top and bottom panels (i is top, j is bottom)
                            call top_panels%append(i)
                            call bottom_panels%append(j)

                            ! If this vertex does not already belong to a wake-shedding edge, add it to the list of wake edge vertices
                            if (this%vertices(shared_verts(1))%N_wake_edges == 0) then 

                                ! Add the first time
                                call wake_edge_verts%append(shared_verts(1))
                                this%vertices(shared_verts(1))%index_in_wake_vertices = wake_edge_verts%len()

                            else

                                ! It is in an edge, so it will likely need to be cloned
                                this%vertices(shared_verts(1))%needs_clone = .true.

                            end if
                            
                            ! Update number of wake edges touching this vertex
                            this%vertices(shared_verts(1))%N_wake_edges = this%vertices(shared_verts(1))%N_wake_edges + 1

                            ! Do the same for the other vertex
                            if (this%vertices(shared_verts(2))%N_wake_edges == 0) then 

                                ! Add the first time
                                call wake_edge_verts%append(shared_verts(2))
                                this%vertices(shared_verts(2))%index_in_wake_vertices = wake_edge_verts%len()

                            else

                                ! It is in an edge, so it will likely need to be cloned
                                this%vertices(shared_verts(2))%needs_clone = .true.

                            end if
                            
                            ! Update number of wake edges touching this vertex
                            this%vertices(shared_verts(2))%N_wake_edges = this%vertices(shared_verts(2))%N_wake_edges + 1

                            ! Store the fact that the panels have a Kutta edge
                            this%panels(i)%on_wake_edge = .true.
                            this%panels(j)%on_wake_edge = .true.

                            ! Store opposing panels
                            call this%panels(i)%opposing_panels%append(j)
                            call this%panels(j)%opposing_panels%append(i)

                        end if
                    end if

                    ! If abutting but not a wake-shedding edge
                    if (.not. is_wake_edge) then

                        ! Add to each others' list
                        call this%panels(i)%abutting_panels%append(j)
                        call this%panels(j)%abutting_panels%append(i)

                    end if

                end if

            end do
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

        ! Create array of wake-shedding edges
        allocate(this%wake_edges(N_wake_edges))
        do i=1,N_wake_edges

            ! Get indices of starting and ending vertices
            call wake_edge_starts%get(i, m)
            call wake_edge_stops%get(i, n)
            call top_panels%get(i, top_panel)
            call bottom_panels%get(i, bottom_panel)

            ! Store
            call this%wake_edges(i)%init(m, n, top_panel, bottom_panel)

        end do

        write(*,*) "Done. Found", N_wake_edges, "wake-shedding edges."

    end subroutine surface_mesh_locate_wake_shedding_edges


    subroutine surface_mesh_set_up_mirroring(this)
        ! Sets up information necessary for mirroring the mesh

        implicit none

        class(surface_mesh),intent(inout) :: this

        integer :: i, m, n

        write(*,*)
        write(*,'(a)',advance='no') "     Setting up mesh mirror..."

        ! Search for vertices lying on mirror plane
        do i=1,this%N_verts

            ! Check coordinate normal to mirror plane
            if (abs(this%vertices(i)%loc(this%mirror_plane))<1e-12) then

                ! The vertex is on the mirror plane
                this%vertices(i)%on_mirror_plane = .true.

                ! If the vertex doesn't belong to a wake edge, then its mirror will not be unique
                if (this%vertices(i)%N_wake_edges == 0) then
                    this%vertices(i)%mirrored_is_unique = .false.
                end if

            end if
            
        end do

        ! Check if any wake-shedding edges touch the mirror plane
        do i=1,size(this%wake_edges)

            ! Get endpoint indices
            m = this%wake_edges(i)%i1
            n = this%wake_edges(i)%i2

            ! If a given wake edge has only one of its endpoints lying on the mirror plane, then that endpoint has another
            ! adjacent edge, since the edge will be mirrored across that plane. This endpoint will need a clone, but it's
            ! mirrored vertex will be the same

            ! Check for edge touching mirror plane at only one end
            if (this%vertices(m)%on_mirror_plane .neqv. this%vertices(n)%on_mirror_plane) then

                ! Only add wake edge if the endpoint is on the mirror plane
                if (this%vertices(m)%on_mirror_plane) then

                    this%vertices(m)%N_wake_edges = this%vertices(m)%N_wake_edges + 1
                    this%vertices(m)%needs_clone = .true.
                    this%vertices(m)%mirrored_is_unique = .false.

                else

                    this%vertices(n)%N_wake_edges = this%vertices(n)%N_wake_edges + 1
                    this%vertices(n)%needs_clone = .true.
                    this%vertices(n)%mirrored_is_unique = .false.

                end if

            ! If the wake edge has both endpoints lying on the mirror plane, then these vertices need no clones, but the mirrored
            ! vertices will still be unique
            else if (this%vertices(m)%on_mirror_plane .and. this%vertices(n)%on_mirror_plane) then

                ! The mirrored vertices will function as clones
                this%vertices(m)%needs_clone = .false.
                this%vertices(n)%needs_clone = .false.

                ! Store that this edge lies on the mirror plane
                this%wake_edges(i)%on_mirror_plane = .true.

            end if

        end do
        write(*,*) "Done."

    end subroutine surface_mesh_set_up_mirroring


    subroutine surface_mesh_clone_wake_shedding_vertices(this)
        ! Takes vertices which lie within wake-shedding edges and splits them into two vertices.
        ! Handles rearranging of necessary dependencies.

        implicit none

        class(surface_mesh),intent(inout),target :: this

        integer :: i, j, k, m, n, N_clones, ind, new_ind, N_wake_verts, bottom_panel_ind, abutting_panel_ind, adj_vert_ind
        type(vertex),dimension(:),allocatable :: cloned_vertices, temp_vertices

        write(*,*)
        write(*,'(a)',advance='no') "     Cloning vertices on wake-shedding edges..."

        ! Allocate array which will store which wake-shedding vertices need to be cloned
        N_wake_verts = size(this%wake_edge_top_verts)

        ! Determine number of vertices which need to be cloned
        N_clones = 0
        do i=1,N_wake_verts

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
        do i=1,N_wake_verts

            ! Get index of vertex to be cloned
            ind = this%wake_edge_top_verts(i)

            ! Check if this vertex needs to be cloned
            if (this%vertices(ind)%needs_clone) then

                ! Clone a vertex for each attached edge minus one (i.e. if three wake edges meet at one point, two clones will need to be made)
                do n=1,this%vertices(ind)%N_wake_edges-1

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
                    do k=1,size(this%wake_edges)

                        ! Check if this vertex belongs to this wake-shedding edge
                        if (this%wake_edges(k)%i1 == ind .or. this%wake_edges(k)%i2 == ind) then

                            ! Get bottom panel index
                            bottom_panel_ind = this%wake_edges(k)%bottom_panel

                            ! Remove bottom panel index from original vertex
                            call this%vertices(ind)%panels_not_across_wake_edge%delete(bottom_panel_ind)

                            ! Add to cloned vertex
                            if (.not. this%vertices(new_ind)%panels_not_across_wake_edge%is_in(bottom_panel_ind)) then
                                call this%vertices(new_ind)%panels_not_across_wake_edge%append(bottom_panel_ind)
                            end if

                            ! If there are any panels attached to this vertex and abutting the bottom panel, shift them over as well
                            do m=1,this%panels(bottom_panel_ind)%abutting_panels%len()

                                ! Get the index of the panel abutting this bottom panel
                                call this%panels(bottom_panel_ind)%abutting_panels%get(m, abutting_panel_ind)

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

                end do

            else

                ! If this vertex did not need to be cloned, but it is on the mirror plane and its mirror is unique
                ! then the wake strength is partially determined by its mirror. But this is only in the case of an asymmetric flow.
                if (this%mirrored_and_asym .and. &
                    this%vertices(ind)%on_mirror_plane .and. &
                    this%vertices(ind)%mirrored_is_unique) then

                    this%wake_edge_bot_verts(i) = ind + this%N_verts

                end if

            end if

        end do

        ! Calculate average edge lengths for each vertex
        do i=1,this%N_verts
            call this%vertices(i)%calc_average_edge_length(this%vertices)
        end do

        write(*,*) "Done. Cloned", N_clones, "vertices. Mesh now has", this%N_verts, "vertices."

    end subroutine surface_mesh_clone_wake_shedding_vertices


    subroutine surface_mesh_calc_vertex_normals(this)
        ! Initializes the normal vectors associated with each vertex.
        ! Must be called only once wake-shedding edge vertices have been cloned.

        implicit none

        class(surface_mesh),intent(inout) :: this
        real,dimension(3) :: vec_sum, normal
        integer :: i, j, N, ind

        write(*,*)
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

        end do

        write(*,*) "Done."

    end subroutine surface_mesh_calc_vertex_normals


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
            offset_ratio = 0.5*sqrt(0.5*(1.0+this%C_min_wake_shedding_angle))

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

    end subroutine surface_mesh_place_interior_control_points


    subroutine surface_mesh_output_results(this, body_file, wake_file, control_point_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: body_file, wake_file, control_point_file

        ! Write out data for body
        call write_surface_vtk(body_file, this%vertices, this%panels, this%sigma, this%mu, .false., this%V, this%C_p)
        write(*,*) "        Surface results written to: ", body_file
        
        ! Write out data for wake
        if (this%wake%N_panels > 0) then
            call write_surface_vtk(wake_file, this%wake%vertices, this%wake%panels, this%sigma, this%mu, .true.)
            write(*,*) "        Wake results written to: ", wake_file
        else
            write(*,*) "        No wake to export."
        end if
        
        ! Write out data for control points
        call write_point_vtk(control_point_file, this%control_points, this%phi_cp, this%phi_cp_mu, this%phi_cp_sigma)
        write(*,*) "        Control point results written to: ", control_point_file
    
    end subroutine surface_mesh_output_results

end module surface_mesh_mod