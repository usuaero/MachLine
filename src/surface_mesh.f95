! A surface mesh type encompassing a body, wakes, and shocks
module surface_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use vtk_mod
    use vertex_mod
    use panel_mod
    use adt_mod
    use flow_mod
    use math_mod
    use wake_edge_mod

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels, N_wake_verts, N_wake_panels
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        type(vertex),allocatable,dimension(:) :: wake_vertices
        type(panel),allocatable,dimension(:) :: wake_panels
        type(list) :: wake_edge_vertices
        type(wake_edge),allocatable,dimension(:) :: wake_edges
        character(len=:),allocatable :: mesh_file
        type(alternating_digital_tree) :: vertex_tree
        real :: wake_shedding_angle, C_wake_shedding_angle, trefftz_distance
        integer :: N_wake_panels_streamwise, N_wake_edges
        real,dimension(:,:),allocatable :: control_points
        real,dimension(:),allocatable :: phi_cp
        real :: control_point_offset
        logical :: xy_sym, xz_sym, yz_sym ! Whether the mesh is to be mirrored about any planes

        contains

            procedure :: init => surface_mesh_init
            procedure :: load_adt => surface_mesh_load_adt
            procedure :: init_with_flow => surface_mesh_init_with_flow
            procedure :: output_results => surface_mesh_output_results
            procedure :: locate_wake_shedding_edges => surface_mesh_locate_wake_shedding_edges
            procedure :: clone_wake_shedding_vertices => surface_mesh_clone_wake_shedding_vertices
            procedure :: initialize_wake => surface_mesh_initialize_wake
            procedure :: calc_vertex_normals => surface_mesh_calc_vertex_normals
            procedure :: place_control_points => surface_mesh_place_control_points
            procedure :: update_wake => surface_mesh_update_wake

    end type surface_mesh
    

contains


    subroutine surface_mesh_init(this, settings)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        character(len=:),allocatable :: extension
        integer :: loc

        ! Get mesh file
        call json_get(settings, 'file', this%mesh_file)
        this%mesh_file = trim(this%mesh_file)
        write(*,*)
        write(*,*) "    Reading surface mesh in from file: ", this%mesh_file

        ! Determine the type of mesh file
        loc = index(this%mesh_file, '.')
        extension = this%mesh_file(loc:len(this%mesh_file))

        ! Load vtk
        if (extension .eq. '.vtk') then
            call load_surface_vtk(this%mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        end if

        ! Display mesh info
        write(*,*)
        write(*,*) "    Surface mesh has", this%N_verts, "vertices and", this%N_panels, "panels."

        ! Load into adt
        !call this%load_adt()

        ! Get symmetry settings
        call json_xtnsn_get(settings, 'symmetry.xy', this%xy_sym, .false.)
        call json_xtnsn_get(settings, 'symmetry.xz', this%xz_sym, .false.)
        call json_xtnsn_get(settings, 'symmetry.yz', this%yz_sym, .false.)

        ! Store settings for wake models
        call json_xtnsn_get(settings, 'wake_model.wake_shedding_angle', this%wake_shedding_angle, 90.0)
        call json_xtnsn_get(settings, 'wake_model.trefftz_distance', this%trefftz_distance, 100.0)
        call json_xtnsn_get(settings, 'wake_model.N_panels', this%N_wake_panels_streamwise, 20)
        this%C_wake_shedding_angle = cos(this%wake_shedding_angle*pi/180.0)

        ! Store boundary condition settings
        call json_get(settings, 'boundary_conditions.control_point_offset', this%control_point_offset)
    
    end subroutine surface_mesh_init


    subroutine surface_mesh_load_adt(this)

        implicit none

        class(surface_mesh),intent(inout) :: this
        real,dimension(3) :: p_min, p_max
        integer :: i

        write(*,*)
        write(*,'(a)',advance='no') "     Loading vertices into ADT..."

        ! Determine bounds of alternating digital tree
        p_min = this%vertices(1)%loc
        p_max = this%vertices(1)%loc
        do i=2,this%N_verts

            ! Check mins
            p_min(1) = min(this%vertices(i)%loc(1), p_min(1))
            p_min(2) = min(this%vertices(i)%loc(2), p_min(2))
            p_min(3) = min(this%vertices(i)%loc(3), p_min(3))

            ! Check maxs
            p_max(1) = max(this%vertices(i)%loc(1), p_max(1))
            p_max(2) = max(this%vertices(i)%loc(2), p_max(2))
            p_max(3) = max(this%vertices(i)%loc(3), p_max(3))

        end do
        
        ! Store
        this%vertex_tree%p_min = p_min
        this%vertex_tree%p_max = p_max

        ! Load vertices into alternating digital tree
        do i=1,this%N_verts
            call this%vertex_tree%add(this%vertices(i))
        end do
        write(*,*) "Done."


    end subroutine surface_mesh_load_adt


    subroutine surface_mesh_init_with_flow(this, freestream_flow)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow

        ! Call subroutines which initialize flow-dependent properties
        call this%locate_wake_shedding_edges(freestream_flow)
        call this%clone_wake_shedding_vertices()
        call this%calc_vertex_normals()
        call this%initialize_wake(freestream_flow)
        call this%place_control_points()
    
    end subroutine surface_mesh_init_with_flow


    subroutine surface_mesh_locate_wake_shedding_edges(this, freestream_flow)
        ! Locates wake-shedding edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        integer :: i, j, m, n, mm, temp, top_panel, bottom_panel
        real,dimension(3) :: d
        integer,dimension(2) :: shared_verts
        type(list) :: wake_edge_starts, wake_edge_stops, top_panels, bottom_panels
        logical :: abutting, already_found_shared, is_wake_edge
        real :: distance

        write(*,*)
        write(*,'(a)',advance='no') "     Locating wake-shedding edges..."

        ! Loop through each pair of panels
        this%N_wake_edges = 0
        do i=1,this%N_panels
            do j=i+1,this%N_panels

                ! Initialize for this panel pair
                already_found_shared = .false.
                abutting = .false.

                ! Check if the panels are abutting
                abutting_loop: do m=1,this%panels(i)%N
                    do n=1,this%panels(j)%N

                        ! Get distance between vertices
                        d = this%panels(i)%get_vertex_loc(m)-this%panels(j)%get_vertex_loc(n) ! More robust than checking vertex indices
                        distance = norm(d)

                        ! Check distance
                        if (distance < 1e-10) then

                            ! First shared vertex
                            if (.not. already_found_shared) then
                                already_found_shared = .true.
                                shared_verts(1) = this%panels(i)%get_vertex_index(m)
                                mm = m

                            ! Previously found a shared vertex
                            else
                                abutting = .true.
                                shared_verts(2) = this%panels(i)%get_vertex_index(m)
                                exit abutting_loop
                            end if

                        end if

                    end do
                end do abutting_loop

                if (abutting) then

                    is_wake_edge = .false.

                    ! Check angle between panels
                    if (inner(this%panels(i)%normal, this%panels(j)%normal) < this%C_wake_shedding_angle) then

                        ! Check angle of panel normal with freestream
                        if (inner(this%panels(i)%normal, freestream_flow%V_inf) > 0.0 .or. &
                            inner(this%panels(j)%normal, freestream_flow%V_inf) > 0.0) then

                            ! Check order the vertices were stored in
                            if (mm == 1 .and. m == this%panels(i)%N) then

                                ! Rearrange so the wake-shedding edge proceeds in the counterclockwise direction around the "top" panel
                                ! I'll use "top" and "bottom" to refer to panels neighboring a wake-shedding edge. These terms are arbitrary but consistent.
                                temp = shared_verts(2)
                                shared_verts(2) = shared_verts(1)
                                shared_verts(1) = temp

                            end if

                            ! Update number of wake-shedding edges
                            this%N_wake_edges = this%N_wake_edges + 1
                            is_wake_edge = .true.

                            ! Store in starts and stops list
                            call wake_edge_starts%append(shared_verts(1))
                            call wake_edge_stops%append(shared_verts(2))

                            ! Store top and bottom panels (i is top, j is bottom)
                            call top_panels%append(i)
                            call bottom_panels%append(j)

                            ! Store the fact that these vertices belong to a wake-shedding edge
                            if (this%vertices(shared_verts(1))%on_wake_edge) then

                                ! If it's already on a wake-shedding edge, then this means it's also in one
                                this%vertices(shared_verts(1))%in_wake_edge = .true.

                            else

                                ! Add the first time
                                this%vertices(shared_verts(1))%on_wake_edge = .true.
                                call this%wake_edge_vertices%append(shared_verts(1))
                                this%vertices(shared_verts(1))%index_in_wake_vertices = this%wake_edge_vertices%len()

                            end if

                            ! Do the same for the other vertex
                            if (this%vertices(shared_verts(2))%on_wake_edge) then

                                this%vertices(shared_verts(2))%in_wake_edge = .true.

                            else

                                this%vertices(shared_verts(2))%on_wake_edge = .true.
                                call this%wake_edge_vertices%append(shared_verts(2))
                                this%vertices(shared_verts(2))%index_in_wake_vertices = this%wake_edge_vertices%len()

                            end if

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

        ! Create array of wake-shedding edges
        allocate(this%wake_edges(this%N_wake_edges))
        do i=1,this%N_wake_edges

            ! Get indices of starting and ending vertices
            call wake_edge_starts%get(i, m)
            call wake_edge_stops%get(i, n)
            call top_panels%get(i, top_panel)
            call bottom_panels%get(i, bottom_panel)

            ! Store
            call this%wake_edges(i)%init(m, n, top_panel, bottom_panel)

        end do

        write(*,*) "Done. Found", this%N_wake_edges, "wake-shedding edges."

    end subroutine surface_mesh_locate_wake_shedding_edges


    subroutine surface_mesh_clone_wake_shedding_vertices(this)
        ! Takes vertices which lie within wake-shedding edges and splits them into two.
        ! Handles rearranging of necessary dependencies.

        implicit none

        class(surface_mesh),intent(inout),target :: this
        integer :: i, j, k, m, N_clones, ind, new_ind, N_wake_verts, bottom_panel_ind, abutting_panel_ind
        type(vertex),dimension(:),allocatable :: cloned_vertices, temp_vertices
        logical,dimension(:),allocatable :: need_cloned

        write(*,*)
        write(*,'(a)',advance='no') "     Cloning vertices on wake-shedding edges..."

        ! Allocate array which will store which wake-shedding vertices need to be cloned
        allocate(need_cloned(this%wake_edge_vertices%len()))

        ! Determine number of vertices which need to be cloned
        N_wake_verts = this%wake_edge_vertices%len()
        N_clones = 0
        do i=1,N_wake_verts

            ! Get the vertex index
            call this%wake_edge_vertices%get(i, ind)

            ! Check if it is *in* a wake-shedding edge
            if (this%vertices(ind)%in_wake_edge) then
                N_clones = N_clones + 1
                need_cloned(i) = .true.
            end if
        end do

        ! Extend allocation of vertex array
        allocate(temp_vertices, source=this%vertices)
        deallocate(this%vertices)
        allocate(this%vertices(this%N_verts + N_clones))
        this%vertices(1:this%N_verts) = temp_vertices
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

            ! Check if this vertex needs to be cloned
            if (need_cloned(i)) then

                ! Get information for the vertex clone
                call this%wake_edge_vertices%get(i, ind)
                new_ind = this%N_verts-N_clones+j ! Will be at position N_verts-N_clones+j in the new vertex array

                ! Initialize new vertex
                call this%vertices(new_ind)%init(this%vertices(ind)%loc, new_ind)

                ! Store that it is on and in a wake-shedding edge (probably unecessary at this point, but let's be consistent)
                this%vertices(new_ind)%on_wake_edge = .true.
                this%vertices(new_ind)%in_wake_edge = .true.

                ! Remove bottom panels from top vertex and give them to the bottom vertex
                do k=1,this%N_wake_edges

                    ! Check if this vertex belongs to this wake-shedding edge
                    if (this%wake_edges(k)%i1 == ind .or. this%wake_edges(k)%i2 == ind) then

                        ! Get bottom panel index
                        bottom_panel_ind = this%wake_edges(k)%bottom_panel

                        ! Remove bottom panel index from original vertex
                        call this%vertices(ind)%panels%delete(bottom_panel_ind)

                        ! Add to cloned vertex
                        if (.not. this%vertices(new_ind)%panels%is_in(bottom_panel_ind)) then
                            call this%vertices(new_ind)%panels%append(bottom_panel_ind)
                        end if

                        ! If there are any panels attached to this vertex and abutting the bottom panel, shift them over as well
                        do m=1,this%panels(bottom_panel_ind)%abutting_panels%len()

                            ! Get the index of the panel abutting this bottom panel
                            call this%panels(bottom_panel_ind)%abutting_panels%get(m, abutting_panel_ind)

                            ! Check if it touches the current index
                            if (this%panels(abutting_panel_ind)%touches_vertex(ind)) then

                                ! Remove from original vertex
                                call this%vertices(ind)%panels%delete(abutting_panel_ind)

                                ! Add to cloned vertex
                                if (.not. this%vertices(new_ind)%panels%is_in(abutting_panel_ind)) then
                                    call this%vertices(new_ind)%panels%append(abutting_panel_ind)
                                end if

                            end if
                        end do

                    end if

                end do

                ! Update bottom panels to point to cloned vertex
                do k=1,this%vertices(new_ind)%panels%len()

                    ! Get panel index
                    call this%vertices(new_ind)%panels%get(k, bottom_panel_ind)

                    ! Update
                    call this%panels(bottom_panel_ind)%point_to_vertex_clone(this%vertices(new_ind))

                end do

                ! Update clone index
                j = j + 1

            end if

        end do

        write(*,*) "Done. Cloned", N_clones, "vertices. Mesh now has", this%N_verts, "vertices."


    end subroutine surface_mesh_clone_wake_shedding_vertices


    subroutine surface_mesh_calc_vertex_normals(this)
        ! Initializes the normal vectors associated with each vertex.
        ! Must be called only once wake-shedding edge vertices have been cloned.

        implicit none

        class(surface_mesh),intent(inout) :: this
        real,dimension(3) :: sum
        integer :: i, j, N, ind

        write(*,*)
        write(*,'(a)',advance='no') "     Calculating vertex normals..."

        ! Loop through vertices
        do j=1,this%N_verts

            ! Loop through neighboring panels and compute the average of their normal vectors
            N = this%vertices(j)%panels%len()
            sum = 0
            do i=1,N
                call this%vertices(j)%panels%get(i, ind)
                sum = sum + this%panels(ind)%normal
            end do

            ! Store
            this%vertices(j)%normal = sum/norm(sum)

        end do

        write(*,*) "Done."

    end subroutine surface_mesh_calc_vertex_normals


    subroutine surface_mesh_initialize_wake(this, freestream_flow)
        ! Creates wake panels and sets their initial shape. Handles vertex association.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        real :: distance, vertex_separation
        real,dimension(3) :: loc, start
        integer :: i, j, ind, kutta_vert_ind, i_start, i_stop, i1, i2, i3, i4
        integer :: N_wake_edge_verts, N_wake_verts, N_wake_panels

        ! Initialize wake
        write(*,*)
        write(*,'(a)',advance='no') "     Initializing wake..."

        ! Determine sizes
        N_wake_edge_verts = this%wake_edge_vertices%len()
        this%N_wake_verts = N_wake_edge_verts*(this%N_wake_panels_streamwise+1)
        this%N_wake_panels = this%N_wake_edges*this%N_wake_panels_streamwise*2

        ! Allocate storage
        allocate(this%wake_vertices(this%N_wake_verts))
        allocate(this%wake_panels(this%N_wake_panels))

        ! Determine vertex placement
        do i=1,N_wake_edge_verts

            ! Determine distance from origin to wake-shedding vertex in direction of the flow
            call this%wake_edge_vertices%get(i, kutta_vert_ind)
            start = this%vertices(kutta_vert_ind)%loc
            distance = this%trefftz_distance-inner(start, freestream_flow%u_inf)

            ! Determine vertex separation
            vertex_separation = distance/this%N_wake_panels_streamwise

            ! Place vertices
            do j=1,this%N_wake_panels_streamwise+1

                ! Determine location
                ind = (i-1)*(this%N_wake_panels_streamwise+1)+j
                loc = start+vertex_separation*(j-1)*freestream_flow%u_inf

                ! Initialize vertex
                call this%wake_vertices(ind)%init(loc, ind)

                ! Set parent index
                this%wake_vertices(ind)%parent = kutta_vert_ind

            end do
        end do

        ! Initialize wake panels
        do i=1,this%N_wake_edges

            ! Determine which Kutta vertices this panel lies between
            i_start = this%vertices(this%wake_edges(i)%i1)%index_in_wake_vertices
            i_stop = this%vertices(this%wake_edges(i)%i2)%index_in_wake_vertices

            ! Create panels heading downstream
            do j=1,this%N_wake_panels_streamwise

                ! Determine index of first triangular panel
                ind = (i-1)*this%N_wake_panels_streamwise*2+2*j-1

                ! Determine vertex indices
                i1 = (i_start-1)*(this%N_wake_panels_streamwise+1)+j
                i2 = (i_start-1)*(this%N_wake_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(this%N_wake_panels_streamwise+1)+j+1

                ! Initialize
                call this%wake_panels(ind)%init(this%wake_vertices(i1),&
                                                this%wake_vertices(i2),&
                                                this%wake_vertices(i3),&
                                                i1, i2, i3, ind)

                ! Determine index of second triangular panel
                ind = (i-1)*this%N_wake_panels_streamwise*2+2*j

                ! Determine vertex indices
                i1 = (i_start-1)*(this%N_wake_panels_streamwise+1)+j
                i2 = (i_stop-1)*(this%N_wake_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(this%N_wake_panels_streamwise+1)+j

                ! Initialize
                call this%wake_panels(ind)%init(this%wake_vertices(i1),&
                                                this%wake_vertices(i2),&
                                                this%wake_vertices(i3),&
                                                i1, i2, i3, ind)

            end do
        end do

        write(*,*) "Done. Created", this%N_wake_verts, "wake vertices and", this%N_wake_panels, "wake panels."

    end subroutine surface_mesh_initialize_wake


    subroutine surface_mesh_place_control_points(this)

        implicit none

        class(surface_mesh),intent(inout) :: this
        integer :: i

        ! Allocate memory
        allocate(this%control_points(this%N_verts,3))
        allocate(this%phi_cp(this%N_verts), source=0.0)

        ! Loop through vertices
        do i=1,this%N_verts

            ! If it's not in a wake-shedding edge (i.e. has no clone), then placement simply follows the normal vector
            if (.not. this%vertices(i)%in_wake_edge) then

                this%control_points(i,:) = this%vertices(i)%loc-this%control_point_offset*this%vertices(i)%normal

            ! Otherwise, we have to make sure control points stay within the body
            else

                this%control_points(i,:) = this%vertices(i)%loc-this%control_point_offset*this%vertices(i)%normal

            end if

        end do


    end subroutine surface_mesh_place_control_points


    subroutine surface_mesh_update_wake(this)

        implicit none

        class(surface_mesh),intent(inout) :: this
        integer :: i

        ! Update panel properties since vertices were moved
        do i=1,this%N_wake_panels
            call this%wake_panels(i)%calc_derived_properties
        end do
    
    end subroutine surface_mesh_update_wake


    subroutine surface_mesh_output_results(this, body_file, wake_file, control_point_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: body_file, wake_file, control_point_file

        ! Write out data for body
        call write_surface_vtk(body_file, this%vertices, this%panels)
        
        ! Write out data for wake
        call write_surface_vtk(wake_file, this%wake_vertices, this%wake_panels)
        
        ! Write out data for control points
        call write_point_vtk(control_point_file, this%control_points, this%phi_cp)
    
    end subroutine surface_mesh_output_results

end module surface_mesh_mod