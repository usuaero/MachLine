! Types and subroutines for meshes
module surface_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use vtk_mod
    use vertex_mod
    use panel_mod
    use adt_mod
    use flow_mod
    use math_mod
    use kutta_edge_mod

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels, N_wake_verts, N_wake_panels
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        type(vertex),allocatable,dimension(:) :: wake_vertices
        type(panel),allocatable,dimension(:) :: wake_panels
        type(list) :: kutta_vertices
        type(kutta_edge),allocatable,dimension(:) :: kutta_edges
        character(len=:),allocatable :: mesh_file
        type(alternating_digital_tree) :: vertex_tree
        real :: kutta_angle, C_kutta_angle, trefftz_distance
        integer :: N_wake_panels_streamwise, N_kutta_edges

        contains

            procedure :: init => surface_mesh_init
            procedure :: output_results => surface_mesh_output_results
            procedure :: locate_kutta_edges => surface_mesh_locate_kutta_edges
            procedure :: initialize_wake => surface_mesh_initialize_wake

    end type surface_mesh


    type cart_volume_mesh

    end type cart_volume_mesh

    
contains


    subroutine surface_mesh_init(this, settings)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        character(len=:),allocatable :: extension
        real,dimension(3) :: p_min, p_max
        integer :: loc, i

        ! Get mesh file
        call json_get(settings, 'file', this%mesh_file)
        this%mesh_file = trim(this%mesh_file)
        write(*,*) "    Initializing surface mesh from file: ", this%mesh_file

        ! Determine the type of mesh file
        loc = index(this%mesh_file, '.')
        extension = this%mesh_file(loc:len(this%mesh_file))

        ! Load vtk
        if (extension .eq. '.vtk') then
            call load_surface_vtk(this%mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        end if

        ! Display mesh info
        write(*,*) "    Surface mesh has", this%N_verts, "vertices and", this%N_panels, "panels."

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
        write(*,*)
        write(*,'(a)',advance='no') "     Loading vertices into ADT..."
        do i=1,this%N_verts
            call this%vertex_tree%add(this%vertices(i))
        end do
        write(*,*) "Done."

        ! Store other settings for wake models
        call json_get(settings, 'wake_model.wake_shedding_angle', this%kutta_angle)
        call json_get(settings, 'wake_model.trefftz_distance', this%trefftz_distance)
        call json_get(settings, 'wake_model.N_panels', this%N_wake_panels_streamwise)
        this%C_kutta_angle = cos(this%kutta_angle*pi/180.0)
    
    end subroutine surface_mesh_init


    subroutine surface_mesh_output_results(this, body_file, wake_file)

        implicit none

        class(surface_mesh),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: body_file, wake_file

        ! Write out data for body
        call write_surface_vtk(body_file, this%vertices, this%panels)
        
        ! Write out data for wake
        call write_surface_vtk(wake_file, this%wake_vertices, this%wake_panels)
    
    end subroutine surface_mesh_output_results


    subroutine surface_mesh_locate_kutta_edges(this, freestream_flow)
        ! Locates wake-shedding edges on the mesh based on the flow conditions.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        integer :: i, j, m, n, N_shared_verts
        real,dimension(3) :: d
        integer,dimension(2) :: shared_verts
        type(list) :: kutta_edge_starts, kutta_edge_stops
        logical :: abutting
        real :: distance

        write(*,*)
        write(*,'(a)',advance='no') "     Locating wake-shedding edges..."

        ! Loop through each pair of panels
        this%N_kutta_edges = 0
        do i=1,this%N_panels
            do j=i+1,this%N_panels

                ! Initialize for this panel pair
                n_shared_verts = 0
                abutting = .false.

                ! Check if the panels are abutting
                abutting_loop: do m=1,this%panels(i)%N
                    do n=1,this%panels(j)%N

                        ! Get distance between vertices
                        d = this%panels(i)%get_vertex_loc(m)-this%panels(j)%get_vertex_loc(n)
                        distance = norm(d)

                        ! Check distance
                        if (distance < 1e-10) then

                            ! Previously found a shared vertex
                            if (n_shared_verts == 1) then
                                abutting = .true.
                                shared_verts(2) = this%panels(i)%get_vertex_index(m)
                                exit abutting_loop

                            ! First shared vertex
                            else
                                n_shared_verts = 1
                                shared_verts(1) = this%panels(i)%get_vertex_index(m)
                            end if

                        end if

                    end do
                end do abutting_loop

                if (abutting) then

                    ! Check angle between panels
                    if (inner(this%panels(i)%normal, this%panels(j)%normal) < this%C_kutta_angle) then

                        ! Check angle with freestream
                        if (inner(this%panels(i)%normal, freestream_flow%V_inf) > 0.0 .or. &
                            inner(this%panels(j)%normal, freestream_flow%V_inf) > 0.0) then

                            ! Update number of Kutta edges
                            this%N_kutta_edges = this%N_kutta_edges + 1

                            ! Store in starts and stops list
                            call kutta_edge_starts%append(shared_verts(1))
                            call kutta_edge_stops%append(shared_verts(2))

                            ! Set toggle for vertices
                            if (this%vertices(shared_verts(1))%on_kutta_edge) then

                                ! If it's already on one, then this means it's also in one
                                this%vertices(shared_verts(1))%in_kutta_edge = .true.

                            else

                                ! Add the first time
                                this%vertices(shared_verts(1))%on_kutta_edge = .true.
                                call this%kutta_vertices%append(shared_verts(1))
                                this%vertices(shared_verts(1))%index_in_kutta_vertices = this%kutta_vertices%len()

                            end if

                            ! Do the same for the other vertex
                            if (this%vertices(shared_verts(2))%on_kutta_edge) then

                                this%vertices(shared_verts(2))%in_kutta_edge = .true.

                            else

                                this%vertices(shared_verts(2))%on_kutta_edge = .true.
                                call this%kutta_vertices%append(shared_verts(2))
                                this%vertices(shared_verts(2))%index_in_kutta_vertices = this%kutta_vertices%len()

                            end if

                            ! Set toggle for panels
                            this%panels(i)%on_kutta_edge = .true.
                            this%panels(j)%on_kutta_edge = .true.

                        end if
                    end if
                end if

            end do
        end do

        ! Create array of Kutta edges
        allocate(this%kutta_edges(this%N_kutta_edges))
        do i=1,this%N_kutta_edges

            ! Get indices of starting and ending vertices
            call get_item(kutta_edge_starts, i, m)
            call get_item(kutta_edge_stops, i, n)

            ! Store
            call this%kutta_edges(i)%init(this%vertices(m), this%vertices(n), m, n)

        end do

        write(*,*) "Done. Found", this%N_kutta_edges, "wake-shedding edges."

    end subroutine surface_mesh_locate_kutta_edges


    subroutine surface_mesh_initialize_wake(this, freestream_flow)
        ! Creates wake panels and sets their initial shape. Handles vertex association.

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        real :: distance, vertex_separation
        real,dimension(3) :: loc, start
        integer :: i, j, ind, kutta_vert_ind, i_start, i_stop, i1, i2, i3, i4
        integer :: N_kutta_verts, N_wake_verts, N_wake_panels

        ! Initialize wake
        write(*,*)
        write(*,'(a)',advance='no') "     Initializing wake..."

        ! Determine sizes
        N_kutta_verts = this%kutta_vertices%len()
        this%N_wake_verts = N_kutta_verts*(this%N_wake_panels_streamwise+1)
        this%N_wake_panels = this%N_kutta_edges*this%N_wake_panels_streamwise

        ! Allocate storage
        allocate(this%wake_vertices(this%N_wake_verts))
        allocate(this%wake_panels(this%N_wake_panels))

        ! Determine vertex placement
        do i=1,N_kutta_verts

            ! Determine distance from origin to Kutta vertex in direction of the flow
            call get_item(this%kutta_vertices, i, kutta_vert_ind)
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
        write(*,*) this%N_kutta_edges
        do i=1,this%N_kutta_edges

            ! Determine which Kutta vertices this panel lies between
            i_start = this%kutta_edges(i)%v1%index_in_kutta_vertices
            i_stop = this%kutta_edges(i)%v2%index_in_kutta_vertices

            ! Create panels heading downstream
            do j=1,this%N_wake_panels_streamwise

                ! Determine index of panel
                ind = (i-1)*this%N_wake_panels_streamwise+j

                ! Determine vertex indices
                i1 = (i_start-1)*(this%N_wake_panels_streamwise+1)+j
                i2 = (i_start-1)*(this%N_wake_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(this%N_wake_panels_streamwise+1)+j+1
                i4 = (i_stop-1)*(this%N_wake_panels_streamwise+1)+j

                ! Initialize
                call this%wake_panels(ind)%init(this%wake_vertices(i1),&
                                                this%wake_vertices(i2),&
                                                this%wake_vertices(i3),&
                                                this%wake_vertices(i4))

            end do
        end do

        write(*,*) "Done."

    end subroutine surface_mesh_initialize_wake

end module surface_mesh_mod