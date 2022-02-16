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
            procedure :: update => wake_mesh_update

    end type wake_mesh


contains


    subroutine wake_mesh_init(this, freestream, top_edge_verts, bot_edge_verts, mesh_edges, wake_edge_indices, &
                              N_panels_streamwise, mesh_vertices, trefftz_distance, mirrored_and_asym, mirror_plane)
        ! Creates the vertices and panels. Handles vertex association.

        implicit none

        class(wake_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream
        integer,allocatable,dimension(:),intent(in) :: top_edge_verts, bot_edge_verts
        type(edge),allocatable,dimension(:),intent(in) :: mesh_edges
        integer,allocatable,dimension(:),intent(in) :: wake_edge_indices
        integer,intent(in) :: N_panels_streamwise
        type(vertex),allocatable,dimension(:),intent(in) :: mesh_vertices
        real,intent(in) :: trefftz_distance
        logical,intent(in) :: mirrored_and_asym
        integer,intent(in) :: mirror_plane

        real :: distance, vertex_separation, mirrored_distance, mirrored_vertex_separation
        real,dimension(3) :: loc, start, mirrored_start
        integer :: i, j, k, i_vert, i_top_parent, i_bot_parent, i_start, i_stop, i1, i2, i3, i4
        integer :: N_wake_edge_verts, N_wake_edges

        write(*,'(a e10.4 a)',advance='no') "     Initializing wake with a Trefftz distance of ", trefftz_distance, "..."

        ! Determine necessary number of vertices
        N_wake_edge_verts = size(top_edge_verts)
        if (mirrored_and_asym) then
            this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)*2
        else
            this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)
        end if

        ! Allocate vertex storage
        allocate(this%vertices(this%N_verts))

        ! Determine vertex placement
        do i=1,N_wake_edge_verts

            ! Get indices
            i_top_parent = top_edge_verts(i)
            i_bot_parent = bot_edge_verts(i)

            ! Determine distance from origin to wake-shedding vertex in the direction of the freestream flow
            start = mesh_vertices(i_top_parent)%loc
            distance = trefftz_distance-inner(start, freestream%c_hat_g)

            ! Double-check
            if (dist(start, mesh_vertices(i_bot_parent)%loc) > 1.e-12) then
                write(*,*) "!!! Wake edge vertices are not identical. Quitting..."
                stop
            end if

            ! Determine vertex separation
            vertex_separation = distance/N_panels_streamwise

            ! Same for mirror
            if (mirrored_and_asym) then

                ! Determine start location
                mirrored_start = mirror_about_plane(start, mirror_plane)
                mirrored_distance = trefftz_distance-inner(mirrored_start, freestream%c_hat_g)

                ! Determine vertex separation
                mirrored_vertex_separation = mirrored_distance/N_panels_streamwise

            end if

            ! Place vertices
            do j=1,N_panels_streamwise+1

                ! Determine location
                i_vert = (i-1)*(N_panels_streamwise+1)+j
                loc = start + vertex_separation*(j-1)*freestream%c_hat_g

                ! Initialize vertex
                call this%vertices(i_vert)%init(loc, i_vert)

                ! Set parent index
                this%vertices(i_vert)%top_parent = i_top_parent
                this%vertices(i_vert)%bot_parent = i_bot_parent

                ! Initialize mirror
                if (mirrored_and_asym) then

                    ! Determine location
                    i_vert = i_vert + this%N_verts/2
                    mirrored_start = mirror_about_plane(start, mirror_plane)
                    loc = mirrored_start + mirrored_vertex_separation*(j-1)*freestream%c_hat_g

                    ! Initialize vertex
                    call this%vertices(i_vert)%init(loc, i_vert)
                    
                    ! Set parent index
                    this%vertices(i_vert)%top_parent = i_top_parent + size(mesh_vertices)
                    this%vertices(i_vert)%bot_parent = i_bot_parent + size(mesh_vertices)

                end if

            end do
        end do

        ! Determine necessary number of panels
        N_wake_edges = size(wake_edge_indices)
        if (mirrored_and_asym) then
            this%N_panels = N_wake_edges*N_panels_streamwise*4
        else
            this%N_panels = N_wake_edges*N_panels_streamwise*2
        end if
        allocate(this%panels(this%N_panels))

        ! Initialize panels
        do k=1,N_wake_edges

            i = wake_edge_indices(k)

            ! Determine which wake-shedding vertices this panel lies between
            i_start = mesh_vertices(mesh_edges(i)%verts(1))%index_in_wake_vertices
            i_stop = mesh_vertices(mesh_edges(i)%verts(2))%index_in_wake_vertices

            ! Create panels heading downstream
            do j=1,N_panels_streamwise

                ! Determine index of first triangular panel
                i_vert = (k-1)*N_panels_streamwise*2+2*j-1

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_start-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1

                ! Initialize
                call this%panels(i_vert)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_vert)

                ! Specify this panel is in the wake
                this%panels(i_vert)%in_wake = .true.

                ! Create mirror
                if (mirrored_and_asym) then

                    ! Determine index
                    i_vert = i_vert + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_start-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2

                    ! Initialize (order of vertices is reversed to maintain panel orientation through mirror)
                    call this%panels(i_vert)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_vert)

                    ! Specify this panel is in the wake
                    this%panels(i_vert)%in_wake = .true.

                end if

                ! Determine index of second triangular panel
                i_vert = (k-1)*N_panels_streamwise*2+2*j

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j

                ! Initialize
                call this%panels(i_vert)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_vert)

                ! Specify this panel is in the wake
                this%panels(i_vert)%in_wake = .true.

                ! Create mirror
                if (mirrored_and_asym) then

                    ! Determine index
                    i_vert = i_vert + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+this%N_verts/2

                    ! Initialize (again, order is reversed)
                    call this%panels(i_vert)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_vert)

                    ! Specify this panel is in the wake
                    this%panels(i_vert)%in_wake = .true.

                end if

            end do
        end do

        ! Intialize panel coordinate transforms
        do i=1,this%N_panels
            call this%panels(i)%calc_transforms(freestream)
        end do

        write(*,'(a, i7, a, i7, a)') "Done. Created ", this%N_verts, " wake vertices and ", this%N_panels, " wake panels."

    end subroutine wake_mesh_init


    subroutine wake_mesh_update(this)

        implicit none

        class(wake_mesh),intent(inout) :: this
        integer :: i

        ! Update panel properties since vertices were moved
        do i=1,this%N_panels
            call this%panels(i)%calc_derived_properties
        end do
    
    end subroutine wake_mesh_update


end module wake_mesh_mod