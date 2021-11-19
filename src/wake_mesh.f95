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
    use wake_edge_mod

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


    subroutine wake_mesh_init(this, freestream_flow, top_edge_verts, bot_edge_verts, wake_edges, &
                              N_panels_streamwise, mesh_vertices, trefftz_distance, mirrored_and_asym, mirror_plane)
        ! Creates the vertices and panels. Handles vertex association.

        implicit none

        class(wake_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        integer,allocatable,dimension(:),intent(in) :: top_edge_verts, bot_edge_verts
        type(wake_edge),allocatable,dimension(:),intent(in) :: wake_edges
        integer,intent(in) :: N_panels_streamwise
        type(vertex),allocatable,dimension(:),intent(in) :: mesh_vertices
        real,intent(in) :: trefftz_distance
        logical,intent(in) :: mirrored_and_asym
        integer,intent(in) :: mirror_plane

        real :: distance, vertex_separation, mirrored_distance, mirrored_vertex_separation
        real,dimension(3) :: loc, start, mirrored_start
        integer :: i, j, ind, top_parent_ind, bot_parent_ind, i_start, i_stop, i1, i2, i3, i4
        integer :: N_wake_edge_verts, N_wakes

        write(*,*)
        write(*,'(a)',advance='no') "     Initializing wake..."

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

            ! Determine distance from origin to wake-shedding vertex in the direction of the freestream flow
            top_parent_ind = top_edge_verts(i)
            bot_parent_ind = bot_edge_verts(i)
            start = mesh_vertices(top_parent_ind)%loc
            distance = trefftz_distance-inner(start, freestream_flow%c0)

            ! Determine vertex separation
            vertex_separation = distance/N_panels_streamwise

            ! Same for mirror
            if (mirrored_and_asym) then

                ! Determine start location
                mirrored_start = mirror_about_plane(start, mirror_plane)
                mirrored_distance = trefftz_distance-inner(mirrored_start, freestream_flow%c0)

                ! Determine vertex separation
                mirrored_vertex_separation = mirrored_distance/N_panels_streamwise

            end if

            ! Place vertices
            do j=1,N_panels_streamwise+1

                ! Determine location
                ind = (i-1)*(N_panels_streamwise+1)+j
                loc = start + vertex_separation*(j-1)*freestream_flow%c0

                ! Initialize vertex
                call this%vertices(ind)%init(loc, ind)

                ! Set parent index
                this%vertices(ind)%top_parent = top_parent_ind
                this%vertices(ind)%bot_parent = bot_parent_ind

                ! Initialize mirror
                if (mirrored_and_asym) then

                    ! Determine location
                    ind = ind + this%N_verts/2
                    mirrored_start = mirror_about_plane(start, mirror_plane)
                    loc = mirrored_start + mirrored_vertex_separation*(j-1)*freestream_flow%c0

                    ! Initialize vertex
                    call this%vertices(ind)%init(loc, ind)
                    
                    ! Set parent index
                    this%vertices(ind)%top_parent = top_parent_ind + size(mesh_vertices)
                    this%vertices(ind)%bot_parent = bot_parent_ind + size(mesh_vertices)

                end if

            end do
        end do

        ! Determine necessary number of panels
        if (mirrored_and_asym) then
            this%N_panels = size(wake_edges)*N_panels_streamwise*4
        else
            this%N_panels = size(wake_edges)*N_panels_streamwise*2
        end if
        allocate(this%panels(this%N_panels))

        ! Initialize panels
        do i=1,size(wake_edges)

            ! Determine which wake-shedding vertices this panel lies between
            i_start = mesh_vertices(wake_edges(i)%i1)%index_in_wake_vertices
            i_stop = mesh_vertices(wake_edges(i)%i2)%index_in_wake_vertices

            ! Create panels heading downstream
            do j=1,N_panels_streamwise

                ! Determine index of first triangular panel
                ind = (i-1)*N_panels_streamwise*2+2*j-1

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_start-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1

                ! Initialize
                call this%panels(ind)%init(this%vertices(i1),&
                                           this%vertices(i2),&
                                           this%vertices(i3),&
                                           i1, i2, i3, ind)

                ! Specify this panel is in the wake
                this%panels(ind)%in_wake = .true.

                ! Create mirror
                if (mirrored_and_asym) then

                    ! Determine index
                    ind = ind + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_start-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2

                    ! Initialize (order of vertices is reversed to maintain panel orientation through mirror)
                    call this%panels(ind)%init(this%vertices(i3),&
                                               this%vertices(i2),&
                                               this%vertices(i1),&
                                               i3, i2, i1, ind)

                    ! Specify this panel is in the wake
                    this%panels(ind)%in_wake = .true.

                end if

                ! Determine index of second triangular panel
                ind = (i-1)*N_panels_streamwise*2+2*j

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j

                ! Initialize
                call this%panels(ind)%init(this%vertices(i1),&
                                           this%vertices(i2),&
                                           this%vertices(i3),&
                                           i1, i2, i3, ind)

                ! Specify this panel is in the wake
                this%panels(ind)%in_wake = .true.

                ! Create mirror
                if (mirrored_and_asym) then

                    ! Determine index
                    ind = ind + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+this%N_verts/2

                    ! Initialize (again, order is reversed)
                    call this%panels(ind)%init(this%vertices(i3),&
                                               this%vertices(i2),&
                                               this%vertices(i1),&
                                               i3, i2, i1, ind)

                    ! Specify this panel is in the wake
                    this%panels(ind)%in_wake = .true.

                end if

            end do
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