! Class for modeling wake meshes
module wake_mesh_mod

    use json_mod
    use json_xtnsn_mod
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


    subroutine wake_mesh_init(this, freestream_flow, top_edge_verts, bot_edge_verts,&
                              wake_edges, N_panels_streamwise, mesh_vertices, trefftz_distance)
        ! Creates the vertices and panels. Handles vertex association.

        implicit none

        class(wake_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        integer,allocatable,dimension(:),intent(in) :: top_edge_verts, bot_edge_verts
        type(wake_edge),allocatable,dimension(:),intent(in) :: wake_edges
        integer,intent(in) :: N_panels_streamwise
        type(vertex),allocatable,dimension(:),intent(in) :: mesh_vertices
        real,intent(in) :: trefftz_distance

        real :: distance, vertex_separation
        real,dimension(3) :: loc, start
        integer :: i, j, ind, top_parent_ind, bot_parent_ind, i_start, i_stop, i1, i2, i3, i4
        integer :: N_wake_edge_verts

        write(*,*)
        write(*,'(a)',advance='no') "     Initializing wake..."

        ! Determine sizes
        N_wake_edge_verts = size(top_edge_verts)
        this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)
        this%N_panels = size(wake_edges)*N_panels_streamwise*2

        ! Allocate storage
        allocate(this%vertices(this%N_verts))
        allocate(this%panels(this%N_panels))

        ! Determine vertex placement
        do i=1,N_wake_edge_verts

            ! Determine distance from origin to wake-shedding vertex in the direction of the freestream flow
            top_parent_ind = top_edge_verts(i)
            bot_parent_ind = bot_edge_verts(i)
            start = mesh_vertices(top_parent_ind)%loc
            distance = trefftz_distance-inner(start, freestream_flow%u_inf)

            ! Determine vertex separation
            vertex_separation = distance/N_panels_streamwise

            ! Place vertices
            do j=1,N_panels_streamwise+1

                ! Determine location
                ind = (i-1)*(N_panels_streamwise+1)+j
                loc = start+vertex_separation*(j-1)*freestream_flow%u_inf

                ! Initialize vertex
                call this%vertices(ind)%init(loc, ind)

                ! Set parent index
                this%vertices(ind)%top_parent = top_parent_ind
                this%vertices(ind)%bot_parent = bot_parent_ind

            end do
        end do

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

            end do
        end do

        write(*,*) "Done. Created", this%N_verts, "wake vertices and", this%N_panels, "wake panels."

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