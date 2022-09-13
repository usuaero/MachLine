! Types for geometric objects
module vertex_mod

    use linked_list_mod
    use math_mod

    implicit none


    type vertex
        ! A vertex in 3-space

        integer :: vert_type ! Whether this is a 1) true vertex or 2) vertex representing an edge midpoint
        real,dimension(3) :: loc ! Location
        real,dimension(3) :: n_g, n_g_mir ! Normal vector associated with this control point
        real :: l_avg ! Average of the edge lengths adjacent to this vertex
        real :: l_min ! Minimum of the edge lengths adjacent to this vertex
        type(list) :: adjacent_vertices ! List of indices for the vertices which share an edge with this vertex
        type(list) :: adjacent_edges ! List of indices for the edges which touch this vertex
        type(list) :: panels ! List of indices for the panels which connect to this vertex
        type(list) :: panels_not_across_wake_edge ! List of indices for the panels which connect to this vertex not across a wake-shedding edge
        integer :: N_wake_edges ! Number of wake edges this vertex belongs to
        integer :: index ! Index of this vertex in the mesh
        integer :: index_in_wake_vertices ! Index of this vertex in the list of wake-shedding vertices
        integer :: top_parent, bot_parent ! Indices of the top and bottom vertices this vertex's strength is determined by (for a wake vertex)
        logical :: on_mirror_plane ! Whether this vertex lies in the mirroring plane
        logical :: clone ! Whether this vertex needs a clone depending on whether it's in a wake-shedding edge
        logical :: mirrored_is_unique ! Whether this vertice's mirror image will be the same for an asymmetric freestream condition
        integer :: i_wake_partner ! Index of the vertex, which along with this one, will determine wake strength

        contains

            procedure :: init => vertex_init

            ! Initializer-setters
            procedure :: set_whether_on_mirror_plane => vertex_set_whether_on_mirror_plane
            procedure :: set_average_edge_length => vertex_set_average_edge_length

            ! Getters
            procedure :: get_N_needed_clones => vertex_get_N_needed_clones

    end type vertex


    type vertex_pointer
        ! A pointer to a vertex, for creating vertex arrays

        type(vertex),pointer :: ptr

    end type vertex_pointer

    
contains


    subroutine vertex_init(this, loc, index, vert_type)
        ! Initializes a vertex

        implicit none

        class(vertex),intent(inout) :: this
        real,dimension(3),intent(in) :: loc
        integer,intent(in) :: index, vert_type

        ! Store info
        this%loc = loc
        this%index = index
        this%vert_type = vert_type

        ! Intitialize some data
        this%top_parent = 0
        this%bot_parent = 0

        ! Default cases
        this%mirrored_is_unique = .true.
        this%clone = .false.
        this%on_mirror_plane = .false.
        this%i_wake_partner = index

    end subroutine vertex_init


    subroutine vertex_set_average_edge_length(this, vertices)
        ! Calculates the average edge length of edges adjacent to this vertex

        implicit none

        class(vertex),intent(inout) :: this
        type(vertex),dimension(:),allocatable,intent(in) :: vertices

        integer :: i, adj_ind, N
        real :: l_i

        ! Loop through adjacent vertices
        this%l_avg = 0.
        this%l_min = huge(this%l_min)
        N = 0
        do i=1,this%adjacent_vertices%len()
            
            ! Get index of adjacent vertex
            call this%adjacent_vertices%get(i, adj_ind)

            ! Calculate edge length
            l_i = dist(this%loc, vertices(adj_ind)%loc)

            ! Get minimum
            this%l_min = min(this%l_min, l_i)

            ! For a vertex on the mirror plane where the adjacent vertex is not on the mirror plane
            ! that length will need to be added twice
            if (this%on_mirror_plane .and. .not. vertices(adj_ind)%on_mirror_plane) then

                ! Add twice
                this%l_avg = this%l_avg + 2*l_i
                N = N + 2
                
            else

                ! Add once
                this%l_avg = this%l_avg + l_i
                N = N + 1

            end if

        end do
        
        ! Compute average
        if (N > 0) then
            this%l_avg = this%l_avg/N
        else
            this%l_avg = 1.
        end if
    
    end subroutine vertex_set_average_edge_length


    subroutine vertex_set_whether_on_mirror_plane(this, mirror_plane)
        ! Sets the member variable telling whether this vertex is on the mirror plane

        implicit none
        
        class(vertex), intent(inout) :: this
        integer, intent(in) :: mirror_plane

        ! Check distance from mirror plane
        if (abs(this%loc(mirror_plane))<1e-12) then

            ! The vertex is on the mirror plane
            this%on_mirror_plane = .true.

        end if
    
        
    end subroutine vertex_set_whether_on_mirror_plane


    function vertex_get_N_needed_clones(this) result(N_clones)
        ! Returns the number of clones this vertex needs

        implicit none
        
        class(vertex),intent(in) :: this

        integer :: N_clones

        ! Check if clones are needed at all
        if (this%clone) then

            ! Regular vertices need one less than the number of adjacent wake edges
            if (this%vert_type == 1) then
                N_clones = this%N_wake_edges - 1

            ! Midpoints only ever need one
            else
                N_clones = 1
            end if

        else
            N_clones = 0
        end if
        
    end function vertex_get_N_needed_clones

    
end module vertex_mod