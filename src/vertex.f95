! Types for geometric objects
module vertex_mod

    use linked_list_mod
    use math_mod

    implicit none


    type vertex
        ! A vertex in 3-space

        real,dimension(3) :: loc, cp ! Location and associated control point
        logical :: on_wake_edge ! Whether this vertex is on a wake-shedding edge
        logical :: in_wake_edge ! Whether this vertex is not an endpoint for a chain of wake-shedding edges (i.e. has two attached edges)
        integer :: index ! Index of this vertex in the mesh
        integer :: index_in_wake_vertices ! Index of this vertex in the list of wake-shedding vertices
        integer :: top_parent ! Index of the top vertex this vertex's strength is determined by (for a wake vertex)
        integer :: bot_parent ! Index of the bottom vertex this vertex's strength is determined by (for a wake vertex)
        type(list) :: adjacent_vertices ! List of vertices which share a panel with this vertex
        real :: l_avg ! Average of the edge lengths adjacent to this vertex
        type(list) :: panels ! List of indices for the panels which connect to this vertex
        type(list) :: panels_not_across_wake_edge ! List of indices for the panels which connect to this vertex not across a wake-shedding edge
        real,dimension(3) :: normal ! Normal vector associated with this control point

        contains

            procedure :: init => vertex_init
            procedure :: calc_average_edge_length => vertex_calc_average_edge_length

    end type vertex


    type vertex_pointer
        ! A pointer to a vertex, for creating vertex arrays

        type(vertex),pointer :: ptr

    end type vertex_pointer

    
contains


    subroutine vertex_init(this, loc, index)
        ! Initializes a vertex

        implicit none

        class(vertex),intent(inout) :: this
        real,dimension(3),intent(in) :: loc
        integer,intent(in) :: index

        ! Store info
        this%loc = loc
        this%index = index

        ! Intitialize some data
        this%top_parent = 0
        this%bot_parent = 0

    end subroutine vertex_init


    subroutine vertex_calc_average_edge_length(this, vertices)
        ! Calculates the average edge length of edges adjacent to this vertex

        implicit none

        class(vertex),intent(inout) :: this
        type(vertex),dimension(:),allocatable,intent(in) :: vertices

        integer :: i, ind

        ! Loop through adjacent vertices
        this%l_avg = 0.
        do i=1,this%adjacent_vertices%len()
            
            ! Get index of adjacent vertex
            call this%adjacent_vertices%get(i, ind)

            ! Add length
            this%l_avg = this%l_avg + dist(this%loc, vertices(ind)%loc)

        end do
        
        ! Compute average
        this%l_avg = this%l_avg/this%adjacent_vertices%len()
    
    end subroutine vertex_calc_average_edge_length

    
end module vertex_mod