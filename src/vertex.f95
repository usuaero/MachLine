! Types for geometric objects
module vertex_mod

    use linked_list_mod

    implicit none


    type vertex
        ! A vertex in 3-space

        real,dimension(3) :: loc, cp ! Location and associated control point
        logical :: on_kutta_edge ! Whether this vertex is on a Kutta edge
        logical :: in_kutta_edge ! Whether this vertex is not an endpoint for a chain of Kutta edges (i.e. has two attached edges)
        integer :: index ! Index of this vertex in the mesh
        integer :: index_in_kutta_vertices ! Index of this vertex in the list of Kutta vertices
        integer :: parent = 0 ! Index of the vertex this vertex's strength is determined by (for a wake vertex)
        type(list) :: panels ! List of indices for the panels which connect to this vertex
        type(list) :: panels_not_across_kutta_edge ! List of indices for the panels which connect to this vertex not across a Kutta edge
        real :: phi = 0 ! Perturbation doublet strength
        real,dimension(3) :: normal ! Normal vector associated with this control point

        contains

            procedure :: init => vertex_init

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

    end subroutine vertex_init

    
end module vertex_mod