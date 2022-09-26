module node_mod

    use vertex_mod

    implicit none

    type :: node

        integer :: index
        type(vertex),pointer :: vert
        logical :: mirror_of_vert

        contains

            procedure :: init => node_init
            procedure :: get_loc => node_get_loc
    
    end type node
    
contains

    subroutine node_init(this, i, v, mirror)
        ! Initializes the nod

        implicit none
        
        class(node),intent(inout) :: this
        integer,intent(in) :: i
        type(vertex),intent(in),target :: v
        logical,intent(in) :: mirror
    
        ! Store
        this%index = i
        this%vert => v
        this%mirror_of_vert = mirror
        
    end subroutine node_init


    function node_get_loc(this) result(loc)
        ! Returns the location of the vertex this node is tied to

        implicit none
        
        class(vertex),intent(in) :: this

        real,dimension(3) :: loc

        ! Get location
        loc = this%vertex%loc
        
    end function node_get_loc

    
end module node_mod