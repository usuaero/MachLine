module kutta_edge_mod

    use vertex_mod

    implicit none

    type kutta_edge

        type(vertex),pointer :: v1, v2
        integer :: i1, i2 ! Indices of the end vertices in the mesh vertex array
        integer :: i1_kutta_vert, i2_kutta_vert ! Indices of the end vertices in the Kutta vertex list

        contains

            procedure :: init => kutta_edge_init

    end type kutta_edge

    
contains


    subroutine kutta_edge_init(this, v1, v2, i1, i2)

        implicit none

        class(kutta_edge),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2
        integer,intent(in) :: i1, i2

        ! Point to vertices
        this%v1 => v1
        this%v2 => v2

        ! Store indices
        this%i1 = i1
        this%i2 = i2
    
    end subroutine kutta_edge_init
    
end module kutta_edge_mod