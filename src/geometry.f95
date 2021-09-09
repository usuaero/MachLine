! Types for geometric objects
module geometry

    implicit none


    type panel

        integer :: N, i1, i2, i3 ! Number of sides (will always be 3) and indices of each vertex
        real,dimension(3) :: v1, v2, v3 ! Vertices

    end type panel

    type vertex

        real,dimension(3) :: p, cp ! Location and associated control point
        logical :: on_kutta_edge
        integer,dimension(20) :: neighboring_panels = -1

    end type vertex
    
    contains
    
end module geometry