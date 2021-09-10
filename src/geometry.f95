! Types for geometric objects
module geometry

    implicit none


    type vertex

        real,dimension(3) :: loc, cp ! Location and associated control point
        logical :: on_kutta_edge
        integer,dimension(20) :: neighboring_panels = -1
        integer :: index ! Index of this vertex in the mesh

        contains

            procedure :: init => vertex_init

    end type vertex


    type vertex_pointer

        type(vertex),pointer :: ptr

    end type vertex_pointer


    type base_panel

        integer :: N ! Number of sides/vertices
        type(vertex_pointer),dimension(:),allocatable :: vertices
        real,dimension(3) :: n_hat ! Normal vector
        real :: A ! Surface area

    end type base_panel


    type,extends(base_panel) :: tri_panel

        contains

            procedure :: init => tri_panel_init

    end type tri_panel


    type,extends(base_panel) :: quad_panel

        contains

            procedure :: init => quad_panel_init

    end type quad_panel


    type panel_pointer

        class(base_panel),pointer :: ptr

    end type panel_pointer

    
contains


    subroutine vertex_init(this, loc, index)

        implicit none

        class(vertex),intent(inout) :: this
        real,dimension(3),intent(in) :: loc
        integer,intent(in) :: index

        ! Store info
        this%loc = loc
        this%index = index

    end subroutine vertex_init


    subroutine tri_panel_init(this, v1, v2, v3)

        implicit none

        class(tri_panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3

        ! Set number of sides
        this%N = 3

        ! Allocate vertex array
        allocate(this%vertices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3

        ! Calculate normal vec

        ! Calculate area

    end subroutine tri_panel_init


    subroutine quad_panel_init(this, v1, v2, v3, v4)

        implicit none

        class(quad_panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3, v4
        
        ! Set number of sides
        this%N = 4

        ! Allocate vertex array
        allocate(this%vertices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3
        this%vertices(4)%ptr => v4

        ! Calculate normal vector

        ! Calculate area

    end subroutine quad_panel_init
    
end module geometry