! Panel type
module panel_mod

    use linked_list_mod
    use vertex_mod
    use math_mod

    implicit none


    type panel
        ! A panel with an arbitrary number of sides

        integer :: N ! Number of sides/vertices
        type(vertex_pointer),dimension(:),allocatable :: vertices
        real,dimension(3) :: normal ! Normal vector
        real,dimension(:,:),allocatable :: midpoints
        real :: A ! Surface area
        real :: phi_n = 0 ! Perturbation source strength
        logical :: on_kutta_edge
        logical :: wake_panel ! Whether this panel belongs to a wake
        logical :: shock_panel ! Whether this panel belongs to a shock

        contains

            procedure :: panel_init_3
            procedure :: panel_init_4
            generic :: init => panel_init_3, panel_init_4
            procedure :: init_common =>panel_init_common
            procedure :: calc_area => panel_calc_area
            procedure :: calc_normal => panel_calc_normal
            procedure :: get_vertex_loc => panel_get_vertex_loc
            procedure :: get_vertex_index => panel_get_vertex_index

    end type panel

    
contains


    subroutine panel_init_3(this, v1, v2, v3)
        ! Initializes a 3-panel

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3

        ! Set number of sides
        this%N = 3

        ! Allocate vertex array
        allocate(this%vertices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3

        call this%init_common()

    end subroutine panel_init_3


    subroutine panel_init_4(this, v1, v2, v3, v4)
        ! Initializes a panel with 4 sides

        implicit none

        class(panel),intent(inout) :: this
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

        call this%init_common()

    end subroutine panel_init_4


    subroutine panel_init_common(this)
        ! Initializes attributes common to both types of panels

        implicit none

        class(panel),intent(inout) :: this
        integer :: i

        ! Determine midpoints
        allocate(this%midpoints(this%N,3))
        do i=1,this%N-1
            this%midpoints(i,:) = 0.5*(this%get_vertex_loc(i)+this%get_vertex_loc(i+1))
        end do
        this%midpoints(this%N,:) = 0.5*(this%get_vertex_loc(1)+this%get_vertex_loc(this%N))

        ! Calculate normal vec
        call this%calc_normal()

        ! Calculate area
        call this%calc_area()
    
    end subroutine panel_init_common


    subroutine panel_calc_area(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: d1, d2, d3, d4

        ! 3-sided panel
        if (this%N == 3) then

            ! Get side vectors
            d1 = this%get_vertex_loc(2)-this%get_vertex_loc(1)
            d2 = this%get_vertex_loc(3)-this%get_vertex_loc(2)

            ! Calculate area from cross product
            this%A = 0.5*norm(cross(d1, d2))

        ! 4-sided panel
        else

            ! Get side vectors
            d1 = this%get_vertex_loc(2)-this%get_vertex_loc(1)
            d2 = this%get_vertex_loc(3)-this%get_vertex_loc(2)
            d3 = this%get_vertex_loc(4)-this%get_vertex_loc(3)
            d4 = this%get_vertex_loc(1)-this%get_vertex_loc(4)

            ! Calculate area from cross product
            this%A = 0.5*(norm(cross(d1, d2))+norm(cross(d3, d4)))

        end if

    end subroutine panel_calc_area


    subroutine panel_calc_normal(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: d1, d2

        ! Get two chord vectors from midpoints
        d1 = this%midpoints(2,:)-this%midpoints(1,:)
        d2 = this%midpoints(3,:)-this%midpoints(2,:)

        ! Find normal
        this%normal = cross(d1, d2)
        this%normal = this%normal/norm(this%normal)

    end subroutine panel_calc_normal


    function panel_get_vertex_loc(this, i) result(loc)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        real,dimension(3) :: loc

        loc = this%vertices(i)%ptr%loc

    end function panel_get_vertex_loc


    function panel_get_vertex_index(this, i) result(index)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        integer :: index

        index = this%vertices(i)%ptr%index

    end function panel_get_vertex_index
    
end module panel_mod