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
        real,dimension(3) :: centroid
        real,dimension(3,3) :: A_t ! Local coordinate transform matrix
        real,dimension(:,:),allocatable :: vertices_local ! Location of the vertices described in local coords
        real :: A ! Surface area
        real :: phi_n = 0 ! Perturbation source strength
        logical :: on_kutta_edge ! Whether this panel belongs to a Kutta edge
        logical :: wake_panel ! Whether this panel belongs to a wake
        logical :: shock_panel ! Whether this panel belongs to a shock
        integer,dimension(:),allocatable :: vertex_indices ! Indices of this panel's vertices in the mesh vertex array
        type(list) :: opposing_kutta_panels ! Indices of panels opposite this one on the Kutta edge(s)
        type(list) :: abutting_panels ! Indices of panels abutting this one (not across Kutta edge)

        contains

            procedure :: panel_init_3
            procedure :: panel_init_4
            generic :: init => panel_init_3, panel_init_4
            procedure :: calc_derived_properties =>panel_calc_derived_properties
            procedure :: calc_area => panel_calc_area
            procedure :: calc_normal => panel_calc_normal
            procedure :: calc_centroid => panel_calc_centroid
            procedure :: calc_coord_transform => panel_calc_coord_transform
            procedure :: get_vertex_loc => panel_get_vertex_loc
            procedure :: get_vertex_index => panel_get_vertex_index
            procedure :: touches_vertex => panel_touches_vertex
            procedure :: point_to_vertex_clone => panel_point_to_vertex_clone

    end type panel

    
contains


    subroutine panel_init_3(this, v1, v2, v3, i1, i2, i3)
        ! Initializes a 3-panel

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3
        integer,intent(in) :: i1, i2, i3

        ! Set number of sides
        this%N = 3

        ! Allocate vertex array
        allocate(this%vertices(this%N))
        allocate(this%vertices_local(this%N,3))
        allocate(this%vertex_indices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3
        this%vertex_indices(1) = i1
        this%vertex_indices(2) = i2
        this%vertex_indices(3) = i3

        call this%calc_derived_properties()

    end subroutine panel_init_3


    subroutine panel_init_4(this, v1, v2, v3, v4, i1, i2, i3, i4)
        ! Initializes a panel with 4 sides

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3, v4
        integer,intent(in) :: i1, i2, i3, i4
        
        ! Set number of sides
        this%N = 4

        ! Allocate vertex array
        allocate(this%vertices(this%N))
        allocate(this%vertices_local(this%N,3))
        allocate(this%vertex_indices(this%N))

        ! Store info
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3
        this%vertices(4)%ptr => v4
        this%vertex_indices(1) = i1
        this%vertex_indices(2) = i2
        this%vertex_indices(3) = i3
        this%vertex_indices(4) = i4

        call this%calc_derived_properties()

    end subroutine panel_init_4


    subroutine panel_calc_derived_properties(this)
        ! Initializes properties based on the location of the vertices.
        ! Should be called when panel geometry is updated.

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

        ! Calculate centroid
        call this%calc_centroid()

        ! Calculate coordinate transform
        call this%calc_coord_transform()

    end subroutine panel_calc_derived_properties


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


    subroutine panel_calc_centroid(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: sum = 0
        integer :: i

        ! Get average of corner points
        do i=1,this%N
            sum = sum + this%get_vertex_loc(i)
        end do

        ! Set centroid
        this%centroid = sum/this%N

    end subroutine panel_calc_centroid


    subroutine panel_calc_coord_transform(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: d
        integer :: i

        ! Choose first edge tangent as local xi-axis
        ! (will need to be projected for quadrilateral panel)
        d = this%get_vertex_loc(2)-this%get_vertex_loc(1)
        this%A_t(1,:) = d/norm(d)

        ! Panel normal is the zeta axis
        this%A_t(3,:) = this%normal

        ! Calculate eta axis from the other two
        this%A_t(2,:) = cross(this%normal, this%A_t(1,:))

        ! Transform vertex coords
        do i=1,this%N
            this%vertices_local(i,:) = matmul(this%A_t, this%get_vertex_loc(i))
        end do

    end subroutine panel_calc_coord_transform


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


    function panel_touches_vertex(this, i) result(touches)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        logical :: touches
        integer :: j

        touches = .false.

        ! Loop through vertices
        do j=1,this%N

            ! Check index
            if (this%vertex_indices(j) == i) then
                touches = .true.
                return
            end if

        end do

    end function panel_touches_vertex


    subroutine panel_point_to_vertex_clone(this, clone)
        ! Updates the panel to point to this new vertex (assumed to be clone of a current vertex)

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: clone
        integer :: i

        ! Loop through vertices
        do i=1,this%N

            ! Check which vertex this will replace
            if (norm(this%vertices(i)%ptr%loc-clone%loc) < 1e-10) then

                ! Update pointer
                this%vertices(i)%ptr => clone

                ! Update index
                this%vertex_indices(i) = clone%index

                return

            end if

        end do
    
    end subroutine panel_point_to_vertex_clone
    
end module panel_mod