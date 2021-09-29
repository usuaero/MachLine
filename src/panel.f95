! Panel type
module panel_mod

    use linked_list_mod
    use vertex_mod
    use math_mod

    implicit none


    type panel
        ! A panel with an arbitrary number of sides

        integer :: N = 3 ! Number of sides/vertices
        integer :: index ! Index of this panel in the mesh array
        type(vertex_pointer),dimension(:),allocatable :: vertices
        real,dimension(3) :: normal ! Normal vector
        real,dimension(:,:),allocatable :: midpoints
        real,dimension(3) :: centroid
        real,dimension(3,3) :: A_t ! Local coordinate transform matrix
        real,dimension(:,:),allocatable :: vertices_local ! Location of the vertices described in local coords
        real,dimension(:,:),allocatable :: t_hat, t_hat_local ! Edge unit tangents
        real,dimension(:,:),allocatable :: n_hat, n_hat_local ! Edge unit outward normals
        real,dimension(:),allocatable :: l ! Edge lengths
        real :: A ! Surface area
        real :: phi_n = 0 ! Perturbation source strength
        real :: mu_0_1, mu_x_1, mu_y_1 ! Influence of vertex 1 on the doublet integral parameters
        real :: mu_0_2, mu_x_2, mu_y_2 ! Influence of vertex 2 on the doublet integral parameters
        real :: mu_0_3, mu_x_3, mu_y_3 ! Influence of vertex 3 on the doublet integral parameters
        logical :: on_kutta_edge ! Whether this panel belongs to a Kutta edge
        logical :: wake_panel ! Whether this panel belongs to a wake
        logical :: shock_panel ! Whether this panel belongs to a shock
        integer,dimension(:),allocatable :: vertex_indices ! Indices of this panel's vertices in the mesh vertex array
        type(list) :: opposing_kutta_panels ! Indices of panels opposite this one on the Kutta edge(s)
        type(list) :: abutting_panels ! Indices of panels abutting this one (not across Kutta edge)

        contains

            procedure :: panel_init_3
            procedure :: panel_init_4
            generic :: init => panel_init_3 !, panel_init_4 ! I only want to deal with 3-sided panels for now
            procedure :: calc_derived_properties =>panel_calc_derived_properties
            procedure :: calc_area => panel_calc_area
            procedure :: calc_normal => panel_calc_normal
            procedure :: calc_centroid => panel_calc_centroid
            procedure :: calc_coord_transform => panel_calc_coord_transform
            procedure :: calc_edge_tangents => panel_calc_edge_tangents
            procedure :: get_vertex_loc => panel_get_vertex_loc
            procedure :: get_vertex_index => panel_get_vertex_index
            procedure :: touches_vertex => panel_touches_vertex
            procedure :: point_to_vertex_clone => panel_point_to_vertex_clone
            procedure :: get_source_potential => panel_get_source_potential
            procedure :: get_source_velocity => panel_get_source_velocity
            procedure :: get_doublet_potential_influence => panel_get_doublet_potential_influence
            procedure :: calc_hH_1_1_3 => panel_calc_hH_1_1_3

    end type panel

    
contains


    subroutine panel_init_3(this, v1, v2, v3, i1, i2, i3, index)
        ! Initializes a 3-panel

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3
        integer,intent(in) :: i1, i2, i3, index

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
        this%index = index

        call this%calc_derived_properties()

    end subroutine panel_init_3


    subroutine panel_init_4(this, v1, v2, v3, v4, i1, i2, i3, i4, index)
        ! Initializes a panel with 4 sides

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3, v4
        integer,intent(in) :: i1, i2, i3, i4, index
        
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
        this%index = index

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

        ! Calculate edge tangents
        call this%calc_edge_tangents()

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
        real :: dx1, dx2, dy1, dy2, det

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

        ! Determine influence of vertex doublet strengths on integral parameters
        ! Preliminaries
        dx1 = this%vertices_local(1,1)-this%centroid(1)
        dx2 = this%vertices_local(2,1)-this%centroid(1)
        dy1 = this%vertices_local(1,2)-this%centroid(2)
        dy2 = this%vertices_local(2,2)-this%centroid(2)
        det = dx1*dy2-dx2*dy1

        ! Influence on mu_0
        this%mu_0_1 = 1/3
        this%mu_0_2 = 1/3
        this%mu_0_3 = 1/3

        ! Influence on mu_x
        this%mu_x_1 = 1/(3*det)*(2*dy2+dy1)
        this%mu_x_2 = 1/(3*det)*(-dy2-2*dy1)
        this%mu_x_3 = 1/(3*det)*(-dy2+dy1)

        ! Influence on mu_y
        this%mu_y_1 = 1/(3*det)*(-2*dx2-dx1)
        this%mu_y_2 = 1/(3*det)*(dx2+2*dx1)
        this%mu_y_3 = 1/(3*det)*(dx2-dx1)

    end subroutine panel_calc_coord_transform


    subroutine panel_calc_edge_tangents(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: d
        integer :: i

        ! Allocate memory
        allocate(this%l(this%N))
        allocate(this%t_hat(this%N,3))
        allocate(this%t_hat_local(this%N,3))
        allocate(this%n_hat(this%N,3))
        allocate(this%n_hat_local(this%N,3))

        ! Calculate tangents
        do i=1,this%N

            ! Calculate difference based on index
            ! Edge tangent i starts at vertex i
            if (i==this%N) then
                d = this%get_vertex_loc(1)-this%get_vertex_loc(i)
            else
                d = this%get_vertex_loc(i+1)-this%get_vertex_loc(i)
            end if

            ! Calculate edge length
            this%l(i) = norm(d)

            ! Calculate tangent
            this%t_hat(i,:) = d/this%l(i)
            this%t_hat_local(i,:) = matmul(this%A_t, this%t_hat(i,:))

        end do

        ! Calculate outward normals
        do i=1,this%N

            this%n_hat(i,:) = cross(this%t_hat(i,:), this%normal)
            this%n_hat_local(i,:) = matmul(this%A_t, this%n_hat(i,:))

        end do
    
    end subroutine panel_calc_edge_tangents


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
            if (norm(this%get_vertex_loc(i)-clone%loc) < 1e-10) then

                ! Update pointer
                this%vertices(i)%ptr => clone

                ! Update index
                this%vertex_indices(i) = clone%index

                return

            end if

        end do
    
    end subroutine panel_point_to_vertex_clone


    function panel_get_source_potential(this, eval_point) result(phi)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        real :: phi

        ! Get H(1,1,1)
    
    end function panel_get_source_potential


    function panel_get_source_velocity(this, eval_point) result(vel)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        real,dimension(3) :: vel
        real :: hH113

        ! Get H(2,1,3), H(1,2,3), and hH(1,1,3)
        hH113 = this%calc_hH_1_1_3(eval_point)
    
    end function panel_get_source_velocity


    function panel_get_doublet_potential_influence(this, eval_point) result(phi)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        real :: phi, hH113

        ! Get fundamental integrals
        hH113 = this%calc_hH_1_1_3(eval_point)
    
    end function panel_get_doublet_potential_influence


    function panel_calc_hH_1_1_3(this, eval_point) result(val)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        real,dimension(3) :: r, r_in_plane, d, d1, d2
        real :: h, phi, S_beta, C_beta
        real,dimension(3) :: a, g2, l1, l2, l1l2, s1, s2, c1, c2
        real :: val
        integer :: i

        ! Transform to panel coordinates
        r = eval_point-this%centroid
        r = matmul(this%A_t, r)
        r_in_plane = r
        r_in_plane(3) = 0
        h = r(3)

        ! Calculate intermediate quantities
        do i=1,this%N

            ! Perpendicular distance in plane from evaluation point to edge
            d = this%vertices_local(i,:)-r_in_plane
            a(i) = inner(d, this%n_hat_local(i,:))

            ! Integration lengths on edges
            l1(i) = inner(d, this%t_hat_local(i,:))
            l2(i) = l1(i)+this%l(i)

        end do

        ! Calculate perpendicular distance to edges
        g2 = a**2+h**2

        ! Other intermediate quantities
        l1l2 = l1*l2
        s1 = sqrt(l1**2+g2)
        s2 = sqrt(l2**2+g2)
        c1 = g2+abs(h)*s1
        c2 = g2+abs(h)*s2

        ! Check location of point relative to the panel
        val = 0.0
        if (abs(h) < 1e-10) then

            ! Check if point is colinear with edge
            if (any(g2 == 0.0)) then

                ! Check if point is on edge
                if (any(l1l2 < 0.0)) then
                    write(*,*) "Error: Evaluation point", eval_point, "is on an edge of panel", this%index
                    stop
                
                ! Not on edge
                else

                    ! Calculate H(1,1,3) using complicated formulation
                    ! Loop through edges
                    do i=1,this%N

                        ! Calculate sine and cosine
                        S_beta = a(i)*( l2(i) - l1(i) + ( abs(h) * ( l2(i)**2 - l1(i)**2 ) / ( l2(i) * s1(i) + l1(i) * s2(i) )))
                        C_beta = g2(i) + abs(h)*( s1(i) + s2(i) ) + l1(i)*l2(i) + &
                                 (h**2 * ( g2(i) + l1(i)**2 + l2(i)**2 )/( s1(i)*s2(i) + l1(i)*l2(i)))

                        ! Add to influence
                        val = val + atan2(S_beta, C_beta)

                    end do

                    val = sign(val, h)

                end if

            ! Point is coplanar but not colinear
            else

                ! Loop through vertices to find total angle swept
                phi = 0
                do i=1,this%N

                    ! Get displacement vectors
                    d1 = this%vertices_local(i,:)-r_in_plane
                    if (i == this%N) then
                        d2 = this%vertices_local(1,:)-r_in_plane
                    else
                        d2 = this%vertices_local(i+1,:)-r_in_plane
                    end if

                    ! Calculate angle swept
                    d = cross(d1, d2)
                    phi = phi + asin(norm(d)/(norm(d1)*norm(d2)))

                end do

                ! Check if point is in panel
                ! Analytically, phi will either be +/-2pi or 0
                if (phi > 3.0 .or. phi < -3.0) then
                    val = sign(2.0*pi, h)
                else
                    val = 0.0
                end if

            end if

        ! Calculate H(1,1,3) at point not coplanar with panel (simple formulation)
        else

            ! Loop through edges
            do i=1,this%N

                ! Calculate sine and cosine
                S_beta = a(i) * ( l2(i)*c1(i) - l1(i)*c2(i) )
                C_beta = c1(i)*c2(i) + a(i)**2 * l1(i)*l2(i)

                ! Add to influence
                val = val + atan2(S_beta, C_beta)
            end do

            val = sign(val, h)

        end if

    end function panel_calc_hH_1_1_3

    
end module panel_mod