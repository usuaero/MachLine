! Panel type
module panel_mod

    use linked_list_mod
    use vertex_mod
    use math_mod

    implicit none

    integer :: doublet_order
    integer :: source_order

    type eval_point_geom
        ! Container type for the geometric parameters necessary for calculating a panel's influence on a given field point

        real,dimension(3) :: r, r_local
        real,dimension(2) :: r_in_plane
        real :: h
        real,dimension(3) :: a, g2, l1, l2, s1, s2, c1, c2

    end type eval_point_geom

    type panel
        ! A panel with an arbitrary number of sides

        integer :: N = 3 ! Number of sides/vertices
        integer :: index ! Index of this panel in the mesh array
        type(vertex_pointer),dimension(:),allocatable :: vertices
        real,dimension(3) :: normal ! Normal vector
        real,dimension(:,:),allocatable :: midpoints
        real,dimension(3) :: centroid
        real,dimension(3,3) :: A_t ! Local coordinate transform matrix
        real,dimension(:,:),allocatable :: vertices_local, midpoints_local ! Location of the vertices and edge midpoints described in local coords
        real,dimension(:,:),allocatable :: t_hat, t_hat_local ! Edge unit tangents
        real,dimension(:,:),allocatable :: n_hat, n_hat_local ! Edge unit outward normals
        real,dimension(:),allocatable :: l ! Edge lengths
        real :: A ! Surface area
        real,dimension(:,:),allocatable :: S_mu, S_mu_inv ! Matrices relating doublet strengths to doublet influence parameters
        real,dimension(:,:),allocatable :: S_sigma, S_sigma_inv ! Matrices relating source strengths to source influence parameters
        logical :: on_wake_edge ! Whether this panel belongs to a wake-shedding edge (on the body)
        integer,dimension(:),allocatable :: vertex_indices ! Indices of this panel's vertices in the mesh vertex array
        type(list) :: opposing_panels ! Indices of panels opposite this one on the wake-shedding edge
        type(list) :: abutting_panels ! Indices of panels abutting this one not across wake-shedding edge
        logical :: xy_sym, xz_sym, yz_sym ! Whether this panel is reflected about any planes

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
            procedure :: R_i => panel_R_i
            procedure :: E_i_M_N_K => panel_E_i_M_N_K
            procedure :: F_i_1_1_1 => panel_F_i_1_1_1
            procedure :: calc_F_integrals => panel_calc_F_integrals
            procedure :: calc_H_integrals => panel_calc_H_integrals
            procedure :: calc_integrals => panel_calc_integrals
            procedure :: get_source_potential => panel_get_source_potential
            procedure :: get_field_point_geometry => panel_get_field_point_geometry

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
        allocate(this%vertices_local(this%N,2))
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
        real,dimension(3) :: sum
        integer :: i

        ! Get average of corner points
        sum = 0.
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

        ! Allocate memory
        allocate(this%vertices_local(this%N,2))
        allocate(this%midpoints_local(this%N,2))

        ! Choose first edge tangent as local xi-axis
        ! (will need to be projected for quadrilateral panel)
        d = this%get_vertex_loc(2)-this%get_vertex_loc(1)
        this%A_t(1,:) = d/norm(d)

        ! Panel normal is the zeta axis
        this%A_t(3,:) = this%normal

        ! Calculate eta axis from the other two
        this%A_t(2,:) = cross(this%normal, this%A_t(1,:))

        ! Transform vertex and midpoint coords
        do i=1,this%N

            ! Vertices
            this%vertices_local(i,:) = matmul(this%A_t(1:2,:), this%get_vertex_loc(i)-this%centroid)

            ! Midpoints
            this%midpoints_local(i,:) = matmul(this%A_t(1:2,:), this%midpoints(i,:)-this%centroid)

        end do

        ! Determine influence of vertex doublet strengths on integral parameters
        if (.not. doublet_order .eq. 0) then

            ! Linear distribution
            if (doublet_order .eq. 1) then

                ! Allocate influence matrices
                allocate(this%S_mu(3,3))
                allocate(this%S_mu_inv(3,3))

                ! Set values
                this%S_mu_inv(:,1) = 1.
                this%S_mu_inv(:,2) = this%vertices_local(:,1)
                this%S_mu_inv(:,3) = this%vertices_local(:,2)

                ! Invert
                call matinv(3, this%S_mu_inv, this%S_mu)

            else if (doublet_order .eq. 2) then

                ! Allocate influence matrix
                allocate(this%S_mu(6,6))
                allocate(this%S_mu_inv(6,6))

                ! Set values
                this%S_mu_inv(:,1) = 1.

                this%S_mu_inv(1:3,2) = this%vertices_local(:,1)
                this%S_mu_inv(1:3,3) = this%vertices_local(:,2)
                this%S_mu_inv(1:3,4) = this%vertices_local(:,1)**2
                this%S_mu_inv(1:3,5) = this%vertices_local(:,1)*this%vertices_local(:,2)
                this%S_mu_inv(1:3,6) = this%vertices_local(:,2)**2
                
                this%S_mu_inv(4:6,2) = this%midpoints_local(:,1)
                this%S_mu_inv(4:6,3) = this%midpoints_local(:,2)
                this%S_mu_inv(4:6,4) = this%midpoints_local(:,1)**2
                this%S_mu_inv(4:6,5) = this%midpoints_local(:,1)*this%midpoints_local(:,2)
                this%S_mu_inv(4:6,6) = this%midpoints_local(:,2)**2

                ! Invert
                call matinv(6, this%S_mu_inv, this%S_mu)

            end if
        end if

        ! Determine influence of vertex source strengths on integral parameters
        if (.not. source_order .eq. 0) then

            ! Linear distribution
            if (source_order .eq. 1) then

                ! Allocate influence matrices
                allocate(this%S_sigma(3,3))
                allocate(this%S_sigma_inv(3,3))

                ! Set values
                this%S_sigma_inv(:,1) = 1.
                this%S_sigma_inv(:,2) = this%vertices_local(:,1)
                this%S_sigma_inv(:,3) = this%vertices_local(:,2)

                ! Invert
                call matinv(3, this%S_sigma_inv, this%S_sigma)

            else if (source_order .eq. 2) then

                ! Allocate influence matrix
                allocate(this%S_sigma(6,6))
                allocate(this%S_sigma_inv(6,6))

                ! Set values
                this%S_sigma_inv(:,1) = 1.

                this%S_sigma_inv(1:3,2) = this%vertices_local(:,1)
                this%S_sigma_inv(1:3,3) = this%vertices_local(:,2)
                this%S_sigma_inv(1:3,4) = this%vertices_local(:,1)**2
                this%S_sigma_inv(1:3,5) = this%vertices_local(:,1)*this%vertices_local(:,2)
                this%S_sigma_inv(1:3,6) = this%vertices_local(:,2)**2
                
                this%S_sigma_inv(4:6,2) = this%midpoints_local(:,1)
                this%S_sigma_inv(4:6,3) = this%midpoints_local(:,2)
                this%S_sigma_inv(4:6,4) = this%midpoints_local(:,1)**2
                this%S_sigma_inv(4:6,5) = this%midpoints_local(:,1)*this%midpoints_local(:,2)
                this%S_sigma_inv(4:6,6) = this%midpoints_local(:,2)**2

                ! Invert
                call matinv(6, this%S_sigma_inv, this%S_sigma)

            end if
        end if

    end subroutine panel_calc_coord_transform


    subroutine panel_calc_edge_tangents(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: d
        integer :: i

        ! Allocate memory
        allocate(this%l(this%N))
        allocate(this%t_hat(this%N,3))
        allocate(this%t_hat_local(this%N,2))
        allocate(this%n_hat(this%N,3))
        allocate(this%n_hat_local(this%N,2))

        ! Loop through edges
        do i=1,this%N

            ! Calculate difference based on index
            ! Edge i starts at vertex i
            if (i==this%N) then
                d = this%get_vertex_loc(1)-this%get_vertex_loc(i)
            else
                d = this%get_vertex_loc(i+1)-this%get_vertex_loc(i)
            end if

            ! Calculate edge length
            this%l(i) = norm(d)

            ! Calculate tangent
            this%t_hat(i,:) = d/this%l(i)
            this%t_hat_local(i,:) = matmul(this%A_t(1:2,:), this%t_hat(i,:))

            ! Calculate outward normal
            this%n_hat(i,:) = cross(this%t_hat(i,:), this%normal)
            this%n_hat_local(i,:) = matmul(this%A_t(1:2,:), this%n_hat(i,:))

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


    function panel_R_i(this, eval_point, i) result(R)
        ! Calculates the distance from the given point (in local coords) to the i-th vertex

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        integer,intent(in) :: i
        real :: R

        R = sqrt((this%vertices_local(i,1)-eval_point(1))**2 + (this%vertices_local(i,2)-eval_point(2))**2+eval_point(3)**2)

    end function panel_R_i


    function panel_get_field_point_geometry(this, eval_point) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(eval_point_geom) :: geom

        real,dimension(2) :: d
        integer :: i

        ! Store point
        geom%r = eval_point

        ! Transform to panel coordinates
        geom%r_local = matmul(this%A_t, geom%r-this%centroid)
        geom%r_in_plane = geom%r_local(1:2)
        geom%h = geom%r_local(3)

        ! Calculate intermediate quantities
        do i=1,this%N

            ! Perpendicular distance in plane from evaluation point to edge
            d = this%vertices_local(i,:)-geom%r_in_plane
            geom%a(i) = inner2(d, this%n_hat_local(i,:))

            ! Integration lengths on edges
            geom%l1(i) = inner2(d, this%t_hat_local(i,:))
            geom%l2(i) = geom%l1(i)+this%l(i)

        end do

        ! Calculate perpendicular distance to edges
        geom%g2 = geom%a**2+geom%h**2

        ! Other intermediate quantities
        geom%s1 = sqrt(geom%l1**2+geom%g2) ! Distance to first vertex (at this point, this method is more efficient than taking the norm)
        geom%s2 = sqrt(geom%l2**2+geom%g2) ! Distance to second vertex
        geom%c1 = geom%g2+abs(geom%h)*geom%s1
        geom%c2 = geom%g2+abs(geom%h)*geom%s2

    end function panel_get_field_point_geometry


    function panel_E_i_M_N_K(this, geom, i, M, N, K) result(E)
        ! Calculates E_i(M,N,K)

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: i, M, N, K

        real :: E, E_1, E_2

        ! Evaluate at start vertex
        E_1 = ((this%vertices_local(i,1)-geom%r_local(1))**(M-1)*(this%vertices_local(i,2)-geom%r_local(2))**(N-1))&
              /this%R_i(geom%r_local, i)**K

        ! Evaluate at end vertex
        if (i .eq. this%N) then
            E_2 = ((this%vertices_local(1,1)-geom%r_local(1))**(M-1)*(this%vertices_local(1,2)-geom%r_local(2))**(N-1))&
                  /this%R_i(geom%r_local, 1)**K
        else
            E_2 = ((this%vertices_local(i,1)-geom%r_local(1))**(M-1)*(this%vertices_local(i,2)-geom%r_local(2))**(N-1))&
                  /this%R_i(geom%r_local, i)**K
        end if

        ! Calculate difference
        E = E_2-E_1

    end function panel_E_i_M_N_K


    function panel_F_i_1_1_1(this, geom, i) result(F)
        ! Calculates F_i(1,1,1)

        implicit none

        class(panel) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: i
        real :: F

        real :: x1, x2

        ! Calculate intermediate quantities
        x1 = sqrt(geom%l1(i)**2+geom%g2(i))
        x2 = sqrt(geom%l2(i)**2+geom%g2(i))

        ! Calculate F(1,1,1)
        ! Below edge
        if (geom%l1(i) >= 0. .and. geom%l2(i) >= 0.) then
            F = log((x2+geom%l2(i))/(x1+geom%l1(i)))
        
        ! Above edge
        else if (geom%l1(i) < 0. .and. geom%l2(i) < 0.) then
            F = log((x1-geom%l2(i))/(x2-geom%l2(i)))

        ! Within edge
        else
            F = log(((x1-geom%l1(i))*(x2+geom%l2(i)))/geom%g2(i))
        end if
        
    end function panel_F_i_1_1_1


    function panel_calc_F_integrals(this, geom, proc_H, MXK, MXQ, NHK, dH) result (F)
        ! Calculates the F integrals necessary

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: proc_H, MXK, MXQ, NHK
        real,intent(in) :: dH
        real,dimension(:,:,:,:),allocatable :: F

        real :: E1, E2, v_xi, v_eta
        real,dimension(3) :: d
        integer :: i, MXFK, NFK, proc_F, k, m, n

        ! Determine which F integrals are needed
        NFK = 16
        if (proc_H .eq. 1) then
            MXFK = MXK - 2
        else
            MXFK = NHK+MXK-2
        end if
        
        ! Allocate integral storage
        allocate(F(this%N,1:MXQ,1:MXQ,1:MXFK+NFK))

        ! Loop through edges
        do i=1,this%N

            ! Store edge derivs
            v_xi = this%n_hat_local(i,1)
            v_eta = this%n_hat_local(i,2)

            ! Determine which procedure is necessary for this edge
            if (sqrt(geom%g2(i)) >= 0.01*dH) then
                proc_F = 4
            else
                proc_F = 5
            end if

            ! Calculate F(1,1,1)
            F(i,1,1,1) = this%F_i_1_1_1(geom, i)

            ! Procedure 4: not close to perimeter
            if (proc_F .eq. 4) then
                
                ! Calculate F(1,1,K) integrals
                do k=3,MXFK,2

                    ! Get necessary E
                    E1 = this%E_i_M_N_K(geom, i, 2, 1, k-2)
                    E2 = this%E_i_M_N_K(geom, i, 1, 2, k-2)

                    ! Calculate F
                    F(i,1,1,k) = 1./(geom%g2(i)*(k-2))*((k-3)*F(i,1,1,k-2)-v_eta*E1+v_xi*E2)
                end do

            ! Procedure 5: close to perimeter
            else

                ! Initialize
                F(i,1,1,MXFK+NFK) = 0.

                ! Calculate other F(1,1,K) integrals
                do k=MXFK+NFK,5,2

                    ! Get necessary E
                    E1 = this%E_i_M_N_K(geom, i, 2, 1, k-2)
                    E2 = this%E_i_M_N_K(geom, i, 1, 2, k-2)

                    ! Calculate F
                    F(i,1,1,k-2) = 1./(k-3)*(geom%g2(i)*(k-2)*F(i,1,1,k)+v_eta*E1-v_xi*E2)

                end do
            end if

            ! Calculate other F integrals (same for both procedures)
            if (v_eta <= v_xi) then ! Case a
                
                ! Calculate F(1,N,1) integrals
                do n=2,MXQ

                    ! Get E
                    E1 = this%E_i_M_N_K(geom, i, 1, n-1, -1)

                    if (.not. n .eq. 2) then
                        F(i,1,N,1) = 1./(n-1)*((2*n-3)*geom%a(i)*v_eta*F(i,1,n-1,1) &
                                     - (n-2)*(geom%a(i)**2+v_xi**2*geom%h**2)*F(i,1,n-2,1) &
                                     + v_xi*E1)
                    else
                        F(i,1,N,1) = 1./(n-1)*((2*n-3)*geom%a(i)*v_eta*F(i,1,n-1,1) + v_xi*E1)
                    end if
                end do

                ! Calculate F(M,N,1) inegrals
                do m=2,MXQ
                    do n=1,MXQ-m+1
                        F(i,m,n,1) = -v_eta/v_xi*F(i,m-1,n+1,1) + geom%a(i)/v_xi*F(i,m-1,n,1)
                    end do
                end do

            else ! Case b

                ! Calculate F(M,1,1) integrals
                do m=2,MXQ

                    ! Get E
                    E1 = this%E_i_M_N_K(geom, i, m-1, 1, -1)

                    if (.not. m .eq. 2) then
                        F(i,m,1,1) = 1./(m-1)*((2*m-3)*geom%a(i)*v_xi*F(i,m-1,1,1) &
                                     - (m-2)*(geom%a(i)**2+v_eta**2*geom%h**2)*F(i,m-2,1,1) &
                                     - v_eta*E1)
                    else
                        F(i,m,1,1) = 1./(m-1)*((2*m-3)*geom%a(i)*v_xi*F(i,m-1,1,1) - v_eta*E1)
                    end if
                end do

                ! Calculate F(M,N,1) integrals
                do n=2,MXQ
                    do m=1,MXQ-n+1
                        F(i,m,n,1) = -v_xi/v_eta*F(i,m+1,n-1,1)+geom%a(i)/v_eta*F(i,m,n-1,1)
                    end do
                end do

            end if
        end do

        ! Calculate F(1,2,K) integrals
        do k=3,MXK-2,2
            F(i,1,2,k) = v_eta*geom%a(i)*F(i,1,1,k)-v_xi/(k-2)*this%E_i_M_N_K(geom, i, 1, 1, k-2)
        end do

        ! Calculate F(1,N,K) integrals
        do n=3,MXQ
            do k=3,MXK-2,2
                F(i,1,n,k) = 2.*geom%a(i)*v_eta*F(i,1,n-1,k) &
                             -(geom%a(i)**2+v_xi**2*geom%h**2)*F(i,1,n-2,k) &
                             +v_xi**2*F(i,1,n-2,k-2)
            end do
        end do


    end function panel_calc_F_integrals


    function panel_calc_H_integrals(this, geom, proc_H, MXK, MXQ, NHK, F) result(H)
        ! Calculates the necessary H integrals

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: proc_H, MXK, MXQ, NHK
        real,dimension(:,:,:,:),allocatable,intent(in) :: F
        real,dimension(:,:,:),allocatable :: H

        real :: S, C
        integer :: i, m, n, k

        ! Allocate integral storage
        allocate(H(1:MXQ,1:MXQ,1:MXK+NHK))

        ! Procedure 1: not close to panel plane
        if (proc_H .eq. 1) then

            ! Calculate H(1,1,1)
            H(1,1,1) = 0.
            do i=1,this%N
        
                ! Add surface integral
                S = geom%a(i)*(geom%l2(i)*geom%c1(i) - geom%l1(i)*geom%c2(i))
                C = geom%c1(i)*geom%c2(i) + geom%a(i)**2*geom%l1(i)*geom%l2(i)
                H(1,1,1) = H(1,1,1) - abs(geom%h)*atan2(S, C)
        
                ! Add line integral
                H(1,1,1) = H(1,1,1) + geom%a(i)*F(i,1,1,1)
        
            end do

            ! Calculate H(1,1,K) integrals
            do k=3,MXK,2
                H(1,1,k) = 1./((k-2)*geom%h**2)*((k-4)*H(1,1,k-2)+sum(geom%a*F(:,1,1,k-2)))
            end do

        ! Procedures 2 and 3: close to panel plane
        else if (proc_H .eq. 2 .or. proc_H .eq. 3) then

            ! Initialize
            H(1,1,NHK+MXK) = 0.

            ! Calculate H(1,1,K) integrals
            do k=NHK+MXK,3,-2
                H(1,1,k-2) = 1./(k-4)*(geom%h**2*(k-2)*h(1,1,k)-sum(geom%a*F(:,1,1,k-2)))
            end do

        end if

        ! Calculate H(2,N,1) integrals
        do n=1,MXQ-1
            H(2,n,1) = 1./(n+1)*(geom%h**2*sum(this%n_hat_local(:,1)*F(:,1,n,1)) &
                       + sum(geom%a*F(:,2,n,1)))
        end do

        ! Calculate H(1,N,1) integrals
        do n=2,MXQ
            if (.not. n .eq. 2) then
                H(2,n,1) = 1./n*(-geom%h**2*(n-2)*H(1,n-2,1) & 
                           + geom%h**2*sum(this%n_hat_local(:,2)*F(i,1,n-1,1)) &
                           + sum(geom%a*F(:,1,n,1)))
            else
                H(2,n,1) = 1./n*(geom%h**2*sum(this%n_hat_local(:,2)*F(i,1,n-1,1)) &
                           + sum(geom%a*F(:,1,n,1)))
            end if
        end do

        ! Calculate H(M,N,1) integrals
        do m=3,MXQ
            do n=1,MXQ-m+1
                if (.not. m .eq. 2) then
                    H(m,n,1) = 1./(m+n-1)*(-geom%h**2*(m-2)*H(m-2,n,1) &
                               + geom%h**2*sum(this%n_hat_local(:,1)*F(:,m-1,n,1)) &
                               + sum(geom%a*F(:,m,n,1)))
                else
                    H(m,n,1) = 1./(m+n-1)*(geom%h**2*sum(this%n_hat_local(:,1)*F(:,m-1,n,1)) &
                               + sum(geom%a*F(:,m,n,1)))
                end if
            end do
        end do

        ! Calculate H(1,N,K) integrals
        do n=2,MXQ
            do k=3,MXK,2
                if (.not. n .eq. 2) then
                    H(1,n,k) = 1./(k-2)*((n-2)*H(1,n-2,k-2) &
                               -sum(this%n_hat_local(:,2)*F(:,1,n-1,k-2)))
                else
                    H(1,n,k) = -1./(k-2)*sum(this%n_hat_local(:,2)*F(:,1,n-1,k-2))
                end if
            end do
        end do

        ! Calculate H(2,N,K) integrals
        do n=1,MXQ-1
            do k=3,MXK,2
                H(2,n,k) = -1./(k-2)*sum(this%n_hat_local(:,1)*F(:,1,n,k-2))
            end do
        end do

        ! Calculate remaining H(M,N,K) integrals
        do m=3,MXQ
            do n=1,MXQ-m+1
                do k=3,MXK,2
                    H(m,n,k) = -h(M-2,N+2,k)-geom%h**2*H(m-2,n,k)+H(m-2,n,k-2)
                end do
            end do
        end do

        ! Convert H* to H in case of Procedure 3
        if (proc_H .eq. 3) then
            do m=1,MXQ
                do n=1,MXQ
                    do k=1,MXK+NHK,2
                        H(m,n,k) = H(m,n,k)+2.*pi*nu_M_N_K(m, n, k)*geom%h**(m+n-k)
                    end do
                end do
            end do
        end if

    end function panel_calc_H_integrals


    function nu_M_N_K(M, N, K) result(nu)
        ! Calculates nu(M,N,K) based on Eq. (D.51) in Johnson 1980

        implicit none

        integer,intent(in) :: M, N, K
        real :: nu

        integer :: mm, nn, kk, i

        ! Check for even M or N
        if (mod(M, 2) .eq. 0 .or. mod(N, 2) .eq. 0) then
            nu = 0
        else
             
            ! Initialize
            mm = 1
            nn = 1
            kk = 1

            ! Run factorials (ish. not sure what you'd call these...)
            if (.not. M .eq. 1) then
                do i=1,M-2,2
                    mm = mm*i
                end do
            end if

            if (.not. N .eq. 1) then
                do i=1,N-2,2
                    nn = nn*i
                end do
            end if

            do i=K-2,K-M-N,2
                kk = kk*i
            end do

            ! Calculate nu
            nu = mm*nn/kk

        end if

    end function nu_M_N_K


    subroutine panel_calc_integrals(this, geom, influence_type, singularity_type, H, F)
        ! Calculates the H and F integrals necessary for the given influence

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        character(len=*),intent(in) :: influence_type, singularity_type
        real,dimension(:,:,:,:),allocatable,intent(out) :: F
        real,dimension(:,:,:),allocatable,intent(out) :: H

        real :: dH
        real,dimension(3) :: d
        integer :: MXQ, MXK, NHK, proc_H

        ! Determine which H integrals are needed based on distribution and type of influence
        if (singularity_type .eq. "source") then
            if (influence_type .eq. "potential") then
                if (source_order .eq. 0) then
                    MXQ = 1
                    MXK = 1
                else if (source_order .eq. 1) then
                    MXQ = 2
                    MXK = 1
                end if
            else if (influence_type .eq. "velocity") then
                if (source_order .eq. 0) then
                    MXQ = 1
                    MXK = 3
                else if (source_order .eq. 1) then
                    MXQ = 3
                    MXK = 3
                end if
            end if
        else if (singularity_type .eq. "doublet") then
            if (influence_type .eq. "potential") then
                if (doublet_order .eq. 0) then
                    MXQ = 1
                    MXK = 3
                else if (doublet_order .eq. 1) then
                    MXQ = 2
                    MXK = 3
                else if (doublet_order .eq. 2) then
                    MXQ = 3
                    MXK = 3
                end if
            else if (influence_type .eq. "velocity") then
                if (doublet_order .eq. 0) then
                    MXQ = 1
                    MXK = 3
                else if (doublet_order .eq. 1) then
                    MXQ = 1
                    MXK = 3
                else if (doublet_order .eq. 2) then
                    MXQ = 4
                    MXK = 5
                end if
            end if
        end if

        ! Check distance to panel perimeter (minimum perpendicular distance to edge)
        dH = minval(sqrt(geom%g2))

        ! Determine which procedure needs to be used
        if (abs(geom%h) >= 0.01*dH) then
            proc_H = 1 ! Not near plane of panel
            NHK = 0

        else

            ! Check if the projected point falls inside the panel
            d(1) = inner2(geom%r_in_plane, this%n_hat_local(1,:))
            d(2) = inner2(geom%r_in_plane, this%n_hat_local(2,:))
            d(3) = inner2(geom%r_in_plane, this%n_hat_local(3,:))

            ! Outside panel (Procedure 2)
            if (all(d > 0)) then
                proc_H = 2

            ! Inside panel (Procedure 3)
            else
                proc_H = 3

            end if
            NHK = 16
        end if

        ! Calculate F integrals
        F = this%calc_F_integrals(geom, proc_H, MXK, MXQ, NHK, dH)

        ! Calculate H integrals
        H = this%calc_H_integrals(geom, proc_H, MXK, MXQ, NHK, F)

    end subroutine panel_calc_integrals


    function panel_get_source_potential(this, eval_point) result(phi)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        real :: phi

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        ! Calculate geometric parameters
        geom = this%get_field_point_geometry(eval_point)

        ! Get H integrals
        call this%calc_integrals(geom, "potential", "source", H, F)

        ! Compute induced potential
        if (source_order .eq. 0) then
            phi = -1./(4.*pi)*H(1,1,1)
        end if
    
    end function panel_get_source_potential

    
end module panel_mod
