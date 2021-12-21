! Panel type
module panel_mod

    use helpers_mod
    use linked_list_mod
    use vertex_mod
    use math_mod
    use flow_mod

    implicit none

    integer :: doublet_order
    integer :: source_order
    integer :: eval_count ! Developer counter for optimization purposes
    logical :: debug = .false. ! Developer toggle
    character(len=:),allocatable :: influence_calc_type ! Either 'johnson-ehlers' or 'gauss quad'

    type eval_point_geom
        ! Container type for the geometric parameters necessary for calculating a panel's influence on a given field point

        real,dimension(3) :: r
        real,dimension(2) :: r_in_plane
        real :: h
        real,dimension(3) :: a, g2, l1, l2, s1, s2, c1, c2, g

    end type eval_point_geom


    type dod
        ! Container type for parameters of whether a panel lies in a point's domain of dependence

        logical :: in_dod
        logical,dimension(3) :: verts_in_dod, edges_in_dod

    end type dod


    type panel
        ! A three-sided panel

        integer :: N = 3 ! Number of sides/vertices
        integer :: index ! Index of this panel in the mesh array
        type(vertex_pointer),dimension(:),allocatable :: vertices
        real :: radius
        real,dimension(3) :: normal, conormal ! Normal and conormal vectors
        real,dimension(:,:),allocatable :: midpoints
        real,dimension(3) :: centroid
        real,dimension(3,3) :: A_g_to_l, A_s_to_ls, A_g_to_ls, A_ls_to_g ! Coordinate transformation matrices
        real,dimension(:,:),allocatable :: vertices_l, midpoints_l ! Location of the vertices and edge midpoints described in local coords
        real,dimension(:,:),allocatable :: t_hat_g, t_hat_l ! Edge unit tangents
        real,dimension(:,:),allocatable :: n_hat_g, n_hat_l ! Edge unit outward normals
        real,dimension(:),allocatable :: l ! Edge lengths
        real :: A ! Surface area
        real,dimension(:,:),allocatable :: S_mu_inv, S_sigma_inv ! Matrix relating doublet/source strengths to doublet/source influence parameters
        logical :: on_wake_edge ! Whether this panel belongs to a wake-shedding edge (on the body)
        integer,dimension(:),allocatable :: vertex_indices ! Indices of this panel's vertices in the mesh vertex array
        logical :: in_wake ! Whether this panel belongs to a wake mesh
        integer,dimension(3) :: abutting_panels ! Indices of panels abutting this one
        real :: r ! Panel inclination indicator; r=-1 -> superinclined, r=1 -> subinclined
        logical,dimension(3) :: edge_subsonic ! Whether each edge is subsonic

        contains

            procedure :: init => panel_init_3
            procedure :: calc_derived_properties =>panel_calc_derived_properties
            procedure :: calc_area => panel_calc_area
            procedure :: calc_normal => panel_calc_normal
            procedure :: calc_centroid => panel_calc_centroid
            procedure :: calc_transforms => panel_calc_transforms
            procedure :: calc_g_to_l_transform => panel_calc_g_to_l_transform
            procedure :: calc_g_to_ls_transform => panel_calc_g_to_ls_transform
            procedure :: calc_edge_tangents => panel_calc_edge_tangents
            procedure :: calc_singularity_matrices => panel_calc_singularity_matrices
            procedure :: add_abutting_panel => panel_add_abutting_panel
            procedure :: get_vertex_loc => panel_get_vertex_loc
            procedure :: get_vertex_index => panel_get_vertex_index
            procedure :: touches_vertex => panel_touches_vertex
            procedure :: point_to_vertex_clone => panel_point_to_vertex_clone
            procedure :: check_dod => panel_check_dod
            procedure :: get_field_point_geometry => panel_get_field_point_geometry
            procedure :: E_i_M_N_K => panel_E_i_M_N_K
            procedure :: F_i_1_1_1 => panel_F_i_1_1_1
            procedure :: calc_F_integrals => panel_calc_F_integrals
            procedure :: calc_H_integrals => panel_calc_H_integrals
            procedure :: calc_integrals => panel_calc_integrals
            procedure :: get_source_potential => panel_get_source_potential
            procedure :: get_source_velocity => panel_get_source_velocity
            procedure :: get_doublet_potential => panel_get_doublet_potential
            procedure :: get_doublet_velocity => panel_get_doublet_velocity
            procedure :: get_velocity_jump => panel_get_velocity_jump

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

        ! Initialize a few things
        this%abutting_panels = 0

        call this%calc_derived_properties()

    end subroutine panel_init_3


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
        real :: x

        ! Get average of corner points
        sum = 0.
        do i=1,this%N
            sum = sum + this%get_vertex_loc(i)
        end do

        ! Set centroid
        this%centroid = sum/this%N

        ! Calculate radius
        this%radius = 0.
        do i=1,this%N

            ! Distance to i-th vertex
            x = dist(this%get_vertex_loc(i), this%centroid)

            ! Check max
            if (x > this%radius) then
                this%radius = x
            end if

        end do

    end subroutine panel_calc_centroid


    subroutine panel_calc_transforms(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        integer :: i

        ! Calculate transforms
        call this%calc_g_to_l_transform(freestream)
        call this%calc_g_to_ls_transform(freestream)

        ! Transform vertex and midpoint coords to l
        allocate(this%vertices_l(this%N,2))
        allocate(this%midpoints_l(this%N,2))
        do i=1,this%N

            ! Vertices
            this%vertices_l(i,:) = matmul(this%A_g_to_l(1:2,:), this%get_vertex_loc(i)-this%centroid)

            ! Midpoints
            this%midpoints_l(i,:) = matmul(this%A_g_to_l(1:2,:), this%midpoints(i,:)-this%centroid)

        end do

        ! Calculate properties dependent on the transforms
        call this%calc_edge_tangents()
        call this%calc_singularity_matrices()

        ! Check character of edges
        do i=1,this%N
            this%edge_subsonic(i) = freestream%C0_inner(this%t_hat_g(i,:), this%t_hat_g(i,:)) > 0.
        end do

    end subroutine panel_calc_transforms


    subroutine panel_calc_g_to_l_transform(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        ! Calculate local eta axis
        this%A_g_to_l(2,:) = cross(this%normal, freestream%c0)
        this%A_g_to_l(2,:) = this%A_g_to_l(2,:)/norm(this%A_g_to_l(2,:))

        ! Calculate local xi axis
        this%A_g_to_l(1,:) = cross(this%A_g_to_l(2,:), this%normal)
        this%A_g_to_l(1,:) = this%A_g_to_l(1,:)/norm(this%A_g_to_l(1,:))

        ! Store local zeta axis
        this%A_g_to_l(3,:) = this%normal

    end subroutine panel_calc_g_to_l_transform


    subroutine panel_calc_g_to_ls_transform(this, freestream)
        ! Calculates the necessary transformations to move from global to local, scaled coordinates

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(3) :: u0, v0
        real,dimension(3,3) :: C0, B0, I
        real :: x
        integer :: j

        ! Get basis vectors
        u0 = this%A_g_to_l(1,:)
        v0 = this%A_g_to_l(2,:)

        ! Calculate compressible parameters
        this%conormal = matmul(freestream%psi, this%normal)
        x = inner(this%normal, this%conormal)
        this%r = sign(1., x) ! r=-1 -> superinclined, r=1 -> subinclined

        ! Check for Mach-inclined panel
        if (inner(this%normal, this%conormal) == 0.) then
            write(*,*) "    !!! Mach-inclined panels are not allowed. Panel", this%index, "is Mach-inclined. Quitting..."
            stop
        end if

        ! Calculate intermediate matrices
        C0 = 0.
        B0 = 0.
        I = 0.

        ! Construct identity matrix
        do j=1,3
            I(j,j) = 1.
        end do

        ! Construct other matrices
        B0 = outer(freestream%c0, freestream%c0)
        C0 = freestream%s*freestream%B**2*I + freestream%M_inf**2*B0
        B0 = I - freestream%M_inf**2*B0

        ! Calculate transformation
        this%A_g_to_ls(1,:) = 1./sqrt(abs(x))*matmul(C0, u0)
        this%A_g_to_ls(2,:) = this%r*freestream%s/freestream%B*matmul(C0, v0)
        this%A_g_to_ls(3,:) = freestream%B/sqrt(abs(x))*this%normal

        ! Calculate inverse
        call matinv(3, this%A_g_to_ls, this%A_ls_to_g)
    
    end subroutine panel_calc_g_to_ls_transform


    subroutine panel_calc_edge_tangents(this)

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3) :: d
        integer :: i

        ! Allocate memory
        allocate(this%l(this%N))
        allocate(this%t_hat_g(this%N,3))
        allocate(this%t_hat_l(this%N,2))
        allocate(this%n_hat_g(this%N,3))
        allocate(this%n_hat_l(this%N,2))

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
            this%t_hat_g(i,:) = d/this%l(i)
            this%t_hat_l(i,:) = matmul(this%A_g_to_l(1:2,:), this%t_hat_g(i,:))

            ! Calculate outward normal
            this%n_hat_g(i,:) = cross(this%t_hat_g(i,:), this%normal)
            this%n_hat_l(i,:) = matmul(this%A_g_to_l(1:2,:), this%n_hat_g(i,:))

        end do
    
    end subroutine panel_calc_edge_tangents


    subroutine panel_calc_singularity_matrices(this)
        ! Calculates the matrices which relate the singularity strengths to the singularity parameters

        implicit none

        class(panel),intent(inout) :: this

        real,dimension(:,:),allocatable :: S_mu, S_sigma

        ! Determine influence of vertex doublet strengths on integral parameters
        if (.not. doublet_order .eq. 0) then

            ! Linear distribution
            if (doublet_order .eq. 1) then

                ! Allocate influence matrices
                allocate(S_mu(3,3))
                allocate(this%S_mu_inv(3,3))

                ! Set values
                S_mu(:,1) = 1.
                S_mu(:,2) = this%vertices_l(:,1)
                S_mu(:,3) = this%vertices_l(:,2)

                ! Invert
                call matinv(3, S_mu, this%S_mu_inv)

            else if (doublet_order .eq. 2) then

                ! Allocate influence matrix
                allocate(S_mu(6,6))
                allocate(this%S_mu_inv(6,6))

                ! Set values
                S_mu(:,1) = 1.

                S_mu(1:3,2) = this%vertices_l(:,1)
                S_mu(1:3,3) = this%vertices_l(:,2)
                S_mu(1:3,4) = this%vertices_l(:,1)**2
                S_mu(1:3,5) = this%vertices_l(:,1)*this%vertices_l(:,2)
                S_mu(1:3,6) = this%vertices_l(:,2)**2
                
                S_mu(4:6,2) = this%midpoints_l(:,1)
                S_mu(4:6,3) = this%midpoints_l(:,2)
                S_mu(4:6,4) = this%midpoints_l(:,1)**2
                S_mu(4:6,5) = this%midpoints_l(:,1)*this%midpoints_l(:,2)
                S_mu(4:6,6) = this%midpoints_l(:,2)**2

                ! Invert
                call matinv(6, S_mu, this%S_mu_inv)

            end if
            
            deallocate(S_mu)

        end if

        ! Determine influence of vertex source strengths on integral parameters
        if (.not. source_order .eq. 0) then

            ! Linear distribution
            if (source_order .eq. 1) then

                ! Allocate influence matrices
                allocate(S_sigma(3,3))
                allocate(this%S_sigma_inv(3,3))

                ! Set values
                S_sigma(:,1) = 1.
                S_sigma(:,2) = this%vertices_l(:,1)
                S_sigma(:,3) = this%vertices_l(:,2)

                ! Invert
                call matinv(3, S_sigma, this%S_sigma_inv)

            else if (source_order .eq. 2) then

                ! Allocate influence matrix
                allocate(S_sigma(6,6))
                allocate(this%S_sigma_inv(6,6))

                ! Set values
                S_sigma(:,1) = 1.

                S_sigma(1:3,2) = this%vertices_l(:,1)
                S_sigma(1:3,3) = this%vertices_l(:,2)
                S_sigma(1:3,4) = this%vertices_l(:,1)**2
                S_sigma(1:3,5) = this%vertices_l(:,1)*this%vertices_l(:,2)
                S_sigma(1:3,6) = this%vertices_l(:,2)**2
                
                S_sigma(4:6,2) = this%midpoints_l(:,1)
                S_sigma(4:6,3) = this%midpoints_l(:,2)
                S_sigma(4:6,4) = this%midpoints_l(:,1)**2
                S_sigma(4:6,5) = this%midpoints_l(:,1)*this%midpoints_l(:,2)
                S_sigma(4:6,6) = this%midpoints_l(:,2)**2

                ! Invert
                call matinv(6, S_sigma, this%S_sigma_inv)

            end if

            deallocate(S_sigma)

        end if
    
    end subroutine panel_calc_singularity_matrices


    subroutine panel_add_abutting_panel(this, panel_ind)
        ! Adds a panel index to this panel's list of abutting panels

        implicit none

        class(panel),intent(inout) :: this
        integer,intent(in) :: panel_ind

        integer :: i

        do i=1,3
            if (this%abutting_panels(i) == 0) then
                this%abutting_panels(i) = panel_ind
                return
            end if
        end do
    
    end subroutine panel_add_abutting_panel


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
            if (norm(this%get_vertex_loc(i)-clone%loc) < 1e-12) then

                ! Update pointer
                this%vertices(i)%ptr => clone

                ! Update index
                this%vertex_indices(i) = clone%index

                return

            end if

        end do
    
    end subroutine panel_point_to_vertex_clone


    function panel_check_dod(this, eval_point, freestream) result(dod_info)
        ! Determines how (if) this panel lies within the domain of dependence of the evaluation point

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod) :: dod_info

        real,dimension(3) :: d
        integer :: i, N_verts_in_dod, N_edges_in_dod
        real :: C_theta, x, y
        logical :: totally_out, totally_in, centroid_in, radius_smaller

        ! Fast check based on centroid and radius (Epton & Magnus p. J.3-1)

        ! First, check if the centroid is in the dod
        centroid_in = freestream%point_in_dod(this%centroid, eval_point)

        ! Check if the centroid is further from the Mach cone than the panel radius
        d = this%centroid-eval_point
        x = inner(d, freestream%c0)
        y = sqrt(norm(d)**2-x**2)
        if (y >= freestream%B*x .and. sqrt((x+freestream%B*y)**2/(1.+freestream%B**2)) >= this%radius) then
            radius_smaller = .true.
        else if (y < freestream%B*x .and. norm(d) >= this%radius) then
            radius_smaller = .true.
        else
            radius_smaller = .false.
        end if

        ! Determine condition
        totally_out = .not. centroid_in .and. radius_smaller
        totally_in = centroid_in .and. radius_smaller

        ! Set parameters
        if (totally_out) then
            
            dod_info%in_dod = .false.
            dod_info%verts_in_dod = .false.
            dod_info%edges_in_dod = .false.

        else if (totally_in) then
            
            dod_info%in_dod = .true.
            dod_info%verts_in_dod = .true.
            dod_info%edges_in_dod = .true.

        ! If it is not guaranteed to be totally out or in, then check all the vertices and edges
        else

            ! Check each of the vertices
            N_verts_in_dod = 0
            do i=1,this%N

                ! Check
                dod_info%verts_in_dod(i) = freestream%point_in_dod(this%get_vertex_loc(i), eval_point)
                if (dod_info%verts_in_dod(i)) then
                    N_verts_in_dod = N_verts_in_dod + 1
                end if

            end do

            ! Check edges (if 2 or more vertices are in the dod, then how the edges fall is known)
            if (N_verts_in_dod < 2) then

                ! Initialize count
                N_edges_in_dod = 0

                ! Loop through edges
                do i=1,this%N

                    ! For subsonic edges, this is entirely dependent on the most downstream endpoint

                end do

                ! Check if the dod is encompassed by the panel (for superinclined panels)
                if (this%r == -1. .and. N_edges_in_dod == 0) then

                end if

            else

                ! Store which edges are in based on which vertices are in
                do i=1,this%N-1
                    dod_info%edges_in_dod(i) = dod_info%verts_in_dod(i) .or. dod_info%verts_in_dod(i+1)
                end do
                dod_info%edges_in_dod(this%N) = dod_info%verts_in_dod(this%N) .or. dod_info%verts_in_dod(1)

            end if

        end if

    
    end function panel_check_dod


    function panel_get_field_point_geometry(this, eval_point) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(eval_point_geom) :: geom

        real,dimension(2) :: d
        real,dimension(3) :: r_l
        integer :: i

        ! Store point
        geom%r = eval_point

        ! Transform to panel coordinates
        r_l = matmul(this%A_g_to_l, geom%r-this%centroid)
        geom%r_in_plane = r_l(1:2)
        geom%h = r_l(3)

        ! Calculate intermediate quantities
        do i=1,this%N

            ! Perpendicular distance in plane from evaluation point to edge
            d = this%vertices_l(i,:)-geom%r_in_plane
            geom%a(i) = inner2(d, this%n_hat_l(i,:)) ! a is positive when the evaluation point is towards the interior of the panel

            ! Integration lengths on edges
            geom%l1(i) = inner2(d, this%t_hat_l(i,:))
            geom%l2(i) = geom%l1(i)+this%l(i)

            ! Distance from evaluation point to start vertex
            ! This is not the definition given by Johnson, but this is equivalent.
            ! I believe this method suffers from less numerical error.
            geom%s1(i) = norm(this%get_vertex_loc(i)-geom%r)

        end do

        ! Distance from evaluation point to end vertices
        geom%s2 = cshift(geom%s1, 1)

        ! Calculate perpendicular distance to edges
        geom%g2 = geom%a**2+geom%h**2
        geom%g = sqrt(geom%g2)

        ! Other intermediate quantities
        !geom%s1 = sqrt(geom%l1**2+geom%g2) ! Distance to first vertex (I'm leaving these methods here, as they are Johnson's, but they seem to suffer from buildup of numerical error)
        !geom%s2 = sqrt(geom%l2**2+geom%g2) ! Distance to second vertex
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
        E_1 = ((this%vertices_l(i,1)-geom%r_in_plane(1))**(M-1) &
              *(this%vertices_l(i,2)-geom%r_in_plane(2))**(N-1)) &
              /geom%s1(i)**K

        ! Evaluate at end vertex
        if (i .eq. this%N) then
            E_2 = ((this%vertices_l(1,1)-geom%r_in_plane(1))**(M-1) &
                   *(this%vertices_l(1,2)-geom%r_in_plane(2))**(N-1)) &
                   /geom%s2(i)**K
        else
            E_2 = ((this%vertices_l(i,1)-geom%r_in_plane(1))**(M-1) &
                   *(this%vertices_l(i,2)-geom%r_in_plane(2))**(N-1)) &
                   /geom%s2(i)**K
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

        ! Calculate F(1,1,1)
        ! Below edge
        if (geom%l1(i) >= 0. .and. geom%l2(i) >= 0.) then
            F = log((geom%s2(i)+geom%l2(i))/(geom%s1(i)+geom%l1(i)))
        
        ! Above edge
        else if (geom%l1(i) < 0. .and. geom%l2(i) < 0.) then
            F = log((geom%s1(i)-geom%l1(i))/(geom%s2(i)-geom%l2(i)))

        ! Within edge
        else
            F = log(((geom%s1(i)-geom%l1(i))*(geom%s2(i)+geom%l2(i)))/geom%g2(i))
        end if
        
    end function panel_F_i_1_1_1


    function panel_calc_F_integrals(this, geom, proc_H, MXK, MXQ, NHK) result (F)
        ! Calculates the F integrals necessary

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: proc_H, MXK, MXQ, NHK
        real,dimension(:,:,:,:),allocatable :: F

        real :: E1, E2, v_xi, v_eta, dF
        real,dimension(3) :: d
        real,dimension(this%N) :: min_dist_to_edge
        integer :: i, MXFK, NFK, k, m, n

        ! Determine which F integrals are needed
        NFK = 16
        if (proc_H .eq. 1) then
            MXFK = MXK - 2
        else
            MXFK = NHK+MXK-2
        end if
        
        ! Allocate integral storage
        allocate(F(this%N,1:MXQ,1:MXQ,1:MXFK+NFK), source=0.)

        ! Calculate minimum distance to perimeter of S
        do i=1,this%N

            ! Within edge, the minimum distance is the perpendicular distance
            if (geom%l1(i) < 0. .and. geom%l2(i) >= 0.) then
                min_dist_to_edge(i) = geom%g(i)
        
            ! Otherwise, it is the minimum of the distances to the corners
            else
                min_dist_to_edge(i) = min(geom%s1(i), geom%s2(i))
            end if
        end do
        dF = minval(min_dist_to_edge)

        ! Check for point on perimeter
        if (abs(dF) < 1e-12) then
            write(*,*) "Detected point on perimeter of panel. Quitting..."
            stop
        end if

        ! Loop through edges
        do i=1,this%N

            ! Store edge derivs
            v_xi = this%n_hat_l(i,1)
            v_eta = this%n_hat_l(i,2)

            ! Calculate F(1,1,1)
            F(i,1,1,1) = this%F_i_1_1_1(geom, i)

            ! Procedure 4: not close to perimeter
            if (geom%g(i) >= 0.01*dF) then
                
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
            if (abs(v_eta) <= abs(v_xi)) then ! Case a
                
                ! Calculate F(1,N,1) integrals
                do n=2,MXQ

                    ! Get E
                    E1 = this%E_i_M_N_K(geom, i, 1, n-1, -1)

                    if (n .eq. 2) then
                        F(i,1,N,1) = 1./(n-1)*((2*n-3)*geom%a(i)*v_eta*F(i,1,n-1,1) + v_xi*E1)
                    else
                        F(i,1,N,1) = 1./(n-1)*((2*n-3)*geom%a(i)*v_eta*F(i,1,n-1,1) &
                                     - (n-2)*(geom%a(i)**2+v_xi**2*geom%h**2)*F(i,1,n-2,1) &
                                     + v_xi*E1)
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

                    if (m .eq. 2) then
                        F(i,m,1,1) = 1./(m-1)*((2*m-3)*geom%a(i)*v_xi*F(i,m-1,1,1) - v_eta*E1)
                    else
                        F(i,m,1,1) = 1./(m-1)*((2*m-3)*geom%a(i)*v_xi*F(i,m-1,1,1) &
                                     - (m-2)*(geom%a(i)**2+v_eta**2*geom%h**2)*F(i,m-2,1,1) &
                                     - v_eta*E1)
                    end if
                end do

                ! Calculate F(M,N,1) integrals
                do n=2,MXQ
                    do m=1,MXQ-n+1
                        F(i,m,n,1) = -v_xi/v_eta*F(i,m+1,n-1,1)+geom%a(i)/v_eta*F(i,m,n-1,1)
                    end do
                end do

            end if

            ! Calculate F(1,2,K) integrals
            do k=3,MXK-2,2
                F(i,1,2,k) = v_eta*geom%a(i)*F(i,1,1,k)-v_xi/(k-2)*this%E_i_M_N_K(geom, i, 1, 1, k-2)
            end do

            ! Calculate F(1,N,K) integrals
            do n=3,MXQ
                do k=3,MXK-2,2
                    F(i,1,n,k) = 2.*geom%a(i)*v_eta*F(i,1,n-1,k) &
                                 - (geom%a(i)**2+v_xi**2*geom%h**2)*F(i,1,n-2,k) &
                                 + v_xi**2*F(i,1,n-2,k-2)
                end do
            end do
        end do
        if (debug .and. any(isnan(F))) then
            write(*,*)
            write(*,*) "NaN found in F"
        end if

    end function panel_calc_F_integrals


    function panel_calc_H_integrals(this, geom, proc_H, MXK, MXQ, NHK, F) result(H)
        ! Calculates the necessary H integrals

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: proc_H, MXK, MXQ, NHK
        real,dimension(:,:,:,:),allocatable,intent(in) :: F
        real,dimension(:,:,:),allocatable :: H

        real :: S, C, nu
        integer :: i, m, n, k
        real,dimension(:),allocatable :: v_xi, v_eta

        ! Get edge normal derivatives
        allocate(v_xi(this%N), source=this%n_hat_l(:,1))
        allocate(v_eta(this%N), source=this%n_hat_l(:,2))

        ! Allocate integral storage
        allocate(H(1:MXQ,1:MXQ,1:MXK+NHK), source=0.)

        ! Procedure 1: not close to panel plane
        if (proc_H == 1) then

            ! Calculate H(1,1,1)
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
        else if (proc_H == 2 .or. proc_H == 3) then

            ! Initialize
            H(1,1,NHK+MXK) = 0.

            ! Calculate H(1,1,K) integrals
            do k=NHK+MXK,3,-2
                H(1,1,k-2) = 1./(k-4)*(geom%h**2*(k-2)*H(1,1,k)-sum(geom%a*F(:,1,1,k-2)))
            end do

        end if

        ! Calculate H(2,N,1) integrals
        do n=1,MXQ-1
            H(2,n,1) = 1./(n+1)*(geom%h**2*sum(v_xi*F(:,1,n,1)) &
                       + sum(geom%a*F(:,2,n,1)))
        end do

        ! Calculate H(1,N,1) integrals
        do n=2,MXQ
            if (n .eq. 2) then
                H(2,n,1) = 1./n*(geom%h**2*sum(v_eta*F(:,1,n-1,1)) &
                           + sum(geom%a*F(:,1,n,1)))
            else
                H(2,n,1) = 1./n*(-geom%h**2*(n-2)*H(1,n-2,1) & 
                           + geom%h**2*sum(v_eta*F(:,1,n-1,1)) &
                           + sum(geom%a*F(:,1,n,1)))
            end if
        end do

        ! Calculate H(M,N,1) integrals
        do m=3,MXQ
            do n=1,MXQ-m+1
                if (m .eq. 2) then
                    H(m,n,1) = 1./(m+n-1)*(geom%h**2*sum(v_xi*F(:,m-1,n,1)) &
                               + sum(geom%a*F(:,m,n,1)))
                else
                    H(m,n,1) = 1./(m+n-1)*(-geom%h**2*(m-2)*H(m-2,n,1) &
                               + geom%h**2*sum(v_xi*F(:,m-1,n,1)) &
                               + sum(geom%a*F(:,m,n,1)))
                end if
            end do
        end do

        ! Calculate H(1,N,K) integrals
        do n=2,MXQ
            do k=3,MXK,2
                if (n .eq. 2) then
                    H(1,n,k) = -1./(k-2)*sum(v_eta*F(:,1,n-1,k-2))
                else
                    H(1,n,k) = 1./(k-2)*((n-2)*H(1,n-2,k-2) &
                               -sum(v_eta*F(:,1,n-1,k-2)))
                end if
            end do
        end do

        ! Calculate H(2,N,K) integrals
        do n=1,MXQ-1
            do k=3,MXK,2
                H(2,n,k) = -1./(k-2)*sum(v_xi*F(:,1,n,k-2))
            end do
        end do

        ! Calculate remaining H(M,N,K) integrals
        do k=3,MXK,2
            do n=1,MXQ-m+1
                do m=3,MXQ
                    H(m,n,k) = -H(M-2,N+2,k)-geom%h**2*H(m-2,n,k)+H(m-2,n,k-2)
                end do
            end do
        end do

        ! Convert H* to H in case of Procedure 3
        if (proc_H .eq. 3) then
            do m=1,MXQ
                do n=1,MXQ
                    do k=1,MXK+NHK,2
                         
                        ! Get factorial dingas
                        nu = nu_M_N_K(m, n, k)

                        ! Convert H* to H
                        ! We need to make this check because h is sometimes zero, which can cause issues if the exponent is negative. If nu is zero, just don't bother.
                        if (abs(nu) <= 1e-12) then
                            H(m,n,k) = H(m,n,k)+2.*pi*nu*abs(geom%h)**(m+n-k)
                        end if
                    end do
                end do
            end do
        end if

        ! Clean up
        deallocate(v_xi)
        deallocate(v_eta)

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
        real,dimension(this%N) :: min_dist_to_edge
        integer :: MXQ, MXK, NHK, proc_H, i

        ! Determine which H integrals are needed based on distribution and type of influence
        if (singularity_type .eq. "source") then
            if (influence_type .eq. "potential") then
                if (source_order .eq. 0) then
                    MXQ = 1
                    MXK = 1
                end if
            else if (influence_type .eq. "velocity") then
                if (source_order .eq. 0) then
                    MXQ = 1
                    MXK = 3
                end if
            end if
        else if (singularity_type .eq. "doublet") then
            if (influence_type .eq. "potential") then
                if (doublet_order .eq. 1) then
                    MXQ = 2
                    MXK = 3
                end if
            else if (influence_type .eq. "velocity") then
                if (doublet_order .eq. 1) then
                    MXQ = 1
                    MXK = 3
                end if
            end if
        end if

        ! Calculate minimum distance to perimeter of S (in plane)
        do i=1,this%N

            ! Within edge, the minimum distance is the perpendicular distance
            if (geom%l1(i) < 0. .and. geom%l2(i) >= 0.) then
                min_dist_to_edge(i) = abs(geom%a(i))
        
            ! Otherwise, it is the minimum of the distances to the corners
            else
                min_dist_to_edge(i) = min(sqrt(geom%l1(i)**2+geom%a(i)**2), sqrt(geom%l2(i)**2+geom%a(i)**2))
            end if
        end do
        dH = minval(min_dist_to_edge)

        ! Determine which procedure needs to be used
        if (abs(geom%h) > 1e-12) then ! The nonzero h check seems to be more reliable than that proposed by Johnson
            proc_H = 1 ! Not near plane of panel
            NHK = 0

        else

            ! Check if the projected point falls inside the panel
            ! Outside panel (Procedure 2)
            if (all(geom%a < 0.)) then
                proc_H = 2

            ! Inside panel (Procedure 3)
            else
                proc_H = 3

            end if
            NHK = 16
        end if

        ! Calculate F integrals
        F = this%calc_F_integrals(geom, proc_H, MXK, MXQ, NHK)

        ! Calculate H integrals
        H = this%calc_H_integrals(geom, proc_H, MXK, MXQ, NHK, F)

    end subroutine panel_calc_integrals


    function panel_get_source_potential(this, eval_point, freestream, vertex_indices) result(phi)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        real,dimension(:),allocatable :: phi

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        if (influence_calc_type == 'johnson-ehlers') then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point)

            ! Get integrals
            call this%calc_integrals(geom, "potential", "source", H, F)

            if (source_order == 0) then

                ! Compute induced potential
                allocate(phi(1))
                phi = -0.25/pi*H(1,1,1)

            end if

            ! Clean up
            deallocate(H)
            deallocate(F)
        end if
    
    end function panel_get_source_potential


    function panel_get_source_velocity(this, eval_point, freestream, vertex_indices) result(v)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        real,dimension(:,:),allocatable :: v

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        if (influence_calc_type == 'johnson-ehlers') then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point)

            ! Get integrals
            call this%calc_integrals(geom, "velocity", "source", H, F)

            if (source_order .eq. 0) then

                ! Specify influencing vertices
                allocate(vertex_indices(1), source=0)

                ! Calculate velocity
                allocate(v(1,3))
                v(1,1) = 0.25/pi*sum(this%n_hat_l(:,1)*F(:,1,1,1))
                v(1,2) = 0.25/pi*sum(this%n_hat_l(:,2)*F(:,1,1,1))
                v(1,3) = 0.25/pi*geom%h*H(1,1,3)

            end if

            ! Clean up
            deallocate(H)
            deallocate(F)
        end if

    end function panel_get_source_velocity


    function panel_get_doublet_potential(this, eval_point, freestream, vertex_indices) result(phi)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        real,dimension(:),allocatable :: phi

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F
        integer,dimension(:),allocatable :: H_shape
        integer :: i, j, k

        if (influence_calc_type == 'johnson-ehlers') then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point)

            ! Get integrals
            call this%calc_integrals(geom, "potential", "doublet", H, F)

            ! Calculate influence
            if (doublet_order == 1) then

                ! Specify influencing vertices
                if (this%in_wake) then

                    ! Wake panels are influenced by two sets of vertices
                    allocate(vertex_indices(6))
                    allocate(phi(6), source=0.)
                    vertex_indices(1) = this%vertices(1)%ptr%top_parent
                    vertex_indices(2) = this%vertices(2)%ptr%top_parent
                    vertex_indices(3) = this%vertices(3)%ptr%top_parent
                    vertex_indices(4) = this%vertices(1)%ptr%bot_parent
                    vertex_indices(5) = this%vertices(2)%ptr%bot_parent
                    vertex_indices(6) = this%vertices(3)%ptr%bot_parent

                else

                    ! Body panels are influenced by only one set of vertices
                    allocate(vertex_indices, source=this%vertex_indices)
                    allocate(phi(3), source=0.)

                end if

                ! Compute induced potential
                phi(1) = geom%h*H(1,1,3)
                phi(2) = geom%r_in_plane(1)*geom%h*H(1,1,3)+geom%h*H(2,1,3)
                phi(3) = geom%r_in_plane(2)*geom%h*H(1,1,3)+geom%h*H(1,2,3)

                ! Convert to vertex influences
                phi(1:3) = 0.25/pi*matmul(phi(1:3), this%S_mu_inv)

                ! Wake bottom influence is opposite the top influence
                if (this%in_wake) then
                    phi(4:6) = -phi(1:3)
                end if

            end if

            ! Clean up
            deallocate(H)
            deallocate(F)

        end if
    
    end function panel_get_doublet_potential


    function panel_get_doublet_velocity(this, eval_point, freestream, vertex_indices) result(v)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        real,dimension(:,:),allocatable :: v

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        if (influence_calc_type == 'johnson-ehlers') then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point)

            ! Get integrals
            call this%calc_integrals(geom, "velocity", "doublet", H, F)

            if (doublet_order .eq. 1) then

                ! Specify influencing vertices
                allocate(vertex_indices(3), source=this%vertex_indices)

                ! Calculate velocity
                allocate(v(3,3))

            end if

            ! Clean up
            deallocate(H)
            deallocate(F)

        end if

    end function panel_get_doublet_velocity


    function panel_get_velocity_jump(this, mu, sigma, mirrored, mirror_plane) result(dv)
        ! Calculates the jump in perturbation velocity across this panel in global coordinates

        implicit none

        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: mu, sigma
        logical,intent(in) :: mirrored
        integer,intent(in) :: mirror_plane
        real,dimension(3) :: dv

        real,dimension(3) :: mu_verts, mu_params
        integer :: i

        if (doublet_order /= 1 .or. source_order /= 0) then
            write(*,*) "Velocity jump calculation has only been implemented for linear doublet and constant source distributions."
            stop
        end if

        ! Set up array of doublet strengths to calculate doublet parameters
        if (mirrored) then
            do i=1,this%N
                mu_verts(i) = mu(this%vertex_indices(i)+size(mu)/2)
            end do
        else
            do i=1,this%N
                mu_verts(i) = mu(this%vertex_indices(i))
            end do
        end if

        ! Calculate doublet parameters
        mu_params = matmul(this%S_mu_inv, mu_verts)

        ! Calculate velocity jump in panel coordinates
        dv(1) = mu_params(2)
        dv(2) = mu_params(3)
        if (mirrored) then
            dv(3) = sigma(this%index+size(sigma)/2)
        else
            dv(3) = sigma(this%index)
        end if

        ! Transform to global coordinates
        dv = matmul(transpose(this%A_g_to_l), dv)

        ! Mirror if necessary
        if (mirrored) then
            dv = mirror_about_plane(dv, mirror_plane)
        end if

    end function panel_get_velocity_jump
    
end module panel_mod