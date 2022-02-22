! Panel type
module panel_mod

    use helpers_mod
    use linked_list_mod
    use vertex_mod
    use math_mod
    use flow_mod
    use linalg_mod

    implicit none

    integer :: doublet_order
    integer :: source_order
    integer :: eval_count ! Developer counter for optimization purposes
    logical :: debug = .true. ! Developer toggle
    character(len=:),allocatable :: influence_calc_type ! Either 'johnson' or 'epton-magnus'


    type eval_point_geom
        ! Container type for the geometric parameters necessary for calculating a panel's influence on a given field point

        real,dimension(3) :: P_g ! Point position in global coords
        real,dimension(2) :: P_ls ! Transformed point in panel plane
        real :: h ! Transformed height above panel
        real,dimension(3) :: a, g2, l1, l2, R1, R2, g ! Edge parameters

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
        real,dimension(3) :: n_g, nu_g ! Normal and conormal vectors
        real,dimension(:,:),allocatable :: midpoints
        real,dimension(3) :: centroid
        real,dimension(3,3) :: A_g_to_l, A_s_to_ls, A_g_to_ls, A_ls_to_g ! Coordinate transformation matrices
        real,dimension(3,3) :: B_mat_ls, C_mat_ls ! Local scaled metric matrices
        real,dimension(:,:),allocatable :: vertices_ls, midpoints_ls ! Location of the vertices and edge midpoints described in local scaled coords
        real,dimension(:,:),allocatable :: t_hat_g, t_hat_l, t_hat_ls ! Edge unit tangents
        real,dimension(:,:),allocatable :: n_hat_g, n_hat_l, n_hat_ls ! Edge unit outward normals
        real,dimension(:,:),allocatable :: nu_hat_ls ! Edge unit outward conormals
        real,dimension(:),allocatable :: l ! Edge lengths
        real :: A ! Surface area
        real,dimension(:,:),allocatable :: S_mu_inv, S_sigma_inv ! Matrix relating doublet/source strengths to doublet/source influence parameters
        integer,dimension(:),allocatable :: vertex_indices ! Indices of this panel's vertices in the mesh vertex array
        logical :: in_wake ! Whether this panel belongs to a wake mesh
        integer,dimension(3) :: abutting_panels ! Indices of panels abutting this one
        integer,dimension(3) :: edges ! Indices of this panel's edges
        integer :: top_parent, bot_parent ! Indices of the top and bottom panels this panel's strength is determined by (for a wake panel w/ constant doublet strength)
        real :: r ! Panel inclination indicator; r=-1 -> superinclined, r=1 -> subinclined
        real,dimension(3) :: tau ! Edge inclination parameter
        integer,dimension(3) :: q ! Edge type indicator; q=1 -> subsonic, q=-1 -> supersonic
        real :: J ! Local scaled transformation Jacobian
        real :: rs ! Product of the inclination indicator and the flow type indicator
        real :: iota ! iota = |{n_g, n_g}|^(1/2); this quantity is often reused

        contains

            procedure :: init => panel_init_3
            procedure :: calc_derived_properties =>panel_calc_derived_properties
            procedure :: calc_area => panel_calc_area
            procedure :: calc_normal => panel_calc_normal
            procedure :: calc_centroid => panel_calc_centroid
            procedure :: calc_transforms => panel_calc_transforms
            procedure :: calc_g_to_l_transform => panel_calc_g_to_l_transform
            procedure :: calc_g_to_ls_transform => panel_calc_g_to_ls_transform
            procedure :: calc_edge_vectors => panel_calc_edge_vectors
            procedure :: calc_singularity_matrices => panel_calc_singularity_matrices
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
            procedure :: calc_source_potential => panel_calc_source_potential
            procedure :: calc_source_velocity => panel_calc_source_velocity
            procedure :: calc_doublet_potential => panel_calc_doublet_potential
            procedure :: calc_doublet_velocity => panel_calc_doublet_velocity
            procedure :: calc_potentials => panel_calc_potentials
            procedure :: calc_velocities => panel_calc_velocities
            procedure :: get_velocity_jump => panel_get_velocity_jump

    end type panel

    
contains


    subroutine panel_init_3(this, v1, v2, v3, index)
        ! Initializes a 3-panel

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: v1, v2, v3
        integer,intent(in) :: index

        ! Set number of sides
        this%N = 3

        ! Allocate vertex array
        allocate(this%vertices(this%N))
        allocate(this%vertex_indices(this%N))

        ! Assign vertex pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3

        ! Store vertex indices
        this%vertex_indices(1) = v1%index
        this%vertex_indices(2) = v2%index
        this%vertex_indices(3) = v3%index

        ! Store the index of the panel
        this%index = index

        ! Initialize a few things
        this%abutting_panels = 0
        this%top_parent = 0
        this%bot_parent = 0

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

            ! Check for zero area
            if (this%A < 1.e-12) then
                write(*,*) "!!! Panel", this%index, "has zero area. Quitting..."
                stop
            end if

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
        this%n_g = cross(d1, d2)
        this%n_g = this%n_g/norm(this%n_g)

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

        ! Calculate transforms
        call this%calc_g_to_l_transform(freestream)
        call this%calc_g_to_ls_transform(freestream)

        ! Calculate properties dependent on the transforms
        call this%calc_edge_vectors(freestream)
        call this%calc_singularity_matrices()

    end subroutine panel_calc_transforms


    subroutine panel_calc_g_to_l_transform(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        ! Calculate local eta axis
        this%A_g_to_l(2,:) = cross(this%n_g, freestream%c_hat_g)

        ! If the normal and freestream are aligned, then we have to pick a different vector.
        ! We pick the first edge, as we already know this is perpendicular to the normal, so
        ! orthogonality is automatically satisfied.
        if (norm(this%A_g_to_l(2,:)) < 1e-12) then
            this%A_g_to_l(2,:) = this%get_vertex_loc(2)-this%get_vertex_loc(1)
        end if

        ! Normalize
        this%A_g_to_l(2,:) = this%A_g_to_l(2,:)/norm(this%A_g_to_l(2,:))

        ! Calculate local xi axis
        this%A_g_to_l(1,:) = cross(this%A_g_to_l(2,:), this%n_g)
        this%A_g_to_l(1,:) = this%A_g_to_l(1,:)/norm(this%A_g_to_l(1,:))

        ! Store local zeta axis
        this%A_g_to_l(3,:) = this%n_g

    end subroutine panel_calc_g_to_l_transform


    subroutine panel_calc_g_to_ls_transform(this, freestream)
        ! Calculates the necessary transformations to move from global to local, scaled coordinates (Eq. (E.0.1) in Epton and Magnus)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(3) :: u0, v0
        real :: x, y
        integer :: i

        ! Get in-panel basis vectors
        u0 = this%A_g_to_l(1,:)
        v0 = this%A_g_to_l(2,:)

        ! Calculate compressible parameters
        this%nu_g = matmul(freestream%B_mat_g, this%n_g)
        x = inner(this%n_g, this%nu_g)

        ! Check for Mach-inclined panels
        if (freestream%supersonic .and. abs(x) < 1e-12) then
            write(*,*) "!!! Mach-inclined panels are not allowed in supersonic flow. Panel", &
                       this%index, "is Mach-inclined. Quitting..."
            stop
        end if

        ! Calculate panel inclination indicator (E&M Eq. (E.3.16b))
        this%r = sign(x) ! r=-1 -> superinclined, r=1 -> subinclined
        this%rs = this%r*freestream%s
        this%iota = sqrt(abs(x))

        ! Calculate transformation
        y = 1./this%iota
        this%A_g_to_ls(1,:) = y*matmul(freestream%C_mat_g, u0)
        this%A_g_to_ls(2,:) = this%rs/freestream%B*matmul(freestream%C_mat_g, v0)
        this%A_g_to_ls(3,:) = freestream%B*y*this%n_g

        ! Check determinant
        x = det3(this%A_g_to_ls)
        if (abs(x-freestream%B**2) >= 1e-12) then
            write(*,*) "!!! Calculation of local scaled coordinate transform failed. Quitting..."
            stop
        end if

        ! Calculate inverse
        if (freestream%M_inf == 0.) then
            this%A_ls_to_g = transpose(this%A_g_to_ls)
        else
            call matinv(3, this%A_g_to_ls, this%A_ls_to_g)
        end if

        ! Calculate Jacobian
        u0 = matmul(freestream%A_g_to_c, this%n_g)
        this%J = 1./(freestream%B*sqrt(abs(1.-freestream%M_inf**2*u0(1))))

        ! Transform vertex and midpoint coords to ls
        allocate(this%vertices_ls(this%N,2))
        allocate(this%midpoints_ls(this%N,2))
        do i=1,this%N

            ! Vertices
            this%vertices_ls(i,:) = matmul(this%A_g_to_l(1:2,:), this%get_vertex_loc(i)-this%centroid)

            ! Midpoints
            this%midpoints_ls(i,:) = matmul(this%A_g_to_l(1:2,:), this%midpoints(i,:)-this%centroid)

        end do

        ! Calculate local scaled metric matrices
        this%B_mat_ls = 0.
        this%B_mat_ls(1,1) = freestream%B**2*this%rs
        this%B_mat_ls(2,2) = freestream%B**2
        this%B_mat_ls(3,3) = freestream%B**2*this%r

        this%C_mat_ls = 0.
        this%C_mat_ls(1,1) = this%r
        this%C_mat_ls(2,2) = freestream%s
        this%C_mat_ls(3,3) = this%rs

        ! Check calculation (E&M Eq. (E.2.19))
        if (any(abs(this%B_mat_ls - matmul(this%A_g_to_ls, matmul(freestream%B_mat_g, transpose(this%A_g_to_ls)))) > 1e-12)) then
            write(*,*) "!!! Calculation of local scaled coordinate transform failed. Quitting..."
            stop
        end if
    
    end subroutine panel_calc_g_to_ls_transform


    subroutine panel_calc_edge_vectors(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(3) :: d
        real,dimension(2) :: e
        integer :: i

        ! Allocate memory
        allocate(this%l(this%N))
        allocate(this%t_hat_g(this%N,3))
        allocate(this%t_hat_l(this%N,2))
        allocate(this%t_hat_ls(this%N,2))
        allocate(this%n_hat_g(this%N,3))
        allocate(this%n_hat_l(this%N,2))
        allocate(this%n_hat_ls(this%N,2))
        allocate(this%nu_hat_ls(this%N,2))

        ! Loop through edges
        do i=1,this%N

            ! Calculate edge vector based on index
            d = this%get_vertex_loc(mod(i, this%N)+1)-this%get_vertex_loc(i)

            ! Calculate edge length (unscaled)
            this%l(i) = norm(d)

            ! Calculate tangent in global and local coords
            this%t_hat_g(i,:) = d/this%l(i)
            this%t_hat_l(i,:) = matmul(this%A_g_to_l(1:2,:), this%t_hat_g(i,:))

            ! Calculate tangent in local scaled coords
            e = matmul(this%A_g_to_ls(1:2,:), d)
            this%t_hat_ls(i,:) = e/sqrt(abs(this%rs*e(1)**2 + e(2)**2))

            ! Calculate edge outward normal
            this%n_hat_g(i,:) = cross(this%t_hat_g(i,:), this%n_g)
            this%n_hat_l(i,1) = this%t_hat_l(i,2)
            this%n_hat_l(i,2) = -this%t_hat_l(i,1)

            ! Calculate edge normal in local scaled coords
            this%n_hat_ls(i,1) = this%t_hat_ls(i,2)
            this%n_hat_ls(i,2) = -this%t_hat_ls(i,1)

            ! Calculate edge conormal
            this%nu_hat_ls(i,1) = this%rs*this%n_hat_ls(i,1)
            this%nu_hat_ls(i,2) = this%n_hat_ls(i,2)

            ! Calculate the edge type parameter (E&M Eq. (J.3.28) or Eq. (J.7.51))
            this%tau(i) = sqrt(abs(freestream%C_g_inner(this%t_hat_g(i,:), this%t_hat_g(i,:))))

            ! Calculate edge type indicator (E&M Eq. (J.6.48); this is the sign of b in Ehlers Eq. (E14))
            this%q(i) = sign(this%r*this%t_hat_ls(i,1)**2 + freestream%s*this%t_hat_ls(i,2)**2)

        end do
    
    end subroutine panel_calc_edge_vectors


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
                S_mu(:,2) = this%vertices_ls(:,1)
                S_mu(:,3) = this%vertices_ls(:,2)

                ! Invert
                call matinv(3, S_mu, this%S_mu_inv)

            else if (doublet_order .eq. 2) then

                ! Allocate influence matrix
                allocate(S_mu(6,6))
                allocate(this%S_mu_inv(6,6))

                ! Set values
                S_mu(:,1) = 1.

                S_mu(1:3,2) = this%vertices_ls(:,1)
                S_mu(1:3,3) = this%vertices_ls(:,2)
                S_mu(1:3,4) = this%vertices_ls(:,1)**2
                S_mu(1:3,5) = this%vertices_ls(:,1)*this%vertices_ls(:,2)
                S_mu(1:3,6) = this%vertices_ls(:,2)**2
                
                S_mu(4:6,2) = this%midpoints_ls(:,1)
                S_mu(4:6,3) = this%midpoints_ls(:,2)
                S_mu(4:6,4) = this%midpoints_ls(:,1)**2
                S_mu(4:6,5) = this%midpoints_ls(:,1)*this%midpoints_ls(:,2)
                S_mu(4:6,6) = this%midpoints_ls(:,2)**2

                ! Invert
                call matinv(6, S_mu, this%S_mu_inv)

            end if
            
            deallocate(S_mu)

        end if

        ! Determine influence of vertex source strengths on integral parameters
        ! Linear distribution
        if (source_order .eq. 1) then

            ! Allocate influence matrices
            allocate(S_sigma(3,3))
            allocate(this%S_sigma_inv(3,3))

            ! Set values
            S_sigma(:,1) = 1.
            S_sigma(:,2) = this%vertices_ls(:,1)
            S_sigma(:,3) = this%vertices_ls(:,2)

            ! Invert
            call matinv(3, S_sigma, this%S_sigma_inv)

            deallocate(S_sigma)

        else if (source_order .eq. 2) then

            ! Allocate influence matrix
            allocate(S_sigma(6,6))
            allocate(this%S_sigma_inv(6,6))

            ! Set values
            S_sigma(:,1) = 1.

            S_sigma(1:3,2) = this%vertices_ls(:,1)
            S_sigma(1:3,3) = this%vertices_ls(:,2)
            S_sigma(1:3,4) = this%vertices_ls(:,1)**2
            S_sigma(1:3,5) = this%vertices_ls(:,1)*this%vertices_ls(:,2)
            S_sigma(1:3,6) = this%vertices_ls(:,2)**2
            
            S_sigma(4:6,2) = this%midpoints_ls(:,1)
            S_sigma(4:6,3) = this%midpoints_ls(:,2)
            S_sigma(4:6,4) = this%midpoints_ls(:,1)**2
            S_sigma(4:6,5) = this%midpoints_ls(:,1)*this%midpoints_ls(:,2)
            S_sigma(4:6,6) = this%midpoints_ls(:,2)**2

            ! Invert
            call matinv(6, S_sigma, this%S_sigma_inv)

            deallocate(S_sigma)

        end if
    
    end subroutine panel_calc_singularity_matrices


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
            if (dist(this%get_vertex_loc(i), clone%loc) < 1e-12) then

                ! Update pointer
                this%vertices(i)%ptr => clone

                ! Update index
                this%vertex_indices(i) = clone%index

                return

            end if

        end do
    
    end subroutine panel_point_to_vertex_clone


    function panel_check_dod(this, eval_point, freestream, verts_in_dod, panel_mirrored, mirror_plane) result(dod_info)
        ! Determines how (if) this panel lies within the domain of dependence of the evaluation point

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        logical,dimension(:),intent(in) :: verts_in_dod
        logical,intent(in),optional :: panel_mirrored
        integer,intent(in),optional :: mirror_plane

        type(dod) :: dod_info

        real,dimension(3) :: d, point, a, b, R_star, Q_end
        integer :: i, i_next
        real :: x, s_star
        logical :: mirrored, in_panel

        ! Set default mirroring
        if (present(panel_mirrored)) then
            mirrored = panel_mirrored
        else
            mirrored = .false.
        end if

        ! First check the flow is supersonic (check_dod is always called)
        if (freestream%supersonic) then

            ! Read in vertex information
            do i=1,this%N
                if (mirrored) then
                    dod_info%verts_in_dod(i) = verts_in_dod(this%get_vertex_index(i)+size(verts_in_dod)/2)
                else
                    dod_info%verts_in_dod(i) = verts_in_dod(this%get_vertex_index(i))
                end if
            end do

            ! If all the vertices are in, then the panel is totally in and we can be done
            if (all(dod_info%verts_in_dod)) then
                dod_info%in_dod = .true.
                dod_info%edges_in_dod = .true.

            ! If it is not guaranteed to be totally in, then check all the edges
            else

                ! Check edges
                do i=1,this%N
                    
                    i_next = mod(i, this%N)+1

                    ! If if at least on endpoint is in, then the edge is in (you gotta love convex subspaces)
                    if (dod_info%verts_in_dod(i) .or. dod_info%verts_in_dod(i_next)) then
                        dod_info%edges_in_dod(i) = .true.

                    ! If both aren't in, then the intersection will depend on the edge type
                    else

                        ! For a subsonic edge, both being out means the edge is out
                        if (this%q(i) == 1) then
                            dod_info%edges_in_dod(i) = .false.

                        ! For a supersonic edge, the edge can still intersect the DoD, so calculate the point of closest approach
                        else

                            ! Get end vertex and vector describing edge
                            if (mirrored) then
                                Q_end = mirror_about_plane(this%get_vertex_loc(i_next), mirror_plane)
                                d = Q_end - mirror_about_plane(this%get_vertex_loc(i), mirror_plane)
                            else
                                Q_end = this%get_vertex_loc(i_next)
                                d = Q_end - this%get_vertex_loc(i)
                            end if
                        
                            ! Calculate nondimensional location of the point of closest approach (E&M Eq. (J.3.39))
                            a = cross(freestream%c_hat_g, d)
                            b = cross(freestream%c_hat_g, Q_end - eval_point)
                            s_star = inner(a, b)/abs(inner(a, a))

                            ! Calculate point of closest approach
                            R_star = Q_end - s_star*d

                            ! Check if the point of closest approach is in the edge and in the DoD
                            if (s_star > 0. .and. s_star < 1. .and. freestream%point_in_dod(R_star, eval_point)) then
                                dod_info%edges_in_dod(i) = .true.

                            ! If not, this edge is not in the DoD
                            else
                                dod_info%edges_in_dod(i) = .false.

                            end if
                        end if
                    end if

                end do

                ! If any edge or vertex is in the DoD, then the panel is in
                if (any(dod_info%verts_in_dod) .or. any(dod_info%edges_in_dod)) then
                    dod_info%in_dod = .true.

                ! If a supersonic panel has no edges or vertices in the DoD, check if the DoD is encompassed by the panel
                else if (this%r == -1) then

                    ! Get the projection of the evaluation point onto the panel in the direction of c_hat
                    s_star = inner(this%get_vertex_loc(1)-eval_point, this%n_g)/inner(freestream%c_hat_g, this%n_g)
                    R_star = eval_point + freestream%c_hat_g*s_star

                    ! See if the projected point is in the panel
                    in_panel = .true.
                    do i=1,this%N

                        ! Get edge displacement
                        x = inner(R_star, this%n_hat_g(:,i))

                        ! Check sign (should be negative if interior to the panel)
                        if (x >= 0.) then
                            in_panel = .false.
                            exit ! Don't need to check any more
                        end if

                    end do

                    ! Store information
                    dod_info%in_dod = in_panel

                ! Not supersonic and no edges or vertices in. Not in.
                else
                    dod_info%in_dod = .false.

                end if
            end if

        else

            ! Subsonic flow. DoD is everywhere. Life is easy.
            dod_info%in_dod = .true.
            dod_info%verts_in_dod = .true.
            dod_info%edges_in_dod = .true.

        end if
    
    end function panel_check_dod


    function panel_get_field_point_geometry(this, eval_point, freestream, dod_info) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        type(eval_point_geom) :: geom

        real,dimension(2) :: d_ls
        real,dimension(3) :: d_g, x
        integer :: i
        real :: val

        ! Store point
        geom%P_g = eval_point

        ! Transform to local scaled coordinates
        x = matmul(this%A_g_to_ls, geom%P_g-this%centroid)
        geom%P_ls = x(1:2)
        geom%h = x(3) ! Equivalent to E&M Eq. (J.7.41)

        ! Calculate intermediate quantities
        do i=1,this%N
            if (dod_info%edges_in_dod(i)) then

                ! Calculate displacements
                d_ls = this%vertices_ls(i,:) - geom%P_ls
                d_g = this%get_vertex_loc(i) - geom%P_g

                ! Get vector perpendicular to edge
                x = cross(d_g, this%t_hat_g(i,:))

                ! Perpendicular distance in plane from evaluation point to edge E&M Eq. (J.6.46) and (J.7.53)
                geom%a(i) = inner2(d_ls, this%n_hat_ls(i,:))
                val = this%r*freestream%B/(this%tau(i)*this%iota) * freestream%B_g_inner(this%n_hat_g(i,:), x)
                !write(*,*) val-geom%a(i)

                ! Integration length on edge to start vertex (E&M Eq. (J.6.47))
                geom%l1(i) = this%rs*d_ls(1)*this%t_hat_ls(i,1) + d_ls(2)*this%t_hat_ls(i,2)

                ! Projected point displacement from end vertex
                d_ls = this%vertices_ls(mod(i, this%N)+1,:)-geom%P_ls

                ! Integration length on edge to end vertex
                geom%l2(i) = this%rs*d_ls(1)*this%t_hat_ls(i,1) + d_ls(2)*this%t_hat_ls(i,2)

                ! Distance from evaluation point to start vertex E&M Eq. (J.8.8)
                ! The distance should be zero in the case of a negative squared distance, as the flow is supersonic and the point lies outside the DoD
                val = freestream%C_g_inner(-d_g, -d_g)
                if (val > 0.) then
                    geom%R1(i) = sqrt(val)
                else
                    geom%R1(i) = 0.
                end if

                ! Calculate square of the perpendicular distance to edge
                geom%g2(i) = (freestream%B/this%tau(i))**2*freestream%B_g_inner(x, x) ! E&M Eq. (J.8.23) or (J.7.70)

                ! Calculate the perpendicular distance to edge
                if (geom%g2(i) >= 0.) then
                    geom%g(i) = sqrt(geom%g2(i))
                else
                    geom%g(i) = 0.
                end if

            else

                geom%R1(i) = 0.
                geom%a(i) = 0.

            end if
        end do

        ! Distance from evaluation point to end vertices
        geom%R2 = cshift(geom%R1, 1)

    end function panel_get_field_point_geometry


    function panel_E_i_M_N_K(this, geom, i, M, N, K, dod_info, freestream) result(E)
        ! Calculates E_i(M,N,K)

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: i, M, N, K
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream

        real :: E

        real :: E1, E2
        integer :: i_next

        i_next = mod(i, this%N) + 1

        ! Evaluate at start vertex
        if (dod_info%verts_in_dod(i)) then
            E1 = ((this%vertices_ls(i,1)-geom%P_ls(1))**(M-1)*(this%vertices_ls(i,2)-geom%P_ls(2))**(N-1))/geom%R1(i)**K
        else
            E1 = 0.
        end if

        ! Evaluate at end vertex
        if (dod_info%verts_in_dod(i_next)) then
            E2 = ((this%vertices_ls(i_next,1)-geom%P_ls(1))**(M-1)*(this%vertices_ls(i_next,2)-geom%P_ls(2))**(N-1))/geom%R2(i)**K
        else
            E2 = 0.
        end if

        ! Calculate difference
        E = E2-E1

    end function panel_E_i_M_N_K


    function panel_F_i_1_1_1(this, geom, i, dod_info, freestream) result(F)
        ! Calculates F_i(1,1,1)

        implicit none

        class(panel) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: i
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream

        real :: F, F1, F2
        integer :: i_next

        i_next = mod(i, this%N) + 1

        ! Check DoD
        if (dod_info%edges_in_dod(i)) then

            ! Supersonic edge
            if (this%q(i) == -1) then

                ! Neither endpoint in DoD (Ehlers Eq. (E22))
                if (geom%R1(i) == 0. .and. geom%R2(i) == 0.) then
                    F = pi

                ! At least one endpoint in the DoD (Ehlers Eq. (E22))
                else

                    ! Calculate preliminary quantities
                    F1 = geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)
                    F2 = geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)

                    ! Calculate F
                    F = -atan2(F1, F2)

                end if

            ! Subsonic edge (all edges in subsonic flow because the definition of a subsonic edge is one that lies inside its own DoD)
            else

                ! Within edge (Johnson Eq. (D.60))
                if (sign(geom%l1(i)) /= sign(geom%l2(i))) then
                    F = log(((geom%R1(i)-geom%l1(i))*(geom%R2(i)+geom%l2(i)))/geom%g2(i))

                ! Above or below edge; this is a unified form of Johnson Eq. (D.60) and should be equivalent to Ehlers Eq. (E.22)
                else
                    F = sign(geom%l1(i))*log((geom%R2(i) + abs(geom%l2(i))) / (geom%R1(i) + abs(geom%l1(i))))

                end if
            end if

        ! Not in DoD
        else
            F = 0.
        end if
        
    end function panel_F_i_1_1_1


    function panel_calc_F_integrals(this, geom, proc_H, MXK, MXQ, NHK, dod_info, freestream) result (F)
        ! Calculates the F integrals necessary

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: proc_H, MXK, MXQ, NHK
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        real,dimension(:,:,:,:),allocatable :: F

        real :: E1, E2, v_xi, v_eta, dF
        real,dimension(3) :: d
        real,dimension(this%N) :: min_dist_to_edge
        integer :: i, MXFK, NFK, k, m, n, i_next

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

            ! Check this edge is in the DoD
            if (dod_info%edges_in_dod(i)) then

                i_next = mod(i, this%N) + 1

                ! Within edge, the minimum distance is the perpendicular distance
                if (sign(geom%l1(i)) /= sign(geom%l2(i))) then
                    min_dist_to_edge(i) = geom%g(i)
        
                ! Otherwise, it is the minimum of the distances to the corners
                else
                    if (dod_info%verts_in_dod(i) .and. dod_info%verts_in_dod(i_next)) then
                        min_dist_to_edge(i) = min(geom%R1(i), geom%R2(i))
                    else if (dod_info%verts_in_dod(i)) then
                        min_dist_to_edge(i) = geom%R1(i)
                    else if (dod_info%verts_in_dod(i_next)) then
                        min_dist_to_edge(i) = geom%R2(i)

                    ! If neither is in the DoD, we again go back to the perpendicular distance
                    else
                        min_dist_to_edge(i) = geom%g(i)
                    end if
                end if

            ! Otherwise, set the minimum distance arbitrarily high
            else
                min_dist_to_edge(i) = 100000.
            end if

        end do
        dF = minval(min_dist_to_edge)

        ! Check for point on perimeter
        if (abs(dF) < 1e-12) then
            write(*,*) "Detected control point on perimeter of panel. Quitting..."
            stop
        end if

        ! Loop through edges
        do i=1,this%N

            if (dod_info%edges_in_dod(i)) then

                ! Store edge derivs
                v_xi = this%n_hat_ls(i,1)
                v_eta = this%n_hat_ls(i,2)

                ! Calculate F(1,1,1)
                F(i,1,1,1) = this%F_i_1_1_1(geom, i, dod_info, freestream)

                ! Procedure 4: not close to perimeter
                if (geom%g(i) >= 0.01*dF) then

                    ! Calculate F(1,1,K) integrals
                    do k=3,MXFK,2

                        ! Get necessary E
                        E1 = this%E_i_M_N_K(geom, i, 2, 1, k-2, dod_info, freestream)
                        E2 = this%E_i_M_N_K(geom, i, 1, 2, k-2, dod_info, freestream)

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
                        E1 = this%E_i_M_N_K(geom, i, 2, 1, k-2, dod_info, freestream)
                        E2 = this%E_i_M_N_K(geom, i, 1, 2, k-2, dod_info, freestream)

                        ! Calculate F
                        F(i,1,1,k-2) = 1./(k-3)*(geom%g2(i)*(k-2)*F(i,1,1,k)+v_eta*E1-v_xi*E2)

                    end do
                end if

                ! Calculate other F integrals (same for both procedures)
                if (abs(v_eta) <= abs(v_xi)) then ! Case a

                    ! Calculate F(1,N,1) integrals (Johnson Eq. (D.62))
                    do n=2,MXQ

                        ! Get E
                        E1 = this%E_i_M_N_K(geom, i, 1, n-1, -1, dod_info, freestream)

                        if (n .eq. 2) then
                            F(i,1,N,1) = 1./(n-1)*((2*n-3)*geom%a(i)*v_eta*F(i,1,n-1,1) + v_xi*E1)
                        else
                            F(i,1,N,1) = 1./(n-1)*((2*n-3)*geom%a(i)*v_eta*F(i,1,n-1,1) &
                                         - (n-2)*(geom%a(i)**2+v_xi**2*geom%h**2)*F(i,1,n-2,1) &
                                         + v_xi*E1)
                        end if
                    end do

                    ! Calculate F(M,N,1) inegrals (Johnsons Eq. (D.63))
                    do m=2,MXQ
                        do n=1,MXQ-m+1
                            F(i,m,n,1) = -v_eta/v_xi*F(i,m-1,n+1,1) + geom%a(i)/v_xi*F(i,m-1,n,1)
                        end do
                    end do

                else ! Case b

                    ! Calculate F(M,1,1) integrals (Johnson Eq. (D.64))
                    do m=2,MXQ

                        ! Get E
                        E1 = this%E_i_M_N_K(geom, i, m-1, 1, -1, dod_info, freestream)

                        if (m .eq. 2) then
                            F(i,m,1,1) = 1./(m-1)*((2*m-3)*geom%a(i)*v_xi*F(i,m-1,1,1) - v_eta*E1)
                        else
                            F(i,m,1,1) = 1./(m-1)*((2*m-3)*geom%a(i)*v_xi*F(i,m-1,1,1) &
                                         - (m-2)*(geom%a(i)**2+v_eta**2*geom%h**2)*F(i,m-2,1,1) &
                                         - v_eta*E1)
                        end if
                    end do

                    ! Calculate F(M,N,1) integrals (Johnson Eq. (D.65))
                    do n=2,MXQ
                        do m=1,MXQ-n+1
                            F(i,m,n,1) = -v_xi/v_eta*F(i,m+1,n-1,1)+geom%a(i)/v_eta*F(i,m,n-1,1)
                        end do
                    end do

                end if

                ! Calculate F(1,2,K) integrals (Johnson Eq. (D.66))
                do k=3,MXK-2,2
                    F(i,1,2,k) = v_eta*geom%a(i)*F(i,1,1,k)-v_xi/(k-2)*this%E_i_M_N_K(geom, i, 1, 1, k-2, dod_info, freestream)
                end do

                ! Calculate F(1,N,K) integrals (Johnson Eq. (D.67))
                do n=3,MXQ
                    do k=3,MXK-2,2
                        F(i,1,n,k) = 2.*geom%a(i)*v_eta*F(i,1,n-1,k) &
                                     - (geom%a(i)**2+v_xi**2*geom%h**2)*F(i,1,n-2,k) &
                                     + v_xi**2*F(i,1,n-2,k-2)
                    end do
                end do
            end if
        end do

        if (debug .and. any(isnan(F))) then
            write(*,*)
            write(*,*) "NaN found in F"
        end if

    end function panel_calc_F_integrals


    function panel_calc_H_integrals(this, geom, proc_H, MXK, MXQ, NHK, dod_info, freestream, F) result(H)
        ! Calculates the necessary H integrals

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: proc_H, MXK, MXQ, NHK
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        real,dimension(:,:,:,:),allocatable,intent(in) :: F
        real,dimension(:,:,:),allocatable :: H

        real :: S, C, nu, c1, c2, F1, F2
        integer :: i, m, n, k
        real,dimension(:),allocatable :: v_xi, v_eta

        ! Get edge normal derivatives
        allocate(v_xi(this%N), source=this%n_hat_ls(:,1))
        allocate(v_eta(this%N), source=this%n_hat_ls(:,2))

        ! Allocate integral storage
        allocate(H(1:MXQ,1:MXQ,1:MXK+NHK), source=0.)

        ! Procedure 1: not close to panel plane
        if (proc_H == 1) then

            ! Calculate H(1,1,1) (Johnson Eq. (D.41))
            do i=1,this%N

                if (dod_info%edges_in_dod(i)) then

                    ! Supersonic edge
                    if (this%q(i) == -1) then

                        ! Check for neither endpoint in
                        if (.not. dod_info%verts_in_dod(i) .and. .not. dod_info%verts_in_dod(mod(i, this%N)+1)) then

                            ! Add influence of this edge
                            H(1,1,1) = H(1,1,1) - abs(geom%h)*sign(v_xi(i))*pi ! Adapted from Ehlers Eq. (E18) to match the form of Johnson Eq. (D.41)

                        ! At least one in
                        else

                            ! Calculate intermediate quantities Ehlers Eq. (E19) and (E20)
                            F1 = (geom%l1(i)*geom%R2(i)-geom%l2(i)*geom%R1(i))/geom%g2(i)
                            F2 = (geom%R1(i)*geom%R2(i)+geom%l1(i)*geom%l2(i))/geom%g2(i)

                            ! Add to surface integral; adapted from Ehlers Eq. (E18) to match the form of Johnson Eq. (D.41)
                            H(1,1,1) = H(1,1,1) - geom%h*atan2(geom%h*geom%a(i)*F1, geom%R1(i)*geom%R2(i)+geom%h**2*F2)

                        end if

                    ! Subsonic edge
                    else

                        ! Calculate intermediate quantities (Johnson Eq. (D.41))
                        c1 = geom%g2(i)+abs(geom%h)*geom%R1(i)
                        c2 = geom%g2(i)+abs(geom%h)*geom%R2(i)
        
                        ! Add surface integral
                        S = geom%a(i)*(geom%l2(i)*c1 - geom%l1(i)*c2)
                        C = c1*c2 + geom%a(i)**2*geom%l1(i)*geom%l2(i)
                        H(1,1,1) = H(1,1,1) - atan2(S, C)

                    end if

                end if
        
            end do
        
            ! Add line integrals (rs factor added to match Ehlers Eq. (E9))
            H(1,1,1) = abs(geom%h)*H(1,1,1) + this%rs*sum(geom%a*F(:,1,1,1))

            ! Calculate H(1,1,K) integrals (Johnson Eq. (D.42) altered to match Ehlers Eq. (E9))
            do k=3,MXK,2
                if (this%rs == 1.) then
                    H(1,1,k) = 1./((k-2)*geom%h**2)*((k-4)*H(1,1,k-2) + sum(geom%a*F(:,1,1,k-2)))
                else
                    H(1,1,k) = 1./((k-2)*geom%h**2)*(-k*H(1,1,k-2) - sum(geom%a*F(:,1,1,k-2)))
                end if
            end do

        ! Procedures 2 and 3: close to panel plane
        else if (proc_H == 2 .or. proc_H == 3) then

            ! Initialize
            H(1,1,NHK+MXK) = 0.

            ! Calculate H(1,1,K) integrals (Johnson Eq. (D.50) altered to match Ehlers Eq. (E9))
            do k=NHK+MXK,3,-2
                H(1,1,k-2) = 1./(k-4)*(geom%h**2*(k-2)*H(1,1,k) - this%rs*sum(geom%a*F(:,1,1,k-2)))
            end do

        end if

        ! Step 3
        ! Calculate H(2,N,1) integrals (Johnson Eq. (D.43) altered to match Ehlers Eq. (E7))
        do n=1,MXQ-1
            H(2,n,1) = 1./(this%rs*n+1)*(this%rs*geom%h**2*sum(v_xi*F(:,1,n,1)) + sum(geom%a*F(:,2,n,1)))
        end do

        ! Step 4
        ! Calculate H(1,N,1) integrals (Johnson Eq. (D.44) altered to match Ehlers Eq. (E8))
        do n=2,MXQ
            if (n .eq. 2) then
                H(2,n,1) = this%rs/n*(geom%h**2*sum(v_eta*F(:,1,n-1,1)) &
                           + sum(geom%a*F(:,1,n,1)))
            else
                H(2,n,1) = this%rs/n*(-geom%h**2*(n-2)*H(1,n-2,1) & 
                           + geom%h**2*sum(v_eta*F(:,1,n-1,1)) &
                           + sum(geom%a*F(:,1,n,1)))
            end if
        end do

        ! Step 5
        ! Calculate H(M,N,1) integrals (Johnson Eq. (D.45) altered to match Ehlers Eq. (E7))
        do m=3,MXQ
            do n=1,MXQ-m+1
                if (m .eq. 2) then
                    H(m,n,1) = 1./(m+n-1)*(this%rs*geom%h**2*sum(v_xi*F(:,m-1,n,1)) &
                               + sum(geom%a*F(:,m,n,1)))
                else
                    H(m,n,1) = 1./(m+n-1)*(-geom%h**2*(m-2)*H(m-2,n,1) &
                               + this%rs*geom%h**2*sum(v_xi*F(:,m-1,n,1)) &
                               + sum(geom%a*F(:,m,n,1)))
                end if
            end do
        end do

        ! Step 6
        ! Calculate H(1,N,K) integrals (Johnson Eq. (D.46) altered to match Ehlers Eq. (E6))
        do n=2,MXQ
            do k=3,MXK,2
                if (n .eq. 2) then
                    H(1,n,k) = -this%rs/(k-2)*sum(v_eta*F(:,1,n-1,k-2))
                else
                    H(1,n,k) = this%rs/(k-2)*((n-2)*H(1,n-2,k-2) &
                               -sum(v_eta*F(:,1,n-1,k-2)))
                end if
            end do
        end do

        ! Step 7
        ! Calculate H(2,N,K) integrals (Johnson Eq. (D.47); already matches Ehlers Eq. (E5))
        do n=1,MXQ-1
            do k=3,MXK,2
                H(2,n,k) = -1./(k-2)*sum(v_xi*F(:,1,n,k-2))
            end do
        end do

        ! Step 8
        ! Calculate remaining H(M,N,K) integrals (Johnson Eq. (D.48) altered to match Ehlers Eq. (E4))
        do k=3,MXK,2
            do n=1,MXQ-m+1
                do m=3,MXQ
                    H(m,n,k) = this%rs*(-H(M-2,N+2,k)-geom%h**2*H(m-2,n,k))+H(m-2,n,k-2)
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
                        ! If nu or h is zero, this can just be skipped.
                        if (abs(geom%h) > 1e-12 .and. abs(nu) > 1e-12) then
                            H(m,n,k) = H(m,n,k) + 2.*pi*nu*abs(geom%h)**(m+n-k)
                        end if
                    end do
                end do
            end do
        end if

        ! Clean up
        deallocate(v_xi)
        deallocate(v_eta)

        if (debug .and. any(isnan(H))) then
            write(*,*)
            write(*,*) "NaN found in H"
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


    subroutine panel_calc_integrals(this, geom, influence_type, singularity_type, dod_info, freestream, H, F)
        ! Calculates the H and F integrals necessary for the given influence

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        character(len=*),intent(in) :: influence_type, singularity_type
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
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
        F = this%calc_F_integrals(geom, proc_H, MXK, MXQ, NHK, dod_info, freestream)

        ! Calculate H integrals
        H = this%calc_H_integrals(geom, proc_H, MXK, MXQ, NHK, dod_info, freestream, F)

    end subroutine panel_calc_integrals


    subroutine panel_calc_potentials(this, P, freestream, dod_info, mirror_panel, phi_s, phi_d, i_vert_s, i_vert_d)
        ! Calculates the source- and doublet-induced potentials at the given point P

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: P
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        logical,intent(in) :: mirror_panel
        real,dimension(:),allocatable,intent(out) :: phi_s, phi_d
        integer,dimension(:),allocatable,intent(out) :: i_vert_s, i_vert_d

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        ! Specify influencing vertices (also sets zero default influence)

        ! Source
        if (source_order == 0) then
            allocate(phi_s(1), source=0.)
        end if

        ! Doublet
        if (doublet_order == 1) then

            ! Check if this panel belongs to the wake
            if (this%in_wake) then

                ! Wake panels are influenced by two sets of vertices
                allocate(i_vert_d(6))
                i_vert_d(1) = this%vertices(1)%ptr%top_parent
                i_vert_d(2) = this%vertices(2)%ptr%top_parent
                i_vert_d(3) = this%vertices(3)%ptr%top_parent
                i_vert_d(4) = this%vertices(1)%ptr%bot_parent
                i_vert_d(5) = this%vertices(2)%ptr%bot_parent
                i_vert_d(6) = this%vertices(3)%ptr%bot_parent

                ! Set default influence
                allocate(phi_d(6), source=0.)

            else

                ! Body panels are influenced by only one set of vertices
                allocate(i_vert_d, source=this%vertex_indices)

                ! Set default influence
                allocate(phi_d(3), source=0.)

            end if
        end if

        ! Check DoD
        if (dod_info%in_dod) then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(P, freestream, dod_info)

            if (influence_calc_type == 'johnson') then

                ! Get integrals
                call this%calc_integrals(geom, "potential", "doublet", dod_info, freestream, H, F)

                ! Source potential
                if (source_order == 0) then
                    phi_s = -this%J*freestream%K_inv*H(1,1,1)
                end if

                ! Doublet potential
                if (doublet_order == 1) then

                    ! Compute induced potential (Johnson Eq. (D.30); Ehlers Eq. (5.17))
                    phi_d(1) = geom%h*H(1,1,3)
                    phi_d(2) = geom%h*H(1,1,3)*geom%P_ls(1) + geom%h*H(2,1,3)
                    phi_d(3) = geom%h*H(1,1,3)*geom%P_ls(2) + geom%h*H(1,2,3)

                    ! Convert to vertex influences (Davis Eq. (4.41))
                    phi_d(1:3) = freestream%K_inv*matmul(phi_d(1:3), this%S_mu_inv)

                    ! Wake bottom influence is opposite the top influence
                    if (this%in_wake) then
                        phi_d(4:6) = -phi_d(1:3)
                    end if
                end if

                ! Clean up
                deallocate(H)
                deallocate(F)

            end if
        end if
    
    end subroutine panel_calc_potentials


    subroutine panel_calc_velocities(this, P, freestream, dod_info, mirror_panel, v_s, v_d, i_vert_s, i_vert_d)
        ! Calculates the source- and doublet-induced potentials at the given point P

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: P
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        logical,intent(in) :: mirror_panel
        real,dimension(:,:),allocatable,intent(out) :: v_s, v_d
        integer,dimension(:),allocatable,intent(out) :: i_vert_s, i_vert_d

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        ! Specify influencing vertices (also sets zero default influence)

        ! Source
        if (source_order == 0) then
            allocate(v_s(1,3), source=0.)
        end if

        ! Doublet
        if (doublet_order == 1) then

            ! Check if this panel belongs to the wake
            if (this%in_wake) then

                ! Wake panels are influenced by two sets of vertices
                allocate(i_vert_d(6))
                i_vert_d(1) = this%vertices(1)%ptr%top_parent
                i_vert_d(2) = this%vertices(2)%ptr%top_parent
                i_vert_d(3) = this%vertices(3)%ptr%top_parent
                i_vert_d(4) = this%vertices(1)%ptr%bot_parent
                i_vert_d(5) = this%vertices(2)%ptr%bot_parent
                i_vert_d(6) = this%vertices(3)%ptr%bot_parent

                ! Set default influence
                allocate(v_d(6,3), source=0.)

            else

                ! Body panels are influenced by only one set of vertices
                allocate(i_vert_d, source=this%vertex_indices)

                ! Set default influence
                allocate(v_d(3,3), source=0.)

            end if
        end if

        ! Check DoD
        if (dod_info%in_dod) then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(P, freestream, dod_info)

            if (influence_calc_type == 'johnson') then

                ! Get integrals
                call this%calc_integrals(geom, "velocity", "doublet", dod_info, freestream, H, F)

                ! Source velocity
                if (source_order == 0) then
                end if

                ! Doublet velocity
                if (doublet_order == 1) then
                end if

                ! Clean up
                deallocate(H)
                deallocate(F)

            end if
        end if
    
    end subroutine panel_calc_velocities


    function panel_calc_source_potential(this, eval_point, freestream, dod_info, vertex_indices, panel_mirrored) result(phi)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        logical,intent(in) :: panel_mirrored
        real,dimension(:),allocatable :: phi

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        ! Allocate potential
        if (source_order == 0) then
            allocate(phi(1), source=0.)
        end if

        ! In dod
        if (dod_info%in_dod) then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point, freestream, dod_info)

            if (influence_calc_type == 'johnson') then

                ! Get integrals
                call this%calc_integrals(geom, "potential", "source", dod_info, freestream, H, F)

                if (source_order == 0) then

                    ! Compute induced potential
                    phi = -this%J*freestream%K_inv*H(1,1,1)

                end if

                ! Clean up
                deallocate(H)
                deallocate(F)

            end if
        end if
    
    end function panel_calc_source_potential


    function panel_calc_source_velocity(this, eval_point, freestream, dod_info, vertex_indices, panel_mirrored) result(v)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        logical,intent(in) :: panel_mirrored
        real,dimension(:,:),allocatable :: v

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        ! Allocate velocity
        if (source_order == 0) then
            allocate(v(1,3), source=0.)
        end if

        ! In dod
        if (dod_info%in_dod) then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point, freestream, dod_info)

            if (influence_calc_type == 'johnson') then

                ! Get integrals
                call this%calc_integrals(geom, "velocity", "source", dod_info, freestream, H, F)

                if (source_order == 0) then

                    ! Calculate velocity
                    v(1,1) = this%J*freestream%K_inv*sum(this%n_hat_ls(:,1)*F(:,1,1,1))
                    v(1,2) = this%J*freestream%K_inv*sum(this%n_hat_ls(:,2)*F(:,1,1,1))
                    v(1,3) = this%J*freestream%K_inv*geom%h*H(1,1,3)

                end if

                ! Clean up
                deallocate(H)
                deallocate(F)

            end if
        end if

    end function panel_calc_source_velocity


    function panel_calc_doublet_potential(this, eval_point, freestream, dod_info, vertex_indices, panel_mirrored) result(phi)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        logical,intent(in) :: panel_mirrored
        real,dimension(:),allocatable :: phi

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F
        integer,dimension(:),allocatable :: H_shape
        integer :: i, j, k
        real :: a
        real,dimension(2) :: a_bar, x

        ! Specify influencing vertices (also sets zero default influence)
        if (doublet_order == 1) then
            if (this%in_wake) then

                ! Wake panels are influenced by two sets of vertices
                allocate(vertex_indices(6))
                vertex_indices(1) = this%vertices(1)%ptr%top_parent
                vertex_indices(2) = this%vertices(2)%ptr%top_parent
                vertex_indices(3) = this%vertices(3)%ptr%top_parent
                vertex_indices(4) = this%vertices(1)%ptr%bot_parent
                vertex_indices(5) = this%vertices(2)%ptr%bot_parent
                vertex_indices(6) = this%vertices(3)%ptr%bot_parent

                ! Set default influence
                allocate(phi(6), source=0.)

            else

                ! Body panels are influenced by only one set of vertices
                allocate(vertex_indices, source=this%vertex_indices)

                ! Set default influence
                allocate(phi(3), source=0.)

            end if
        end if

        ! Check DoD
        if (dod_info%in_dod) then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point, freestream, dod_info)

            if (influence_calc_type == 'johnson') then

                ! Get integrals
                call this%calc_integrals(geom, "potential", "doublet", dod_info, freestream, H, F)

                ! Calculate influence
                if (doublet_order == 1) then

                    ! Compute induced potential (Johnson Eq. (D.30); Ehlers Eq. (5.17))
                    phi(1) = geom%h*H(1,1,3)
                    phi(2) = geom%h*H(1,1,3)*geom%P_ls(1) + geom%h*H(2,1,3)
                    phi(3) = geom%h*H(1,1,3)*geom%P_ls(2) + geom%h*H(1,2,3)

                    ! Convert to vertex influences (Davis Eq. (4.41))
                    phi(1:3) = freestream%K_inv*matmul(phi(1:3), this%S_mu_inv)

                    ! Wake bottom influence is opposite the top influence
                    if (this%in_wake) then
                        phi(4:6) = -phi(1:3)
                    end if

                end if

                ! Clean up
                deallocate(H)
                deallocate(F)

            end if
        end if
    
    end function panel_calc_doublet_potential


    function panel_calc_doublet_velocity(this, eval_point, freestream, dod_info, vertex_indices, panel_mirrored) result(v)

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        integer,dimension(:),allocatable,intent(out) :: vertex_indices
        logical,intent(in) :: panel_mirrored
        real,dimension(:,:),allocatable :: v

        type(eval_point_geom) :: geom
        real,dimension(:,:,:),allocatable :: H
        real,dimension(:,:,:,:),allocatable :: F

        ! Allocate velocity and vertices
        if (doublet_order .eq. 1) then
            if (this%in_wake) then

                ! Wake panels are influenced by two sets of vertices
                allocate(vertex_indices(6))
                vertex_indices(1) = this%vertices(1)%ptr%top_parent
                vertex_indices(2) = this%vertices(2)%ptr%top_parent
                vertex_indices(3) = this%vertices(3)%ptr%top_parent
                vertex_indices(4) = this%vertices(1)%ptr%bot_parent
                vertex_indices(5) = this%vertices(2)%ptr%bot_parent
                vertex_indices(6) = this%vertices(3)%ptr%bot_parent

                ! Set default influence
                allocate(v(6,3), source=0.)

            else

                ! Body panels are influenced by only one set of vertices
                allocate(vertex_indices, source=this%vertex_indices)

                ! Set default influence
                allocate(v(3,3), source=0.)

            end if
        end if

        ! Check DoD
        if (dod_info%in_dod) then

            ! Calculate geometric parameters
            geom = this%get_field_point_geometry(eval_point, freestream, dod_info)

            if (influence_calc_type == 'johnson') then

                ! Get integrals
                call this%calc_integrals(geom, "velocity", "doublet", dod_info, freestream, H, F)

                if (doublet_order .eq. 1) then

                    ! Specify influencing vertices
                    allocate(vertex_indices(3), source=this%vertex_indices)

                    ! Calculate velocity
                    allocate(v(3,3), source=0.)

                end if

                ! Clean up
                deallocate(H)
                deallocate(F)

            end if
        end if

    end function panel_calc_doublet_velocity


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
        real :: s

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

        ! Calculate doublet parameters (derivatives)
        mu_params = matmul(this%S_mu_inv, mu_verts)

        ! Calculate tangential velocity jump in panel coordinates E&M Eq. (N.1.11b)
        dv(1) = mu_params(2)
        dv(2) = mu_params(3)
        dv(3) = 0.

        ! Transform to global coordinates
        dv = matmul(transpose(this%A_g_to_ls), dv)

        ! Get source strength
        if (mirrored) then
            s = sigma(this%index+size(sigma)/2)
        else
            s = sigma(this%index)
        end if

        ! Add normal velocity jump in global coords E&M Eq. (N.1.11b)
        dv = dv + s*this%n_g/inner(this%nu_g, this%n_g)

        ! Mirror if necessary
        if (mirrored) then
            dv = mirror_about_plane(dv, mirror_plane)
        end if

    end function panel_get_velocity_jump

    
end module panel_mod