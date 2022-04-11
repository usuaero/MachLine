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


    type eval_point_geom
        ! Container type for the geometric parameters necessary for calculating a panel's influence on a given field point

        real,dimension(3) :: P_g ! Point position in global coords
        real,dimension(2) :: P_ls ! Transformed point in panel plane
        real :: h, h2 ! Transformed height above panel
        real,dimension(3) :: a, g2, l1, l2, R1, R2 ! Edge integration parameters (Johnson's method)
        real,dimension(3) :: xm, ym1, ym2, sm1, sm2, s1, s2 ! Edge integration parameters (Ehler's method)

        contains

            procedure :: init => eval_point_geom_init

    end type eval_point_geom


    type integrals
        ! Container type for the fundamental integrals used to calculate influence coefficients

        real :: H111, H211, H121 ! Source integrals
        real :: hH113, H213, H123, H313, H223, H133 ! Doublet integrals; we use hH(1,1,3) because it can be reliably calculated, unlike H(1,1,3)
        real,dimension(:),allocatable :: F111, F211, F121 ! Necessary line integrals
        real :: Q1 ! Davis-Ehlers integral
        real,dimension(:),allocatable :: w0 ! Davis-Ehlers integral

    end type integrals


    type dod
        ! Container type for parameters of whether a panel lies in a point's domain of dependence

        logical :: in_dod = .true.
        logical,dimension(3) :: verts_in_dod = .true.
        logical,dimension(3) :: edges_in_dod = .true.

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
        real,dimension(3,3) :: A_g_to_ls, A_ls_to_g ! Coordinate transformation matrices
        real,dimension(:,:),allocatable :: vertices_ls, midpoints_ls ! Location of the vertices and edge midpoints described in local scaled coords
        real,dimension(:,:),allocatable :: t_hat_g, t_hat_ls ! Edge unit tangents
        real,dimension(:,:),allocatable :: n_hat_g, n_hat_ls ! Edge unit outward normals
        real,dimension(:),allocatable :: b, sqrt_b ! Edge parameter
        real :: A ! Surface area
        real,dimension(:,:),allocatable :: S_mu_inv, S_sigma_inv ! Matrix relating doublet/source strengths to doublet/source influence parameters
        integer,dimension(:),allocatable :: vertex_indices ! Indices of this panel's vertices in the mesh vertex array
        logical :: in_wake ! Whether this panel belongs to a wake mesh
        integer,dimension(3) :: abutting_panels ! Indices of panels abutting this one
        integer,dimension(3) :: edges ! Indices of this panel's edges
        integer :: top_parent, bot_parent ! Indices of the top and bottom panels this panel's strength is determined by (for a wake panel w/ constant doublet strength)
        integer :: r ! Panel inclination indicator; r=-1 -> superinclined, r=1 -> subinclined
        real,dimension(3) :: m, l ! Edge slope and its inverse
        integer,dimension(3) :: q ! Edge type indicator; q=1 -> subsonic, q=-1 -> supersonic
        real :: J ! Local scaled transformation Jacobian

        contains

            procedure :: init => panel_init_3
            procedure :: calc_derived_properties =>panel_calc_derived_properties
            procedure :: calc_area => panel_calc_area
            procedure :: calc_normal => panel_calc_normal
            procedure :: calc_centroid => panel_calc_centroid
            procedure :: calc_transforms => panel_calc_transforms
            procedure :: calc_g_to_ls_transform => panel_calc_g_to_ls_transform
            procedure :: calc_edge_vectors => panel_calc_edge_vectors
            procedure :: calc_singularity_matrices => panel_calc_singularity_matrices
            procedure :: get_vertex_loc => panel_get_vertex_loc
            procedure :: get_vertex_index => panel_get_vertex_index
            procedure :: touches_vertex => panel_touches_vertex
            procedure :: point_to_vertex_clone => panel_point_to_vertex_clone
            procedure :: check_dod => panel_check_dod
            procedure :: calc_subsonic_geom => panel_calc_subsonic_geom
            procedure :: calc_supersonic_subinc_geom => panel_calc_supersonic_subinc_geom
            procedure :: E_i_M_N_K => panel_E_i_M_N_K
            procedure :: calc_subsonic_edge_integrals => panel_calc_subsonic_edge_integrals
            procedure :: calc_subsonic_panel_integrals => panel_calc_subsonic_panel_integrals
            procedure :: calc_supersonic_subinc_edge_integrals => panel_calc_supersonic_subinc_edge_integrals
            procedure :: calc_supersonic_subinc_panel_integrals => panel_calc_supersonic_subinc_panel_integrals
            procedure :: calc_integrals => panel_calc_integrals
            procedure :: calc_potentials => panel_calc_potentials
            procedure :: calc_velocities => panel_calc_velocities
            procedure :: get_velocity_jump => panel_get_velocity_jump

    end type panel

    
contains


    subroutine eval_point_geom_init(this, P, A_g_to_ls, centroid)

        implicit none

        class(eval_point_geom),intent(out) :: this
        real,dimension(3),intent(in) :: P, centroid
        real,dimension(3,3),intent(in) :: A_g_to_ls

        real,dimension(3) :: x

        ! Store point
        this%P_g = P

        ! Transform to local scaled coordinates
        x = matmul(A_g_to_ls, this%P_g-centroid)
        this%P_ls = x(1:2)
        this%h = x(3) ! Equivalent to E&M Eq. (J.7.41)
        this%h2 = this%h**2

        ! These are sometimes accessed when the DoD is not checked, so they need to be set to zero
        this%R1 = 0.
        this%R2 = 0.
        this%a = 0.
    
    end subroutine eval_point_geom_init


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
        allocate(this%midpoints(3,this%N))
        do i=1,this%N-1
            this%midpoints(:,i) = 0.5*(this%get_vertex_loc(i)+this%get_vertex_loc(i+1))
        end do
        this%midpoints(:,this%N) = 0.5*(this%get_vertex_loc(1)+this%get_vertex_loc(this%N))

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
            this%A = 0.5*norm2(cross(d1, d2))

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
        d1 = this%midpoints(:,2)-this%midpoints(:,1)
        d2 = this%midpoints(:,3)-this%midpoints(:,2)

        ! Find normal
        this%n_g = cross(d1, d2)
        this%n_g = this%n_g/norm2(this%n_g)

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
        call this%calc_g_to_ls_transform(freestream)

        ! Calculate properties dependent on the transforms
        call this%calc_edge_vectors(freestream)
        call this%calc_singularity_matrices()

    end subroutine panel_calc_transforms


    subroutine panel_calc_g_to_ls_transform(this, freestream)
        ! Calculates the necessary transformations to move from global to local, scaled coordinates (Eq. (E.0.1) in Epton and Magnus)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(3) :: u0, v0
        real,dimension(3,3) :: B_mat_ls
        real :: x, y
        integer :: i, rs

        ! Get in-panel basis vectors
        if (abs(inner(this%n_g, freestream%c_hat_g) - 1.) < 1e-12) then ! Check the freestream isn't aligned with the normal vector
            v0 = this%get_vertex_loc(2)-this%get_vertex_loc(1)
        else
            v0 = cross(this%n_g, freestream%c_hat_g)
        end if
        v0 = v0/norm2(v0)
        u0 = cross(v0, this%n_g)
        u0 = u0/norm2(u0)

        ! Calculate compressible parameters
        this%nu_g = matmul(freestream%B_mat_g, this%n_g)
        x = inner(this%n_g, this%nu_g)

        ! Check for Mach-inclined panels
        if (freestream%supersonic .and. abs(x) < 1e-12) then
            write(*,*) "!!! Panel", this%index, "is Mach-inclined, which is not allowed. Quitting..."
            stop
        end if

        ! Calculate panel inclination indicator (E&M Eq. (E.3.16b))
        this%r = sign(1., x) ! r=-1 -> superinclined, r=1 -> subinclined

        ! Check for superinclined panels
        if (this%r < 0) then
            write(*,*) "!!! Panel", this%index, "is superinclined, which is not allowed. Quitting..."
            stop
        end if

        ! Other inclination parameters
        rs = this%r*freestream%s

        ! Calculate transformation
        y = 1./sqrt(abs(x))
        this%A_g_to_ls(1,:) = y*matmul(freestream%C_mat_g, u0)
        this%A_g_to_ls(2,:) = rs/freestream%B*matmul(freestream%C_mat_g, v0)
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
        this%J = 1./(freestream%B*sqrt(abs(1.-freestream%M_inf**2*inner(freestream%c_hat_g, this%n_g)**2)))

        ! Transform vertex and midpoint coords to ls
        allocate(this%vertices_ls(2,this%N))
        allocate(this%midpoints_ls(2,this%N))
        do i=1,this%N

            ! Vertices
            this%vertices_ls(:,i) = matmul(this%A_g_to_ls(1:2,:), this%get_vertex_loc(i)-this%centroid)

            ! Midpoints
            this%midpoints_ls(:,i) = matmul(this%A_g_to_ls(1:2,:), this%midpoints(:,i)-this%centroid)

        end do

        ! Calculate local scaled metric matrices
        B_mat_ls = 0.
        B_mat_ls(1,1) = freestream%B**2*rs
        B_mat_ls(2,2) = freestream%B**2
        B_mat_ls(3,3) = freestream%B**2*this%r

        ! Check calculation (E&M Eq. (E.2.19))
        if (any(abs(B_mat_ls - matmul(this%A_g_to_ls, matmul(freestream%B_mat_g, transpose(this%A_g_to_ls)))) > 1e-12)) then
            write(*,*) "!!! Calculation of local scaled coordinate transform failed. Quitting..."
            stop
        end if
    
    end subroutine panel_calc_g_to_ls_transform


    subroutine panel_calc_edge_vectors(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(3) :: d_g
        real,dimension(2) :: d_ls
        integer :: i, i_next

        ! Allocate memory
        allocate(this%t_hat_g(3,this%N))
        allocate(this%n_hat_g(3,this%N))
        allocate(this%t_hat_ls(2,this%N))
        allocate(this%n_hat_ls(2,this%N))
        allocate(this%b(this%N))
        allocate(this%sqrt_b(this%N))

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N)+1

            ! Calculate edge vector based on index
            d_g = this%get_vertex_loc(i_next)-this%get_vertex_loc(i)

            ! Calculate tangent in global coords
            this%t_hat_g(:,i) = d_g/norm2(d_g)

            ! Calculate tangent in local scaled coords 
            ! This purposefully does not match E&M Eq. (J.6.43);
            ! This is the formula used in the PAN AIR source code (for subinclined panels)
            d_ls = this%vertices_ls(:,i_next) - this%vertices_ls(:,i)
            this%t_hat_ls(:,i) = d_ls/norm2(d_ls)
            !this%t_hat_ls(:,i) = e/sqrt(abs(rs*e(1)**2 + e(2)**2)) ! E&M Eq. (J.6.43)

            ! Calculate edge outward normal
            this%n_hat_g(:,i) = cross(this%t_hat_g(:,i), this%n_g)

            ! Calculate edge type indicator (E&M Eq. (J.6.48)
            this%q(i) = sign(1., this%r*this%t_hat_ls(1,i)**2 + freestream%s*this%t_hat_ls(2,i)**2)

        end do

        ! Calculate edge normal in local scaled coords E&M Eq. (J.6.45)
        this%n_hat_ls(1,:) = this%t_hat_ls(2,:)
        this%n_hat_ls(2,:) = -this%t_hat_ls(1,:)

        ! Calculate edge parameter (Ehlers Eq. (E14))
        this%b = (this%n_hat_ls(1,:) - this%n_hat_ls(2,:))*(this%n_hat_ls(1,:) + this%n_hat_ls(2,:))
        this%sqrt_b = sqrt(abs(this%b))

        ! Calculate edge slope
        this%m = this%t_hat_ls(2,:)/this%t_hat_ls(1,:)
        this%l = 1./this%m
    
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
                S_mu(:,2) = this%vertices_ls(1,:)
                S_mu(:,3) = this%vertices_ls(2,:)

                ! Invert
                call matinv(3, S_mu, this%S_mu_inv)

            else if (doublet_order .eq. 2) then

                ! Allocate influence matrix
                allocate(S_mu(6,6))
                allocate(this%S_mu_inv(6,6))

                ! Set values
                S_mu(:,1) = 1.

                S_mu(1:3,2) = this%vertices_ls(1,:)
                S_mu(1:3,3) = this%vertices_ls(2,:)
                S_mu(1:3,4) = this%vertices_ls(1,:)**2
                S_mu(1:3,5) = this%vertices_ls(1,:)*this%vertices_ls(2,:)
                S_mu(1:3,6) = this%vertices_ls(2,:)**2
                
                S_mu(4:6,2) = this%midpoints_ls(1,:)
                S_mu(4:6,3) = this%midpoints_ls(2,:)
                S_mu(4:6,4) = this%midpoints_ls(1,:)**2
                S_mu(4:6,5) = this%midpoints_ls(1,:)*this%midpoints_ls(2,:)
                S_mu(4:6,6) = this%midpoints_ls(2,:)**2

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
            S_sigma(:,2) = this%vertices_ls(1,:)
            S_sigma(:,3) = this%vertices_ls(2,:)

            ! Invert
            call matinv(3, S_sigma, this%S_sigma_inv)

            deallocate(S_sigma)

        else if (source_order .eq. 2) then

            ! Allocate influence matrix
            allocate(S_sigma(6,6))
            allocate(this%S_sigma_inv(6,6))

            ! Set values
            S_sigma(:,1) = 1.

            S_sigma(1:3,2) = this%vertices_ls(1,:)
            S_sigma(1:3,3) = this%vertices_ls(2,:)
            S_sigma(1:3,4) = this%vertices_ls(1,:)**2
            S_sigma(1:3,5) = this%vertices_ls(1,:)*this%vertices_ls(2,:)
            S_sigma(1:3,6) = this%vertices_ls(2,:)**2
            
            S_sigma(4:6,2) = this%midpoints_ls(1,:)
            S_sigma(4:6,3) = this%midpoints_ls(2,:)
            S_sigma(4:6,4) = this%midpoints_ls(1,:)**2
            S_sigma(4:6,5) = this%midpoints_ls(1,:)*this%midpoints_ls(2,:)
            S_sigma(4:6,6) = this%midpoints_ls(2,:)**2

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
            ! This shouldn't be necessary, but I'll keep it here for now.
            dod_info%in_dod = .true.
            dod_info%verts_in_dod = .true.
            dod_info%edges_in_dod = .true.

        end if
    
    end function panel_check_dod


    function panel_calc_subsonic_geom(this, eval_point, freestream) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point in subsonic flow.

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(eval_point_geom) :: geom

        real,dimension(2) :: d_ls
        real,dimension(3) :: d_g
        integer :: i, i_next
        real :: val

        ! Initialize
        call geom%init(eval_point, this%A_g_to_ls, this%centroid)

        ! Calculate edge quantities
        do i=1,this%N

            i_next = mod(i, this%N) + 1

            ! Calculate displacements
            d_ls = this%vertices_ls(:,i) - geom%P_ls

            ! Perpendicular distance in plane from evaluation point to edge E&M Eq. (J.6.46) and (J.7.53)
            geom%a(i) = inner2(d_ls, this%n_hat_ls(:,i)) ! Definition

            ! Integration length on edge to start vertex (E&M Eq. (J.6.47))
            geom%l1(i) = inner2(d_ls, this%t_hat_ls(:,i))

            ! Distance from evaluation point to start vertex E&M Eq. (J.8.8)
            geom%R1(i) = sqrt(d_ls(1)**2 + d_ls(2)**2 + geom%h2)

            ! Calculate square of the perpendicular distance to edge
            geom%g2(i) = geom%a(i)**2 + geom%h2 ! Subsonic version of Ehlers Eq. (E14) and what is used in PAN AIR

            ! Displacement from end vertex
            d_ls = this%vertices_ls(:,i_next) - geom%P_ls

            ! Integration length on edge to end vertex
            geom%l2(i) = inner2(d_ls, this%t_hat_ls(:,i)) ! Definition; same as used in PAN AIR

        end do

        ! Distance from evaluation point to end vertices
        geom%R2 = cshift(geom%R1, 1)

    end function panel_calc_subsonic_geom


    function panel_calc_supersonic_subinc_geom(this, eval_point, freestream, dod_info) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        type(eval_point_geom) :: geom

        real,dimension(2) :: d_ls, d
        real :: x
        integer :: i, i_next

        ! Initialize
        call geom%init(eval_point, this%A_g_to_ls, this%centroid)

        ! Calculate edge quantities
        do i=1,this%N
            if (dod_info%edges_in_dod(i)) then

                i_next = mod(i, this%N) + 1

                ! Calculate displacements
                d_ls = this%vertices_ls(:,i) - geom%P_ls

                ! Perpendicular distance in plane from evaluation point to edge E&M Eq. (J.6.46) and (J.7.53)
                geom%a(i) = inner2(d_ls, this%n_hat_ls(:,i)) ! Definition

                ! Calculate square of the perpendicular distance to edge
                geom%g2(i) = geom%a(i)**2 - this%b(i)*geom%h2 ! Ehlers Eq. (E14) and what is used in supsbi in PAN AIR.

                if (dod_info%verts_in_dod(i)) then
        
                    ! Integration length on edge to start vertex Ehlers Eq. (E14)
                    geom%l1(i) = -d_ls(1)*this%t_hat_ls(1,i) + d_ls(2)*this%t_hat_ls(2,i)

                    ! Distance from evaluation point to start vertex Ehlers Eq. (E15)
                    geom%R1(i) = sqrt(d_ls(1)**2 - d_ls(2)**2 - geom%h2)

                else

                    ! Comes from PAN AIR (subroutine supsbi) to enforce Ehlers Eq. (E15)
                    geom%l1(i) = sqrt(abs(geom%g2(i)))

                end if

                if (dod_info%verts_in_dod(i_next)) then

                    ! Displacement from end vertex
                    d_ls = this%vertices_ls(:,i_next) - geom%P_ls

                    ! Integration length on edge to end vertex
                    geom%l2(i) = -d_ls(1)*this%t_hat_ls(1,i) + d_ls(2)*this%t_hat_ls(2,i) ! Definition; same as used in supsbi in PAN AIR

                    ! Distance from evaluation point to end vertex
                    geom%R2(i) = sqrt(d_ls(1)**2 - d_ls(2)**2 - geom%h2) ! Ehlers Eq. (E15)

                else

                    ! Comes from PAN AIR (subroutine supsbi) to enforce Ehlers Eq. (E15)
                    geom%l2(i) = -sqrt(abs(geom%g2(i)))

                end if

                ! Get vector describing edge
                d = this%vertices_ls(:,i_next) - this%vertices_ls(:,i)

                ! Check integration direction
                x = d(1)*this%m(i) - d(2)

                ! Calculate s and sm
                ! First vertex
                d_ls = geom%P_ls - this%vertices_ls(:,i)
                geom%sm1(i) = d_ls(1) ! May be derived from Ehlers Eq. (E25) or (5.13)
                geom%s1(i) = d_ls(2)

                ! Second vertex
                d_ls = geom%P_ls - this%vertices_ls(:,i_next)
                geom%sm2(i) = d_ls(1) ! May be derived from Ehlers Eq. (E25) or (5.13)
                geom%s2(i) = d_ls(2)

                ! Subsonic or sonic edge (these are the hatted quantities in Ehlers and Davis)
                if (abs(this%m(i)) <= 1.) then

                    ! Calculate xhatm
                    geom%xm(i) = -geom%sm1(i)*this%m(i) + geom%s1(i) ! Opposite of Ehlers Eq. (A25), but lines up with Ehlers Eq. (E35)

                    ! Calculate yhatm for first vertex
                    geom%ym1(i) = geom%sm1(i) - geom%s1(i)*this%m(i) ! Ehlers p. 109

                    ! Calculate yhatm for second vertex
                    geom%ym2(i) = geom%sm2(i) - geom%s2(i)*this%m(i) ! Ehlers p. 109

                ! Supersonic edge
                else

                    ! Calculate xm
                    geom%xm(i) = geom%sm1(i) - geom%s1(i)*this%l(i) ! Ehlers Eq. (5.13) or (A2)

                    ! Calculate ym for first vertex
                    geom%ym1(i) = geom%s1(i) - this%l(i)*geom%sm1(i) ! Ehlers p. 104; the definition of this given on p. 108 is opposite

                    ! Calculate ym for second vertex
                    geom%ym2(i) = geom%s2(i) - this%l(i)*geom%sm2(i) ! Ehlers p. 104

                end if

            end if
        end do

    end function panel_calc_supersonic_subinc_geom


    function panel_E_i_M_N_K(this, geom, i, M, N, K, freestream) result(E)
        ! Calculates E_i(M,N,K)

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: i, M, N, K
        type(flow),intent(in) :: freestream

        real :: E

        real :: E1, E2
        integer :: i_next

        i_next = mod(i, this%N) + 1

        ! Evaluate at start vertex
        if (geom%R1(i) /= 0.) then
            E1 = ((this%vertices_ls(1,i)-geom%P_ls(1))**(M-1)*(this%vertices_ls(2,i)-geom%P_ls(2))**(N-1))/geom%R1(i)**K
        else
            E1 = 0.
        end if

        ! Evaluate at end vertex
        if (geom%R1(i_next) /= 0.) then
            E2 = ((this%vertices_ls(1,i_next)-geom%P_ls(1))**(M-1)*(this%vertices_ls(2,i_next)-geom%P_ls(2))**(N-1))/geom%R2(i)**K
        else
            E2 = 0.
        end if

        ! Calculate difference
        E = E2 - E1

    end function panel_E_i_M_N_K


    subroutine panel_calc_subsonic_edge_integrals(this, geom, freestream, int)
        ! Calculates the F integrals necessary for determining the influence of a triangular panel in subsonic flow.
        ! This is a pared-down version of the algorithm presented by Johnson (1980) Appendix D.3.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(flow),intent(in) :: freestream
        type(integrals),intent(inout) :: int

        integer :: i
        
        ! Allocate storage
        allocate(int%F111(this%N))

        ! Loop through edges
        do i=1,this%N

            ! Calculate F(1,1,1)
            ! Within edge (Johnson Eq. (D.60))
            if (sign(1., geom%l1(i)) /= sign(1., geom%l2(i))) then

                ! Check for point on perimeter
                if(sqrt(geom%g2(i)) < 1e-12) then
                    write(*,*) "Detected control point on perimeter of panel. Quitting..."
                    stop
                end if

                ! Calculate
                int%F111(i) = log( ( (geom%R1(i) - geom%l1(i)) * (geom%R2(i) + geom%l2(i)) ) / geom%g2(i) )

            ! Above or below edge; this is a unified form of Johnson Eq. (D.60)
            else

                ! Check for point on perimeter
                if (min(geom%R1(i), geom%R2(i)) < 1e-12) then
                    write(*,*) "Detected control point on perimeter of panel. Quitting..."
                    stop
                end if

                ! Calculate
                int%F111(i) = sign(1., geom%l1(i)) * log( (geom%R2(i) + abs(geom%l2(i))) / (geom%R1(i) + abs(geom%l1(i))) )

            end if

        end do

    end subroutine panel_calc_subsonic_edge_integrals


    subroutine panel_calc_supersonic_subinc_edge_integrals(this, geom, dod_info, freestream, int)
        ! Calculates the F integrals necessary to determine the influence of a subinclined triangular panel in supersonic flow.
        ! Taken from Ehlers et al. (1979) Appendix E.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        type(integrals),intent(inout) :: int

        real :: v_xi, v_eta, F1, F2, eps, eps2, series, x, x2, x_inv, zr, val
        integer :: i
        
        ! Allocate integral storage
        allocate(int%F111(this%N), source=0.)
        allocate(int%w0(this%N), source=0.)

        ! Loop through edges
        do i=1,this%N

            if (dod_info%edges_in_dod(i)) then

                ! Calculate F(1,1,1) (Ehlers Eq. (E22))

                ! Neither endpoint in DoD (Ehlers Eq. (E22))
                if (geom%R1(i) == 0. .and. geom%R2(i) == 0.) then
                    int%F111(i) = pi/this%sqrt_b(i)

                else

                    ! Calculate preliminary quantities (Ehlers Eqs. (E19-20))
                    if (this%b(i) >= 0) then
                        F1 = (geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)) / geom%g2(i)
                        F2 = (this%b(i)*geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)) / geom%g2(i)

                    else
                        F1 = (geom%R2(i) + geom%R1(i))*(geom%R2(i) - geom%R1(i)) / (geom%l1(i)*geom%R2(i) + geom%l2(i)*geom%R1(i))
                        F2 = (geom%g2(i) - geom%l1(i)**2 - geom%l2(i)**2) / &
                             (this%b(i)*geom%R1(i)*geom%R2(i) - geom%l1(i)*geom%l2(i))

                    end if

                    ! Check for poorly-conditioned edge
                    if (abs(F2) > 100.*this%sqrt_b(i)*abs(F1)) then

                        ! This mimics what is done in PAN AIR, for efficiency
                        !int%F111(i) = -eps*(1 - this%b(i)*eps**2/3. + this%b(i)**2*eps**4/5. - this%b(i)**3*eps**6/7.)
                        eps = F1/F2
                        eps2 = eps*eps
                        series = eps*eps2*( 1./3. - 0.2*this%b(i)*eps2 + this%b(i)**2*eps2**2/7. )
                        int%F111(i) = -eps + this%b(i)*series

                    ! Supersonic edge
                    else if (this%b(i) > 0) then

                        int%F111(i) = -atan2(this%sqrt_b(i)*F1, F2) / this%sqrt_b(i)

                    ! Subsonic edge
                    else

                        ! Calculate preliminary quantities
                        F1 = this%sqrt_b(i)*geom%R1(i) + abs(geom%l1(i))
                        F2 = this%sqrt_b(i)*geom%R2(i) + abs(geom%l2(i))

                        ! Calculate F(1,1,1)
                        int%F111(i) = -sign(1., this%n_hat_ls(2,i)) / this%sqrt_b(i) * log(F1/F2)

                    end if
                end if

                ! Calculate w0
                
                ! Subsonic edge
                if (abs(this%m(i)) < 1.) then

                    x2 = 1.-this%m(i)**2
                    x = sqrt(x2)
                    x_inv = 1./x

                    ! Both endpoints in DoD (Davis Eq. (A.17) has a typo)
                    if (geom%R1(i) /= 0. .and. geom%R2(i) /= 0.) then

                        F1 = geom%ym2(i) + geom%R2(i)*x
                        F2 = geom%ym1(i) + geom%R1(i)*x
                        int%w0(i) = x_inv*log(F1/F2)

                    ! First endpoint in DoD
                    else if (geom%R1(i) /= 0.) then

                        F1 = geom%ym1(i) + geom%R1(i)*x
                        F2 = geom%ym1(i) - geom%R1(i)*x
                        int%w0(i) = -0.5*x_inv*log(F1/F2)

                    ! Second endpoint in DoD
                    else if (geom%R2(i) /= 0.) then

                        F1 = geom%ym2(i) + geom%R2(i)*x
                        F2 = geom%ym2(i) - geom%R2(i)*x
                        int%w0(i) = 0.5*x_inv*log(F1/F2)

                    end if

                ! Supersonic edge
                else if (abs(this%m(i)) > 1.) then

                    x2 = 1.-this%l(i)**2
                    x = sqrt(x2)
                    x_inv = 1./x

                    ! Both endpoints in DoD (Ehlers p. 108)
                    if (geom%R1(i) /= 0. .and. geom%R2(i) /= 0.) then

                        F1 = x*(geom%ym1(i)*geom%R2(i) - geom%ym2(i)*geom%R1(i))
                        F2 = geom%ym1(i)*geom%ym2(i) + x2*geom%R1(i)*geom%R2(i)
                        int%w0(i) = x_inv*atan2(F1, F2)

                    ! First endpoint in DoD
                    else if (geom%R1(i) /= 0.) then

                        int%w0(i) = ( sign(pi2, -geom%ym2(i)) - atan2(-geom%ym1(i), geom%R1(i)*x) ) * x_inv

                    ! Second endpoint in DoD
                    else if (geom%R2(i) /= 0.) then

                        int%w0(i) = ( atan2(-geom%ym2(i), geom%R2(i)*x) - sign(pi2, -geom%ym1(i)) ) * x_inv

                    ! Neither endpoint in DoD (Ehlers p. 108)
                    else
                        int%w0(i) = ( sign(pi2, -geom%ym2(i)) - sign(pi2, -geom%ym1(i)) ) * x_inv

                    end if

                ! Sonic edge
                else

                    ! Combination of w0 equation in Ehlers A5 and Davis Eqs. (A29-30) and equation for w0 on Ehlers p. 109
                    int%w0(i) = geom%R2(i)/geom%ym2(i) - geom%R1(i)/geom%ym1(i)

                end if

                ! Check computation for nearly-sonic edges
                if (abs(this%m(i)-1.) < 1e-8) then
                    x = 1.-this%l(i)**2
                    zr = (geom%ym1(i)*geom%R2(i) - geom%ym2(i)*geom%R1(i))
                    zr = zr/(geom%ym1(i)*geom%ym2(i) + x*geom%R1(i)*geom%R2(i))
                    val = zr*(1-x*zr**2/3. + x**2*zr**4/5. - x**3*zr**6/7.)
                    write(*,*)
                    write(*,*) abs(this%m(i))
                    write(*,*) int%w0(i)
                    write(*,*) val
                    write(*,*) abs(int%w0(i)) - val
                end if

            end if
        end do

    end subroutine panel_calc_supersonic_subinc_edge_integrals


    subroutine panel_calc_subsonic_panel_integrals(this, geom, freestream, int)
        ! Calculates the necessary H integrals to determine the influence of a panel in subsonic flow.
        ! Taken from Johnson (1980) Appendix D.3.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(flow),intent(in) :: freestream
        type(integrals),intent(inout) :: int

        real :: S, C, nu, c1, c2, x
        integer :: i, m, n, k
        real,dimension(:),allocatable :: v_xi, v_eta

        ! Get edge normal derivatives
        allocate(v_xi(this%N), source=this%n_hat_ls(1,:))
        allocate(v_eta(this%N), source=this%n_hat_ls(2,:))

        ! Calculate hH(1,1,3)
        ! Not close to panel plane
        if (abs(geom%h) > 1e-12) then ! The nonzero h check seems to be more reliable than that proposed by Johnson

            ! Calculate and hH(1,1,3) (Johnson Eqs. (D.41) and (G.24))
            int%hH113 = 0.
            do i=1,this%N

                ! Calculate intermediate quantities
                c1 = geom%g2(i)+abs(geom%h)*geom%R1(i)
                c2 = geom%g2(i)+abs(geom%h)*geom%R2(i)
        
                ! Add surface integral
                S = geom%a(i)*(geom%l2(i)*c1 - geom%l1(i)*c2)
                C = c1*c2 + geom%a(i)**2*geom%l1(i)*geom%l2(i)
                x = atan2(S, C)
                int%hH113 = int%hH113 + x
        
            end do
        
            ! Calculate hH(1,1,3) (Johnson Eq. (D.42)
            int%hH113 = sign(int%hH113, geom%h)

        else

            ! Close to panel plane but outside Sigma
            if (all(geom%a < 0.)) then
                int%hH113 = 0.

            ! Close to panel plane but inside Sigma
            else
                int%hH113 = sign(2.*pi, geom%h)

            end if
        end if

        ! Calculate H(1,1,1)
        int%H111 = -geom%h*int%hH113 + sum(geom%a*int%F111)

        ! Calculate H(2,1,3) and H(1,2,3)
        int%H213 = -sum(v_xi*int%F111)
        int%H123 = -sum(v_eta*int%F111)

        ! Clean up
        deallocate(v_xi)
        deallocate(v_eta)

    end subroutine panel_calc_subsonic_panel_integrals


    subroutine panel_calc_supersonic_subinc_panel_integrals(this, geom, dod_info, freestream, int)
        ! Calculates the necessary H integrals to determine the influence of a subinclined panel in supersonic flow.
        ! Taken from Ehlers et al. (1979) Appendix E.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        type(integrals),intent(inout) :: int

        real :: F1, F2
        integer :: i, i_next
        real,dimension(:),allocatable :: v_xi, v_eta

        ! Get edge normal derivatives
        allocate(v_xi(this%N), source=this%n_hat_ls(1,:))
        allocate(v_eta(this%N), source=this%n_hat_ls(2,:))

        ! Calculate hH(1,1,3) (Ehlers Eq. (E18))
        int%hH113 = 0.
        int%Q1 = 0.
        if (abs(geom%h) > 1e-12) then ! Check the point is off the panel plane

            ! Add influence of each edge
            do i=1,this%N

                if (dod_info%edges_in_dod(i)) then

                    i_next = mod(i, this%N) + 1

                    ! Check for at least one endpoint in the DoD
                    if (dod_info%verts_in_dod(i) .or. dod_info%verts_in_dod(i_next)) then

                        ! Calculate intermediate quantities Ehlers Eq. (E19) and (E20)
                        if (this%b(i) >= 0 ) then
                            F1 = (geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)) / geom%g2(i)
                            F2 = (this%b(i)*geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)) / geom%g2(i)

                        else
                            F1 = (geom%R2(i) + geom%R1(i))*(geom%R2(i) - geom%R1(i)) / &
                                 (geom%l1(i)*geom%R2(i) + geom%l2(i)*geom%R1(i))
                            F2 = (geom%g2(i) - geom%l1(i)**2 - geom%l2(i)**2) / &
                                 (this%b(i)*geom%R1(i)*geom%R2(i) - geom%l1(i)*geom%l2(i))
                        end if

                        ! Add to surface integral
                        int%hH113 = int%hH113 + atan2(geom%h*geom%a(i)*F1, geom%R1(i)*geom%R2(i) + geom%h2*F2)

                    ! Neither endpoint is in
                    else
                        int%hH113 = int%hH113 + sign(pi, geom%h*v_xi(i))

                    end if

                    ! Calculate Q1

                    ! Subsonic or sonic edge
                    if (abs(this%m(i)) <= 1.) then

                        ! Both endpoints in
                        if (geom%R1(i) /= 0. .and. geom%R2(i) /= 0.) then

                            F1 = geom%h*geom%xm(i)*(geom%ym1(i)*geom%R2(i) - geom%ym2(i)*geom%R1(i))
                            F2 = geom%xm(i)**2*geom%R1(i)*geom%R2(i) + geom%h2*geom%ym1(i)*geom%ym2(i)
                            int%Q1 = int%Q1 + atan2(F1, F2)

                        ! First endpoint in
                        else if (geom%R1(i) /= 0.) then

                            int%Q1 = int%Q1 - sign(atan2(geom%xm(i)*geom%R1(i), abs(geom%h)*geom%ym1(i)), geom%h)

                        ! Second endpoint in
                        else if (geom%R2(i) /= 0.) then

                            int%Q1 = int%Q1 + sign(atan2(geom%xm(i)*geom%R2(i), abs(geom%h)*geom%ym2(i)), geom%h)

                        end if

                    ! Supersonic edge (Ehlers p. 108)
                    else if (abs(this%m(i)) > 1.) then

                        ! Both endpoints in
                        if (geom%R1(i) /= 0. .and. geom%R2(i) /= 0.) then

                            F1 = geom%h*geom%xm(i)*(geom%ym1(i)*geom%R2(i) - geom%ym2(i)*geom%R1(i))
                            F2 = geom%xm(i)**2*geom%R1(i)*geom%R2(i) + geom%h2*geom%ym1(i)*geom%ym2(i)
                            int%Q1 = int%Q1 + atan2(F1, F2)

                        ! First endpoint in (Davis Eq. (A.28))
                        else if (geom%R1(i) /= 0.) then

                            int%Q1 = int%Q1 + sign(pi2, -geom%h*geom%ym2(i)) - atan2(-geom%h*geom%ym1(i), geom%xm(i)*geom%R1(i))

                        ! Second endpoint in (Davis Eq. (A.28))
                        else if (geom%R2(i) /= 0.) then

                            int%Q1 = int%Q1 + atan2(-geom%h*geom%ym2(i), geom%xm(i)*geom%R2(i)) - sign(pi2, -geom%h*geom%ym1(i))

                        ! Neither endpoint in (Davis Eq. (A.28))
                        else

                            int%Q1 = int%Q1 + sign(pi2, -geom%h*geom%ym2(i)) - sign(pi2, -geom%h*geom%ym1(i))

                        end if

                    end if

                end if
        
            end do
        end if

        ! Calculate H(1,1,1) (Ehlers Eq. (E9) has a typo)
        ! What exactly that typo is, I don't know
        int%H111 = -geom%h*int%hH113 - sum(geom%a*int%F111)

        ! Calculate H(2,1,3) (Ehlers Eq. (E5) has a typo)
        int%H213 = sum(v_xi*int%F111)

        ! Calculate H(1,2,3) (Ehlers Eq. (E6) has a typo)
        int%H123 = -sum(v_eta*int%F111)

        ! Clean up
        deallocate(v_xi)
        deallocate(v_eta)

    end subroutine panel_calc_supersonic_subinc_panel_integrals


    function panel_calc_integrals(this, geom, influence_type, singularity_type, dod_info, freestream) result(int)
        ! Calculates the H and F integrals necessary for the given influence

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        character(len=*),intent(in) :: influence_type, singularity_type
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream

        type(integrals) :: int

        ! Calculate necessary integrals based on the flow condition and panel type
        if (freestream%supersonic) then
            call this%calc_supersonic_subinc_edge_integrals(geom, dod_info, freestream, int)
            call this%calc_supersonic_subinc_panel_integrals(geom, dod_info, freestream, int)
        else
            call this%calc_subsonic_edge_integrals(geom, freestream, int)
            call this%calc_subsonic_panel_integrals(geom, freestream, int)
        end if

    end function panel_calc_integrals


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
        type(integrals) :: int
        logical :: ehlers_calc
        integer :: i

        ehlers_calc = .false.

        ! Specify influencing vertices (also sets zero default influence)

        ! Source
        if (source_order == 0) then
            allocate(phi_s(1), source=0.)
            allocate(i_vert_s(1), source=this%index)
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
            if (freestream%supersonic) then
                geom = this%calc_supersonic_subinc_geom(P, freestream, dod_info)
            else
                geom = this%calc_subsonic_geom(P, freestream)
            end if

            ! Get integrals
            int = this%calc_integrals(geom, 'potential', 'doublet', dod_info, freestream)

            ! Source potential
            if (source_order == 0) then

                if (freestream%supersonic) then

                    if (ehlers_calc) then
                        
                        ! Sum up w0 terms
                        do i=1,this%N
                            if (abs(this%m(i)) < 1.) then
                                phi_s = phi_s + geom%xm(i)*int%w0(i)*this%m(i)
                            else
                                phi_s = phi_s + geom%xm(i)*int%w0(i)
                            end if
                        end do

                        ! Add Q1 term
                        phi_s = this%J*freestream%K_inv*(phi_s - geom%h*int%Q1)

                    else

                        ! Equivalent to Ehlers Eq. (8.6)
                        phi_s = this%J*freestream%K_inv*(geom%h*int%hH113 + sum(geom%a*int%F111))

                    end if
                else

                    ! Johnson Eq. (D21) including the area factor discussed by Ehlers in Sec. 10.3
                    phi_s = this%J*freestream%K_inv*(geom%h*int%hH113 - sum(geom%a*int%F111))

                end if

            end if

            ! Doublet potential
            if (doublet_order == 1) then

                if (freestream%supersonic) then

                    if (ehlers_calc) then

                        ! Sum up w0 terms (Ehlers Eq. (5.17))
                        do i=1,this%N
                            if (abs(this%m(i)) < 1.) then
                                phi_d(2) = phi_d(2) + int%w0(i)*this%m(i)
                                phi_d(3) = phi_d(3) + int%w0(i)
                            else
                                phi_d(2) = phi_d(2) + int%w0(i)
                                phi_d(3) = phi_d(3) + int%w0(i)*this%l(i)
                            end if
                        end do

                        ! Add Q1 terms (Ehlers Eq. (5.17))
                        phi_d(1) = -int%Q1
                        phi_d(2) = geom%h*phi_d(2) - int%Q1*geom%P_ls(1)
                        phi_d(3) = geom%h*phi_d(3) - int%Q1*geom%P_ls(2)

                    else

                        ! Equivalent to Ehlers Eq. (5.17))
                        phi_d(1) = int%hH113
                        phi_d(2) = int%hH113*geom%P_ls(1) - geom%h*sum(this%n_hat_ls(1,:)*int%F111)
                        phi_d(3) = int%hH113*geom%P_ls(1) + geom%h*sum(this%n_hat_ls(2,:)*int%F111)

                    end if
                else

                    ! Johnson Eq. (D.30)
                    phi_d(1) = int%hH113
                    phi_d(2) = int%hH113*geom%P_ls(1) + geom%h*int%H213
                    phi_d(3) = int%hH113*geom%P_ls(2) + geom%h*int%H123
                end if

                ! Convert to vertex influences (Davis Eq. (4.41))
                phi_d(1:3) = freestream%K_inv*matmul(phi_d(1:3), this%S_mu_inv)

                ! Wake bottom influence is opposite the top influence
                if (this%in_wake) then
                    phi_d(4:6) = -phi_d(1:3)
                end if
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
        type(integrals) :: int

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
            if (freestream%supersonic) then
                geom = this%calc_supersonic_subinc_geom(P, freestream, dod_info)
            else
                geom = this%calc_subsonic_geom(P, freestream)
            end if

            ! Get integrals
            int = this%calc_integrals(geom, "velocity", "doublet", dod_info, freestream)

            ! Source velocity
            if (source_order == 0) then
                v_s(1,1) = this%J*freestream%K_inv*sum(this%n_hat_ls(1,:)*int%F111(:))
                v_s(1,2) = this%J*freestream%K_inv*sum(this%n_hat_ls(2,:)*int%F111(:))
                v_s(1,3) = this%J*freestream%K_inv*int%hH113
            end if

            ! Doublet velocity
            if (doublet_order == 1) then
            end if
        end if
    
    end subroutine panel_calc_velocities


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