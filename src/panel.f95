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


    type eval_point_geom
        ! Container type for the geometric parameters necessary for calculating a panel's influence on a given field point

        real,dimension(3) :: P_g ! Point position in global coords
        real,dimension(2) :: P_ls ! Transformed point in panel plane
        real :: h, h2 ! Transformed height above panel
        real,dimension(3) :: a, g2, l1, l2, R1, R2, dR ! Edge integration parameters for the Ehlers-Johnson method
        real,dimension(3) :: v_xi, v_eta ! Edge in-plane normal vectors

        contains

            procedure :: init => eval_point_geom_init

    end type eval_point_geom


    type integrals
        ! Container type for the fundamental integrals used to calculate influence coefficients

        real :: H111, H211, H121 ! Source integrals
        real :: hH113, H213, H123, H313, H223, H133 ! Doublet integrals; we use hH(1,1,3) because it can be reliably calculated, unlike H(1,1,3)
        real,dimension(:),allocatable :: F111, F211, F121 ! Necessary line integrals

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
        type(vertex_pointer),dimension(:),allocatable :: midpoints
        real,dimension(3) :: n_g, nu_g ! Normal and conormal vectors
        real,dimension(3) :: n_g_mir, nu_g_mir ! Mirrored normal and conormal vectors
        real,dimension(3) :: centr, centr_mir ! Centroid
        real,dimension(3,3) :: A_g_to_ls, A_ls_to_g ! Coordinate transformation matrices
        real,dimension(3,3) :: A_g_to_ls_mir, A_ls_to_g_mir
        real,dimension(:,:),allocatable :: vertices_ls, midpoints_ls ! Location of the vertices and edge midpoints described in local scaled coords
        real,dimension(:,:),allocatable :: vertices_ls_mir, midpoints_ls_mir
        real,dimension(:,:),allocatable :: t_hat_ls ! Edge unit tangents
        real,dimension(:,:),allocatable :: t_hat_ls_mir
        real,dimension(:,:),allocatable :: n_hat_g, n_hat_ls ! Edge unit outward normals
        real,dimension(:,:),allocatable :: n_hat_g_mir, n_hat_ls_mir
        real,dimension(:),allocatable :: b, sqrt_b ! Edge parameter
        real,dimension(:),allocatable :: b_mir, sqrt_b_mir
        real :: A ! Surface area (same for mirror, in global coordinates at least)
        real,dimension(:,:),allocatable :: S_mu_inv, S_sigma_inv ! Matrix relating doublet/source strengths to doublet/source influence parameters
        real,dimension(:,:),allocatable :: S_mu_inv_mir, S_sigma_inv_mir
        logical :: in_wake ! Whether this panel belongs to a wake mesh
        integer,dimension(3) :: abutting_panels ! Indices of panels abutting this one
        integer :: r, r_mir ! Panel inclination indicator; r=-1 -> superinclined, r=1 -> subinclined
        real :: J, J_mir ! Local scaled transformation Jacobian
        integer,dimension(:),allocatable :: i_vert_d, i_vert_s

        contains

            ! Initialization procedures
            procedure :: init => panel_init_3
            procedure :: calc_derived_geom => panel_calc_derived_geom
            procedure :: calc_normal => panel_calc_normal
            procedure :: calc_area => panel_calc_area
            procedure :: calc_centroid => panel_calc_centroid

            ! Flow-dependent initialization procedures
            procedure :: init_with_flow => panel_init_with_flow
            procedure :: calc_g_to_ls_transform => panel_calc_g_to_ls_transform
            procedure :: calc_edge_vectors => panel_calc_edge_vectors
            procedure :: calc_singularity_matrices => panel_calc_singularity_matrices
            procedure :: set_influencing_verts => panel_set_influencing_verts

            ! Mirror initialization
            procedure :: init_mirror => panel_init_mirror
            procedure :: calc_mirrored_g_to_ls_transform => panel_calc_mirrored_g_to_ls_transform
            procedure :: calc_mirrored_edge_vectors => panel_calc_mirrored_edge_vectors
            procedure :: calc_mirrored_singularity_matrices => panel_calc_mirrored_singularity_matrices

            ! Getters
            procedure :: get_vertex_loc => panel_get_vertex_loc
            procedure :: get_midpoint_loc => panel_get_midpoint_loc
            procedure :: get_vertex_index => panel_get_vertex_index
            procedure :: get_midpoint_index => panel_get_midpoint_index
            procedure :: touches_vertex => panel_touches_vertex

            ! Update information
            procedure :: point_to_new_vertex => panel_point_to_new_vertex

            ! Influence calculations
            procedure :: check_dod => panel_check_dod
            procedure :: calc_subsonic_geom => panel_calc_subsonic_geom
            procedure :: calc_supersonic_subinc_geom => panel_calc_supersonic_subinc_geom
            procedure :: E_i_M_N_K => panel_E_i_M_N_K
            procedure :: calc_subsonic_edge_integrals => panel_calc_subsonic_edge_integrals
            procedure :: calc_subsonic_panel_integrals => panel_calc_subsonic_panel_integrals
            procedure :: calc_supersonic_subinc_edge_integrals => panel_calc_supersonic_subinc_edge_integrals
            procedure :: calc_supersonic_subinc_panel_integrals => panel_calc_supersonic_subinc_panel_integrals
            procedure :: calc_integrals => panel_calc_integrals
            procedure :: allocate_potential_influences => panel_allocate_potential_influences
            procedure :: calc_potentials => panel_calc_potentials
            procedure :: calc_velocities => panel_calc_velocities
            procedure :: get_velocity_jump => panel_get_velocity_jump

    end type panel

    
contains


    subroutine eval_point_geom_init(this, P, A_g_to_ls, centr)

        implicit none

        class(eval_point_geom),intent(out) :: this
        real,dimension(3),intent(in) :: P, centr
        real,dimension(3,3),intent(in) :: A_g_to_ls

        real,dimension(3) :: x

        ! Store point
        this%P_g = P

        ! Transform to local scaled coordinates
        x = matmul(A_g_to_ls, this%P_g-centr)
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

        ! Allocate vertex arrays
        allocate(this%vertices(this%N))

        ! Assign vertex pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3

        ! Store the index of the panel
        this%index = index

        ! Initialize a few things
        this%abutting_panels = 0

        ! Allocate midpoint arrays
        if (doublet_order == 2) then
            allocate(this%midpoints(this%N))
        end if

        call this%calc_derived_geom()

    end subroutine panel_init_3


    subroutine  panel_calc_derived_geom(this)
        ! Initializes geometry based on the location of the vertices.
        ! Should be called when panel geometry is updated.

        implicit none

        class(panel),intent(inout) :: this

        ! Calculate midpoints

        ! Calculate normal vec
        call this%calc_normal()

        ! Calculate area
        call this%calc_area()

        ! Calculate centroid
        call this%calc_centroid()

    end subroutine  panel_calc_derived_geom


    subroutine panel_calc_normal(this)

        implicit none

        class(panel),intent(inout) :: this

        real,dimension(3) :: d1, d2

        ! Get two edge vectors
        d1 = this%get_vertex_loc(2)-this%get_vertex_loc(1)
        d2 = this%get_vertex_loc(3)-this%get_vertex_loc(2)

        ! Find normal
        this%n_g = cross(d1, d2)
        this%n_g = this%n_g/norm2(this%n_g)

    end subroutine panel_calc_normal


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
                write(*,*) "!!! Vertex 1: ", this%get_vertex_loc(1)
                write(*,*) "!!! Vertex 2: ", this%get_vertex_loc(2)
                write(*,*) "!!! Vertex 3: ", this%get_vertex_loc(3)
                stop
            end if

        end if

    end subroutine panel_calc_area


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
        this%centr = sum/this%N

    end subroutine panel_calc_centroid


    subroutine panel_init_with_flow(this, freestream, initialize_mirror, mirror_plane)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream
        logical,intent(in) :: initialize_mirror
        integer,intent(in) :: mirror_plane

        ! Calculate transforms
        call this%calc_g_to_ls_transform(freestream)

        ! Calculate properties dependent on the transforms
        call this%calc_edge_vectors(freestream)
        call this%calc_singularity_matrices()

        ! Calculate mirrored properties
        if (initialize_mirror) then
            call this%init_mirror(freestream, mirror_plane)
        end if

    end subroutine panel_init_with_flow


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
        if (doublet_order == 2) allocate(this%midpoints_ls(2,this%N))
        do i=1,this%N

            ! Vertices
            this%vertices_ls(:,i) = matmul(this%A_g_to_ls(1:2,:), this%get_vertex_loc(i)-this%centr)

            ! Midpoints
            if (doublet_order == 2) then
                this%midpoints_ls(:,i) = matmul(this%A_g_to_ls(1:2,:), this%get_midpoint_loc(i)-this%centr)
            end if

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

        real,dimension(3) :: d_g, t_hat_g
        real,dimension(2) :: d_ls
        integer :: i, i_next

        ! Allocate memory
        allocate(this%n_hat_g(3,this%N))
        allocate(this%t_hat_ls(2,this%N))
        allocate(this%n_hat_ls(2,this%N))
        allocate(this%b(this%N))
        allocate(this%b_mir(this%N)) ! This needs to be initialized here because of some DoD checks. It will have no effect.
        allocate(this%sqrt_b(this%N))

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N)+1

            ! Calculate edge vector based on index
            d_g = this%get_vertex_loc(i_next)-this%get_vertex_loc(i)

            ! Calculate tangent in global coords
            t_hat_g = d_g/norm2(d_g)

            ! Calculate tangent in local scaled coords 
            d_ls = this%vertices_ls(:,i_next) - this%vertices_ls(:,i)
            this%t_hat_ls(:,i) = d_ls/norm2(d_ls)

            ! Calculate edge outward normal
            this%n_hat_g(:,i) = cross(t_hat_g, this%n_g)

        end do

        ! Calculate edge normal in local scaled coords E&M Eq. (J.6.45)
        this%n_hat_ls(1,:) = this%t_hat_ls(2,:)
        this%n_hat_ls(2,:) = -this%t_hat_ls(1,:)

        ! Calculate edge parameter (Ehlers Eq. (E14))
        this%b = (this%n_hat_ls(1,:) - this%n_hat_ls(2,:))*(this%n_hat_ls(1,:) + this%n_hat_ls(2,:))
        this%sqrt_b = sqrt(abs(this%b))
    
    end subroutine panel_calc_edge_vectors


    subroutine panel_calc_singularity_matrices(this)
        ! Calculates the matrices which relate the singularity strengths to the singularity parameters

        implicit none

        class(panel),intent(inout) :: this

        real,dimension(:,:),allocatable :: S_mu, S_sigma

        ! Determine influence of vertex doublet strengths on integral parameters

        ! Linear distribution
        if (doublet_order == 1) then

            ! Allocate influence matrices
            allocate(S_mu(3,3))
            allocate(this%S_mu_inv(3,3))

            ! Set values
            S_mu(:,1) = 1.
            S_mu(:,2) = this%vertices_ls(1,:)
            S_mu(:,3) = this%vertices_ls(2,:)

            ! Invert
            call matinv(3, S_mu, this%S_mu_inv)

        else if (doublet_order == 2) then

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

        ! Determine influence of vertex source strengths on integral parameters
        ! Linear distribution
        if (source_order == 1) then

            ! Allocate influence matrices
            allocate(S_sigma(3,3))
            allocate(this%S_sigma_inv(3,3))

            ! Set values
            S_sigma(:,1) = 1.
            S_sigma(:,2) = this%vertices_ls(1,:)
            S_sigma(:,3) = this%vertices_ls(2,:)

            ! Invert
            call matinv(3, S_sigma, this%S_sigma_inv)

        end if
    
    end subroutine panel_calc_singularity_matrices


    subroutine panel_set_influencing_verts(this)
        ! Sets up the arrays of vertex indices which set the influence for this panel

        class(panel),intent(inout) :: this

        ! Source
        if (source_order == 0) then
            allocate(this%i_vert_s(1), source=this%index)
        else if (source_order == 1) then
            allocate(this%i_vert_s(3))
            this%i_vert_s(1) = this%get_vertex_index(1)
            this%i_vert_s(2) = this%get_vertex_index(2)
            this%i_vert_s(3) = this%get_vertex_index(3)
        end if

        ! Doublet
        if (doublet_order == 1) then

            ! Check if this panel belongs to the wake
            if (this%in_wake) then

                ! Wake panels are influenced by two sets of vertices
                allocate(this%i_vert_d(6))
                this%i_vert_d(1) = this%vertices(1)%ptr%top_parent
                this%i_vert_d(2) = this%vertices(2)%ptr%top_parent
                this%i_vert_d(3) = this%vertices(3)%ptr%top_parent
                this%i_vert_d(4) = this%vertices(1)%ptr%bot_parent
                this%i_vert_d(5) = this%vertices(2)%ptr%bot_parent
                this%i_vert_d(6) = this%vertices(3)%ptr%bot_parent

            else

                ! Body panels are influenced by only one set of vertices
                allocate(this%i_vert_d(3))
                this%i_vert_d(1) = this%get_vertex_index(1)
                this%i_vert_d(2) = this%get_vertex_index(2)
                this%i_vert_d(3) = this%get_vertex_index(3)

            end if

        else if (doublet_order == 2) then

            ! Check if this panel belongs to the wake
            if (this%in_wake) then

                ! Wake panels are influenced by two sets of vertices
                allocate(this%i_vert_d(12))
                this%i_vert_d(1) = this%vertices(1)%ptr%top_parent
                this%i_vert_d(2) = this%vertices(2)%ptr%top_parent
                this%i_vert_d(3) = this%vertices(3)%ptr%top_parent
                this%i_vert_d(4) = this%midpoints(1)%ptr%top_parent
                this%i_vert_d(5) = this%midpoints(2)%ptr%top_parent
                this%i_vert_d(6) = this%midpoints(3)%ptr%top_parent
                this%i_vert_d(7) = this%vertices(1)%ptr%bot_parent
                this%i_vert_d(8) = this%vertices(2)%ptr%bot_parent
                this%i_vert_d(9) = this%vertices(3)%ptr%bot_parent
                this%i_vert_d(10) = this%midpoints(1)%ptr%bot_parent
                this%i_vert_d(11) = this%midpoints(2)%ptr%bot_parent
                this%i_vert_d(12) = this%midpoints(3)%ptr%bot_parent

            else

                ! Body panels are influenced by only one set of vertices
                allocate(this%i_vert_d(6))
                this%i_vert_d(1) = this%get_vertex_index(1)
                this%i_vert_d(2) = this%get_vertex_index(2)
                this%i_vert_d(3) = this%get_vertex_index(3)
                this%i_vert_d(4) = this%get_midpoint_index(1)
                this%i_vert_d(5) = this%get_midpoint_index(2)
                this%i_vert_d(6) = this%get_midpoint_index(3)

            end if

        end if
    
        
    end subroutine panel_set_influencing_verts


    subroutine panel_init_mirror(this, freestream, mirror_plane)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream
        integer,intent(in) :: mirror_plane

        integer :: i

        ! Calculate mirrored normal vector
        this%n_g_mir = mirror_across_plane(this%n_g, mirror_plane)

        ! Calculate mirrored centroid
        this%centr_mir = mirror_across_plane(this%centr, mirror_plane)

        ! Calculate mirrored g to ls transform
        call this%calc_mirrored_g_to_ls_transform(freestream, mirror_plane)

        ! Calculate mirrored edge vectors
        ! Global
        allocate(this%n_hat_g_mir(3,this%N))
        do i=1,this%N
            this%n_hat_g_mir(:,i) = mirror_across_plane(this%n_hat_g(:,i), mirror_plane)
        end do

        ! Local-scaled
        call this%calc_mirrored_edge_vectors()

        ! Calculate mirrored singularity matrices
        call this%calc_mirrored_singularity_matrices()

    end subroutine panel_init_mirror


    subroutine panel_calc_mirrored_g_to_ls_transform(this, freestream, mirror_plane)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream
        integer,intent(in) :: mirror_plane

        real,dimension(3) :: u0, v0
        real,dimension(3,3) :: B_mat_ls
        real :: x, y
        integer :: i, rs

        ! Get in-panel basis vectors
        if (abs(inner(this%n_g_mir, freestream%c_hat_g) - 1.) < 1e-12) then ! Check the freestream isn't aligned with the normal vector
            v0 = mirror_across_plane(this%get_vertex_loc(2) - this%get_vertex_loc(1), mirror_plane)
        else
            v0 = cross(this%n_g_mir, freestream%c_hat_g)
        end if
        v0 = v0/norm2(v0)
        u0 = cross(v0, this%n_g_mir)
        u0 = u0/norm2(u0)

        ! Calculate compressible parameters
        this%nu_g_mir = matmul(freestream%B_mat_g, this%n_g_mir)
        x = inner(this%n_g_mir, this%nu_g_mir)

        ! Check for Mach-inclined panels
        if (freestream%supersonic .and. abs(x) < 1e-12) then
            write(*,*) "!!! Mirror of panel", this%index, "is Mach-inclined, which is not allowed. Quitting..."
            stop
        end if

        ! Calculate panel inclination indicator (E&M Eq. (E.3.16b))
        this%r_mir = sign(1., x) ! r=-1 -> superinclined, r=1 -> subinclined

        ! Check for superinclined panels
        if (this%r_mir < 0) then
            write(*,*) "!!! Mirror of panel", this%index, "is superinclined, which is not allowed. Quitting..."
            stop
        end if

        ! Other inclination parameters
        rs = this%r*freestream%s

        ! Calculate transformation
        y = 1./sqrt(abs(x))
        this%A_g_to_ls_mir(1,:) = y*matmul(freestream%C_mat_g, u0)
        this%A_g_to_ls_mir(2,:) = rs/freestream%B*matmul(freestream%C_mat_g, v0)
        this%A_g_to_ls_mir(3,:) = freestream%B*y*this%n_g_mir

        ! Check determinant
        x = det3(this%A_g_to_ls_mir)
        if (abs(x-freestream%B**2) >= 1e-12) then
            write(*,*) "!!! Calculation of mirrored local scaled coordinate transform failed. Quitting..."
            stop
        end if

        ! Calculate inverse
        if (freestream%M_inf == 0.) then
            this%A_ls_to_g_mir = transpose(this%A_g_to_ls_mir)
        else
            call matinv(3, this%A_g_to_ls_mir, this%A_ls_to_g_mir)
        end if

        ! Calculate Jacobian
        this%J_mir = 1./(freestream%B*sqrt(abs(1.-freestream%M_inf**2*inner(freestream%c_hat_g, this%n_g_mir)**2)))

        ! Transform vertex and midpoint coords to ls
        allocate(this%vertices_ls_mir(2,this%N))
        if (doublet_order == 2) allocate(this%midpoints_ls_mir(2,this%N))
        do i=1,this%N

            ! Vertices
            this%vertices_ls_mir(:,i) = matmul(this%A_g_to_ls_mir(1:2,:), &
                                               mirror_across_plane(this%get_vertex_loc(i), mirror_plane)-this%centr_mir)

            ! Midpoints
            if (doublet_order == 2) then
                this%midpoints_ls_mir(:,i) = matmul(this%A_g_to_ls_mir(1:2,:), &
                                                    mirror_across_plane(this%get_midpoint_loc(i), mirror_plane)-this%centr_mir)
            end if

        end do
    
    end subroutine panel_calc_mirrored_g_to_ls_transform


    subroutine panel_calc_mirrored_edge_vectors(this)

        implicit none

        class(panel),intent(inout) :: this

        real,dimension(2) :: d_ls
        integer :: i, i_next

        ! Allocate memory
        allocate(this%t_hat_ls_mir(2,this%N))
        allocate(this%n_hat_ls_mir(2,this%N))
        allocate(this%sqrt_b_mir(this%N))

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N)+1

            ! Calculate tangent in local scaled coords 
            d_ls = this%vertices_ls_mir(:,i_next) - this%vertices_ls_mir(:,i)
            this%t_hat_ls_mir(:,i) = d_ls/norm2(d_ls)

        end do

        ! Calculate edge normal in local scaled coords E&M Eq. (J.6.45)
        this%n_hat_ls_mir(1,:) = this%t_hat_ls_mir(2,:)
        this%n_hat_ls_mir(2,:) = -this%t_hat_ls_mir(1,:)

        ! Calculate edge parameter (Ehlers Eq. (E14))
        this%b_mir = (this%n_hat_ls_mir(1,:) - this%n_hat_ls_mir(2,:))*(this%n_hat_ls_mir(1,:) + this%n_hat_ls_mir(2,:))
        this%sqrt_b_mir = sqrt(abs(this%b_mir))
    
    end subroutine panel_calc_mirrored_edge_vectors


    subroutine panel_calc_mirrored_singularity_matrices(this)

        implicit none

        class(panel),intent(inout) :: this

        real,dimension(:,:),allocatable :: S_mu, S_sigma

        ! Linear distribution
        if (doublet_order == 1) then

            ! Allocate influence matrices
            allocate(S_mu(3,3))
            allocate(this%S_mu_inv_mir(3,3))

            ! Set values
            S_mu(:,1) = 1.
            S_mu(:,2) = this%vertices_ls_mir(1,:)
            S_mu(:,3) = this%vertices_ls_mir(2,:)

            ! Invert
            call matinv(3, S_mu, this%S_mu_inv_mir)

        else if (doublet_order == 2) then

            ! Allocate influence matrix
            allocate(S_mu(6,6))
            allocate(this%S_mu_inv_mir(6,6))

            ! Set values
            S_mu(:,1) = 1.

            S_mu(1:3,2) = this%vertices_ls_mir(1,:)
            S_mu(1:3,3) = this%vertices_ls_mir(2,:)
            S_mu(1:3,4) = this%vertices_ls_mir(1,:)**2
            S_mu(1:3,5) = this%vertices_ls_mir(1,:)*this%vertices_ls_mir(2,:)
            S_mu(1:3,6) = this%vertices_ls_mir(2,:)**2
            
            S_mu(4:6,2) = this%midpoints_ls_mir(1,:)
            S_mu(4:6,3) = this%midpoints_ls_mir(2,:)
            S_mu(4:6,4) = this%midpoints_ls_mir(1,:)**2
            S_mu(4:6,5) = this%midpoints_ls_mir(1,:)*this%midpoints_ls_mir(2,:)
            S_mu(4:6,6) = this%midpoints_ls_mir(2,:)**2

            ! Invert
            call matinv(6, S_mu, this%S_mu_inv_mir)

        end if

        ! Determine influence of vertex source strengths on integral parameters
        ! Linear distribution
        if (source_order .eq. 1) then

            ! Allocate influence matrices
            allocate(S_sigma(3,3))
            allocate(this%S_sigma_inv_mir(3,3))

            ! Set values
            S_sigma(:,1) = 1.
            S_sigma(:,2) = this%vertices_ls_mir(1,:)
            S_sigma(:,3) = this%vertices_ls_mir(2,:)

            ! Invert
            call matinv(3, S_sigma, this%S_sigma_inv_mir)

        else if (source_order .eq. 2) then

            ! Allocate influence matrix
            allocate(S_sigma(6,6))
            allocate(this%S_sigma_inv_mir(6,6))

            ! Set values
            S_sigma(:,1) = 1.

            S_sigma(1:3,2) = this%vertices_ls_mir(1,:)
            S_sigma(1:3,3) = this%vertices_ls_mir(2,:)
            S_sigma(1:3,4) = this%vertices_ls_mir(1,:)**2
            S_sigma(1:3,5) = this%vertices_ls_mir(1,:)*this%vertices_ls_mir(2,:)
            S_sigma(1:3,6) = this%vertices_ls_mir(2,:)**2
            
            S_sigma(4:6,2) = this%midpoints_ls_mir(1,:)
            S_sigma(4:6,3) = this%midpoints_ls_mir(2,:)
            S_sigma(4:6,4) = this%midpoints_ls_mir(1,:)**2
            S_sigma(4:6,5) = this%midpoints_ls_mir(1,:)*this%midpoints_ls_mir(2,:)
            S_sigma(4:6,6) = this%midpoints_ls_mir(2,:)**2

            ! Invert
            call matinv(6, S_sigma, this%S_sigma_inv_mir)

        end if
    
    end subroutine panel_calc_mirrored_singularity_matrices


    function panel_get_vertex_loc(this, i) result(loc)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        real,dimension(3) :: loc

        loc = this%vertices(i)%ptr%loc

    end function panel_get_vertex_loc


    function panel_get_midpoint_loc(this, i) result(loc)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        real,dimension(3) :: loc

        loc = this%midpoints(i)%ptr%loc

    end function panel_get_midpoint_loc


    function panel_get_vertex_index(this, i) result(index)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        integer :: index

        index = this%vertices(i)%ptr%index

    end function panel_get_vertex_index


    function panel_get_midpoint_index(this, i) result(index)

        implicit none

        class(panel),intent(in) :: this
        integer,intent(in) :: i
        integer :: index

        index = this%midpoints(i)%ptr%index

    end function panel_get_midpoint_index


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
            if (this%get_vertex_index(j) == i) then
                touches = .true.
                return
            end if

        end do

    end function panel_touches_vertex


    subroutine panel_point_to_new_vertex(this, new_vertex)
        ! Updates the panel to point to this new vertex (assumed to be a copy of a current vertex)

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(in),target :: new_vertex
        integer :: i

        ! Loop through vertices/midpoints
        do i=1,this%N

            ! It's a vertex, so check the vertex locations
            if (new_vertex%vert_type == 1) then

                if (dist(this%get_vertex_loc(i), new_vertex%loc) < 1e-12) then

                    ! Update pointer
                    this%vertices(i)%ptr => new_vertex

                    return

                end if

            ! Check midpoint locations
            else

                if (dist(this%get_midpoint_loc(i), new_vertex%loc) < 1e-12) then

                    ! Update pointer
                    this%midpoints(i)%ptr => new_vertex

                    return

                end if

            end if

        end do
    
    end subroutine panel_point_to_new_vertex


    function panel_check_dod(this, eval_point, freestream, verts_in_dod, mirror_panel, mirror_plane) result(dod_info)
        ! Determines how (if) this panel lies within the domain of dependence of the evaluation point

        implicit none

        class(panel),intent(inout) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        logical,dimension(:),intent(in) :: verts_in_dod
        logical,intent(in),optional :: mirror_panel
        integer,intent(in),optional :: mirror_plane

        type(dod) :: dod_info

        real,dimension(3) :: d, point, a, b, R_star, Q_end
        integer :: i, i_next
        real :: x, s_star
        logical :: mirrored, in_panel

        ! Set default mirroring
        if (present(mirror_panel)) then
            mirrored = mirror_panel
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

                        ! For a subsonic or sonic edge, both being out means the edge is out
                        if ((.not. mirrored .and. this%b(i) <= 0.) .or. (mirrored .and. this%b_mir(i) <= 0.)) then
                            dod_info%edges_in_dod(i) = .false.

                        ! For a supersonic edge, the edge can still intersect the DoD, so calculate the point of closest approach
                        else

                            ! Get end vertex and vector describing edge
                            if (mirrored) then
                                Q_end = mirror_across_plane(this%get_vertex_loc(i_next), mirror_plane)
                                d = Q_end - mirror_across_plane(this%get_vertex_loc(i), mirror_plane)
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
                    if (mirrored) then
                        s_star = inner(mirror_across_plane(this%get_vertex_loc(1), mirror_plane) - eval_point, this%n_g) &
                                 / inner(freestream%c_hat_g, this%n_g)
                    else
                        s_star = inner(this%get_vertex_loc(1)-eval_point, this%n_g)/inner(freestream%c_hat_g, this%n_g)
                    end if
                    R_star = eval_point + freestream%c_hat_g*s_star

                    ! See if the projected point is in the panel
                    in_panel = .true.
                    do i=1,this%N

                        ! Get edge displacement
                        if (mirrored) then
                            x = inner(R_star, this%n_hat_g_mir(:,i))
                        else
                            x = inner(R_star, this%n_hat_g(:,i))
                        end if

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


    function panel_calc_subsonic_geom(this, eval_point, freestream, mirror_panel) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point in subsonic flow.

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(eval_point_geom) :: geom

        real,dimension(2) :: d_ls
        real,dimension(3) :: d_g
        integer :: i, i_next
        real :: val
        real,dimension(this%N) :: dummy

        ! Initialize
        if (mirror_panel) then
            call geom%init(eval_point, this%A_g_to_ls_mir, this%centr_mir)
        else
            call geom%init(eval_point, this%A_g_to_ls, this%centr)
        end if

        ! Store edge in-plane normal vectors
        if (mirror_panel) then
            geom%v_xi = this%n_hat_ls_mir(1,:)
            geom%v_eta = this%n_hat_ls_mir(2,:)
        else
            geom%v_xi = this%n_hat_ls(1,:)
            geom%v_eta = this%n_hat_ls(2,:)
        end if

        ! Calculate edge quantities
        do i=1,this%N

            i_next = mod(i, this%N) + 1

            ! Calculate displacement vector from start vertex
            if (mirror_panel) then
                d_ls = this%vertices_ls_mir(:,i) - geom%P_ls
            else
                d_ls = this%vertices_ls(:,i) - geom%P_ls
            end if

            ! Perpendicular distance in plane from evaluation point to edge
            geom%a(i) = d_ls(1)*geom%v_xi(i) + d_ls(2)*geom%v_eta(i)

            ! Integration length on edge to start vertex
            geom%l1(i) = -d_ls(1)*geom%v_eta(i) + d_ls(2)*geom%v_xi(i)

            ! Distance from evaluation point to start vertex
            geom%R1(i) = sqrt(d_ls(1)**2 + d_ls(2)**2 + geom%h2)

            ! Calculate displacement vector from end vertex
            if (mirror_panel) then
                d_ls = this%vertices_ls_mir(:,i_next) - geom%P_ls
            else
                d_ls = this%vertices_ls(:,i_next) - geom%P_ls
            end if

            ! Integration length on edge to start vertex
            geom%l2(i) = -d_ls(1)*geom%v_eta(i) + d_ls(2)*geom%v_xi(i)

        end do

        ! Square of the perpendicular distance to edge
        geom%g2 = geom%a**2 + geom%h2

        ! Distance from evaluation point to end vertices
        geom%R2 = cshift(geom%R1, 1)

        ! Swap directions for mirror
        if (mirror_panel) then
            dummy = geom%l1
            geom%l1 = geom%l2
            geom%l2 = dummy

            dummy = geom%R1
            geom%R1 = geom%R2
            geom%R2 = dummy
        end if

        ! Difference in R
        geom%dR = geom%R2 - geom%R1

    end function panel_calc_subsonic_geom


    function panel_calc_supersonic_subinc_geom(this, eval_point, freestream, mirror_panel, dod_info) result(geom)
        ! Calculates the geometric parameters necessary for calculating the influence of the panel at the given evaluation point

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(dod),intent(in) :: dod_info
        type(eval_point_geom) :: geom

        real,dimension(2) :: d_ls, d
        real :: x
        integer :: i, i_next
        real,dimension(this%N) :: dummy

        ! Initialize
        if (mirror_panel) then
            call geom%init(eval_point, this%A_g_to_ls_mir, this%centr_mir)
            geom%v_xi = this%n_hat_ls_mir(1,:)
            geom%v_eta = this%n_hat_ls_mir(2,:)
        else
            call geom%init(eval_point, this%A_g_to_ls, this%centr)
            geom%v_xi = this%n_hat_ls(1,:)
            geom%v_eta = this%n_hat_ls(2,:)
        end if

        ! Loop through edges
        do i=1,this%N

            ! Check DoD
            if (dod_info%edges_in_dod(i)) then

                ! Calculate displacement from first vertex
                if (mirror_panel) then
                    d_ls = this%vertices_ls_mir(:,i) - geom%P_ls
                else
                    d_ls = this%vertices_ls(:,i) - geom%P_ls
                end if

                ! Edge integration length
                geom%l1(i) = geom%v_eta(i)*d_ls(1) + geom%v_xi(i)*d_ls(2)

                ! Perpendicular in-plane distance
                geom%a(i) = geom%v_xi(i)*d_ls(1) + geom%v_eta(i)*d_ls(2)

                ! Perpendicular hyperbolic distance
                if (mirror_panel) then
                    geom%g2(i) = geom%a(i)**2 - this%b_mir(i)*geom%h2
                else
                    geom%g2(i) = geom%a(i)**2 - this%b(i)*geom%h2
                end if

                ! Hyperbolic radius to first vertex
                x = d_ls(1)*d_ls(1) - d_ls(2)*d_ls(2) - geom%h2
                if (x > 0. .and. d_ls(1) < 0.) then
                    geom%R1(i) = sqrt(x)
                else
                    geom%l1(i) = -sqrt(abs(geom%g2(i)))
                    geom%R1(i) = 0.
                end if

                ! Get index of end vertex
                i_next = mod(i, this%N)+1

                ! Calculate displacement from second vertex
                if (mirror_panel) then
                    d_ls = this%vertices_ls_mir(:,i_next) - geom%P_ls
                else
                    d_ls = this%vertices_ls(:,i_next) - geom%P_ls
                end if

                ! Edge integration length
                geom%l2(i) = geom%v_eta(i)*d_ls(1) + geom%v_xi(i)*d_ls(2)

                ! Hyperbolic radius to first vertex
                x = d_ls(1)*d_ls(1) - d_ls(2)*d_ls(2) - geom%h2
                if (x > 0. .and. d_ls(1) < 0.) then
                    geom%R2(i) = sqrt(x)
                else
                    geom%l2(i) = sqrt(abs(geom%g2(i)))
                    geom%R2(i) = 0.
                end if

            end if

        end do

        ! Swap directions for mirror
        if (mirror_panel) then
            dummy = geom%l1
            geom%l1 = geom%l2
            geom%l2 = dummy

            dummy = geom%R1
            geom%R1 = geom%R2
            geom%R2 = dummy
        end if

        ! Difference in R
        geom%dR = geom%R2 - geom%R1

    end function panel_calc_supersonic_subinc_geom


    function panel_E_i_M_N_K(this, geom, i, M, N, K, freestream, mirror_panel) result(E)
        ! Calculates E_i(M,N,K)

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        integer,intent(in) :: i, M, N, K
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel

        real :: E

        real :: E1, E2
        integer :: i_next

        i_next = mod(i, this%N) + 1

        if (.not. mirror_panel) then

            ! Evaluate at start vertex
            E1 = ((this%vertices_ls(1,i)-geom%P_ls(1))**(M-1)*(this%vertices_ls(2,i)-geom%P_ls(2))**(N-1))/geom%R1(i)**K

            ! Evaluate at end vertex
            E2 = ((this%vertices_ls(1,i_next)-geom%P_ls(1))**(M-1)*(this%vertices_ls(2,i_next)-geom%P_ls(2))**(N-1))/geom%R2(i)**K

        else

            ! Evaluate at start vertex (order is swapped due to mirroring)
            E2 = ((this%vertices_ls_mir(1,i)-geom%P_ls(1))**(M-1)*(this%vertices_ls_mir(2,i)-geom%P_ls(2))**(N-1))/geom%R1(i)**K

            ! Evaluate at end vertex
            E1 = ((this%vertices_ls_mir(1,i_next)-geom%P_ls(1))**(M-1)*(this%vertices_ls_mir(2,i_next)-geom%P_ls(2))**(N-1)) &
                 /geom%R2(i)**K

        end if

        ! Calculate difference
        E = E2 - E1

    end function panel_E_i_M_N_K


    subroutine panel_calc_subsonic_edge_integrals(this, geom, freestream, mirror_panel, int)
        ! Calculates the F integrals necessary for determining the influence of a triangular panel in subsonic flow.
        ! This is a pared-down version of the algorithm presented by Johnson (1980) Appendix D.3.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(integrals),intent(inout) :: int

        integer :: i

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

        ! Calculate F(1,2,1) and F(2,1,1)
        if (source_order == 1 .or. doublet_order == 2) then

            ! Calculate (these formulas come from PAN AIR and are equivalent to Johnson, but simplified)
            int%F121 = geom%a*geom%v_eta*int%F111 + geom%v_xi*geom%dR
            int%F211 = geom%a*geom%v_xi*int%F111 - geom%v_eta*geom%dR

        end if

    end subroutine panel_calc_subsonic_edge_integrals


    subroutine panel_calc_supersonic_subinc_edge_integrals(this, geom, dod_info, freestream, mirror_panel, int)
        ! Calculates the F integrals necessary to determine the influence of a subinclined triangular panel in supersonic flow.
        ! Taken from Ehlers et al. (1979) Appendix E.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(integrals),intent(inout) :: int

        real :: F1, F2, eps, eps2, series, b, s_b
        integer :: i, i_next

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N) + 1

            ! Check DoD
            if (dod_info%edges_in_dod(i)) then

                ! Get b a and its square root
                if (mirror_panel) then
                    b = this%b_mir(i)
                    s_b = this%sqrt_b_mir(i)
                else
                    b = this%b(i)
                    s_b = this%sqrt_b(i)
                end if

                ! Calculate F factors
                if (b > 0.) then
                    F1 = (geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)) / geom%g2(i)
                    F2 = (b*geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)) / geom%g2(i)
                else
                    F1 = (geom%R2(i) - geom%R1(i))*(geom%R2(i) + geom%R1(i)) / (geom%l1(i)*geom%R2(i) + geom%l2(i)*geom%R1(i))
                    F2 = (geom%g2(i) - geom%l1(i)**2 - geom%l2(i)**2) / (b*geom%R1(i)*geom%R2(i) - geom%l1(i)*geom%l2(i))
                end if

                ! Calculate F(1,1,1) and other higher integrals if necessary

                ! Nearly-sonic edge
                if (abs(F2) > 100.0*abs(s_b*F1)) then

                    ! Calculate series solution
                    eps = F1/F2
                    eps2 = eps*eps
                    series = eps*eps2*(1./3. - b*eps2/5. + (b*eps2)*(b*eps2)/7.)
                    int%F111(i) = -eps + b*series

                    if (source_order == 1 .or. doublet_order == 2) then
                        if (mirror_panel) then
                            int%F121(i) = (-geom%v_xi(i)*geom%dR(i)*geom%R1(i)*geom%R2(i) &
                                           + geom%l2(i)*geom%R1(i)*(this%vertices_ls_mir(2,i_next) - geom%P_ls(2)) &
                                           - geom%l1(i)*geom%R2(i)*(this%vertices_ls_mir(2,i) - geom%P_ls(2))) / (geom%g2(i)*F2) &
                                          - geom%a(i)*geom%v_eta(i)*series
                        else
                            int%F121(i) = (-geom%v_xi(i)*geom%dR(i)*geom%R1(i)*geom%R2(i) &
                                           + geom%l2(i)*geom%R1(i)*(this%vertices_ls(2,i) - geom%P_ls(2)) &
                                           - geom%l1(i)*geom%R2(i)*(this%vertices_ls(2,i_next) - geom%P_ls(2))) / (geom%g2(i)*F2) &
                                          - geom%a(i)*geom%v_eta(i)*series
                        end if
                        int%F211(i) = -geom%v_eta(i)*geom%dR(i) + geom%a(i)*geom%v_xi(i)*int%F111(i) - &
                                      2.*geom%v_xi(i)*geom%v_eta(i)*int%F121(i)
                    end if

                ! Supersonic edge
                else if (b > 0.) then

                    ! Mach wedge
                    if (geom%R1(i) == 0. .and. geom%R2(i) == 0) then
                        int%F111(i) = pi/s_b

                        if (source_order == 1) then
                            int%F121(i) = -geom%a(i)*geom%v_eta(i)*int%F111(i)/b
                            int%F211(i) = geom%a(i)*geom%v_xi(i)*int%F111(i)/b
                        end if

                    ! At least one endpoint in
                    else
                        int%F111(i) = -atan2(s_b*F1, F2) / s_b

                        if (source_order == 1 .or. doublet_order == 2) then
                            int%F121(i) = -(geom%v_xi(i)*geom%dR(i) + geom%a(i)*geom%v_eta(i)*int%F111(i)) / b
                            int%F211(i) = -geom%v_eta(i)*geom%dR(i) + geom%a(i)*geom%v_xi(i)*int%F111(i) - &
                                          2.*geom%v_xi(i)*geom%v_eta(i)*int%F121(i)
                        end if
                    end if

                ! Subsonic edge
                else
                    F1 = s_b*geom%R1(i) + abs(geom%l1(i))
                    F2 = s_b*geom%R2(i) + abs(geom%l2(i))
                    int%F111(i) = -sign(1., geom%v_eta(i))*log(F1/F2)

                    if (source_order == 1 .or. doublet_order == 2) then
                        int%F121(i) = -(geom%v_xi(i)*geom%dR(i) + geom%a(i)*geom%v_eta(i)*int%F111(i)) / b
                        int%F211(i) = -geom%v_eta(i)*geom%dR(i) + geom%a(i)*geom%v_xi(i)*int%F111(i) - &
                                      2.*geom%v_xi(i)*geom%v_eta(i)*int%F121(i)
                    end if
                end if

            end if

        end do

    end subroutine panel_calc_supersonic_subinc_edge_integrals


    subroutine panel_calc_subsonic_panel_integrals(this, geom, freestream, mirror_panel, int)
        ! Calculates the necessary H integrals to determine the influence of a panel in subsonic flow.
        ! Taken from Johnson (1980) Appendix D.3. with alterations made based on PAN AIR.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(integrals),intent(inout) :: int

        real :: S, C, nu, c1, c2, x
        integer :: i

        ! Calculate hH(1,1,3) (Johnson Eqs. (D.41) and (G.24))
        ! No check on the magnitude of h is necessary since we never divide by it
        int%hH113 = 0.
        do i=1,this%N

            ! Calculate intermediate quantities
            c1 = geom%g2(i) + abs(geom%h)*geom%R1(i)
            c2 = geom%g2(i) + abs(geom%h)*geom%R2(i)
        
            ! Calculate integral for edge
            S = geom%a(i)*(geom%l2(i)*c1 - geom%l1(i)*c2)
            C = c1*c2 + geom%a(i)**2*geom%l1(i)*geom%l2(i)
            x = atan2(S, C)

            ! Sum
            int%hH113 = int%hH113 + x

        end do
        
        ! Apply sign factor (Johnson Eq. (D.42)
        int%hH113 = sign(int%hH113, geom%h)

        ! Calculate H(1,1,1)
        int%H111 = -geom%h*int%hH113 + sum(geom%a*int%F111)

        ! Calculate H(2,1,3) and H(1,2,3)
        int%H213 = -sum(geom%v_xi*int%F111)
        int%H123 = -sum(geom%v_eta*int%F111)

        ! Calculate higher-order source integrals
        if (source_order == 1) then
            int%H211 = 0.5*(-geom%h2*int%H213 + sum(geom%a*int%F211))
            int%H121 = 0.5*(-geom%h2*int%H123 + sum(geom%a*int%F121))
        end if

        ! Calculate higher-order doublet integrals
        if (doublet_order == 2) then
            int%H313 = -sum(geom%v_eta*int%F121) + geom%h*int%hH113
            int%H223 = sum(geom%v_xi*int%F121)
            int%H133 = -int%H111 + sum(geom%v_eta*int%F121)
        end if

    end subroutine panel_calc_subsonic_panel_integrals


    subroutine panel_calc_supersonic_subinc_panel_integrals(this, geom, dod_info, freestream, mirror_panel, int)
        ! Calculates the necessary H integrals to determine the influence of a subinclined panel in supersonic flow.
        ! Taken from Ehlers et al. (1979) Appendix E.

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(integrals),intent(inout) :: int

        real :: F1, F2, b, s_b
        integer :: i, i_next

        ! Calculate hH(1,1,3) (Ehlers Eq. (E18))
        int%hH113 = 0.

        ! Loop through edges
        do i=1,this%N

            ! Check DoD
            if (dod_info%edges_in_dod(i)) then

                ! Get b
                if (mirror_panel) then
                    b = this%b_mir(i)
                else
                    b = this%b(i)
                end if

                ! Check not on panel plane
                if (abs(geom%h) > 1.e-12) then

                    ! Calculate F factors
                    if (b > 0.) then
                        F1 = (geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)) / geom%g2(i)
                        F2 = (b*geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)) / geom%g2(i)
                    else
                        F1 = geom%dR(i)*(geom%R2(i) + geom%R1(i)) / (geom%l1(i)*geom%R2(i) + geom%l2(i)*geom%R1(i))
                        F2 = (geom%g2(i) - geom%l1(i)**2 - geom%l2(i)**2) &
                             / (b*geom%R1(i)*geom%R2(i) - geom%l1(i)*geom%l2(i))
                    end if

                    ! Supersonic edge
                    if (b > 0.) then

                        ! Mach wedge
                        if (geom%R1(i) == 0. .and. geom%R2(i) == 0.) then
                            int%hH113 = int%hH113 + pi*sign(1., geom%h*geom%v_xi(i))
                        else
                            int%hH113 = int%hH113 + atan2(geom%h*geom%a(i)*F1, geom%R1(i)*geom%R2(i) + geom%h2*F2)
                        end if

                    ! Subsonic edge
                    else
                        int%hH113 = int%hH113 + atan2(geom%h*geom%a(i)*F1, geom%R1(i)*geom%R2(i) + geom%h2*F2)
                    end if

                end if
            end if
        end do

        ! Calculate H(1,1,1)
        int%H111 = -geom%h*int%hH113 - sum(geom%a*int%F111)

        ! Calculate H(2,1,3) and H(1,2,3)
        int%H213 = -sum(geom%v_xi*int%F111)
        int%H123 = sum(geom%v_eta*int%F111)

        ! Calculate higher-order source integrals
        if (source_order == 1) then
            int%H211 = 0.5*(-geom%h2*int%H213 - sum(geom%a*int%F211))
            int%H121 = 0.5*(-geom%h2*int%H123 - sum(geom%a*int%F121))
        end if

        ! Calculate higher-order doublet integrals
        if (doublet_order == 2) then
            int%H313 = -sum(geom%v_eta*int%F121) - geom%h*int%hH113
            int%H223 = sum(geom%v_xi*int%F121)
            int%H133 = -sum(geom%v_eta*int%F121) - int%H111
        end if

    end subroutine panel_calc_supersonic_subinc_panel_integrals


    function panel_calc_integrals(this, geom, influence_type, singularity_type, freestream, mirror_panel, dod_info) result(int)
        ! Calculates the H and F integrals necessary for the given influence

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        character(len=*),intent(in) :: influence_type, singularity_type
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(dod),intent(in) :: dod_info

        type(integrals) :: int

        ! Allocate space for integrals
        allocate(int%F111(this%N), source=0.)
        if (source_order == 1 .or. doublet_order == 2) then
            allocate(int%F121(this%N))
            allocate(int%F211(this%N))
        end if

        ! Calculate necessary integrals based on the flow condition and panel type
        if (freestream%supersonic) then
            call this%calc_supersonic_subinc_edge_integrals(geom, dod_info, freestream, mirror_panel, int)
            call this%calc_supersonic_subinc_panel_integrals(geom, dod_info, freestream, mirror_panel, int)
        else
            call this%calc_subsonic_edge_integrals(geom, freestream, mirror_panel, int)
            call this%calc_subsonic_panel_integrals(geom, freestream, mirror_panel, int)
        end if

    end function panel_calc_integrals


    subroutine panel_allocate_potential_influences(this, phi_s, phi_d)
        ! Allocates the necessary influence arrays

        implicit none

        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(out) :: phi_s, phi_d

        ! Source
        if (source_order == 0) then
            allocate(phi_s(1), source=0.)
        else if (source_order == 1) then
            allocate(phi_s(3), source=0.)
        end if

        ! Doublet
        if (doublet_order == 1) then

            ! Check if this panel belongs to the wake
            if (this%in_wake) then
                allocate(phi_d(6), source=0.)
            else
                allocate(phi_d(3), source=0.)
            end if

        else if (doublet_order == 2) then

            ! Check if this panel belongs to the wake
            if (this%in_wake) then
                allocate(phi_d(12), source=0.)
            else
                allocate(phi_d(6), source=0.)
            end if

        end if
        
    end subroutine panel_allocate_potential_influences


    subroutine panel_calc_potentials(this, P, freestream, dod_info, mirror_panel, phi_s, phi_d)
        ! Calculates the source- and doublet-induced potentials at the given point P

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: P
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        logical,intent(in) :: mirror_panel
        real,dimension(:),allocatable,intent(out) :: phi_s, phi_d

        type(eval_point_geom) :: geom
        type(integrals) :: int
        integer :: i

        ! Specify influencing vertices (also sets zero default influence)
        call this%allocate_potential_influences(phi_s, phi_d)

        ! Check DoD
        if (dod_info%in_dod .and. this%A > 0.) then

            ! Calculate geometric parameters
            if (freestream%supersonic) then
                geom = this%calc_supersonic_subinc_geom(P, freestream, mirror_panel, dod_info)
            else
                geom = this%calc_subsonic_geom(P, freestream, mirror_panel)
            end if

            ! Get integrals
            int = this%calc_integrals(geom, 'potential', 'doublet', freestream, mirror_panel, dod_info)

            ! Source potential
            phi_s(1) = int%H111
            if (source_order == 1) then

                ! Johnson Eq. (D21)
                ! Equivalent to Ehlers Eq. (8.6)
                phi_s(2) = int%H111*geom%P_ls(1) + int%H211
                phi_s(3) = int%H111*geom%P_ls(2) + int%H121

                ! Convert to vertex influences (Davis Eq. (4.41))
                if (mirror_panel) then
                    phi_s = matmul(phi_s, this%S_sigma_inv_mir)
                else
                    phi_s = matmul(phi_s, this%S_sigma_inv)
                end if

            end if
            
            ! Add area Jacobian and kappa factor
            phi_s = -this%J*freestream%K_inv*phi_s

            ! Doublet potential
            ! Johnson Eq. (D.30)
            ! Equivalent to Ehlers Eq. (5.17))
            phi_d(1) = int%hH113
            phi_d(2) = int%hH113*geom%P_ls(1) + geom%h*int%H213
            phi_d(3) = int%hH113*geom%P_ls(2) + geom%h*int%H123

            if (doublet_order == 1) then

                ! Convert to vertex influences (Davis Eq. (4.41))
                if (mirror_panel) then
                    phi_d(1:3) = freestream%K_inv*matmul(phi_d(1:3), this%S_mu_inv_mir)
                else
                    phi_d(1:3) = freestream%K_inv*matmul(phi_d(1:3), this%S_mu_inv)
                end if

                ! Wake bottom influence is opposite the top influence
                if (this%in_wake) then
                    phi_d(4:6) = -phi_d(1:3)
                end if

            else if (doublet_order == 2) then

                ! Add quadratic terms
                phi_d(4) = 0.5*int%hH113*geom%P_ls(1)**2 + geom%h*(geom%P_ls(1)*int%H213 + 0.5*int%H313)

                phi_d(5) = int%hH113*geom%P_ls(1)*geom%P_ls(2) + geom%h*(geom%P_ls(2)*int%H213 + geom%P_ls(1)*int%H123 + int%H223)

                phi_d(6) = 0.5*int%hH113*geom%P_ls(2)**2 + geom%h*(geom%P_ls(2)*int%H123 + 0.5*int%H133)

                ! Convert to vertex influences (Davis Eq. (4.41))
                if (mirror_panel) then
                    phi_d(1:6) = freestream%K_inv*matmul(phi_d(1:6), this%S_mu_inv_mir)
                else
                    phi_d(1:6) = freestream%K_inv*matmul(phi_d(1:6), this%S_mu_inv)
                end if

                ! Wake bottom influence is opposite the top influence
                if (this%in_wake) then
                    phi_d(7:12) = -phi_d(1:6)
                end if

            end if
        end if
    
    end subroutine panel_calc_potentials


    subroutine panel_calc_velocities(this, P, freestream, dod_info, mirror_panel, v_s, v_d)
        ! Calculates the source- and doublet-induced potentials at the given point P

        implicit none

        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: P
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        logical,intent(in) :: mirror_panel
        real,dimension(:,:),allocatable,intent(out) :: v_s, v_d

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

                ! Set default influence
                allocate(v_d(6,3), source=0.)

            else

                ! Set default influence
                allocate(v_d(3,3), source=0.)

            end if
        end if

        ! Check DoD
        if (dod_info%in_dod .and. this%A > 0.) then

            ! Calculate geometric parameters
            if (freestream%supersonic) then
                geom = this%calc_supersonic_subinc_geom(P, freestream, mirror_panel, dod_info)
            else
                geom = this%calc_subsonic_geom(P, freestream, mirror_panel)
            end if

            ! Get integrals
            int = this%calc_integrals(geom, "velocity", "doublet", freestream, mirror_panel, dod_info)

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

        real,dimension(:),allocatable :: mu_verts, mu_params
        integer :: i
        real :: s

        ! Jump calculations for mirrored panel
        if (mirrored) then

            ! Set up array of doublet strengths to calculate doublet parameters
            if (doublet_order == 1) then

                ! Allocate
                allocate(mu_verts(this%N))

                ! Get doublet values
                do i=1,this%N
                    mu_verts(i) = mu(this%get_vertex_index(i)+size(mu)/2)
                end do

            else if (doublet_order == 2) then

                ! Allocate
                allocate(mu_verts(2*this%N))

                ! Get doublet values
                do i=1,this%N
                    mu_verts(i) = mu(this%get_vertex_index(i)+size(mu)/2)
                    mu_verts(i+this%N) = mu(this%get_midpoint_index(i)+size(mu)/2)
                end do

            end if
        
            ! Calculate doublet parameters (derivatives)
            mu_params = matmul(this%S_mu_inv_mir, mu_verts)

            ! Calculate tangential velocity jump in panel coordinates E&M Eq. (N.1.11b)
            dv(1) = mu_params(2)
            dv(2) = mu_params(3)
            dv(3) = 0.

            ! Transform to global coordinates
            dv = matmul(transpose(this%A_g_to_ls_mir), dv)

            ! Get source strength
            
            ! Constant
            if (source_order == 0) then
                s = sigma(this%index+size(sigma)/2)

            ! For linear distribution, use the average
            else if (source_order == 1) then
                s = 0.
                do i=1,this%N
                    s = s + sigma(this%get_vertex_index(i)+size(sigma)/2)
                end do
                s = s/this%N
            end if

            ! Add normal velocity jump in global coords E&M Eq. (N.1.11b)
            dv = dv + s*this%n_g_mir/inner(this%nu_g_mir, this%n_g_mir)

        ! Jump calculations for original panel
        ! Same steps as above
        else

            if (doublet_order == 1) then

                allocate(mu_verts(this%N))

                do i=1,this%N
                    mu_verts(i) = mu(this%get_vertex_index(i))
                end do

            else if (doublet_order == 2) then

                allocate(mu_verts(2*this%N))

                do i=1,this%N
                    mu_verts(i) = mu(this%get_vertex_index(i))
                    mu_verts(i+this%N) = mu(this%get_midpoint_index(i))
                end do

            end if

            do i=1,this%N
                mu_verts(i) = mu(this%get_vertex_index(i))
            end do
        
            mu_params = matmul(this%S_mu_inv, mu_verts)

            dv(1) = mu_params(2)
            dv(2) = mu_params(3)
            dv(3) = 0.

            dv = matmul(transpose(this%A_g_to_ls), dv)

            if (source_order == 0) then
                s = sigma(this%index)

            else if (source_order == 1) then
                s = 0.
                do i=1,this%N
                    s = s + sigma(this%get_vertex_index(i))
                end do
                s = s/this%N
            end if

            dv = dv + s*this%n_g/inner(this%nu_g, this%n_g)

        end if

    end function panel_get_velocity_jump

    
end module panel_mod