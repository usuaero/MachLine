! Panel type
module panel_mod

    use helpers_mod
    use linked_list_mod
    use base_geom_mod
    use math_mod
    use flow_mod
    use linalg_mod

    implicit none

    integer :: doublet_order
    integer :: source_order


    type integrals
        ! Container type for the fundamental integrals used to calculate influence coefficients

        real :: H111, H211, H121 ! Source integrals
        real :: hH113, H213, H123, H313, H223, H133 ! Doublet integrals; we use hH(1,1,3) because it can be reliably calculated, unlike H(1,1,3)
        real,dimension(:),allocatable :: F111, F211, F121 ! Necessary line integrals

    end type integrals


    type dod
        ! Container type for parameters of whether a panel lies in a point's domain of dependence

        logical :: in_dod = .true.
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
        real,dimension(:,:),allocatable :: n_hat_g, n_hat_ls ! Edge unit outward normals
        real,dimension(:,:),allocatable :: n_hat_g_mir, n_hat_ls_mir
        real,dimension(:),allocatable :: b, sqrt_b ! Edge parameter
        real,dimension(:),allocatable :: b_mir, sqrt_b_mir
        real :: A ! Surface area (same for mirror, in global coordinates at least)
        real,dimension(:,:),allocatable :: S_mu_inv, S_sigma_inv ! Matrix relating doublet/source strengths to doublet/source influence parameters
        real,dimension(:,:),allocatable :: S_mu_inv_mir, S_sigma_inv_mir
        logical :: in_wake ! Whether this panel belongs to a wake mesh
        integer,dimension(3) :: abutting_panels ! Indices of panels abutting this one
        integer,dimension(3) :: edges ! Indices of the edges of this panel
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
            procedure :: calc_g_edge_vectors => panel_calc_g_edge_vectors

            ! Flow-dependent initialization procedures
            procedure :: init_with_flow => panel_init_with_flow
            procedure :: calc_g_to_ls_transform => panel_calc_g_to_ls_transform
            procedure :: calc_ls_edge_vectors => panel_calc_ls_edge_vectors
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
            procedure :: get_subpanel_centroid => panel_get_subpanel_centroid
            procedure :: get_corner_angle => panel_get_corner_angle
            procedure :: get_weighted_normal_at_corner => panel_get_weighted_normal_at_corner
            procedure :: get_projection_onto_surface => panel_get_projection_onto_surface

            ! Checks
            procedure :: touches_vertex => panel_touches_vertex
            procedure :: check_abutting_mirror_plane => panel_check_abutting_mirror_plane
            procedure :: projection_inside => panel_projection_inside
            procedure :: point_outside => panel_point_outside

            ! Update information
            procedure :: point_to_new_vertex => panel_point_to_new_vertex

            ! Domain of dependence checking
            procedure :: check_dod => panel_check_dod

            ! Influence calculations

            ! Geometry
            procedure :: calc_subsonic_geom => panel_calc_subsonic_geom
            procedure :: calc_supersonic_subinc_geom => panel_calc_supersonic_subinc_geom
            procedure :: calc_supersonic_supinc_geom => panel_calc_supersonic_supinc_geom

            ! Fundamental integrals
            procedure :: calc_subsonic_edge_integrals => panel_calc_subsonic_edge_integrals
            procedure :: calc_subsonic_panel_integrals => panel_calc_subsonic_panel_integrals
            procedure :: calc_supersonic_subinc_edge_integrals => panel_calc_supersonic_subinc_edge_integrals
            procedure :: calc_supersonic_subinc_panel_integrals => panel_calc_supersonic_subinc_panel_integrals
            procedure :: calc_supersonic_supinc_edge_integrals => panel_calc_supersonic_supinc_edge_integrals
            procedure :: calc_supersonic_supinc_panel_integrals => panel_calc_supersonic_supinc_panel_integrals
            procedure :: calc_integrals => panel_calc_integrals

            ! Influences
            procedure :: allocate_potential_influences => panel_allocate_potential_influences
            procedure :: calc_potential_influences => panel_calc_potential_influences
            procedure :: calc_potentials => panel_calc_potentials

            ! Results
            procedure :: get_source_strengths => panel_get_source_strengths
            procedure :: get_doublet_strengths => panel_get_doublet_strengths
            procedure :: get_source_dist_parameters => panel_get_source_dist_parameters
            procedure :: get_doublet_dist_parameters => panel_get_doublet_dist_parameters
            procedure :: get_velocity_jump => panel_get_velocity_jump

    end type panel

    
contains


    subroutine panel_init_3(this, v1, v2, v3, index)
        ! Initializes a 3-panel

        implicit none

        class(panel),intent(inout) :: this
        type(vertex),intent(inout),target :: v1, v2, v3
        integer,intent(in) :: index

        integer :: i

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

        ! Store that this panel is attached to its vertices
        do i=1,this%N
            call this%vertices(i)%ptr%panels%append(this%index)
            call this%vertices(i)%ptr%panels_not_across_wake_edge%append(this%index)
        end do

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

        ! Calculate ledge vectors
        call this%calc_g_edge_vectors()

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


    subroutine panel_calc_g_edge_vectors(this)

        implicit none

        class(panel),intent(inout) :: this

        real,dimension(3) :: d_g, t_hat_g
        integer :: i, i_next

        ! Allocate memory
        allocate(this%n_hat_g(3,this%N))

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N)+1

            ! Calculate edge vector based on index
            d_g = this%get_vertex_loc(i_next)-this%get_vertex_loc(i)

            ! Calculate tangent in global coords
            t_hat_g = d_g/norm2(d_g)

            ! Calculate edge outward normal
            this%n_hat_g(:,i) = cross(t_hat_g, this%n_g)

        end do
    
    end subroutine panel_calc_g_edge_vectors


    subroutine panel_init_with_flow(this, freestream, initialize_mirror, mirror_plane)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream
        logical,intent(in) :: initialize_mirror
        integer,intent(in) :: mirror_plane

        ! Calculate transforms
        call this%calc_g_to_ls_transform(freestream)

        ! Calculate properties dependent on the transforms
        call this%calc_ls_edge_vectors(freestream)
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
        if (abs(abs(inner(this%n_g, freestream%c_hat_g)) - 1.) < 1e-12) then ! Check the freestream isn't aligned with the normal vector
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
        if (freestream%supersonic .and. abs(x) < 1.e-12) then
            write(*,*) "!!! Panel", this%index, "is Mach-inclined, which is not allowed. Quitting..."
            stop
        end if

        ! Calculate panel inclination indicator (E&M Eq. (E.3.16b))
        this%r = sign(1., x) ! r = -1 -> superinclined, r = 1 -> subinclined

        ! Other inclination parameters
        rs = this%r*freestream%s

        ! Calculate transformation
        y = 1./sqrt(abs(x))
        this%A_g_to_ls(1,:) = y*matmul(freestream%C_mat_g, u0)
        this%A_g_to_ls(2,:) = rs/freestream%B*matmul(freestream%C_mat_g, v0)
        this%A_g_to_ls(3,:) = freestream%B*y*this%n_g

        ! Check determinant
        x = det3(this%A_g_to_ls)
        if (abs(x - freestream%B*freestream%B) > 1.e-12) then
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


    subroutine panel_calc_ls_edge_vectors(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(2) :: d_ls
        real,dimension(:,:),allocatable :: t_hat_ls
        integer :: i, i_next

        ! Allocate memory
        allocate(t_hat_ls(2,this%N))
        allocate(this%n_hat_ls(2,this%N))
        allocate(this%b(this%N))
        allocate(this%b_mir(this%N)) ! This needs to be initialized here because of some DoD checks. It will have no effect.
        allocate(this%sqrt_b(this%N))

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N)+1

            ! Calculate tangent in local scaled coords 
            d_ls = this%vertices_ls(:,i_next) - this%vertices_ls(:,i)
            t_hat_ls(:,i) = d_ls/norm2(d_ls)

        end do

        ! Calculate edge normal in local scaled coords E&M Eq. (J.6.45)
        this%n_hat_ls(1,:) = t_hat_ls(2,:)
        this%n_hat_ls(2,:) = -t_hat_ls(1,:)

        ! Calculate edge parameter (Ehlers Eq. (E14))
        ! This really only matters for subinclined, supersonic panels
        ! But we set defaults for the other cases to make unified calcs work
        if (freestream%supersonic) then
            if (this%r > 0) then
                this%b = (this%n_hat_ls(1,:) - this%n_hat_ls(2,:))*(this%n_hat_ls(1,:) + this%n_hat_ls(2,:))
                this%sqrt_b = sqrt(abs(this%b))
            else
                this%b = 1.
                this%sqrt_b = 1.
            end if
        else
            this%b = -1.
            this%sqrt_b = 1.
        end if
    
    end subroutine panel_calc_ls_edge_vectors


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

            ! x
            S_mu(1:3,2) = this%vertices_ls(1,:)
            S_mu(4:6,2) = this%midpoints_ls(1,:)

            ! y
            S_mu(1:3,3) = this%vertices_ls(2,:)
            S_mu(4:6,3) = this%midpoints_ls(2,:)

            ! x^2
            S_mu(:,4) = S_mu(:,2)**2*0.5

            ! xy
            S_mu(:,5) = S_mu(:,2)*S_mu(:,3)

            ! y^2
            S_mu(:,6) = S_mu(:,3)**2*0.5

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
        call this%calc_mirrored_edge_vectors(freestream)

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
        if (abs(abs(inner(this%n_g_mir, freestream%c_hat_g)) - 1.) < 1e-12) then ! Check the freestream isn't aligned with the normal vector
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


    subroutine panel_calc_mirrored_edge_vectors(this, freestream)

        implicit none

        class(panel),intent(inout) :: this
        type(flow),intent(in) :: freestream

        real,dimension(2) :: d_ls
        real,dimension(:,:),allocatable :: t_hat_ls_mir
        integer :: i, i_next

        ! Allocate memory
        allocate(t_hat_ls_mir(2,this%N))
        allocate(this%n_hat_ls_mir(2,this%N))
        allocate(this%sqrt_b_mir(this%N))

        ! Loop through edges
        do i=1,this%N

            i_next = mod(i, this%N)+1

            ! Calculate tangent in local scaled coords 
            ! Direction is flipped so that we're still going counter-clockwise about the panel
            d_ls = this%vertices_ls_mir(:,i) - this%vertices_ls_mir(:,i_next)
            t_hat_ls_mir(:,i) = d_ls/norm2(d_ls)

        end do

        ! Calculate edge normal in local scaled coords E&M Eq. (J.6.45)
        this%n_hat_ls_mir(1,:) = t_hat_ls_mir(2,:)
        this%n_hat_ls_mir(2,:) = -t_hat_ls_mir(1,:)

        ! Calculate edge parameter (Ehlers Eq. (E14))
        if (freestream%supersonic) then
            if (this%r_mir > 0) then
                this%b_mir = (this%n_hat_ls_mir(1,:) - this%n_hat_ls_mir(2,:))*(this%n_hat_ls_mir(1,:) + this%n_hat_ls_mir(2,:))
                this%sqrt_b_mir = sqrt(abs(this%b_mir))
            else
                this%b_mir = 1.
                this%sqrt_b_mir = 1.
            end if
        else
            this%b_mir = -1.
            this%sqrt_b_mir = 1.
        end if
    
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

            ! x
            S_mu(1:3,2) = this%vertices_ls_mir(1,:)
            S_mu(4:6,2) = this%midpoints_ls_mir(1,:)

            ! y
            S_mu(1:3,3) = this%vertices_ls_mir(2,:)
            S_mu(4:6,3) = this%midpoints_ls_mir(2,:)

            ! x^2
            S_mu(:,4) = 0.5*S_mu(:,2)**2

            ! xy
            S_mu(:,5) = S_mu(:,2)*S_mu(:,3)

            ! y^2
            S_mu(:,6) = 0.5*S_mu(:,3)**2

            ! Invert
            call matinv(6, S_mu, this%S_mu_inv_mir)

        end if

        ! Determine influence of vertex source strengths on integral parameters
        ! Linear distribution
        if (source_order == 1) then

            ! Allocate influence matrices
            allocate(S_sigma(3,3))
            allocate(this%S_sigma_inv_mir(3,3))

            ! Set values
            S_sigma(:,1) = 1.
            S_sigma(:,2) = this%vertices_ls_mir(1,:)
            S_sigma(:,3) = this%vertices_ls_mir(2,:)

            ! Invert
            call matinv(3, S_sigma, this%S_sigma_inv_mir)

        ! Quadratic distribution; I don't anticipate this will ever be used
        else if (source_order == 2) then

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


    function panel_check_abutting_mirror_plane(this, N_panels, i_endpoints, edge_index) result(abuts)
        ! Tells whether this panel abuts the mirror plane.

        class(panel),intent(inout) :: this
        integer,intent(in) :: N_panels
        integer,dimension(2),intent(out) :: i_endpoints
        integer,intent(out) :: edge_index

        logical :: abuts

        logical :: already_found_vert_on_mirror_plane
        integer :: m, m1, temp

        ! Initialize checks
        already_found_vert_on_mirror_plane = .false.
        abuts = .false.

        ! Loop through vertices
        mirror_loop: do m=1,this%N

            ! Check if vertex is on the mirror plane
            if (this%vertices(m)%ptr%on_mirror_plane) then

                ! Previously found a vertex on mirror plane, so the panels are abutting
                if (already_found_vert_on_mirror_plane) then

                    abuts = .true.

                    ! Store the second shared vertex
                    i_endpoints(2) = this%get_vertex_index(m)

                    ! Check order
                    if (m1 == 1 .and. m == 3) then
                        temp = i_endpoints(1)
                        i_endpoints(1) = i_endpoints(2)
                        i_endpoints(2) = temp
                    end if

                    ! Store adjacent panel
                    if (m-m1 == 1) then
                        this%abutting_panels(m1) = this%index + N_panels
                        edge_index = m1
                    else
                        this%abutting_panels(m) = this%index + N_panels
                        edge_index = m
                    end if

                    return

                ! First vertex on the mirror plane
                else

                    already_found_vert_on_mirror_plane = .true.
                    i_endpoints(1) = this%get_vertex_index(m)
                    m1 = m

                end if
            end if

        end do mirror_loop
        
    end function panel_check_abutting_mirror_plane


    function panel_projection_inside(this, point, mirror_panel, mirror_plane) result(inside)
        ! Checks whether the given point, when projected into the plane of the panel, is inside the panel

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: point
        logical,intent(in) :: mirror_panel
        integer,intent(in) :: mirror_plane

        logical :: inside

        real,dimension(3) :: d
        integer :: i
        real :: x

        ! Loop through edges
        inside = .true.
        do i=1,this%N

            ! Shift origin to the edge
            if (mirror_panel) then
                d = point - mirror_across_plane(this%get_vertex_loc(i), mirror_plane)
            else
                d = point - this%get_vertex_loc(i)
            end if
            
            ! Get dot product with outer normal
            if (mirror_panel) then
                x = inner(d, this%n_hat_g_mir)
            else
                x = inner(d, this%n_hat_g)
            end if

            ! Check
            if (x > 0.) then
                inside = .false.
                return
            end if

        end do
        
    end function panel_projection_inside


    function panel_point_outside(this, point, mirror_panel, mirror_plane) result(outside)
        ! Tells whether the given point is above the panel and its projection is inside the surface of the panel

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: point
        logical,intent(in) :: mirror_panel
        integer,intent(in) :: mirror_plane

        logical :: outside

        real :: h

        ! Get height above panel
        if (mirror_panel) then
            h = inner(point-this%centr_mir, this%n_g_mir)
        else
            h = inner(point-this%centr, this%n_g)
        end if

        ! If height is negative, we know this isn't outside the panel
        if (h < 0.) then
            outside = .false.

        ! Otherwise, it's dependent upon whether the projection is inside the panel surface
        else
            outside = this%projection_inside(point, mirror_panel, mirror_plane)
        end if
        
    end function panel_point_outside


    function panel_get_subpanel_centroid(this, j) result(cent)
        ! Calculates the centroid of the subpanel on edges j and j+1

        implicit none
        
        class(panel),intent(in) :: this
        integer,intent(in) :: j

        real,dimension(3) :: cent

        cent = ( this%get_midpoint_loc(j) &
               + this%get_vertex_loc(modulo(j, this%N)+1) & 
               + this%get_midpoint_loc(modulo(j, this%N)+1) ) / 3.0
        
    end function panel_get_subpanel_centroid


    function panel_get_corner_angle(this, vert_loc) result(angle)
        ! Calculates the angle of the corner at which the given vertex lies

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: vert_loc

        real :: angle

        integer :: i, i_prev

        ! Find the right corner
        do i=1,this%N

            ! Check vertex
            if (dist(this%get_vertex_loc(i), vert_loc) < 1.e-12) then

                ! Get previous edge index
                if (i == 1) then
                    i_prev = this%N
                else
                    i_prev = i-1
                end if

                ! Calculate angle
                angle = acos(inner(-this%n_hat_g(:,i), this%n_hat_g(:,i_prev)))
                return

            end if

            ! Check midpoint
            if (doublet_order == 2) then
                if (dist(this%get_midpoint_loc(i), vert_loc) < 1.e-12) then
                    angle = pi
                    return
                end if
            end if
        end do

        ! This vertex doesn't belong, so it has no angle
        angle = 0.
        
    end function panel_get_corner_angle


    function panel_get_weighted_normal_at_corner(this, vert_loc) result(n_weighted)
        ! Returns the panel normal weighted by the angle of the corner given

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: vert_loc

        real,dimension(3) :: n_weighted

        real :: W

        ! Get angle
        W = this%get_corner_angle(vert_loc)
    
        ! Apply weight
        n_weighted = this%n_g*W

    end function panel_get_weighted_normal_at_corner


    function panel_get_projection_onto_surface(this, v, mirror_panel) result(v_proj)
        ! Projects the given vector into the plane of the panel

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: v
        logical,intent(in) :: mirror_panel

        real,dimension(3) :: v_proj
    
        ! Project
        if (mirror_panel) then
            v_proj = v - this%n_g_mir*inner(v, this%n_g_mir)
        else
            v_proj = v - this%n_g*inner(v, this%n_g)
        end if
        
    end function panel_get_projection_onto_surface


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

        class(panel),intent(in) :: this
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
        logical,dimension(3) :: these_verts_in_dod

        ! Set default mirroring
        if (present(mirror_panel)) then
            mirrored = mirror_panel
        else
            mirrored = .false.
        end if

        ! First check the flow is supersonic
        if (freestream%supersonic) then

            ! Read in vertex information
            do i=1,this%N
                if (mirrored) then
                    these_verts_in_dod(i) = verts_in_dod(this%get_vertex_index(i)+size(verts_in_dod)/2)
                else
                    these_verts_in_dod(i) = verts_in_dod(this%get_vertex_index(i))
                end if
            end do

            ! If all the vertices are in, then the panel is totally in and we can be done
            if (all(these_verts_in_dod)) then
                dod_info%in_dod = .true.
                dod_info%edges_in_dod = .true.

            ! If it is not guaranteed to be totally in, then check all the edges
            else

                ! Check edges
                do i=1,this%N
                    
                    i_next = mod(i, this%N)+1

                    ! If if at least one endpoint is in, then the edge is in (you gotta love convex subspaces)
                    if (these_verts_in_dod(i) .or. these_verts_in_dod(i_next)) then
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
                if (any(these_verts_in_dod) .or. any(dod_info%edges_in_dod)) then
                    dod_info%in_dod = .true.

                ! If a superinclined panel has no edges or vertices in the DoD, check if the DoD is encompassed by the panel
                else if (this%r == -1) then

                    ! Get the projection of the evaluation point onto the panel in the direction of c_hat
                    if (mirrored) then
                        s_star = inner(mirror_across_plane(this%get_vertex_loc(1), mirror_plane) - eval_point, this%n_g_mir) &
                                 / inner(freestream%c_hat_g, this%n_g_mir)
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

                ! Not superinclined and no edges or vertices in. Not in.
                else
                    dod_info%in_dod = .false.

                end if
            end if

        else

            ! Subsonic flow. DoD is everywhere. Life is easy.
            ! This shouldn't be necessary, but I'll keep it here for now.
            dod_info%in_dod = .true.
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
            geom%v_xi = this%n_hat_ls_mir(1,:)
            geom%v_eta = this%n_hat_ls_mir(2,:)
        else
            call geom%init(eval_point, this%A_g_to_ls, this%centr)
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
        real :: dummy

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

                ! Swap directions for mirror
                if (mirror_panel) then

                    ! Swap l1 and l2
                    ! The check is necessary because we set the sign of l1 and l2 in the case of R=0 based on which end each came from
                    dummy = geom%l1(i)
                    if (geom%R2(i) == 0.) then
                        geom%l1(i) = -geom%l2(i)
                    else
                        geom%l1(i) = geom%l2(i)
                    end if
                    if (geom%R1(i) == 0) then
                        geom%l2(i) = -dummy
                    else
                        geom%l2(i) = dummy
                    end if

                    ! Swap R1 and R2
                    dummy = geom%R1(i)
                    geom%R1(i) = geom%R2(i)
                    geom%R2(i) = dummy

                end if

            end if

        end do

        ! Difference in R
        geom%dR = geom%R2 - geom%R1

    end function panel_calc_supersonic_subinc_geom


    function panel_calc_supersonic_supinc_geom(this, eval_point, freestream, mirror_panel, dod_info) result(geom)
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
        real :: dummy

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

                ! Perpendicular in-plane distance
                geom%a(i) = geom%v_xi(i)*d_ls(1) + geom%v_eta(i)*d_ls(2)

                ! Edge integration length
                geom%l1(i) = -geom%v_eta(i)*d_ls(1) + geom%v_xi(i)*d_ls(2)

                ! Hyperbolic radius to first vertex
                x = -d_ls(1)*d_ls(1) - d_ls(2)*d_ls(2) + geom%h2
                if (geom%h > 0. .and. x > 0.) then
                    geom%R1(i) = sqrt(x)
                else
                    geom%l1(i) = -1.
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
                geom%l2(i) = -geom%v_eta(i)*d_ls(1) + geom%v_xi(i)*d_ls(2)

                ! Hyperbolic radius to first vertex
                x = -d_ls(1)*d_ls(1) - d_ls(2)*d_ls(2) + geom%h2
                if (geom%h > 0. .and. x > 0.) then
                    geom%R2(i) = sqrt(x)
                else
                    geom%l2(i) = 1.
                    geom%R2(i) = 0.
                end if

                ! Swap directions for mirror
                if (mirror_panel) then

                    ! Swap l1 and l2
                    ! The check is necessary because we set the sign of l1 and l2 in the case of R=0 based on which end each came from
                    dummy = geom%l1(i)
                    if (geom%R2(i) == 0.) then
                        geom%l1(i) = -geom%l2(i)
                    else
                        geom%l1(i) = geom%l2(i)
                    end if
                    if (geom%R1(i) == 0) then
                        geom%l2(i) = -dummy
                    else
                        geom%l2(i) = dummy
                    end if

                    ! Swap R1 and R2
                    dummy = geom%R1(i)
                    geom%R1(i) = geom%R2(i)
                    geom%R2(i) = dummy

                end if

            end if

        end do

        ! Difference in R
        geom%dR = geom%R2 - geom%R1

    end function panel_calc_supersonic_supinc_geom


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
        logical :: higher_order

        ! Check order
        higher_order = source_order == 1 .or. doublet_order == 2

        ! Loop through edges
        do i=1,this%N

            ! Check DoD
            if (dod_info%edges_in_dod(i)) then

                i_next = mod(i, this%N) + 1

                ! Get b and its square root; doing this removes a lot of mirror checks later
                if (mirror_panel) then
                    b = this%b_mir(i)
                    s_b = this%sqrt_b_mir(i)
                else
                    b = this%b(i)
                    s_b = this%sqrt_b(i)
                end if

                ! Mach wedge
                if (geom%R1(i) == 0. .and. geom%R2(i) == 0) then

                    ! F(1,1,1)
                    int%F111(i) = pi/s_b

                    ! Higher-order
                    if (higher_order) then
                        int%F121(i) = -geom%a(i)*geom%v_eta(i)*int%F111(i)/b
                        int%F211(i) = geom%a(i)*geom%v_xi(i)*int%F111(i)/b
                    end if

                else

                    ! Calculate F factors
                    if (b > 0.) then
                        F1 = (geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)) / geom%g2(i)
                        F2 = (b*geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)) / geom%g2(i)
                    else
                        F1 = (geom%R2(i) - geom%R1(i))*(geom%R2(i) + geom%R1(i)) / (geom%l1(i)*geom%R2(i) + geom%l2(i)*geom%R1(i))
                        F2 = (geom%g2(i) - geom%l1(i)**2 - geom%l2(i)**2) / (b*geom%R1(i)*geom%R2(i) - geom%l1(i)*geom%l2(i))
                    end if

                    ! Nearly-sonic edge
                    if (abs(F2) > 100.0*abs(s_b*F1)) then

                        ! F(1,1,1)
                        eps = F1/F2
                        eps2 = eps*eps
                        series = eps*eps2*(1./3. - b*eps2/5. + (b*eps2)*(b*eps2)/7.)
                        int%F111(i) = -eps + b*series

                        ! Higher-order
                        if (higher_order) then
                            if (mirror_panel) then
                                int%F121(i) = (-geom%v_xi(i)*geom%dR(i)*geom%R1(i)*geom%R2(i) &
                                               + geom%l2(i)*geom%R1(i)*(this%vertices_ls_mir(2,i_next) - geom%P_ls(2)) &
                                               - geom%l1(i)*geom%R2(i)*(this%vertices_ls_mir(2,i) - geom%P_ls(2)) &
                                              ) / (geom%g2(i)*F2) - geom%a(i)*geom%v_eta(i)*series
                            else
                                int%F121(i) = (-geom%v_xi(i)*geom%dR(i)*geom%R1(i)*geom%R2(i) &
                                               + geom%l2(i)*geom%R1(i)*(this%vertices_ls(2,i) - geom%P_ls(2)) &
                                               - geom%l1(i)*geom%R2(i)*(this%vertices_ls(2,i_next) - geom%P_ls(2)) &
                                              ) / (geom%g2(i)*F2) - geom%a(i)*geom%v_eta(i)*series
                            end if
                            int%F211(i) = -geom%v_eta(i)*geom%dR(i) + geom%a(i)*geom%v_xi(i)*int%F111(i) - &
                                          2.*geom%v_xi(i)*geom%v_eta(i)*int%F121(i)
                        end if

                    ! Supersonic edge
                    else if (b > 0.) then

                        ! F(1,1,1)
                        int%F111(i) = -atan2(s_b*F1, F2) / s_b

                        ! Higher-order
                        if (higher_order) then
                            int%F121(i) = -(geom%v_xi(i)*geom%dR(i) + geom%a(i)*geom%v_eta(i)*int%F111(i)) / b
                            int%F211(i) = -geom%v_eta(i)*geom%dR(i) + geom%a(i)*geom%v_xi(i)*int%F111(i) - &
                                          2.*geom%v_xi(i)*geom%v_eta(i)*int%F121(i)
                        end if

                    ! Subsonic edge
                    else
                        
                        ! F(1,1,1)
                        F1 = s_b*geom%R1(i) + abs(geom%l1(i))
                        F2 = s_b*geom%R2(i) + abs(geom%l2(i))
                        if (F1 /= 0. .and. F2 /= 0.) then
                            int%F111(i) = -sign(1., geom%v_eta(i))*log(F1/F2)/s_b
                        else
                            if (verbose) write(*,*) "!!! Detected evaluation point on perimeter of panel. Solution may be affected."
                        end if

                        ! Higher-order
                        if (higher_order) then
                            int%F121(i) = -(geom%v_xi(i)*geom%dR(i) + geom%a(i)*geom%v_eta(i)*int%F111(i)) / b
                            int%F211(i) = -geom%v_eta(i)*geom%dR(i) + geom%a(i)*geom%v_xi(i)*int%F111(i) - &
                                          2.*geom%v_xi(i)*geom%v_eta(i)*int%F121(i)
                        end if
                    end if

                end if

            end if

            ! Check
            if (higher_order) then
                if (abs(geom%v_xi(i)*int%F211(i) + geom%v_eta(i)*int%F121(i) - geom%a(i)*int%F111(i)) > 1.e-12) then
                    write(*,*) "!!! Calculation of F(2,1,1) and F(1,2,1) failed. Quitting..."
                    stop
                end if
            end if

        end do

    end subroutine panel_calc_supersonic_subinc_edge_integrals


    subroutine panel_calc_supersonic_supinc_edge_integrals(this, geom, dod_info, freestream, mirror_panel, int)
        ! Calculates the F integrals necessary to determine the influence of a superinclined triangular panel in supersonic flow.
        ! Taken from Epton and Magnus, but mostly the PAN AIR source code

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(integrals),intent(inout) :: int

        real :: F1, F2
        integer :: i
        logical :: higher_order

        ! Check order
        higher_order = source_order == 1 .or. doublet_order == 2

        ! Loop through edges
        do i=1,this%N

            ! Check DoD
            if (dod_info%edges_in_dod(i)) then

                ! Mach wedge
                if (geom%R1(i) == 0. .and. geom%R2(i) == 0) then

                    ! F(1,1,1)
                    int%F111(i) = -pi

                    ! Higher-order
                    if (higher_order) then
                        int%F211(i) = geom%a(i)*geom%v_xi(i)*int%F111(i)
                        int%F121(i) = geom%a(i)*geom%v_eta(i)*int%F111(i)
                    end if

                else

                    ! Calculate F factors
                    F1 = geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)
                    F2 = geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)

                    ! F(1,1,1)
                    int%F111(i) = -atan2(F1, F2)

                    ! Higher-order
                    if (higher_order) then
                        int%F211(i) = geom%a(i)*geom%v_xi(i)*int%F111(i) - geom%v_eta(i)*geom%dR(i)
                        int%F121(i) = geom%a(i)*geom%v_eta(i)*int%F111(i) + geom%v_xi(i)*geom%dR(i)
                    end if

                end if

            end if

            ! Check
            if (higher_order) then
                if (abs(geom%v_xi(i)*int%F211(i) + geom%v_eta(i)*int%F121(i) - geom%a(i)*int%F111(i)) > 1.e-12) then
                    write(*,*) "!!! Calculation of F(2,1,1) and F(1,2,1) failed. Quitting..."
                    stop
                end if
            end if

        end do

    end subroutine panel_calc_supersonic_supinc_edge_integrals


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
            int%H313 = sum(geom%v_eta*int%F121) - geom%h*int%hH113
            int%H223 = -sum(geom%v_xi*int%F121)
            int%H133 = int%H111 - sum(geom%v_eta*int%F121)

            ! Run checks
            if (abs(sum(geom%v_eta*int%F211) + int%H223) > 1e-12) then
                write(*,*) "!!! Influence calculation failed for H(2,2,3). Quitting..."
                stop
            end if

            if (abs(int%H111 - int%H313 - int%H133 - geom%h*int%hH113) > 1e-12) then
                write(*,*) "!!! Influence calculation failed for H(3,1,3) and H(1,3,3). Quitting..."
                stop
            end if
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

                    ! Mach wedge
                    if (geom%R1(i) == 0. .and. geom%R2(i) == 0.) then
                        int%hH113 = int%hH113 - pi*sign(1., geom%h*geom%v_xi(i))
                    else

                        ! Calculate F factors for supersonic edge
                        if (b > 0.) then
                            F1 = (geom%l1(i)*geom%R2(i) - geom%l2(i)*geom%R1(i)) / geom%g2(i)
                            F2 = (b*geom%R1(i)*geom%R2(i) + geom%l1(i)*geom%l2(i)) / geom%g2(i)

                        ! Calculate F factors for subsonic edge
                        else
                            F1 = geom%dR(i)*(geom%R2(i) + geom%R1(i)) / (geom%l1(i)*geom%R2(i) + geom%l2(i)*geom%R1(i))
                            F2 = (geom%g2(i) - geom%l1(i)**2 - geom%l2(i)**2) &
                                 / (b*geom%R1(i)*geom%R2(i) - geom%l1(i)*geom%l2(i))
                        end if

                        ! Calculate hH113
                        int%hH113 = int%hH113 - atan2(geom%h*geom%a(i)*F1, geom%R1(i)*geom%R2(i) + geom%h2*F2)

                    end if

                end if
            end if
        end do

        ! Calculate H(1,1,1)
        int%H111 = -geom%h*int%hH113 + sum(geom%a*int%F111)

        ! Calculate H(2,1,3) and H(1,2,3)
        int%H213 = sum(geom%v_xi*int%F111)
        int%H123 = -sum(geom%v_eta*int%F111)

        ! Calculate higher-order source integrals
        if (source_order == 1) then
            int%H211 = 0.5*(-geom%h2*int%H213 + sum(geom%a*int%F211))
            int%H121 = 0.5*(-geom%h2*int%H123 + sum(geom%a*int%F121))
        end if

        ! Calculate higher-order doublet integrals
        if (doublet_order == 2) then
            int%H313 = -int%H111 + sum(geom%v_xi*int%F211)
            int%H223 = sum(geom%v_xi*int%F121)
            int%H133 = int%H111 - sum(geom%v_eta*int%F121)

            ! Run checks
            if (abs(sum(geom%v_eta*int%F211) + int%H223) > 1e-12) then
                write(*,*) "!!! Influence calculation failed for H(2,2,3). Quitting..."
                stop
            end if

            if (abs(int%H111 + int%H313 - int%H133 - geom%h*int%hH113) > 1e-12) then
                write(*,*) "!!! Influence calculation failed for H(3,1,3) and H(1,3,3). Quitting..."
                stop
            end if

        end if

    end subroutine panel_calc_supersonic_subinc_panel_integrals


    subroutine panel_calc_supersonic_supinc_panel_integrals(this, geom, dod_info, freestream, mirror_panel, int)
        ! Calculates the necessary H integrals to determine the influence of a superinclined panel in supersonic flow.
        ! Taken from Epton and Magnus, but mostly the PAN AIR source code

        implicit none

        class(panel),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_panel
        type(integrals),intent(inout) :: int

        real :: t_dot, t_cross, X, Y
        integer :: i, i_prev

        ! Calculate hH(1,1,3)
        int%hH113 = 2.*pi

        ! Loop through corners
        do i=1,this%N

            ! Edge influence
            if (dod_info%edges_in_dod(i)) then
                int%hH113 = int%hH113 + pi
            end if

            ! Corner influence
            if (geom%R1(i) > 0.) then

                ! Cancel out edge influence
                int%hH113 = int%hH113 - pi

                ! Get index of previous corner
                if (i == 1) then
                    i_prev = this%N
                else
                    i_prev = i-1
                end if

                ! Get dot and cross products
                t_dot = geom%v_xi(i)*geom%v_xi(i_prev) + geom%v_eta(i)*geom%v_eta(i_prev)
                t_cross = -geom%v_eta(i_prev)*geom%v_xi(i) + geom%v_xi(i_prev)*geom%v_eta(i)

                ! Intermediate values
                X = geom%a(i)*geom%a(i_prev) - geom%h2*t_dot
                Y = geom%h*geom%R1(i)*t_cross

                ! Update hH113
                int%hH113 = int%hH113 - atan2(Y, -X)

            end if
        end do

        ! Calculate H(1,1,1)
        int%H111 = geom%h*int%hH113 - sum(geom%a*int%F111)

        ! Calculate H(2,1,3) and H(1,2,3)
        int%H213 = sum(geom%v_xi*int%F111)
        int%H123 = sum(geom%v_eta*int%F111)

        ! Calculate higher-order source integrals
        if (source_order == 1) then
            int%H211 = 0.5*(geom%h2*int%H213 - sum(geom%a*int%F211))
            int%H121 = 0.5*(geom%h2*int%H123 - sum(geom%a*int%F121))
        end if

        ! Calculate higher-order doublet integrals
        if (doublet_order == 2) then
            int%H313 = int%H111 + sum(geom%v_xi*int%F211)
            int%H223 = sum(geom%v_xi*int%F121)
            int%H133 = int%H111 + sum(geom%v_eta*int%F121)

            ! Run checks
            if (abs(sum(geom%v_eta*int%F211) - int%H223) > 1e-12) then
                write(*,*) "!!! Influence calculation failed for H(2,2,3). Quitting..."
                stop
            end if

            if (abs(int%H111 - int%H313 - int%H133 + geom%h*int%hH113) > 1e-12) then
                write(*,*) "!!! Influence calculation failed for H(3,1,3) and H(1,3,3). Quitting..."
                stop
            end if

        end if

    end subroutine panel_calc_supersonic_supinc_panel_integrals


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

        ! Allocate space for edge integrals
        allocate(int%F111(this%N), source=0.)
        if (source_order == 1 .or. doublet_order == 2) then
            allocate(int%F121(this%N), source=0.)
            allocate(int%F211(this%N), source=0.)
        end if

        ! Calculate necessary integrals based on the flow condition and panel type
        if (freestream%supersonic) then
            if ((mirror_panel .and. this%r_mir < 0.) .or. (.not. mirror_panel .and. this%r < 0.)) then
                call this%calc_supersonic_supinc_edge_integrals(geom, dod_info, freestream, mirror_panel, int)
                call this%calc_supersonic_supinc_panel_integrals(geom, dod_info, freestream, mirror_panel, int)
            else
                call this%calc_supersonic_subinc_edge_integrals(geom, dod_info, freestream, mirror_panel, int)
                call this%calc_supersonic_subinc_panel_integrals(geom, dod_info, freestream, mirror_panel, int)
            end if
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


    subroutine panel_calc_potential_influences(this, P, freestream, dod_info, mirror_panel, phi_s, phi_d)
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
                if ((mirror_panel .and. this%r_mir < 0.) .or. (.not. mirror_panel .and. this%r < 0.)) then
                    geom = this%calc_supersonic_supinc_geom(P, freestream, mirror_panel, dod_info)
                else
                    geom = this%calc_supersonic_subinc_geom(P, freestream, mirror_panel, dod_info)
                end if
            else
                geom = this%calc_subsonic_geom(P, freestream, mirror_panel)
            end if

            ! Get integrals
            int = this%calc_integrals(geom, 'potential', 'doublet', freestream, mirror_panel, dod_info)

            ! Source potential
            if (.not. this%in_wake) then
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
            end if

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
    
    end subroutine panel_calc_potential_influences


    subroutine panel_calc_potentials(this, P, freestream, dod_info, mirror_panel, sigma, mu, phi_s, phi_d)
        ! Calculates the potentials induced at the given point

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(3),intent(in) :: P
        type(flow),intent(in) :: freestream
        type(dod),intent(in) :: dod_info
        logical,intent(in) :: mirror_panel
        real,dimension(:),allocatable,intent(in) :: sigma, mu
        real,intent(out) :: phi_d, phi_s

        real,dimension(:),allocatable :: source_inf, doublet_inf
        real,dimension(6) :: mu_verts
        real,dimension(3) :: sigma_verts

        ! Get influences
        call this%calc_potential_influences(P, freestream, dod_info, mirror_panel, source_inf, doublet_inf)

        ! Get strengths
        sigma_verts = this%get_source_strengths(sigma, mirror_panel)
        mu_verts = this%get_doublet_strengths(mu, mirror_panel)

        ! Apply strengths to calculate potentials
        ! Source
        if (source_order == 1) then
            phi_s = sum(source_inf*sigma_verts)
        else
            phi_s = source_inf(1)*sigma_verts(1)
        end if

        ! Doublet
        if (doublet_order == 2) then
            phi_d = sum(doublet_inf*mu_verts)
        else
            phi_d = sum(doublet_inf(1:3)*mu_verts(1:3))
        end if
    
    end subroutine panel_calc_potentials


    function panel_get_source_strengths(this, sigma, mirror) result(sigma_strengths)
        ! Returns a vector of the relevant source strengths for this panel

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: sigma
        logical,intent(in) :: mirror

        real,dimension(3) :: sigma_strengths

        integer :: i, shift

        ! Check we're not in the wake
        if (.not. this%in_wake) then

            ! Determine shift to use
            if (mirror) then
                shift = size(sigma)/2
            else
                shift = 0
            end if

            ! Constant sources
            if (source_order == 0) then
                sigma_strengths(1) = sigma(this%index+shift)
                sigma_strengths(2:3) = 0.

            ! Linear sources
            else
                do i=1,this%N
                    sigma_strengths(i) = sigma(this%get_vertex_index(i)+shift)
                end do
            end if

        ! Wakes no not have source distributions
        else
            sigma_strengths = 0.
        end if
    
        
    end function panel_get_source_strengths


    function panel_get_source_dist_parameters(this, sigma, mirror) result(sigma_params)
        ! Returns a vector describing the distribution of source strength across the panel surface

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: sigma
        logical,intent(in) :: mirror
        
        real,dimension(3) :: sigma_params

        real,dimension(3) :: sigma_verts
        real,dimension(3,3) :: S_sigma_inv

        ! Get strengths
        sigma_verts = this%get_source_strengths(sigma, mirror)

        ! Constant sources
        if (source_order == 0) then
            sigma_params = sigma_verts

        ! Linear sources
        else
            if (mirror) then
                sigma_params = matmul(this%S_sigma_inv_mir, sigma_verts)
            else
                sigma_params = matmul(this%S_sigma_inv, sigma_verts)
            end if
        end if
            
    end function panel_get_source_dist_parameters


    function panel_get_doublet_strengths(this, mu, mirror) result(mu_strengths)
        ! Returns the relevant doublet strengths for this panel

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: mu
        logical,intent(in) :: mirror

        real,dimension(6) :: mu_strengths

        integer :: shift, i

        ! Determine shift
        if (mirror) then
            shift = size(mu)/2
        else
            shift = 0
        end if

        ! Get doublet strengths based on parents
        if (this%in_wake) then
            do i=1,this%N

                ! Vertices
                mu_strengths(i) = mu(this%vertices(i)%ptr%top_parent+shift) - mu(this%vertices(i)%ptr%bot_parent+shift)

                ! Midpoints
                if (doublet_order == 2) then
                    mu_strengths(i+3) = mu(this%midpoints(i)%ptr%top_parent+shift) - mu(this%midpoints(i)%ptr%bot_parent+shift)
                end if

            end do

        ! Get doublet strengths at vertices
        else
            do i=1,this%N

                ! Vertices
                mu_strengths(i) = mu(this%get_vertex_index(i)+shift)

                ! Midpoints
                if (doublet_order == 2) mu_strengths(i+3) = mu(this%get_midpoint_index(i)+shift)

            end do
        end if
        
    end function panel_get_doublet_strengths


    function panel_get_doublet_dist_parameters(this, mu, mirror) result(mu_params)
        ! Returns a vector describing the distribution of doublet strength across the panel surface

        implicit none
        
        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: mu
        logical,intent(in) :: mirror

        real,dimension(6) :: mu_params

        real,dimension(6) :: mu_verts

        ! Get doublet strengths at vertices
        mu_verts = this%get_doublet_strengths(mu, mirror)

        ! Calculate doublet parameters (derivatives)
        if (doublet_order == 2) then
            if (mirror) then
                mu_params = matmul(this%S_mu_inv_mir, mu_verts)
            else
                mu_params = matmul(this%S_mu_inv, mu_verts)
            end if
        else
            if (mirror) then
                mu_params(1:3) = matmul(this%S_mu_inv_mir, mu_verts(1:3))
            else
                mu_params(1:3) = matmul(this%S_mu_inv, mu_verts(1:3))
            end if
            mu_params(4:6) = 0.
        end if
        
    end function panel_get_doublet_dist_parameters


    function panel_get_velocity_jump(this, mu, sigma, mirrored, point) result(dv)
        ! Calculates the jump in perturbation velocity across this panel in global coordinates at the given location
        ! If a location is not given, this will default to the centroid

        implicit none

        class(panel),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: mu, sigma
        logical,intent(in) :: mirrored
        real,dimension(3),intent(in),optional :: point

        real,dimension(3) :: dv

        real,dimension(6) :: mu_params
        real,dimension(3) :: s_dir, sigma_params
        real,dimension(2) :: Q_ls
        integer :: i
        real :: s

        ! Get point
        if (present(point)) then
            if (mirrored) then
                Q_ls = matmul(this%A_g_to_ls_mir(1:2,:), point-this%centr_mir)
            else
                Q_ls = matmul(this%A_g_to_ls(1:2,:), point-this%centr)
            end if
        else
            Q_ls = 0.
        end if

        ! Calculate doublet parameters (derivatives)
        mu_params = this%get_doublet_dist_parameters(mu, mirrored)

        ! Calculate tangential velocity jump in panel coordinates E&M Eq. (N.1.11b)
        dv(1) = mu_params(2) + mu_params(4)*Q_ls(1) + mu_params(5)*Q_ls(2)
        dv(2) = mu_params(3) + mu_params(5)*Q_ls(1) + mu_params(6)*Q_ls(2)
        dv(3) = 0.

        ! Transform to global coordinates
        if (mirrored) then
            dv = matmul(transpose(this%A_g_to_ls_mir), dv)
        else
            dv = matmul(transpose(this%A_g_to_ls), dv)
        end if

        ! Get source direction
        if (mirrored) then
            s_dir = this%n_g_mir/inner(this%nu_g_mir, this%n_g_mir)
        else
            s_dir = this%n_g/inner(this%nu_g, this%n_g)
        end if

        ! Source strength
        sigma_params = this%get_source_dist_parameters(sigma, mirrored)
        s = sigma_params(1) + sigma_params(2)*Q_ls(1) + sigma_params(3)*Q_ls(2)

        ! Add normal velocity jump in global coords E&M Eq. (N.1.11b)
        dv = dv + s*s_dir

    end function panel_get_velocity_jump

    
end module panel_mod