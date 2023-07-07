module panel_solver_mod

    use helpers_mod
    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use base_geom_mod
    use surface_mesh_mod
    use flow_mod
    use math_mod
    use linalg_mod
    use sort_mod

    implicit none


    type panel_solver

        real :: M_inf_corr 
        character(len=:),allocatable :: formulation, pressure_for_forces, matrix_solver, preconditioner, iteration_file
        logical :: incompressible_rule, isentropic_rule, second_order_rule, slender_rule, linear_rule
        logical :: morino, write_A_and_b, sort_system
        logical :: compressible_correction, prandtl_glauert, karman_tsien, laitone
        type(dod),dimension(:,:),allocatable :: dod_info
        type(dod),dimension(:,:,:),allocatable :: wake_dod_info
        type(flow) :: freestream
        real :: norm_res, max_res, tol, rel
        real :: sort_time, prec_time, solver_time
        real,dimension(3) :: C_F, C_M, inner_flow
        real,dimension(:,:),allocatable :: A
        real,dimension(:),allocatable :: b, I_known, BC
        integer,dimension(:),allocatable :: P
        integer :: N_cells, block_size, max_iterations, N_unknown, N_d_unknown, N_s_unknown, solver_iterations, N_sigma
        integer :: restart_iterations, B_l_system
        logical,dimension(:),allocatable :: sigma_known
        integer,dimension(:),allocatable :: i_sigma_in_sys, i_sys_sigma_in_body

        contains

            ! Initialization
            procedure :: init => panel_solver_init
            procedure :: parse_solver_settings => panel_solver_parse_solver_settings
            procedure :: parse_processing_settings => panel_solver_parse_processing_settings

            ! Setup
            procedure :: init_dirichlet => panel_solver_init_dirichlet
            procedure :: set_panel_sources => panel_solver_set_panel_sources
            procedure :: determine_dirichlet_unknowns => panel_solver_determine_dirichlet_unknowns
            procedure :: init_dirichlet_boundary_conditions => panel_solver_init_dirichlet_boundary_conditions
            procedure :: calc_domains_of_dependence => panel_solver_calc_domains_of_dependence
            procedure :: set_permutation => panel_solver_set_permutation

            ! Solve
            procedure :: solve => panel_solver_solve
            procedure :: assemble_BC_vector => panel_solver_assemble_BC_vector
            procedure :: calc_source_strengths => panel_solver_calc_source_strengths
            procedure :: update_system_row => panel_solver_update_system_row
            procedure :: calc_body_influences => panel_solver_calc_body_influences
            procedure :: calc_wake_influences => panel_solver_calc_wake_influences
            procedure :: check_system => panel_solver_check_system
            procedure :: write_system => panel_solver_write_system
            procedure :: solve_system => panel_solver_solve_system

            ! Surface properties
            procedure :: calc_cell_velocities => panel_solver_calc_cell_velocities
            procedure :: calc_surface_potentials => panel_solver_calc_surface_potentials

            ! Pressures
            procedure :: allocate_pressure_storage => panel_solver_allocate_pressure_storage
            procedure :: calc_avg_pressure_on_panel => panel_solver_calc_avg_pressure_on_panel
            procedure :: calc_pressures => panel_solver_calc_pressures
            procedure :: output_pressures_to_terminal => panel_solver_output_pressures_to_terminal
            procedure :: calc_crit_mach => panel_solver_calc_crit_mach

            ! Force and moment integration
            procedure :: calc_forces => panel_solver_calc_forces
            procedure :: calc_moments => panel_solver_calc_moments
            procedure :: calc_forces_with_pressure => panel_solver_calc_forces_with_pressure
            procedure :: calc_moment_about_centroid_of_panel => panel_solver_calc_moment_about_centroid_of_panel

            ! Results export
            procedure :: update_report => panel_solver_update_report
            procedure :: add_pressure_to_report => panel_solver_add_pressure_to_report
            procedure :: export_off_body_points => panel_solver_export_off_body_points

    end type panel_solver


contains


    subroutine panel_solver_init(this, solver_settings, processing_settings, body, freestream, control_point_file)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: solver_settings, processing_settings
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(inout) :: freestream
        character(len=:),allocatable,intent(in) :: control_point_file

        ! Store
        this%freestream = freestream

        ! Get solver settings
        call this%parse_solver_settings(solver_settings)

        ! Get post-processing settings
        call this%parse_processing_settings(processing_settings)

        ! Initialize based on formulation
        if (this%morino .or. this%formulation == 'source-free') then
            call this%init_dirichlet(solver_settings, body)
        end if
        
        ! Write out control point geometry
        if (control_point_file /= 'none') then
            call body%write_control_points(control_point_file, solved=.false.)
        end if

        ! Calculate domains of dependence
        call this%calc_domains_of_dependence(body)

        ! Set up permutation for linear system
        call this%set_permutation(body)

    end subroutine panel_solver_init


    subroutine panel_solver_parse_solver_settings(this, solver_settings)
        ! Parses the solver settings from the input

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value),pointer, intent(in) :: solver_settings
        
        ! Get formulation
        call json_xtnsn_get(solver_settings, 'formulation', this%formulation, 'morino')        
        if (this%formulation /= 'morino' .and. this%formulation /= 'source-free') then
            write(*,*) "!!! ", this%formulation, " is not a recognized boundary condition formulation. Quitting..."
            stop
        end if
        this%morino = this%formulation == 'morino'

        ! Get matrix solver settings
        if (this%freestream%supersonic) then
            call json_xtnsn_get(solver_settings, 'matrix_solver', this%matrix_solver, 'GMRES')
        else
            call json_xtnsn_get(solver_settings, 'matrix_solver', this%matrix_solver, 'GMRES')
        end if
        call json_xtnsn_get(solver_settings, 'block_size', this%block_size, -1)
        call json_xtnsn_get(solver_settings, 'tolerance', this%tol, 1.e-12)
        call json_xtnsn_get(solver_settings, 'relaxation', this%rel, 0.8)
        call json_xtnsn_get(solver_settings, 'max_iterations', this%max_iterations, 1000)
        call json_xtnsn_get(solver_settings, 'restart_iterations', this%restart_iterations, 20)
        call json_xtnsn_get(solver_settings, 'preconditioner', this%preconditioner, 'DIAG')
        call json_xtnsn_get(solver_settings, 'iterative_solver_output', this%iteration_file, 'none')
        call json_xtnsn_get(solver_settings, 'sort_system', this%sort_system, this%freestream%supersonic)
        this%solver_iterations = -1

        ! Whether to write the linear system to file
        call json_xtnsn_get(solver_settings, 'write_A_and_b', this%write_A_and_b, .false.)
        
    end subroutine panel_solver_parse_solver_settings


    subroutine panel_solver_parse_processing_settings(this, processing_settings)
        ! Parses the post-processing settings from the input

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value),pointer, intent(in) :: processing_settings
        
        ! Get incompressible/isentropic pressure rules
        if (this%freestream%M_inf > 0.) then

            ! Isentropic is default for M > 0
            call json_xtnsn_get(processing_settings, 'pressure_rules.isentropic', this%isentropic_rule, .true.)

            ! Check for incompressible rule
            call json_xtnsn_get(processing_settings, 'pressure_rules.incompressible', this%incompressible_rule, .false.)

            ! Notify user if we're throwing out the incompressible rule
            if (this%incompressible_rule) then
                write(*,*) "!!! The incompressible pressure rule cannot be used for M > 0. Switching to the isentropic rule."
                this%incompressible_rule = .false.
                this%isentropic_rule = .true.
            end if

        else if (this%freestream%M_inf == 0.) then

            ! Incompressible is deafult
            call json_xtnsn_get(processing_settings, 'pressure_rules.incompressible', this%incompressible_rule, .true.)

            ! Check for isentropic
            call json_xtnsn_get(processing_settings, 'pressure_rules.isentropic', this%isentropic_rule, .false.)

            ! Notify user if pressure rule applied is changed based on selected freestream Mach number
            if (this%isentropic_rule) then
                write(*,*) "!!! The isentropic pressure rule cannot be used for M = 0. Switching to the incompressible rule."
                this%isentropic_rule = .false.
                this%incompressible_rule = .true.
            end if

        else

            ! MachLine does not support a negative freestream Mach number
            write(*,*) "!!! A negative freestream Mach number is not allowed. Quitting..."
            stop

        end if

        ! Get other pressure rules (these are applicable for all Mach numbers)
        call json_xtnsn_get(processing_settings, 'pressure_rules.second-order', this%second_order_rule, .false.)
        call json_xtnsn_get(processing_settings, 'pressure_rules.slender-body', this%slender_rule, .false.)
        call json_xtnsn_get(processing_settings, 'pressure_rules.linear', this%linear_rule, .false.)
        
        ! Get information for pressure corrections
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.correction_mach_number', this%M_inf_corr, 0.0)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.prandtl-glauert', this%prandtl_glauert, .false.)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.karman-tsien', this%karman_tsien, .false.)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.laitone', this%laitone, .false.)
        
        ! Check the correction Mach number is subsonic
        if (this%M_inf_corr < 0.0 .or. this%M_inf_corr >= 1.0) then
            write(*,*) "!!! The pressure correction Mach number must be between zero and one. Quitting..."
            stop
        end if

        ! Check freestream Mach number is set to 0 if pressure correction is selected
        if (((this%prandtl_glauert) .or. (this%karman_tsien) .or. (this%laitone)) .and. (this%freestream%M_inf /= 0.0)) then
            write(*,*) "!!! In order to apply a subsonic pressure correction, the freestream Mach number must be set to '0'."
            write(*,*) "!!! Include the desired Mach number for the correction calculations as 'correction_mach_number' under"
            write(*,*) "!!! 'post-processing'->'subsonic_pressure_correction'. Quitting..."
            stop
        end if

        ! Check to see if user selected incompressible rule when applying the compressible pressure corrections
        if (((this%prandtl_glauert) .or. (this%karman_tsien) .or. (this%laitone)) .and. .not. this%incompressible_rule) then
            this%incompressible_rule = .true.
        end if 

        ! Get which pressure rule will be used for force calculation
        if (this%incompressible_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'incompressible')
        else if (this%isentropic_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'isentropic')
        else if (this%second_order_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'second-order')
        else if (this%linear_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'linear')
        else if (this%slender_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'slender-body')
        else if (this%prandtl_glauert) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'prandtl-glauert')
        else if (this%karman_tsien) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'karman-tsien')
        else if (this%laitone) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'laitone')
        end if
        
    end subroutine panel_solver_parse_processing_settings


    subroutine panel_solver_init_dirichlet(this, solver_settings, body)
        ! Initializes the solver to use one of the Dirichlet formulations

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: solver_settings
        type(surface_mesh),intent(inout) :: body

        real :: offset
        character(len=:),allocatable :: offset_type

        ! Get offset
        call json_xtnsn_get(solver_settings, 'control_point_offset', offset, 1.e-7)
        call json_xtnsn_get(solver_settings, 'control_point_offset_type', offset_type, 'direct')
        if (offset <= 0.) then
            write(*,*) "!!! Control point offset must be greater than 0. Defaulting to 1e-7."
            offset = 1.e-7
        end if
        
        ! Place control points
        if (verbose) then
            if (offset_type == 'direct') then
                write(*,'(a ES10.4 a)',advance='no') "     Placing control points using a direct offset of ", offset, "..."
            else
                write(*,'(a a a ES10.4 a)',advance='no') "     Placing control points using a ", offset_type, &
                    " offset ratio of ", offset, "..."
            end if
        end if

        ! Place control points inside the body
        if (this%morino .or. this%formulation == 'source-free') then
            call body%place_interior_control_points(offset, offset_type, this%freestream)
        end if

        ! Set needed sources
        call this%set_panel_sources(body)

        ! Determine unknowns
        call this%determine_dirichlet_unknowns(body)
        if (this%N_unknown /= body%N_cp) then
            write(*,*) "!!! The number of unknowns is not the same as the number of control points. Quitting..."
            stop
        end if

        ! Set boundary conditions
        call this%init_dirichlet_boundary_conditions(body)

        ! Calculate target inner flow
        this%inner_flow = this%freestream%c_hat_g
        if (this%formulation == 'source-free') then
            this%inner_flow = this%inner_flow - matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        end if

        if (verbose) write(*,'(a, i6, a)') "Done. Placed", body%N_cp, " control points."
    
    end subroutine panel_solver_init_dirichlet


    subroutine panel_solver_set_panel_sources(this, body)
        ! Tells the panels whether they have sources or not

        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i

        ! Loop through panels
        do i=1,body%N_panels
            body%panels(i)%has_sources = (this%formulation == 'morino') .or. (body%panels(i)%r < 0.)
        end do
        
    end subroutine panel_solver_set_panel_sources


    subroutine panel_solver_determine_dirichlet_unknowns(this, body)
        ! Determines the number and type of unknown parameters to be solved for

        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(in) :: body

        integer :: i, j

        ! Determine number of source strengths
        if (body%asym_flow) then
            this%N_sigma = body%N_panels*2
        else
            this%N_sigma = body%N_panels
        end if

        ! The number of doublet unknowns is equal to the number of mesh vertices
        ! In the case of an asymmetric flow over a mirrored mesh, this must be doubled
        if (body%asym_flow) then
            this%N_d_unknown = body%N_verts*2
        else
            this%N_d_unknown = body%N_verts
        end if

        ! The number of source unknowns is simply the number of superinclined panels
        allocate(this%sigma_known(this%N_sigma), source=.true.)
        this%N_s_unknown = body%N_supinc

        ! Allocate vectors mapping unknown sigmas into the (unpermuted) linear system and back
        allocate(this%i_sigma_in_sys(this%N_sigma), source=0)
        allocate(this%i_sys_sigma_in_body(this%N_s_unknown), source=0)

        ! Create mapping
        j = this%N_d_unknown

        ! Loop through original panels
        do i=1,body%N_panels

            ! Check for superinclined
            if (body%panels(i)%r < 0) then
                j = j + 1
                this%i_sigma_in_sys(i) = j
                this%i_sys_sigma_in_body(j-this%N_d_unknown) = i
                this%sigma_known(i) = .false.
            end if
        end do

        ! Loop through mirrored panels
        if (body%asym_flow) then
            do i=1,body%N_panels

                ! Check for superinclined
                if (body%panels(i)%r_mir < 0) then
                    j = j + 1
                    this%i_sigma_in_sys(i+body%N_panels) = j
                    this%i_sys_sigma_in_body(j-this%N_d_unknown) = i+body%N_panels
                    this%sigma_known(i+body%N_panels) = .false.
                end if
            end do
        end if

        ! Total number of unknowns
        this%N_unknown = this%N_d_unknown + this%N_s_unknown
        
    end subroutine panel_solver_determine_dirichlet_unknowns


    subroutine panel_solver_init_dirichlet_boundary_conditions(this, body)
        ! Sets up the desired boundary conditions on the control points

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body

        integer :: i
    
        ! Loop through control points
        !$OMP parallel do
        do i=1,body%N_cp

            ! Check if we need to enforce strength matching
            if (body%cp(i)%is_mirror) then
                if (body%cp(i)%tied_to_type == 1) then
                    if (.not. body%vertices(body%cp(i)%tied_to_index)%mirrored_is_unique) then
                        call body%cp(i)%set_bc(4)
                        cycle
                    end if
                end if
            end if

            select case (body%cp(i)%cp_type)

            ! Interior control points
            case (1)

                if (this%morino) then
                    call body%cp(i)%set_bc(1)
                else
                    call body%cp(i)%set_bc(2)
                end if

            ! Exterior control points
            case (2)

                write(*,*) "!!! MachLine cannot handle exterior control points. Quitting..."
                stop

            ! Surface control points
            case (3)

                ! Need to get correct normal vector
                if (body%cp(i)%tied_to_type ==1 ) then
                    if (body%cp(i)%is_mirror) then
                        call body%cp(i)%set_bc(3, body%vertices(body%cp(i)%tied_to_index)%n_g_mir)
                    else
                        call body%cp(i)%set_bc(3, body%vertices(body%cp(i)%tied_to_index)%n_g)
                    end if
                else
                    if (body%cp(i)%is_mirror) then
                        call body%cp(i)%set_bc(3, body%panels(body%cp(i)%tied_to_index)%n_g_mir)
                    else
                        call body%cp(i)%set_bc(3, body%panels(body%cp(i)%tied_to_index)%n_g)
                    end if
                end if
            
            end select

        end do
        
    end subroutine panel_solver_init_dirichlet_boundary_conditions


    subroutine panel_solver_calc_domains_of_dependence(this, body)
        ! Determines the domains of dependence for each control point based on the freestream condition

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k, stat, N_panels, N_strip_panels
        real,dimension(3) :: vert_loc, mirrored_vert_loc

        if (this%freestream%supersonic .and. verbose) write(*,'(a)',advance='no') "     Calculating domains of dependence..."

        ! Figure out how many verts/panels we're going to consider
        if (body%mirrored) then
            if (body%asym_flow) then
                N_strip_panels = body%wake%N_max_strip_panels
            else
                N_strip_panels = 2*body%wake%N_max_strip_panels
            end if
            N_panels = 2*body%N_panels
        else
            N_panels = body%N_panels
            N_strip_panels = body%wake%N_max_strip_panels
        end if

        ! Allocate

        ! DoD info for panels
        allocate(this%dod_info(N_panels, this%N_unknown), stat=stat)
        call check_allocation(stat, "domain of dependence storage")

        ! DoD info for panels
        allocate(this%wake_dod_info(N_strip_panels, body%wake%N_strips, this%N_unknown), stat=stat)
        call check_allocation(stat, "wake domain of dependence storage")

        ! If the freestream is subsonic, these don't need to be checked
        if (this%freestream%supersonic) then

            ! Loop through control points
            !$OMP parallel do private(i, k, vert_loc, mirrored_vert_loc)
            do j=1,body%N_cp

                ! Loop through body panels
                do i=1,body%N_panels

                    ! Check DoD for original panel
                    this%dod_info(i,j) = body%panels(i)%check_dod(body%cp(j)%loc, this%freestream)

                    if (body%mirrored) then

                        ! Check DoD for mirrored panel
                        this%dod_info(i+body%N_panels,j) = body%panels(i)%check_dod(body%cp(j)%loc, this%freestream, &
                                                                                    .true.)
                    end if
                end do

                ! Loop through wake strip panels
                do i=1,body%wake%N_strips
                    do k=1,body%wake%strips(i)%N_panels

                        ! Check DoD for panel
                        this%wake_dod_info(k,i,j) = body%wake%strips(i)%panels(k)%check_dod(body%cp(j)%loc, this%freestream)

                        if (body%mirrored .and. .not. body%asym_flow) then

                            ! Check DoD for mirrored panel
                            this%wake_dod_info(k+body%wake%N_max_strip_panels,i,j) = &
                                body%wake%strips(i)%panels(k)%check_dod(body%cp(j)%loc, this%freestream, &
                                                                        .true.)

                        end if
                    end do
                end do

            end do

            if (verbose) write(*,*) "Done."
        end if
    
    end subroutine panel_solver_calc_domains_of_dependence


    subroutine panel_solver_set_permutation(this, body)
        ! Creates a vertex/control point permutation which should lead to a more efficient solution of the matrix equation

        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(in) :: body

        real,dimension(:),allocatable :: x
        real,dimension(3) :: loc
        integer :: i, j, i_neighbor, i_cp, i_vert, i_panel, i_vert_for_panel, i2, i3, i_panel_abutting, i_opp_edge
        integer,dimension(:),allocatable :: P_inv_1, P_inv_2
        integer(8) :: start_count, end_count
        real(16) :: count_rate

        ! Sort control points in the compressibility direction
        ! We do this using the location of the vertex/panel centroid tied to each control point
        ! We proceed from most downstream to most upstream so as to get an upper-pentagonal matrix
        ! This sorting sometimes improves the performance of iterative solvers for supersonic flow as well
        if (this%sort_system) then

            if (verbose) write(*,'(a)',advance='no') "     Permuting linear system for efficient solution..."

            ! Sort by actual vertex location first
            call system_clock(start_count, count_rate)
            allocate(x(this%N_unknown))

            ! Add compressibility distance of each vertex/panel centroid
            do i=1,body%N_cp

                ! Get location
                if (body%cp(i)%is_mirror) then
                    if (body%cp(i)%tied_to_type == 1) then
                        loc = mirror_across_plane(body%vertices(body%cp(i)%tied_to_index)%loc, body%mirror_plane)
                    else
                        loc = body%panels(body%cp(i)%tied_to_index)%centr_mir
                    end if
                else
                    if (body%cp(i)%tied_to_type == 1) then
                        loc = body%vertices(body%cp(i)%tied_to_index)%loc
                    else
                        loc = body%panels(body%cp(i)%tied_to_index)%centr
                    end if
                end if

                ! Get compressibility distance
                x(i) = -inner(this%freestream%c_hat_g, loc)

            end do

            ! Get inverse permutation based on vertex location
            call insertion_arg_sort(x, P_inv_1)

            ! Now sort again based on location of most-downstream neighboring vertex
            ! For vertices sharing a downstream vertex, this will preserve the previous order, since insertion sort is stable
            do i=1,size(P_inv_1)

                ! Get control point index
                i_cp = P_inv_1(i)

                ! Mirrored control point
                if (body%cp(i_cp)%is_mirror) then

                    ! Tied to vertex
                    if (body%cp(i_cp)%tied_to_type == 1) then

                        ! Initialize
                        i_vert = body%cp(i_cp)%tied_to_index ! Index of original, not mirrored, vertex
                        x(i) = huge(x(i))

                        ! Loop through neighboring vertices to find the furthest back
                        do j=1,body%vertices(i_vert)%adjacent_vertices%len()

                            ! Get neighbor index
                            call body%vertices(i_vert)%adjacent_vertices%get(j, i_neighbor)

                            ! Update sorting parameter
                            x(i) = min(x(i), -inner(this%freestream%c_hat_g, mirror_across_plane(body%vertices(i_neighbor)%loc, &
                                                                                                 body%mirror_plane)))
                        end do

                        ! For higher-order distributions, loop through panels abutting neighboring panels as well (not across discontinuous edges)
                        abutting_panel_loop_mir: do j=1,body%vertices(i_vert)%panels_not_across_wake_edge%len()

                            ! Get panel index
                            call body%vertices(i_vert)%panels_not_across_wake_edge%get(j, i_panel)

                            ! Check if this panel has a higher-order distribution
                            if (body%panels(i_panel)%order == 2) then

                                ! Get index of the edge opposite this vertex
                                i_opp_edge = body%panels(i_panel)%get_opposite_edge(i_vert)
                            
                                ! Check whether this edge is discontinuous
                                if (body%edges(i_opp_edge)%discontinuous) cycle abutting_panel_loop_mir

                                ! Check whether this edge is on the mirror plane; if so, then the vertex we're looking for is the original one, and we can skip this
                                if (body%edges(i_opp_edge)%on_mirror_plane) cycle abutting_panel_loop_mir

                                ! Get abutting panel
                                if (body%edges(i_opp_edge)%panels(1) == i_panel) then
                                    i_panel_abutting = body%edges(i_opp_edge)%panels(2)
                                else
                                    i_panel_abutting = body%edges(i_opp_edge)%panels(1)
                                end if

                                ! Get index of vertex opposite this edge on the abutting panel
                                i_neighbor = &
                                    body%panels(i_panel_abutting)%get_opposite_vertex(body%edges(i_opp_edge)%top_verts(1), &
                                                                                      body%edges(i_opp_edge)%top_verts(2)) ! We can always use top_verts here because no wake-shedding edge will be continuous

                                ! Update sorting parameter
                                x(i) = min(x(i), -inner(this%freestream%c_hat_g, mirror_across_plane(body%vertices(i_neighbor)%loc,&
                                                                                                     body%mirror_plane)))
                            end if

                        end do abutting_panel_loop_mir

                    ! Tied to panel
                    else

                        ! Loop through this panel's vertices to find the furthest back
                        i_panel = body%cp(i_cp)%tied_to_index
                        x(i) = huge(x(i))
                        do j=1,body%panels(i_panel)%N
                            loc = mirror_across_plane(body%panels(i_panel)%get_vertex_loc(j), body%mirror_plane)
                            x(i) = min(x(i), -inner(this%freestream%c_hat_g, loc))
                        end do

                        ! TODO: For higher-order distributions, loop through abutting panels as well (not across discontinuous edges)

                    end if

                else

                    ! Tied to vertex
                    if (body%cp(i_cp)%tied_to_type == 1) then

                        ! Initialize
                        i_vert = body%cp(i_cp)%tied_to_index
                        x(i) = huge(x(i))

                        ! Loop through neighboring vertices to find the furthest back
                        do j=1,body%vertices(i_vert)%adjacent_vertices%len()

                            ! Get index of neighbor
                            call body%vertices(i_vert)%adjacent_vertices%get(j, i_neighbor)

                            ! Update sorting parameter
                            x(i) = min(x(i), -inner(this%freestream%c_hat_g, body%vertices(i_neighbor)%loc))

                        end do

                        ! For higher-order distributions, loop through panels abutting neighboring panels as well (not across discontinuous edges)
                        abutting_panel_loop: do j=1,body%vertices(i_vert)%panels_not_across_wake_edge%len()

                            ! Get panel index
                            call body%vertices(i_vert)%panels_not_across_wake_edge%get(j, i_panel)

                            ! Check if this panel has a higher-order distribution
                            if (body%panels(i_panel)%order == 2) then

                                ! Get index of the edge opposite this vertex
                                i_opp_edge = body%panels(i_panel)%get_opposite_edge(i_vert)
                                
                                ! Check whether this edge is discontinuous
                                if (body%edges(i_opp_edge)%discontinuous) cycle abutting_panel_loop

                                ! Check whether this edge is on the mirror plane; if so, then the vertex we're looking for is the original one, and we can skip this
                                if (body%edges(i_opp_edge)%on_mirror_plane) cycle abutting_panel_loop

                                ! Get abutting panel
                                if (body%edges(i_opp_edge)%panels(1) == i_panel) then
                                    i_panel_abutting = body%edges(i_opp_edge)%panels(2)
                                else
                                    i_panel_abutting = body%edges(i_opp_edge)%panels(1)
                                end if

                                ! Get index of vertex opposite this edge on the abutting panel
                                i_neighbor = &
                                    body%panels(i_panel_abutting)%get_opposite_vertex(body%edges(i_opp_edge)%top_verts(1), &
                                                                                      body%edges(i_opp_edge)%top_verts(2)) ! We can always use top_verts here because no wake-shedding edge will be continuous

                                ! Update sorting parameter
                                x(i) = min(x(i), -inner(this%freestream%c_hat_g, body%vertices(i_neighbor)%loc))

                            end if

                        end do abutting_panel_loop

                    ! Tied to panel
                    else

                        ! Loop through this panel's vertices to find the furthest back
                        i_panel = body%cp(i_cp)%tied_to_index
                        x(i) = huge(x(i))
                        do j=1,body%panels(i_panel)%N
                            loc = body%panels(i_panel)%get_vertex_loc(j)
                            x(i) = min(x(i), -inner(this%freestream%c_hat_g, loc))
                        end do

                        ! TODO: For higher-order distributions, loop through abutting panels as well (not across discontinuous edges)

                    end if

                end if

            end do

            ! Get inverse permutation
            ! Sorts into increasing order
            call insertion_arg_sort(x, P_inv_2)

            ! Get overall permuation
            allocate(this%P(this%N_unknown))
            do i=1,this%N_unknown
                this%P(P_inv_1(P_inv_2(i))) = i
            end do

            ! Get timing
            call system_clock(end_count)

            if (verbose) write(*,*) "Done."

        else

            ! Only permute as needed for placing vertex clones next to each other
            call system_clock(start_count, count_rate)
            allocate(this%P(this%N_unknown))

            ! Copy vertex ordering
            this%P(1:body%N_verts) = body%vertex_ordering

            ! Copy identity for unknown sources
            do i=1,body%N_supinc
                this%P(body%N_verts+i) = body%N_verts + i
            end do

            if (body%asym_flow) then

                ! For a mirrored, asymmetric condition, copy the vertex ordering again
                this%P(body%N_verts+body%N_supinc+1:2*body%N_verts+body%N_supinc) = body%vertex_ordering + body%N_verts

                ! Copy identity for unknown sources
                do i=1,body%N_supinc
                    this%P(2*body%N_verts+body%N_supinc+i) = 2*body%N_verts+body%N_supinc + i
                end do
            end if
            call system_clock(end_count)

        end if

        ! Calculate elapsed time
        this%sort_time = real(end_count - start_count)/count_rate
    
    end subroutine panel_solver_set_permutation


    subroutine panel_solver_solve(this, body, solver_stat)
        ! Calls the relevant subroutine to solve the case based on the selected formulation
        ! We are solving the equation
        !
        !    [A] | mu    | + I_known = BC
        !        | sigma |
        !
        ! Thus, b = BC - I_known
        !
        ! I_known represents known singularity influences
        ! BC represents target values for each control point based on the boundary condition enforced there

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        integer,intent(out) :: solver_stat

        integer :: stat

        ! Set default status
        solver_stat = 0

        ! Allocate known influence storage
        allocate(this%I_known(this%N_unknown), source=0., stat=stat)
        call check_allocation(stat, "known influence vector")

        ! Allocate AIC matrix
        allocate(this%A(this%N_unknown, this%N_unknown), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")

        ! Allocate b vector
        allocate(this%b(this%N_unknown), source=0., stat=stat)
        call check_allocation(stat, "b vector")

        ! Calculate source strengths
        call this%calc_source_strengths(body)

        ! Calculate body influences
        call this%calc_body_influences(body)

        ! Calculate wake influences
        if (body%wake%N_panels > 0) call this%calc_wake_influences(body)

        ! Assemble boundary condition vector
        call this%assemble_BC_vector(body)

        ! Solve the linear system
        call this%solve_system(body, solver_stat)
        
        ! Check for errors
        if (solver_stat /= 0) return

        ! Calculate velocities
        call this%calc_cell_velocities(body)
        
        ! Calculate potentials
        call this%calc_surface_potentials(body)

        ! Calculate pressures
        call this%calc_pressures(body)

        ! Calculate forces
        call this%calc_forces(body)

        ! Calculate moments
        call this%calc_moments(body)

    end subroutine panel_solver_solve


    subroutine panel_solver_assemble_BC_vector(this, body)
        ! Sets up the {BC} vector used in the linear system synthesis

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(in) :: body

        integer :: i
        real,dimension(3) :: x

        ! Allocate
        allocate(this%BC(this%N_unknown))

        ! Vector for source-free target potential calculation
        x = matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)

        ! Loop through control points
        do i=1,body%N_cp
            
            ! Check boundary condition type
            select case (body%cp(i)%bc)

            ! Source-free formulation
            case (2)
                this%BC(this%P(i)) = -inner(x, body%cp(i)%loc)

            ! All other cases
            ! 1: For Morino, desired BC is zero inner potential
            ! 3: Velocity BC not yet implemented
            ! 4: For strength-matching, desired BC is zero difference between strengths
            case default
                this%BC(this%P(i)) = 0.

            end select
        end do
        
    end subroutine panel_solver_assemble_BC_vector


    subroutine panel_solver_calc_source_strengths(this, body)
        ! Calculates the necessary source strengths for subinclined panels

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat

        ! Allocate source strength array
        allocate(body%sigma(this%N_sigma), source=0., stat=stat)
        call check_allocation(stat, "source strength array")

        ! Set source strengths
        if (this%morino) then

            if (verbose) write(*,'(a)',advance='no') "     Calculating source strengths..."

            ! Loop through panels
            do i=1,body%N_panels

                ! Existing panels
                if (this%sigma_known(i)) then
                    body%sigma(i) = -inner(body%panels(i)%n_g, this%freestream%c_hat_g)
                end if

                ! Mirrored panels for asymmetric flow
                if (body%asym_flow) then
                    if (this%sigma_known(i+body%N_panels)) then
                        body%sigma(i+body%N_panels) = -inner(body%panels(i)%n_g_mir, this%freestream%c_hat_g)
                    end if
                end if
            end do

            if (verbose) write(*,*) "Done."
        end if
    
    end subroutine panel_solver_calc_source_strengths


    subroutine panel_solver_update_system_row(this, body, cp, A_row, I_known_i, i_panel, source_inf, doublet_inf, mirrored_panel)
        ! Updates the linear system with the source and doublet influences

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(control_point),intent(in) :: cp
        real,dimension(this%N_unknown),intent(inout) :: A_row
        real,intent(inout) :: I_known_i
        integer,intent(in) :: i_panel
        real,dimension(:),allocatable,intent(in) :: source_inf, doublet_inf
        logical,intent(in) :: mirrored_panel

        integer :: k, index

        ! Add source influence depending on boundary condition
        if (cp%bc == 1 .or. cp%bc == 3) then

            ! Loop through influencing panels
            do k=1,body%panels(i_panel)%S_dim

                ! Need to shift indices for mirrored panels
                if (mirrored_panel) then

                    ! Check for dependence across mirror plane
                    if (body%panels(i_panel)%i_panel_s(k) > body%N_panels) then
                        index = body%panels(i_panel)%i_panel_s(k) - body%N_panels
                    else
                        index = body%panels(i_panel)%i_panel_s(k) + body%N_panels
                    end if

                else
                    ! Check for dependence across mirror plane
                    if (body%panels(i_panel)%i_panel_s(k) > body%N_panels) then
                        index = body%panels(i_panel)%i_panel_s(k) - body%N_panels
                    else
                        index = body%panels(i_panel)%i_panel_s(k)
                    end if
                end if

                ! Add to known influences if sigma is known
                if (this%sigma_known(index)) then
                    I_known_i = I_known_i + source_inf(k)*body%sigma(index)

                ! Add to A matrix if not
                else
                    A_row(this%P(this%i_sigma_in_sys(index))) = A_row(this%P(this%i_sigma_in_sys(index))) + source_inf(k)
                end if

            end do

        end if

        ! Add doublet influence

        ! Loop through influencing vertices
        do k=1,body%panels(i_panel)%M_dim

            ! Need to shift indices for mirrored panels
            if (mirrored_panel) then

                ! Check for dependence across mirror plane
                if (body%panels(i_panel)%i_vert_d(k) > body%N_verts) then
                    index = body%panels(i_panel)%i_vert_d(k) - body%N_verts
                else
                    index = body%panels(i_panel)%i_vert_d(k) + body%N_verts
                end if

            else

                ! Check for dependence across mirror plane
                if (body%panels(i_panel)%i_vert_d(k) > body%N_verts) then
                    index = body%panels(i_panel)%i_vert_d(k) - body%N_verts
                else
                    index = body%panels(i_panel)%i_vert_d(k)
                end if
            end if

            ! Update
            A_row(this%P(index)) = A_row(this%P(index)) + doublet_inf(k)

        end do
    
    end subroutine panel_solver_update_system_row


    subroutine panel_solver_calc_body_influences(this, body)
        ! Calculates the influence of the body on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j
        real,dimension(:),allocatable :: source_inf, doublet_inf
        real,dimension(this%N_unknown) :: A_i
        real :: I_known_i

        if (verbose) write(*,'(a)',advance='no') "     Calculating body influences..."

        ! Calculate source and doublet influences from body on each control point
        !$OMP parallel do private(j, source_inf, doublet_inf, A_i, I_known_i) schedule(dynamic)
        do i=1,body%N_cp

            ! Initialize
            A_i = 0.
            I_known_i = 0.

            ! Determine the type of boundary condition on this control point
            select case (body%cp(i)%bc)

            case (4) ! Strength matching

                A_i(this%P(i)) = 1.
                A_i(this%P(i-body%N_cp/2)) = -1.

            case (3) ! Calculate velocity influences

                ! TODO: Do this

            case default ! Calculate potential influences

                ! Loop through panels
                do j=1,body%N_panels

                    ! Influence of existing panel on control point
                    if (this%dod_info(j,i)%in_dod) then

                        ! Calculate influence
                        call body%panels(j)%calc_potential_influences(body%cp(i)%loc, this%freestream, this%dod_info(j,i), &
                                                                      .false., source_inf, doublet_inf)

                        ! Add influence
                        call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, .false.)

                    end if

                    if (body%mirrored) then

                        ! Calculate influence of mirrored panel on control point
                        if (this%dod_info(j+body%N_panels,i)%in_dod) then

                            ! Calculate influence
                            call body%panels(j)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
                                                                          this%dod_info(j+body%N_panels,i), &
                                                                          .true., source_inf, doublet_inf)

                            ! Add influence
                            call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, &
                                                        source_inf, doublet_inf, body%asym_flow)
                        end if

                    end if

                end do

            end select

            ! Update system
            !$OMP critical

            ! Update A matrix
            this%A(this%P(i),:) = A_i

            ! Update I_known
            this%I_known(this%P(i)) = I_known_i

            !$OMP end critical

        end do

        ! Clean up memory
        deallocate(this%dod_info)

        if (verbose) write(*,*) "Done."
    
    end subroutine panel_solver_calc_body_influences


    subroutine panel_solver_calc_wake_influences(this, body)
        ! Calculates the influence of the wake on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k, l
        real,dimension(:),allocatable ::  doublet_inf, source_inf
        real,dimension(this%N_unknown) :: A_i

        ! Calculate influence of wake
        if (verbose) write(*,'(a)',advance='no') "     Calculating wake influences..."

        ! Loop through control points
        !$OMP parallel do private(j, k, l, source_inf, doublet_inf, A_i) schedule(dynamic)
        do i=1,body%N_cp

            ! Check boundary condition
            select case (body%cp(i)%bc)

            case (4) ! Strength matching
                cycle

            case (3) ! Calculate velocity influences

            case default ! Calculate potential influences

                ! Initialize
                A_i = 0.

                ! Get doublet influence from wake strips
                do j=1,body%wake%N_strips
                    do l=1,body%wake%strips(j)%N_panels

                        ! Caclulate influence of existing panel on control point
                        call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
                                                                                     this%wake_dod_info(l,j,i), &
                                                                                     .false., source_inf, doublet_inf)

                        ! Add influence
                        do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                        end do

                        ! Get influence of mirrored panel
                        if (body%wake%strips(j)%mirrored) then

                            ! Calculate influence of mirrored panel on control point
                            call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
                                                                     this%wake_dod_info(l+body%wake%N_max_strip_panels,j,i), & ! No, this is not the DoD for this computation; yes, it is equivalent
                                                                     .true., source_inf, doublet_inf)

                            ! Add influence
                            do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                            end do

                        end if
                    end do
                end do

            end select

            ! Update row of A
            !$OMP critical
            this%A(this%P(i),:) = this%A(this%P(i),:) + A_i
            !$OMP end critical

        end do

        ! Clean up memory
        deallocate(this%wake_dod_info)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_wake_influences


    subroutine panel_solver_check_system(this, solver_stat)
        ! Checks the validity of the linear system

        implicit none
        
        class(panel_solver),intent(in) :: this
        integer,intent(inout) :: solver_stat

        integer :: i

        if (verbose) write(*,'(a)',advance='no') "     Checking validity of linear system..."

        ! Check for NaNs
        if (any(isnan(this%A))) then
            write(*,*) "!!! Invalid value detected in A matrix."
            solver_stat = 1
            return
        end if
        if (any(isnan(this%b))) then
            write(*,*) "!!! Invalid value detected in b vector."
            solver_stat = 1
            return
        end if

        ! Check for uninfluenced/ing points
        do i=1,this%N_unknown
            if (all(this%A(this%P(i),:) == 0.)) then
                write(*,*) "!!! Control point ", i, " is not influenced."
                solver_stat = 2
                return
            end if
            if (all(this%A(:,this%P(i)) == 0.)) then
                write(*,*) "!!! Vertex ", i, " exerts no influence."
                solver_stat = 2
                return
            end if
        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_check_system


    subroutine panel_solver_write_system(this)
        ! Writes the linear system to file

        implicit none
        
        class(panel_solver),intent(in) :: this

        integer :: i, j, unit

        if (verbose) write(*,'(a)',advance='no') "     Writing linear system to file..."

        ! A
        open(newunit=unit, file="A_mat.txt")
        do i=1,this%N_unknown
            do j=1,this%N_unknown
                write(unit,'(e20.12)',advance='no') this%A(i,j)
            end do
            write(unit,*)
        end do
        close(unit)

        ! b
        open(newunit=unit, file="b_vec.txt")
        do i=1,this%N_unknown
            write(unit,*) this%b(i)
        end do
        close(unit)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_write_system


    subroutine panel_solver_solve_system(this, body, solver_stat)
        ! Solves the linear system for the singularity strengths

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        integer,intent(inout) :: solver_stat

        real,dimension(:,:),allocatable :: A_p
        real,dimension(:),allocatable :: b_p, x
        integer :: stat, i
        integer(8) :: start_count, end_count
        real(16) :: count_rate

        ! Set b vector
        this%b = this%BC - this%I_known

        ! Run checks
        if (run_checks) call this%check_system(solver_stat)

        ! Calculate lower bandwidth
        if (this%sort_system) this%B_l_system = get_lower_bandwidth(this%N_unknown, this%A)

        ! Write to file
        if (this%write_A_and_b) call this%write_system()

        if (verbose) write(*,'(a, a, a)',advance='no') "     Solving linear system (method: ", this%matrix_solver, ")..."

        ! Precondition
        select case(this%preconditioner)

        ! Diagonal preconditioning
        case ('DIAG')
            call system_clock(start_count, count_rate)
            call diagonal_preconditioner(this%N_unknown, this%A, this%b, A_p, b_p)
            call system_clock(end_count)

        ! No preconditioning
        case default
            call system_clock(start_count, count_rate)
            allocate(A_p, source=this%A, stat=stat)
            call check_allocation(stat, "solver copy of AIC matrix")
            allocate(b_p, source=this%b, stat=stat)
            call check_allocation(stat, "solver copy of b vector")
            call system_clock(end_count)

        end select

        ! Calculate how much time preconditioning took
        this%prec_time = real(end_count-start_count)/count_rate

        ! Check block size
        if (this%block_size <= 0) then
            this%block_size = this%N_unknown / 5
        end if

        ! Solve
        select case(this%matrix_solver)

        ! LU decomposition
        case ('LU')
            call system_clock(start_count, count_rate)
            call lu_solve(this%N_unknown, A_p, b_p, x)
            call system_clock(end_count)

        ! QR via Givens rotations for upper-pentagonal
        case ('QRUP')
            call system_clock(start_count, count_rate)
            call QR_givens_solve_UP(this%N_unknown, A_p, b_p, x)
            call system_clock(end_count)

        ! QR via fast Givens rotations for upper-pentagonal
        case ('FQRUP')
            call system_clock(start_count, count_rate)
            call QR_fast_givens_solve_upper_pentagonal(this%N_unknown, A_p, b_p, x)
            call system_clock(end_count)

        ! GMRES
        case ('GMRES')
            call system_clock(start_count, count_rate)
            call GMRES(this%N_unknown, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        ! Restarted GMRES
        case ('RGMRES')
            call system_clock(start_count, count_rate)
            call restarted_GMRES(this%N_unknown, A_p, b_p, this%tol, this%max_iterations, this%restart_iterations, &
                                 this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        ! Purcell's method
        case ('PURC')
            call system_clock(start_count, count_rate)
            call purcell_solve(this%N_unknown, A_p, b_p, x)
            call system_clock(end_count)

        ! Block symmetric successive over-relaxation
        case ('BSSOR')
            call system_clock(start_count, count_rate)
            call block_ssor_solve(this%N_unknown, A_p, b_p, this%block_size, this%tol, this%rel, &
                                  this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)
        
        ! Block Jacobi
        case ('BJAC')
            call system_clock(start_count, count_rate)
            call block_jacobi_solve(this%N_unknown, A_p, b_p, this%block_size, this%tol, this%rel, &
                                    this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        ! Improper specification
        case default
            write(*,*) "!!! ", this%matrix_solver, " is not a valid solver option. Defaulting to GMRES."
            call system_clock(start_count, count_rate)
            call GMRES(this%N_unknown, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        end select
        if (verbose) write(*,*) "Done."

        ! Calculate how much time the matrix solver took
        this%solver_time = real(end_count-start_count)/count_rate

        ! Clean up memory
        deallocate(A_p)
        deallocate(b_p)

        ! Get residual vector
        allocate(body%R_cp(this%N_unknown))
        body%R_cp = matmul(this%A, x) - this%b
        deallocate(this%A)
        deallocate(this%b)

        ! Calculate residual parameters
        this%max_res = maxval(abs(body%R_cp))
        this%norm_res = sqrt(sum(body%R_cp*body%R_cp))
        if (verbose) then
            write(*,*) "        Maximum residual:", this%max_res
            write(*,*) "        Norm of residual:", this%norm_res
            if (this%sort_system) write(*,*) "        Lower bandwidth of A matrix:", this%B_l_system
        end if

        ! Check
        if (isnan(this%norm_res)) then
            write(*,*) "!!! Linear system failed to produce a valid solution."
            solver_stat = 4
            return
        end if

        ! Transfer solved doublet strengths to body storage
        if (body%asym_flow) then
            allocate(body%mu(body%N_verts*2))
        else
            allocate(body%mu(body%N_verts))
        end if
        do i=1,this%N_d_unknown
            body%mu(i) = x(this%P(i))
        end do

        ! Transfer solved source strengths to body storage
        do i=1,this%N_s_unknown
            body%sigma(this%i_sys_sigma_in_body(i)) = x(this%P(this%N_d_unknown+i))
        end do

    end subroutine panel_solver_solve_system


    subroutine panel_solver_calc_cell_velocities(this, body)
        ! Calculates the surface velocities on the mesh panels

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat, j
        real,dimension(3) :: v_jump, cent

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface velocities..."

        ! Determine number of surface cells
        this%N_cells = body%N_panels
        if (body%asym_flow) then
            this%N_cells = this%N_cells*2
        end if

        ! Allocate cell velocity storage
        allocate(body%V_cells(3,this%N_cells), stat=stat)
        call check_allocation(stat, "surface velocity vectors")

        ! Get the surface velocity on each existing panel
        !$OMP parallel do private(v_jump, j, cent) schedule(static)
        do i=1,body%N_panels
            body%V_cells(:,i) = body%panels(i)%get_velocity(body%mu, body%sigma, .false., body%N_panels, &
                                                            body%N_verts, body%asym_flow, this%freestream, this%inner_flow)

            ! Get surface velocity on each mirrored panel
            if (body%asym_flow) then
                body%V_cells(:,i+body%N_panels) = body%panels(i)%get_velocity(body%mu, body%sigma, .true., body%N_panels, &
                                                                              body%N_verts, body%asym_flow, this%freestream, &
                                                                              this%inner_flow)

            end if
        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_cell_velocities


    subroutine panel_solver_calc_surface_potentials(this, body)
        ! Calculates the total outer potential at mesh vertices

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: loc_mir

        if (verbose) write(*,'(a)',advance='no') "     Calculating outer surface potentials..."

        ! Allocate storage
        allocate(body%Phi_u, source=body%mu, stat=stat)
        call check_allocation(stat, "total outer surface potential")

        ! Calculate total potential on outside of mesh and quadratic-doublet velocity discrepancies
        !$OMP parallel do private(loc_mir) schedule(dynamic)
        do i=1,body%N_verts

            ! Existing points
            body%Phi_u(i) = body%Phi_u(i) + inner(this%inner_flow, body%vertices(i)%loc)

            ! Mirrored points
            if (body%asym_flow) then
                loc_mir = mirror_across_plane(body%vertices(i)%loc, body%mirror_plane)
                body%Phi_u(i+body%N_verts) = body%Phi_u(i+body%N_verts) + inner(this%inner_flow, loc_mir)
            end if
        end do

        ! Factor in freestream velocity magnitude
        body%Phi_u = body%Phi_u*this%freestream%U

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_surface_potentials


    subroutine panel_solver_allocate_pressure_storage(this, body)
        ! Allocates vectors for storing pressure results

        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: stat

        if (this%incompressible_rule) then
            allocate(body%C_p_inc(this%N_cells), stat=stat)
            call check_allocation(stat, "incompressible surface pressures")
        end if
        
        if (this%isentropic_rule) then
            allocate(body%C_p_ise(this%N_cells), stat=stat)
            call check_allocation(stat, "isentropic surface pressures")
        end if
        
        if (this%second_order_rule) then
            allocate(body%C_p_2nd(this%N_cells), stat=stat)
            call check_allocation(stat, "second-order surface pressures")
        end if
        
        if (this%slender_rule) then
            allocate(body%C_p_sln(this%N_cells), stat=stat)
            call check_allocation(stat, "slender-body surface pressures")
        end if
        
        if (this%linear_rule) then
            allocate(body%C_p_lin(this%N_cells), stat=stat)
            call check_allocation(stat, "linear surface pressures")
        end if
        
        if (this%prandtl_glauert) then
            allocate(body%C_p_pg(this%N_cells), stat=stat)
            call check_allocation(stat, "Prandtl-Glauert corrected surface pressures")
        end if
        
        if (this%laitone) then
            allocate(body%C_p_lai(this%N_cells), stat=stat)
            call check_allocation(stat, "Laitone corrected surface pressures")
        end if
        
        if (this%karman_tsien) then
            allocate(body%C_p_kt(this%N_cells), stat=stat)
            call check_allocation(stat, "Karman-Tsien corrected surface pressures")
        end if
        
    end subroutine panel_solver_allocate_pressure_storage


    function panel_solver_calc_avg_pressure_on_panel(this, i_panel, body, mirrored, rule) result(C_P_avg)
        ! Calls the panel function with the needed variables to get the average pressure coefficient on that panel

        implicit none
        
        class(panel_solver),intent(in) :: this
        integer,intent(in) :: i_panel
        type(surface_mesh),intent(in) :: body
        logical,intent(in) :: mirrored
        character(len=*),intent(in) :: rule

        real :: C_P_avg

        C_P_avg = body%panels(i_panel)%get_avg_pressure_coef(body%mu, body%sigma, mirrored, body%N_panels, body%N_verts, &
                                                             body%asym_flow, this%freestream, this%inner_flow, &
                                                             rule, M_corr=this%M_inf_corr)
        
    end function panel_solver_calc_avg_pressure_on_panel


    subroutine panel_solver_calc_pressures(this, body)
        ! Calculates the surface pressures

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: V_pert
        real :: a, b, c, lin, sln

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface pressures..."

        ! Allocate storage
        call this%allocate_pressure_storage(body)

        ! Calculate pressures
        !$OMP parallel do private(V_pert, lin, sln) schedule(static)
        do i=1,body%N_panels

            ! Incompressible rule
            if (this%incompressible_rule) then
                body%C_p_inc(i) = this%calc_avg_pressure_on_panel(i, body, .false., "incompressible")
                if (body%asym_flow) then
                    body%C_p_inc(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "incompressible")
                end if
            end if
        
            ! Isentropic rule
            if (this%isentropic_rule) then
                body%C_p_ise(i) = this%calc_avg_pressure_on_panel(i, body, .false., "isentropic")
                if (body%asym_flow) then
                    body%C_p_ise(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "isentropic")
                end if
            end if

            ! Second-order rule
            if (this%second_order_rule) then
                body%C_p_2nd(i) = this%calc_avg_pressure_on_panel(i, body, .false., "second-order")
                if (body%asym_flow) then
                    body%C_p_2nd(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "second-order")
                end if
            end if
        
            ! Slender-body rule
            if (this%slender_rule) then
                body%C_p_sln(i) = this%calc_avg_pressure_on_panel(i, body, .false., "slender-body")
                if (body%asym_flow) then
                    body%C_p_sln(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "slender-body")
                end if
            end if
        
            ! Linear rule
            if (this%linear_rule) then
                body%C_p_lin(i) = this%calc_avg_pressure_on_panel(i, body, .false., "linear")
                if (body%asym_flow) then
                    body%C_p_lin(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "linear")
                end if
            end if

            ! Prandtl-Glauert correction
            if (this%prandtl_glauert) then
                body%C_p_pg(i) = this%calc_avg_pressure_on_panel(i, body, .false., "prandtl-glauert")
                if (body%asym_flow) then
                    body%C_p_pg(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "prandtl-glauert")
                end if
            end if

            ! Karman-Tsien correction
            if (this%karman_tsien) then
                body%C_p_kt(i) = this%calc_avg_pressure_on_panel(i, body, .false., "karman-tsien")
                if (body%asym_flow) then
                    body%C_p_kt(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "karman-tsien")
                end if
            end if

            ! Laitone correction
            if (this%laitone) then
                body%C_p_lai(i) = this%calc_avg_pressure_on_panel(i, body, .false., "laitone")
                if (body%asym_flow) then
                    body%C_p_lai(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "laitone")
                end if
            end if

        end do

        if (verbose) then
            write(*,*) "Done."
            call this%output_pressures_to_terminal(body)
        end if
        
        ! Check if critical Mach number has been exceeded
        if (this%incompressible_rule) then
            if ((this%M_inf_corr < 1.) .and. (this%M_inf_corr /= 0.)) then
                call this%calc_crit_mach(body)
            end if
        else
            if ((this%freestream%M_inf < 1.) .and. (this%freestream%M_inf /= 0.)) then
                call this%calc_crit_mach(body)
            end if
        end if
        
    end subroutine panel_solver_calc_pressures


    subroutine panel_solver_output_pressures_to_terminal(this, body)
        ! Writes pertinent pressure values to the terminal

        implicit none
        
        class(panel_solver), intent(in) :: this
        type(surface_mesh), intent(in) :: body
        
        ! Report min and max pressure coefficients
        if (this%incompressible_rule) then
            write(*,*) "        Maximum incompressible pressure coefficient:", maxval(body%C_p_inc)
            write(*,*) "        Minimum incompressible pressure coefficient:", minval(body%C_p_inc)
        end if
        
        if (this%isentropic_rule) then
            write(*,*) "        Maximum isentropic pressure coefficient:", maxval(body%C_p_ise)
            write(*,*) "        Minimum isentropic pressure coefficient:", minval(body%C_p_ise)
        end if
        
        if (this%second_order_rule) then
            write(*,*) "        Maximum second-order pressure coefficient:", maxval(body%C_p_2nd)
            write(*,*) "        Minimum second-order pressure coefficient:", minval(body%C_p_2nd)
        end if
        
        if (this%slender_rule) then
            write(*,*) "        Maximum slender-body pressure coefficient:", maxval(body%C_p_sln)
            write(*,*) "        Minimum slender-body pressure coefficient:", minval(body%C_p_sln)
        end if
        
        if (this%linear_rule) then
            write(*,*) "        Maximum linear pressure coefficient:", maxval(body%C_p_lin)
            write(*,*) "        Minimum linear pressure coefficient:", minval(body%C_p_lin)
        end if
        
        ! Report vacuum pressure coefficient
        if (this%freestream%M_inf > 0.) then
            write(*,*) "        Vacuum pressure coefficient:", this%freestream%C_P_vac
        end if
    
    end subroutine panel_solver_output_pressures_to_terminal
    

    subroutine panel_solver_calc_crit_mach(this, body)
        ! Calculates the critical Mach number based on the freestream condition and the pressure solution

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        integer,dimension(1) :: min_loc
        real :: C_p_crit, C_p_min, numerator, denominator, M_inf_selected

        ! Locate minimum pressure location on body based on pressure for forces selection
        select case (this%pressure_for_forces)
            case ("incompressible")
                min_loc = MINLOC(body%C_p_inc)
                C_p_min = MINVAL(body%C_p_inc)
                M_inf_selected = this%M_inf_corr  
            
            case ("prandtl-glauert")
                min_loc = MINLOC(body%C_p_pg)
                C_p_min = MINVAL(body%C_p_pg)
                M_inf_selected = this%M_inf_corr

            case ("karman-tsien")
                min_loc = MINLOC(body%C_p_kt)
                C_p_min = MINVAL(body%C_p_kt)
                M_inf_selected = this%M_inf_corr

            case ("laitone")
                min_loc = MINLOC(body%C_p_lai)
                C_p_min = MINVAL(body%C_p_lai)
                M_inf_selected = this%M_inf_corr

            case ("isentropic")
                min_loc = MINLOC(body%C_p_ise)
                C_p_min = MINVAL(body%C_p_ise)
                M_inf_selected = this%freestream%M_inf

            case ("second-order")
                min_loc = MINLOC(body%C_p_2nd)
                C_p_min = MINVAL(body%C_p_2nd)
                M_inf_selected = this%freestream%M_inf

            case ("linear")
                min_loc = MINLOC(body%C_p_lin)
                C_p_min = MINVAL(body%C_p_lin)
                M_inf_selected = this%freestream%M_inf

            case ("slender-body")
                min_loc = MINLOC(body%C_p_sln)
                C_p_min = MINVAL(body%C_p_sln)
                M_inf_selected = this%freestream%M_inf
                
            case default
                write(*,*) "!!! ", this%pressure_for_forces," pressure for forces is not available. Quitting..."
                stop
        end select

        ! Calculate critical pressure coefficient for the selected Mach number (Modern Compressible Flow by John Anderson EQ 9.55)
        C_p_crit = this%freestream%get_C_P_crit(M_inf_selected)
    
        ! Report the critical mach pressure coefficient and the minimum pressure coefficients to the terminal
        if (verbose) then
            write(*,*) "        Critical Mach C_p = ", C_p_crit
            write(*,*) "        Minimum C_p       = ", C_p_min
        end if

        ! Throw warning message indicating body is in transonic flow
        if ((C_p_min <= C_p_crit) .and. (verbose)) then
            write(*,*) "!!! The critical Mach number has been exceeded over panel ", min_loc, ". Results may not be reliable."
        end if

    end subroutine panel_solver_calc_crit_mach
    

    subroutine panel_solver_calc_forces(this, body)
        ! Calculates the forces
        
        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        
        if (verbose) write(*,'(a, a, a)',advance='no') "     Calculating forces using the ", & 
                                                       this%pressure_for_forces, " pressure rule..."

        ! Calculate forces
        select case (this%pressure_for_forces)

        case ('incompressible')
            call this%calc_forces_with_pressure(body, body%C_p_inc)

        case ('isentropic')
            call this%calc_forces_with_pressure(body, body%C_p_ise)

        case ('second-order')
            call this%calc_forces_with_pressure(body, body%C_p_2nd)

        case ('slender-body')
            call this%calc_forces_with_pressure(body, body%C_p_sln)

        case ('linear')
            call this%calc_forces_with_pressure(body, body%C_p_lin)

        case ('prandtl-glauert')
            call this%calc_forces_with_pressure(body, body%C_p_pg)

        case ('karman-tsien')
            call this%calc_forces_with_pressure(body, body%C_p_kt)

        case ('laitone')
            call this%calc_forces_with_pressure(body, body%C_p_lai)

        end select

        ! Sum discrete forces
        this%C_F(:) = sum(body%dC_f, dim=2)/body%S_ref

        ! Add contributions from mirrored half
        if (body%mirrored .and. .not. body%asym_flow) then
            this%C_F = 2.*this%C_F
            this%C_F(body%mirror_plane) = 0. ! We know this
        end if

        if (verbose) then
            write(*,*) "Done."
            write(*,*) "        Cx:", this%C_F(1)
            write(*,*) "        Cy:", this%C_F(2)
            write(*,*) "        Cz:", this%C_F(3)
        end if
    
    end subroutine panel_solver_calc_forces
    

    subroutine panel_solver_calc_forces_with_pressure(this, body, pressures)
        ! Calculates the forces
        
        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:),allocatable :: pressures
        
        integer :: i, stat

        ! Allocate force storage
        allocate(body%dC_f(3,this%N_cells), stat=stat)
        call check_allocation(stat, "cell forces")

        ! Calculate total forces
        !$OMP parallel do schedule(static)
        do i=1,body%N_panels

            ! Discrete force coefficient acting on panel
            body%dC_f(:,i) = -pressures(i)*body%panels(i)%A*body%panels(i)%n_g

            ! Mirror
            if (body%asym_flow) then
                body%dC_f(:,i+body%N_panels) = -pressures(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
            end if

        end do
    
    end subroutine panel_solver_calc_forces_with_pressure


    function panel_solver_calc_moment_about_centroid_of_panel(this, i_panel, body, mirrored, rule) result(C_M)
        ! Calls the panel function with the needed variables to get the moment coefficient about the panel centroid

        implicit none
        
        class(panel_solver),intent(in) :: this
        integer,intent(in) :: i_panel
        type(surface_mesh),intent(in) :: body
        logical,intent(in) :: mirrored
        character(len=*),intent(in) :: rule

        real,dimension(3) :: C_M

        C_M = body%panels(i_panel)%get_moment_about_centroid(body%mu, body%sigma, mirrored, body%N_panels, body%N_verts, &
                                                             body%asym_flow, this%freestream, this%inner_flow, &
                                                             rule, M_corr=this%M_inf_corr)
        
    end function panel_solver_calc_moment_about_centroid_of_panel


    subroutine panel_solver_calc_moments(this, body)
        ! Calculates the moments acting on the body

        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer i, stat
        real,dimension(:,:),allocatable :: dC_m

        if (verbose) write(*,'(a)',advance='no') "     Calculating moments..."

        ! Allocate moment storage
        allocate(dC_m(3,this%N_cells), stat=stat)
        call check_allocation(stat, "cell moments")
    
        ! Calculate moment induced by each panel
        !$OMP parallel do schedule(static)
        do i=1,body%N_panels

            ! Discrete moment acting on each panel
            dC_m(:,i) = cross(body%panels(i)%centr-body%CG, body%dC_f(:,i))

            ! Add in contribution of varying pressure
            if (body%panels(i)%order == 2) then
                dC_m(:,i) = dC_m(:,i) + this%calc_moment_about_centroid_of_panel(i, body, .false., this%pressure_for_forces)
            end if

            ! Mirror
            if (body%asym_flow) then

                dC_m(:,i+body%N_panels) = cross(body%panels(i)%centr_mir-body%CG, body%dC_f(:,i))

                ! Add in contribution of varying pressure
                if (body%panels(i)%order == 2) then
                    dC_m(:,i+body%N_panels) = dC_m(:,i+body%N_panels) + this%calc_moment_about_centroid_of_panel(i, body, &
                                                                            .true., this%pressure_for_forces)
                end if

            end if
        end do

        ! Sum
        this%C_M = sum(dC_m, dim=2)/body%l_ref

        ! Add contributions from mirrored half
        if (body%mirrored .and. .not. body%asym_flow) then
            do i=1,3
                if (i == body%mirror_plane) then
                    this%C_M(i) = 2.*this%C_M(i)
                else
                    this%C_M(i) = 0.
                end if
            end do
        end if

        if (verbose) then
            write(*,*) "Done."
            write(*,*) "        CMx:", this%C_M(1)
            write(*,*) "        CMy:", this%C_M(2)
            write(*,*) "        CMz:", this%C_M(3)
        end if
        
    end subroutine panel_solver_calc_moments


    subroutine panel_solver_update_report(this, p_json, body, solver_stat)
        ! Updates the report JSON with the information relevant to the solver

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(inout) :: p_json
        type(surface_mesh),intent(in) :: body
        integer,intent(in) :: solver_stat

        type(json_value),pointer :: p_parent, p_child
        real :: max_flow_turning_angle

        ! Write mesh info
        call json_value_create(p_parent)
        call to_object(p_parent, 'mesh_info')
        call json_value_add(p_json, p_parent)
        call json_value_add(p_parent, 'N_body_panels', body%N_panels)
        call json_value_add(p_parent, 'N_body_vertices', body%N_verts)
        call json_value_add(p_parent, 'N_wake_panels', body%wake%N_panels)
        call json_value_add(p_parent, 'average_characteristic_length', body%get_avg_characteristic_panel_length())
        max_flow_turning_angle = acos(body%C_min_panel_angle)*180./pi
        call json_value_add(p_parent, 'max_flow_turning_angle', max_flow_turning_angle)

        ! Write solver results
        call json_value_create(p_parent)
        call to_object(p_parent, 'solver_results')
        call json_value_add(p_json, p_parent)

        ! Solver results
        call json_value_add(p_parent, 'solver_status_code', solver_stat)
        call json_value_add(p_parent, 'system_dimension', this%N_unknown)

        ! Timing
        call json_value_create(p_child)
        call to_object(p_child, 'timing')
        call json_value_add(p_parent, p_child)
        call json_value_add(p_child, 'system_sorting', this%sort_time)
        call json_value_add(p_child, 'preconditioner', this%prec_time)
        call json_value_add(p_child, 'matrix_solver', this%solver_time)
        nullify(p_child)

        ! Check there wasn't an error
        if (solver_stat == 0) then

            ! Iterations
            if (this%solver_iterations > -1) call json_value_add(p_parent, 'iterations', this%solver_iterations)

            ! System lower bandwidth
            if (this%sort_system) call json_value_add(p_parent, "system_lower_bandwidth", this%B_l_system)

            ! Residuals
            call json_value_create(p_child)
            call to_object(p_child, 'residual')
            call json_value_add(p_parent, p_child)
            call json_value_add(p_child, 'max', this%max_res)
            call json_value_add(p_child, 'norm', this%norm_res)

            ! Clean up pointers
            nullify(p_parent)
            nullify(p_child)

            ! Write pressure results
            call json_value_create(p_parent)
            call to_object(p_parent, 'pressure_calculations')
            call json_value_add(p_json, p_parent)

            ! Incompressible rule
            if (this%incompressible_rule) then
                call this%add_pressure_to_report(p_parent, 'incompressible_rule', body%C_p_inc)
            end if

            ! Isentropic rule
            if (this%isentropic_rule) then
                call this%add_pressure_to_report(p_parent, 'isentropic_rule', body%C_p_ise)
            end if

            ! Second-order rule
            if (this%second_order_rule) then
                call this%add_pressure_to_report(p_parent, 'second_order_rule', body%C_p_2nd)
            end if

            ! Slender-body rule
            if (this%slender_rule) then
                call this%add_pressure_to_report(p_parent, 'slender_body_rule', body%C_p_sln)
            end if

            ! Linear rule
            if (this%linear_rule) then
                call this%add_pressure_to_report(p_parent, 'linear_rule', body%C_p_lin)
            end if

            ! Prandtl-Galuert rule
            if (this%prandtl_glauert) then
                call this%add_pressure_to_report(p_parent, 'prandtl_glauert', body%C_p_pg)
            end if

            ! Second-order rule
            if (this%karman_tsien) then
                call this%add_pressure_to_report(p_parent, 'karman_tsien', body%C_p_kt)
            end if

            ! Second-order rule
            if (this%laitone) then
                call this%add_pressure_to_report(p_parent, 'laitone', body%C_p_lai)
            end if
            nullify(p_parent)

            ! Write forces
            call json_value_create(p_parent)
            call to_object(p_parent, 'total_forces')
            call json_value_add(p_json, p_parent)
            call json_value_add(p_parent, 'Cx', this%C_F(1))
            call json_value_add(p_parent, 'Cy', this%C_F(2))
            call json_value_add(p_parent, 'Cz', this%C_F(3))
            nullify(p_parent)

            ! Write moments
            call json_value_create(p_parent)
            call to_object(p_parent, 'total_moments')
            call json_value_add(p_json, p_parent)
            call json_value_add(p_parent, 'CMx', this%C_M(1))
            call json_value_add(p_parent, 'CMy', this%C_M(2))
            call json_value_add(p_parent, 'CMz', this%C_M(3))
            nullify(p_parent)

        end if

    end subroutine panel_solver_update_report


    subroutine panel_solver_add_pressure_to_report(this, p_parent, pressure_label, pressure_values)
        ! Adds the results for a given pressure to the report file

        implicit none
        
        class(panel_solver), intent(in) :: this
        type(json_value),pointer,intent(inout) :: p_parent
        character(len=*),intent(in) :: pressure_label
        real,dimension(:),allocatable,intent(in) :: pressure_values

        type(json_value),pointer :: p_child

        call json_value_create(p_child)
        call to_object(p_child, pressure_label)
        call json_value_add(p_parent, p_child)
        call json_value_add(p_child, 'max', maxval(pressure_values))
        call json_value_add(p_child, 'min', minval(pressure_values))
        nullify(p_child)
        
    end subroutine panel_solver_add_pressure_to_report


    subroutine panel_solver_export_off_body_points(this, points_file, points_output_file, body)
        ! Writes out a csv file of potentials on a user-specified slice of the flow

        implicit none

        class(panel_solver),intent(in) :: this
        character(len=:),allocatable,intent(in) :: points_file, points_output_file
        type(surface_mesh),intent(in) :: body

        integer :: i, unit, N_points, stat
        real :: phi_inf
        real,dimension(3) :: v_inf
        real,dimension(:,:),allocatable :: points
        character(len=200) :: dummy_read
        real,dimension(:),allocatable :: phi_s, phi_d
        real,dimension(:,:),allocatable :: v_s, v_d

        if (verbose) write(*,'(a)',advance='no') "    Calculating potentials at off-body points "

        ! Get number of points from file
        N_points = 0
        open(newunit=unit, file=points_file)

        ! Skip header
        read(unit,*)

        ! Loop through lines
        do
            read(unit,*,iostat=stat) dummy_read

            ! If a line was actually read, increment the number of points
            if (stat == 0) then
                N_points = N_points + 1
            else
                exit ! No more lines
            end if
        end do

        close(unit)

        if (verbose) write(*,'(a, i7, a)',advance='no') "(got ", N_points, " points)..."

        ! Get points
        allocate(points(3,N_points))
        open(newunit=unit, file=points_file)

        ! Skip header
        read(unit,*)

        ! Loop through lines
        do i=1,N_points
            read(unit,*) points(1,i), points(2,i), points(3,i)
        end do

        close(unit)

        ! Allocate potential storage
        allocate(phi_s(N_points), stat=stat)
        call check_allocation(stat, "off-body source potentials")
        allocate(phi_d(N_points), stat=stat)
        call check_allocation(stat, "off-body doublet potentials")

        ! Allocate Velocity storage
        allocate(v_s(3,N_points), stat=stat)
        call check_allocation(stat, "off-body source velocities")
        allocate(v_d(3,N_points), stat=stat)
        call check_allocation(stat, "off-body doublet velocities")

        ! Calculate potentials
        !$OMP parallel do
        do i=1,N_points

            ! Get induced potentials
            call body%get_induced_potentials_at_point(points(:,i), this%freestream, phi_d(i), phi_s(i))
            
            ! Get induced velocities
            call body%get_induced_velocities_at_point(points(:,i), this%freestream, v_d(:,i), v_s(:,i))

        end do

        phi_d = phi_d*this%freestream%U
        phi_s = phi_s*this%freestream%U
        v_d = v_d*this%freestream%U
        v_s = v_s*this%freestream%U
        

        ! Delete old output file
        call delete_file(points_output_file)

        ! Open output file
        100 format(e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, &
                   ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', &
                   e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12)
        open(newunit=unit, file=points_output_file)

        ! Write header
        write(unit,*) 'x,y,z,phi_inf,phi_d,phi_s,phi,Phi,v_inf_x,v_inf_y,v_inf_z,  &
                       v_d_x,v_d_y,v_d_z,v_s_x,v_s_y,v_s_z,v_x,v_y,v_z,V_x,V_y,V_z'

        ! Write potentials out to file
        do i=1,N_points

            ! Calculate freestream potential
            phi_inf = this%freestream%U*inner(points(:,i), this%freestream%c_hat_g)
            ! Calculate freestream velocity
            v_inf = this%freestream%v_inf

            ! Write to file
            write(unit,100) points(1,i), points(2,i), points(3,i), phi_inf, phi_d(i), phi_s(i), phi_d(i) + phi_s(i), &
                            phi_inf + phi_d(i) + phi_s(i), v_inf(1), v_inf(2), v_inf(3), v_d(1,i), v_d(2,i), v_d(3,i), &
                            v_s(1,i), v_s(2,i), v_s(3,i), v_d(1,i) + v_s(1,i), v_d(2,i) + v_s(2,i), v_d(3,i) + v_s(3,i), &
                            v_inf(1) + v_d(1,i) + v_s(1,i), v_inf(2) + v_d(2,i) + v_s(2,i), v_inf(3) + v_d(3,i) + v_s(3,i)

        end do

        close(unit)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_export_off_body_points


end module panel_solver_mod