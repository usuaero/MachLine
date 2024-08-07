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


    ! Actual
    character(len=*),parameter :: D_MORINO = "dirichlet-morino"
    character(len=*),parameter :: D_SOURCE_FREE = "dirichlet-source-free"
    character(len=*),parameter :: N_MF_D_LS = "neumann-mass-flux"
    character(len=*),parameter :: N_V_D_LS = "neumann-velocity"

    ! Experimental
    character(len=*),parameter :: N_MF_D_VCP = "neumann-mass-flux-VCP" !!!! wake version
    character(len=*),parameter :: N_MF_DS_LS = "neumann-doublet-source-mass-flux-ls"
    character(len=*),parameter :: N_MF_D_IF = "neumann-mass-flux-inner-flow"
    character(len=*),parameter :: N_MF_D = "neumann-doublet-only-mass-flux"


    type panel_solver

        real :: M_inf_corr 
        character(len=:),allocatable :: formulation, pressure_for_forces, matrix_solver, preconditioner, iteration_file
        logical :: incompressible_rule, isentropic_rule, second_order_rule, slender_rule, linear_rule
        logical :: write_A_and_b, sort_system, use_sort_for_cp, overdetermined_ls, underdetermined_ls, dirichlet
        logical :: compressible_correction, prandtl_glauert, karman_tsien, laitone
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

        !!!!!! ADJOINT variables !!!!!!
        type(sparse_vector),dimension(:,:),allocatable :: d_A_matrix
        type(sparse_vector),dimension(:),allocatable :: d_b_vector
        type(sparse_matrix) :: d_C_F_wrt_vars, d_C_F_wrt_mu
        ! real,dimension(:),allocatable :: CFx_sensitivities, CFy_sensitivities, CFz_sensitivities
        real,dimension(:,:),allocatable :: CF_sensitivities
        real,dimension(3) :: norms_of_CF_sensitivities

        !!!!! END ADJOINT variables !!!

        contains

            ! Initialization
            procedure :: init => panel_solver_init
            procedure :: parse_solver_settings => panel_solver_parse_solver_settings
            procedure :: parse_processing_settings => panel_solver_parse_processing_settings

            ! Dirichlet setup
            procedure :: init_dirichlet => panel_solver_init_dirichlet
            procedure :: determine_dirichlet_unknowns => panel_solver_determine_dirichlet_unknowns

            ! Neumann setup
            procedure :: init_neumann => panel_solver_init_neumann
            procedure :: determine_neumann_unknowns => panel_solver_determine_neumann_unknowns

            ! Both
            procedure :: set_panel_sources => panel_solver_set_panel_sources

            ! Control points
            procedure :: init_cp_neumann_condition => panel_solver_init_cp_neumann_condition
            procedure :: init_control_point_boundary_conditions => panel_solver_init_control_point_boundary_conditions

            ! Misc setup
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


            !!!!!!!!! ADJOINT PROCEDURES !!!!!!!!!!

            ! update adjoint system A matrix row
            procedure :: update_adjoint_A_row => panel_solver_update_adjoint_A_row
            procedure :: assemble_adjoint_b_vector=>panel_solver_assemble_adjoint_b_vector
            
            procedure :: calc_cell_velocities_adjoint => panel_solver_calc_cell_velocities_adjoint

            procedure :: calc_d_avg_pressure_on_panel_wrt_vars => panel_solver_calc_d_avg_pressure_on_panel_wrt_vars
            procedure :: calc_d_avg_pressure_on_panel_wrt_mu => panel_solver_calc_d_avg_pressure_on_panel_wrt_mu

            ! calc forces
            procedure :: calc_forces_adjoint => panel_solver_calc_forces_adjoint
            procedure :: calc_d_forces_wrt_vars => panel_solver_calc_d_forces_wrt_vars
            procedure :: calc_d_forces_wrt_mu => panel_solver_calc_d_forces_wrt_mu

            ! solve for adjoint v^T terms
            procedure :: solve_for_adjoint_vT_terms => panel_solver_solve_for_adjoint_vT_terms

            ! solve for total derivative (sensitivities) this is the end goal
            procedure :: calc_total_derivative => panel_solver_calc_total_derivative

            !!!!!!!!! END ADJOINT !!!!!!!!!!!!

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
        if (this%formulation == D_MORINO .or. this%formulation == D_SOURCE_FREE) then
            this%dirichlet = .true.
            call this%init_dirichlet(solver_settings, body)
        else if (this%formulation == N_MF_D_LS .or. this%formulation == N_MF_D .or. this%formulation == N_MF_DS_LS &
            .or. this%formulation == N_V_D_LS .or. this%formulation == N_MF_D_IF .or. this%formulation == N_MF_D_VCP) then
                this%dirichlet = .false.
                call this%init_neumann(solver_settings, body)
        else
            write(*,*) "!!! '", this%formulation, "' is not a valid formulation. Quitting..."
            stop
        end if
            
        ! Set boundary conditions
        call this%init_control_point_boundary_conditions(body,freestream)
        
        ! Write out control point geometry
        if (control_point_file /= 'none') then
            call body%write_control_points(control_point_file, solved=.false.)
        end if

        ! Set up permutation for linear system
        call this%set_permutation(body)

    end subroutine panel_solver_init


    subroutine panel_solver_parse_solver_settings(this, solver_settings)
        ! Parses the solver settings from the input

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value),pointer, intent(in) :: solver_settings
        
        ! Get formulation
        call json_xtnsn_get(solver_settings, 'formulation', this%formulation, 'dirichlet-morino')        

        ! Get matrix solver settings
        if (this%formulation == N_MF_D_LS) then
            call json_xtnsn_get(solver_settings, 'matrix_solver', this%matrix_solver, 'HHLS')
            if (this%matrix_solver /= 'HHLS' .and. verbose) write(*,*) "!!! The HHLS solver is recommended for use &
                with the neumann-mass-flux formulation."
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

        ! Special settings for the Neumann formulation
        if (this%formulation == N_MF_D_LS .or. this%formulation == N_V_D_LS) then
            this%sort_system = .false.
            this%use_sort_for_cp = .false.
            this%overdetermined_ls = .true.
            this%underdetermined_ls = .false.
        else if (this%formulation == N_MF_DS_LS) then
            this%sort_system = .false.
            this%use_sort_for_cp = .false.
            this%underdetermined_ls = .true.
            this%overdetermined_ls = .false.
        else if (this%formulation == N_MF_D_VCP) then
            this%sort_system = .false.
            this%use_sort_for_cp = .false.
            this%overdetermined_ls = .false.
            this%underdetermined_ls = .false.
        else 
            this%use_sort_for_cp = .true.
            this%overdetermined_ls = .false.
            this%underdetermined_ls = .false.
        end if

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
        call json_xtnsn_get(solver_settings, 'control_point_offset', offset, 1.0e-7)
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
        call body%place_internal_vertex_control_points(offset, offset_type, this%freestream)
        
        ! Set needed sources
        call this%set_panel_sources(body)

        ! Determine unknowns
        call this%determine_dirichlet_unknowns(body)
        if (this%N_unknown /= body%N_cp) then
            write(*,*) "!!! The number of unknowns is not the same as the number of control points. Quitting..."
            stop
        end if

        ! Calculate target inner flow
        this%inner_flow = this%freestream%c_hat_g
        if (this%formulation == D_SOURCE_FREE) then
            this%inner_flow = this%inner_flow - matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        end if

        if (verbose) write(*,'(a, i6, a)') "Done. Placed", body%N_cp, " control points."
    
    end subroutine panel_solver_init_dirichlet


    subroutine panel_solver_init_neumann(this, solver_settings, body)
        ! Initializes things for a Neumann formulation run

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
        !!!! this should probably be a case at some point
        ! Place control points
        if (this%formulation == N_MF_D_LS) then
            if (verbose) write(*,'(a)',advance='no') "     Placing control points at panel centroids..."
            call body%place_centroid_control_points(.false., 0.)

        else if (this%formulation == N_V_D_LS) then
            if (verbose) write(*,'(a)',advance='no') "     Placing control points above panel centroids..."
            call body%place_centroid_control_points(.false., offset)

        else if (this%formulation == N_MF_D) then
            if (verbose) write(*,'(a)',advance='no') "     Placing control points above select panel centroids..."
            call body%place_sparse_centroid_control_points(0.)
            !call body%place_internal_vertex_control_points(offset, offset_type, this%freestream)
            !call body%place_centroid_vertex_avg_control_points(.false.)

        else if (this%formulation == N_MF_DS_LS) then
            if (verbose) write(*,'(a)',advance='no') "     Placing control points above panel centroids..."
            call body%place_centroid_control_points(.false., offset)

        else if (this%formulation == N_MF_D_IF) then
            if (verbose) write(*,'(a a a ES10.4 a)',advance='no') "     Placing control points control points using a ", &
                                offset_type," offset ratio of ", offset, "..."
            call body%place_internal_vertex_control_points(offset, offset_type, this%freestream)
            
        else if (this%formulation == N_MF_D_VCP) then
            if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Placing control points using a direct offset of ",offset,"..."
            call body%place_vertex_control_points(offset, this%freestream)        
        end if

        ! Sources
        call this%set_panel_sources(body)

        ! Determine unknowns
        call this%determine_neumann_unknowns(body)

        if (verbose) write(*,'(a, i6, a)') "Done. Placed", body%N_cp, " control points."
        
    end subroutine panel_solver_init_neumann


    subroutine panel_solver_set_panel_sources(this, body)
        ! Tells the panels whether they have sources or not

        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i

        ! Loop through panels
        if (this%formulation == D_MORINO .or. this%formulation == N_MF_DS_LS) then
            do i=1,body%N_panels
                body%panels(i)%has_sources = .true.
            end do
        else
            do i=1,body%N_panels
                body%panels(i)%has_sources = body%panels(i)%r < 0.
            end do
        end if
        
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


    subroutine panel_solver_determine_neumann_unknowns(this, body)
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

        if (this%formulation == N_MF_DS_LS) then

            ! All sources are unknown here
            allocate(this%sigma_known(this%N_sigma), source=.false.)
            this%N_s_unknown = this%N_sigma
            ! Allocate vectors mapping unknown sigmas into the (unpermuted) linear system and back
            allocate(this%i_sigma_in_sys(this%N_sigma))
            allocate(this%i_sys_sigma_in_body(this%N_s_unknown))

            ! Create mapping
            do i=1,this%N_sigma
                this%i_sigma_in_sys(i) = i + this%N_d_unknown
                this%i_sys_sigma_in_body(i) = i
            end do

        else

            ! The number of source unknowns is simply the number of superinclined panels
            ! Not yet implemented
            allocate(this%sigma_known(this%N_sigma), source=.true.)
            this%N_s_unknown = 0
        end if

        ! Total number of unknowns
        this%N_unknown = this%N_d_unknown + this%N_s_unknown
        
    end subroutine panel_solver_determine_neumann_unknowns


    subroutine panel_solver_init_cp_neumann_condition(this, cp, freestream, bc_type, body)
        ! Initializes the necessary Neumann condition at the given control point
        ! Handles selecting the proper normal vector for you

        implicit none
        
        class(panel_solver),intent(in) :: this
        type(control_point),intent(inout) :: cp
        type(flow),intent(inout) :: freestream !!!! changed here 
        integer,intent(in) :: bc_type
        type(surface_mesh),intent(inout) :: body
        real, dimension(3) :: average_edge, n_g_new
        integer :: i, i_edge,vert1,vert2,counter,num_edges
        integer, dimension(:), allocatable :: points

        ! wake_normal = cross(freestream%c_hat_g, body%vertices(cp%tied_to_index)%edge) !!!! changed here 
        ! Get vertex normal
        if (cp%tied_to_type == TT_VERTEX .and. cp%cp_type == SURFACE) then
            
            if (body%vertices(cp%tied_to_index)%N_wake_edges>0) then

                ! num_edges = body%vertices(cp%tied_to_index)%adjacent_edges%len()
                ! allocate(points(num_edges))
                ! counter = 0
                ! do i=1,num_edges
                !     ! Check if it's a wake-shedding edge
                !     call body%vertices(cp%tied_to_index)%adjacent_edges%get(i, i_edge)
                !     if (body%edges(i_edge)%sheds_wake) then
                !         !!!! make the average vector here, and then do the cross product
                !         ! get the index of the two points attacehd to the edge 
                !         vert1 = body%edges(i_edge)%top_verts(1)
                !         vert2 = body%edges(i_edge)%top_verts(2)
                !         ! pick which one isn't the current point
                !         if (vert1 == cp%tied_to_index) then
                !             counter = counter + 1
                !             points(counter) = vert2
                !         else
                !             counter = counter + 1
                !             points(counter) = vert1
                !         end if
                !     end if
                ! end do
                ! if (counter == 1) then
                !     ! only use the single wake sheading edge as the vector
                !     average_edge = body%vertices(points(1))%loc - body%vertices(cp%tied_to_index)%loc
                !     ! get the new normal vector
                !     n_g_new = cross(freestream%c_hat_g, average_edge)
                !     ! check to see if this should be an upper or lower normal vector
                !     if (inner(n_g_new,body%vertices(cp%tied_to_index)%n_g_wake) < 0) then
                !         n_g_new = cross(average_edge,freestream%c_hat_g)
                !     end if
                ! else if (counter == 2) then
                !     ! make a vector using the two points
                !     average_edge = body%vertices(points(1))%loc - body%vertices(points(2))%loc
                !     ! get the new normal vector
                !     n_g_new = cross(freestream%c_hat_g, average_edge)
                !     ! check to see if this should be an upper or lower normal vector
                !     if (inner(n_g_new,body%vertices(cp%tied_to_index)%n_g_wake) < 0) then
                !         n_g_new = cross(average_edge,freestream%c_hat_g)
                !     end if
                    
                
                ! else 
                !     write(*,*) "!!!! ERROR wake contains a triple point. These are currently not allowed by this formulation"
                ! end if
                if (body%vertices(cp%tied_to_index)%N_wake_edges==1) then
                    n_g_new = body%vertices(cp%tied_to_index)%n_g
                else
                    n_g_new = body%vertices(cp%tied_to_index)%n_g_wake - freestream%c_hat_g&
                             * inner(freestream%c_hat_g,body%vertices(cp%tied_to_index)%n_g_wake)
                end if 


                n_g_new = n_g_new / norm2(n_g_new) 
                
                if (cp%is_mirror) then
                    call cp%set_bc(bc_type, n_g_new) !!!! this will need to be fixed the rest of the way before mirroring will work
                else
                    call cp%set_bc(bc_type, n_g_new)
                end if
                write(*,*)"vertex_loc",body%vertices(cp%tied_to_index)%loc
                write(*,*)"N_wake_edges",body%vertices(cp%tied_to_index)%N_wake_edges
                write(*,*)"n_g_new",n_g_new
                write(*,*)"n_g_wake",body%vertices(cp%tied_to_index)%n_g_wake
                write(*,*)"n_g",body%vertices(cp%tied_to_index)%n_g



            else 
                if (cp%is_mirror) then
                    call cp%set_bc(bc_type, body%vertices(cp%tied_to_index)%n_g_mir)
                else
                    call cp%set_bc(bc_type, body%vertices(cp%tied_to_index)%n_g)
                    
                    ! if adjoint
                    if (body%calc_adjoint) then
                        call cp%set_bc_adjoint(body%vertices(cp%tied_to_index)%d_n_g)
                    end if

                end if
            end if
        else if (cp%tied_to_type == TT_VERTEX) then
        ! Get panel normal
            if (cp%is_mirror) then
                call cp%set_bc(bc_type, body%vertices(cp%tied_to_index)%n_g_mir)
            else
                call cp%set_bc(bc_type, body%vertices(cp%tied_to_index)%n_g)

            end if
        else
            
            if (cp%is_mirror) then
                call cp%set_bc(bc_type, body%panels(cp%tied_to_index)%n_g_mir)
            else
                call cp%set_bc(bc_type, body%panels(cp%tied_to_index)%n_g)
            end if
        end if        
        
    end subroutine panel_solver_init_cp_neumann_condition


    subroutine panel_solver_init_control_point_boundary_conditions(this, body, freestream)
        ! Sets up the desired boundary conditions on the control points

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body
        type(flow), intent(inout) :: freestream

        integer :: i
    
        ! Loop through control points
        !$OMP parallel do
        do i=1,body%N_cp

            ! If we already know this is a strength-matching control point, move on
            if (body%cp(i)%bc == STRENGTH_MATCHING) cycle

            ! Mirrored control points that are tied to a vertex on the mirror plane whose mirror is not unique are used for strength matching
            if (body%cp(i)%is_mirror .and. body%cp(i)%tied_to_type == TT_VERTEX) then
                if (.not. body%vertices(body%cp(i)%tied_to_index)%mirrored_is_unique) then

                    call body%cp(i)%set_bc(STRENGTH_MATCHING)
                    cycle ! That's all we need to do

                end if
            end if

            select case (this%formulation)

            case (D_MORINO)
                call body%cp(i)%set_bc(ZERO_POTENTIAL)

            case (D_SOURCE_FREE)
                call body%cp(i)%set_bc(SF_POTENTIAL)

            case (N_MF_D_IF)
                call this%init_cp_neumann_condition(body%cp(i), freestream, MF_INNER_FLOW, body)
            case (N_V_D_LS)
                call this%init_cp_neumann_condition(body%cp(i), freestream, ZERO_NORMAL_VEL, body)
            case (N_MF_D_VCP) !!!! this is ours
                call this%init_cp_neumann_condition(body%cp(i), freestream, ZERO_NORMAL_MF, body)
            case default 
                call this%init_cp_neumann_condition(body%cp(i), freestream, ZERO_NORMAL_MF, body)
            
            end select

        end do
        
    end subroutine panel_solver_init_control_point_boundary_conditions


    subroutine panel_solver_set_permutation(this, body)
        ! Creates a vertex/control point permutation which should lead to a more efficient solution of the matrix equation

        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(in) :: body

        real,dimension(:),allocatable :: x
        real,dimension(3) :: loc
        integer :: i, j, k, i_neighbor, i_cp, i_vert, i_panel, i_vert_for_panel, i2, i3, i_panel_abutting, i_opp_edge, source_start
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
            allocate(this%P(body%N_cp))
            do i=1,body%N_cp
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
            source_start = body%N_verts

            ! For a mirrored, asymmetric condition, copy the vertex ordering again
            if (body%asym_flow) then
                this%P(body%N_verts+1:2*body%N_verts) = body%vertex_ordering + body%N_verts
                source_start = 2*body%N_verts
            end if

            ! Copy identity for unknown sources
            do i=1,this%N_s_unknown
                this%P(source_start + i) = source_start + i
            end do

            call system_clock(end_count)

        end if

        ! Calculate elapsed time
        this%sort_time = real(end_count - start_count)/count_rate
    
    end subroutine panel_solver_set_permutation


    subroutine panel_solver_solve(this, body, solver_stat, formulation,freestream) !!!! changed so formulation is in solve
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
        type(flow),intent(inout)::freestream
        integer,intent(out) :: solver_stat
        character(len=:),allocatable,intent(in) :: formulation !!!! changed to get wake influence working 

        integer :: stat
        type(sparse_vector) :: zeros
        real,dimension(:,:),allocatable :: vT_terms
        integer(8) :: start_count, end_count
        real(16) :: count_rate, adjoint_solver_time

        ! Set default status
        solver_stat = 0

        ! Allocate known influence storage
        allocate(this%I_known(body%N_cp), source=0., stat=stat)
        call check_allocation(stat, "known influence vector")
        
        ! Allocate AIC matrix
        allocate(this%A(body%N_cp, this%N_unknown), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")
        
        ! Allocate b vector
        allocate(this%b(body%N_cp), source=0., stat=stat)
        call check_allocation(stat, "b vector")
        
        ! if adjoint
        if (body%calc_adjoint) then
            ! init empty sparse vector
            call zeros%init(body%adjoint_size)
            
            ! allocate AIC sensitivity matrix
            allocate(this%d_A_matrix(body%N_cp, this%N_unknown), source=zeros, stat=stat)
            call check_allocation(stat, "Adjoint AIC sensitivity matrix")
            
            ! allocate b vector 
            allocate(this%d_b_vector(body%N_cp), source=zeros, stat=stat)
            ! call check_allocation(stat, "Adjoint b sensitivity vector")
            
            ! assemble b sensitivity vector
            call this%assemble_adjoint_b_vector(body)
            
            write(*,*) "d_A and d_b allocated"
        end if

        ! Calculate source strengths
        call this%calc_source_strengths(body)
        
        ! Calculate body influences
        call this%calc_body_influences(body)
        
        ! Calculate wake influences
        if (body%wake%N_panels > 0 .or. body%filament_wake%N_filaments > 0)&
        call this%calc_wake_influences(body, formulation,freestream) !!!! formulation part is a change
        !write(*,*) "ROCK YOU LIKE A HURRICANE"
        ! Assemble boundary condition vector
        
        call this%assemble_BC_vector(body)
        
        ! Solve the linear system
        call this%solve_system(body, solver_stat)
        
        ! Check for errors
        if (solver_stat /= 0) return
        
        
        ! Calculate velocities
        call this%calc_cell_velocities(body)
        
        ! if adjoint, calc d_V_inner
        if (body%calc_adjoint) then
            write(*,'(a)',advance='no') "     Calculating adjoint cell velocities..."
            call this%calc_cell_velocities_adjoint(body)
            write(*,*) "Done"
        end if 
        
        ! Calculate potentials
        call this%calc_surface_potentials(body)

        ! Calculate pressures
        call this%calc_pressures(body)

        ! Calculate forces
        call this%calc_forces(body)
        
        ! Calculate moments
        call this%calc_moments(body)
        
        ! if adjoint finish sensitivity calcs
        if (body%calc_adjoint) then
            ! if adjoint, calc d_C_F
            call this%calc_forces_adjoint(body)
            ! call this%calc_moments_adjoint(body)

            call system_clock(start_count, count_rate)

            ! ! if adjoint, solve for v^T terms, then calc total derivatives (sensitivities)
            vT_terms = this%solve_for_adjoint_vT_terms(body)
            
            ! calc the total derivative. this is the end game!
            call this%calc_total_derivative(body, vT_terms)

            call system_clock(end_count)

            adjoint_solver_time = real(end_count-start_count)/count_rate
            write(*,'(A, f10.5, A)') "adjoint solver time = ", adjoint_solver_time, " seconds"
        end if

    end subroutine panel_solver_solve


    subroutine panel_solver_assemble_BC_vector(this, body)
        ! Sets up the {BC} vector used in the linear system synthesis

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(in) :: body

        integer :: i, ind
        real,dimension(3) :: x

        ! Allocate
        allocate(this%BC(body%N_cp))

        ! Vector for source-free target potential calculation
        x = matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)

        ! Loop through control points
        do i=1,body%N_cp

            ! Get index
            if (this%use_sort_for_cp) then
                ind = this%P(i)
            else
                ind = i
            end if
            
            ! Check boundary condition type
            select case (body%cp(i)%bc)

            ! Source-free formulation
            case (SF_POTENTIAL)
                this%BC(ind) = -inner(x, body%cp(i)%loc)
            ! Neumann (mass flux)
            case (ZERO_NORMAL_MF)
                this%BC(ind) = -inner(body%cp(i)%n_g,this%freestream%c_hat_g)
            ! Neumann (mass-flux-inner-flow)
            case (MF_INNER_FLOW)
                this%BC(ind) = -inner(x, body%cp(i)%loc)

            ! Neumann (velocity)
            case (ZERO_NORMAL_VEL)
                this%BC(ind) = -inner(this%freestream%c_hat_g, body%cp(i)%n_g)

            ! All other cases
            ! 1: For Morino, desired BC is zero inner potential
            ! 4: For strength-matching, desired BC is zero difference between strengths
            case default
                this%BC(ind) = 0.

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
        if (this%formulation == D_MORINO) then

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

        ! Add source influence depending whether the panel has sources
        if (body%panels(i_panel)%has_sources) then

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

        integer :: i, j, k
        real,dimension(:),allocatable :: source_inf, doublet_inf
        real,dimension(:,:),allocatable :: v_s, v_d
        real,dimension(this%N_unknown) :: A_i
        real :: I_known_i

        type(sparse_vector), dimension(3,3) ::  d_v_d
        type(sparse_vector), dimension(:),allocatable :: inf_adjoint
        type(sparse_vector) :: zeros
        type(sparse_vector), dimension(this%N_unknown) :: d_AIC_row

        ! if adjoint, set up zero sparse vector for initialization
        if (body%calc_adjoint) then
            call zeros%init(body%adjoint_size)
        end if
        
        if (body%calc_adjoint) then
            if (verbose) write(*,'(a)',advance='no') "     Calculating body influences and adjoint body influences..."
        else 
            if (verbose) write(*,'(a)',advance='no') "     Calculating body influences..."
        end if

        ! Calculate source and doublet influences from body on each control point
        !$OMP parallel do private(j, source_inf, doublet_inf, v_s, v_d, A_i, I_known_i, inf_adjoint,d_AIC_row) schedule(dynamic)
        do i=1,body%N_cp

            ! Initialize
            A_i = 0.
            I_known_i = 0.

            ! if adjoint, set d_AIC_row elements to zero vectors
            if (body%calc_adjoint) then
                d_AIC_row = zeros
            end if

            ! Determine the type of boundary condition on this control point
            select case (body%cp(i)%bc)

            case (STRENGTH_MATCHING) ! Strength matching

                A_i(this%P(i)) = 1.
                A_i(this%P(i-body%N_cp/2)) = -1.

            case (ZERO_NORMAL_MF) ! Calculate normal mass flux influences
                ! Loop through panels
                
                do j=1,body%N_panels
                    ! Influence of existing panel on control point
                    call body%panels(j)%calc_velocity_influences(body%cp(i)%loc, this%freestream,.false.,v_s, v_d)
                    
                    !!!!!!!!!!!! ADJOINT CALCS HAPPEN HERE !!!!!!!!!!!!!!!!!!!!

                    ! ! if calc_adjoint is specified, do adjoint calcs (METHOD 1)
                    ! if (body%calc_adjoint) then
                    !     call body%panels(j)%calc_velocity_influences_adjoint(body%cp(i)%loc, body%cp(i)%d_loc, &
                    !         this%freestream, d_v_d)
                    !     inf_adjoint = body%panels(j)%calc_doublet_inf_adjoint(body%cp(i), this%freestream, d_v_d)
                    !     call this%update_adjoint_A_row(body, body%cp(i), d_AIC_row, j, inf_adjoint, .false.)
                    ! end if

                    ! if calc_adjoint is specified, do adjoint calcs (METHOD 2)
                    if (body%calc_adjoint) then
                        inf_adjoint = body%panels(j)%calc_doublet_inf_adjoint2(body%cp(i),this%freestream)
                        call this%update_adjoint_A_row(body, body%cp(i), d_AIC_row, j, inf_adjoint, .false.)
                    end if

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    source_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_s))
                    doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))
                    
                    ! Add influence
                    call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, .false.)         
 
                    if (body%mirrored) then

                        ! Calculate influence of mirrored panel on control point
                        call body%panels(j)%calc_velocity_influences(body%cp(i)%loc, this%freestream, .true., v_s, v_d)
                        source_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_s))
                        doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))

                        ! Add influence
                        call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, body%asym_flow)

                    end if

                end do
                !Write(*,*) body%cp(i)%n_g , body%cp(i)%loc
            case (MF_INNER_FLOW) ! Calculate inner flow influences

                ! Loop through panels
                do j=1,body%N_panels

                    ! Calculate influence of existing panel on control point
                    call body%panels(j)%calc_velocity_influences(body%cp(i)%loc, this%freestream, .false., v_s, v_d)
                    source_inf = matmul(body%cp(i)%loc, v_s)
                    doublet_inf = matmul(body%cp(i)%loc, v_d)
                    
                    ! Add influence
                    call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, .false.)

                    if (body%mirrored) then

                        ! Calculate influence of mirrored panel on control point
                        call body%panels(j)%calc_velocity_influences(body%cp(i)%loc,this%freestream, .true.,v_s, v_d)
                        source_inf = matmul(body%cp(i)%loc, v_s)
                        doublet_inf = matmul(body%cp(i)%loc, v_d)

                        ! Add influence
                        call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, body%asym_flow)

                    end if

                end do

            case (ZERO_NORMAL_VEL) ! Calculate normal velocity influences

                ! Loop through panels
                do j=1,body%N_panels

                    ! Influence of existing panel on control point
                    call body%panels(j)%calc_velocity_influences(body%cp(i)%loc, this%freestream, .false., v_s, v_d)
                    source_inf = matmul(body%cp(i)%n_g, v_s)
                    doublet_inf = matmul(body%cp(i)%n_g, v_d)

                    ! Add influence
                    call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, .false.)
                    if (body%mirrored) then

                        ! Calculate influence of mirrored panel on control point
                        call body%panels(j)%calc_velocity_influences(body%cp(i)%loc, this%freestream, .true., v_s, v_d)
                        source_inf = matmul(body%cp(i)%n_g, v_s)
                        doublet_inf = matmul(body%cp(i)%n_g, v_d)

                        ! Add influence
                        call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, body%asym_flow)

                    end if

                end do

            case default ! Calculate potential influences
                ! Loop through panels
                do j=1,body%N_panels

                    ! Influence of existing panel on control point
                    call body%panels(j)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
                                                                  .false., source_inf, doublet_inf)

                    !!!!!!!!!!!! ADJOINT CALCS HAPPEN HERE !!!!!!!!!!!!!!!!!!!!

                    ! if calc_adjoint is specified, do adjoint calcs (METHOD 2)
                    if (body%calc_adjoint) then
                        inf_adjoint = body%panels(j)%calc_adjoint_potential_influences(body%cp(i),this%freestream, .false.)
                        call this%update_adjoint_A_row(body, body%cp(i), d_AIC_row, j, inf_adjoint, .false.)
                    end if
                                            
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    ! Add influence
                    call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, .false.)

                    if (body%mirrored) then

                        ! Calculate influence of mirrored panel on control point
                        call body%panels(j)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
                                                                      .true., source_inf, doublet_inf)

                        ! Add influence
                        call this%update_system_row(body, body%cp(i), A_i, I_known_i, j, source_inf, doublet_inf, body%asym_flow)

                    end if

                end do

            end select

            ! Update system
            !$OMP critical

            ! Update A matrix and I_known
            if (this%use_sort_for_cp) then
                this%A(this%P(i),:) = A_i
                this%I_known(this%P(i)) = I_known_i

                ! if adjoint, assemble d_A_matrix (sort used for Dirichlet adjoint)
                if (body%calc_adjoint) then
                    this%d_A_matrix(this%P(i),:) = d_AIC_row
                end if
                
            else
                this%A(i,:) = A_i
                this%I_known(i) = I_known_i

                ! if adjoint, assemble d_A_matrix (sort not used for Neumann adjoint)
                if (body%calc_adjoint) then
                    this%d_A_matrix(i,:) = d_AIC_row
                end if
            end if
            !$OMP end critical

        end do
        if (verbose) write(*,*) "Done."
    
    end subroutine panel_solver_calc_body_influences


    subroutine panel_solver_calc_wake_influences(this, body, formulation,freestream)
        ! Calculates the influence of the wake on the control points

        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(inout)::freestream
        character(len=:),allocatable,intent(in) :: formulation
        integer :: i, j, k, l, m
        real,dimension(:),allocatable ::  doublet_inf, source_inf
        real,dimension(this%N_unknown) :: A_i
        real,dimension(:,:),allocatable :: v_s, v_d
        real :: s_star,x !!!! might remove
        real,dimension(3) :: d_from_te
        logical :: downstream,passes_through

        type(sparse_vector), dimension(:),allocatable :: inf_adjoint
        type(sparse_vector) :: zeros
        type(sparse_vector), dimension(this%N_unknown) :: d_AIC_row

        ! if adjoint, set up zero sparse vector for initialization
        if (body%calc_adjoint) then
            call zeros%init(body%adjoint_size)
        end if
        
        if (body%calc_adjoint) then
            if (verbose) write(*,'(a)',advance='no') "     Calculating wake influences and adjoint wake influences..."
        else 
            if (verbose) write(*,'(a)',advance='no') "     Calculating wake influences..."
        end if
        ! Calculate influence of wake

        ! Loop through control points
        !$OMP parallel do private(j, k, l, source_inf, doublet_inf, v_d, v_s, A_i, s_star)&
        !$OMP & private(passes_through,downstream,d_from_te,x) schedule(dynamic)
        do i=1,body%N_cp

            ! if adjoint, set d_AIC_row elements to zero vectors
            if (body%calc_adjoint) then
                d_AIC_row = zeros
            end if

            ! Check boundary condition
            select case (body%cp(i)%bc)

            case (STRENGTH_MATCHING) ! Strength matching
                cycle

            case (ZERO_NORMAL_MF) ! Calculate normal mass flux influences !!!! add our stuff with logic to choose filaments or panels
                if (body%wake_has_filaments(formulation)) then !!!! start wake_dev, might need to adjust cases instead of doing an if statement like this - SA
                    ! Initialize
                    !!!! need to see if we need to ignore a vertex here -jjh
                    A_i = 0.
                
                    ! Get doublet influence from wake filaments 
                    do j=1,body%filament_wake%N_filaments 
                        do l=1,body%filament_wake%filaments(j)%N_segments
                            
                            ! check if the filament passes through the panel of the control point
                            !! first check if it is downstream
                            d_from_te = body%cp(i)%loc - body%filament_wake%filaments(j)%segments(l)&
                                                                                %vertices(1)%ptr%loc
                            x = inner(d_from_te, freestream%c_hat_g)    
                            downstream = x > 0.     
                            ! write(*,*) "c_hat_g",freestream%c_hat_g                                         
                            ! write(*,*) "d_from_te",d_from_te
                            ! write(*,*) "X",X
                            ! write(*,*) "DS?",downstream

                            !! now check if the filament passes through 
                            passes_through = body%panels(i)%filament_passes_through(body%filament_wake%filaments(j)&
                            %segments(l)%vertices(1)%ptr%loc, &
                            body%filament_wake%filaments(j)%segments(l)%vertices(2)%ptr%loc, &
                            .false., s_star)
                            if (passes_through.and.downstream .and. body%cp(i)%tied_to_type==2)  then !!!! tied to type=2 means the cp is on a panel
                                doublet_inf = (/0.0, 0.0, 0.0, 0.0/) !!!! is this actually a scalar? 
                                write(*,*) "Panel:", body%cp(i)%tied_to_index-1,"ignores filament segment", &
                                ((j-1)*body%filament_wake%filaments(j)%N_segments)+&
                                body%filament_wake%filaments(j)%segments(l)%index-1, "due to impingement"
                            else
                                call body%filament_wake%filaments(j)%segments(l)%calc_velocity_influences(& 
                                body%cp(i)%loc, this%freestream,.false., v_d) 
                                doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))
                            end if

                            !write(*,*) "doublet_inf", doublet_inf
                            !write(*,*) ""
                            ! call body%filament_wake%filaments(j)%segments(l)%calc_velocity_influences(& 
                            !                                     body%cp(i)%loc, this%freestream,.false., v_d) 
                            ! doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d)) 
                
                            ! Add influence
                            do k=1,size(body%filament_wake%filaments(j)%segments(l)%parents) !!!! will need to change i_vert_d. Look at panel_set_doublet_verts for where that's set
                                A_i(this%P(body%filament_wake%filaments(j)%segments(l)%parents(k))) = &
                                    A_i(this%P(body%filament_wake%filaments(j)%segments(l)%parents(k))) + doublet_inf(k)
                            end do
                
                            ! Get influence of mirrored panel
                            if (body%filament_wake%filaments(j)%mirrored) then
                                    ! check if the filament passes through the panel of the control point
                                if (body%panels(i)%filament_passes_through(body%filament_wake%filaments(j)&
                                    %segments(l)%vertices(1)%ptr%loc, &
                                    body%filament_wake%filaments(j)%segments(l)%vertices(2)%ptr%loc, &
                                    .true., s_star)) then ! vertices list might be called vertex
                                    doublet_inf = (/0.0, 0.0, 0.0, 0.0/) !!!! is this actually a scalar? 
                                else
                                    ! Calculate influence of mirrored panel on control point
                                    call body%filament_wake%filaments(j)%segments(l)%calc_velocity_influences(& 
                                    body%cp(i)%loc, this%freestream,.false., v_d) 
                                    doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))      
                                end if
                                ! Calculate influence of mirrored panel on control point
                                ! call body%filament_wake%filaments(j)%segments(l)%calc_velocity_influences(& 
                                !                                      body%cp(i)%loc, this%freestream,.true., v_d)
                                ! doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))
                
                                ! Add influence
                                do k=1,size(body%filament_wake%filaments(j)%segments(l)%parents) !!!! will need to change i_vert_d. Look at panel_set_doublet_verts for where that's set
                                    A_i(this%P(body%filament_wake%filaments(j)%segments(l)%parents(k))) = &
                                              A_i(this%P(body%filament_wake%filaments(j)%segments(l)%parents(k))) &
                                               + doublet_inf(k)
                                end do
                
                            end if
                        end do
                    end do
                else 
                    !!!! if it isn't a filament wake, then do panel wake solver stuff. - SA
                    ! Initialize
                    A_i = 0.

                    ! Get doublet influence from wake strips
                    do j=1,body%wake%N_strips
                        do l=1,body%wake%strips(j)%N_panels
                            ! Caclulate influence of existing panel on control point
                            call body%wake%strips(j)%panels(l)%calc_velocity_influences(body%cp(i)%loc, this%freestream, &
                                                                                        .false., v_s, v_d)
                            doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))

                            ! Add influence
                            do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d) 
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                    A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                            end do

                            ! Get influence of mirrored panel
                            if (body%wake%strips(j)%mirrored) then

                                ! Calculate influence of mirrored panel on control point
                                call body%wake%strips(j)%panels(l)%calc_velocity_influences(body%cp(i)%loc, this%freestream, &
                                                                        .true., v_s, v_d)
                                doublet_inf = matmul(body%cp(i)%n_g, matmul(this%freestream%B_mat_g, v_d))

                                ! Add influence
                                do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                                    A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                                end do

                            end if
                        end do
                    end do
                end if 
            case (MF_INNER_FLOW) ! Calculate inner flow wake influences

                ! Initialize
                A_i = 0.

                ! Get doublet influence from wake strips
                do j=1,body%wake%N_strips
                    do l=1,body%wake%strips(j)%N_panels

                        ! Caclulate influence of existing panel on control point
                        call body%wake%strips(j)%panels(l)%calc_velocity_influences(body%cp(i)%loc, this%freestream, &
                                                                                     .false., v_s, v_d)
                        doublet_inf = matmul(body%cp(i)%loc, v_d)

                        ! Add influence
                        do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                        end do

                        ! Get influence of mirrored panel
                        if (body%wake%strips(j)%mirrored) then

                            ! Calculate influence of mirrored panel on control point
                            call body%wake%strips(j)%panels(l)%calc_velocity_influences(body%cp(i)%loc, this%freestream, &
                                                                     .true.,  v_s, v_d)
                            doublet_inf = matmul(body%cp(i)%loc, v_d)

                            ! Add influence
                            do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                            end do

                        end if
                    end do
                end do
                
            case (ZERO_NORMAL_VEL) ! Calculate velocity flux influences !!!! add our stuff

                ! Initialize
                A_i = 0.

                ! Get doublet influence from wake strips
                do j=1,body%wake%N_strips
                    do l=1,body%wake%strips(j)%N_panels

                        ! Caclulate influence of existing panel on control point
                        call body%wake%strips(j)%panels(l)%calc_velocity_influences(body%cp(i)%loc, this%freestream, &
                                                                                     .false.,  v_s, v_d)
                        doublet_inf = matmul(body%cp(i)%n_g, v_d)

                        ! Add influence
                        do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                        end do

                        ! Get influence of mirrored panel
                        if (body%wake%strips(j)%mirrored) then

                            ! Calculate influence of mirrored panel on control point
                            call body%wake%strips(j)%panels(l)%calc_velocity_influences(body%cp(i)%loc, this%freestream, &
                                                                     .true.,  v_s, v_d)
                            doublet_inf = matmul(body%cp(i)%n_g, v_d)

                            ! Add influence
                            do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                            end do

                        end if
                    end do
                end do

            case default ! Calculate potential influences

                ! Initialize
                A_i = 0.

                ! Get doublet influence from wake strips
                do j=1,body%wake%N_strips
                    do l=1,body%wake%strips(j)%N_panels

                        ! Caclulate influence of existing panel on control point
                        call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
                                                                                     .false., source_inf, doublet_inf)

                        !!!!!!!!!!!! ADJOINT CALCS HAPPEN HERE !!!!!!!!!!!!!!!!!!!!
                        ! if calc_adjoint is specified, do adjoint calcs (METHOD 2)
                        if (body%calc_adjoint) then
                            inf_adjoint = body%wake%strips(j)%panels(l)%calc_adjoint_potential_influences(body%cp(i),&
                            this%freestream, .false.)
                            ! call this%update_adjoint_A_row(body, body%cp(i), d_AIC_row, j, inf_adjoint, .false.)
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                        ! Add influence
                        do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                            
                            ! add inf_adjoint(k) to the appropriate sparse_vector in d_AIC_row
                            if (body%calc_adjoint) then
                                call d_AIC_row(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k)))%sparse_add(&
                                inf_adjoint(k))
                            end if

                        end do

                        ! Get influence of mirrored panel
                        if (body%wake%strips(j)%mirrored) then

                            ! Calculate influence of mirrored panel on control point
                            call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp(i)%loc, this%freestream, &
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
            if (this%use_sort_for_cp) then !!!! is this the same in filaments as it is in panels? 
                this%A(this%P(i),:) = this%A(this%P(i),:) + A_i

                ! if adjoint, assemble d_A_matrix (sort used for Dirichlet adjoint)
                if (body%calc_adjoint) then
                    do m = 1, body%N_cp
                        call this%d_A_matrix(this%P(i),m)%sparse_add(d_AIC_row(m))
                    end do
                end if

            else
                this%A(i,:) = this%A(i,:) + A_i

                ! if adjoint, assemble d_A_matrix (sort not used for Neumann adjoint)
                if (body%calc_adjoint) then
                    do m = 1, body%N_cp
                        call this%d_A_matrix(i,m)%sparse_add(d_AIC_row(m))
                    end do
                end if

            end if
            !$OMP end critical

        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_wake_influences


    subroutine panel_solver_check_system(this, solver_stat, body)
        ! Checks the validity of the linear system

        implicit none
        
        class(panel_solver),intent(in) :: this
        integer,intent(inout) :: solver_stat
        type(surface_mesh),intent(in) :: body

        integer :: i

        if (verbose) write(*,'(a)',advance='no') "     Checking validity of linear system..."

        ! Check for NaNs
        if (any(isnan(this%A))) then
            write(*,*) "!!! Invalid value detected in A matrix."
            solver_stat = 1
        end if
        if (any(isnan(this%b))) then
            write(*,*) "!!! Invalid value detected in b vector."
            solver_stat = 1
        end if

        if (solver_stat /= 0) return

        ! Check for uninfluenced/ing points
        if (this%use_sort_for_cp) then
            do i=1,body%N_cp
                if (all(this%A(this%P(i),:) == 0.)) then
                    write(*,*) "!!! Control point ", i, " is not influenced."
                    solver_stat = 2
                end if
            end do
            do i=1,this%N_unknown
                if (all(this%A(:,this%P(i)) == 0.)) then
                    write(*,*) "!!! Vertex ", i, " exerts no influence."
                    solver_stat = 2
                end if
            end do
        else
            do i=1,body%N_cp
                if (all(this%A(i,:) == 0.)) then
                    write(*,*) "!!! Control point ", i, " is not influenced."
                    solver_stat = 2
                end if
            end do
            do i=1,this%N_unknown
                if (all(this%A(:,i) == 0.)) then
                    write(*,*) "!!! Vertex ", i, " exerts no influence."
                    solver_stat = 2
                end if
            end do
        end if

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_check_system


    subroutine panel_solver_write_system(this, body)
        ! Writes the linear system to file

        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(in) :: body

        integer :: i, j, unit

        if (verbose) write(*,'(a)',advance='no') "     Writing linear system to file..."

        ! A
        open(newunit=unit, file="A_mat.txt")
        do i=1,body%N_cp
            do j=1,this%N_unknown
                write(unit,'(e20.12)',advance='no') this%A(i,j)
            end do
            write(unit,*)
        end do
        close(unit)

        ! b
        open(newunit=unit, file="b_vec.txt")
        do i=1,body%N_cp
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

        real,dimension(:,:),allocatable :: A_p, A_T
        real,dimension(:),allocatable :: b_p, x, x_temp
        integer :: stat, i, N
        integer(8) :: start_count, end_count
        real(16) :: count_rate

        ! Set b vector
        this%b = this%BC - this%I_known

        !! Force arbitrary doublet constant
        !if (this%formulation == N_MF_DS_LS) then
        !    this%A(body%N_cp,:) = 0.
        !    this%A(body%N_cp,this%N_d_unknown) = 1.
        !    this%b(body%N_cp) = 0.
        !end if

        ! Run checks
        if (run_checks) call this%check_system(solver_stat, body)

        ! Calculate lower bandwidth
        if (this%sort_system) this%B_l_system = get_lower_bandwidth(this%N_unknown, this%A)
        
        !!!! scale by area 
                ! do i = 1,body%N_panels
        !     this%A(i,:) = this%A(i,:) * body%panels(i)%A
        !     this%b(i) = this%b(i) * body%panels(i)%A
        ! end do

        ! Write to file
        if (this%write_A_and_b) call this%write_system(body)

        if (verbose) write(*,'(a, a, a)',advance='no') "     Solving linear system (method: ", this%matrix_solver, ")..."

        ! Precondition
        select case(this%preconditioner)

        ! Diagonal preconditioning
        case ('DIAG')

            ! Overdetermined least-squares
            if (this%overdetermined_ls) then
                if (this%matrix_solver == 'HHLS') then
                    call system_clock(start_count, count_rate)
                    allocate(A_p, source=this%A, stat=stat)
                    call check_allocation(stat, "solver copy of AIC matrix")
                    allocate(b_p, source=this%b, stat=stat)
                    call check_allocation(stat, "solver copy of b vector")
                    call system_clock(end_count)
                else
                    call system_clock(start_count, count_rate)
                    call diagonal_preconditioner(this%N_unknown, matmul(transpose(this%A), this%A), &
                                                 matmul(transpose(this%A), this%b), A_p, b_p)
                    call system_clock(end_count)
                end if

            ! Underdetermined least-squares
            else if (this%underdetermined_ls) then
                call system_clock(start_count, count_rate)
                !call diagonal_preconditioner(body%N_cp, matmul(this%A, transpose(this%A)), this%b, A_p, b_p)
                A_T = transpose(this%A)
                call diagonal_preconditioner(body%N_cp, parallel_matmul(this%A, A_T), this%b, A_p, b_p)
                call system_clock(end_count)

            ! Standard
            else
                call system_clock(start_count, count_rate)
                call diagonal_preconditioner(this%N_unknown, this%A, this%b, A_p, b_p)
                call system_clock(end_count)
            end if

        ! No preconditioning
        case default

            ! Overdetermined least-squares
            if (this%overdetermined_ls) then
                if (this%matrix_solver == 'HHLS') then
                    call system_clock(start_count, count_rate)
                    allocate(A_p, source=this%A, stat=stat)
                    call check_allocation(stat, "solver copy of AIC matrix")
                    allocate(b_p, source=this%b, stat=stat)
                    call check_allocation(stat, "solver copy of b vector")
                    call system_clock(end_count)
                else
                    call system_clock(start_count, count_rate)
                    allocate(A_p, source=matmul(transpose(this%A), this%A), stat=stat)
                    call check_allocation(stat, "solver copy of AIC matrix")
                    allocate(b_p, source=matmul(transpose(this%A), this%b), stat=stat)
                    call check_allocation(stat, "solver copy of b vector")
                    call system_clock(end_count)
                end if

            ! Underdetermined least-squares
            else if (this%underdetermined_ls) then
                call system_clock(start_count, count_rate)
                allocate(A_p, source=matmul(this%A, transpose(this%A)), stat=stat)
                call check_allocation(stat, "solver copy of AIC matrix")
                allocate(b_p, source=this%b, stat=stat)
                call check_allocation(stat, "solver copy of b vector")
                call system_clock(end_count)

            ! Standard
            else
                call system_clock(start_count, count_rate)
                allocate(A_p, source=this%A, stat=stat)
                call check_allocation(stat, "solver copy of AIC matrix")
                allocate(b_p, source=this%b, stat=stat)
                call check_allocation(stat, "solver copy of b vector")
                call system_clock(end_count)
            end if

        end select

        ! Determine size of system to be solved
        if (this%underdetermined_ls) then
            N = body%N_cp
        else
            N = this%N_unknown
        end if

        ! Calculate how much time preconditioning took
        this%prec_time = real(end_count-start_count)/count_rate

        ! Check block size
        if (this%block_size <= 0) then
            this%block_size = N / 5
        end if

        ! Solve
        select case(this%matrix_solver)

        ! LU decomposition
        case ('LU')
            call system_clock(start_count, count_rate)
            call lu_solve(N, A_p, b_p, x)
            call system_clock(end_count)

        ! QR via Givens rotations for upper-pentagonal
        case ('QRUP')
            call system_clock(start_count, count_rate)
            call QR_givens_solve_UP(N, A_p, b_p, x)
            call system_clock(end_count)

        ! QR via fast Givens rotations for upper-pentagonal
        case ('FQRUP')
            call system_clock(start_count, count_rate)
            call QR_fast_givens_solve_upper_pentagonal(N, A_p, b_p, x)
            call system_clock(end_count)

        ! GMRES
        case ('GMRES')
            call system_clock(start_count, count_rate)
            call GMRES(N, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        ! Restarted GMRES
        case ('RGMRES')
            call system_clock(start_count, count_rate)
            call restarted_GMRES(N, A_p, b_p, this%tol, this%max_iterations, this%restart_iterations, &
                                 this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        ! Purcell's method
        case ('PURC')
            call system_clock(start_count, count_rate)
            call purcell_solve(N, A_p, b_p, x)
            call system_clock(end_count)

        ! Block symmetric successive over-relaxation
        case ('BSSOR')
            call system_clock(start_count, count_rate)
            call block_ssor_solve(N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                  this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)
        
        ! Block Jacobi
        case ('BJAC')
            call system_clock(start_count, count_rate)
            call block_jacobi_solve(N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                    this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)
        
        ! Household least-squares
        case ('HHLS')
            call system_clock(start_count, count_rate)
            call householder_ls_solve(body%N_cp, N, A_p, b_p, x)
            call system_clock(end_count)

        ! Improper specification
        case default
            write(*,*) "!!! ", this%matrix_solver, " is not a valid solver option. Defaulting to GMRES."
            call system_clock(start_count, count_rate)
            call GMRES(N, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, this%solver_iterations, x)
            call system_clock(end_count)

        end select
        if (verbose) write(*,*) "Done."

        ! Calculate how much time the matrix solver took
        this%solver_time = real(end_count-start_count)/count_rate

        ! Clean up memory
        deallocate(A_p)
        deallocate(b_p)

        ! Parse out underdetermined solution
        if (this%underdetermined_ls) then
            x_temp = matmul(transpose(this%A), x)
            call move_alloc(x_temp, x)
        end if

        ! Get residual vector
        body%R_cp = matmul(this%A, x) - this%b
        
        if (body%calc_adjoint) then
            ! don't deallocate
        else
            ! write(*,*) " A matrix"
            ! do i = 1,N 
            !     write(*,*) this%A(i,:)
            ! end do
            
            ! write(*,*) " b vector"
            ! do i = 1,N 
            !     write(*,*) this%b(i)
            ! end do
            deallocate(this%A)
        end if 

        deallocate(this%b)

        ! Calculate residual parameters
        this%max_res = maxval(abs(body%R_cp))
        this%norm_res = sqrt(sum(body%R_cp*body%R_cp))
        if (verbose) then
            write(*,*) "        Maximum residual:", this%max_res
            write(*,*) "        Norm of residual:", this%norm_res
            if (this%sort_system) write(*,*) "        Lower bandwidth of A matrix:", this%B_l_system
        end if
        
        !!! TESTING !!! 
        write(*,*) "        Maximum residual:", this%max_res
        write(*,*) "        Norm of residual:", this%norm_res
        !!!!!!!

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

        integer :: i, stat
        real,dimension(3) :: v_inner, v_d, v_s, P

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface velocities..."

        ! Determine number of surface cells
        this%N_cells = body%N_panels
        if (body%asym_flow) then
            this%N_cells = this%N_cells*2
        end if

        ! Allocate cell velocity storage
        allocate(body%V_cells_inner(3,this%N_cells), stat=stat)
        call check_allocation(stat, "inner velocity vectors")
        allocate(body%V_cells(3,this%N_cells), stat=stat)
        call check_allocation(stat, "surface velocity vectors")

        ! Get the surface velocity on each existing panel
        !$OMP parallel do private(v_d, v_s, P) schedule(static)
        do i=1,body%N_panels

            ! Get inner flow
            if (this%dirichlet) then
                body%V_cells_inner(:,i) = this%inner_flow*this%freestream%U
            else
                P = body%panels(i)%centr - 1.e-10*body%panels(i)%n_g
                call body%get_induced_velocities_at_point(P, this%freestream, v_d, v_s)
                body%V_cells_inner(:,i) = this%freestream%v_inf + this%freestream%U*(v_d + v_s)
                ! body%V_cells_inner(:,i) = 0
            end if

            ! Get surface velocity on each panel
            body%V_cells(:,i) = body%panels(i)%get_velocity(body%mu, body%sigma, .false., body%N_panels, &
                                                            body%N_verts, body%asym_flow, this%freestream, &
                                                            body%V_cells_inner(:,i)/this%freestream%U)

            ! Get surface velocity on each mirrored panel
            if (body%asym_flow) then

                ! Get inner flow
                if (this%dirichlet) then
                    body%V_cells_inner(:,i+body%N_panels) = this%inner_flow*this%freestream%U
                else
                    P = mirror_across_plane(P, body%mirror_plane)
                    call body%get_induced_velocities_at_point(P, this%freestream, v_d, v_s)
                    body%V_cells_inner(:,i+body%N_panels) = this%freestream%v_inf + this%freestream%U*(v_d + v_s)
                end if

                ! Get mirrored surface velocity
                body%V_cells(:,i+body%N_panels) = body%panels(i)%get_velocity(body%mu, body%sigma, .true., body%N_panels, &
                                                                            body%N_verts, body%asym_flow, this%freestream, &
                                                                            body%V_cells_inner(:,i+body%N_panels)/this%freestream%U)

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
            if (body%calc_adjoint) then
                allocate(body%d_C_p_inc_wrt_vars(this%N_cells), stat=stat)
                call check_allocation(stat, "adjoint incompressible pressures (wrt vars)")
                allocate(body%d_C_p_inc_wrt_mu(this%N_cells), stat=stat)
                call check_allocation(stat, "adjoint incompressible pressures (wrt mu)")
            end if
        end if
        
        if (this%isentropic_rule) then
            allocate(body%C_p_ise(this%N_cells), stat=stat)
            call check_allocation(stat, "isentropic surface pressures")
            if (body%calc_adjoint) then
                allocate(body%d_C_p_ise_wrt_vars(this%N_cells), stat=stat)
                call check_allocation(stat, "adjoint isentropic pressures (wrt vars)")
                allocate(body%d_C_p_ise_wrt_mu(this%N_cells), stat=stat)
                call check_allocation(stat, "adjoint isentropic pressures (wrt mu)")
            end if
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

        if (mirrored) then
            C_P_avg = body%panels(i_panel)%get_avg_pressure_coef(body%mu, body%sigma, mirrored, body%N_panels, body%N_verts, &
                                                                 body%asym_flow, this%freestream, &
                                                                 body%V_cells_inner(:,i_panel+body%N_panels)/this%freestream%U, &
                                                                 rule, M_corr=this%M_inf_corr)
        else
            C_P_avg = body%panels(i_panel)%get_avg_pressure_coef(body%mu, body%sigma, mirrored, body%N_panels, body%N_verts, &
                                                                 body%asym_flow, this%freestream, &
                                                                 body%V_cells_inner(:,i_panel)/this%freestream%U, &
                                                                 rule, M_corr=this%M_inf_corr)
        end if
        
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

                if (body%calc_adjoint) then
                    body%d_C_p_inc_wrt_vars(i) = this%calc_d_avg_pressure_on_panel_wrt_vars(i,body,.false.,"incompressible") 
                    body%d_C_p_inc_wrt_mu(i) = this%calc_d_avg_pressure_on_panel_wrt_mu(i,body,.false.,"incompressible") 
                end if
            end if
        
            ! Isentropic rule
            if (this%isentropic_rule) then
                body%C_p_ise(i) = this%calc_avg_pressure_on_panel(i, body, .false., "isentropic")
                if (body%asym_flow) then
                    body%C_p_ise(i+body%N_panels) = this%calc_avg_pressure_on_panel(i, body, .true., "isentropic")
                end if

                if (body%calc_adjoint) then
                    body%d_C_p_ise_wrt_vars(i) = this%calc_d_avg_pressure_on_panel_wrt_vars(i,body,.false.,"isentropic") 
                    body%d_C_p_ise_wrt_mu(i) = this%calc_d_avg_pressure_on_panel_wrt_mu(i,body,.false.,"isentropic") 
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
        
        integer :: i
        
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
        call json_value_add(p_parent, "N_wake_filaments", body%filament_wake%N_filaments)
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

            if (body%calc_adjoint) then
                ! Write CF sensitivities
                call json_value_create(p_parent)
                call to_object(p_parent, 'CF_sensitivities')
                call json_value_add(p_json, p_parent)
                call json_value_add(p_parent, 'd_CFx', this%CF_sensitivities(:,1))
                call json_value_add(p_parent, 'd_CFy', this%CF_sensitivities(:,2))
                call json_value_add(p_parent, 'd_CFz', this%CF_sensitivities(:,3))
                nullify(p_parent)

                ! Write norms of CF sensitivities
                call json_value_create(p_parent)
                call to_object(p_parent, 'norms_of_CF_sensitivities')
                call json_value_add(p_json, p_parent)
                call json_value_add(p_parent, 'norm_of_d_CFx', this%norms_of_CF_sensitivities(1))
                call json_value_add(p_parent, 'norm_of_d_CFy', this%norms_of_CF_sensitivities(2))
                call json_value_add(p_parent, 'norm_of_d_CFz', this%norms_of_CF_sensitivities(3))
                nullify(p_parent)
            end if


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

        if (verbose) write(*,'(a)',advance='no') "    Calculating flow properties at off-body points "

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

        ! Scale by the freestream magnitude
        phi_d = phi_d*this%freestream%U
        phi_s = phi_s*this%freestream%U
        v_d = v_d*this%freestream%U
        v_s = v_s*this%freestream%U

        ! Delete old output file
        call delete_file(points_output_file)

        ! Open output file
        open(newunit=unit, file=points_output_file)

        ! Write header
        write(unit,*) 'x,y,z,phi_inf,phi_d,phi_s,phi,Phi,v_inf_x,v_inf_y,v_inf_z,v_d_x,v_d_y,v_d_z,v_s_x,v_s_y,v_s_z,&
                       v_x,v_y,v_z,V_x,V_y,V_z,V'

        ! Write potentials out to file
        do i=1,N_points

            ! Get freestream properties
            phi_inf = this%freestream%U*inner(points(:,i), this%freestream%c_hat_g)
            v_inf = this%freestream%v_inf

            ! Write to file
            write(unit," (e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, &
                        a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, &
                        e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, e20.13, a, &
                        e20.13, a)", advance="no") &
                            points(1,i),',', points(2,i),',', points(3,i),',', &
                            phi_inf,',', phi_d(i),',', phi_s(i),',', phi_d(i) + phi_s(i),',', &
                            phi_inf + phi_d(i) + phi_s(i),',', v_inf(1),',', &
                            v_inf(2),',', v_inf(3),',', v_d(1,i),',', v_d(2,i),',', v_d(3,i),',', &
                            v_s(1,i),',', v_s(2,i),',', v_s(3,i),',', &
                            v_d(1,i) + v_s(1,i),',', v_d(2,i) + v_s(2,i),',', v_d(3,i) + v_s(3,i),',', &
                            v_inf(1) + v_d(1,i) + v_s(1,i),',',&
                            v_inf(2) + v_d(2,i) + v_s(2,i),',', v_inf(3) + v_d(3,i) + v_s(3,i), ',', &
                            norm2(v_inf + v_d(:,i) + v_s(:,i)), achar(10)

        end do

        close(unit)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_export_off_body_points


    subroutine panel_solver_update_adjoint_A_row(this, body, cp, d_AIC_row, i_panel, inf_adjoint, mirrored_panel)
        ! places the doublet influences sensitivities (sparse vectors) in the appropriate place in the 
        ! d_A matrix
        
        !subroutine panel_solver_update_system_row(this, body, cp, A_row, I_known_i, i_panel, source_inf, doublet_inf, mirrored_panel)
        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(control_point),intent(in) :: cp
        type(sparse_vector),dimension(this%N_unknown),intent(inout) :: d_AIC_row
        integer,intent(in) :: i_panel
        type(sparse_vector),dimension(3),intent(in) :: inf_adjoint
        logical,intent(in) :: mirrored_panel

        integer :: k, index

        ! ! Add source influence depending whether the panel has sources
        ! if (body%panels(i_panel)%has_sources) then

        !     ! Loop through influencing panels
        !     do k=1,body%panels(i_panel)%S_dim

        !         ! Need to shift indices for mirrored panels
        !         if (mirrored_panel) then

        !             ! Check for dependence across mirror plane
        !             if (body%panels(i_panel)%i_panel_s(k) > body%N_panels) then
        !                 index = body%panels(i_panel)%i_panel_s(k) - body%N_panels
        !             else
        !                 index = body%panels(i_panel)%i_panel_s(k) + body%N_panels
        !             end if

        !         else
        !             ! Check for dependence across mirror plane
        !             if (body%panels(i_panel)%i_panel_s(k) > body%N_panels) then
        !                 index = body%panels(i_panel)%i_panel_s(k) - body%N_panels
        !             else
        !                 index = body%panels(i_panel)%i_panel_s(k)
        !             end if
        !         end if

        !         ! Add to known influences if sigma is known
        !         if (this%sigma_known(index)) then
        !             I_known_i = I_known_i + source_inf(k)*body%sigma(index)

        !         ! Add to A matrix if not
        !         else
        !             A_row(this%P(this%i_sigma_in_sys(index))) = A_row(this%P(this%i_sigma_in_sys(index))) + source_inf(k)
        !         end if

        !     end do

        ! end if


        

        ! Add doublet influence sensitivities
        ! Loop through influencing vertices
        do k=1,body%panels(i_panel)%M_dim

            ! Need to shift indices for mirrored panels
            if (mirrored_panel) then
                write(*,*) "!!! Cannot calculate adjoint for mirrored mesh. Quitting..."
                stop
                ! ! Check for dependence across mirror plane
                ! if (body%panels(i_panel)%i_vert_d(k) > body%N_verts) then
                !     index = body%panels(i_panel)%i_vert_d(k) - body%N_verts
                ! else
                !     index = body%panels(i_panel)%i_vert_d(k) + body%N_verts
                ! end if

            else

                ! Check for dependence across mirror plane
                if (body%panels(i_panel)%i_vert_d(k) > body%N_verts) then
                    index = body%panels(i_panel)%i_vert_d(k) - body%N_verts
                else
                    index = body%panels(i_panel)%i_vert_d(k)
                end if
            end if

            ! Update
            call d_AIC_row(this%P(index))%sparse_add(inf_adjoint(k))

        end do

    end subroutine panel_solver_update_adjoint_A_row


    subroutine panel_solver_assemble_adjoint_b_vector(this, body)
        ! Sets up the {BC} sensitivity vector used in adjoint calculation.

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body

        integer :: i, ind
        real,dimension(3) :: x


        ! Vector for source-free target potential calculation
        x = matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        

        !!!!!! NOTE:  in the non-adjoint assemble_BC_vector subroutine, we are calculating this%BC, 
        !            but here we calculate this%d_b instead of this%d_BC because 
        !            this%b = this%BC - this%I_known    and in the adjoint formulations I_known is all zeros 
        !!!!!!!

        ! Loop through control points
        do i=1,body%N_cp

            if (this%use_sort_for_cp) then
                ind = this%P(i)
            else
                ind = i
            end if

            ! Check boundary condition type
            select case (body%cp(i)%bc)

            ! Neumann (mass flux)
            case (ZERO_NORMAL_MF)
                this%d_b_vector(ind) = body%cp(i)%d_n_g%broadcast_vector_dot_element(this%freestream%c_hat_g)
                call this%d_b_vector(ind)%broadcast_element_times_scalar(-1.)
            ! Source-free formulation
            case (SF_POTENTIAL)
                !this%BC(ind) = -inner(x, body%cp(i)%loc)
                this%d_b_vector(ind) = body%cp(i)%d_loc%broadcast_vector_dot_element(-x)

            case default
                write(*,*) "!!! '", body%cp(i)%bc, "' is not a valid boundary condition for &
                adjoint calculation. Must use Neumann Zero Normal Mass Flux with control &
                points at vertices. Quitting..."
                stop

            end select
        end do

    end subroutine panel_solver_assemble_adjoint_b_vector


    subroutine panel_solver_calc_cell_velocities_adjoint(this, body)
        ! Calculates the surface velocities sensitivities on the mesh panels

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: P
        type(sparse_matrix) :: d_P, d_P_term2, zeros_vars, zeros_mu, d_inner_flow_wrt_vars, d_inner_flow_wrt_mu

        
        ! Determine number of surface cells
        this%N_cells = body%N_panels
        ! if (body%asym_flow) then
        !     this%N_cells = this%N_cells*2
        ! end if
        
        ! Allocate cell velocity storage
        allocate(body%d_V_cells_inner_wrt_vars(this%N_cells), stat=stat)
        call check_allocation(stat, "inner velocity senstivitviy wrt design variables") 

        allocate(body%d_V_cells_inner_wrt_mu(this%N_cells), stat=stat)
        call check_allocation(stat, "inner velocity senstivitviy wrt mu")
        
        allocate(body%d_V_cells_wrt_vars(this%N_cells), stat=stat)
        call check_allocation(stat, "velocity senstivitviy wrt design variables")
        
        allocate(body%d_V_cells_wrt_mu(this%N_cells), stat=stat)
        call check_allocation(stat, "velocity senstivitviy wrt mu")
        

        ! Allocate inner cell velocity storage for Neumann
        if (this%dirichlet) then

            ! d_inner_flow_wrt_vars and d_inner_flow_wrt_mu are both zeros
            call zeros_vars%init(body%adjoint_size)
            call zeros_mu%init(size(body%mu))

        end if


        ! Get the surface velocity on each existing panel
        !$OMP parallel do private(P, d_P, d_P_term2, d_inner_flow_wrt_vars, d_inner_flow_wrt_mu) schedule(static)
        do i=1,body%N_panels
            

            ! Get inner flow
            if (this%dirichlet) then
                !body%V_cells_inner(:,i) = this%inner_flow*this%freestream%U

                !!! dont need to do store this info because the inner flow is the same everywhere in dirichlet
                ! call body%d_V_cells_inner_wrt_vars(i)%init_from_sparse_matrix(this%d_inner_flow_wrt_vars)
                ! call body%d_V_cells_inner_wrt_vars(i)%broadcast_element_times_scalar(this%freestream%U)

                ! call body%d_V_cells_inner_wrt_mu(i)%init_from_sparse_matrix(this%d_inner_flow_wrt_mu)
                ! call body%d_V_cells_inner_wrt_mu(i)%broadcast_element_times_scalar(this%freestream%U)
                
                ! future subroutines access the panel_solver%d_V_cells_inner_wrt_vars/mu attributes, so populate them
                call body%d_V_cells_inner_wrt_vars(i)%init(body%adjoint_size)
                call body%d_V_cells_inner_wrt_mu(i)%init(size(body%mu))

                !!!!!!!!!!!!!!!!!! d_V_cells with respect to design variables(vars)   !!!!!!!!!!!!
                body%d_V_cells_wrt_vars(i) = body%panels(i)%get_d_velocity_wrt_vars(body%mu,.false., body%N_panels, &
                                                                        body%N_verts, body%asym_flow, &
                                                                        this%freestream, zeros_vars)

                !!!!!!!!!!!!!!! end d_V_cells wrt vars !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                

                !!!!!!!!!!!!!!!!!! d_V_cells with respect to mu   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                body%d_V_cells_wrt_mu(i) = body%panels(i)%get_d_velocity_wrt_mu(body%mu,.false., body%N_panels, &
                body%N_verts, body%asym_flow, &
                this%freestream, zeros_mu)

                !!!!!!!!!!!!!!!!!   end d_V_cells with respect to mu  !!!!!!!!!!!!!!!!!!!!!!!!!!

                

            else
                P = body%panels(i)%centr - 1.e-10*body%panels(i)%n_g
                
                ! calc d_P
                call d_P%init_from_sparse_matrix(body%panels(i)%d_centr) 
                call d_P_term2%init_from_sparse_matrix(body%panels(i)%d_n_g)
                call d_P_term2%broadcast_element_times_scalar(-1.e-10)
                call d_P%sparse_add(d_P_term2)
                

                !!!!!!!!!!!!!!!!!! d_V_inner with respect to design variables(vars)   !!!!!!!!!!!!!!!!!!!!!!
                body%d_V_cells_inner_wrt_vars(i) = body%get_d_v_inner_at_point_wrt_vars(P, d_P, this%freestream)
                
                call d_inner_flow_wrt_vars%init_from_sparse_matrix(body%d_V_cells_inner_wrt_vars(i))
                call d_inner_flow_wrt_vars%broadcast_element_times_scalar(1./this%freestream%U)
                !!!!!!!!!!!!!!! end d_V_inner with respect to design variables(vars) !!!!!!!!!!!!!!!!!!!!!!!
                


                !!!!!!!!!!!!!!!!!! d_V_inner with respect to mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                body%d_V_cells_inner_wrt_mu(i) = body%get_d_v_inner_at_point_wrt_mu(P, this%freestream)
                
                call d_inner_flow_wrt_mu%init_from_sparse_matrix(body%d_V_cells_inner_wrt_mu(i))
                call d_inner_flow_wrt_mu%broadcast_element_times_scalar(1./this%freestream%U)
                !!!!!!!!!!!!!!!!!! end d_V_inner with respect to design mu !!!!!!!!!!!!!!!!!!!!!!


                !!!!!!!!!!!!!!!!!! d_V_cells with respect to design variables(vars)   !!!!!!!!!!!!
                body%d_V_cells_wrt_vars(i) = body%panels(i)%get_d_velocity_wrt_vars(body%mu,.false., body%N_panels, &
                                                                        body%N_verts, body%asym_flow, &
                                                                        this%freestream, d_inner_flow_wrt_vars)
                
                ! deallocate stuff for next loop
                deallocate(d_P%columns, d_P_term2%columns, d_inner_flow_wrt_vars%columns)
                !!!!!!!!!!!!!!! end d_V_cells wrt vars !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       



                !!!!!!!!!!!!!!!!!! d_V_cells with respect to mu   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                body%d_V_cells_wrt_mu(i) = body%panels(i)%get_d_velocity_wrt_mu(body%mu,.false., body%N_panels, &
                body%N_verts, body%asym_flow, &
                this%freestream, d_inner_flow_wrt_mu)
                
                ! deallocate stuff for next loop
                deallocate(d_inner_flow_wrt_mu%columns)
                !!!!!!!!!!!!!!!!!   end d_V_cells with respect to mu  !!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
            
        end do
            

    end subroutine panel_solver_calc_cell_velocities_adjoint


    function panel_solver_calc_d_avg_pressure_on_panel_wrt_vars(this, i_panel, body, mirrored, rule) &
                                                        result(d_C_P_avg_wrt_vars)
        ! Calls the panel function with the needed variables to get the average pressure coefficient on that panel

        implicit none
        
        class(panel_solver),intent(in) :: this
        integer,intent(in) :: i_panel
        type(surface_mesh),intent(in) :: body
        logical,intent(in) :: mirrored
        character(len=*),intent(in) :: rule

        type(sparse_matrix) :: d_inner_flow_wrt_vars
        type(sparse_vector) :: d_C_P_avg_wrt_vars

        if (mirrored) then
            write(*,*) "!!! Cannot calculate adjoint avg pressure on panel for mirrored mesh yet. Quitting..."
            stop
        else
            ! sensitivity of C_P_avg with respect to vars (points)
            call d_inner_flow_wrt_vars%init_from_sparse_matrix(body%d_V_cells_inner_wrt_vars(i_panel))
            call d_inner_flow_wrt_vars%broadcast_element_times_scalar(1./this%freestream%U)
            
            d_C_P_avg_wrt_vars = body%panels(i_panel)%get_d_avg_pressure_coef_wrt_vars(body%mu, body%sigma,mirrored, body%N_panels,&
                body%N_verts, body%asym_flow, this%freestream, &
                body%V_cells_inner(:,i_panel)/this%freestream%U,&
                d_inner_flow_wrt_vars, rule, M_corr=this%M_inf_corr)
        end if
        
    end function panel_solver_calc_d_avg_pressure_on_panel_wrt_vars


    function panel_solver_calc_d_avg_pressure_on_panel_wrt_mu(this, i_panel, body, mirrored, rule) &
                                                            result(d_C_P_avg_wrt_mu)
        ! Calls the panel function with the needed variables to get the average pressure coefficient on that panel

        implicit none
        
        class(panel_solver),intent(in) :: this
        integer,intent(in) :: i_panel
        type(surface_mesh),intent(in) :: body
        logical,intent(in) :: mirrored
        character(len=*),intent(in) :: rule

        type(sparse_matrix) :: d_inner_flow_wrt_mu
        type(sparse_vector) :: d_C_P_avg_wrt_mu

        if (mirrored) then
            write(*,*) "!!! Cannot calculate adjoint avg pressure on panel for mirrored mesh yet. Quitting..."
            stop
        else

            ! sensitivity of C_P_avg with respect to mu
            call d_inner_flow_wrt_mu%init_from_sparse_matrix(body%d_V_cells_inner_wrt_mu(i_panel))
            call d_inner_flow_wrt_mu%broadcast_element_times_scalar(1./this%freestream%U)

            d_C_P_avg_wrt_mu = body%panels(i_panel)%get_d_avg_pressure_coef_wrt_mu(body%mu, body%sigma,mirrored, body%N_panels,&
                                                body%N_verts, body%asym_flow, this%freestream, &
                                                body%V_cells_inner(:,i_panel)/this%freestream%U,&
                                                d_inner_flow_wrt_mu, rule, M_corr=this%M_inf_corr)
        end if
        
    end function panel_solver_calc_d_avg_pressure_on_panel_wrt_mu


    subroutine panel_solver_calc_forces_adjoint(this, body)
        ! Calculates the forces
        
        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i
        
        if (verbose) write(*,'(a, a, a)',advance='no') "     Calculating adjoint force sensitivities &
                                            using the ", this%pressure_for_forces, " pressure rule..."

        ! Calculate forces
        select case (this%pressure_for_forces)

        case ('incompressible')
            call this%calc_d_forces_wrt_vars(body, body%C_p_inc, body%d_C_p_inc_wrt_vars)
            call this%calc_d_forces_wrt_mu(body, body%C_p_inc, body%d_C_p_inc_wrt_mu)

        case ('isentropic')
            call this%calc_d_forces_wrt_vars(body, body%C_p_ise, body%d_C_p_ise_wrt_vars)
            call this%calc_d_forces_wrt_mu(body, body%C_p_ise, body%d_C_p_ise_wrt_mu)

        ! case ('second-order')
        !     call this%calc_forces_with_pressure(body, body%C_p_2nd)

        ! case ('slender-body')
        !     call this%calc_forces_with_pressure(body, body%C_p_sln)

        ! case ('linear')
        !     call this%calc_forces_with_pressure(body, body%C_p_lin)

        ! case ('prandtl-glauert')
        !     call this%calc_forces_with_pressure(body, body%C_p_pg)

        ! case ('karman-tsien')
        !     call this%calc_forces_with_pressure(body, body%C_p_kt)

        ! case ('laitone')
        !     call this%calc_forces_with_pressure(body, body%C_p_lai)

        end select

        !!!!!!!!!!!!!!!!!!!!!!!!!!! begin wrt vars !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! init d_C_F_wrt_variables
        call this%d_C_F_wrt_vars%init(body%adjoint_size)
        
        ! Sum discrete force sensitivities
        do i = 1,body%N_panels
            call this%d_C_F_wrt_vars%sparse_add(body%d_cell_forces_wrt_vars(i))
        end do

        ! divide by S_ref
        call this%d_C_F_wrt_vars%broadcast_element_times_scalar(1./body%S_ref)
        !!!!!!!!!!!!!!!!!!!!!!!!!! end wrt vars !!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!! begin wrt mu !!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! init d_C_F_wrt_mu
        call this%d_C_F_wrt_mu%init(body%N_verts)
        
        ! Sum discrete force sensitivities
        do i = 1,body%N_panels
            call this%d_C_F_wrt_mu%sparse_add(body%d_cell_forces_wrt_mu(i))
        end do

        ! divide by S_ref
        call this%d_C_F_wrt_mu%broadcast_element_times_scalar(1./body%S_ref)
        !!!!!!!!!!!!!!!!!!!!!!!!!! end wrt mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        ! Add contributions from mirrored half
        if (body%mirrored .and. .not. body%asym_flow) then
            ! this%C_F = 2.*this%C_F
            ! this%C_F(body%mirror_plane) = 0. ! We know this
            write(*,*) "!!! Cannot calculate adjoint d_C_F for mirrored mesh yet. Quitting..."
            stop
        end if

       
    
    end subroutine panel_solver_calc_forces_adjoint


    subroutine panel_solver_calc_d_forces_wrt_vars(this, body, pressures, d_pressures_wrt_vars)
        ! Calculates the force sensitivities with respect to the design variables
        
        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:),allocatable :: pressures
        type(sparse_vector),dimension(:),allocatable,intent(inout) :: d_pressures_wrt_vars

        type(sparse_matrix) :: term2, term3
        
        integer :: i, stat

        ! Allocate force storage
        allocate(body%d_cell_forces_wrt_vars(this%N_cells), stat=stat)
        call check_allocation(stat, "adjoint cell forces wrt design variables")

        ! Calculate total forces
        !$OMP parallel do private(term2, term3) schedule(static)
        do i=1,body%N_panels
            ! this is the first term
            body%d_cell_forces_wrt_vars(i) = d_pressures_wrt_vars(i)%broadcast_element_times_vector&
                                                                (-body%panels(i)%A*body%panels(i)%n_g)
            ! calc the second term to be added
            term2 = body%panels(i)%d_A%broadcast_element_times_vector(-pressures(i)*body%panels(i)%n_g)

            ! calc the third term to be added
            call term3%init_from_sparse_matrix(body%panels(i)%d_n_g)
            call term3%broadcast_element_times_scalar(-pressures(i)*body%panels(i)%A)

            ! add terms to first term
            call body%d_cell_forces_wrt_vars(i)%sparse_add(term2)
            call body%d_cell_forces_wrt_vars(i)%sparse_add(term3)

            ! Mirror (havent done mirroring for adjoints yet)
            if (body%asym_flow) then
                ! body%dC_f(:,i+body%N_panels) = -pressures(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                write(*,*) "!!! Cannot calculate adjoint d_dC_f wrt to design variables for mirrored mesh yet. Quitting..."
                stop
            end if

            deallocate(term2%columns,term3%columns)

        end do
    
    end subroutine panel_solver_calc_d_forces_wrt_vars


    subroutine panel_solver_calc_d_forces_wrt_mu(this, body, pressures, d_pressures_wrt_mu)
        ! Calculates the force sensitivities with respect to the mu
        
        implicit none
        
        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:),allocatable :: pressures
        type(sparse_vector),dimension(:),allocatable,intent(inout) :: d_pressures_wrt_mu
        
        integer :: i, stat

        ! Allocate force storage
        allocate(body%d_cell_forces_wrt_mu(body%N_panels), stat=stat)
        call check_allocation(stat, "adjoint cell forces wrt doublet strengths (mu)")

        ! Calculate total forces
        !$OMP parallel do schedule(static)
        do i=1,body%N_panels
            ! this is the first term
            body%d_cell_forces_wrt_mu(i) = d_pressures_wrt_mu(i)%broadcast_element_times_vector&
                                                                (-body%panels(i)%A*body%panels(i)%n_g)

            ! Mirror (havent done mirroring for adjoints yet)
            if (body%asym_flow) then
                ! body%dC_f(:,i+body%N_panels) = -pressures(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                write(*,*) "!!! Cannot calculate adjoint d_dC_f wrt to mu for mirrored mesh yet. Quitting..."
                stop
            end if

        end do
    
    end subroutine panel_solver_calc_d_forces_wrt_mu


    function panel_solver_solve_for_adjoint_vT_terms(this,body) result(vT_forces)
        implicit none
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        ! integer,intent(inout) :: solver_stat

        real,dimension(:,:),allocatable :: A_p
        real,dimension(:),allocatable :: b_p, x, R_cp
        integer :: stat, i, j, N
        integer(8) :: start_count, end_count
        real(16) :: count_rate, norm_res, max_res

        real,dimension(:,:),allocatable :: d_forces_wrt_mu, vT_forces, d_moments_wrt_mu, expanded
        logical :: tall

        !!!!!!!!!!!!!!!!!!!!!!!!!!!! Forces v^T terms !!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! expand sparse matrix to a full matrix N by 3
        tall = .true.
        expanded = this%d_C_F_wrt_mu%expand(tall)

        
        N = this%N_unknown

        ! allocate d_forces_wrt_mu 
        allocate(d_forces_wrt_mu(N,3))

        ! if use_sort_for_cp, need to re-order the values of d_forces_wrt_mu
        if (this%use_sort_for_cp) then
            do j = 1,N 
                ! the sensitivities of CFx, CFy, and CFz wrt to mu need to be re-ordered to 
                ! the sorted order because the mu values were solved for in a sorted order 
                ! with a sorted A matrix
                d_forces_wrt_mu(this%P(j),:) = expanded(j,:)

            end do

        else
            d_forces_wrt_mu = expanded
        end if

        ! allocate vT_forces (N by 3)
        allocate(vT_forces(N,3))

        
        write(*,*) ""
        write(*,*) "        Solving  [A]^T v = g  for d_Cx, d_Cy, d_Cz"
        
        do i=1,3
            allocate(A_p, source=transpose(this%A), stat=stat)
            call check_allocation(stat, "solver copy of  A transpose")
            allocate(b_p, source=d_forces_wrt_mu(:,i), stat=stat)
            call check_allocation(stat, "solver copy of  d_forces_wrt_mu(:,i)")

            call system_clock(start_count, count_rate)
            call GMRES(N, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, &
            this%solver_iterations, x)
            ! call lu_solve(N, A_p, b_p, x)
            ! call householder_ls_solve(body%N_cp, N, A_p, b_p, x)
            call system_clock(end_count)

            ! write(*,*) " x vector"
            ! do j = 1,N 
            !     write(*,*) x(j)
            ! end do

            ! Get residual vector
            R_cp = matmul(transpose(this%A), x) - d_forces_wrt_mu(:,i)
            
            ! Calculate residual parameters
            max_res = maxval(abs(R_cp))
            norm_res = sqrt(sum(R_cp*R_cp))

            
            
            write(*,*) "        Maximum residual:", max_res
            write(*,*) "        Norm of residual:", norm_res
        
            ! write(*,*) " A^T matrix"
            ! do j = 1,N 
            !     write(*,*) this%A(:,j)
            ! end do
            
            ! write(*,*) " g vector"
            ! do j = 1,N 
            !     write(*,*) d_forces_wrt_mu(j,i)
            ! end do

            ! Check
            if (isnan(norm_res)) then
                write(*,*) "!!! Linear system failed to solve [A]^T v = g ."
                ! return
            end if


            ! ! store forces v^T term:if use_sort_for_cp, need to re-order the values of vT_forces
            ! if (this%use_sort_for_cp) then
            !     do j = 1,N 
            !         ! undo the sorting that was needed to solve for x
            !         vT_forces(j,i) = x(this%P(j))

            !     end do

            ! else
            !     vT_forces(:,i) = x
            ! end if

            vT_forces(:,i) = x

            deallocate(b_p, x)
            deallocate(A_p)
    
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!! end Forces v^T terms !!!!!!!!!!!!!!!!!!!!!!!!!!!
        

    end function panel_solver_solve_for_adjoint_vT_terms


    subroutine panel_solver_calc_total_derivative(this,body, vT_forces)
        ! this is the end game. calc total force and moment sensitivities
        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:,:),allocatable,intent(inout) :: vT_forces

        real,dimension(:,:),allocatable :: d_AIC_i, d_CF_wrt_vars, adjoint_CF
        real,dimension(:),allocatable :: f_i, d_b_i, vT_n , d_AIC_times_mu, mu_vector
        integer :: i,j,k,m, N_original_verts, N_unknown, p
        
        N_original_verts = body%N_original_verts
        
        ! this may need to be adjusted for appended wakes
        N_unknown = this%N_unknown
        
        ! allocations
        allocate(d_AIC_i(N_unknown,N_unknown))
        allocate(d_b_i(N_unknown))
        allocate(mu_vector(N_unknown))
        allocate(this%CF_sensitivities(3*N_original_verts,3))

        ! expand d_C_F_wrt_vars into a full matrix. tall = true (N by 3)
        d_CF_wrt_vars = this%d_C_F_wrt_vars%expand(.true.)

        ! if use_sort_for_cp, need to use the sorted values of mu
        if (this%use_sort_for_cp) then

            ! get the sorted values of mu
            do i = 1,N_unknown
                mu_vector(this%P(i)) = body%mu(i)
            end do

        else
            mu_vector = body%mu
        end if

        
        ! for CF_x, CF_y, and CF_z
        do m=1,3
           
            ! for each design variable
            do i=1,3*N_original_verts

                ! k = Number of columns of A matrix
                do k=1,N_unknown 

                    ! populate i^th slice of d_b_vector
                    d_b_i(k) = this%d_b_vector(k)%get_value(i)

                    ! j = Number of rows of A matrix
                    do j=1,N_unknown
                        ! populate i^th slice of the d_A tensor
                        d_AIC_i(j,k) = this%d_A_matrix(j,k)%get_value(i)

                    end do

                end do
                ! for each design variable, multiply d_A_matrix i^th slice by mu
                
                d_AIC_times_mu = matmul(-d_AIC_i, mu_vector)

                ! add d_b_vector i^th slice
                f_i = d_AIC_times_mu + d_b_i
                
                ! NOTE: the result of the dot product below doesn't depend on order as long as the inputs
                ! are ordered the same. 
                this%CF_sensitivities(i,m) = dot_product(vT_forces(:,m),f_i) + d_CF_wrt_vars(i,m)
                
            end do ! end i loop
                
            ! calc and store L2 norms
            this%norms_of_CF_sensitivities(m) = sqrt(sum(this%CF_sensitivities(:,m)*this%CF_sensitivities(:,m)))
            

        end do

        if (verbose) then
            write(*,*) ""
            write(*,*) "       d_Cx:                   d_Cy:                   d_Cz:"
            do i=1,body%adjoint_size
                write(*, '(3(f20.10, 4x))') this%CF_sensitivities(i,:)
                end do
            write(*,*) ""
            write(*,*) ""
            write(*,'((A), f20.10)') "Norm of adjoint d_CFx = ", this%norms_of_CF_sensitivities(1)
            write(*,'((A), f20.10)') "Norm of adjoint d_CFy = ", this%norms_of_CF_sensitivities(2)
            write(*,'((A), f20.10)') "Norm of adjoint d_CFz = ", this%norms_of_CF_sensitivities(3)
            
            write(*,*) ""
        end if

    

    end subroutine panel_solver_calc_total_derivative
    

    



end module panel_solver_mod