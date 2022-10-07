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
    use preconditioners_mod
    use sparsity_mod
    use sort_mod

    implicit none


    type panel_solver

        real :: corrected_M_inf 
        character(len=:),allocatable :: formulation, pressure_for_forces, matrix_solver, preconditioner, iteration_file
        logical :: incompressible_rule, isentropic_rule, second_order_rule, slender_rule, linear_rule
        logical :: morino, write_A_and_b, sort_system
        logical :: compressible_correction, prandtl_glauert, karman_tsien, laitone
        type(dod),dimension(:,:),allocatable :: dod_info
        type(dod),dimension(:,:,:),allocatable :: wake_dod_info
        type(flow) :: freestream
        real :: norm_res, max_res, tol, rel
        real,dimension(3) :: C_F, inner_flow
        real,dimension(:,:),allocatable :: A
        real,dimension(:),allocatable :: b
        integer,dimension(:),allocatable :: P
        integer :: N, N_cells, block_size, max_iterations, N_unknown, N_d_unknown, N_s_unknown, solver_iterations

        contains

            ! Initialization
            procedure :: init => panel_solver_init
            procedure :: parse_solver_settings => panel_solver_parse_solver_settings
            procedure :: parse_processing_settings => panel_solver_parse_processing_settings

            ! Setup
            procedure :: init_dirichlet => panel_solver_init_dirichlet
            procedure :: calc_domains_of_dependence => panel_solver_calc_domains_of_dependence
            procedure :: set_permutation => panel_solver_set_permutation

            ! Solve
            procedure :: solve => panel_solver_solve
            procedure :: calc_source_strengths => panel_solver_calc_source_strengths
            procedure :: update_system_row => panel_solver_update_system_row
            procedure :: calc_body_influences => panel_solver_calc_body_influences
            procedure :: calc_wake_influences => panel_solver_calc_wake_influences
            procedure :: check_system => panel_solver_check_system
            procedure :: write_system => panel_solver_write_system
            procedure :: solve_system => panel_solver_solve_system

            ! Post-processing
            procedure :: calc_cell_velocities => panel_solver_calc_cell_velocities
            procedure :: calc_point_velocities => panel_solver_calc_point_velocities
            procedure :: calc_surface_potentials => panel_solver_calc_surface_potentials
            procedure :: calc_pressures => panel_solver_calc_pressures
            procedure :: subsonic_pressure_correction => panel_solver_subsonic_pressure_correction
            procedure :: calc_crit_mach => panel_solver_calc_crit_mach
            procedure :: calc_forces => panel_solver_calc_forces
            procedure :: calc_forces_with_pressure => panel_solver_calc_forces_with_pressure

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
        this%morino = this%formulation == 'morino'

        ! Get matrix solver settings
        if (this%freestream%supersonic) then
            call json_xtnsn_get(solver_settings, 'matrix_solver', this%matrix_solver, 'GMRES')
        else
            call json_xtnsn_get(solver_settings, 'matrix_solver', this%matrix_solver, 'GMRES')
        end if
        call json_xtnsn_get(solver_settings, 'block_size', this%block_size, 400)
        call json_xtnsn_get(solver_settings, 'tolerance', this%tol, 1.e-12)
        call json_xtnsn_get(solver_settings, 'relaxation', this%rel, 0.8)
        call json_xtnsn_get(solver_settings, 'max_iterations', this%max_iterations, 1000)
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
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.correction_mach_number', this%corrected_M_inf, 0.0)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.prandtl-glauert', this%prandtl_glauert, .false.)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.karman-tsien', this%karman_tsien, .false.)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.laitone', this%laitone, .false.)
        
        ! Check the correction Mach number is subsonic
        if (this%corrected_M_inf < 0.0 .or. this%corrected_M_inf >= 1.0) then
            write(*,*) "!!! The pressure correction Mach number must be greater than zero and subsonic. Quitting..."
            stop
        end if

        ! Check freestream Mach number is set to 0 if pressure correction is selected
        if (((this%prandtl_glauert) .or. (this%karman_tsien) .or. (this%laitone)) .and. (this%freestream%M_inf /= 0.0)) then
            write(*,*) "!!! In order to apply a subsonic pressure correction, the freestream Mach number must be set to '0'."
            write(*,*) "!!! Include the desired Mach number for the correction calculations as 'correction_mach_number' under"
            write(*,*) "!!! 'post-processing'-'subsonic_pressure_correction'. Quitting..."
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

        ! Get offset
        call json_xtnsn_get(solver_settings, 'control_point_offset', offset, 1.e-6)
        if (offset <= 0.) then
            write(*,*) "!!! Control point offset must be greater than 0. Defaulting to 1e-6."
            offset = 1.e-6
        end if
        
        ! Place control points
        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Placing control points using offset of ", offset, "..."

        ! Place control points inside the body
        if (this%morino .or. this%formulation == 'source-free') then
            call body%place_interior_control_points(offset)
        end if

        ! Determine size of linear system
        if (body%asym_flow) then
            this%N = body%N_cp*2
        else
            this%N = body%N_cp
        end if

        ! Calculate target inner flow
        this%inner_flow = this%freestream%c_hat_g
        if (this%formulation == 'source-free') then
            this%inner_flow = this%inner_flow - matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        end if

        if (verbose) write(*,'(a, i6, a)') "Done. Placed", body%N_cp, " control points."
    
    end subroutine panel_solver_init_dirichlet


    subroutine panel_solver_calc_domains_of_dependence(this, body)
        ! Determines the domains of dependence for each control point based on the freestream condition

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k, stat
        real,dimension(3) :: vert_loc, mirrored_vert_loc
        logical,dimension(:),allocatable :: verts_in_dod, verts_in_dod_of_mirrored
        logical,dimension(:,:),allocatable :: wake_verts_in_dod, wake_verts_in_dod_of_mirrored

        ! For asymmetric flow on a mirrored mesh, all domains of dependence must be calculated. There are no shortcuts.
        ! For symmetric flow on a mirrored mesh, domains of dependence will be the same between mirrored panels and mirrored
        ! control points. So, we just need to calculate the DoD for mirrored control points, and then we're good.

        if (this%freestream%supersonic .and. verbose) write(*,'(a)',advance='no') "     Calculating domains of dependence..."

        ! Allocate arrays for domain of dependence information for the body
        if (body%mirrored) then

            if (body%asym_flow) then ! This is the only case when we have mirrored control points

                ! Whether vertices are in the DoD of the mirrored control point
                allocate(verts_in_dod_of_mirrored(2*body%N_verts), stat=stat)
                call check_allocation(stat, "vertex domain of dependence storage")

            end if

            ! DoD info for panels
            allocate(this%dod_info(2*body%N_panels, this%N), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            ! Whether vertices are in the DoD of the original control point
            allocate(verts_in_dod(2*body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")

        else

            ! DoD info for panels
            allocate(this%dod_info(body%N_panels, this%N), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            ! Whether vertices are in the DoD of the original control point
            allocate(verts_in_dod(body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")

        end if

        ! Allocate arrays for domain of dependence information for the wake
        if (body%mirrored) then

            if (body%asym_flow) then ! Wake is not mirrored but we do have more control points

                ! DoD info for panels
                allocate(this%wake_dod_info(body%wake%N_max_strip_panels, body%wake%N_strips, this%N), stat=stat)
                call check_allocation(stat, "wake strip domain of dependence storage")

                ! Whether vertices are in the DoD of the original control point
                allocate(wake_verts_in_dod(body%wake%N_max_strip_verts, body%wake%N_strips), stat=stat)
                call check_allocation(stat, "vertex domain of dependence storage")

                ! Whether vertices are in the DoD of the mirrored control point
                allocate(wake_verts_in_dod_of_mirrored(body%wake%N_max_strip_verts, body%wake%N_strips), stat=stat)
                call check_allocation(stat, "vertex domain of dependence storage")

            else ! Wake is mirrored and we have more control points

                ! DoD info for panels
                allocate(this%wake_dod_info(2*body%wake%N_max_strip_panels, body%wake%N_strips, this%N), stat=stat)
                call check_allocation(stat, "wake strip domain of dependence storage")

                ! Whether vertices are in the DoD of the original control point
                allocate(wake_verts_in_dod(2*body%wake%N_max_strip_verts, body%wake%N_strips), stat=stat)
                call check_allocation(stat, "vertex domain of dependence storage")

            end if

        else

            ! DoD info for panels
            allocate(this%wake_dod_info(body%wake%N_max_strip_panels, body%wake%N_strips, this%N), stat=stat)
            call check_allocation(stat, "wake strip domain of dependence storage")

            ! Whether vertices are in the DoD of the original control point
            allocate(wake_verts_in_dod(body%wake%N_max_strip_verts, body%wake%N_strips), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")

        end if

        ! If the freestream is subsonic, these don't need to be checked
        if (this%freestream%supersonic) then

            ! Loop through control points
            !$OMP parallel do private(i, k, vert_loc, mirrored_vert_loc, verts_in_dod, verts_in_dod_of_mirrored) &
            !$OMP & private(wake_verts_in_dod_of_mirrored, wake_verts_in_dod)
            do j=1,body%N_cp

                ! Get whether body vertices are in the DoD of this control point
                verts_in_dod(1:body%N_verts) = body%get_verts_in_dod_of_point(body%cp(:,j), this%freestream, .false.)

                if (body%mirrored) then

                    ! Mirrored vertices and original control point
                    verts_in_dod(body%N_verts+1:) = body%get_verts_in_dod_of_point(body%cp(:,j), this%freestream, .true.)

                    if (body%asym_flow) then

                        ! Get whether original body vertices are in the DoD of the mirrored control point
                        verts_in_dod_of_mirrored(1:body%N_verts) = &
                            body%get_verts_in_dod_of_point(body%cp_mir(:,j), this%freestream, .false.)

                        ! Mirrored vertices and mirrored control point
                        verts_in_dod_of_mirrored(body%N_verts+1:) = &
                            body%get_verts_in_dod_of_point(body%cp_mir(:,j), this%freestream, .true.)
                    end if

                end if

                ! Loop through wake strips
                do i=1,body%wake%N_strips

                    ! Get whether wake vertices are in the DoD of this control point
                    wake_verts_in_dod(1:body%wake%strips(i)%N_verts,i) = &
                        body%wake%strips(i)%get_verts_in_dod_of_point(body%cp(:,j), this%freestream, .false.)

                    if (body%mirrored) then

                        if (body%asym_flow) then

                            ! Original vertex and mirrored control point
                            wake_verts_in_dod_of_mirrored(1:body%wake%strips(i)%N_verts,i) = &
                                body%wake%strips(i)%get_verts_in_dod_of_point(body%cp_mir(:,j), this%freestream, .false.)

                        else

                            ! Mirrored vertex and original control point
                            wake_verts_in_dod(body%wake%strips(i)%N_verts+1:,i) = &
                                body%wake%strips(i)%get_verts_in_dod_of_point(body%cp(:,j), this%freestream, .true.)

                        end if
                    end if
                end do

                ! Loop through body panels
                do i=1,body%N_panels

                    ! Original panel and original control point
                    this%dod_info(i,j) = body%panels(i)%check_dod(body%cp(:,j), this%freestream, verts_in_dod)

                    if (body%mirrored) then

                        if (body%asym_flow) then

                            ! Check DoD for mirrored panel and original control point
                            this%dod_info(i+body%N_panels,j) = body%panels(i)%check_dod(body%cp(:,j), this%freestream, &
                                                                                        verts_in_dod, &
                                                                                        .true., body%mirror_plane)

                            ! Check DoD for original panel and mirrored control point
                            this%dod_info(i,j+body%N_cp) = body%panels(i)%check_dod(body%cp_mir(:,j), this%freestream, &
                                                                                    verts_in_dod_of_mirrored)

                            ! Check DoD for mirrored panel and mirrored control point
                            this%dod_info(i+body%N_panels,j+body%N_cp) = body%panels(i)%check_dod(body%cp_mir(:,j), &
                                                                                                  this%freestream, &
                                                                                                  verts_in_dod_of_mirrored, &
                                                                                                  .true., body%mirror_plane)

                        else

                            ! Check DoD for mirrored panel and original control point
                            this%dod_info(i+body%N_panels,j) = body%panels(i)%check_dod(body%cp(:,j), this%freestream, &
                                                                                        verts_in_dod, &
                                                                                        .true., body%mirror_plane)

                        end if
                    end if
                end do

                ! Loop through wake strip panels
                do i=1,body%wake%N_strips
                    do k=1,body%wake%strips(i)%N_panels

                        ! Check DoD for panel and original control point
                        this%wake_dod_info(k,i,j) = body%wake%strips(i)%panels(k)%check_dod(body%cp(:,j), this%freestream, &
                                                                                                  wake_verts_in_dod(:,i))

                        if (body%mirrored) then

                            if (body%asym_flow) then

                                ! Check DoD for panel and mirrored control point
                                this%wake_dod_info(k,i,j+body%N_cp) = &
                                    body%wake%strips(i)%panels(k)%check_dod(body%cp_mir(:,j), &
                                                                            this%freestream, &
                                                                            wake_verts_in_dod_of_mirrored(:,i))

                            else

                                ! Check DoD for mirrored panel and original control point
                                this%wake_dod_info(k+body%wake%N_max_strip_panels,i,j) = &
                                    body%wake%strips(i)%panels(k)%check_dod(body%cp(:,j), this%freestream, &
                                                                            wake_verts_in_dod(:,i), &
                                                                            .true., body%mirror_plane)

                            end if
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
        integer :: i, j, i_neighbor, i_vert
        integer,dimension(:),allocatable :: P_inv_1, P_inv_2

        ! Sort vertices in compressibility direction
        ! We proceed from most downstream to most upstream so as to get an upper-pentagonal matrix
        ! This sorting seems to improve the performance of iterative solvers for supersonic flow as well
        if (this%sort_system) then

            if (verbose) write(*,'(a)',advance='no') "     Permuting vertices for efficient system solution..."

            ! Sort by actual vertex location first
            allocate(x(this%N))

            ! Add compressibility distance of each vertex
            do i=1,body%N_verts
                x(i) = -inner(this%freestream%c_hat_g, body%vertices(i)%loc)

                ! Mirrored vertex
                if (body%asym_flow) then
                    x(i+body%N_verts) = -inner(this%freestream%c_hat_g, &
                                               mirror_across_plane(body%vertices(i)%loc, body%mirror_plane))
                end if

            end do

            ! Get inverse permutation based on vertex location
            call insertion_arg_sort(x, P_inv_1)

            ! Now sort again based on location of most-downstream neighboring vertex
            ! For vertices sharing a downstream vertex, this will preserve the previous order, since insertion sort is stable
            do i=1,size(P_inv_1)

                ! Get vertex index
                i_vert = P_inv_1(i)

                ! Mirrored vertex
                if (i_vert > body%N_verts) then

                    ! Loop through neighboring vertices to find the furthest back
                    x(i) = huge(x(i))
                    do j=1,body%vertices(i_vert-body%N_verts)%adjacent_vertices%len()
                        call body%vertices(i_vert-body%N_verts)%adjacent_vertices%get(j, i_neighbor)
                        x(i) = min(x(i), -inner(this%freestream%c_hat_g, mirror_across_plane(body%vertices(i_neighbor)%loc, &
                                                                                             body%mirror_plane)))
                    end do

                else

                    ! Loop through neighboring vertices to find the furthest back
                    x(i) = huge(x(i))
                    do j=1,body%vertices(i_vert)%adjacent_vertices%len()
                        call body%vertices(i_vert)%adjacent_vertices%get(j, i_neighbor)
                        x(i) = min(x(i), -inner(this%freestream%c_hat_g, body%vertices(i_neighbor)%loc))
                    end do

                end if

            end do

            ! Get inverse permutation
            ! Sorts into increasing order
            call insertion_arg_sort(x, P_inv_2)

            ! Get overall permuation
            allocate(this%P(this%N))
            do i=1,this%N
                this%P(P_inv_1(P_inv_2(i))) = i
            end do

            if (verbose) write(*,*) "Done."

        else

            ! Identity permutation
            allocate(this%P(this%N))
            this%P(1:body%N_verts) = body%vertex_ordering
            if (this%N > body%N_verts) then
                this%P(body%N_verts+1:) = body%vertex_ordering + body%N_verts
            end if

        end if
    
    end subroutine panel_solver_set_permutation


    subroutine panel_solver_solve(this, body)
        ! Calls the relevant subroutine to solve the case based on the selected formulation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        ! Calculate source strengths
        call this%calc_source_strengths(body)

        ! Calculate body influences
        call this%calc_body_influences(body)

        ! Calculate wake influences
        if (body%wake%N_panels > 0) call this%calc_wake_influences(body)

        ! Solve the linear system
        call this%solve_system(body)

        ! Calculate velocities
        call this%calc_cell_velocities(body)
        if (doublet_order == 2) call this%calc_point_velocities(body)
        
        ! Calculate potentials
        call this%calc_surface_potentials(body)

        ! Calculate pressures
        call this%calc_pressures(body)

        ! Calculate forces
        call this%calc_forces(body)

    end subroutine panel_solver_solve


    subroutine panel_solver_calc_source_strengths(this, body)
        ! Calculates the necessary source strengths

        implicit none

        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: N_sigma, i, stat

        ! Determine number of source strengths
        if (source_order == 0) then

            if (body%asym_flow) then
                N_sigma = body%N_panels*2
            else
                N_sigma = body%N_panels
            end if

        else if (source_order == 1) then

            if (body%asym_flow) then
                N_sigma = body%N_verts*2
            else
                N_sigma = body%N_verts
            end if

        end if

        ! Allocate source strength array (yes this does need to be allocated for the source-free formulation)
        allocate(body%sigma(N_sigma), source=0., stat=stat)
        call check_allocation(stat, "source strength array")

        ! Set source strengths
        if (this%morino) then

            if (verbose) write(*,'(a)',advance='no') "     Calculating source strengths..."

            ! Use panel normals
            if (source_order == 0) then

                ! Loop through panels
                do i=1,body%N_panels

                    ! Existing panels
                    body%sigma(i) = -inner(body%panels(i)%n_g, this%freestream%c_hat_g)

                    ! Mirrored panels for asymmetric flow
                    if (body%asym_flow) then
                        body%sigma(i+body%N_panels) = -inner(body%panels(i)%n_g_mir, this%freestream%c_hat_g)
                    end if
                end do

            ! Use vertex normals
            else if (source_order == 1) then

                ! Loop through vertices
                do i=1,body%N_verts

                    ! Existing vertices
                    body%sigma(i) = -inner(body%vertices(i)%n_g, this%freestream%c_hat_g)

                    ! Mirrored panels for asymmetric flow
                    if (body%asym_flow) then
                        body%sigma(i+body%N_verts) = -inner(body%vertices(i)%n_g_mir, this%freestream%c_hat_g)
                    end if
                end do

            end if

            if (verbose) write(*,*) "Done."
        end if
    
    end subroutine panel_solver_calc_source_strengths


    subroutine panel_solver_update_system_row(this, body, A_row, phi_cp_s, i_panel, source_inf, doublet_inf, mirrored_panel)
        ! Updates the linear system with the source and doublet influences

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:),allocatable,intent(inout) :: A_row
        real,intent(inout) :: phi_cp_s
        integer,intent(in) :: i_panel
        real,dimension(:),allocatable,intent(in) :: source_inf, doublet_inf
        logical,intent(in) :: mirrored_panel

        integer :: k

        ! Add source influence (if sources are present)
        if (this%morino) then

            ! Constant
            if (source_order == 0) then
                if (mirrored_panel) then
                    phi_cp_s = phi_cp_s + source_inf(1)*body%sigma(i_panel+body%N_panels)
                else
                    phi_cp_s = phi_cp_s + source_inf(1)*body%sigma(i_panel)
                end if

            ! Linear
            else
                do k=1,size(body%panels(i_panel)%i_vert_s)
                    if (mirrored_panel) then
                        phi_cp_s = phi_cp_s + source_inf(k)*body%sigma(body%panels(i_panel)%i_vert_s(k)+body%N_cp)
                    else
                        phi_cp_s = phi_cp_s + source_inf(k)*body%sigma(body%panels(i_panel)%i_vert_s(k))
                    end if
                end do
            end if
        end if

        ! Add doublet influence
        ! This method is the same for linear and quadratic doublets
        ! Loop through panel vertices
        do k=1,size(body%panels(i_panel)%i_vert_d)
            if (mirrored_panel) then
                A_row(this%P(body%panels(i_panel)%i_vert_d(k)+body%N_cp)) = &
                        A_row(this%P(body%panels(i_panel)%i_vert_d(k)+body%N_cp)) + doublet_inf(k)
            else
                A_row(this%P(body%panels(i_panel)%i_vert_d(k))) = A_row(this%P(body%panels(i_panel)%i_vert_d(k))) + doublet_inf(k)
            end if
        end do
    
    end subroutine panel_solver_update_system_row


    subroutine panel_solver_calc_body_influences(this, body)
        ! Calculates the influence of the body on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k, stat
        real,dimension(:),allocatable :: source_inf, doublet_inf, A_i, A_i_mir
        real :: phi_cp_s, phi_cp_s_mir
        real,dimension(3) :: x

        ! Allocate space for inner potential calculations
        allocate(body%phi_cp_sigma(this%N), source=0., stat=stat)
        call check_allocation(stat, "induced potential vector")

        ! Allocate AIC matrix
        allocate(this%A(this%N, this%N), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")

        ! Allocate b vector
        allocate(this%b(this%N), source=0., stat=stat)
        call check_allocation(stat, "b vector")

        ! Allocate row of A
        allocate(A_i(this%N))
        if (body%asym_flow) then
            allocate(A_i_mir(this%N))
        end if

        if (verbose) write(*,'(a)',advance='no') "     Calculating body influences..."
        
        ! Parameter used for calculating inner potential for the source-free formulation
        if (this%formulation == 'source-free') then
            x = matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        end if

        ! Calculate source and doublet influences from body on each control point
        !$OMP parallel do private(j, source_inf, doublet_inf, k, A_i, A_i_mir, phi_cp_s, phi_cp_s_mir) &
        !$OMP schedule(dynamic)
        do i=1,body%N_cp

            ! Initialize
            A_i = 0.
            phi_cp_s = 0.
            if (body%asym_flow) then
                A_i_mir = 0.
                phi_cp_s_mir = 0.
            end if

            ! Loop through panels
            do j=1,body%N_panels

                ! Existing panel on existing control point
                if (this%dod_info(j,i)%in_dod) then
                    
                    ! Calculate influence
                    call body%panels(j)%calc_potential_influences(body%cp(:,i), this%freestream, this%dod_info(j,i), .false., &
                                                        source_inf, doublet_inf)

                    ! Add influence
                    call this%update_system_row(body, A_i, phi_cp_s, j, source_inf, doublet_inf, .false.)

                end if

                ! Calculate mirrored influences
                if (body%mirrored) then

                    if (body%asym_flow) then

                        if (body%vertices(i)%mirrored_is_unique) then

                            ! Recalculate influence of mirrored panel on mirrored control point for compressible flow
                            ! If the flow is symmetric or incompressible, this will be the same as already calculated
                            if (this%dod_info(j+body%N_panels,i+body%N_cp)%in_dod) then
                                if (.not. this%freestream%incompressible) then
                                    call body%panels(j)%calc_potential_influences(body%cp_mir(:,i), this%freestream, &
                                                                        this%dod_info(j+body%N_panels,i+body%N_cp), .true., &
                                                                        source_inf, doublet_inf)
                                end if

                                ! Add influence of mirrored panel on mirrored control point
                                call this%update_system_row(body, A_i_mir, phi_cp_s_mir, j, source_inf, doublet_inf, .true.)
                            end if

                        end if

                        ! Calculate influence of existing panel on mirrored control point
                        if (this%dod_info(j,i+body%N_cp)%in_dod) then
                            if (body%vertices(i)%mirrored_is_unique .or. this%freestream%incompressible) then
                                call body%panels(j)%calc_potential_influences(body%cp_mir(:,i), this%freestream, &
                                                                    this%dod_info(j,i+body%N_cp), .false., &
                                                                    source_inf, doublet_inf)
                            end if

                            ! Add influence of existing panel on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                call this%update_system_row(body, A_i_mir, phi_cp_s_mir, j, source_inf, doublet_inf, .false.)
                            end if
                        end if

                        ! Recalculate mirrored->existing influences for compressible flow
                        if (this%dod_info(j+body%N_panels,i)%in_dod) then
                            if (.not. this%freestream%incompressible) then
                                call body%panels(j)%calc_potential_influences(body%cp(:,i), this%freestream, &
                                                                    this%dod_info(j+body%N_panels,i), .true., source_inf, &
                                                                    doublet_inf)
                            end if

                            ! Add influence of mirrored panel on existing control point
                            call this%update_system_row(body, A_i, phi_cp_s, j, source_inf, doublet_inf, .true.)
                        end if

                    else

                        ! Calculate influence of existing panel on mirrored control point
                        ! This is the same as the influence of the mirrored panel on the existing control point,
                        ! even for compressible flow, since we know the flow is symmetric here
                        if (this%dod_info(j+body%N_panels,i)%in_dod) then
                            call body%panels(j)%calc_potential_influences(body%cp_mir(:,i), this%freestream, &
                                                                this%dod_info(j+body%N_panels,i), &
                                                                .false., source_inf, doublet_inf)

                            ! Add influence of mirrored panel on existing control point
                            call this%update_system_row(body, A_i, phi_cp_s, j, source_inf, doublet_inf, .false.)
                        end if

                    end if
                end if

            end do

            ! Update A matrix with rows
            !$OMP critical
            this%A(this%P(i),:) = A_i

            ! Update for mirrored points
            if (body%asym_flow) then

                ! Enforce doublet strength matching (i.e. for non-unique, mirrored control points, the
                ! doublet strengths must be the same). The RHS for these rows should still be zero.
                if (.not. body%vertices(i)%mirrored_is_unique) then
                    this%A(this%P(i+body%N_cp),this%P(i)) = 1.
                    this%A(this%P(i+body%N_cp),this%P(i+body%N_cp)) = -1.

                ! Otherwise, add influences
                else
                    this%A(this%P(i+body%N_cp),:) = A_i_mir
                end if
            end if

            ! Set potential for morino formulation
            if (this%morino) then
                body%phi_cp_sigma(i) = phi_cp_s

                ! Set for unique mirrored control  points
                if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then
                    body%phi_cp_sigma(i+body%N_cp) = phi_cp_s_mir
                end if

            ! Set target potential for source-free formulation
            else
                this%b(this%P(i)) = -inner(x, body%cp(:,i))

                ! Set for unique mirrored control  points
                if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then
                    this%b(this%P(i+body%N_cp)) = -inner(x, body%cp_mir(:,i))
                end if
            end if
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
        real,dimension(:),allocatable ::  doublet_inf, source_inf, A_i, A_i_mir

        ! Calculate influence of wake
        if (verbose) write(*,'(a)',advance='no') "     Calculating wake influences..."

        ! Allocate A rows
        allocate(A_i(this%N))
        if (body%asym_flow) then
            allocate(A_i_mir(this%N))
        end if

        ! Loop through control points
        !$OMP parallel do private(j, source_inf, doublet_inf, k, A_i, A_i_mir) schedule(dynamic)
        do i=1,body%N_cp

            ! Initialize
            A_i = 0.
            if (body%asym_flow) then
                A_i_mir = 0.
            end if

            ! Get doublet influence from wake strips
            do j=1,body%wake%N_strips
                do l=1,body%wake%strips(j)%N_panels

                    ! Caclulate influence of existing panel on existing control point
                    call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp(:,i), this%freestream, &
                                                                                 this%wake_dod_info(l,j,i), &
                                                                                 .false., source_inf, doublet_inf)

                    ! Add influence
                    do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                        A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                    end do

                    ! Get influence on mirrored control point
                    if (body%mirrored) then

                        if (body%asym_flow) then

                            ! Calculate influence of existing panel on mirrored point
                            if (body%vertices(i)%mirrored_is_unique) then
                                call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp_mir(:,i), this%freestream, &
                                                                         this%wake_dod_info(l,j,i+body%N_cp), .false., &
                                                                         source_inf, doublet_inf)

                                ! Add influence
                                do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                                    A_i_mir(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                            A_i_mir(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                                end do
                            end if

                        else

                            ! Calculate influence of existing panel on mirrored control point
                            ! This is the same as the influence of a mirrored panel on an existing control point
                            ! even for compressible flow, since we know the flow is symmetric here
                            call body%wake%strips(j)%panels(l)%calc_potential_influences(body%cp_mir(:,i), this%freestream, &
                                                                     this%wake_dod_info(l+body%wake%N_max_strip_panels,j,i), & ! No, this is not the DoD for this computation; yes, it is equivalent
                                                                     .false., source_inf, doublet_inf)

                            ! Add influence
                            do k=1,size(body%wake%strips(j)%panels(l)%i_vert_d)
                                A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) = &
                                            A_i(this%P(body%wake%strips(j)%panels(l)%i_vert_d(k))) + doublet_inf(k)
                            end do

                        end if

                    end if
                end do
            end do

            ! Update rows of A
            !$OMP critical
            this%A(this%P(i),:) = this%A(this%P(i),:) + A_i
            if (body%asym_flow) then
                this%A(this%P(i+body%N_cp),:) = this%A(this%P(i+body%N_cp),:) + A_i_mir
            end if
            !$OMP end critical

        end do

        ! Clean up memory
        deallocate(this%wake_dod_info)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_wake_influences


    subroutine panel_solver_check_system(this)
        ! Checks the validity of the linear system

        implicit none
        
        class(paneL_solver),intent(in) :: this

        integer :: i

        if (verbose) write(*,'(a)',advance='no') "     Checking validity of linear system..."

        ! Check for NaNs
        if (any(isnan(this%A))) then
            write(*,*) "!!! Invalid value detected in A matrix. Quitting..."
            stop
        end if
        if (any(isnan(this%b))) then
            write(*,*) "!!! Invalid value detected in b vector. Quitting..."
            stop
        end if

        ! Check for uninfluenced/ing points
        do i=1,this%N
            if (all(this%A(this%P(i),:) == 0.)) then
                write(*,*) "!!! Control point ", i, " is not influenced. Quitting..."
                stop
            end if
            if (all(this%A(:,this%P(i)) == 0.)) then
                write(*,*) "!!! Vertex ", i, " exerts no influence. Quitting..."
                stop
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
        do i=1,this%N
            do j=1,this%N
                write(unit,'(e20.12)',advance='no') this%A(i,j)
            end do
            write(unit,*)
        end do
        close(unit)

        ! b
        open(newunit=unit, file="b_vec.txt")
        do i=1,this%N
            write(unit,*) this%b(i)
        end do
        close(unit)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_write_system


    subroutine panel_solver_solve_system(this, body)
        ! Solves the linear system for the singularity strengths

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        real,dimension(:,:),allocatable :: A_p
        real,dimension(:),allocatable :: b_p, x, R
        integer :: stat, i, j

        ! Set b vector for Morino formulation
        if (this%formulation == "morino") then
            do i=1,this%N
                this%b(this%P(i)) = -body%phi_cp_sigma(i)
            end do
        end if

        ! Run checks
        if (run_checks) call this%check_system()

        ! Write to file
        if (this%write_A_and_b) call this%write_system()

        if (verbose) write(*,'(a, a, a)',advance='no') "     Solving linear system (method: ", this%matrix_solver, ")..."

        ! Precondition
        select case(this%preconditioner)

        ! Diagonal preconditioning
        case ('DIAG')
            call diagonal_preconditioner(this%N, this%A, this%b, A_p, b_p)

        ! No preconditioning
        case default
            allocate(A_p, source=this%A, stat=stat)
            call check_allocation(stat, "solver copy of AIC matrix")
            allocate(b_p, source=this%b, stat=stat)
            call check_allocation(stat, "solver copy of b vector")

        end select

        ! Solve
        select case(this%matrix_solver)

        ! LU decomposition
        case ('LU')
            call lu_solve(this%N, A_p, b_p, x)

        ! QR via Givens rotations for upper-pentagonal
        case ('QRUP')
            call QR_row_givens_solve_UP(this%N, A_p, b_p, x)

        ! QR via fast Givens rotations for upper-pentagonal
        case ('FQRUP')
            call QR_fast_givens_solve_upper_pentagonal(this%N, A_p, b_p, x)

        ! GMRES
        case ('GMRES')
            call GMRES(this%N, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, this%solver_iterations, x)

        ! Purcell's method
        case ('PURC')
            call purcell_solve(this%N, A_p, b_p, x)

        ! Block successive over-relaxation
        case ('BSOR')
            call block_sor_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                 this%max_iterations, this%iteration_file, this%solver_iterations, x)

        ! Adaptive block SOR
        case ('ABSOR')
            this%rel = -1.
            call block_sor_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                 this%max_iterations, this%iteration_file, this%solver_iterations, x)
        
        ! Block Jacobi
        case ('BJAC')
            call block_jacobi_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                    this%max_iterations, this%iteration_file, this%solver_iterations, x)

        ! Optimally relaxed block Jacobi
        case ('ORBJ')
            this%rel = -1.
            call block_jacobi_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                    this%max_iterations, this%iteration_file, this%solver_iterations, x)

        ! Improper specification
        case default
            write(*,*) "!!! ", this%matrix_solver, " is not a valid solver option. Defaulting to GMRES."
            call GMRES(this%N, A_p, b_p, this%tol, this%max_iterations, this%iteration_file, this%solver_iterations, x)

        end select
        if (verbose) write(*,*) "Done."

        ! Clean up memory
        deallocate(A_p)
        deallocate(b_p)

        ! Get residual vector
        allocate(R(this%N))
        R = matmul(this%A, x) - this%b
        deallocate(this%A)

        ! Calculate residual parameters
        this%max_res = maxval(abs(R))
        this%norm_res = sqrt(sum(R*R))
        if (verbose) then
            write(*,*) "        Maximum residual:", this%max_res
            write(*,*) "        Norm of residual:", this%norm_res
        end if

        ! Check
        if (isnan(this%norm_res)) then
            write(*,*) "!!! Linear system failed to produce a valid solution. Quitting..."
            stop
        end if

        ! Transfer solution to body storage
        allocate(body%mu(this%N))
        do i=1,this%N
            body%mu(i) = x(this%P(i))
        end do
        deallocate(x)

        ! Get potentials at control points
        allocate(body%phi_cp_mu(this%N))
        do i=1,this%N
            body%phi_cp_mu(i) = R(this%P(i)) + this%b(this%P(i))
        end do
        body%phi_cp = body%phi_cp_mu + body%phi_cp_sigma
        deallocate(this%b)

    end subroutine panel_solver_solve_system


    subroutine panel_solver_calc_cell_velocities(this, body)
        ! Calculates the surface velocities on the mesh cells

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat, j, N_cycle
        real,dimension(3) :: v_jump, cent

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface velocities..."

        ! Determine number of surface cells
        this%N_cells = body%N_panels
        N_cycle = 1
        if (body%asym_flow) then
            this%N_cells = this%N_cells*2
        end if
        if (doublet_order == 2) then
            this%N_cells = this%N_cells*4
            N_cycle = 4
        end if

        ! Allocate cell velocity storage
        allocate(body%V_cells(3,this%N_cells), stat=stat)
        call check_allocation(stat, "surface velocity vectors")

        ! Calculate the surface velocity on each existing panel
        !$OMP parallel do private(v_jump, j, cent) schedule(static)
        do i=1,body%N_panels

            ! Get velocity jump at centroid
            v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false.)

            ! Calculate velocity
            body%V_cells(:,(i-1)*N_cycle+1) = this%freestream%U*(this%inner_flow + v_jump)

            ! Calculate for subpanels
            if (doublet_order == 2) then
                do j=1,body%panels(i)%N

                    ! Get centroid
                    cent = body%panels(i)%get_subpanel_centroid(j)

                    ! Calculate velocity
                    v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., cent)
                    body%V_cells(:,(i-1)*N_cycle+1+j) = this%freestream%U*(this%inner_flow + v_jump)

                end do
            end if

            ! Calculate surface velocity on each mirrored panel
            if (body%asym_flow) then

                ! Get velocity jump
                v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .true.)

                ! Calculate velocity
                body%V_cells(:,(i-1)*N_cycle+1+this%N_cells/2) = this%freestream%U*(this%inner_flow + v_jump)

                ! Calculate for subpanels
                if (doublet_order == 2) then
                    do j=1,body%panels(i)%N

                        ! Get centroid
                        cent = mirror_across_plane(body%panels(i)%get_subpanel_centroid(j), body%mirror_plane)

                        ! Calculate velocity
                        v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .true., cent)
                        body%V_cells(:,(i-1)*N_cycle+1+j+this%N_cells/2) = this%freestream%U*(this%inner_flow + v_jump)

                    end do
                end if

            end if
        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_cell_velocities


    subroutine panel_solver_calc_point_velocities(this, body)
        ! Calculates the surface velocities the mesh vertices

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat, N_neighbors, j, i_panel
        real,dimension(3) :: v_jump, loc_mir
        real,dimension(:,:),allocatable :: Vs

        if (verbose) write(*,'(a)',advance='no') "     Calculating point velocities..."

        ! Allocate vertex velocity storage
        allocate(body%V_verts_avg(3,this%N), stat=stat)
        call check_allocation(stat, "average point velocity vectors")
        allocate(body%V_verts_std(3,this%N), source=0., stat=stat)
        call check_allocation(stat, "standard deviation of point velocity vectors")

        ! Calculate total potential on outside of mesh and quadratic-doublet velocity discrepancies
        !$OMP parallel do private(Vs, N_neighbors, v_jump, j, i_panel, loc_mir) schedule(dynamic)
        do i=1,body%N_verts

            ! Allocate storage
            N_neighbors = body%vertices(i)%panels_not_across_wake_edge%len()
            allocate(Vs(3,N_neighbors))

            ! Get velocities
            do j=1,N_neighbors

                ! Get panel index
                call body%vertices(i)%panels_not_across_wake_edge%get(j, i_panel)

                ! Get velocity jump at the vertex location
                v_jump = body%panels(i_panel)%get_velocity_jump(body%mu, body%sigma, .false., body%vertices(i)%loc)

                ! Store
                Vs(:,j) = this%freestream%U*(this%inner_flow + v_jump)
            end do

            ! Calculate stats
            body%V_verts_avg(:,i) = sum(Vs, dim=2)/N_neighbors
            if (N_neighbors > 1) then
                do j=1,N_neighbors
                    body%V_verts_std(:,i) = body%V_verts_std(:,i) + (Vs(:,j) - body%V_verts_avg(:,i))**2
                end do
                body%V_verts_std(:,i) = body%V_verts_std(:,i)/(N_neighbors - 1)
            end if

            ! Mirrored vertices
            if (body%asym_flow) then

                ! Get mirrored location
                loc_mir = mirror_across_plane(body%vertices(i)%loc, body%mirror_plane)

                ! Get velocities
                do j=1,N_neighbors

                    ! Get panel index
                    call body%vertices(i)%panels_not_across_wake_edge%get(j, i_panel)

                    ! Get velocity jump at the vertex location
                    v_jump = body%panels(i_panel)%get_velocity_jump(body%mu, body%sigma, .true., loc_mir)

                    ! Store
                    Vs(:,j) = this%freestream%U*(this%inner_flow + v_jump)
                end do

                ! Calculate stats
                body%V_verts_avg(:,i+body%N_verts) = sum(Vs, dim=2)/N_neighbors
                if (N_neighbors > 1) then
                    do j=1,N_neighbors
                        body%V_verts_std(:,i+body%N_verts) = body%V_verts_std(:,i+body%N_verts) &
                                                             + (Vs(:,j) - body%V_verts_avg(:,i))**2
                    end do
                    body%V_verts_std(:,i+body%N_verts) = body%V_verts_std(:,i+body%N_verts)/(N_neighbors - 1)
                end if

            end if

            deallocate(Vs)

        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_point_velocities


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

        ! Factor in freestream velocity
        body%Phi_u = body%Phi_u*this%freestream%U

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_surface_potentials


    subroutine panel_solver_calc_pressures(this, body)
        ! Calculates the surface pressures

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: V_pert
        real :: a, b, c, C_p_vac, lin, sln

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface pressures..."

        ! Calculate vacuum pressure coefficient
        C_p_vac = -2./(this%freestream%gamma*this%freestream%M_inf**2)

        ! Allocate storage
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

        ! Calculate reusable terms for the isentropic rule
        if (this%isentropic_rule) then
            a = 2./(this%freestream%gamma*this%freestream%M_inf**2)
            b = 0.5*(this%freestream%gamma-1.)*this%freestream%M_inf**2
            c = this%freestream%gamma/(this%freestream%gamma-1.)
        end if

        ! Calculate pressures
        !$OMP parallel do private(V_pert, lin, sln) schedule(static)
        do i=1,this%N_cells

            ! Incompressible rule
            if (this%incompressible_rule) then
                body%C_p_inc(i) = 1.-inner(body%V_cells(:,i), body%V_cells(:,i))*this%freestream%U_inv**2
            end if
        
            ! Isentropic rule
            if (this%isentropic_rule) then
                
                body%C_p_ise(i) = 1. - inner(body%V_cells(:,i), body%V_cells(:,i))*this%freestream%U_inv**2
                body%C_p_ise(i) = a*( (1. + b*body%C_p_ise(i))**c - 1.)

                ! Check for NaN and replace with vacuum pressure
                if (isnan(body%C_p_ise(i))) then
                    body%C_p_ise(i) = C_p_vac
                end if
            end if

            ! Get perturbation velocity in the compressible frame
            V_pert = matmul(this%freestream%A_g_to_c, body%V_cells(:,i)-this%freestream%v_inf)

            ! Linear term
            if (this%linear_rule .or. this%slender_rule .or. this%second_order_rule) then
                lin = -2.*V_pert(1)*this%freestream%U_inv
            end if

            ! Slender-body term
            if (this%slender_rule .or. this%second_order_rule) then
                sln = lin - (V_pert(2)**2 + V_pert(3)**2)*this%freestream%U_inv**2
            end if
        
            ! Second-order rule
            if (this%second_order_rule) then
                body%C_p_2nd(i) = sln - (1.-this%freestream%M_inf**2)*V_pert(1)**2*this%freestream%U_inv**2
            end if
        
            ! Slender-body rule
            if (this%slender_rule) then
                body%C_p_sln(i) = sln
            end if
        
            ! Linear rule
            if (this%linear_rule) then
                body%C_p_lin(i) = lin
            end if

        end do

        if (verbose) write(*,*) "Done."
        
        ! Report min and max pressure coefficients
        if (this%incompressible_rule .and. verbose) then
            write(*,*) "        Maximum incompressible pressure coefficient:", maxval(body%C_p_inc)
            write(*,*) "        Minimum incompressible pressure coefficient:", minval(body%C_p_inc)
        end if
        
        if (this%isentropic_rule .and. verbose) then
            write(*,*) "        Maximum isentropic pressure coefficient:", maxval(body%C_p_ise)
            write(*,*) "        Minimum isentropic pressure coefficient:", minval(body%C_p_ise)
        end if
        
        if (this%second_order_rule .and. verbose) then
            write(*,*) "        Maximum second-order pressure coefficient:", maxval(body%C_p_2nd)
            write(*,*) "        Minimum second-order pressure coefficient:", minval(body%C_p_2nd)
        end if
        
        if (this%slender_rule .and. verbose) then
            write(*,*) "        Maximum slender-body pressure coefficient:", maxval(body%C_p_sln)
            write(*,*) "        Minimum slender-body pressure coefficient:", minval(body%C_p_sln)
        end if
        
        if (this%linear_rule .and. verbose) then
            write(*,*) "        Maximum linear pressure coefficient:", maxval(body%C_p_lin)
            write(*,*) "        Minimum linear pressure coefficient:", minval(body%C_p_lin)
        end if
        
        ! Report vacuum pressure coefficient
        if (this%freestream%M_inf > 0.) then
            if (verbose) write(*,*) "        Vacuum pressure coefficient:", C_p_vac
        end if
        
        ! Apply subsonic pressure corrections
        if (this%prandtl_glauert .or. this%karman_tsien .or. this%laitone) then
            call this%subsonic_pressure_correction(body)
        end if

        ! Check if critical Mach number has been exceeded
        if (this%incompressible_rule) then
            if ((this%corrected_M_inf < 1) .and. (this%corrected_M_inf /= 0.)) then
                call this%calc_crit_mach(body)
            end if
        else
            if ((this%freestream%M_inf < 1) .and. (this%freestream%M_inf /= 0.)) then
                call this%calc_crit_mach(body)
            end if
        end if
        
    end subroutine panel_solver_calc_pressures

    
    subroutine panel_solver_subsonic_pressure_correction(this, body)
        ! Apply selected method of correcting subsonic pressures
        
        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        integer :: i,stat
        real :: val_holder_L, val_holder_KT

        if (verbose) write(*,'(a)',advance='no') "     Calculating compressibility pressure corrections..."        
        
        ! Prandtl-Glauert rule
        if (this%prandtl_glauert) then
            ! Allocate storage
            allocate(body%C_p_pg(this%N_cells), stat=stat)
            call check_allocation(stat, "Prandtl-Glauert corrected surface pressures")
            
            ! Perform calculations for Prandtl-Glauert Rule (Modern Compressible Flow by John Anderson EQ 9.36)
            body%C_p_pg = body%C_p_inc / (sqrt(1 - (this%corrected_M_inf**2)))
        end if
        
        ! Laitone rule
        if (this%laitone) then
            ! Allocate storage
            allocate(body%C_p_lai(this%N_cells), stat=stat)
            call check_allocation(stat, "Laitone corrected surface pressures")
            
            ! Perform calculations (Modern Compressible Flow by John Anderson EQ 9.39)
            val_holder_L = this%corrected_M_inf**2 * (1 + (0.5 * (this%freestream%gamma - 1) * this%corrected_M_inf**2)) &
            / (2 * sqrt(1 - this%corrected_M_inf**2))
            
            body%C_p_lai = body%C_p_inc / &
            (sqrt(1 - this%corrected_M_inf**2) + (val_holder_L * body%C_p_inc))
        end if
        
        ! Karman-Tsien rule
        if (this%karman_tsien) then
            ! Allocate storage
            allocate(body%C_p_kt(this%N_cells), stat=stat)
            call check_allocation(stat, "Karman-Tsien corrected surface pressures")
            
            ! Perform calculations (Modern Compressible Flow by John Anderson EQ 9.40)
            val_holder_KT = this%corrected_M_inf**2 / (1 + sqrt(1 - this%corrected_M_inf**2))
            
            body%C_p_kt = body%C_p_inc / &
            (sqrt(1 - this%corrected_M_inf**2) + val_holder_KT * (0.5 * body%C_p_inc))
        end if

        if (verbose) write(*,*) "Done. "  
        
    end subroutine panel_solver_subsonic_pressure_correction
    

    subroutine panel_solver_calc_crit_mach(this, body)
        ! Warn the user if body is in transonic flow

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
                M_inf_selected = this%corrected_M_inf  
            
            case ("prandtl-glauert")
                min_loc = MINLOC(body%C_p_pg)
                C_p_min = MINVAL(body%C_p_pg)
                M_inf_selected = this%corrected_M_inf

            case ("karman-tsien")
                min_loc = MINLOC(body%C_p_kt)
                C_p_min = MINVAL(body%C_p_kt)
                M_inf_selected = this%corrected_M_inf

            case ("laitone")
                min_loc = MINLOC(body%C_p_lai)
                C_p_min = MINVAL(body%C_p_lai)
                M_inf_selected = this%corrected_M_inf

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
        numerator = 1 + ((this%freestream%gamma - 1) / 2) * M_inf_selected**2 
        denominator = 1 + (this%freestream%gamma -1) / 2

        C_p_crit = (2 / (this%freestream%gamma * M_inf_selected**2)) &
                    * (((numerator/denominator) ** (this%freestream%gamma / (this%freestream%gamma - 1))) - 1)
    
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
        
        integer :: i, stat

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
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:),allocatable :: pressures
        
        integer :: i, j, stat

        ! Allocate force storage
        allocate(body%dC_f(3,this%N_cells), stat=stat)
        call check_allocation(stat, "cell forces")

        ! Calculate total forces
        !$OMP parallel do schedule(static)
        do i=1,body%N_panels

            ! Discrete force coefficient acting on panel
            if (doublet_order == 2) then
                do j=1,4
                    body%dC_f(:,4*(i-1)+j) = -pressures(4*(i-1)+j)*0.25*body%panels(i)%A*body%panels(i)%n_g
                end do
            else
                body%dC_f(:,i) = -pressures(i)*body%panels(i)%A*body%panels(i)%n_g
            end if

            ! Mirror
            if (body%asym_flow) then
                if (doublet_order == 2) then
                    do j=1,4
                        body%dC_f(:,4*(i-1)+j+this%N_cells/2) = &
                            -pressures(4*(i-1)+j+this%N_cells/2)*0.25*body%panels(i)%A*body%panels(i)%n_g_mir
                    end do
                else
                    body%dC_f(:,i+this%N_cells/2) = -pressures(i+this%N_cells/2)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if
            end if

        end do
    
    end subroutine panel_solver_calc_forces_with_pressure


    subroutine panel_solver_update_report(this, p_json, body)
        ! Updates the report JSON with the information relevant to the solver

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(inout) :: p_json
        type(surface_mesh),intent(in) :: body

        type(json_value),pointer :: p_parent, p_child
        integer :: i_unit

        ! Write solver results
        call json_value_create(p_parent)
        call to_object(p_parent, 'solver_results')
        call json_value_add(p_json, p_parent)
        if (this%solver_iterations > -1) call json_value_add(p_parent, 'iterations', this%solver_iterations)
        call json_value_create(p_child)
        call to_object(p_child, 'residual')
        call json_value_add(p_parent, p_child)
        call json_value_add(p_child, 'max', this%max_res)
        call json_value_add(p_child, 'norm', this%norm_res)
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
        real,dimension(:,:),allocatable :: points
        character(len=200) :: dummy_read
        real,dimension(:),allocatable :: phi_s, phi_d

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

        ! Calculate potentials
        !$OMP parallel do
        do i=1,N_points

            ! Get induced potentials
            call body%get_induced_potentials_at_point(points(:,i), this%freestream, phi_d(i), phi_s(i))

        end do
        phi_d = phi_d*this%freestream%U
        phi_s = phi_s*this%freestream%U

        ! Delete old output file
        call delete_file(points_output_file)

        ! Open output file
        100 format(e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12, ',', e20.12)
        open(newunit=unit, file=points_output_file)

        ! Write header
        write(unit,*) 'x,y,z,phi_inf,phi_d,phi_s,phi,Phi'

        ! Write potentials out to file
        do i=1,N_points

            ! Calculate freestream potential
            phi_inf = this%freestream%U*inner(points(:,i), this%freestream%c_hat_g)

            ! Write to file
            write(unit,100) points(1,i), points(2,i), points(3,i), phi_inf, phi_d(i), phi_s(i), phi_d(i) + phi_s(i), &
                            phi_inf + phi_d(i) + phi_s(i)

        end do

        close(unit)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_export_off_body_points


end module panel_solver_mod