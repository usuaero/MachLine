module panel_solver_mod

    use helpers_mod
    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use vertex_mod
    use surface_mesh_mod
    use flow_mod
    use math_mod
    use linalg_mod
    use sort_mod

    implicit none


    type panel_solver


        character(len=:),allocatable :: formulation, pressure_for_forces
        logical :: incompressible_rule, isentropic_rule, second_order_rule
        logical,dimension(:,:),allocatable :: verts_in_dod, wake_verts_in_dod
        logical,dimension(:,:),allocatable :: edges_in_dod, wake_edges_in_dod
        type(dod),dimension(:,:),allocatable :: dod_info, wake_dod_info
        type(flow) :: freestream
        real :: norm_res, max_res
        real,dimension(3) :: C_F
        real,dimension(:,:),allocatable :: A
        real,dimension(:), allocatable :: b
        integer :: N, wake_start, N_pressures
        integer,dimension(:),allocatable :: i_cp_sorted

        contains

            procedure :: init => panel_solver_init
            procedure :: init_dirichlet => panel_solver_init_dirichlet
            procedure :: calc_domains_of_dependence => panel_solver_calc_domains_of_dependence
            procedure :: sort_control_points => panel_solver_sort_control_points
            procedure :: solve => panel_solver_solve
            procedure :: calc_source_strengths => panel_solver_calc_source_strengths
            procedure :: calc_body_influences => panel_solver_calc_body_influences
            procedure :: calc_wake_influences => panel_solver_calc_wake_influences
            procedure :: solve_system => panel_solver_solve_system
            procedure :: calc_velocities => panel_solver_calc_velocities
            procedure :: calc_pressures => panel_solver_calc_pressures
            procedure :: calc_forces => panel_solver_calc_forces
            procedure :: write_report => panel_solver_write_report

    end type panel_solver


contains


    subroutine panel_solver_init(this, solver_settings, processing_settings, body, freestream, control_point_file)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: solver_settings, processing_settings
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(inout) :: freestream
        character(len=:),allocatable,intent(in) :: control_point_file

        integer :: i, j
        real,dimension(:),allocatable :: cp_indices
        type(vtk_out) :: cp_vtk

        ! Get solver_settings
        call json_xtnsn_get(solver_settings, 'formulation', this%formulation, 'morino')
        call json_xtnsn_get(solver_settings, 'influence_calculations', influence_calc_type, 'johnson')

        ! Get pressure rules
        if (freestream%M_inf > 0.) then
            call json_xtnsn_get(processing_settings, 'pressure_rules.incompressible', this%incompressible_rule, .false.)
            call json_xtnsn_get(processing_settings, 'pressure_rules.isentropic', this%isentropic_rule, .true.)
        else
            call json_xtnsn_get(processing_settings, 'pressure_rules.incompressible', this%incompressible_rule, .true.)
            this%isentropic_rule = .false.
        end if
        call json_xtnsn_get(processing_settings, 'pressure_rules.second-order', this%second_order_rule, .false.)

        ! Get which pressure rule will be used for force calculation
        if (this%incompressible_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'incompressible')
        else if (this%isentropic_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'isentropic')
        else if (this%second_order_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'second-order')
        end if

        ! Check the force calculation pressure rule is available
        if (this%pressure_for_forces == 'incompressible' .and. .not. this%incompressible_rule) then
            write(*,*) "!!! Incompressible pressure rule is not available for force calculation. Quitting..."
            stop
        end if
        if (this%pressure_for_forces == 'isentropic' .and. .not. this%isentropic_rule) then
            write(*,*) "!!! Isentropic pressure rule is not available for force calculation. Quitting..."
            stop
        end if
        if (this%pressure_for_forces == 'second-order' .and. .not. this%second_order_rule) then
            write(*,*) "!!! Second-order pressure rule is not available for force calculation. Quitting..."
            stop
        end if

        write(*,*) "    The ", this%pressure_for_forces, " pressure rule will be used for force calculation."

        ! Store
        this%freestream = freestream

        ! Initialize based on formulation
        if (this%formulation == 'morino' .or. this%formulation == 'source-free') then
            call this%init_dirichlet(solver_settings, body)
        end if
        
        ! Write out control point geometry
        if (control_point_file /= 'none') then

            ! Clear old file
            call delete_file(control_point_file)

            ! Get indices
            allocate(cp_indices(size(body%cp)/3))
            do i=1,size(cp_indices)
                cp_indices(i) = i
            end do

            ! Write out points
            call cp_vtk%begin(control_point_file)
            call cp_vtk%write_points(body%cp)
            call cp_vtk%write_vertices(body%cp)
            call cp_vtk%write_point_scalars(cp_indices, 'index')
            call cp_vtk%finish()

        end if

        ! Calculate domains of dependence
        call this%calc_domains_of_dependence(body)

    end subroutine panel_solver_init


    subroutine panel_solver_calc_domains_of_dependence(this, body)
        ! Determines the domains of dependence for each control point based on the freestream condition

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j
        real,dimension(3) :: vert_loc

        ! For asymmetric flow on a mirrored mesh, all domains of dependence must be calculated. There are no shortcuts.
        ! For symmetric flow on a mirrored mesh, domains of dependence will be the same between mirrored panels and mirrored
        ! control points. So, we just need to calculate the DoD for mirrored control points, and then we're good.

        write(*,'(a)',advance='no') "     Calculating domains of dependence..."

        ! Allocate arrays for domain of dependence information for the body
        if (body%mirrored) then
            allocate(this%dod_info(2*body%N_panels, this%N))
            allocate(this%verts_in_dod(2*body%N_verts, this%N))
            allocate(this%edges_in_dod(2*body%N_edges, this%N))
        else
            allocate(this%dod_info(body%N_panels, this%N))
            allocate(this%verts_in_dod(body%N_verts, this%N))
            allocate(this%edges_in_dod(body%N_edges, this%N))
        end if

        ! Allocate arrays for domain of dependence information for the wake
        if (body%mirrored .and. .not. body%asym_flow) then
            allocate(this%wake_dod_info(2*body%wake%N_panels, this%N))
            allocate(this%wake_verts_in_dod(2*body%wake%N_verts, this%N))
        else
            allocate(this%wake_dod_info(body%wake%N_panels, this%N))
            allocate(this%wake_verts_in_dod(body%wake%N_verts, this%N))
        end if

        ! Loop through control points
        do j=1,body%N_cp

            ! Loop through vertices
            do i=1,body%N_verts

                vert_loc = body%vertices(i)%loc

                ! Original vertex and original control point
                this%verts_in_dod(i,j) = this%freestream%point_in_dod(vert_loc, body%cp(j,:))

                if (body%mirrored) then

                    ! Mirrored vertex and original control point
                    this%verts_in_dod(i+body%N_verts,j) = this%freestream%point_in_dod( &
                                                          mirror_about_plane(vert_loc, body%mirror_plane), &
                                                          body%cp(j,:))

                    if (body%asym_flow) then

                        ! Original vertex and mirrored control point
                        this%verts_in_dod(i,j+body%N_cp) = this%freestream%point_in_dod(vert_loc, body%cp_mirrored(j,:))

                        ! Mirrored vertex and mirrored control point
                        this%verts_in_dod(i+body%N_verts,j+body%N_cp) = this%freestream%point_in_dod( &
                                                                        mirror_about_plane(vert_loc, body%mirror_plane), &
                                                                        body%cp_mirrored(j,:))

                    end if
                end if
            end do

            ! Loop through edges

            ! Loop through body panels
            do i=1,body%N_panels

                ! Original panel and original control point
                this%dod_info(i,j) = body%panels(i)%check_dod(body%cp(j,:), this%freestream)

                if (body%mirrored) then

                    ! Check DoD for mirrored panel and original control point
                    this%dod_info(i+body%N_panels,j) = body%panels(i)%check_dod(body%cp(j,:), this%freestream, &
                                                                                .true., body%mirror_plane)
                    
                    if (body%asym_flow) then

                        ! Check DoD for original panel and mirrored control point
                        this%dod_info(i,j+body%N_cp) = body%panels(i)%check_dod(body%cp_mirrored(j,:), this%freestream)

                        ! Check DoD for mirrored panel and mirrored control point
                        this%dod_info(i+body%N_panels,j+body%N_cp) = body%panels(i)%check_dod(body%cp_mirrored(j,:), &
                                                                                              this%freestream, .true., &
                                                                                              body%mirror_plane)

                    end if
                end if
            end do

            ! Loop through wake panels
            do i=1,body%wake%N_panels

                ! Check DoD for panel and original control point
                this%wake_dod_info(i,j) = body%wake%panels(i)%check_dod(body%cp(j,:), this%freestream)

                if (body%mirrored) then

                    if (body%asym_flow) then

                        ! Check DoD for panel and mirrored control point
                        this%wake_dod_info(i,j+body%N_verts) = body%wake%panels(i)%check_dod(body%cp_mirrored(j,:), this%freestream)

                    else

                        ! Check DoD for mirrored panel and original control point
                        this%wake_dod_info(i+body%wake%N_panels,j) = body%wake%panels(i)%check_dod(body%cp(j,:), this%freestream, &
                                                                                                   .true., body%mirror_plane)

                    end if
                end if
            end do

        end do

        write(*,*) "Done"
    
    end subroutine panel_solver_calc_domains_of_dependence


    subroutine panel_solver_init_dirichlet(this, solver_settings, body)
        ! Initializes the solver to use one of the Dirichlet formulations

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: solver_settings
        type(surface_mesh),intent(inout) :: body

        real :: offset
        
        ! Place control points
        write(*,'(a)',advance='no') "     Placing control points..."

        ! Get offset
        call json_xtnsn_get(solver_settings, 'control_point_offset', offset, 1e-5)

        ! Place control points inside the body
        if (this%formulation == 'morino' .or. this%formulation == 'source-free') then
            call body%place_interior_control_points(offset)
        end if

        ! Determine size of linear system
        if (body%mirrored .and. body%asym_flow) then
            this%N = body%N_cp*2
        else
            this%N = body%N_cp
        end if

        ! Sort control points
        if (this%freestream%supersonic) then
            call this%sort_control_points(body)
        end if

        write(*,'(a, i5, a)') "Done. Placed", body%N_cp, " control points."
    
    end subroutine panel_solver_init_dirichlet


    subroutine panel_solver_sort_control_points(this, body)
        ! Sorts the control points in the compressibility direction to allow for fast matrix solving

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        real,dimension(:),allocatable :: x
        integer :: i

        ! Allocate the compressibility distance array
        allocate(x(this%N))

        ! Add original control points
        do i=1,body%N_cp
            x(i) = -inner(this%freestream%c_hat_g, body%cp(i,:))
        end do

        ! Add mirrored control points
        if (body%mirrored .and. body%asym_flow) then
            do i=1,body%N_cp
                x(i+body%N_cp) = -inner(this%freestream%c_hat_g, body%cp_mirrored(i,:))
            end do
        end if

        ! Get sorted indices
        call insertion_sort_indices(x, this%i_cp_sorted)
    
    end subroutine panel_solver_sort_control_points


    subroutine panel_solver_solve(this, body, report_file)
        ! Calls the relevant subroutine to solve the case based on the selected formulation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        character(len=:),allocatable :: report_file

        ! Calculate source strengths
        call this%calc_source_strengths(body)

        ! Calculate body influences
        call this%calc_body_influences(body)

        ! Calculate wake influences
        call this%calc_wake_influences(body)

        ! Solve the linear system
        call this%solve_system(body)

        ! Calculate velocities
        call this%calc_velocities(body)

        ! Calculate pressures
        call this%calc_pressures(body)

        ! Calculate forces
        call this%calc_forces(body)

        ! Write report
        if (report_file /= 'none') then
            call this%write_report(body, report_file)
        end if

    end subroutine panel_solver_solve


    subroutine panel_solver_calc_source_strengths(this, body)
        ! Calculates the necessary source strengths

        implicit none

        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: N_sigma, i, stat
        real,dimension(3) :: n_mirrored

        ! Set source strengths
        if (source_order == 0) then

            ! Determine necessary number of source strengths
            if (body%mirrored .and. body%asym_flow) then
                N_sigma = body%N_panels*2
            else
                N_sigma = body%N_panels
            end if

            ! Allocate source strength array
            allocate(body%sigma(N_sigma), source=0., stat=stat)
            call check_allocation(stat, "source strength array")

            ! Morino formulation
            if (this%formulation == "morino") then

                write(*,'(a)',advance='no') "     Calculating source strengths..."

                ! Loop through panels
                do i=1,body%N_panels

                    ! Existing panels
                    body%sigma(i) = -inner(body%panels(i)%normal, this%freestream%c_hat_g)

                    ! Mirrored panels for asymmetric flow
                    if (body%mirrored .and. body%asym_flow) then

                        ! Get mirrored normal vector
                        n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)

                        ! Calculate source strength
                        body%sigma(i+body%N_panels) = -inner(n_mirrored, this%freestream%c_hat_g)

                    end if
                end do

                write(*,*) "Done."

            end if
        end if
    
    end subroutine panel_solver_calc_source_strengths


    subroutine panel_solver_calc_body_influences(this, body)
        ! Calculates the influence of the body on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k, stat
        real,dimension(:),allocatable :: source_inf, doublet_inf
        integer,dimension(:),allocatable :: i_vert_s, i_vert_d
        logical :: morino

        ! Determine formulation
        morino = this%formulation == 'morino'

        ! Allocate space for inner potential calculations
        allocate(body%phi_cp_sigma(this%N), source=0., stat=stat)
        call check_allocation(stat, "induced potential vector")

        ! Allocate AIC matrix
        allocate(this%A(this%N, this%N), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")

        ! Allocate b vector
        allocate(this%b(this%N), source=0., stat=stat)
        call check_allocation(stat, "b vector")

        write(*,'(a)',advance='no') "     Calculating body influences..."

        ! Calculate source and doublet influences from body
        do i=1,body%N_cp
            do j=1,body%N_panels

                ! Calculate influence of existing panel on existing control point
                call body%panels(j)%calc_potentials(body%cp(i,:), this%freestream, this%dod_info(j,i), .false., &
                                                    source_inf, doublet_inf, i_vert_s, i_vert_d)

                ! Add influence of existing panel on existing control point
                if (morino) then
                    if (source_order == 0) then
                        body%phi_cp_sigma(i) = body%phi_cp_sigma(i) + source_inf(1)*body%sigma(j)
                    end if
                end if

                if (doublet_order == 1) then
                    do k=1,size(i_vert_d)
                        this%A(i,i_vert_d(k)) = this%A(i,i_vert_d(k)) + doublet_inf(k)
                    end do
                end if

                if (body%mirrored) then

                    if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then

                        ! Recalculate influence of mirrored panel on mirrored control point for compressible flow
                        ! If the flow is symmetric or incompressible, this will be the same as already calculated
                        if (.not. this%freestream%incompressible) then
                            call body%panels(j)%calc_potentials(body%cp_mirrored(i,:), this%freestream, &
                                                                this%dod_info(j+body%N_panels,i+body%N_cp), .true., &
                                                                source_inf, doublet_inf, i_vert_s, i_vert_d)
                        end if

                        ! Add influence of mirrored panel on mirrored control point
                        if (morino) then
                            if (source_order == 0) then
                                body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) &
                                                                 + source_inf(1)*body%sigma(j+body%N_panels)
                            end if
                        end if

                        if (doublet_order == 1) then
                            do k=1,size(i_vert_d)
                                this%A(i+body%N_cp,i_vert_d(k)+body%N_cp) = this%A(i+body%N_cp,i_vert_d(k)+body%N_cp) &
                                                                                 + doublet_inf(k)
                            end do
                        end if

                    end if

                    if (body%asym_flow) then

                        ! Calculate influence of existing panel on mirrored control point
                        if (body%vertices(i)%mirrored_is_unique .or. this%freestream%incompressible) then
                            call body%panels(j)%calc_potentials(body%cp_mirrored(i,:), this%freestream, &
                                                                this%dod_info(j,i+body%N_cp), .false., &
                                                                source_inf, doublet_inf, i_vert_s, i_vert_d)
                        end if

                        ! Add influence of existing panel on mirrored control point
                        if (body%vertices(i)%mirrored_is_unique) then

                            if (morino) then
                                if (source_order == 0) then
                                    body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) + source_inf(1)*body%sigma(j)
                                end if
                            end if

                            if (doublet_order == 1) then
                                do k=1,size(i_vert_d)
                                    this%A(i+body%N_cp,i_vert_d(k)) = this%A(i+body%N_cp,i_vert_d(k)) + doublet_inf(k)
                                end do
                            end if

                        end if

                        ! Recalculate mirrored->existing influences for compressible flow
                        if (.not. this%freestream%incompressible) then
                            call body%panels(j)%calc_potentials(body%cp(i,:), this%freestream, &
                                                                this%dod_info(j+body%N_panels,i), .true., source_inf, &
                                                                doublet_inf, i_vert_s, i_vert_d)
                        end if

                        ! Add influence of mirrored panel on existing control point
                        if (morino) then
                            if (source_order == 0) then
                                body%phi_cp_sigma(i) = body%phi_cp_sigma(i) + source_inf(1)*body%sigma(j+body%N_panels)
                            end if
                        end if

                        if (doublet_order == 1) then
                            do k=1,size(i_vert_d)
                                this%A(i,i_vert_d(k)+body%N_cp) = this%A(i,i_vert_d(k)+body%N_cp) + doublet_inf(k)
                            end do
                        end if

                    else

                        ! Calculate influence of existing panel on mirrored control point
                        ! This is the same as the influence of the mirrored panel on the existing control point,
                        ! even for compressible flow, since we know the flow is symmetric here
                        call body%panels(j)%calc_potentials(body%cp_mirrored(i,:), this%freestream, this%dod_info(j,i), &
                                                            .false., source_inf, doublet_inf, i_vert_s, i_vert_d)

                        ! Add influence of mirrored panel on existing control point
                        if (morino) then
                            if (source_order == 0) then
                                body%phi_cp_sigma(i) = body%phi_cp_sigma(i) + source_inf(1)*body%sigma(j)
                            end if
                        end if

                        if (doublet_order == 1) then
                            do k=1,size(i_vert_d)
                                this%A(i,i_vert_d(k)) = this%A(i,i_vert_d(k)) + doublet_inf(k)
                            end do
                        end if

                    end if

                end if

            end do

            ! Enforce doublet strength matching (i.e. for non-unique, mirrored control points, the
            ! doublet strengths must be the same). The RHS for these rows should still be zero.
            if (body%mirrored .and. body%asym_flow) then
                if (.not. body%vertices(i)%mirrored_is_unique) then
                    this%A(i+body%N_cp,i) = 1.
                    this%A(i+body%N_cp,i+body%N_cp) = -1.

                ! If the control point is unique, it's target potential will need to be set for the source-free formulation
                else if (.not. morino) then
                    this%b(i+body%N_cp) = -inner(body%cp_mirrored(i,:), this%freestream%c_hat_g)
                end if
            end if

            ! Set target potential for source-free formulation
            if (.not. morino) then
                this%b(i) = -inner(body%cp(i,:), this%freestream%c_hat_g)
            end if

        end do

        write(*,*) "Done."
        if (any(isnan(this%A))) then
            write(*,*)
            write(*,*) "!!! NaN in A after body computations. Quitting..."
            stop
        end if
    
    end subroutine panel_solver_calc_body_influences


    subroutine panel_solver_calc_wake_influences(this, body)
        ! Calculates the influence of the wake on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k
        real,dimension(:),allocatable ::  doublet_inf, source_inf
        integer,dimension(:),allocatable :: i_vert_d, i_vert_s

        ! Calculate influence of wake
        if (body%wake%N_panels > 0) then
            write(*,'(a)',advance='no') "     Calculating wake influences..."

            ! Loop through control points
            do i=1,body%N_cp

                ! Get doublet influence from wake
                ! Note that for the wake, in the case of mirrored meshes with asymmetric flow, the mirrored wake panels have actually been created.
                ! In this case, there are technically no mirrored panels, and this loop will cycle through both existing and mirrored panels.
                ! For symmetric flow, mirrored panels still need to be added as before.
                do j=1,body%wake%N_panels

                    ! Caclulate influence of existing panel on existing control point
                    call body%wake%panels(j)%calc_potentials(body%cp(i,:), this%freestream, &
                                                             this%wake_dod_info(j,i), .false., &
                                                             source_inf, doublet_inf, i_vert_s, i_vert_d)

                    ! Add doublet influence
                    if (doublet_order == 1) then
                        do k=1,size(i_vert_d)
                            this%A(i,i_vert_d(k)) = this%A(i,i_vert_d(k)) + doublet_inf(k)
                        end do
                    end if

                    ! Get influence on mirrored control point
                    if (body%mirrored) then

                        if (body%asym_flow) then

                            ! Calculate influence of existing panel on mirrored point
                            call body%wake%panels(j)%calc_potentials(body%cp_mirrored(i,:), this%freestream, &
                                                                     this%wake_dod_info(j,i+body%N_cp), .false., &
                                                                     source_inf, doublet_inf, i_vert_s, i_vert_d)

                            ! Add doublet influence
                            if (body%vertices(i)%mirrored_is_unique) then
                                if (doublet_order == 1) then
                                    do k=1,size(i_vert_d)
                                        this%A(i+body%N_cp,i_vert_d(k)) = this%A(i+body%N_cp,i_vert_d(k)) + doublet_inf(k)
                                    end do
                                end if
                            end if

                        else

                            ! Calculate influence of mirrored panel on existing control point
                            call body%wake%panels(j)%calc_potentials(body%cp_mirrored(i,:), this%freestream, &
                                                                     this%wake_dod_info(j+body%wake%N_panels,i), .true., &
                                                                     source_inf, doublet_inf, i_vert_s, i_vert_d)

                            ! Add doublet influence
                            if (doublet_order == 1) then
                                do k=1,size(i_vert_d)
                                    this%A(i,i_vert_d(k)) = this%A(i,i_vert_d(k)) + doublet_inf(k)
                                end do
                            end if

                        end if

                    end if
                end do
            end do

            write(*,*) "Done."

        end if

    end subroutine panel_solver_calc_wake_influences


    subroutine panel_solver_solve_system(this, body)
        ! Solves the linear system for the singularity strengths

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        real,dimension(:,:),allocatable :: A_copy
        integer :: stat, i, j

        write(*,'(a)',advance='no') "     Solving linear system..."

        ! Check for NaNs; I'd rather have it fail here than give the user garbage results
        if (any(isnan(this%A))) then
            write(*,*) "!!! Invalid value detected in A matrix. Quitting..."
            stop
        end if
        if (any(isnan(this%b))) then
            write(*,*) "!!! Invalid value detected in b vector. Quitting..."
            stop
        end if

        ! Make a copy of A (lu_solve replaces A with its decomposition)
        allocate(A_copy, source=this%A, stat=stat)
        call check_allocation(stat, "solver copy of AIC matrix")

        ! Check for uninfluenced/ing points
        do i=1,this%N
            if (all(this%A(i,:) == 0.)) then
                write(*,*) "!!! Control point ", i, " is not influenced. Quitting"
                stop
            end if
            if (all(this%A(:,i) == 0.)) then
                write(*,*) "!!! Vertex ", i, " exerts no influence. Quitting"
                stop
            end if
        end do

        ! Write A and b to file
        !open(34, file="./dev/A_mat.txt")
        !do i=1,this%N
        !    write(34,*) this%A(i,:)
        !end do
        !close(34)
        !open(34, file="./dev/b_vec.txt")
        !do i=1,this%N
        !    write(34,*) this%b(i)
        !end do
        !close(34)

        ! Set b vector for Morino formulation
        if (this%formulation == "morino") then
            this%b = -body%phi_cp_sigma
        end if

        ! Solve
        call lu_solve(this%N, A_copy, this%b, body%mu)
        write(*,*) "Done."

        ! Clean up memory
        deallocate(A_copy)

        ! Calculate potential at control points
        body%phi_cp_mu = matmul(this%A, body%mu)
        body%phi_cp = body%phi_cp_mu+body%phi_cp_sigma

        ! Calculate residual parameters
        this%max_res = maxval(abs(body%phi_cp_mu-this%b))
        this%norm_res = sqrt(sum((body%phi_cp_mu-this%b)**2))
        write(*,*) "        Maximum residual:", this%max_res
        write(*,*) "        Norm of residual:", this%norm_res

    end subroutine panel_solver_solve_system


    subroutine panel_solver_calc_velocities(this, body)
        ! Calculates the surface velocities

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, N_vels, stat

        write(*,'(a)',advance='no') "     Calculating surface velocities..."

        ! Determine surface velocities
        if (body%mirrored .and. body%asym_flow) then

            ! Allocate velocity storage
            N_vels = body%N_panels*2
            allocate(body%V(N_vels,3), stat=stat)
            call check_allocation(stat, "surface velocity vectors")

            ! Calculate the surface velocity on each panel
            do i=1,body%N_panels

                if (this%formulation == "morino") then

                    ! Original panel
                    body%V(i,:) = this%freestream%U*(this%freestream%c_hat_g + body%panels(i)%get_velocity_jump(body%mu, &
                                  body%sigma, .false., body%mirror_plane))

                    ! Mirror
                    body%V(i+body%N_panels,:) = this%freestream%U*(this%freestream%c_hat_g + &
                                                body%panels(i)%get_velocity_jump(body%mu, body%sigma, .true., body%mirror_plane))
                
                else

                    ! Original panel
                    body%V(i,:) = this%freestream%U*body%panels(i)%get_velocity_jump(body%mu, body%sigma, &
                                                                                     .false., body%mirror_plane)

                    ! Mirror
                    body%V(i+body%N_panels,:) = this%freestream%U*body%panels(i)%get_velocity_jump(body%mu, body%sigma, &
                                                                                              .true., body%mirror_plane)

                end if

            end do

        else

            ! Allocate velocity storage
            N_vels = body%N_panels
            allocate(body%V(N_vels,3), stat=stat)
            call check_allocation(stat, "surface velocity vectors")

            ! Calculate the surface velocity on each panel
            if (this%formulation == "morino") then
                do i=1,body%N_panels
                    body%V(i,:) = this%freestream%U*(this%freestream%c_hat_g &
                                  + body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., 0))
                end do
            else
                do i=1,body%N_panels
                    body%V(i,:) = this%freestream%U*body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., 0)
                end do
            end if

        end if
        write(*,*) "Done."

    end subroutine panel_solver_calc_velocities


    subroutine panel_solver_calc_pressures(this, body)
        ! Calculates the surface pressures

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: V_pert
        real :: a, b, c, C_p_vac

        write(*,'(a)',advance='no') "     Calculating surface pressures..."
        this%N_pressures = size(body%V)/3

        ! Calculate vacuum pressure coefficient
        C_p_vac = -2./(this%freestream%gamma*this%freestream%M_inf**2)

        ! Incompressible rule
        if (this%incompressible_rule) then

            ! Allocate storage
            allocate(body%C_p_inc(this%N_pressures), stat=stat)
            call check_allocation(stat, "incompressible surface pressures")

            ! Calculate
            do i=1,this%N_pressures
                body%C_p_inc(i) = 1.-inner(body%V(i,:), body%V(i,:))*this%freestream%U_inv**2
            end do

        end if
        
        ! Isentropic rule
        if (this%isentropic_rule) then

            ! Calculate isentropic pressure correction terms
            a = 2./(this%freestream%gamma*this%freestream%M_inf**2)
            b = 0.5*(this%freestream%gamma-1.)*this%freestream%M_inf**2
            c = this%freestream%gamma/(this%freestream%gamma-1.)

            ! Allocate storage
            allocate(body%C_p_ise(this%N_pressures), stat=stat)
            call check_allocation(stat, "isentropic surface pressures")

            ! Calculate
            do i=1,this%N_pressures
                
                ! Incompressible first
                body%C_p_ise(i) = 1.-inner(body%V(i,:), body%V(i,:))*this%freestream%U_inv**2

                ! Apply compressible correction
                body%C_p_ise(i) = a*((1.+b*body%C_p_ise(i))**c - 1.)

                ! Check for NaN and replace with vacuum pressure
                if (isnan(body%C_p_ise(i))) then
                    body%C_p_ise(i) = C_p_vac
                end if

            end do

        end if
        
        ! Second-order rule
        if (this%second_order_rule) then

            ! Allocate storage
            allocate(body%C_p_2nd(this%N_pressures), stat=stat)
            call check_allocation(stat, "second-order surface pressures")

            ! Calculate (E&M Eq. (N..2.43))
            do i=1,this%N_pressures

                ! Get perturbation velocity in the compressible frame
                V_pert = matmul(this%freestream%A_g_to_c, body%V(i,:)-this%freestream%v_inf)

                ! Calculate first term
                body%C_p_2nd(i) = -2.*V_pert(1)*this%freestream%U_inv

                ! Calculate second term
                body%C_p_2nd(i) = body%C_p_2nd(i) &
                                  - ((1.-this%freestream%M_inf**2)*V_pert(1)**2 + V_pert(2)**2 + V_pert(3)**2) &
                                  *this%freestream%U_inv**2
            end do

        end if

        ! TODO: IMPLEMENT SUBSONIC PRESSURE CORRECTIONS

        write(*,*) "Done."
        
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

        if (this%freestream%M_inf > 0.) then
            write(*,*) "        Vacuum pressure coefficient:", C_p_vac
        end if

    end subroutine panel_solver_calc_pressures


    subroutine panel_solver_calc_forces(this, body)
        ! Calculates the forces

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: n_mirrored

        write(*,'(a)',advance='no') "     Calculating forces..."

        ! Allocate force storage
        allocate(body%dC_f(this%N_pressures,3), stat=stat)
        call check_allocation(stat, "forces")

        ! Calculate total forces
        do i=1,body%N_panels

            select case (this%pressure_for_forces)

            case ('incompressible')

                ! Discrete force coefficient acting on panel
                body%dC_f(i,:) = -body%C_p_inc(i)*body%panels(i)%A*body%panels(i)%normal

                ! Mirror
                if (body%mirrored .and. body%asym_flow) then
                    n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)
                    body%dC_f(i+body%N_panels,:) = -body%C_p_inc(i+body%N_panels)*body%panels(i)%A*n_mirrored
                end if

            case ('isentropic')

                ! Discrete force coefficient acting on panel
                body%dC_f(i,:) = -body%C_p_ise(i)*body%panels(i)%A*body%panels(i)%normal

                ! Mirror
                if (body%mirrored .and. body%asym_flow) then
                    n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)
                    body%dC_f(i+body%N_panels,:) = -body%C_p_ise(i+body%N_panels)*body%panels(i)%A*n_mirrored
                end if

            case ('second-order')

                ! Discrete force coefficient acting on panel
                body%dC_f(i,:) = -body%C_p_2nd(i)*body%panels(i)%A*body%panels(i)%normal

                ! Mirror
                if (body%mirrored .and. body%asym_flow) then
                    n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)
                    body%dC_f(i+body%N_panels,:) = -body%C_p_2nd(i+body%N_panels)*body%panels(i)%A*n_mirrored
                end if

            end select

        end do

        ! Sum discrete forces
        this%C_F(:) = sum(body%dC_f, dim=1)/body%S_ref

        ! Add contributions from mirrored half
        if (body%mirrored .and. .not. body%asym_flow) then
            this%C_F = 2.*this%C_F
            this%C_F(body%mirror_plane) = 0. ! We know this
        end if

        write(*,*) "Done."
        write(*,*) "        Cx:", this%C_F(1)
        write(*,*) "        Cy:", this%C_F(2)
        write(*,*) "        Cz:", this%C_F(3)
    
    end subroutine panel_solver_calc_forces


    subroutine panel_solver_write_report(this, body, report_file)
        ! Writes the report file

        implicit none

        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body
        character(len=:),allocatable :: report_file

        open(12, file=report_file)

        ! Header
        write(12,'(a)') "MachLine Report (c) 2022 USU AeroLab"

        ! Solver results
        write(12,*) "Maximum residual:", this%max_res
        write(12,*) "Norm of residual:", this%norm_res
        
        if (this%incompressible_rule) then
            write(12,*) "Maximum incompressible pressure coefficient:", maxval(body%C_p_inc)
            write(12,*) "Minimum incompressible pressure coefficient:", minval(body%C_p_inc)
        end if
        
        if (this%isentropic_rule) then
            write(12,*) "Maximum isentropic pressure coefficient:", maxval(body%C_p_ise)
            write(12,*) "Minimum isentropic pressure coefficient:", minval(body%C_p_ise)
        end if
        
        if (this%second_order_rule) then
            write(12,*) "Maximum second-order pressure coefficient:", maxval(body%C_p_2nd)
            write(12,*) "Minimum second-order pressure coefficient:", minval(body%C_p_2nd)
        end if

        write(12,*) "Cx:", this%C_F(1)
        write(12,*) "Cy:", this%C_F(2)
        write(12,*) "Cz:", this%C_F(3)

        close(12)
   
    end subroutine panel_solver_write_report


end module panel_solver_mod