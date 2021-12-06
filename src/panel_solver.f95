module panel_solver_mod

    use helpers_mod
    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use vertex_mod
    use surface_mesh_mod
    use flow_mod
    use math_mod

    implicit none


    type panel_solver


        character(len=:),allocatable :: formulation

        contains

            procedure :: init => panel_solver_init
            procedure :: init_morino => panel_solver_init_morino
            procedure :: init_source_free => panel_solver_init_source_free
            procedure :: solve => panel_solver_solve
            procedure :: solve_morino => panel_solver_solve_morino
            procedure :: solve_source_free => panel_solver_solve_source_free

    end type panel_solver


contains


    subroutine panel_solver_init(this, settings, body)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body

        ! Get settings
        call json_xtnsn_get(settings, 'formulation', this%formulation, 'morino')

        ! Initialize based on formulation
        if (this%formulation == 'morino') then
            call this%init_morino(settings, body)
        else if (this%formulation == 'source-free') then
            call this%init_source_free(settings, body)
        end if

    end subroutine panel_solver_init


    subroutine panel_solver_init_morino(this, settings, body)
        ! Initializes the solver to use the Morino formulation (i.e. zero inner perturbation potential)

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body

        real :: offset
        
        write(*,'(a)',advance='no') "     Placing control points..."

        ! Get offset
        call json_xtnsn_get(settings, 'control_point_offset', offset, 1e-5)

        ! Place control points inside the body
        call body%place_interior_control_points(offset)
        write(*,*) "Done."
    
    end subroutine panel_solver_init_morino


    subroutine panel_solver_init_source_free(this, settings, body)
        ! Initializes the solver to use the source-free formulation (i.e. zero inner total potential)

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body

        real :: offset
        
        write(*,'(a)',advance='no') "     Placing control points..."

        ! Get offset
        call json_xtnsn_get(settings, 'control_point_offset', offset, 1e-5)

        ! Place control points inside the body
        call body%place_interior_control_points(offset)
        write(*,*) "Done."
    
    end subroutine panel_solver_init_source_free


    subroutine panel_solver_solve(this, body, freestream, report_file)
        ! Calls the relevant subroutine to solve the case based on the selected formulation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(in) :: freestream
        character(len=:),allocatable :: report_file

        ! Morino formulation
        if (this%formulation == 'morino') then
            call this%solve_morino(body, freestream, report_file)
        else if (this%formulation == 'source-free') then
            call this%solve_source_free(body, freestream, report_file)
        end if

    end subroutine panel_solver_solve


    subroutine panel_solver_solve_morino(this, body, freestream, report_file)
        ! Solves the Morino formulation for the given conditions

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(in) :: freestream
        character(len=:),allocatable :: report_file

        integer :: i, j, k
        real,dimension(:),allocatable :: source_inf, doublet_inf
        integer,dimension(:),allocatable :: source_verts, doublet_verts
        real,dimension(:,:),allocatable :: A, A_copy
        real,dimension(:),allocatable :: b
        integer :: stat, N_sigma, N_mu, N_pressures
        real,dimension(3) :: n_mirrored, cp_mirrored, C_F
        real,dimension(:,:),allocatable :: dC_F

        write(*,'(a)',advance='no') "     Calculating source strengths..."

        ! Set source strengths
        if (source_order == 0) then

            ! Determine necessary number of source strengths
            if (body%mirrored .and. body%asym_flow) then
                N_sigma = body%N_panels*2
            else
                N_sigma = body%N_panels
            end if

            ! Allocate source strength array
            allocate(body%sigma(N_sigma))

            ! Loop through panels
            do i=1,body%N_panels

                ! Existing panels
                body%sigma(i) = -inner(body%panels(i)%normal, freestream%c0)

                ! Mirrored panels for asymmetric flow
                if (body%mirrored .and. body%asym_flow) then

                    ! Get mirrored normal vector
                    n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)

                    ! Calculate source strength
                    body%sigma(i+body%N_panels) = -inner(n_mirrored, freestream%c0)

                end if
            end do

        end if
        write(*,*) "Done."

        ! Determine number of doublet strengths (some will be repeats for mirrored vertices)
        if (body%mirrored .and. body%asym_flow) then
            N_mu = body%N_cp*2
        else
            N_mu = body%N_cp
        end if

        ! Allocate space for inner potential calculations
        allocate(body%phi_cp_sigma(N_mu), source=0., stat=stat)
        call check_allocation(stat, "induced potential vector")

        ! Allocate AIC matrix
        allocate(A(N_mu, N_mu), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")

        write(*,'(a)',advance='no') "     Calculating body influences..."

        ! Calculate source and doublet influences from body
        do i=1,body%N_cp
            do j=1,body%N_panels

                ! Get influence for existing->existing and mirrored->mirrored
                source_inf = body%panels(j)%get_source_potential(body%control_points(i,:), source_verts)
                doublet_inf = body%panels(j)%get_doublet_potential(body%control_points(i,:), doublet_verts)

                ! Influence of existing panel on existing control point
                if (source_order == 0) then
                    body%phi_cp_sigma(i) = body%phi_cp_sigma(i) + source_inf(1)*body%sigma(j)
                end if

                if (doublet_order == 1) then
                    do k=1,size(doublet_verts)
                        A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                    end do
                end if

                ! Get influences for mirroring
                if (body%mirrored) then

                    ! Influence of mirrored panels on mirrored control points (uses the influences already calculated)
                    if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then

                        if (source_order == 0) then
                            body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) &
                                                             + source_inf(1)*body%sigma(j+body%N_panels)
                        end if

                        if (doublet_order == 1) then
                            do k=1,size(doublet_verts)
                                A(i+body%N_cp,doublet_verts(k)+body%N_cp) = A(i+body%N_cp,doublet_verts(k)+body%N_cp) &
                                                                            + doublet_inf(k)
                            end do
                        end if

                    end if

                    ! Get mirrored point (recall the Green's function is reciprocal)
                    cp_mirrored = mirror_about_plane(body%control_points(i,:), body%mirror_plane)

                    ! Calculate influences
                    source_inf = body%panels(j)%get_source_potential(cp_mirrored, source_verts)
                    doublet_inf = body%panels(j)%get_doublet_potential(cp_mirrored, doublet_verts)

                    if (body%asym_flow) then

                        ! Influence of mirrored panel on existing control point
                        if (source_order == 0) then
                            body%phi_cp_sigma(i) = body%phi_cp_sigma(i) + source_inf(1)*body%sigma(j+body%N_panels)
                        end if

                        if (doublet_order == 1) then
                            do k=1,size(doublet_verts)
                                A(i,doublet_verts(k)+body%N_cp) = A(i,doublet_verts(k)+body%N_cp) + doublet_inf(k)
                            end do
                        end if

                        ! Influence of existing panel on mirrored control point
                        if (body%vertices(i)%mirrored_is_unique) then

                            if (source_order == 0) then
                                body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) + source_inf(1)*body%sigma(j)
                            end if

                            if (doublet_order == 1) then
                                do k=1,size(doublet_verts)
                                    A(i+body%N_cp,doublet_verts(k)) = A(i+body%N_cp,doublet_verts(k)) + doublet_inf(k)
                                end do
                            end if

                        end if

                    else

                        ! Influence of mirrored panel on existing control point
                        if (source_order == 0) then
                            body%phi_cp_sigma(i) = body%phi_cp_sigma(i) + source_inf(1)*body%sigma(j)
                        end if

                        if (doublet_order == 1) then
                            do k=1,size(doublet_verts)
                                A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                            end do
                        end if

                    end if

                end if

            end do

            ! Enforce doublet strength matching (i.e. for non-unique, mirrored control points, the
            ! doublet strengths must be the same). The RHS for these rows should still be zero.
            if (body%mirrored .and. body%asym_flow .and. .not. body%vertices(i)%mirrored_is_unique) then
                A(i+body%N_cp,i) = 1.
                A(i+body%N_cp,i+body%N_cp) = -1.
            end if

        end do
        write(*,*) "Done."

        ! Calculate influence of wake
        if (body%wake%N_panels > 0) then
            write(*,'(a)',advance='no') "     Calculating wake influences..."

            ! Loop through control points
            do i=1,body%N_cp

                ! Get doublet influence from wake
                ! Note that for the wake, in the case of mirrored mesh with asymmetric flow, the mirrored wake panels have actually been created.
                ! In this case, there are technically no mirrored panels, and this loop will cycle through both existing and mirrored panels.
                ! For symmetric flow, mirrored panels still need to be added as before.
                do j=1,body%wake%N_panels

                    ! Caclulate influence
                    doublet_inf = body%wake%panels(j)%get_doublet_potential(body%control_points(i,:), doublet_verts)

                    ! Influence on existing control point
                    if (doublet_order == 1) then
                        do k=1,size(doublet_verts)
                            A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                        end do
                    end if

                    ! Get influence on mirrored control point
                    if (body%mirrored) then

                        ! Get mirrored point
                        cp_mirrored = mirror_about_plane(body%control_points(i,:), body%mirror_plane)

                        ! Calculate influences
                        doublet_inf = body%wake%panels(j)%get_doublet_potential(cp_mirrored, doublet_verts)

                        if (body%asym_flow) then

                            ! Influence on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                if (doublet_order == 1) then
                                    do k=1,size(doublet_verts)
                                        A(i+body%N_cp,doublet_verts(k)) = A(i+body%N_cp,doublet_verts(k)) + doublet_inf(k)
                                    end do
                                end if
                            end if

                        else

                            ! Influence of mirrored panel on existing control point
                            if (doublet_order == 1) then
                                do k=1,size(doublet_verts)
                                    A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                                end do
                            end if

                        end if

                    end if
                end do
            end do
            write(*,*) "Done."
        end if

        write(*,'(a)',advance='no') "     Solving linear system..."

        ! Make a copy of A (lu_solve replaces A with its decomposition)
        allocate(A_copy, source=A, stat=stat)
        call check_allocation(stat, "solver copy of AIC matrix")

        ! Solve
        b = -body%phi_cp_sigma
        call lu_solve(N_mu, A_copy, b, body%mu)
        write(*,*) "Done."

        ! Clean up memory
        deallocate(A_copy)
        deallocate(b)

        ! Calculate potential at control points
        body%phi_cp_mu = matmul(A, body%mu)
        body%phi_cp = body%phi_cp_mu+body%phi_cp_sigma
        write(*,*) "        Maximum residual inner potential:", maxval(abs(body%phi_cp))
        write(*,*) "        Norm of residual innner potential:", sqrt(sum(body%phi_cp**2))

        write(*,'(a)',advance='no') "     Calculating surface velocities and pressures..."

        ! Determine surface velocities
        if (body%mirrored .and. body%asym_flow) then

            ! Allocate velocity storage
            N_pressures = body%N_panels*2
            allocate(body%V(N_pressures,3), stat=stat)
            call check_allocation(stat, "surface velocity vectors")

            ! Calculate the surface velocity on each panel
            do i=1,body%N_panels

                ! Original panel
                body%V(i,:) = freestream%U*(freestream%c0 + body%panels(i)%get_velocity_jump(body%mu, &
                              body%sigma, .false., body%mirror_plane))

                ! Mirror
                body%V(i+body%N_panels,:) = freestream%U*(freestream%c0 + body%panels(i)%get_velocity_jump(body%mu, &
                                            body%sigma, .true., body%mirror_plane))

            end do

        else

            ! Allocate velocity storage
            N_pressures = body%N_panels
            allocate(body%V(N_pressures,3), stat=stat)
            call check_allocation(stat, "surface velocity vectors")

            ! Calculate the surface velocity on each panel
            do i=1,body%N_panels
                body%V(i,:) = freestream%U*(freestream%c0 + body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., 0))
            end do

        end if

        ! Calculate coefficients of pressure
        allocate(body%C_p(N_pressures), stat=stat)
        call check_allocation(stat, "surface pressures")
        do i=1,N_pressures
            body%C_p(i) = 1.-(norm(body%V(i,:))*freestream%U_inv)**2
        end do

        write(*,*) "Done."
        write(*,*) "        Maximum pressure coefficient:", maxval(body%C_p)
        write(*,*) "        Minimum pressure coefficient:", minval(body%C_p)

        ! Calculate total forces
        write(*,'(a)',advance='no') "     Calculating forces..."
        allocate(dC_F(N_pressures,3), stat=stat)
        call check_allocation(stat, "forces")
        do i=1,body%N_panels

            ! Discrete force coefficient acting on panel
            dC_F(i,:) = body%C_p(i)*body%panels(i)%A*body%panels(i)%normal

            ! Mirror
            if (body%mirrored .and. body%asym_flow) then
                n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)
                dC_F(i+body%N_panels,:) = body%C_p(i+body%N_panels)*body%panels(i)%A*n_mirrored
            end if
        end do

        ! Sum discrete forces
        C_F(:) = sum(dC_F, dim=1)/body%S_ref

        write(*,*) "Done."
        write(*,*) "        Cx:", C_F(1)
        write(*,*) "        Cy:", C_F(2)
        write(*,*) "        Cz:", C_F(3)

        ! Write report file
        if (report_file /= 'none') then

            open(1, file=report_file)

            ! Header
            write(1,'(a)') "TriPan Report (c) 2021 USU AeroLab"

            ! Solver results
            write(1,*) "Maximum residual inner potential:", maxval(abs(body%phi_cp))
            write(1,*) "Norm of residual innner potential:", sqrt(sum(body%phi_cp**2))
            write(1,*) "Maximum pressure coefficient:", maxval(body%C_p)
            write(1,*) "Minimum pressure coefficient:", minval(body%C_p)
            write(1,*) "Cx:", C_F(1)
            write(1,*) "Cy:", C_F(2)
            write(1,*) "Cz:", C_F(3)

            close(1)

        end if
    
    end subroutine panel_solver_solve_morino


    subroutine panel_solver_solve_source_free(this, body, freestream, report_file)
        ! Solves the source_free formulation for the given conditions

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(in) :: freestream
        character(len=:),allocatable :: report_file

        integer :: i, j, k
        real,dimension(:),allocatable :: doublet_inf
        integer,dimension(:),allocatable :: doublet_verts
        real,dimension(:,:),allocatable :: A, A_copy
        real,dimension(:),allocatable :: b
        integer :: stat, N_sigma, N_mu, N_pressures
        real,dimension(3) :: n_mirrored, cp_mirrored, C_F
        real,dimension(:,:),allocatable :: dC_F

        ! Set source strengths
        if (source_order == 0) then

            ! Determine necessary number of source strengths
            if (body%mirrored .and. body%asym_flow) then
                N_sigma = body%N_panels*2
            else
                N_sigma = body%N_panels
            end if

            ! Allocate source strength array
            allocate(body%sigma(N_sigma), source=0.)

        end if

        ! Determine number of doublet strengths (some will be repeats for mirrored vertices)
        if (body%mirrored .and. body%asym_flow) then
            N_mu = body%N_cp*2
        else
            N_mu = body%N_cp
        end if

        ! Allocate space for inner potential calculations
        allocate(body%phi_cp_sigma(N_mu), source=0., stat=stat)
        call check_allocation(stat, "induced potential vector")

        ! Allocate AIC matrix
        allocate(A(N_mu, N_mu), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")

        ! Allocate target inner potential
        allocate(b(N_mu), source=0., stat=stat)
        call check_allocation(stat, "solver copy of AIC matrix")

        write(*,'(a)',advance='no') "     Calculating body influences..."

        ! Calculate source and doublet influences from body
        do i=1,body%N_cp

            ! Get mirrored control point
            if (body%mirrored) then
                cp_mirrored = mirror_about_plane(body%control_points(i,:), body%mirror_plane)
            end if

            do j=1,body%N_panels

                ! Get influence for existing->existing and mirrored->mirrored
                doublet_inf = body%panels(j)%get_doublet_potential(body%control_points(i,:), doublet_verts)

                ! Influence of existing panel on existing control point
                if (doublet_order == 1) then
                    do k=1,size(doublet_verts)
                        A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                    end do
                end if

                ! Get influences for mirroring
                if (body%mirrored) then

                    ! Influence of mirrored panels on mirrored control points (uses the influences already calculated)
                    if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then

                        if (doublet_order == 1) then
                            do k=1,size(doublet_verts)
                                A(i+body%N_cp,doublet_verts(k)+body%N_cp) = A(i+body%N_cp,doublet_verts(k)+body%N_cp) &
                                                                            + doublet_inf(k)
                            end do
                        end if

                    end if

                    ! Calculate influence on mirrored control points (recall the Green's function is reciprocal)
                    doublet_inf = body%panels(j)%get_doublet_potential(cp_mirrored, doublet_verts)

                    if (body%asym_flow) then

                        ! Influence of mirrored panel on existing control point
                        if (doublet_order == 1) then
                            do k=1,size(doublet_verts)
                                A(i,doublet_verts(k)+body%N_cp) = A(i,doublet_verts(k)+body%N_cp) + doublet_inf(k)
                            end do
                        end if

                        ! Influence of existing panel on mirrored control point
                        if (body%vertices(i)%mirrored_is_unique) then

                            if (doublet_order == 1) then
                                do k=1,size(doublet_verts)
                                    A(i+body%N_cp,doublet_verts(k)) = A(i+body%N_cp,doublet_verts(k)) + doublet_inf(k)
                                end do
                            end if

                        end if

                    else

                        ! Influence of mirrored panel on existing control point
                        if (doublet_order == 1) then
                            do k=1,size(doublet_verts)
                                A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                            end do
                        end if

                    end if

                end if

            end do

            ! Enforce doublet strength matching (i.e. for non-unique, mirrored control points, the
            ! doublet strengths must be the same). The RHS for these rows must still be zero.
            if (body%mirrored .and. body%asym_flow) then
                if (.not. body%vertices(i)%mirrored_is_unique) then
                    A(i+body%N_cp,i) = 1.
                    A(i+body%N_cp,i+body%N_cp) = -1.

                ! If the control point is unique, it's target potential will need to be set
                else
                    b(i+body%N_cp) = -inner(cp_mirrored, freestream%c0)

                end if
            end if

            ! Set target inner potential for existing control points
            b(i) = -inner(body%control_points(i,:), freestream%c0)

        end do
        write(*,*) "Done."

        ! Calculate influence of wake
        if (body%wake%N_panels > 0) then
            write(*,'(a)',advance='no') "     Calculating wake influences..."

            ! Loop through control points
            do i=1,body%N_cp

                ! Get doublet influence from wake
                ! Note that for the wake, in the case of mirrored mesh with asymmetric flow, the mirrored wake panels have actually been created.
                ! In this case, there are technically no mirrored panels, and this loop will cycle through both existing and mirrored panels.
                ! For symmetric flow, mirrored panels still need to be added as before.
                do j=1,body%wake%N_panels

                    ! Caclulate influence
                    doublet_inf = body%wake%panels(j)%get_doublet_potential(body%control_points(i,:), doublet_verts)

                    ! Influence on existing control point
                    if (doublet_order == 1) then
                        do k=1,size(doublet_verts)
                            A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                        end do
                    end if

                    ! Get influence on mirrored control point
                    if (body%mirrored) then

                        ! Get mirrored point
                        cp_mirrored = mirror_about_plane(body%control_points(i,:), body%mirror_plane)

                        ! Calculate influences
                        doublet_inf = body%wake%panels(j)%get_doublet_potential(cp_mirrored, doublet_verts)

                        if (body%asym_flow) then

                            ! Influence on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                if (doublet_order == 1) then
                                    do k=1,size(doublet_verts)
                                        A(i+body%N_cp,doublet_verts(k)) = A(i+body%N_cp,doublet_verts(k)) + doublet_inf(k)
                                    end do
                                end if
                            end if

                        else

                            ! Influence of mirrored panel on existing control point
                            if (doublet_order == 1) then
                                do k=1,size(doublet_verts)
                                    A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                                end do
                            end if

                        end if

                    end if
                end do
            end do
            write(*,*) "Done."
        end if

        write(*,'(a)',advance='no') "     Solving linear system..."

        ! Make a copy of A (lu_solve replaces A with its decomposition)
        allocate(A_copy, source=A, stat=stat)
        call check_allocation(stat, "solver copy of AIC matrix")

        ! Solve
        call lu_solve(N_mu, A_copy, b, body%mu)
        write(*,*) "Done."

        ! Clean up memory
        deallocate(A_copy)

        ! Calculate total potential at control points
        body%phi_cp_mu = matmul(A, body%mu)
        body%phi_cp = body%phi_cp_mu
        write(*,*) "        Maximum residual inner total potential:", maxval(abs(body%phi_cp-b))
        write(*,*) "        Norm of residual innner total potential:", sqrt(sum((body%phi_cp-b)**2))
        deallocate(b)

        write(*,'(a)',advance='no') "     Calculating surface velocities and pressures..."

        ! Determine surface velocities
        if (body%mirrored .and. body%asym_flow) then

            ! Allocate velocity storage
            N_pressures = body%N_panels*2
            allocate(body%V(N_pressures,3), stat=stat)
            call check_allocation(stat, "surface velocity vectors")

            ! Calculate the surface velocity on each panel
            do i=1,body%N_panels

                ! Original panel
                body%V(i,:) = freestream%U*body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., body%mirror_plane)

                ! Mirror
                body%V(i+body%N_panels,:) = freestream%U*body%panels(i)%get_velocity_jump(body%mu, body%sigma, &
                                                                                          .true., body%mirror_plane)

            end do

        else

            ! Allocate velocity storage
            N_pressures = body%N_panels
            allocate(body%V(N_pressures,3), stat=stat)
            call check_allocation(stat, "surface velocity vectors")

            ! Calculate the surface velocity on each panel
            do i=1,body%N_panels
                body%V(i,:) = freestream%U*body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., 0)
            end do

        end if

        ! Calculate coefficients of pressure
        allocate(body%C_p(N_pressures), stat=stat)
        call check_allocation(stat, "surface pressures")
        do i=1,N_pressures
            body%C_p(i) = 1.-(norm(body%V(i,:))*freestream%U_inv)**2
        end do

        write(*,*) "Done."
        write(*,*) "        Maximum pressure coefficient:", maxval(body%C_p)
        write(*,*) "        Minimum pressure coefficient:", minval(body%C_p)

        ! Calculate total forces
        write(*,'(a)',advance='no') "     Calculating forces..."
        allocate(dC_F(N_pressures,3), stat=stat)
        call check_allocation(stat, "forces")
        do i=1,body%N_panels

            ! Discrete force coefficient acting on panel
            dC_F(i,:) = body%C_p(i)*body%panels(i)%A*body%panels(i)%normal

            ! Mirror
            if (body%mirrored .and. body%asym_flow) then
                n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)
                dC_F(i+body%N_panels,:) = body%C_p(i+body%N_panels)*body%panels(i)%A*n_mirrored
            end if
        end do

        ! Sum discrete forces
        C_F(:) = sum(dC_F, dim=1)/body%S_ref

        write(*,*) "Done."
        write(*,*) "        Cx:", C_F(1)
        write(*,*) "        Cy:", C_F(2)
        write(*,*) "        Cz:", C_F(3)

        ! Write report file
        if (report_file /= 'none') then

            open(1, file=report_file)

            ! Header
            write(1,'(a)') "TriPan Report (c) 2021 USU AeroLab"

            ! Solver results
            write(1,*) "Maximum residual inner total potential:"
            write(1,*) "Norm of residual innner total potential:"
            write(1,*) "Maximum pressure coefficient:", maxval(body%C_p)
            write(1,*) "Minimum pressure coefficient:", minval(body%C_p)
            write(1,*) "Cx:", C_F(1)
            write(1,*) "Cy:", C_F(2)
            write(1,*) "Cz:", C_F(3)

            close(1)

        end if
    
    end subroutine panel_solver_solve_source_free


end module panel_solver_mod