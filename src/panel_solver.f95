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
            procedure :: solve => panel_solver_solve
            procedure :: solve_morino => panel_solver_solve_morino

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
        end if

    end subroutine panel_solver_init


    subroutine panel_solver_init_morino(this, settings, body)
        ! Initializes the solver to use the Morino formulation (i.e. zero inner perturbation potential)

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body

        real :: offset

        write(*,*)
        write(*,'(a)',advance='no') "     Placing control points..."

        ! Place control points inside the body
        call json_xtnsn_get(settings, 'control_point_offset', offset, 1e-5)
        call body%place_interior_control_points(offset)
        write(*,*) "Done."
    
    end subroutine panel_solver_init_morino


    subroutine panel_solver_solve(this, body, freestream)
        ! Calls the relevant subroutine to solve the case based on the selected formulation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(in) :: freestream

        ! Morino formulation
        if (this%formulation == 'morino') then
            call this%solve_morino(body, freestream)
        end if

    end subroutine panel_solver_solve


    subroutine panel_solver_solve_morino(this, body, freestream)
        ! Solves the Morino formulation for the given conditions

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(in) :: freestream

        integer :: i, j, k
        real,dimension(:),allocatable :: source_inf, doublet_inf
        integer,dimension(:),allocatable :: source_verts, doublet_verts
        real,dimension(:,:),allocatable :: A, A_copy
        real,dimension(:),allocatable :: b
        integer :: stat, N_sigma, N_mu
        real,dimension(3) :: n_mirrored, cp_mirrored

        write(*,'(a)',advance='no') "     Calculating source strengths..."

        ! Set source strengths
        if (source_order == 0) then

            ! Determine necessary number of source strengths
            if (body%mirrored_and_asym) then
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
                if (body%mirrored_and_asym) then

                    ! Get mirrored normal vector
                    n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)

                    ! Calculate source strength
                    body%sigma(i+body%N_panels) = -inner(n_mirrored, freestream%c0)

                end if
            end do

        end if
        write(*,*) "Done."

        ! Determine number of doublet strengths (some will be repeats for mirrored vertices)
        if (body%mirrored_and_asym) then
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

                ! Influence of mirrored panels on mirrored control points, if the mirrored control point is meant to be unique
                if (body%mirrored_and_asym .and. body%vertices(i)%mirrored_is_unique) then

                    if (source_order == 0) then
                        body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) &
                                                         + source_inf(1)*body%sigma(j+body%N_panels)
                    end if

                    if (doublet_order == 1) then
                        do k=1,size(doublet_verts)
                            A(i+body%N_cp,doublet_verts(k)+body%N_cp) = A(i+body%N_cp,doublet_verts(k)+body%N_cp) + doublet_inf(k)
                        end do
                    end if
                end if

                ! Get influence for existing->mirrored and mirrored->existing
                if (body%mirrored) then

                    ! Get mirrored point (recall the Green's function is reciprocal)
                    cp_mirrored = mirror_about_plane(body%control_points(i,:), body%mirror_plane)

                    ! Calculate influences
                    source_inf = body%panels(j)%get_source_potential(cp_mirrored, source_verts)
                    doublet_inf = body%panels(j)%get_doublet_potential(cp_mirrored, doublet_verts)

                    if (body%mirrored_and_asym) then

                        ! Source influences
                        if (source_order == 0) then

                            ! Influence of mirrored panel on existing control point
                            body%phi_cp_sigma(i) = body%phi_cp_sigma(i) &
                                                        + source_inf(1)*body%sigma(j+body%N_panels)

                            ! Influence of existing panel on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) + source_inf(1)*body%sigma(j)
                            end if
                        end if

                        ! Doublet influences
                        if (doublet_order == 1) then

                            ! Influence of mirrored panel on existing control point
                            do k=1,size(doublet_verts)
                                A(i,doublet_verts(k)+body%N_cp) = A(i,doublet_verts(k)+body%N_cp) + doublet_inf(k)
                            end do

                            ! Influence of existing panel on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
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
            if (body%mirrored_and_asym .and. .not. body%vertices(i)%mirrored_is_unique) then
                A(i+body%N_cp,i) = 1.
                A(i+body%N_cp,i+body%N_cp) = -1.
            end if

        end do
        write(*,*) "Done."

        ! Calculate influence of wake
        write(*,'(a)',advance='no') "     Calculating wake influences..."

        ! Loop through control points
        do i=1,body%N_cp

            ! Get doublet influence from wake
            ! Note that for the wake, in the case of mirrored mesh with asymmetric flow, the mirrored wake panels have actually been created.
            ! In this case, there are technically no mirrored panels, and this loop will cycle through both existing and mirrored panels.
            ! For symmetric flow, mirrored panels still need to be added as before.
            do j=1,body%wake%N_panels

                ! Caclulate influence on existing control points
                doublet_inf = body%wake%panels(j)%get_doublet_potential(body%control_points(i,:), doublet_verts)

                ! Influence of existing wake panels on existing control points
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

                    if (body%mirrored_and_asym) then

                        ! Influence of mirrored panel on mirrored control point
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

        write(*,'(a)',advance='no') "     Solving linear system..."

        ! Make a copy of A (lu_solve replaces A with its decomposition)
        allocate(A_copy, source=A, stat=stat)
        call check_allocation(stat, "solver copy of AIC matrix")

        ! Solve
        b = -body%phi_cp_sigma
        call lu_solve(size(b), A_copy, b, body%mu)
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
        allocate(body%V(body%N_panels,3), stat=stat)
        call check_allocation(stat, "surface velocity vectors")
        do i=1,body%N_panels
            body%V(i,:) = freestream%U*(freestream%c0 + body%panels(i)%get_velocity_jump(body%mu, body%sigma))
        end do

        ! Calculate coefficients of pressure
        allocate(body%C_p(body%N_panels), stat=stat)
        call check_allocation(stat, "surface pressures")
        do i=1,body%N_panels
            body%C_p(i) = 1.0-(norm(body%V(i,:))/freestream%U)**2
        end do

        write(*,*) "Done."
        write(*,*) "        Maximum pressure coefficient:", maxval(body%C_p)
        write(*,*) "        Minimum pressure coefficient:", minval(body%C_p)
    
    end subroutine panel_solver_solve_morino


end module panel_solver_mod