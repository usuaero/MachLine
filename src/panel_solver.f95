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


    subroutine panel_solver_init(this, settings, body_mesh)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body_mesh
        real :: control_point_offset

        ! Get settings
        call json_xtnsn_get(settings, 'formulation', this%formulation, 'morino')

        ! Initialize based on formulation
        if (this%formulation == 'morino') then
            call this%init_morino(settings, body_mesh)
        end if

    end subroutine panel_solver_init


    subroutine panel_solver_init_morino(this, settings, body_mesh)
        ! Initializes the solver to use the Morino formulation (i.e. zero inner perturbation potential)

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body_mesh

        real :: control_point_offset

        write(*,*)
        write(*,'(a)',advance='no') "     Placing control points..."

        ! Place control points inside the body
        call json_xtnsn_get(settings, 'control_point_offset', control_point_offset, 1e-5)
        call body_mesh%place_interior_control_points(control_point_offset)
        write(*,*) "Done."
    
    end subroutine panel_solver_init_morino


    subroutine panel_solver_solve(this, body_mesh, freestream_flow)
        ! Calls the relevant subroutine to solve the case based on the selected formulation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body_mesh
        type(flow),intent(in) :: freestream_flow

        ! Morino formulation
        if (this%formulation == 'morino') then
            call this%solve_morino(body_mesh, freestream_flow)
        end if

    end subroutine panel_solver_solve


    subroutine panel_solver_solve_morino(this, body, freestream_flow)
        ! Solves the Morino formulation for the given conditions

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(in) :: freestream_flow

        integer :: i, j, k
        real,dimension(:),allocatable :: source_inf, doublet_inf
        integer,dimension(:),allocatable :: source_verts, doublet_verts
        real,dimension(:,:),allocatable :: A, A_copy
        real,dimension(:),allocatable :: b
        integer :: stat, N_sigma, N_mu
        real,dimension(3) :: n_mirrored, cp_mirrored
        logical :: mirrored_and_asym

        ! Check mirroring and flow symmetry condition
        mirrored_and_asym = .false.
        if (body%mirrored) then
            if (.not. freestream_flow%symmetric_plane(body%mirror_plane)) then
                mirrored_and_asym = .true.
            end if
        end if

        write(*,*)
        write(*,'(a)',advance='no') "     Calculating source strengths..."

        ! Set source strengths
        if (source_order == 0) then

            ! Allocate source strength array
            if (mirrored_and_asym) then
                N_sigma = body%N_panels*2
            else
                N_sigma = body%N_panels
            end if

            allocate(body%sigma(N_sigma))

            ! Loop through panels
            do i=1,body%N_panels

                ! Existing panels
                body%sigma(i) = -inner(body%panels(i)%normal, freestream_flow%c0)

                ! Mirrored panels for asymmetric flow
                if (mirrored_and_asym) then

                    ! Get mirrored normal vector
                    n_mirrored = mirror_about_plane(body%panels(i)%normal, body%mirror_plane)

                    ! Calculate source strength
                    body%sigma(i+body%N_panels) = -inner(n_mirrored, freestream_flow%c0)

                end if
            end do

        end if
        write(*,*) "Done."

        ! Allocate space for inner potential calculations
        if (mirrored_and_asym) then
            allocate(body%phi_cp_sigma(body%N_cp*2), source=0., stat=stat)
        else
            allocate(body%phi_cp_sigma(body%N_cp), source=0., stat=stat)
        end if
        call check_allocation(stat, "induced potential vector")

        ! Allocate AIC matrix
        if (mirrored_and_asym) then
            allocate(A(body%N_cp*2, body%N_cp*2), source=0., stat=stat)
        else
            allocate(A(body%N_cp, body%N_cp), source=0., stat=stat)
        end if
        call check_allocation(stat, "AIC matrix")

        write(*,*)
        write(*,'(a)',advance='no') "     Calculating body influences..."

        ! Calculate source and doublet influences from body
        do i=1,body%N_cp
            do j=1,body%N_panels

                ! Get influence for existing->existing and mirrored->mirrored
                source_inf = body%panels(j)%get_source_potential(body%control_points(i,:), &
                                                                            source_verts)
                doublet_inf = body%panels(j)%get_doublet_potential(body%control_points(i,:), &
                                                                              doublet_verts)


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
                if (mirrored_and_asym .and. body%vertices(i)%mirrored_is_unique) then
                    if (source_order == 0) then
                            body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) &
                                                                       + source_inf(1)*body%sigma(j+body%N_panels)
                    end if
                    if (doublet_order == 1) then
                        do k=1,size(doublet_verts)
                            if (body%vertices(k)%mirrored_is_unique) then
                                A(i+body%N_cp,doublet_verts(k)+body%N_cp) = &
                                A(i+body%N_cp,doublet_verts(k)+body%N_cp) + doublet_inf(k)
                            end if
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

                    if (mirrored_and_asym) then

                        ! Source influences
                        if (source_order == 0) then

                            ! Influence of mirrored panel on existing control point
                            body%phi_cp_sigma(i) = body%phi_cp_sigma(i) &
                                                        + source_inf(1)*body%sigma(j+body%N_panels)

                            ! Influence of existing panel on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                body%phi_cp_sigma(i+body%N_cp) = body%phi_cp_sigma(i+body%N_cp) &
                                                                           + source_inf(1)*body%sigma(j)
                            end if
                        end if

                        ! Doublet influences
                        if (doublet_order == 1) then

                            ! Influence of mirrored panel on existing control point
                            do k=1,size(doublet_verts)
                                A(i,doublet_verts(k)+body%N_cp) = &
                                A(i,doublet_verts(k)+body%N_cp) + doublet_inf(k)
                            end do

                            ! Influence of existing panel on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                do k=1,size(doublet_verts)
                                    A(i+body%N_cp,doublet_verts(k)) = &
                                    A(i+body%N_cp,doublet_verts(k)) + doublet_inf(k)
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
        end do
        write(*,*) "Done."

        write(*,*)
        write(*,'(a)',advance='no') "     Calculating wake influences..."

        ! Loop through control points
        do i=1,body%N_cp

            ! Get doublet influence from wake
            do j=1,body%wake%N_panels

                doublet_inf = body%wake%panels(j)%get_doublet_potential(body%control_points(i,:),&
                                                                                   doublet_verts)

                ! Add to LHS
                if (doublet_order == 1) then
                    do k=1,size(doublet_verts)
                        A(i,doublet_verts(k)) = A(i,doublet_verts(k)) + doublet_inf(k)
                    end do
                end if
            end do

            ! For non-unique control points, set that the doublet strengths must be the same
            ! The RHS for these rows should still be zero
            if (mirrored_and_asym) then
                if (.not. body%vertices(i)%mirrored_is_unique) then
                    A(i+body%N_cp,i) = 1.
                    A(i+body%N_cp,i+body%N_cp) = -1.
                end if
            end if

        end do
        write(*,*) "Done."

        ! Make a copy of A (lu_solve replaces A with its decomposition)
        allocate(A_copy, source=A, stat=stat)
        call check_allocation(stat, "solver copy of AIC matrix")

        write(*,*)
        write(*,'(a)',advance='no') "     Solving linear system..."

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

        write(*,*)
        write(*,'(a)',advance='no') "     Calculating surface velocities and pressures..."

        ! Determine surface velocities
        allocate(body%V(body%N_panels,3), stat=stat)
        call check_allocation(stat, "surface velocity vectors")
        do i=1,body%N_panels
            body%V(i,:) = freestream_flow%U*(freestream_flow%c0 &
                               + body%panels(i)%get_velocity_jump(body%mu, body%sigma))
        end do

        ! Calculate coefficients of pressure
        allocate(body%C_p(body%N_panels), stat=stat)
        call check_allocation(stat, "surface pressures")
        do i=1,body%N_panels
            body%C_p(i) = 1.0-(norm(body%V(i,:))/freestream_flow%U)**2
        end do

        write(*,*) "Done."
        write(*,*) "        Maximum pressure coefficient:", maxval(body%C_p)
        write(*,*) "        Minimum pressure coefficient:", minval(body%C_p)
    
    end subroutine panel_solver_solve_morino


end module panel_solver_mod