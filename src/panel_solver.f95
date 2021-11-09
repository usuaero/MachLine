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


    subroutine panel_solver_solve_morino(this, body_mesh, freestream_flow)
        ! Solves the Morino formulation for the given conditions

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body_mesh
        type(flow),intent(in) :: freestream_flow

        integer :: i, j, k
        real,dimension(:),allocatable :: influence
        integer,dimension(:),allocatable :: vertex_indices
        real,dimension(:,:),allocatable :: A, A_copy
        real,dimension(:),allocatable :: b
        integer :: stat

        ! Set source strengths
        write(*,*)
        write(*,'(a)',advance='no') "     Calculating source strengths..."
        if (source_order == 0) then

            ! Allocate source strength array
            allocate(body_mesh%sigma(body_mesh%N_panels))

            ! Loop through panels
            do i=1,body_mesh%N_panels
                body_mesh%sigma(i) = -inner(body_mesh%panels(i)%normal, freestream_flow%c0)
            end do

        end if
        write(*,*) "Done."

        ! Calculate source influences
        write(*,*)
        write(*,'(a)',advance='no') "     Calculating source influences..."
        allocate(body_mesh%phi_cp_sigma(body_mesh%N_verts), source=0., stat=stat)
        call check_allocation(stat)
        do i=1,body_mesh%N_verts
            do j=1,body_mesh%N_panels

                ! Get source influence
                influence = body_mesh%panels(j)%get_source_potential(body_mesh%control_points(i,:), vertex_indices)

                ! Add to RHS
                if (source_order == 0) then
                    body_mesh%phi_cp_sigma(i) = body_mesh%phi_cp_sigma(i) + influence(1)*body_mesh%sigma(j)
                end if

            end do
        end do
        write(*,*) "Done."

        ! Calculate doublet influences
        write(*,*)
        write(*,'(a)',advance='no') "     Calculating doublet influences..."

        ! Allocate AIC matrix
        allocate(A(body_mesh%N_verts, body_mesh%N_verts), source=0., stat=stat)
        call check_allocation(stat)

        ! Loop through control points
        do i=1,body_mesh%N_verts

            ! Get doublet influence from body
            do j=1,body_mesh%N_panels

                influence = body_mesh%panels(j)%get_doublet_potential(body_mesh%control_points(i,:), vertex_indices)

                ! Add to LHS
                if (doublet_order == 1) then
                    do k=1,size(vertex_indices)
                        A(i,vertex_indices(k)) = A(i,vertex_indices(k)) + influence(k)
                    end do
                end if
            end do

            ! Get doublet influence from wake
            do j=1,body_mesh%wake%N_panels

                influence = body_mesh%wake%panels(j)%get_doublet_potential(body_mesh%control_points(i,:), vertex_indices)

                ! Add to LHS
                if (doublet_order == 1) then
                    do k=1,size(vertex_indices)
                        A(i,vertex_indices(k)) = A(i,vertex_indices(k)) + influence(k)
                    end do
                end if
            end do

        end do
        write(*,*) "Done."

        ! Allocate solution memory
        allocate(body_mesh%mu(body_mesh%N_verts), stat=stat)
        call check_allocation(stat)
        allocate(A_copy, source=A, stat=stat)
        call check_allocation(stat)

        ! Solve
        write(*,*)
        write(*,'(a)',advance='no') "     Solving linear system..."
        b = -body_mesh%phi_cp_sigma
        call lu_solve(body_mesh%N_verts, A_copy, b, body_mesh%mu)
        write(*,*) "Done."
        deallocate(A_copy)
        deallocate(b)

        ! Calculate potential at control points
        body_mesh%phi_cp_mu = matmul(A, body_mesh%mu)
        body_mesh%phi_cp = body_mesh%phi_cp_mu+body_mesh%phi_cp_sigma
        write(*,*) "        Maximum residual inner potential:", maxval(abs(body_mesh%phi_cp))
        write(*,*) "        Norm of residual innner potential:", sqrt(sum(body_mesh%phi_cp**2))

        write(*,*)
        write(*,'(a)',advance='no') "     Calculating surface velocities and pressures..."

        ! Determine surface velocities
        allocate(body_mesh%V(body_mesh%N_panels,3), stat=stat)
        call check_allocation(stat)
        do i=1,body_mesh%N_panels
            body_mesh%V(i,:) = freestream_flow%U*(freestream_flow%c0 &
                               + body_mesh%panels(i)%get_velocity_jump(body_mesh%mu, body_mesh%sigma))
        end do

        ! Calculate coefficients of pressure
        allocate(body_mesh%C_p(body_mesh%N_panels), stat=stat)
        call check_allocation(stat)
        do i=1,body_mesh%N_panels
            body_mesh%C_p(i) = 1.0-(norm(body_mesh%V(i,:))/freestream_flow%U)**2
        end do

        write(*,*) "Done."
        write(*,*) "        Maximum pressure coefficient:", maxval(body_mesh%C_p)
        write(*,*) "        Minimum pressure coefficient:", minval(body_mesh%C_p)
    
    end subroutine panel_solver_solve_morino


end module panel_solver_mod