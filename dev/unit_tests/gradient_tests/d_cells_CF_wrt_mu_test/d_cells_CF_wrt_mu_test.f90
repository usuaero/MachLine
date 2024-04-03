program d_cells_CF_wrt_mu_test

    ! tests various intermediate sensitivities 
    use adjoint_mod
    use base_geom_mod
    use panel_mod
    use flow_mod
    use surface_mesh_mod
    use panel_solver_mod
    use json_mod
    use json_xtnsn_mod
    ! use panel_solver_mod
    use helpers_mod
    
    implicit none

    !!!!!!!!!!!!!!!!!!! STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!
    character(100) :: test_input, adjoint_input
    character(len=:),allocatable :: body_file, wake_file, control_point_file, points_file, &
    mirrored_body_file, points_output_file,&
    adjoint_body_file, adjoint_wake_file, adjoint_control_point_file, adjoint_points_file, &
    adjoint_mirrored_body_file, adjoint_points_output_file
    character(len=:),allocatable :: report_file, spanwise_axis, adjoint_spanwise_axis
    character(len=:),allocatable :: formulation, adjoint_formulation

    type(json_file) :: input_json, adjoint_input_json
    type(json_value),pointer :: flow_settings, &
                                geom_settings, &
                                solver_settings, &
                                processing_settings,&
                                output_settings, &
                                adjoint_flow_settings, &
                                adjoint_geom_settings, &
                                adjoint_solver_settings,&
                                adjoint_processing_settings,&
                                adjoint_output_settings
    type(surface_mesh) :: test_mesh, adjoint_mesh
    type(flow) :: freestream_flow, adjoint_freestream_flow
    type(panel_solver) :: test_solver, adjoint_solver
    type(eval_point_geom) :: test_geom, adjoint_geom
    type(dod) :: test_dod_info, adjoint_dod_info
    type(integrals) :: test_int, adjoint_int
    integer :: start_count, end_count, i_unit
    logical :: exists, found 
    integer :: adjoint_solver_stat, test_solver_stat
    type(sparse_vector) :: zeros

    real,dimension(3) :: adjoint_P, test_P, test_v_d, test_v_s
    type(sparse_matrix) :: adjoint_d_P_term2
    type(sparse_matrix) :: adjoint_d_P
    type(sparse_matrix) :: adjoint_d_v_d_panel

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals
    real,dimension(:,:),allocatable ::  residuals3, cells_CF_wrt_mu_up, cells_CF_wrt_mu_dn, d_cells_CF_wrt_mu_FD

    integer :: i,j,k,m,n,p, N_verts, N_panels, vert, index, cp_ind, row,col, stat
    real :: step, error_allowed, cp_offset
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute
    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(50) :: failure_log
    character(len=10) :: m_char

    !!!!!!!!!!!!!!!!!!! END TESTING STUFF !!!!!!!!!!!!!!!!!!!!!11

    test_failed = .true. ! assume test failed, if the test condition is met, test passed
    ! NOTE: on the sparse vector test, I assume the test passes, if it fails a test condition, test fails
    passed_tests = 0
    total_tests = 0

    error_allowed = 1.0e-6
    
    step = 0.00001
    index = 1
    cp_ind = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                             FROM MAIN

    !!!!!!!!!!!!!!! TEST INPUT (calc_adjoint = false) !!!!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()

    test_input = "dev\input_files\adjoint_inputs\test.json"
    test_input = trim(test_input)

    ! Check it exists
    inquire(file=test_input, exist=exists)
    if (.not. exists) then
        write(*,*) "!!! The file ", test_input, " does not exist. Quitting..."
        stop
    end if

    ! Load settings from input file
    call input_json%load_file(filename=test_input)
    call json_check()
    call input_json%get('flow', flow_settings, found)
    call input_json%get('geometry', geom_settings, found)
    call input_json%get('solver', solver_settings, found)
    call input_json%get('post_processing', processing_settings, found)
    call input_json%get('output', output_settings, found)

    ! Initialize surface mesh
    call test_mesh%init(geom_settings)


    ! Initialize flow
    call json_xtnsn_get(geom_settings, 'spanwise_axis', spanwise_axis, '+y')
    call freestream_flow%init(flow_settings, spanwise_axis)
    
    ! Get result files
    call json_xtnsn_get(output_settings, 'body_file', body_file, 'none')
    call json_xtnsn_get(output_settings, 'wake_file', wake_file, 'none')
    call json_xtnsn_get(output_settings, 'control_point_file', control_point_file, 'none')
    call json_xtnsn_get(output_settings, 'mirrored_body_file', mirrored_body_file, 'none')
    call json_xtnsn_get(output_settings, 'offbody_points.points_file', points_file, 'none')
    call json_xtnsn_get(output_settings, 'offbody_points.output_file', points_output_file, 'none')

    !!!!!!!!!!!!!!!!!!!!!! WAKE_DEV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get formulation type                                                  !
    call json_xtnsn_get(solver_settings, 'formulation', formulation, 'none')!
    !!!!!!!!!!!!!!!!!!!!!!! END_WAKE_DEV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Perform flow-dependent initialization on the surface mesh
    call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)

    ! Initialize panel solver
    call test_solver%init(solver_settings, processing_settings, test_mesh, freestream_flow, control_point_file)
    
    ! pull out the cp offset
    call json_xtnsn_get(solver_settings, 'control_point_offset', cp_offset, 1.e-7)
    
    ! Set default status
    test_solver_stat = 0

    ! Allocate known influence storage
    allocate(test_solver%I_known(test_mesh%N_cp), source=0., stat=stat)
    call check_allocation(stat, "known influence vector")

    ! Allocate AIC matrix
    allocate(test_solver%A(test_mesh%N_cp, test_solver%N_unknown), source=0., stat=stat)
    ! call check_allocation(stat, "AIC matrix")

    ! Allocate b vector
    allocate(test_solver%b(test_mesh%N_cp), source=0., stat=stat)
    call check_allocation(stat, "b vector")

    ! Calculate source strengths
    call test_solver%calc_source_strengths(test_mesh)

    ! Calculate body influences
    call test_solver%calc_body_influences(test_mesh)

    call test_solver%assemble_BC_vector(test_mesh)

    ! Solve the linear system
    call test_solver%solve_system(test_mesh, test_solver_stat)
    
    ! Check for errors
    if (test_solver_stat /= 0) return

    ! Calculate velocities
    call test_solver%calc_cell_velocities(test_mesh)

    ! Calculate velocities
    call test_solver%calc_pressures(test_mesh)

    call test_solver%calc_forces(test_mesh)
    
    
    !!!!!!!!!!!!!!!!!!!!! END TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()
    
    adjoint_input = "dev\input_files\adjoint_inputs\adjoint_test.json"
    adjoint_input = trim(adjoint_input)

    ! Check it exists
    inquire(file=adjoint_input, exist=exists)
    if (.not. exists) then
        write(*,*) "!!! The file ", adjoint_input, " does not exist. Quitting..."
        stop
    end if

    ! Load settings from input file
    call adjoint_input_json%load_file(filename=adjoint_input)
    call json_check()
    call adjoint_input_json%get('flow', adjoint_flow_settings, found)
    call adjoint_input_json%get('geometry', adjoint_geom_settings, found)
    call adjoint_input_json%get('solver', adjoint_solver_settings, found)
    call adjoint_input_json%get('post_processing', adjoint_processing_settings, found)
    call adjoint_input_json%get('output', adjoint_output_settings, found)
    
    ! Initialize surface mesh
    call adjoint_mesh%init(adjoint_geom_settings)
    !call adjoint_mesh%init_adjoint()

    ! Initialize flow
    call json_xtnsn_get(adjoint_geom_settings, 'spanwise_axis', adjoint_spanwise_axis, '+y')
    call adjoint_freestream_flow%init(adjoint_flow_settings, adjoint_spanwise_axis)
    
    ! Get result files
    call json_xtnsn_get(adjoint_output_settings, 'body_file', adjoint_body_file, 'none')
    call json_xtnsn_get(adjoint_output_settings, 'wake_file', adjoint_wake_file, 'none')
    call json_xtnsn_get(adjoint_output_settings, 'control_point_file', adjoint_control_point_file, 'none')
    call json_xtnsn_get(adjoint_output_settings, 'mirrored_body_file', adjoint_mirrored_body_file, 'none')
    call json_xtnsn_get(adjoint_output_settings, 'offbody_points.points_file', adjoint_points_file, 'none')
    call json_xtnsn_get(adjoint_output_settings, 'offbody_points.output_file', adjoint_points_output_file, 'none')
    
    !!!!!!!!!!!!!!!!!!!!!! WAKE_DEV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get formulation type                                                  !
    call json_xtnsn_get(adjoint_solver_settings, 'formulation', adjoint_formulation, 'none')!
    !!!!!!!!!!!!!!!!!!!!!!! END_WAKE_DEV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Perform flow-dependent initialization on the surface mesh
    call adjoint_mesh%init_with_flow(adjoint_freestream_flow, adjoint_body_file, adjoint_wake_file, adjoint_formulation)
    
    ! Initialize panel solver
    call adjoint_solver%init(adjoint_solver_settings, adjoint_processing_settings, adjoint_mesh, &
    adjoint_freestream_flow, adjoint_control_point_file)
    write(*,*) "check for size mismatch"
    ! solve
    call adjoint_solver%solve(adjoint_mesh, adjoint_solver_stat, adjoint_formulation,adjoint_freestream_flow)
    

    
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! d_cells_CF_wrt_mu_test !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ d_cells_CF_wrt_mu_test ---&
    ------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""

    ! allocate data holders
    allocate(cells_CF_wrt_mu_up(3,N_verts*3))
    allocate(cells_CF_wrt_mu_dn(3,N_verts*3))
    allocate(d_cells_CF_wrt_mu_FD(3,N_verts*3))

        
    !!!!!!!!! CENTRAL DIFFERENCE  d_cells_CF_wrt_mu_test column ", k, !!!!!!!!!
    write(*,*) ""
    write(*,*) "--------------------------------------------------------------------------------------"
    write(*,*) "                   d_cells_CF_wrt_mu_test "
    write(*,*) "--------------------------------------------------------------------------------------"
    write(*,*) ""

    do k=1,N_panels
        do i=1,test_mesh%N_verts  !!!!!!!!!! THIS IS DIFFERENT because its wrt mu

            ! perturb up the current mu
            test_mesh%mu(i) = test_mesh%mu(i) + step
            ! write(*,*) " perturb up"
            
            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
            ! deallocate stuff
            deallocate(test_mesh%V_cells_inner)
            deallocate(test_mesh%V_cells)
        
            ! Calculate velocities
            call test_solver%calc_cell_velocities(test_mesh)
            
            deallocate(test_mesh%C_P_inc)

            ! Calculate velocities
            call test_solver%calc_pressures(test_mesh)

            deallocate(test_mesh%dC_f)

            call test_solver%calc_forces(test_mesh)
                        
            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
            
            ! get the needed info
            cells_CF_wrt_mu_up(:, i) = test_mesh%dC_F(:,k)
            
            
            ! perturb down the current design variable
            ! write(*,*) " perturb down"
            test_mesh%mu(i) = test_mesh%mu(i) - 2.*step
    
            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
            ! deallocate stuff
            deallocate(test_mesh%V_cells_inner)
            deallocate(test_mesh%V_cells)
        
            ! Calculate velocities
            call test_solver%calc_cell_velocities(test_mesh)
            
            deallocate(test_mesh%C_P_inc)

            ! Calculate velocities
            call test_solver%calc_pressures(test_mesh)

            deallocate(test_mesh%dC_f)

            call test_solver%calc_forces(test_mesh)
                        
                            
            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
            
            ! get the needed info
            cells_CF_wrt_mu_dn(:, i) = test_mesh%dC_F(:,k)

            ! restore geometry
            test_mesh%mu(i) = test_mesh%mu(i) + step
            
        end do

            
            ! central difference 
        d_cells_CF_wrt_mu_FD = (cells_CF_wrt_mu_up - cells_CF_wrt_mu_dn)/(2.*step)
                
        write(*,*) ""
        
        
        ! write results
        write(*,*) ""
        write(*,'(A,I1)') "              Central Diff  d_cells_CF_wrt_mu_test cell ",k
        write(*,*) "  d_cells_CF_wrt_mu_ x         d_cells_CF_wrt_mu_ y         d_cells_CF_wrt_mu_ z "
        do i = 1, N_verts
            write(*, '(3(f16.10, 10x))') d_cells_CF_wrt_mu_FD(:,i)
        end do 
        write(*,*) ""
        
        
        do i=1,N_verts*3
            residuals3(:,i) = adjoint_mesh%d_cell_forces_wrt_mu(k)%get_values(i) - d_cells_CF_wrt_mu_FD(:,i)
        end do
        
        write(*,*) ""
        write(*,'(A,I1)') "          Adjoint      d_cells_CF_wrt_mu cell ",k
        write(*,*) "  d_cells_CF_wrt_mu_x      d_cells_CF_wrt_mu_y       d_cells_CF_wrt_mu_z         &
        residuals"
        do i = 1, N_verts
            write(*, '(3(f16.10, 8x),3x, 3(f16.10, 4x))') &
            adjoint_mesh%d_cell_forces_wrt_mu(k)%get_values(i), residuals3(:,i)
        end do
        write(*,*) ""
        
        ! check if test failed
        do i=1,N_verts
            if (any(abs(residuals3(:,i)) > error_allowed)) then
                test_failed = .true.
                exit
            else 
                test_failed = .false.
            end if
        end do
        if (test_failed) then
            total_tests = total_tests + 1
            
            failure_log(total_tests-passed_tests) = "d_cells_CF_wrt_mu TEST FAILED"
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1,A)') " d_cells_CF_wrt_mu cell ",k," TEST PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""

    end do
    

    
    
    
    
    !!!!!!!!!!!!!! d_cells_CF_wrt_mu_test RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------d_cells_CF_wrt_mu_test RESULTS--------------"
    write(*,*) ""
    write(*,'((A), ES10.1)') "allowed residual = ", error_allowed
    write(*,'((A), ES10.1)') "control point offset = ", cp_offset
    write(*,*) ""

    write(*,'(I15,a14)') total_tests - passed_tests, " tests FAILED"
    write(*,*) ""
    write(*,'(I4,a9,I2,a14)') passed_tests, " out of ", total_tests, " tests PASSED"
    if (passed_tests < total_tests)then
        write(*,*) ""
        write(*,*) "----------------------"
        write(*,*) "Failure Log:"
        write(*,*) ""
        do i=1,total_tests-passed_tests
            write(*,*) failure_log(i)
        end do
    end if
    
    write(*,*) ""
    write(*,*) ""
    write(*,*) "----------------------"
    write(*,*) "Program Complete"
    write(*,*) "----------------------"


end program d_cells_CF_wrt_mu_test