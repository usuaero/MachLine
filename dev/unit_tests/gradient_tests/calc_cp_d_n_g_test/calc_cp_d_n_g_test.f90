program calc_cp_d_n_g_test
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
    type(sparse_matrix),dimension(3) :: d_v_d
    integer :: start_count, end_count, i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals
    real,dimension(:,:),allocatable ::  residuals3, cp_n_g_up, cp_n_g_dn, d_cp_n_g_FD

    type(sparse_vector),dimension(3,3) :: d_v_d_M_adjoint

    integer :: i,j,k,m,n,p, N_verts, N_panels, vert, index, cp_ind
    real,dimension(:,:), allocatable :: v_s, v_d
    real :: step, error_allowed, cp_offset
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(20) :: failure_log
    character(len=10) :: m_char

    !!!!!!!!!!!!!!!!!!! END TESTING STUFF !!!!!!!!!!!!!!!!!!!!!11

    test_failed = .true. ! assume test failed, if the test condition is met, test passed
    ! NOTE: on the sparse vector test, I assume the test passes, if it fails a test condition, test fails
    passed_tests = 0
    total_tests = 0

    error_allowed = 1.0e-8 
    step = 0.00001
    index = 1
    cp_ind = 2
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
    
    ! calc CALC BASIC GEOM geom of relation between cp1 and panel1 
    test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(cp_ind)%loc,freestream_flow,.false.)
    test_dod_info = test_mesh%panels(index)%check_dod(test_mesh%cp(cp_ind)%loc, freestream_flow, .false.)
    test_int = test_mesh%panels(index)%calc_integrals(test_geom, 'velocity', freestream_flow,.false., test_dod_info)
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

    ! !calc CALC BASIC F integral sensitivities cp1 and panel1 
    ! adjoint_geom = adjoint_mesh%panels(index)%calc_subsonic_geom_adjoint(adjoint_mesh%cp(cp_ind),&
    !                                                             adjoint_freestream_flow)
    ! adjoint_dod_info = adjoint_mesh%panels(index)%check_dod(adjoint_mesh%cp(cp_ind)%loc, &
    !     adjoint_freestream_flow, .false.)
    
    ! adjoint_int = adjoint_mesh%panels(index)%calc_integrals(adjoint_geom, 'velocity',&
    !     adjoint_freestream_flow,.false., adjoint_dod_info)
    
    ! call adjoint_mesh%panels(index)%calc_integrals_adjoint(adjoint_geom,adjoint_int,adjoint_freestream_flow&
    !     , .false., adjoint_dod_info)
    
    ! d_v_d_M_adjoint = adjoint_mesh%panels(index)%assemble_v_d_M_space_adjoint(adjoint_int, adjoint_geom, &
    !     adjoint_freestream_flow, .false.)
    
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALC cp d_n_g TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ CALC cp d_n_g TEST ---&
    ------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC cp d_n_g  (cp 1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

    write(*,'(A, I1,A)') "---------------------------------- TEST CALC cp ",cp_ind," d_n_g&
     ---------------------------------"
    write(*,*) ""
    write(*,*) ""


    !!!!!!!!! CENTRAL DIFFERENCE cp 1 d_n_g !!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  CENTRAL DIFFERENCE"
    
    ! perturb x1 up
    allocate(cp_n_g_up(3,N_verts*3))
    allocate(cp_n_g_dn(3,N_verts*3))
    allocate(d_cp_n_g_FD(3,N_verts*3))

    do i=1,3
        do j=1,N_verts

            ! perturb up the current design variable
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            ! write(*,*) " perturb up"
            
            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
           
            ! ! update panel geometry and calcs
            ! deallocate(test_mesh%panels(index)%vertices)
            ! deallocate(test_mesh%panels(index)%n_hat_g)
            ! deallocate(test_mesh%panels(index)%edge_is_discontinuous)
            ! call test_mesh%panels(index)%init(test_mesh%vertices(1),test_mesh%vertices(2),test_mesh%vertices(3),index)
            
            ! update vertex normal
            do m =1,N_panels
                deallocate(test_mesh%panels(m)%n_hat_g)
                call test_mesh%panels(m)%calc_derived_geom()
            end do
            
            call test_mesh%calc_vertex_geometry()
            
            ! ! update with flow
            ! deallocate(test_mesh%panels(index)%vertices_ls)
            ! deallocate(test_mesh%panels(index)%n_hat_ls)
            ! deallocate(test_mesh%panels(index)%b)
            ! deallocate(test_mesh%panels(index)%b_mir)  
            ! deallocate(test_mesh%panels(index)%sqrt_b)
            ! deallocate(test_mesh%panels(index)%i_vert_d)
            ! deallocate(test_mesh%panels(index)%S_mu_inv)
            ! deallocate(test_mesh%panels(index)%T_mu)
            ! deallocate(test_mesh%panels(index)%i_panel_s)
            ! call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
            ! call test_mesh%panels(index)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
            ! test_mesh%vertices,.false.)
            
            ! recalculates cp locations
            deallocate(test_solver%sigma_known)
            deallocate(test_mesh%cp)
            deallocate(test_solver%P)
            call test_solver%init(solver_settings, processing_settings, &
            test_mesh, freestream_flow, control_point_file)
            
            ! update v_d and doublet inf
            ! call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc, freestream_flow,.false.,v_s, v_d)
            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
            
            ! get the needed info
            cp_n_g_up(:,j + (i-1)*N_verts) = test_mesh%cp(cp_ind)%n_g
            
            
            ! perturb down the current design variable
            ! write(*,*) " perturb down"
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
            
            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
            ! ! update panel geometry and calcs
            ! deallocate(test_mesh%panels(index)%vertices)
            ! deallocate(test_mesh%panels(index)%n_hat_g)
            ! deallocate(test_mesh%panels(index)%edge_is_discontinuous)
            ! call test_mesh%panels(index)%init(test_mesh%vertices(1),test_mesh%vertices(2),test_mesh%vertices(3),index)
           
            do m =1,N_panels
                deallocate(test_mesh%panels(m)%n_hat_g)
                call test_mesh%panels(m)%calc_derived_geom()
            end do
            
            ! update vertex normal
            call test_mesh%calc_vertex_geometry()
            
            ! ! update with flow
            ! deallocate(test_mesh%panels(index)%vertices_ls)
            ! deallocate(test_mesh%panels(index)%n_hat_ls)
            ! deallocate(test_mesh%panels(index)%b)
            ! deallocate(test_mesh%panels(index)%b_mir)  
            ! deallocate(test_mesh%panels(index)%sqrt_b)
            ! deallocate(test_mesh%panels(index)%i_vert_d)
            ! deallocate(test_mesh%panels(index)%S_mu_inv)
            ! deallocate(test_mesh%panels(index)%T_mu)
            ! deallocate(test_mesh%panels(index)%i_panel_s)
            ! call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
            ! call test_mesh%panels(index)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
            !     test_mesh%vertices,.false.)
            
            ! recalculates cp locations
            deallocate(test_solver%sigma_known)
            deallocate(test_mesh%cp)
            deallocate(test_solver%P)
            call test_solver%init(solver_settings, processing_settings, &
            test_mesh, freestream_flow, control_point_file)
           
            ! update v_d and doublet inf
            ! call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc, freestream_flow,.false.,v_s, v_d)
            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
            
            ! get the needed info
            cp_n_g_dn(:,j + (i-1)*N_verts) = test_mesh%cp(cp_ind)%n_g

            
            ! restore geometry
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
        end do 
    end do

    ! central difference 
    d_cp_n_g_FD = (cp_n_g_up - cp_n_g_dn)/(2.*step)

    ! write results
    write(*,*) ""
    write(*,*) "                d_n_g_FD cp 1"
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z "
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x))') d_cp_n_g_FD(:,i)
    end do 
    
    !!!!!!!!!! ADJOINT d_n_g!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  ADJOINT d_n_g cp 1"
    write(*,*) ""
    
    !write sparse matrix
    write(*,*) ""
    write(*,*) "         d_n_g cp 1"
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z             sparse_index       full_index"

    do i=1,adjoint_mesh%cp(cp_ind)%d_n_g%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%cp(cp_ind)%d_n_g%columns(i)%vector_values, &
        i, adjoint_mesh%cp(cp_ind)%d_n_g%columns(i)%full_index
    end do
    write(*,*) ""


    ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = adjoint_mesh%cp(cp_ind)%d_n_g%get_values(i) - d_cp_n_g_FD(:,i)
    end do

    
    write(*,*) "         d_n_g cp 1 expanded "
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z                                 residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%cp(cp_ind)%d_n_g%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,N_verts*3
        if (any(abs(residuals3(:,i)) > error_allowed)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        
        failure_log(total_tests-passed_tests) = "CALC cp d_n_g test FAILED"

        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "CALC cp d_n_g test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
        
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""
    





    !!!!!!!!!!!!!! CALC cp d_n_g TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------CALC cp d_n_g TEST RESULTS--------------"
    write(*,*) ""
    write(*,'((A), ES10.2)') "allowed residual = ", error_allowed
    write(*,'((A), ES10.2)') "control point offset = ", cp_offset
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

end program calc_cp_d_n_g_test