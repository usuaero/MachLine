program calc_inf_adjoint_test2
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
    type(sparse_vector),dimension(3,3) :: d_v_d
    integer :: start_count, end_count, i_unit
    logical :: exists, found
    real,dimension(:),allocatable :: doublet_inf
    real,dimension(:,:),allocatable :: v_s, v_d
    type(sparse_vector),dimension(3) :: inf_adjoint, inf_adjoint2

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals, inf_up, inf_dn, d_inf_FD
    real,dimension(:,:),allocatable ::  residuals3 

    integer :: i,j,k,m,n,p, N_verts, N_panels, vert, index, cp_ind
    real,dimension(3,3) :: v_d_M_space
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
    
    call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc, freestream_flow,.false.,v_s, v_d)
    doublet_inf = matmul(test_mesh%cp(cp_ind)%n_g, matmul(freestream_flow%B_mat_g, v_d))
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
    ! , .false., adjoint_dod_info)


    inf_adjoint2 = adjoint_mesh%panels(index)%calc_doublet_inf_adjoint2(adjoint_mesh%cp(cp_ind), &
    adjoint_freestream_flow)

    ! call adjoint_mesh%panels(index)%calc_velocity_influences_adjoint(adjoint_mesh%cp(cp_ind), &
    ! adjoint_freestream_flow, d_v_d)
    
    ! inf_adjoint = adjoint_mesh%panels(index)%calc_doublet_inf_adjoint(adjoint_mesh%cp(cp_ind),&
    ! adjoint_freestream_flow, d_v_d)
    
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALC inf_adjoint SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ CALC inf_adjoint SENSITIVITIES TEST ---&
    ------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""

    ! allocate data holders
    allocate(inf_up(N_verts*3))
    allocate(inf_dn(N_verts*3))
    allocate(d_inf_FD(N_verts*3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST2 CALC inf_adjoint (panel x, cp x) row 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    ! panel vertex 
    do k=1,3

        write(*,'(A, I1,A,I1,A,I1,A)') "---------------------------------- TEST2 CALC inf_adjoint &
        (sensitivity of the influence of panel ",index,", panel vertex ", k,", on cp ",cp_ind,") ------&
        ---------------------------"
        write(*,*) ""
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE (panel 1, cp 1) d_v_d_M row 1, column ", k, !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
        do i=1,3
            do j=1,N_verts

                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                ! write(*,*) " perturb up"
                
                !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
               
                ! update panel geometry and calcs
                deallocate(test_mesh%panels(index)%vertices)
                deallocate(test_mesh%panels(index)%n_hat_g)
                deallocate(test_mesh%panels(index)%edge_is_discontinuous)
                call test_mesh%panels(index)%init(test_mesh%vertices(1),test_mesh%vertices(2),test_mesh%vertices(3),index)
                
                ! update vertex normal
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do
                
                call test_mesh%calc_vertex_geometry()
                
                ! update with flow
                deallocate(test_mesh%panels(index)%vertices_ls)
                deallocate(test_mesh%panels(index)%n_hat_ls)
                deallocate(test_mesh%panels(index)%b)
                deallocate(test_mesh%panels(index)%b_mir)  
                deallocate(test_mesh%panels(index)%sqrt_b)
                deallocate(test_mesh%panels(index)%i_vert_d)
                deallocate(test_mesh%panels(index)%S_mu_inv)
                deallocate(test_mesh%panels(index)%T_mu)
                deallocate(test_mesh%panels(index)%i_panel_s)
                call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                call test_mesh%panels(index)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                test_mesh%vertices,.false.)
                
                ! recalculates cp locations
                deallocate(test_solver%sigma_known)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)
                
                ! update v_d and doublet inf
                call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc, freestream_flow,.false.,v_s, v_d)
                doublet_inf = matmul(test_mesh%cp(cp_ind)%n_g, matmul(freestream_flow%B_mat_g, v_d))
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                inf_up(j + (i-1)*N_verts) = doublet_inf(k)
                
                
                ! perturb down the current design variable
                ! write(*,*) " perturb down"
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
                
                !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
                ! update panel geometry and calcs
                deallocate(test_mesh%panels(index)%vertices)
                deallocate(test_mesh%panels(index)%n_hat_g)
                deallocate(test_mesh%panels(index)%edge_is_discontinuous)
                call test_mesh%panels(index)%init(test_mesh%vertices(1),test_mesh%vertices(2),test_mesh%vertices(3),index)
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do
                
                ! update vertex normal
                call test_mesh%calc_vertex_geometry()
                
                ! update with flow
                deallocate(test_mesh%panels(index)%vertices_ls)
                deallocate(test_mesh%panels(index)%n_hat_ls)
                deallocate(test_mesh%panels(index)%b)
                deallocate(test_mesh%panels(index)%b_mir)  
                deallocate(test_mesh%panels(index)%sqrt_b)
                deallocate(test_mesh%panels(index)%i_vert_d)
                deallocate(test_mesh%panels(index)%S_mu_inv)
                deallocate(test_mesh%panels(index)%T_mu)
                deallocate(test_mesh%panels(index)%i_panel_s)
                call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                call test_mesh%panels(index)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                    test_mesh%vertices,.false.)
                
                ! recalculates cp locations
                deallocate(test_solver%sigma_known)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)
               
                ! update v_d and doublet inf
                call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc, freestream_flow,.false.,v_s, v_d)
                doublet_inf = matmul(test_mesh%cp(cp_ind)%n_g, matmul(freestream_flow%B_mat_g, v_d))
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                inf_dn(j + (i-1)*N_verts) = doublet_inf(k)

                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_inf_FD = (inf_up - inf_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1,A,I1,A,I1,A)') "  CENTRAL DIFFERENCE CALC inf_adjoint2 &
        (panel ",index,", panel vertex ",k,", cp,",cp_ind,")"
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  inf_adjoint2"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_inf_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC inf_adjoint (panel x, cp x) !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------"
        write(*,'(A, I1,A,I1,A,I1,A)') "  ADJOINT  inf_adjoint2, (panel ",index,", panel vertex ",k, &
        ", cp,",cp_ind,")"
        write(*,*) "---------------------------------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1,A,I1,A,I1,A)') "  adjoint CALC inf_adjoint2 &
        (panel ",index,", panel vertex ",k, ", cp,",cp_ind,")  sparse"
        write(*,*) ""
        write(*,*) "  inf_adjoint2              sparse_index       full_index"
        

        do i=1,inf_adjoint2(k)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') inf_adjoint2(k)%elements(i)%value, &
            i, inf_adjoint2(k)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = inf_adjoint2(k)%get_value(i) - d_inf_FD(i)
        end do

        write(*,'(A, I1,A,I1,A,I1,A)') "  adjoint CALC inf_adjoint2 &
        (panel ",index,", panel vertex ",k, ", cp,",cp_ind,"),  expanded"
        write(*,*) ""
        write(*,*) "  inf_adjoint2                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') inf_adjoint2(k)%get_value(i), residuals(i)
        end do
        write(*,*) ""
        write(*,*) ""

        ! check if test failed
        do i=1,N_verts*3
            if (abs(residuals(i)) > error_allowed) then
                test_failed = .true.
                exit
            else 
                test_failed = .false.
            end if
        end do
        if (test_failed) then
            total_tests = total_tests + 1
            if (k ==1) then
                failure_log(total_tests-passed_tests) = "CALC inf_adjoint2 panel vertex &
                1 test FAILED"
            elseif (k ==2) then
                failure_log(total_tests-passed_tests) = "CALC inf_adjoint2 panel vertex &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC inf_adjoint2 panel vertex &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1,A,I1,A,I1,A)') "CALC inf_adjoint2 (panel ",index,", panel vertex ",&
            k, ", cp,",cp_ind,"), test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! k loop end
    end do




    !!!!!!!!!!!!!! CALC inf_adjoint2 RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------CALC inf_adjoint2 TEST RESULTS--------------"
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

end program calc_inf_adjoint_test2