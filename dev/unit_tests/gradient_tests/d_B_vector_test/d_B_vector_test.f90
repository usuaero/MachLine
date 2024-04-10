program d_b_vector_test

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
    type(sparse_vector) :: zeros

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals, b_up, b_dn, d_b_FD
    real,dimension(:,:),allocatable ::  residuals3 

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

    error_allowed = 1.0e-7
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
    
    ! ! Allocate known influence storage
    allocate(test_solver%I_known(test_mesh%N_cp), source=0., stat=stat)
    call check_allocation(stat, "known influence vector")

    ! Allocate b vector
    allocate(test_solver%b(test_mesh%N_cp), source=0., stat=stat)
    call check_allocation(stat, "b vector")

    ! assemble BC vector
    call test_solver%assemble_BC_vector(test_mesh)

    ! Set b vector
    test_solver%b = test_solver%BC - test_solver%I_known
    
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

    ! allocate  d_b_vector
    call zeros%init(adjoint_mesh%adjoint_size)
    allocate(adjoint_solver%d_b_vector(adjoint_mesh%N_cp), source=zeros, stat=stat)
    call check_allocation(stat, "Adjoint b sensitivity vector")

    ! if adjoint, assemble b sensitivity vector
    call adjoint_solver%assemble_adjoint_b_vector(adjoint_mesh)

    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! d_b_vector TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ d_b_vector TEST ---&
    ------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""

    ! allocate data holders
    allocate(b_up(N_verts*3))
    allocate(b_dn(N_verts*3))
    allocate(d_b_FD(N_verts*3))

        
    !!!!!!!!! CENTRAL DIFFERENCE (panel 1, cp 1) d_v_d_M row 1, column ", k, !!!!!!!!!
    write(*,*) ""
    write(*,*) "--------------------------------------------------------------------------------------"
    write(*,*) "                   d_b_vector test "
    write(*,*) "--------------------------------------------------------------------------------------"
    write(*,*) ""


    do k = 1,test_mesh%N_cp
        
        do i=1,3
            do j=1,N_verts

                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                ! write(*,*) " perturb up"
                
                !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
            
                ! update vertex normal
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do
                
                call test_mesh%calc_vertex_geometry()
                
                ! update with flow
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%vertices_ls)
                    deallocate(test_mesh%panels(m)%n_hat_ls)
                    deallocate(test_mesh%panels(m)%b)
                    deallocate(test_mesh%panels(m)%b_mir)  
                    deallocate(test_mesh%panels(m)%sqrt_b)
                    deallocate(test_mesh%panels(m)%i_vert_d)
                    deallocate(test_mesh%panels(m)%S_mu_inv)
                    deallocate(test_mesh%panels(m)%T_mu)
                    ! deallocate(test_mesh%panels(m)%i_panel_s)
                    call test_mesh%panels(m)%init_with_flow(freestream_flow, .false., 0)
                    call test_mesh%panels(m)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                    test_mesh%vertices,.false.)
                end do
                
                ! recalculates cp locations
                deallocate(test_solver%sigma_known)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)
                
                ! deallocate stuff
                deallocate(test_solver%I_known)
                deallocate(test_solver%BC)
                deallocate(test_solver%b)

                ! Allocate known influence storage
                allocate(test_solver%I_known(test_mesh%N_cp), source=0., stat=stat)
                call check_allocation(stat, "known influence vector")
                ! ! Allocate known influence storage

                ! Allocate b vector
                allocate(test_solver%b(test_mesh%N_cp), source=0., stat=stat)
                call check_allocation(stat, "b vector")

                ! assemble BC vector
                call test_solver%assemble_BC_vector(test_mesh)

                ! Set b vector
                test_solver%b = test_solver%BC - test_solver%I_known
                                
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                b_up(j + (i-1)*N_verts) = test_solver%b(k)
                
                
                ! perturb down the current design variable
                ! write(*,*) " perturb down"
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
                
                !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
            
                ! update vertex normal
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do
                
                call test_mesh%calc_vertex_geometry()
                
                ! update with flow
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%vertices_ls)
                    deallocate(test_mesh%panels(m)%n_hat_ls)
                    deallocate(test_mesh%panels(m)%b)
                    deallocate(test_mesh%panels(m)%b_mir)  
                    deallocate(test_mesh%panels(m)%sqrt_b)
                    deallocate(test_mesh%panels(m)%i_vert_d)
                    deallocate(test_mesh%panels(m)%S_mu_inv)
                    deallocate(test_mesh%panels(m)%T_mu)
                    ! deallocate(test_mesh%panels(m)%i_panel_s)
                    call test_mesh%panels(m)%init_with_flow(freestream_flow, .false., 0)
                    call test_mesh%panels(m)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                    test_mesh%vertices,.false.)
                end do
                
                ! recalculates cp locations
                deallocate(test_solver%sigma_known)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)
                
                ! deallocate stuff
                deallocate(test_solver%I_known)
                deallocate(test_solver%BC)
                deallocate(test_solver%b)

                ! Allocate known influence storage
                allocate(test_solver%I_known(test_mesh%N_cp), source=0., stat=stat)
                call check_allocation(stat, "known influence vector")
                ! ! Allocate known influence storage

                ! Allocate b vector
                allocate(test_solver%b(test_mesh%N_cp), source=0., stat=stat)
                call check_allocation(stat, "b vector")

                ! assemble BC vector
                call test_solver%assemble_BC_vector(test_mesh)

                ! Set b vector
                test_solver%b = test_solver%BC - test_solver%I_known
                                
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                b_dn(j + (i-1)*N_verts) = test_solver%b(k)

                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_b_FD = (b_up - b_dn)/(2.*step)
                
        write(*,*) ""
        
        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_solver%d_b_vector(k)%get_value(i) - d_b_FD(i)
        end do
        
        ! write results
        write(*,*) ""
        write(*,'(A, I1,A)') " d_b_vector_test &
        d_b (",k,") "
        write(*,*) ""
        write(*,*) "  d_b        d_b_FD         residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x,(f14.10, 4x),3x, (f14.10, 4x))') &
            adjoint_solver%d_b_vector(k)%get_value(i), d_b_FD(i), residuals(i)
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
                failure_log(total_tests-passed_tests) = "TEST FAILURE in d_b_vector(1)"
            elseif (k ==2) then
                failure_log(total_tests-passed_tests) = "TEST FAILURE in d_b_vector(2)"
            elseif (k ==3) then
                failure_log(total_tests-passed_tests) = "TEST FAILURE in d_b_vector(3)"
            elseif (k ==4) then
                failure_log(total_tests-passed_tests) = "TEST FAILURE in d_b_vector(4)"
            elseif (k ==5) then
                failure_log(total_tests-passed_tests) = "TEST FAILURE in d_b_vector(5)"
            else
                failure_log(total_tests-passed_tests) = "TEST FAILURE in d_b_vector(6)"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1,A,I1,A)') " d_b_vector(",k,") test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""

    ! end k loop
    end do




    !!!!!!!!!!!!!! d_b_vector TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------d_b_vector TEST RESULTS--------------"
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


end program d_b_vector_test