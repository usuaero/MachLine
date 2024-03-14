program control_point_sensitivities

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
    integer :: start_count, end_count, i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, loc_up, loc_dn

    real,dimension(:,:),allocatable :: v, residuals3 , d_loc_FD

    integer :: i,j,k,m,n, N_verts, N_panels, vert, index
    real :: step
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(20) :: failure_log
    character(len=10) :: m_char



    test_failed = .true. ! assume test failed, if the test condition is met, test passed
    ! NOTE: on the sparse vector test, I assume the test passes, if it fails a test condition, test fails
    passed_tests = 0
    total_tests = 0

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
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))
    
    step = 0.00001
    index = 1
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONTROL POINT SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ CONTROL POINT SENSITIVITIES TEST ---------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CP_d_loc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST CP d_loc -----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of control point location WRT each design variable"
    write(*,*) ""

    !!!!!!!!! Finite Difference  d_n_g !!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  CENTRAL DIFFERENCE CP d_loc"
    
    ! perturb x1 up
    allocate(loc_up(N_verts*3))
    allocate(loc_dn(N_verts*3))
    allocate(d_loc_FD(3,N_verts*3))

    ! for each x, y, z of centr 1 
    do k=1,3
        ! do for each design variable
        do i=1,3
            do j=1,N_verts

                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                ! update panel geometry and calc
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do
                call test_mesh%calc_vertex_geometry()
                
                deallocate(test_solver%sigma_known)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)

     
                
                ! put the x y or z component of the vertex of interest (index) in a list
                loc_up(j + (i-1)*N_verts) = test_mesh%cp(index)%loc(k)

                ! perturb down the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                ! update panel geometry and calc
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do
                call test_mesh%calc_vertex_geometry()
                deallocate(test_solver%sigma_known)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)    
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)
                
                ! put the x y or z component of the vertex of interest (index) in a list
                loc_dn(j + (i-1)*N_verts) = test_mesh%cp(index)%loc(k)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_loc_FD(k,:) = (loc_up - loc_dn)/(2.*step)
            
    end do

    ! write results
    write(*,*) ""
    write(*,*) "                d_loc_FD vertex 1"
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z "
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x))') d_loc_FD(:,i)
    end do 
    
    !!!!!!!!!! ADJOINT d_n_g!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  ADJOINT CP d_loc"
    write(*,*) ""
    
    !write sparse matrix
    write(*,*) ""
    write(*,*) "         d_loc CP 1"
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z             sparse_index       full_index"

    do i=1,adjoint_mesh%cp(index)%d_loc%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%cp(index)%d_loc%columns(i)%vector_values, &
        i, adjoint_mesh%cp(index)%d_loc%columns(i)%full_index
    end do
    write(*,*) ""


    ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = adjoint_mesh%cp(index)%d_loc%get_values(i) - d_loc_FD(:,i)
    end do

    
    write(*,*) "         d_loc CP 1 expanded "
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z                                 residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%cp(index)%d_loc%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (any(abs(residuals3(:,i)) > 1.0e-8)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "CP d_loc test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "CP d_loc test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
        
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""


    !!!!!!!!!!!!!! CONTROL POINT SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------CONTROL POINT SENSITIVITIES TEST RESULTS--------------"
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

end program control_point_sensitivities