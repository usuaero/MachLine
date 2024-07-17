program dirichlet_test4

    ! tests various intermediate sensitivities 
    use adjoint_mod
    use base_geom_mod
    use panel_mod
    use flow_mod
    use surface_mesh_mod
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
    ! type(panel_solver) :: linear_solver, adjoint_linear_solver
    integer :: i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, n_g_up, n_g_dn, sum_up, sum_dn

    real,dimension(:,:),allocatable :: v, vertex_locs, residuals3,  d_n_g_FD, d_sum_FD

    ! real,dimension(:,:,:),allocatable ::  d_n_g_FD

    integer :: i,j,k,m,n,z, N_verts, N_panels, vert, index, cp_ind
    real :: step, error_allowed
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(100) :: failure_log
    character(len=10) :: m_char
    integer :: start_count, end_count
    real(16) :: count_rate, time



    test_failed = .false. 
    passed_tests = 0
    total_tests = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                             FROM MAIN

    !!!!!!!!!!!!!!! TEST INPUT (calc_adjoint = false) !!!!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()

    test_input = "dev\input_files\adjoint_inputs\dirichlet_test.json"
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

    !!!!!!!!!!!!!!!!!!!!! END TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call system_clock(start_count, count_rate)
    

    !!!!!!!!!!!!!!!!!!!!!!ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()

    adjoint_input = "dev\input_files\adjoint_inputs\dirichlet_adjoint_test.json"
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

    ! ! Initialize flow
    ! call json_xtnsn_get(adjoint_geom_settings, 'spanwise_axis', adjoint_spanwise_axis, '+y')
    ! call adjoint_freestream_flow%init(adjoint_flow_settings, adjoint_spanwise_axis)
    
    ! ! Get result files
    ! call json_xtnsn_get(adjoint_output_settings, 'body_file', adjoint_body_file, 'none')
    ! call json_xtnsn_get(adjoint_output_settings, 'wake_file', adjoint_wake_file, 'none')
    ! call json_xtnsn_get(adjoint_output_settings, 'control_point_file', adjoint_control_point_file, 'none')
    ! call json_xtnsn_get(adjoint_output_settings, 'mirrored_body_file', adjoint_mirrored_body_file, 'none')
    ! call json_xtnsn_get(adjoint_output_settings, 'offbody_points.points_file', adjoint_points_file, 'none')
    ! call json_xtnsn_get(adjoint_output_settings, 'offbody_points.output_file', adjoint_points_output_file, 'none')

    ! ! Get formulation type                                                  !
    ! call json_xtnsn_get(adjoint_solver_settings, 'formulation', adjoint_formulation, 'none')!

    ! ! Perform flow-dependent initialization on the surface mesh
    ! call adjoint_mesh%init_with_flow(adjoint_freestream_flow, adjoint_body_file, adjoint_wake_file, adjoint_formulation)
    !!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!

    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))
    allocate(n_g_up(N_verts*3))
    allocate(n_g_dn(N_verts*3))
    allocate(d_n_g_FD(3,N_verts*3))
    
    error_allowed = 1.0e-7
    step = 0.000001
    ! index = 1
    cp_ind = 1
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERTEX NORMAL SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) "                    VERTEX NORMAL SENSITIVITIES TEST                    "
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""
    
    do z = 1,N_verts ! for each control point
        cp_ind = z

        write(*,'(A,I5)') "VERTEX NORMAL TEST ", z

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
        
                    
                    ! get desired info
                    n_g_up(j + (i-1)*N_verts) = test_mesh%vertices(cp_ind)%n_g(k)

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    ! update panel geometry and calc
                    do m =1,N_panels
                        deallocate(test_mesh%panels(m)%n_hat_g)
                        call test_mesh%panels(m)%calc_derived_geom()
                    end do
                    
                    call test_mesh%calc_vertex_geometry()
                    
                    ! get desired info
                    n_g_dn(j + (i-1)*N_verts) = test_mesh%vertices(cp_ind)%n_g(k)
                    
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                end do 
            end do 
            
            ! central difference 
            d_n_g_FD(k,:) = (n_g_up - n_g_dn)/(2.*step)
                
        end do


        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%vertices(cp_ind)%d_n_g%get_values(i) - d_n_g_FD(:,i)
        end do

        if (maxval(abs(residuals3(:,:)))>error_allowed) then
            write(*,*) ""
            write(*,*) "     FLAGGED VALUES :"
            do i = 1, N_verts*3
                if (any(abs(residuals3(:,i))>error_allowed)) then
                    write(*,*) ""
                    write(*,*) "                                      d_n_g     "
                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_n_g_FD(:,i)
                    write(*, '(A25,8x,3(f25.10, 4x))') "               adjoint",   &
                    adjoint_mesh%vertices(cp_ind)%d_n_g%get_values(i)
                    write(*, '(A25,8x,3(f25.10, 4x))') "    residuals", residuals3(:,i)
                end if
            end do
        end if

        
        
        ! check if test failed
        do i=1,N_verts*3
            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                do j = 1,3
                    if (abs(d_n_g_FD(j,i))>1000.0) then
                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (1000.0>abs(d_n_g_FD(j,i)) .and. abs(d_n_g_FD(j,i))>100.0) then
                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (100.0>abs(d_n_g_FD(j,i)) .and. abs(d_n_g_FD(j,i))>10.0) then
                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (10.0>abs(d_n_g_FD(j,i)) .and. abs(d_n_g_FD(j,i))>1.0) then
                        if (abs(residuals3(j,i)) > error_allowed*10.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    else
                        if (abs(residuals3(j,i)) > error_allowed) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    end if
                end do
            end if
        end do
        if (test_failed) then
            total_tests = total_tests + 1
            write(*,'(A,I5,A,I5,A)')"                                               &
                               d_n_g vertex ",z," test FAILED"
            failure_log(total_tests-passed_tests) = "d_n_g test FAILED"
        else
            ! write(*,*) "        CALC d_n_g test PASSED"
            ! write(*,*) "" 
            ! write(*,*) ""
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if

        ! reset test failed for the next z loop
        test_failed = .false.


    end do ! z control points


    !!!!!!!!!!!!!! Vertex normal  SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "     VERTEX NORMAL SENSITIVITIES TEST RESULTS "
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) ""
    write(*,'((A), ES10.1)') "allowed residual = ", error_allowed
    write(*,*) ""

    write(*,'(I35,a14)') total_tests - passed_tests, " tests FAILED"
    write(*,*) ""
    write(*,'(I15,a9,I15,a14)') passed_tests, " out of ", total_tests, " tests PASSED"
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
    call system_clock(end_count)
    time = real(end_count - start_count)/(count_rate*60.0)
    write(*,'(A,f16.12, A)') " Total test time = ", time, " minutes"
    write(*,*) ""
    write(*,*) "----------------------"
    write(*,*) "Program Complete"
    write(*,*) "----------------------"



end program dirichlet_test4