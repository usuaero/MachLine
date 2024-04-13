program calc_basic_geom_cp
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
    integer :: i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, P_g_up, P_g_dn, P_ls_up, P_ls_dn, h_up, h_dn, d_h_FD,&
    h2_up, h2_dn, d_h2_FD, v_xi_up, v_xi_dn, d_v_xi_FD, v_eta_up, v_eta_dn, d_v_eta_FD, d_ls_up, d_ls_dn, d_d_ls_FD

    real,dimension(:,:),allocatable :: v, residuals3 , d_P_g_FD, d_P_ls_FD

    integer :: i,j,k,m,n,p,y,z, N_verts, N_panels, vert, index, cp_ind
    real :: step, error_allowed
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(500) :: failure_log
    character(len=10) :: m_char
    real,dimension(:,:),allocatable :: maxRs
    integer(8) :: start_count, end_count
    real(16) :: count_rate, time




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

    ! calc CALC BASIC GEOM geom of relation between cp1 and panel1 
    test_geom = test_mesh%panels(1)%calc_basic_geom(test_mesh%cp(1)%loc,.false.)
    !!!!!!!!!!!!!!!!!!!!! END TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(start_count, count_rate)

    

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

    ! ! calc CALC BASIC GEOM geom sensitivity of relation between cp1 and panel1 
    ! adjoint_geom = adjoint_mesh%panels(index)%calc_basic_geom_adjoint(adjoint_mesh%cp(cp_ind)%loc,&
    !     adjoint_mesh%cp(cp_ind)%d_loc, .false.)
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    allocate(P_g_up(N_verts*3))
    allocate(P_g_dn(N_verts*3))
    allocate(d_P_g_FD(3,N_verts*3))
    allocate(P_ls_up(N_verts*3))
    allocate(P_ls_dn(N_verts*3))
    allocate(d_P_ls_FD(2,N_verts*3))
    allocate(h_up(N_verts*3))
    allocate(h_dn(N_verts*3))
    allocate(d_h_FD(N_verts*3))
    allocate(h2_up(N_verts*3))
    allocate(h2_dn(N_verts*3))
    allocate(d_h2_FD(N_verts*3))
    allocate(v_xi_up(N_verts*3))
    allocate(v_xi_dn(N_verts*3))
    allocate(d_v_xi_FD(N_verts*3))
    allocate(v_eta_up(N_verts*3))
    allocate(v_eta_dn(N_verts*3))
    allocate(d_v_eta_FD(N_verts*3))
    allocate(d_ls_up(N_verts*3))
    allocate(d_ls_dn(N_verts*3))
    allocate(d_d_ls_FD(N_verts*3))

    error_allowed = 1.0e-7
    step = 0.000001
    index = 1
    cp_ind = 1
    allocate(maxRs(500,3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALC BASIC GEOM (CONTROL POINTS) SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "-----------------------------------------------------------------------------------"
    write(*,*) "                CALC BASIC GEOM (CONTROL POINTS) SENSITIVITIES TEST "
    write(*,*) "-----------------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""


    do y=1,N_panels ! loop through panels
        index = y

        do z= 1,N_verts ! loop through control points
            cp_ind  = z

            ! update adjoint calcs
            adjoint_geom = adjoint_mesh%panels(index)%calc_basic_geom_adjoint&
            (adjoint_mesh%cp(cp_ind)%loc, adjoint_mesh%cp(cp_ind)%d_loc, .false.)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_h and d_h2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do i=1,3
                do j=1,N_verts

                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                    !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!

                    ! update panel geometry and calc
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
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    
                    ! recalculates cp locations
                    deallocate(test_solver%sigma_known)
                    deallocate(test_mesh%cp)
                    deallocate(test_solver%P)
                    call test_solver%init(solver_settings, processing_settings, &
                    test_mesh, freestream_flow, control_point_file)

                    ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                    test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)
            
                    !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                    
                    ! get desired info
                    h_up(j + (i-1)*N_verts) = test_geom%h
                    h2_up(j + (i-1)*N_verts) = test_geom%h2

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
                    ! update panel geometry and calc
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
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)

                    ! recalculates cp locations
                    deallocate(test_solver%sigma_known)
                    deallocate(test_mesh%cp)
                    deallocate(test_solver%P) 
                    call test_solver%init(solver_settings, processing_settings, &
                    test_mesh, freestream_flow, control_point_file)

                    ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                    test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)

                    !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                    ! get desired info
                    h_dn(j + (i-1)*N_verts) = test_geom%h
                    h2_dn(j + (i-1)*N_verts) = test_geom%h2
                    
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                end do 
            end do 
            
            ! central difference 
            d_h_FD(:) = (h_up - h_dn)/(2.*step)
            d_h2_FD(:) = (h2_up - h2_dn)/(2.*step)
                    

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_h !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! calculate residuals
            do i =1, N_verts*3
                residuals(i) = adjoint_geom%d_h%get_value(i) - d_h_FD(i)
            end do

            ! write(*,*) ""
            ! write(*,*) "------------------------------------------------------------------"
            ! write(*,*) "------------------------------------------------------------------"
            ! write(*,'(A, I4, A, I4)') "   d_h panel ",index, " cp ", cp_ind
            ! write(*,*) ""

            ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals))
            if (maxval(abs(residuals))>error_allowed) then
                write(*,*) "          Failed values :"
                write(*,*) "          d_h FD             adjoint d_h             residual"
                do i = 1, N_verts*3
                    if (abs(residuals(i))>error_allowed) then
                        write(*, '(8x,(f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))') &
                        d_h_FD(i), adjoint_geom%d_h%get_value(i), residuals(i)
                    end if
                end do
            end if
            
            ! write(*,*) " d_h FD                adjoint d_h                 residual"

            ! do i = 1, N_verts*3
            !     write(*, '((f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))') &
            !     d_h_FD(i), adjoint_geom%d_h%get_value(i), residuals(i)
            ! end do
            ! write(*,*) ""
            ! write(*,*) ""

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
                failure_log(total_tests-passed_tests) = "CALC d_h test FAILED"
                write(*,*) failure_log(total_tests-passed_tests)
            else
                ! write(*,*) "        CALC d_h test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


            
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC   d_h2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = adjoint_geom%d_h2%get_value(i) - d_h2_FD(i)
            end do

            ! write(*,*) ""
            ! write(*,*) "------------------------------------------------------------------"
            ! write(*,*) "------------------------------------------------------------------"
            ! write(*,'(A, I4, A, I4)') "   d_h2 panel ",index, " cp ", cp_ind
            ! write(*,*) ""
            ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals))
            if (maxval(abs(residuals))>error_allowed) then
                write(*,*) "          Failed values :"
                write(*,*) "          d_h2 FD             adjoint d_h2             residual"
                do i = 1, N_verts*3
                    if (abs(residuals(i))>error_allowed) then
                        write(*, '(8x,(f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))') &
                        d_h2_FD(i), adjoint_geom%d_h2%get_value(i), residuals(i)
                    end if
                end do
            end if
            
            ! write(*,*) " d_h2 FD                adjoint d_h2                 residual"

            ! do i = 1, N_verts*3
            !     write(*, '((f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))') &
            !     d_h2_FD(i), adjoint_geom%d_h2%get_value(i), residuals(i)
            ! end do
            ! write(*,*) ""
            ! write(*,*) ""

          

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
                failure_log(total_tests-passed_tests) = "CALC d_h2 test FAILED"
                write(*,*) failure_log(total_tests-passed_tests)
            else
                ! write(*,*) "        CALC d_h2 test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_P_g !!!!!!!!!!!!!!!!!!!!!!!
            ! for each panel vertex
            do k = 1,3
                
                ! do for each design variable
                do i=1,3
                    do j=1,N_verts

                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                        !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!

                        ! update panel geometry and calc
                        do m =1,N_panels
                            deallocate(test_mesh%panels(m)%n_hat_g)
                            call test_mesh%panels(m)%calc_derived_geom()
                        end do

                        call test_mesh%calc_vertex_geometry()
                        deallocate(test_solver%sigma_known)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P)

                        ! recalculates cp locations
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)

                        ! calc CALC BASIC GEOM geom of relation between cp1 and panel1 
                        test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)
            
                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                        
                        ! put the x y or z component of the vertex of interest (index) in a list
                        P_g_up(j + (i-1)*N_verts) = test_geom%P_g(k)
                        
                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                        !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
                        ! update panel geometry and calc
                        do m =1,N_panels
                            deallocate(test_mesh%panels(m)%n_hat_g)
                            call test_mesh%panels(m)%calc_derived_geom()
                        end do

                        call test_mesh%calc_vertex_geometry()
                        deallocate(test_solver%sigma_known)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P) 

                        ! recalculates cp locations
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)

                        ! calc CALC BASIC GEOM geom of relation between cp1 and panel1 
                        test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)

                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                        ! put the x y or z component of the vertex of interest (index) in a list
                        P_g_dn(j + (i-1)*N_verts) = test_geom%P_g(k)
                       
                        
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    end do 
                end do 
                
                ! central difference 
                d_P_g_FD(k,:) = (P_g_up - P_g_dn)/(2.*step)
                
                    
            end do ! end k loop



            ! calculate residuals3
            do i =1, N_verts*3
                residuals3(:,i) = adjoint_geom%d_P_g%get_values(i) - d_P_g_FD(:,i)
            end do


            ! write(*,*) ""
            ! write(*,*) "------------------------------------------------------------------"
            ! write(*,*) "------------------------------------------------------------------"
            ! write(*,'(A,I4,A,I4,A,I4)') "   d_P_g panel ",index," vertex ",&
            ! k," cp ",cp_ind
            ! write(*,*) ""
            ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals3))
            if (abs(maxval(residuals3))>error_allowed) then
                write(*,*) "          Failed values :"
                do i = 1, N_verts*3
                    if (any(abs(residuals3(:,i))>error_allowed)) then
                        write(*,*) "         Central Difference    d_P_g      x, y, and z"
                        write(*, '(8x,3(f14.10, 4x))') d_P_g_FD(:,i)
                        write(*,*) "        adjoint       d_P_g       x, y, and z                  &
                           residuals"
                        write(*, '(8x,3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_geom%d_P_g%get_values(i), residuals3(:,i)
                    end if
                end do
            end if


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
                if (k==1) then
                    failure_log(total_tests-passed_tests) = "Calc d_P_g panel index 1 test FAILED"
                elseif (k==2) then
                    failure_log(total_tests-passed_tests) = "Calc d_P_g panel index 2 test FAILED"
                else 
                    failure_log(total_tests-passed_tests) = "Calc d_P_g panel index 3 test FAILED"
                end if
                write(*,*) failure_log(total_tests-passed_tests)
            else
                ! write(*,*) "        CALC d_P_g test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_v_xi, d_v_eta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

            ! for each panel vertex
            do k = 1,3
                do i=1,3
                    do j=1,N_verts
        
                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
        
                        !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
        
                        ! update panel geometry and calc
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
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
        
                        ! recalculates cp locations
                        deallocate(test_solver%sigma_known)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P)
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)
        
                        ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                        test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(1)%loc,.false.)
        
                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
        
                        ! put the x y or z component of the vertex of interest (index) in a list
                        v_xi_up(j + (i-1)*N_verts) = test_geom%v_xi(k)
                        v_eta_up(j + (i-1)*N_verts) = test_geom%v_eta(k)
        
                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
        
                        !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
                        ! update panel geometry and calc
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
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
        
                        ! recalculates cp locations
                        deallocate(test_solver%sigma_known)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P) 
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)
        
                        ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                        test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(1)%loc,.false.)
        
                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
        
                        ! put the x y or z component of the vertex of interest (index) in a list
                        v_xi_dn(j + (i-1)*N_verts) = test_geom%v_xi(k)
                        v_eta_dn(j + (i-1)*N_verts) = test_geom%v_eta(k)
        
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    end do 
                end do 
                
                ! central difference 
                d_v_xi_FD(:) = (v_xi_up - v_xi_dn)/(2.*step)
                d_v_eta_FD(:) = (v_eta_up - v_eta_dn)/(2.*step)
                    
                
                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_v_xi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                ! calculate residuals
                do i =1, N_verts*3
                    residuals(i) = adjoint_geom%d_v_xi(k)%get_value(i) - d_v_xi_FD(i)
                end do
                
                
                ! write(*,*) ""
                ! write(*,*) "------------------------------------------------------------------"
                ! write(*,*) "------------------------------------------------------------------"
                ! write(*,'(A,I4,A,I4,A,I4)') "   d_v_xi panel ",index,&
                ! " vertex ",k," cp ",cp_ind
                ! write(*,*) ""
                ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals))
                if (maxval(abs(residuals))>error_allowed) then
                    write(*,*) "          Failed values :"
                    write(*,*) "          d_v_xi_FD            adjoint d_v_xi              residuals"
                    do i = 1, N_verts*3
                        if (abs(residuals(i))>error_allowed) then
                            write(*, '(8x,(f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))')&
                            d_v_xi_FD(i), adjoint_geom%d_v_xi(k)%get_value(i), residuals(i)
                        end if
                    end do
                end if

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
                    if (k==1) then
                        failure_log(total_tests-passed_tests) = "calc d_v_xi coordinate 1 test FAILED"
                        elseif (k==2) then
                            failure_log(total_tests-passed_tests) = "calc d_v_xi coordinate 2 test FAILED"
                        else 
                            failure_log(total_tests-passed_tests) = "calc d_v_xi coordinate 3 test FAILED"
                        end if
                        write(*,*) failure_log(total_tests-passed_tests)
                    else
                        ! write(*,*) "        CALC d_v_xi test PASSED"
                        ! write(*,*) "" 
                        ! write(*,*) ""
                        passed_tests = passed_tests + 1
                        total_tests = total_tests + 1
                        
                    end if
                test_failed = .false.
    

                

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_v_eta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                ! calculate residuals
                do i =1, N_verts*3
                    residuals(i) = adjoint_geom%d_v_eta(k)%get_value(i) - d_v_eta_FD(i)
                end do
            
                
                ! write(*,*) ""
                ! write(*,*) "------------------------------------------------------------------"
                ! write(*,*) "------------------------------------------------------------------"
                ! write(*,'(A,I4,A,I4,A,I4)') "   d_v_eta panel ",index,&
                ! " vertex ",k," cp ",cp_ind
                ! write(*,*) ""
                ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals))
                if (maxval(abs(residuals))>error_allowed) then
                    write(*,*) "          Failed values :"
                    write(*,*) "          d_v_eta_FD            adjoint d_v_eta              residuals"
                    do i = 1, N_verts*3
                        if (abs(residuals(i))>error_allowed) then
                            write(*, '(8x,(f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))')&
                            d_v_eta_FD(i), adjoint_geom%d_v_eta(k)%get_value(i), residuals(i)
                        end if
                    end do
                end if


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
                    if (k==1) then
                        failure_log(total_tests-passed_tests) = "calc d_v_eta coordinate 1 test FAILED"
                        elseif (k==2) then
                            failure_log(total_tests-passed_tests) = "calc d_v_eta coordinate 2 test FAILED"
                        else 
                            failure_log(total_tests-passed_tests) = "calc d_v_eta coordinate 3 test FAILED"
                        end if
                        write(*,*) failure_log(total_tests-passed_tests)
                    else
                        ! write(*,*) "        CALC d_v_eta test PASSED"
                        ! write(*,*) "" 
                        ! write(*,*) ""
                        passed_tests = passed_tests + 1
                        total_tests = total_tests + 1
                    
                end if
                test_failed = .false.
        
                
            end do ! end k loop
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_P_ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            

            ! for each xi and eta 
            do k=1,2
               
                do i=1,3
                    do j=1,N_verts

                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                        !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!

                        ! update panel geometry and calc
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
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                        
                        ! recalculates cp locations
                        deallocate(test_solver%sigma_known)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P)
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)

                        ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                        test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)
            
                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                        
                        ! put the x y or z component of the vertex of interest (index) in a list
                        P_ls_up(j + (i-1)*N_verts) = test_geom%P_ls(k)

                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                        !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
                        ! update panel geometry and calc
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
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)

                        ! recalculates cp locations
                        deallocate(test_solver%sigma_known)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P) 
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)

                        ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                        test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)

                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                        ! put the x y or z component of the vertex of interest (index) in a list
                        P_ls_dn(j + (i-1)*N_verts) = test_geom%P_ls(k)
                        
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    end do 
                end do 
                
                ! central difference 
                d_P_ls_FD(k,:) = (P_ls_up - P_ls_dn)/(2.*step)
                        
                ! calculate residuals3
                do i =1, N_verts*3
                    residuals(i) = adjoint_geom%d_P_ls(k)%get_value(i) - d_P_ls_FD(k,i)
                end do
                
                ! write(*,*) ""
                ! ! write results
                ! write(*,*) ""
                ! if (k==1) then
                !     write(*,*) "------------------------------------------------------------------"
                !     write(*,*) "------------------------------------------------------------------"
                !     write(*,'(A,I4,A,I4)') "   d_P_ls xi coord panel ",index,&
                !      " cp ", cp_ind
                ! else
                !     write(*,*) "------------------------------------------------------------------"
                !     write(*,*) "------------------------------------------------------------------"
                !     write(*,'(A,I4,A,I4)') "   d_P_ls eta coord panel ",index,&
                !      " cp ", cp_ind
                ! end if 
                ! write(*,*) ""
                ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals))
                if (maxval(abs(residuals))>error_allowed) then
                    write(*,*) "          Failed values :"
                    if (k==1) then
                        write(*,*) "        d_P_ls xi  FD         adjoint d_P_ls xi            residual"
                    else
                        write(*,*) "        d_P_ls eta  FD        adjoint d_P_ls eta             residual"
                    end if 
                    do i = 1, N_verts*3
                        if (abs(residuals(i))>error_allowed) then
                            write(*, '(8x,(f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))') &
                            d_P_ls_FD(k,i), adjoint_geom%d_P_ls(k)%get_value(i), residuals(i)
                        end if
                    end do
                end if

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
                    if (k==1) then
                        failure_log(total_tests-passed_tests) = "calc d_P_ls xi coordinate test FAILED"
                    else
                        failure_log(total_tests-passed_tests) = "calc d_P_ls eta coordinate test FAILED"
                    end if
                    write(*,*) failure_log(total_tests-passed_tests)
                else
                    ! if (k==1) then
                    !     write(*,*) "        calc d_P_ls xi coordinate test PASSED"
                    ! else
                    !     write(*,*) "        calc d_P_ls xi coordinate test PASSED"
                    ! end if
                    ! write(*,*) "" 
                    ! write(*,*) ""
                    passed_tests = passed_tests + 1
                    total_tests = total_tests + 1
                    
                end if
                test_failed = .false.
                
            
            ! k loop
            end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM  d_d_ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! for each panel side
            do p=1,3

                ! for each xi and eta 
                do k=1,2

                    do i=1,3
                        do j=1,N_verts

                            ! perturb up the current design variable
                            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!

                            ! update panel geometry and calc
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
                            call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                            
                            ! recalculates cp locations
                            deallocate(test_solver%sigma_known)
                            deallocate(test_mesh%cp)
                            deallocate(test_solver%P)
                            call test_solver%init(solver_settings, processing_settings, &
                            test_mesh, freestream_flow, control_point_file)

                            ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                            test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)
                
                            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                            
                            ! put the x y or z component of the vertex of interest (index) in a list
                            d_ls_up(j + (i-1)*N_verts) = test_geom%d_ls(k,p)

                            ! perturb down the current design variable
                            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
                            ! update panel geometry and calc
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
                            call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)

                            ! recalculates cp locations
                            deallocate(test_solver%sigma_known)
                            deallocate(test_mesh%cp)
                            deallocate(test_solver%P) 
                            call test_solver%init(solver_settings, processing_settings, &
                            test_mesh, freestream_flow, control_point_file)

                            ! update CALC BASIC GEOM geom of relation between cp1 and panel1 
                            test_geom = test_mesh%panels(index)%calc_basic_geom(test_mesh%cp(cp_ind)%loc,.false.)

                            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                            ! put the x y or z component of the vertex of interest (index) in a list
                            d_ls_dn(j + (i-1)*N_verts) = test_geom%d_ls(k,p)
                            
                            ! restore geometry
                            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        end do 
                    end do 
                    
                    ! central difference 
                    d_d_ls_FD = (d_ls_up - d_ls_dn)/(2.*step)


                    ! calculate residuals3
                    do i =1, N_verts*3
                        residuals(i) = adjoint_geom%d_d_ls(k,p)%get_value(i) - d_d_ls_FD(i)
                    end do
                       
                    ! if (k==1) then
                    !     write(*,*) "------------------------------------------------------------------"
                    !     write(*,*) "------------------------------------------------------------------"
                    !     write(*,'(A,I3,A,I3,A,I3)') "   d_d_ls panel ", y, ", vert ",p,&
                    !     ", xi coordinate, cp=",z 
                    ! else
                    !     write(*,*) "------------------------------------------------------------------"
                    !     write(*,*) "------------------------------------------------------------------"
                    !     write(*,'(A,I3,A,I3,A,I3)') "   d_d_ls panel ", y, ", vert ",p,&
                    !     ", eta coordinate, cp=",z 
                    ! end if 
                
                    ! write(*,*) ""
                    ! write(*,'(A,f16.10)') "        Max Residual = ", maxval(abs(residuals))
                    if (maxval(abs(residuals))>error_allowed) then
                        write(*,*) "          Failed values :"
                        write(*,*) "          d_d_ls_FD            adjoint d_d_ls              residuals"
                        do i = 1, N_verts*3
                            if (abs(residuals(i))>error_allowed) then
                                write(*, '(8x,(f14.10, 4x),3x, (f14.10, 4x),3x, (f14.10, 4x))')&
                                d_d_ls_FD(i), adjoint_geom%d_d_ls(k,p)%get_value(i), residuals(i)
                            end if
                        end do
                    end if

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
                        if (k==1 .and. p ==1) then
                            failure_log(total_tests-passed_tests) = "calc d_d_ls vert 1, &
                            xi coordinate test FAILED"
                        elseif (k==2 .and. p ==1) then
                            failure_log(total_tests-passed_tests) = "calc d_d_ls vert 1, &
                            eta coordinate test FAILED"
                        elseif (k==1 .and. p ==2) then
                            failure_log(total_tests-passed_tests) = "calc d_d_ls vert 2, &
                            xi coordinate test FAILED"
                        elseif (k==2 .and. p ==2) then
                            failure_log(total_tests-passed_tests) = "calc d_d_ls vert 2, &
                            eta coordinate test FAILED"
                        elseif (k==1 .and. p ==3) then
                            failure_log(total_tests-passed_tests) = "calc d_d_ls vert 3, &
                            xi coordinate test FAILED"
                        elseif (k==2 .and. p ==3) then
                            failure_log(total_tests-passed_tests) = "calc d_d_ls vert 3, &
                            eta coordinate test FAILED"
                        else
                            failure_log(total_tests-passed_tests) = "calc d_d_ls FAILED"
                        end if
                        write(*,*) failure_log(total_tests-passed_tests)
                    else
                        ! if (k==1) then
                        !     write(*,'(A,I3,A,I3,A,I3,A)') "        calc d_d_ls panel ", y, ", vert ",p,&
                        !     ", xi coordinate, cp=",z, " test PASSED"
                        ! else
                        !     write(*,'(A,I3,A,I3,A,I3,A)') "        calc d_d_ls panel ", y, ", vert ",p,&
                        !     ", eta coordinate, cp=",z, " test PASSED"
                        ! end if
                        ! write(*,*) "" 
                        ! write(*,*) ""
                        passed_tests = passed_tests + 1
                        total_tests = total_tests + 1
                        
                    end if
                    test_failed = .false.
                    
                
                end do
            ! k loop
            end do
        
        end do ! z loop (control points)
    
    end do !y loop (panels)

    !!!!!!!!!!!!!! CALC_BASIC_GEOM (CONTROL POINTS) SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "     CALC_BASIC_GEOM (CONTROL POINTS) SENSITIVITIES TEST RESULTS "
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
    write(*,'(A,f12.10, A)') " Total test time = ", time, " minutes"
    write(*,*) ""
    write(*,*) "----------------------"
    write(*,*) "Program Complete"
    write(*,*) "----------------------"

end program calc_basic_geom_cp