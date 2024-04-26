program test11
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
    integer :: i_unit
    logical :: exists, found
    real,dimension(:),allocatable :: doublet_inf
    real,dimension(:,:),allocatable :: v_s, v_d
    type(sparse_vector),dimension(3) :: inf_adjoint, inf_adjoint2

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals
    real,dimension(:,:),allocatable ::  residuals3, inf_up, inf_dn, d_inf_FD

    integer :: i,j,k,m,n,y,z, N_verts, N_panels, vert, index, cp_ind
    real :: step,error_allowed, cp_offset
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute
    
    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(100) :: failure_log
    character(len=10) :: m_char
    integer(8) :: start_count, end_count
    real(16) :: count_rate, time
    
    !!!!!!!!!!!!!!!!!!! END TESTING STUFF !!!!!!!!!!!!!!!!!!!!!11
    
    test_failed = .false. 
    passed_tests = 0
    total_tests = 0
    
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
    
      
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    allocate(inf_up(3,N_verts*3))
    allocate(inf_dn(3,N_verts*3))
    allocate(d_inf_FD(3,N_verts*3))
    

    error_allowed = 1.0e-7
    step = 0.000001
    index = 1
    cp_ind = 1
    

    write(*,*) ""
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) "                           inf_adjoint2 TEST                    "
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""


    
    
    do y =1,N_panels
        index = y
        
        do z = 1,N_verts
            cp_ind = z
            
            write(*,'(A,I5,A,I5)') "inf adjoint2 test: Panel ",y," cp ",z

        
            ! inf adjoint method 2
            inf_adjoint2 = adjoint_mesh%panels(index)%calc_doublet_inf_adjoint2(&
            adjoint_mesh%cp(cp_ind),adjoint_freestream_flow)
            

                
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
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    deallocate(test_mesh%panels(index)%i_vert_d)
                    deallocate(test_mesh%panels(index)%S_mu_inv)
                    deallocate(test_mesh%panels(index)%T_mu)
                    ! deallocate(test_mesh%panels(index)%i_panel_s)
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
                    call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc,&
                     freestream_flow,.false.,v_s, v_d)
                    doublet_inf = matmul(test_mesh%cp(cp_ind)%n_g, matmul(freestream_flow%B_mat_g, v_d))
                    !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                    
                    ! get the needed info
                    inf_up(:, j + (i-1)*N_verts) = doublet_inf(:)
                    
                    
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
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    deallocate(test_mesh%panels(index)%i_vert_d)
                    deallocate(test_mesh%panels(index)%S_mu_inv)
                    deallocate(test_mesh%panels(index)%T_mu)
                    ! deallocate(test_mesh%panels(index)%i_panel_s)
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
                    call test_mesh%panels(index)%calc_velocity_influences(test_mesh%cp(cp_ind)%loc,&
                     freestream_flow,.false.,v_s, v_d)
                    doublet_inf = matmul(test_mesh%cp(cp_ind)%n_g, matmul(freestream_flow%B_mat_g, v_d))
                    !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                    
                    ! get the needed info
                    inf_dn(:, j + (i-1)*N_verts) = doublet_inf(:)

                    
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                end do 
            end do 
            
            ! central difference 
            d_inf_FD(:,:) = (inf_up(:,:) - inf_dn(:,:))/(2.*step)

            ! calculate residuals3
            do i =1, N_verts*3
                residuals3(:,i) = (/inf_adjoint2(1)%get_value(i),&
                                  inf_adjoint2(2)%get_value(i),&
                                  inf_adjoint2(3)%get_value(i)/) - d_inf_FD(:,i)
            end do


            
            if (maxval(abs(residuals3(:,:)))>error_allowed) then
                write(*,*) ""
                write(*,*) "     FLAGGED VALUES :"
                do i = 1, N_verts*3
                    if (any(abs(residuals3(:,i))>error_allowed)) then
                        write(*,*) ""
                        write(*,'(A,I5,A,I5,A)') "                                       inf_adjoint2 &
                        panel ", y,"on cp ",z,"                                              residuals"
                        write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_inf_FD(:,i)
                    
                        write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                inf_adjoint2(1)%get_value(i),&
                                inf_adjoint2(2)%get_value(i),&
                                inf_adjoint2(3)%get_value(i), residuals3(:,i)
                    end if
                end do
            end if

            
            
            ! check if test failed
            do i=1,N_verts*3
                if (any(abs(residuals3(:,i)) > error_allowed)) then 
                    do j = 1,3
                        if (abs(d_inf_FD(j,i))>1000.0) then
                            if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (1000.0>abs(d_inf_FD(j,i)) .and. abs(d_inf_FD(j,i))>100.0) then
                            if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_inf_FD(j,i)) .and. abs(d_inf_FD(j,i))>10.0) then
                            if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_inf_FD(j,i)) .and. abs(d_inf_FD(j,i))>1.0) then
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
                inf_adjoint2 of panel ",y," on cp ",z," test FAILED"
                failure_log(total_tests-passed_tests) = "inf_adjoint2 test FAILED"
            else
                ! write(*,*) "        inf_adjoint2 test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


        ! z loop
        end do

    ! y loop
    end do


    !!!!!!!!!!!!!!  RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "                          inf_adjoint2 TEST RESULTS "
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
    write(*,'(A,f16.10, A)') " Total test time = ", time, " minutes"
    write(*,*) ""
    write(*,*) "----------------------"
    write(*,*) "Program Complete"
    write(*,*) "----------------------"

end program test11