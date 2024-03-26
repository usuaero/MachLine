program calc_inf_point
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
    type(sparse_vector),dimension(3) :: inf_adjoint

    real,dimension(3) :: adjoint_P, test_P, test_v_d, test_v_s
    real,dimension(:,:),allocatable :: test_v_d_inf, test_v_s_inf
    type(sparse_matrix) :: adjoint_d_P_term2
    type(sparse_matrix) :: adjoint_d_P, adjoint_geom_dl1
    type(sparse_matrix),dimension(3) ::  adjoint_d_v_d_inf_row
    type(sparse_vector),dimension(3,3) :: adjoint_d_v_d_M


    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals
    real,dimension(:,:),allocatable ::  residuals3, inf_up, inf_dn, d_inf_FD

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

    test_P = test_mesh%panels(index)%centr - 1.e-10*test_mesh%panels(index)%n_g
    
    test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_P, freestream_flow, .false.)
    write(*,*) " geom a = ", test_geom%a
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

    !!!! inf adjoint method 2
    ! inf_adjoint = adjoint_mesh%panels(index)%calc_doublet_inf_adjoint2(adjoint_mesh%cp(cp_ind), &
    ! adjoint_freestream_flow)
    !!!!!!!!!

    adjoint_P = adjoint_mesh%panels(index)%centr - 1.e-10*adjoint_mesh%panels(index)%n_g

    ! calc d_P
    call adjoint_d_P%init_from_sparse_matrix(adjoint_mesh%panels(index)%d_centr) 
    call adjoint_d_P_term2%init_from_sparse_matrix(adjoint_mesh%panels(index)%d_n_g)
    call adjoint_d_P_term2%broadcast_element_times_scalar(-1.e-10)
    call adjoint_d_P%sparse_add(adjoint_d_P_term2)

    ! ! inf adjoint method 1
    ! call adjoint_mesh%panels(index)%calc_velocity_influences_adjoint(adjoint_P, adjoint_d_P, &
    ! adjoint_freestream_flow, adjoint_d_v_d_M)

    ! do i = 1,3
    !     call adjoint_d_v_d_inf_row(i)%init_from_sparse_vectors(adjoint_d_v_d_M(i,1), adjoint_d_v_d_M(i,2), &
    !     adjoint_d_v_d_M(i,3))
    ! end do
    ! Perform flow-dependent initialization on the surface mesh

    ! calc CALC BASIC GEOM geom sensitivity of relation between cp1 and panel1 
    adjoint_geom = adjoint_mesh%panels(index)%calc_subsonic_geom_adjoint(adjoint_P,&
    adjoint_d_P, adjoint_freestream_flow)

    call adjoint_geom_dl1%init_from_sparse_vectors(adjoint_geom%d_g2(1), adjoint_geom%d_g2(2),adjoint_geom%d_g2(3))
    write(*,*) " adjoint geom a = ", adjoint_geom%a
    
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
    allocate(inf_up(3,N_verts*3))
    allocate(inf_dn(3,N_verts*3))
    allocate(d_inf_FD(3,N_verts*3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC inf_adjoint (panel x, cp x) row 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    ! panel vertex 
    ! do k=1,3

        write(*,'(A, I1,A,I1,A,I1,A)') "---------------------------------- TEST CALC inf_adjoint &
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
                
                test_P = test_mesh%panels(index)%centr - 1.e-10*test_mesh%panels(index)%n_g
                ! update v_d and doublet inf
                ! call test_mesh%panels(index)%calc_velocity_influences(test_P, freestream_flow,.false.,v_s, v_d)
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_P,freestream_flow,.false.)
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                inf_up(:,j + (i-1)*N_verts) = test_geom%a*test_geom%a
                
                
                ! perturb down the current design variable
                ! write(*,*) " perturb down"
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                
                !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
               
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
                
                test_P = test_mesh%panels(index)%centr - 1.e-10*test_mesh%panels(index)%n_g
                ! update v_d and doublet inf
                ! call test_mesh%panels(index)%calc_velocity_influences(test_P, freestream_flow,.false.,v_s, v_d)
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_P,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                inf_dn(:,j + (i-1)*N_verts) = test_geom%a*test_geom%a

                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_inf_FD = (inf_up - inf_dn)/(2.*step)

        write(*,*) ""
        
        ! write results
        write(*,*) "-----------------------------------------------------------------------------------------"
        write(*,*) ""
        write(*,'(A,I1,A,I1)') "                Central Diff d_v_d_inf panel ", index, " row = ",k
        write(*,*) "  d_v_d_inf_x           d_v_d_inf_y           d_v_d_inf_z "
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x))') d_inf_FD(:,i)
        end do 
        
        ! !write sparse matrix
        ! write(*,*) ""
        ! write(*,'(A,I1)') "                d_V_inner_test index = ",index
        ! write(*,*) "  d_V_inner_x           d_V_inner_y           d_V_inner_z     sparse_index       full_index"
        ! write(*,*) "sparse num cols", adjoint_mesh%d_V_cells_inner(k)%sparse_num_cols
        ! do i=1,adjoint_mesh%d_V_cells_inner(k)%sparse_num_cols
        !     write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%d_V_cells_inner(k)%columns(i)%vector_values, &
        !     i, adjoint_mesh%d_V_cells_inner(k)%columns(i)%full_index
        ! end do
        ! write(*,*) ""


        
        do i=1,N_verts*3
            residuals3(:,i) = adjoint_geom_dl1%get_values(i) - d_inf_FD(:,i)
        end do
        
        write(*,'(A,I1,A,I1)') "                adjoint d_v_d_inf panel ", index, " row = ",k
        write(*,*) "  d_v_d_inf_x           d_v_d_inf_y           d_v_d_inf_z                        residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_geom_dl1%get_values(i), residuals3(:,i)
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
            
            failure_log(total_tests-passed_tests) = "TEST FAILURE in calc_vel_inf_test"
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A, I1)') " calc_vel_inf Panel ", index, " row ", k," TEST PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
        
        
    ! end do




    !!!!!!!!!!!!!! calc_vel_inf_test RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------calc_vel_inf_test RESULTS--------------"
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

end program calc_inf_point