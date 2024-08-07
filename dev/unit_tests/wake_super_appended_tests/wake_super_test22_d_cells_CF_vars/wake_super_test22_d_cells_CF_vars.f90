program wake_super_test22
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
    integer :: i_unit, stat
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals
    real,dimension(:,:),allocatable ::  residuals3, cells_CF_up, cells_CF_dn, d_cells_CF_FD

    integer :: i,j,k,m,n,y,z,N_original_verts,N_total_verts, N_panels, vert, index, cp_ind
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

    integer :: adjoint_solver_stat, test_solver_stat
    type(sparse_vector) :: zeros

    real,dimension(3) :: adjoint_P, test_P, test_v_d, test_v_s
    type(sparse_matrix) :: adjoint_d_P_term2
    type(sparse_matrix) :: adjoint_d_P
    type(sparse_matrix) :: adjoint_d_v_d_panel
    
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

    test_input = "dev\input_files\adjoint_inputs\wake_super_test.json"
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
    test_mesh%perturb_point = .true.

    N_original_verts = test_mesh%N_verts

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

    call system_clock(start_count, count_rate)


    !!!!!!!!!!!!!!!!!!!!!!ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()
    
    adjoint_input = "dev\input_files\adjoint_inputs\wake_super_adjoint_test.json"
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
    ! solve
    call adjoint_solver%solve(adjoint_mesh, adjoint_solver_stat, adjoint_formulation,adjoint_freestream_flow)
    

    
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!

    
    N_total_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_original_verts*3))
    allocate(residuals(N_original_verts*3))

    ! allocate data holders
    allocate(cells_CF_up(3,N_original_verts*3))
    allocate(cells_CF_dn(3,N_original_verts*3))
    allocate(d_cells_CF_FD(3,N_original_verts*3))

    

    error_allowed = 1.0e-8
    step = 0.000001
    index = 1
    cp_ind = 1
    

    write(*,*) ""
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) "       SUPERSONIC dirichlet d_cells_CF_wrt_vars TEST (WAKE PRESENT)                   "
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""


    
    
    do z =1,N_panels

        write(*,*) ""
        write(*,*) "--------------------------------------------------------------------------------------"
        write(*,'(A,I5)') "       SUPERSONIC dirichlet d_cells_CF_wrt_vars (WAKE PRESENT) test ",z
        write(*,*) "--------------------------------------------------------------------------------------"
        write(*,*) ""
        
        
        do i=1,3
            do j=1,N_original_verts

                deallocate(test_mesh%vertices, test_mesh%edges, test_mesh%panels, test_mesh%vertex_ordering)
                call test_mesh%init(geom_settings)
                test_mesh%perturb_point = .true.

                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                ! update panel geometry and calc
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do

                call test_mesh%calc_vertex_geometry()
                
                call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)
                
                ! update solver init
                deallocate(test_solver%sigma_known)
                deallocate(test_solver%i_sigma_in_sys)
                deallocate(test_solver%i_sys_sigma_in_body)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)
                
                ! Check for errors
                if (test_solver_stat /= 0) return

                deallocate(test_mesh%V_cells_inner, test_mesh%V_cells)

                ! Calculate velocities
                call test_solver%calc_cell_velocities(test_mesh)

                deallocate(test_mesh%C_P_ise)

                ! Calculate velocities
                call test_solver%calc_pressures(test_mesh)

                deallocate(test_mesh%dC_f)

                call test_solver%calc_forces(test_mesh)
                                
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                cells_CF_up(:,j + (i-1)*N_original_verts) = test_mesh%dC_F(:,z)
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!! UPDATE STEP DOWN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                deallocate(test_mesh%vertices, test_mesh%edges, test_mesh%panels, test_mesh%vertex_ordering)
                call test_mesh%init(geom_settings)
                test_mesh%perturb_point = .true.

                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - step
                
                ! update panel geometry and calc
                do m =1,N_panels
                    deallocate(test_mesh%panels(m)%n_hat_g)
                    call test_mesh%panels(m)%calc_derived_geom()
                end do

                call test_mesh%calc_vertex_geometry()
                
                call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)
                
                ! update solver init
                deallocate(test_solver%sigma_known)
                deallocate(test_solver%i_sigma_in_sys)
                deallocate(test_solver%i_sys_sigma_in_body)
                deallocate(test_mesh%cp)
                deallocate(test_solver%P)
                call test_solver%init(solver_settings, processing_settings, &
                test_mesh, freestream_flow, control_point_file)

                ! Check for errors
                if (test_solver_stat /= 0) return

                deallocate(test_mesh%V_cells_inner, test_mesh%V_cells)

                ! Calculate velocities
                call test_solver%calc_cell_velocities(test_mesh)

                deallocate(test_mesh%C_P_ise)

                ! Calculate velocities
                call test_solver%calc_pressures(test_mesh)

                deallocate(test_mesh%dC_f)

                call test_solver%calc_forces(test_mesh)
                                
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! get the needed info
                cells_CF_dn(:, j + (i-1)*N_original_verts) = test_mesh%dC_F(:,z)

            
            end do 
        end do 

            
            ! central difference 
        d_cells_CF_FD = (cells_CF_up - cells_CF_dn)/(2.*step)
                
        
        do i=1,N_original_verts*3
            residuals3(:,i) = adjoint_mesh%d_cell_forces_wrt_vars(z)%get_values(i) - d_cells_CF_FD(:,i)
        end do
            

        if (maxval(abs(residuals3(:,:)))>error_allowed) then
            write(*,*) ""
            write(*,*) "     FLAGGED VALUES :"
            do i = 1, N_original_verts*3
                if (any(abs(residuals3(:,i))>error_allowed)) then
                    write(*,*) ""
                    write(*,'(A,I5,A)') "                                       d_cells_CF_wrt_vars &
                     ",z,"                                             residuals"
                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_cells_CF_FD(:,i)
                
                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                    adjoint_mesh%d_cell_forces_wrt_vars(z)%get_values(i), residuals3(:,i)
                end if
            end do
        end if

        
        
        ! check if test failed
        do i=1,N_original_verts*3
            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                do j = 1,3
                    if (abs(d_cells_CF_FD(j,i))>1000.0) then
                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (1000.0>abs(d_cells_CF_FD(j,i)) .and. abs(d_cells_CF_FD(j,i))>100.0) then
                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (100.0>abs(d_cells_CF_FD(j,i)) .and. abs(d_cells_CF_FD(j,i))>10.0) then
                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (10.0>abs(d_cells_CF_FD(j,i)) .and. abs(d_cells_CF_FD(j,i))>1.0) then
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
            write(*,'(A,I5,A)')"                                               &
            d_cells_CF_wrt_vars  ",z," test FAILED"
            failure_log(total_tests-passed_tests) = "d_cells_CF_wrt_vars test FAILED"
        else
            ! write(*,*) "        d_cells_CF_wrt_vars test PASSED"
            ! write(*,*) "" 
            ! write(*,*) ""
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if

        ! reset test failed for the next z loop
        test_failed = .false.


        
        

    ! z loop
    end do


    !!!!!!!!!!!!!!  RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "     SUPERSONIC dirichlet d_cells_CF_wrt_vars TEST RESULTS (WAKE PRESENT) "
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

end program wake_super_test22