program dirichlet_super_test9
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

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals, hH113_up, hH113_dn, d_hH113_FD, &
    H213_up, H213_dn, d_H213_FD, H123_up, H123_dn, d_H123_FD
    real,dimension(:,:),allocatable ::  residuals3 

    integer :: i,j,k,m,n,y,z, N_verts, N_panels, vert, index, cp_ind
    real :: step,error_allowed, cp_offset
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed, mirror_panel
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

    test_input = "dev\input_files\adjoint_inputs\dirichlet_supersonic_test.json"
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
    test_int = test_mesh%panels(index)%calc_integrals(test_geom, 'potential', freestream_flow,.false., test_dod_info)
    !!!!!!!!!!!!!!!!!!!!! END TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call system_clock(start_count, count_rate)
    

    !!!!!!!!!!!!!!!!!!!!!!ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()

    adjoint_input = "dev\input_files\adjoint_inputs\dirichlet_supersonic_adjoint_test.json"
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

    allocate(hH113_up(N_verts*3))
    allocate(hH113_dn(N_verts*3))
    allocate(d_hH113_FD(N_verts*3))
    allocate(H213_up(N_verts*3))
    allocate(H213_dn(N_verts*3))
    allocate(d_H213_FD(N_verts*3))
    allocate(H123_up(N_verts*3))
    allocate(H123_dn(N_verts*3))
    allocate(d_H123_FD(N_verts*3))

    
    mirror_panel = .false.
    error_allowed = 1.0e-5
    step = 0.000001
    index = 1
    cp_ind = 1
    

    write(*,*) ""
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) "       Dirichlet Supersonic H Integrals SENSITIVITIES TEST                    "
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""


    
    
    do y =1,N_panels
        index = y
        
        do z = 1,N_verts
            cp_ind = z
            
            test_dod_info = test_mesh%panels(index)%check_dod(test_mesh%cp(cp_ind)%loc, freestream_flow, mirror_panel)
            if (test_dod_info%in_dod .and. test_mesh%panels(index)%A > 0.) then

                ! panel is in domain of dependence of point    
                write(*,'(A,I5,A,I5)') "Supersonic Subinc H integrals test: Panel ",y," cp ",z

                ! Calculate geometric parameters
                if (freestream_flow%supersonic) then
                    if ((mirror_panel .and. test_mesh%panels(index)%r_mir < 0.) .or. &
                                    (.not. mirror_panel .and. test_mesh%panels(index)%r < 0.)) then
                        !geom = this%calc_supersonic_supinc_geom(P, freestream, mirror_panel, dod_info)
                        write(*,*) " can't do super inclined"
                        stop
                    else
                        ! supersonic subinclined
                        adjoint_geom = adjoint_mesh%panels(index)%calc_supersonic_subinc_geom_adjoint(&
                                                    adjoint_mesh%cp(cp_ind)%loc,adjoint_mesh%cp(cp_ind)%d_loc,&
                                                    adjoint_freestream_flow, mirror_panel, test_dod_info)
                    end if
                else
                    !geom = this%calc_subsonic_geom(P, freestream, mirror_panel)
                    write(*,*) " something went wrong, shouldn't be at calc subsonic geom"
                
                end if ! end supersonic if statement

                ! Get integrals
                adjoint_int = adjoint_mesh%panels(index)%calc_integrals(adjoint_geom, 'potential', &
                                            adjoint_freestream_flow, mirror_panel, test_dod_info)

                ! ! integral adjoint
                call adjoint_mesh%panels(index)%calc_integrals_adjoint(&
                adjoint_geom,"potential", adjoint_int,adjoint_freestream_flow, .false., adjoint_dod_info)
        
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
                        deallocate(test_solver%i_sigma_in_sys)
                        deallocate(test_solver%i_sys_sigma_in_body)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P)
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)

                        ! update hH113
                        deallocate(test_int%F111)
                        test_geom = test_mesh%panels(index)%calc_supersonic_subinc_geom(&
                                                test_mesh%cp(cp_ind)%loc,freestream_flow,mirror_panel,test_dod_info)
                        test_int = test_mesh%panels(index)%calc_integrals(test_geom, 'potential', freestream_flow,&
                                                .false., test_dod_info)
                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                        
                        ! get desired info
                        hH113_up(j + (i-1)*N_verts) = test_int%hH113
                        H213_up(j + (i-1)*N_verts) = test_int%H213
                        H123_up(j + (i-1)*N_verts) = test_int%H123

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
                        deallocate(test_solver%i_sigma_in_sys)
                        deallocate(test_solver%i_sys_sigma_in_body)
                        deallocate(test_mesh%cp)
                        deallocate(test_solver%P)
                        call test_solver%init(solver_settings, processing_settings, &
                        test_mesh, freestream_flow, control_point_file)


                        ! update hH113
                        deallocate(test_int%F111)
                        deallocate(test_int%F121)
                        deallocate(test_int%F211)
                        test_geom = test_mesh%panels(index)%calc_supersonic_subinc_geom(&
                                                test_mesh%cp(cp_ind)%loc,freestream_flow,mirror_panel,test_dod_info)
                        test_int = test_mesh%panels(index)%calc_integrals(test_geom, 'potential', freestream_flow,&
                                                .false., test_dod_info)
                        !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                        ! get desired info
                        hH113_dn(j + (i-1)*N_verts) = test_int%hH113
                        H213_dn(j + (i-1)*N_verts) = test_int%H213
                        H123_dn(j + (i-1)*N_verts) = test_int%H123
                        
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    end do 
                end do 
                
                ! central difference 
                d_hH113_FD = (hH113_up - hH113_dn)/(2.*step)
                d_H213_FD = (H213_up - H213_dn)/(2.*step)
                d_H123_FD = (H123_up - H123_dn)/(2.*step)
                
            
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_hH113 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

                ! calculate residuals
                do i =1, N_verts*3
                    residuals(i) = adjoint_int%d_hH113%get_value(i) - d_hH113_FD(i)
                end do

                if (maxval(abs(residuals(:)))>error_allowed) then
                    write(*,*) ""
                    write(*,*) "     FLAGGED VALUES :"
                    do i = 1, N_verts*3
                        if (abs(residuals(i))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "                                                d_hH113 "
                            write(*, '(A25,8x,(f25.10, 4x))') "    Central Difference", d_hH113_FD(i)
                        
                            write(*, '(A25,8x,(f25.10, 4x))') "          adjoint",   &
                            adjoint_int%d_hH113%get_value(i)
                            write(*, '(A25,8x,(f25.10, 4x))') "              residual", residuals(i)
                        end if
                    end do
                end if
        
                
                
                ! check if test failed
                do i=1,N_verts*3
                    if (abs(residuals(i)) > error_allowed) then 
                        
                        if (abs(d_hH113_FD(i))>1000.0) then
                            if (abs(residuals(i)) > error_allowed*10000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (1000.0>abs(d_hH113_FD(i)) .and. abs(d_hH113_FD(i))>100.0) then
                            if (abs(residuals(i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_hH113_FD(i)) .and. abs(d_hH113_FD(i))>10.0) then
                            if (abs(residuals(i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_hH113_FD(i)) .and. abs(d_hH113_FD(i))>1.0) then
                            if (abs(residuals(i)) > error_allowed*10.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        else
                            if (abs(residuals(i)) > error_allowed) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        end if
                    end if
                end do
                if (test_failed) then
                    total_tests = total_tests + 1
                    write(*,'(A,I5,A,I5,A)')"                                               &
                                    d_hH113 panel ",y," cp ",z," test FAILED"
                    failure_log(total_tests-passed_tests) = "d_hH113 test FAILED"
                else
                    ! write(*,*) "        CALC d_hH113 test PASSED"
                    ! write(*,*) "" 
                    ! write(*,*) ""
                    passed_tests = passed_tests + 1
                    total_tests = total_tests + 1
                    
                end if
                test_failed = .false.



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_H213 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

                ! calculate residuals
                do i =1, N_verts*3
                    residuals(i) = adjoint_int%d_H213%get_value(i) - d_H213_FD(i)
                end do

                if (maxval(abs(residuals(:)))>error_allowed) then
                    write(*,*) ""
                    write(*,*) "     FLAGGED VALUES :"
                    do i = 1, N_verts*3
                        if (abs(residuals(i))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "                                             d_H213"
                            write(*, '(A25,8x,(f25.10, 4x))') "    Central Difference", d_H213_FD(i)
                        
                            write(*, '(A25,8x,(f25.10, 4x))') "          adjoint",   &
                            adjoint_int%d_H213%get_value(i)
                            write(*, '(A25,8x,(f25.10, 4x))') "              residual", residuals(i)
                        end if
                    end do
                end if
        
                
                
                ! check if test failed
                do i=1,N_verts*3
                    if (abs(residuals(i)) > error_allowed) then 
                        
                        if (abs(d_H213_FD(i))>1000.0) then
                            if (abs(residuals(i)) > error_allowed*10000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (1000.0>abs(d_H213_FD(i)) .and. abs(d_H213_FD(i))>100.0) then
                            if (abs(residuals(i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_H213_FD(i)) .and. abs(d_H213_FD(i))>10.0) then
                            if (abs(residuals(i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_H213_FD(i)) .and. abs(d_H213_FD(i))>1.0) then
                            if (abs(residuals(i)) > error_allowed*10.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        else
                            if (abs(residuals(i)) > error_allowed) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        end if
                    end if
                end do
                if (test_failed) then
                    total_tests = total_tests + 1
                    write(*,'(A,I5,A,I5,A)')"                                               &
                                    d_H213 panel ",y," cp ",z," test FAILED"
                    failure_log(total_tests-passed_tests) = "d_H213 test FAILED"
                else
                    ! write(*,*) "        CALC d_H213 test PASSED"
                    ! write(*,*) "" 
                    ! write(*,*) ""
                    passed_tests = passed_tests + 1
                    total_tests = total_tests + 1
                    
                end if
                test_failed = .false.




                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_H123 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

                ! calculate residuals
                do i =1, N_verts*3
                    residuals(i) = adjoint_int%d_H123%get_value(i) - d_H123_FD(i)
                end do

                if (maxval(abs(residuals(:)))>error_allowed) then
                    write(*,*) ""
                    write(*,*) "     FLAGGED VALUES :"
                    do i = 1, N_verts*3
                        if (abs(residuals(i))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "                                            d_H123 "
                            write(*, '(A25,8x,(f25.10, 4x))') "    Central Difference", d_H123_FD(i)
                        
                            write(*, '(A25,8x,(f25.10, 4x))') "          adjoint",   &
                            adjoint_int%d_H123%get_value(i)
                            write(*, '(A25,8x,(f25.10, 4x))') "              residual", residuals(i)
                        end if
                    end do
                end if
        
                
                
                ! check if test failed
                do i=1,N_verts*3
                    if (abs(residuals(i)) > error_allowed) then 
                        
                        if (abs(d_H123_FD(i))>1000.0) then
                            if (abs(residuals(i)) > error_allowed*10000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (1000.0>abs(d_H123_FD(i)) .and. abs(d_H123_FD(i))>100.0) then
                            if (abs(residuals(i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_H123_FD(i)) .and. abs(d_H123_FD(i))>10.0) then
                            if (abs(residuals(i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_H123_FD(i)) .and. abs(d_H123_FD(i))>1.0) then
                            if (abs(residuals(i)) > error_allowed*10.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        else
                            if (abs(residuals(i)) > error_allowed) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        end if
                    end if
                end do
                if (test_failed) then
                    total_tests = total_tests + 1
                    write(*,'(A,I5,A,I5,A)')"                                               &
                                    d_H123 panel ",y," cp ",z," test FAILED"
                    failure_log(total_tests-passed_tests) = "d_H123 test FAILED"
                else
                    ! write(*,*) "        CALC d_H123 test PASSED"
                    ! write(*,*) "" 
                    ! write(*,*) ""
                    passed_tests = passed_tests + 1
                    total_tests = total_tests + 1
                    
                end if
                test_failed = .false.

            else ! end check dod if statement    
                write(*,'(A,I5,A,I5)') "Panel ",y," not in dod of cp ",z
            end if

        ! z loop
        end do

    ! y loop
    end do


    !!!!!!!!!!!!!!  RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "      Dirichlet SUPERSONIC H INTEGRAL SENSITIVITIES TEST RESULTS "
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

end program dirichlet_super_test9