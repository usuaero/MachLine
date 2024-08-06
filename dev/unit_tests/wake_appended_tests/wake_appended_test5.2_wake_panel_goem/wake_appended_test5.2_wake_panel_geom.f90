program wake_appended_test5_2

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
    integer :: i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, area_up, area_dn, d_area_FD

    real,dimension(:,:),allocatable :: v, residuals3 , normal_up, normal_dn, d_n_g_FD, centr_up, centr_dn,&
    d_centr_FD, n_hat_g_up, n_hat_g_dn, d_n_hat_g_FD

    integer :: i,j,k,m,n,p,strip,y,z, N_original_verts, N_total_verts, N_panels, vert, &
    num_wake_strips, num_panels_in_strip, index, cp_ind, clones, vert_ind
    real :: step,error_allowed
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!! TEST INPUT (calc_adjoint = false) !!!!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()
    
    test_input = "dev\input_files\adjoint_inputs\wake_appended_test.json"
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
    
    ! set perturb point to true so the clone dir is based on centroid information
    
    !!!!!!!!!!!!!!!!!!!!! END TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call system_clock(start_count, count_rate)
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()
    
    adjoint_input = "dev\input_files\adjoint_inputs\wake_appended_adjoint_test.json"
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
    write(*,*) "made it here"
    call adjoint_mesh%init_with_flow(adjoint_freestream_flow, adjoint_body_file, adjoint_wake_file, adjoint_formulation)
    
    ! Initialize panel solver
    ! call adjoint_solver%init(adjoint_solver_settings, adjoint_processing_settings, adjoint_mesh, &
    ! adjoint_freestream_flow, adjoint_control_point_file)
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!
    
    
    ! do j = 1,188
    !     d_loc_matrix = adjoint_mesh%cp(j)%d_loc%expand(.true.)
    !     allocate(norm(3))
    !     do i =1,3
    !         norm(i) = sqrt(sum(d_loc_matrix(:,i)*d_loc_matrix(:,i)))
    !         write(*,*) "Norm of CP_loc = ", norm(i)
    !     end do
    !     deallocate(norm)
    ! end do

    ! stop

    N_total_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    
    allocate(centr_up(3,N_original_verts*3))
    allocate(centr_dn(3,N_original_verts*3))
    allocate(d_centr_FD(3,N_original_verts*3))
    allocate(normal_up(3,N_original_verts*3))
    allocate(normal_dn(3,N_original_verts*3))
    allocate(d_n_g_FD(3,N_original_verts*3))
    allocate(area_up(N_original_verts*3))
    allocate(area_dn(N_original_verts*3))
    allocate(d_area_FD(N_original_verts*3))
    allocate(n_hat_g_up(3,N_original_verts*3))
    allocate(n_hat_g_dn(3,N_original_verts*3))
    allocate(d_n_hat_g_FD(3,N_original_verts*3))

    allocate(residuals3(3,N_original_verts*3))
    allocate(residuals(N_original_verts*3))

    error_allowed = 1.0e-9
    step = 0.000001
    index = 1
    cp_ind = 1
    

    write(*,*) ""
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) "       SUBsonic Wake Panel Geom Test (WAKE APPENDED)                   "
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""

    num_wake_strips = test_mesh%wake%N_strips

    ! do for each wake strip
    do strip = 1, num_wake_strips

        write(*,'(A,I5)') "Test Wake Strip ", strip
        num_panels_in_strip = test_mesh%wake%strips(strip)%N_panels

        ! do for each panel in the wake strip
        do p = 1, num_panels_in_strip

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CP_d_loc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(*,'(A,I5)') "    Test Panel", p
            write(*,*) ""
        
    

            ! do for each design variable
            do i=1,3
                do j=1,N_original_verts

                    deallocate(test_mesh%vertices, test_mesh%edges, test_mesh%panels, test_mesh%vertex_ordering, &
                    test_mesh%wake%strips)
                    call test_mesh%init(geom_settings)
                    test_mesh%perturb_point = .true.

                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    ! write(*,*) "this vertex is clone? ", test_mesh%vertices(j)%clone 
                    
                    !!!!!!!!!!! update !!!!!!!!!!!!!
                    ! update panel geometry and calc
                    do m =1,N_panels
                        deallocate(test_mesh%panels(m)%n_hat_g)
                        call test_mesh%panels(m)%calc_derived_geom()
                    end do

                    call test_mesh%calc_vertex_geometry()
                        
                    call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)
                    

                    !!!!!!!!!!!!! end update !!!!!!!!!!!!!!!!
                    
                    !!!!!!!!!!!!!!!!!!!!!!! GET INFO !!!!!!!!!!!!!!!!!!!!!!!!     
                    ! put the x y or z component of the vertex of interest (index) in a list
                    normal_up(:,j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%n_g(:)
                    centr_up(:,j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%centr(:)
                    area_up(j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%A


                    !!!! Perturb Down !!!!

                    deallocate(test_mesh%vertices, test_mesh%edges, test_mesh%panels, test_mesh%vertex_ordering, &
                    test_mesh%wake%strips)
                    call test_mesh%init(geom_settings)
                    test_mesh%perturb_point = .true.

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) -step
                        
                    !!!!!!!!!!! update !!!!!!!!!!!!!!!
                    ! update panel geometry and calc
                    do m =1,N_panels
                        deallocate(test_mesh%panels(m)%n_hat_g)
                        call test_mesh%panels(m)%calc_derived_geom()
                    end do

                    call test_mesh%calc_vertex_geometry()
                        
                    call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)
                        
                    !!!!!!!!!!!!!!!! end update !!!!!!!!!!!!!!!!!!!!!

                    ! put the x y or z component of the vertex of interest (cp_ind) in a list
                    normal_dn(:,j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%n_g(:)
                    centr_dn(:, j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%centr(:)
                    area_dn(j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%A

                end do 
            end do 
            
            ! central difference 

            d_n_g_FD = (normal_up - normal_dn)/(2.*step)
            d_centr_FD = (centr_up - centr_dn)/(2.*step)
            d_area_FD = (area_up - area_dn)/(2.*step)


            ! calculate residuals3
            do i =1, N_original_verts*3
                residuals3(:,i) = adjoint_mesh%wake%strips(strip)%panels(p)%d_n_g%get_values(i) - d_n_g_FD(:,i)
            end do

            
            if (maxval(abs(residuals3(:,:)))>error_allowed) then
                write(*,*) ""
                write(*,*) "     FLAGGED VALUES :"
                do i = 1, N_original_verts*3
                    if (any(abs(residuals3(:,i))>error_allowed)) then
                        write(*,*) ""
                        write(*,*) "                          d_n_g      x, y, and z "
                        write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_n_g_FD(:,i)
                    
                        write(*, '(A25,8x,3(f25.10, 4x))') "               adjoint",   &
                        adjoint_mesh%wake%strips(strip)%panels(index)%d_n_g%get_values(i)
                        write(*, '(A25,8x,3(f25.10, 4x))') "             residuals", residuals3(:,i)
                        
                    end if
                end do
            end if


            ! check if test failed
            do i=1,N_original_verts*3
                ! if (any(abs(residuals3(:,i)) > error_allowed) .and. any(abs(d_n_g_FD(:,i))<10.0) .and. &
                ! any(abs(residuals3(:,i)) > error_allowed*10.0)) then
                !     test_failed = .true.
                !     exit
                ! else 
                !     test_failed = .false.
                ! end if
                if (any(abs(residuals3(:,i)) > error_allowed)) then 
                    do j = 1,3
                        if (abs(d_n_g_FD(j,i))>100.0) then
                            if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_n_g_FD(j,i)).and. abs(d_n_g_FD(j,i))>10.0) then
                            if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_n_g_FD(j,i)).and. abs(d_n_g_FD(j,i))>1.0) then
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
                write(*,*)"                              d_n_g panel ",y," test FAILED"
                failure_log(total_tests-passed_tests) = "d_n_g test FAILED"
            else
                ! write(*,*) "        CALC d_n_g test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_centr !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! calculate residuals3
            do i =1, N_original_verts*3
                residuals3(:,i) = adjoint_mesh%wake%strips(strip)%panels(p)%d_centr%get_values(i) - d_centr_FD(:,i)
            end do

            
            if (maxval(abs(residuals3(:,:)))>error_allowed) then
                write(*,*) ""
                write(*,*) "     FLAGGED VALUES :"
                do i = 1, N_original_verts*3
                    if (any(abs(residuals3(:,i))>error_allowed)) then
                        write(*,*) ""
                        write(*,*) "                          d_centr      x, y, and z "
                        write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_centr_FD(:,i)
                    
                        write(*, '(A25,8x,3(f25.10, 4x))') "               adjoint",   &
                        adjoint_mesh%wake%strips(strip)%panels(p)%d_centr%get_values(i)
                        write(*, '(A25,8x,3(f25.10, 4x))') "             residuals", residuals3(:,i)
                    end if
                end do
            end if


            ! check if test failed
            do i=1,N_original_verts*3
                ! if (any(abs(residuals3(:,i)) > error_allowed) .and. any(abs(d_centr_FD(:,i))<10.0) .and. &
                ! any(abs(residuals3(:,i)) > error_allowed*10.0)) then
                !     test_failed = .true.
                !     exit
                ! else 
                !     test_failed = .false.
                ! end if
                if (any(abs(residuals3(:,i)) > error_allowed)) then 
                    do j = 1,3
                        if (abs(d_centr_FD(j,i))>100.0) then
                            if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_centr_FD(j,i)) .and. abs(d_centr_FD(j,i))>10.0) then
                            if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_centr_FD(j,i)) .and. abs(d_centr_FD(j,i))>1.0) then
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
                write(*,*)"                              d_centr panel ",y," test FAILED"
                failure_log(total_tests-passed_tests) = "d_centr test FAILED"
            else
                ! write(*,*) "        CALC d_centr test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


            ! calculate residuals
            do i =1, N_original_verts*3
                residuals(i) = adjoint_mesh%wake%strips(strip)%panels(p)%d_A%get_value(i) - d_area_FD(i)
            end do
            
            if (maxval(abs(residuals))>error_allowed) then
                write(*,*) ""
                write(*,*) "     FLAGGED VALUES :"
                write(*,*) "          d_area FD             adjoint d_area             residual"
                do i = 1, N_original_verts*3
                    if (abs(residuals(i))>error_allowed) then
                        write(*, '(8x,(f25.10, 4x),3x, (f25.10, 4x),3x, (f25.10, 4x))') &
                        d_area_FD(i), adjoint_mesh%wake%strips(strip)%panels(p)%d_A%get_value(i), residuals(i)
                    end if
                end do
            end if
            
            
            ! check if test failed
            do i=1,N_original_verts*3
                if (any(abs(residuals) > error_allowed)) then 
                    if (abs(d_area_FD(i))>1000.0) then
                        if (abs(residuals(i)) > error_allowed*10000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (1000.0>abs(d_area_FD(i)) .and. abs(d_area_FD(i))>100.0) then
                        if (abs(residuals(i)) > error_allowed*1000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (100.0>abs(d_area_FD(i)) .and. abs(d_area_FD(i))>10.0) then
                        if (abs(residuals(i)) > error_allowed*100.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (10.0>abs(d_area_FD(i)) .and. abs(d_area_FD(i))>1.0) then
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
                
                ! if (abs(residuals(i)) > error_allowed .and. abs(d_area_FD(i))<10.0 .and. &
                ! abs(residuals(i)) > error_allowed*10.0) then
                !     test_failed = .true.
                !     exit
                ! else 
                !     test_failed = .false.
                ! end if
            end do
            if (test_failed) then
                total_tests = total_tests + 1
                write(*,'(A,I5,A)')"                              d_area panel ",y," test FAILED"
                failure_log(total_tests-passed_tests) = "d_area test FAILED"
            else
                ! write(*,*) "        d_area test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.


            do k = 1,3
                ! do for each design variable
                do i=1,3
                    do j=1,N_original_verts

                        deallocate(test_mesh%vertices, test_mesh%edges, test_mesh%panels, test_mesh%vertex_ordering, &
                        test_mesh%wake%strips)
                        call test_mesh%init(geom_settings)
                        test_mesh%perturb_point = .true.

                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        ! write(*,*) "this vertex is clone? ", test_mesh%vertices(j)%clone 
                        
                        !!!!!!!!!!! update !!!!!!!!!!!!!
                        ! update panel geometry and calc
                        do m =1,N_panels
                            deallocate(test_mesh%panels(m)%n_hat_g)
                            call test_mesh%panels(m)%calc_derived_geom()
                        end do

                        call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)
                        

                        !!!!!!!!!!!!! end update !!!!!!!!!!!!!!!!
                        
                        !!!!!!!!!!!!!!!!!!!!!!! GET INFO !!!!!!!!!!!!!!!!!!!!!!!!     
                        ! put the x y or z com_verts) = test_mesh%wake%strips(strip)%panels(p)%A
                        n_hat_g_up(:,j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%n_hat_g(:,k)


                        !!!! Perturb Down !!!!

                        deallocate(test_mesh%vertices, test_mesh%edges, test_mesh%panels, test_mesh%vertex_ordering, &
                        test_mesh%wake%strips)
                        call test_mesh%init(geom_settings)
                        test_mesh%perturb_point = .true.

                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) -step
                            
                        !!!!!!!!!!! update !!!!!!!!!!!!!!!
                        ! update panel geometry and calc
                        do m =1,N_panels
                            deallocate(test_mesh%panels(m)%n_hat_g)
                            call test_mesh%panels(m)%calc_derived_geom()
                        end do

                        call test_mesh%init_with_flow(freestream_flow, body_file, wake_file, formulation)
                            
                        !!!!!!!!!!!!!!!! end update !!!!!!!!!!!!!!!!!!!!!

                        ! get info 
                        n_hat_g_dn(:,j + (i-1)*N_original_verts) = test_mesh%wake%strips(strip)%panels(p)%n_hat_g(:,k)

                    end do 
                end do

                d_n_hat_g_FD(:,:) = (n_hat_g_up - n_hat_g_dn)/(2.*step)

                ! calculate residuals3
                do i =1, N_original_verts*3
                    residuals3(:,i) = adjoint_mesh%wake%strips(strip)%panels(p)%d_n_hat_g(k)%get_values(i) - d_n_hat_g_FD(:,i)
                end do


                if (maxval(abs(residuals3(:,:)))>error_allowed) then
                    write(*,*) ""
                    write(*,*) "     FLAGGED VALUES :"
                    do i = 1, N_original_verts*3
                        if (any(abs(residuals3(:,i))>error_allowed)) then
                            write(*,*) ""
                            write(*,'(A,I5,A)') "                d_n_hat_g edge ",k,"     x, y, and z "
                            write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_n_hat_g_FD(:,i)
                        
                            write(*, '(A25,8x,3(f25.10, 4x))') "               adjoint",   &
                            adjoint_mesh%wake%strips(strip)%panels(p)%d_n_hat_g(k)%get_values(i)
                            write(*, '(A25,8x,3(f25.10, 4x))') "             residuals", residuals3(:,i)
                        end if
                    end do
                end if
    
    
                ! check if test failed
                do i=1,N_original_verts*3
                    if (any(abs(residuals3(:,i)) > error_allowed)) then 
                        do j = 1,3
                            if (abs(d_n_hat_g_FD(j,i))>1000.0) then
                                if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                    test_failed = .true.
                                    exit
                                else
                                    test_failed = .false.
                                end if
                            elseif (1000.0>abs(d_n_hat_g_FD(j,i)) .and. abs(d_n_hat_g_FD(j,i))>100.0) then
                                if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                    test_failed = .true.
                                    exit
                                else
                                    test_failed = .false.
                                end if
                            elseif (100.0>abs(d_n_hat_g_FD(j,i)) .and. abs(d_n_hat_g_FD(j,i))>10.0) then
                                if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                    test_failed = .true.
                                    exit
                                else
                                    test_failed = .false.
                                end if
                            elseif (10.0>abs(d_n_hat_g_FD(j,i)) .and. abs(d_n_hat_g_FD(j,i))>1.0) then
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
                    write(*,'(A,I5,A,I5,A)')"                              d_n_hat_g panel ",y," edge ",k," test FAILED"
                    failure_log(total_tests-passed_tests) = "d_n_hat_g test FAILED"
                else
                    ! write(*,*) "        CALC d_n_hat_g test PASSED"
                    ! write(*,*) "" 
                    ! write(*,*) ""
                    passed_tests = passed_tests + 1
                    total_tests = total_tests + 1
                    
                end if

                ! reset test failed for the next loop
                test_failed = .false.
                

            end do


        end do ! end p loop

    end do ! end strip loop




!   !!!!!!!!!!!!!   SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "       SUBsonic Wake PANEL GEOM TEST RESULTS (WAKE APPENDED)"
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
    write(*,'(A,f18.10, A)') " Total test time = ", time, " minutes"
    write(*,*) ""
    write(*,*) "----------------------"
    write(*,*) "Program Complete"
    write(*,*) "----------------------"

end program wake_appended_test5_2