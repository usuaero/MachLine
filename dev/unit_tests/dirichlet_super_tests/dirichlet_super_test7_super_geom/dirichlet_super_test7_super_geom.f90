program dirichlet_super_test7
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
    type(dod) :: dod_info

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals

    real,dimension(:,:),allocatable ::  residuals3, l1_up, l1_dn, d_l1_FD, l2_up, l2_dn, d_l2_FD,&
    a_up, a_dn, d_a_FD,g2_up, g2_dn, d_g2_FD, R1_up, R1_dn, d_R1_FD, R2_up, R2_dn, d_R2_FD,&
    dR_up, dR_dn, d_dR_FD

    integer :: i,j,k,m,n,y,z, N_verts, N_panels, vert, index, cp_ind
    real :: step,error_allowed
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

    ! calc CALC BASIC GEOM geom of relation between cp1 and panel1 
    test_geom = test_mesh%panels(1)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
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
    call adjoint_mesh%init_with_flow(adjoint_freestream_flow, adjoint_body_file, adjoint_wake_file, &
                                                                                    adjoint_formulation)

    ! Initialize panel solver
    call adjoint_solver%init(adjoint_solver_settings, adjoint_processing_settings, adjoint_mesh, &
    adjoint_freestream_flow, adjoint_control_point_file)

    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!
    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    allocate(l1_up(3,N_verts*3))
    allocate(l1_dn(3,N_verts*3))
    allocate(d_l1_FD(3,N_verts*3))
    allocate(l2_up(3,N_verts*3))
    allocate(l2_dn(3,N_verts*3))
    allocate(d_l2_FD(3,N_verts*3))
    allocate(a_up(3,N_verts*3))
    allocate(a_dn(3,N_verts*3))
    allocate(d_a_FD(3,N_verts*3))
    allocate(g2_up(3,N_verts*3))
    allocate(g2_dn(3,N_verts*3))
    allocate(d_g2_FD(3,N_verts*3))
    allocate(R1_up(3,N_verts*3))
    allocate(R1_dn(3,N_verts*3))
    allocate(d_R1_FD(3,N_verts*3))
    allocate(R2_up(3,N_verts*3))
    allocate(R2_dn(3,N_verts*3))
    allocate(d_R2_FD(3,N_verts*3))
    allocate(dR_up(3,N_verts*3))
    allocate(dR_dn(3,N_verts*3))
    allocate(d_dR_FD(3,N_verts*3))


    error_allowed = 1.0e-9
    step = 0.000001
    index = 1
    cp_ind = 1
    
    mirror_panel = .false.

    write(*,*) ""
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) "       Dirichlet SUPERSONIC SUBINC GEOMETRY SENSITIVITIES TEST                    "
    write(*,*) "------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""


    
    
    do y =1,N_panels
        index = y
        
        do z = 1,N_verts
            cp_ind = z
            
            dod_info = test_mesh%panels(index)%check_dod(test_mesh%cp(cp_ind)%loc, freestream_flow, mirror_panel)
            if (dod_info%in_dod .and. test_mesh%panels(index)%A > 0.) then

                ! panel is in domain of dependence of point    
                write(*,'(A,I5,A,I5)') "Supersonic Subinc Geom test Panel ",y," cp ",z

                ! Calculate geometric parameters
                if (freestream_flow%supersonic) then
                    if ((mirror_panel .and. test_mesh%panels(index)%r_mir < 0.) .or. &
                                    (.not. mirror_panel .and. test_mesh%panels(index)%r < 0.)) then
                        !geom = this%calc_supersonic_supinc_geom(P, freestream, mirror_panel, dod_info)
                        write(*,*) " can't do super inclined yet"
                        stop
                    else
                        ! supersonic subinclined
                        adjoint_geom = adjoint_mesh%panels(index)%calc_supersonic_subinc_geom_adjoint(&
                                                    adjoint_mesh%cp(cp_ind)%loc,adjoint_mesh%cp(cp_ind)%d_loc,&
                                                     adjoint_freestream_flow, mirror_panel, dod_info)
                        
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

                                ! update supersonic geom
                                test_geom = test_mesh%panels(index)%calc_supersonic_subinc_geom(&
                                                test_mesh%cp(cp_ind)%loc,freestream_flow,mirror_panel,dod_info)
                    
                                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                                
                                ! get desired info
                                l1_up(:,j + (i-1)*N_verts) = test_geom%l1(:)
                                l2_up(:,j + (i-1)*N_verts) = test_geom%l2(:)
                                a_up(:,j + (i-1)*N_verts) = test_geom%a(:)
                                g2_up(:,j + (i-1)*N_verts) = test_geom%g2(:)
                                R1_up(:,j + (i-1)*N_verts) = test_geom%R1(:)
                                R2_up(:,j + (i-1)*N_verts) = test_geom%R2(:)
                                dR_up(:,j + (i-1)*N_verts) = test_geom%dR(:)
                                
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
                                
                                ! update supersonic geom
                                test_geom = test_mesh%panels(index)%calc_supersonic_subinc_geom(&
                                                test_mesh%cp(cp_ind)%loc,freestream_flow,mirror_panel,dod_info)
                                
                                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                                
                                ! get desired info
                                l1_dn(:,j + (i-1)*N_verts) = test_geom%l1(:)
                                l2_dn(:,j + (i-1)*N_verts) = test_geom%l2(:)
                                a_dn(:,j + (i-1)*N_verts) = test_geom%a(:)
                                g2_dn(:,j + (i-1)*N_verts) = test_geom%g2(:)
                                R1_dn(:,j + (i-1)*N_verts) = test_geom%R1(:)
                                R2_dn(:,j + (i-1)*N_verts) = test_geom%R2(:)
                                dR_dn(:,j + (i-1)*N_verts) = test_geom%dR(:)
                                
                                
                                ! restore geometry
                                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                            end do 
                        end do 
                        
                        ! central difference 
                        d_l1_FD = (l1_up - l1_dn)/(2.*step)
                        d_l2_FD = (l2_up - l2_dn)/(2.*step)
                        d_a_FD = (a_up - a_dn)/(2.*step)
                        d_g2_FD = (g2_up - g2_dn)/(2.*step)
                        d_R1_FD = (R1_up - R1_dn)/(2.*step)
                        d_R2_FD = (R2_up - R2_dn)/(2.*step)
                        d_dR_FD = (dR_up - dR_dn)/(2.*step)
                        
                        
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! d_l1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_l1(1)%get_value(i),&
                                                adjoint_geom%d_l1(2)%get_value(i),&
                                                adjoint_geom%d_l1(3)%get_value(i)/) - d_l1_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_l1   values           & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_l1_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_l1(1)%get_value(i),&
                                    adjoint_geom%d_l1(2)%get_value(i),&
                                    adjoint_geom%d_l1(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_l1_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_l1_FD(j,i)) .and. abs(d_l1_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_l1_FD(j,i)) .and. abs(d_l1_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_l1_FD(j,i)) .and. abs(d_l1_FD(j,i))>1.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    ! elseif (abs(d_l1_FD(j,i))<0.00001) then
                                    !     if (abs(residuals3(j,i)) < error_allowed/0.000001) then
                                    !         test_failed = .true.
                                    !         exit
                                    !     else
                                    !         test_failed = .false.
                                    !     end if
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
                                            d_l1 vales panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_l1 test FAILED"
                        else
                            ! write(*,*) "        CALC d_l1 test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.

                        
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_l2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    

                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_l2(1)%get_value(i),&
                                                adjoint_geom%d_l2(2)%get_value(i),&
                                                adjoint_geom%d_l2(3)%get_value(i)/) - d_l2_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_l2  values            & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_l2_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_l2(1)%get_value(i),&
                                    adjoint_geom%d_l2(2)%get_value(i),&
                                    adjoint_geom%d_l2(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_l2_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_l2_FD(j,i)) .and. abs(d_l2_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_l2_FD(j,i)) .and. abs(d_l2_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_l2_FD(j,i)) .and. abs(d_l2_FD(j,i))>1.0) then
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
                                            d_l2 vales panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_l2 test FAILED"
                        else
                            ! write(*,*) "        CALC d_l2 test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.


                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_a !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    

                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_a(1)%get_value(i),&
                                                adjoint_geom%d_a(2)%get_value(i),&
                                                adjoint_geom%d_a(3)%get_value(i)/) - d_a_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_a values             & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_a_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_a(1)%get_value(i),&
                                    adjoint_geom%d_a(2)%get_value(i),&
                                    adjoint_geom%d_a(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_a_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_a_FD(j,i)) .and. abs(d_a_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_a_FD(j,i)) .and. abs(d_a_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_a_FD(j,i)) .and. abs(d_a_FD(j,i))>1.0) then
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
                                            d_a values panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_a test FAILED"
                        else
                            ! write(*,*) "        CALC d_a test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.



                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_g2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_g2(1)%get_value(i),&
                                                adjoint_geom%d_g2(2)%get_value(i),&
                                                adjoint_geom%d_g2(3)%get_value(i)/) - d_g2_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_g2 values             & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_g2_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_g2(1)%get_value(i),&
                                    adjoint_geom%d_g2(2)%get_value(i),&
                                    adjoint_geom%d_g2(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_g2_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_g2_FD(j,i)) .and. abs(d_g2_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_g2_FD(j,i)) .and. abs(d_g2_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_g2_FD(j,i)) .and. abs(d_g2_FD(j,i))>1.0) then
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
                                            d_g2 vales panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_g2 test FAILED"
                        else
                            ! write(*,*) "        CALC d_g2 test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.


                !!!!!!!!!!!!!!!!!!!!!!!!!!!!! d_R1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_R1(1)%get_value(i),&
                                                adjoint_geom%d_R1(2)%get_value(i),&
                                                adjoint_geom%d_R1(3)%get_value(i)/) - d_R1_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_R1 values             & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_R1_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_R1(1)%get_value(i),&
                                    adjoint_geom%d_R1(2)%get_value(i),&
                                    adjoint_geom%d_R1(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_R1_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_R1_FD(j,i)) .and. abs(d_R1_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_R1_FD(j,i)) .and. abs(d_R1_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_R1_FD(j,i)) .and. abs(d_R1_FD(j,i))>1.0) then
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
                                            d_R1 vales panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_R1 test FAILED"
                        else
                            ! write(*,*) "        CALC d_R1 test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_R2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    


                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_R2(1)%get_value(i),&
                                                adjoint_geom%d_R2(2)%get_value(i),&
                                                adjoint_geom%d_R2(3)%get_value(i)/) - d_R2_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_R2 values             & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_R2_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_R2(1)%get_value(i),&
                                    adjoint_geom%d_R2(2)%get_value(i),&
                                    adjoint_geom%d_R2(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_R2_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_R2_FD(j,i)) .and. abs(d_R2_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_R2_FD(j,i)) .and. abs(d_R2_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_R2_FD(j,i)) .and. abs(d_R2_FD(j,i))>1.0) then
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
                                            d_R2 vales panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_R2 test FAILED"
                        else
                            ! write(*,*) "        CALC d_R2 test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.





                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  d_dR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                        ! calculate residuals3
                        do i =1, N_verts*3
                            residuals3(:,i) = (/adjoint_geom%d_dR(1)%get_value(i),&
                                                adjoint_geom%d_dR(2)%get_value(i),&
                                                adjoint_geom%d_dR(3)%get_value(i)/) - d_dR_FD(:,i)
                        end do

                        if (maxval(abs(residuals3(:,:)))>error_allowed) then
                            write(*,*) ""
                            write(*,*) "     FLAGGED VALUES :"
                            do i = 1, N_verts*3
                                if (any(abs(residuals3(:,i))>error_allowed)) then
                                    write(*,*) ""
                                    write(*,*) "                                  d_dR values             & 
                                                                            residuals"
                                    write(*, '(A25,8x,3(f25.10, 4x))') "    Central Difference", d_dR_FD(:,i)
                                
                                    write(*, '(A25,8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') "          adjoint",   &
                                    adjoint_geom%d_dR(1)%get_value(i),&
                                    adjoint_geom%d_dR(2)%get_value(i),&
                                    adjoint_geom%d_dR(3)%get_value(i), residuals3(:,i)
                                end if
                            end do
                        end if
                
                        
                        
                        ! check if test failed
                        do i=1,N_verts*3
                            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                                do j = 1,3
                                    if (abs(d_dR_FD(j,i))>1000.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (1000.0>abs(d_dR_FD(j,i)) .and. abs(d_dR_FD(j,i))>100.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (100.0>abs(d_dR_FD(j,i)) .and. abs(d_dR_FD(j,i))>10.0) then
                                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                            test_failed = .true.
                                            exit
                                        else
                                            test_failed = .false.
                                        end if
                                    elseif (10.0>abs(d_dR_FD(j,i)) .and. abs(d_dR_FD(j,i))>1.0) then
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
                                            d_dR vales panel ",y," cp ",z," test FAILED"
                            failure_log(total_tests-passed_tests) = "d_dR test FAILED"
                        else
                            ! write(*,*) "        CALC d_dR test PASSED"
                            ! write(*,*) "" 
                            ! write(*,*) ""
                            passed_tests = passed_tests + 1
                            total_tests = total_tests + 1
                            
                        end if
                        test_failed = .false.


                    end if
                else
                    !geom = this%calc_subsonic_geom(P, freestream, mirror_panel)
                    write(*,*) " something went wrong, shouldn't be at calc subsonic geom"
                end if
            else    
                write(*,'(A,I5,A,I5)') "Panel ",y," not in dod of cp ",z
            end if
            
    
    
        ! z loop
        end do
    ! y loop
    end do



!!!!!!!!!!!!!! SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "   Dirichlet SUPERSONIC SUBINC GEOM SENSITIVITIES TEST RESULTS "
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

end program dirichlet_super_test7