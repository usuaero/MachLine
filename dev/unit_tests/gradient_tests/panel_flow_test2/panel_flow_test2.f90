program panel_flow2

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
    integer :: start_count, end_count, i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, A_g_to_ls_up, A_g_to_ls_dn,  A_ls_to_g_up, &
    A_ls_to_g_dn, vertices_ls_up, vertices_ls_dn, d_vertices_ls_FD, n_hat_ls_up, n_hat_ls_dn, &
    d_n_hat_ls_FD, T_mu_up, T_mu_dn

    real,dimension(:,:),allocatable :: v, vertex_locs, residuals3 

    real,dimension(:,:,:),allocatable ::  d_A_g_to_ls_FD, d_A_ls_to_g_FD, d_T_mu_FD

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
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))
    
    step = 0.00001
    index = 1
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FLOW-DEPENDENT PANEL SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ FLOW-DEPENDENT PANEL SENSITIVITIES TEST ---------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_A_g_to_ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "-------------------------------- TEST d_A_g_to_ls -------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the global to ls transformation matrix of panel 1 WRT each design variable"
    write(*,*) ""
    
    ! allocate central difference variables
    allocate(A_g_to_ls_up(N_verts*3))
    allocate(A_g_to_ls_dn(N_verts*3))
    allocate(d_A_g_to_ls_FD(3,N_verts*3,3))
    
    

    
    ! do for each row of A_g_to_ls
    do m=1,3
        
        !!!!!!!!! Finite Difference  d_A_g_to_ls (row m) !!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_A_g_to_ls (row ", m, ")"
    
        
        ! for each x, y, z of A_g_to_ls (row m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    call test_mesh%panels(index)%calc_derived_geom()
                    
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_g_to_ls_up(j + (i-1)*N_verts) = test_mesh%panels(index)%A_g_to_ls(m,k)
                    
                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
                    
                    
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    call test_mesh%panels(index)%calc_derived_geom()
                    
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_g_to_ls_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%A_g_to_ls(m,k)
                    
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    
                    
                end do 
            end do 
            
            ! central difference 
            d_A_g_to_ls_FD(m,:,k) = (A_g_to_ls_up - A_g_to_ls_dn)/(2.*step)

        end do

        ! write results
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_g_to_ls_FD panel 1 (row ", m, ")"
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y            d_A_g_to_ls_z "
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x))') d_A_g_to_ls_FD(m,i, :)
        end do 

        
        !!!!!!!!!! ADJOINT d_A_g_to_ls (edge m)!!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  ADJOINT d_A_g_to_ls (row ", m, ")"
        write(*,*) ""
           
        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_g_to_ls panel 1 (row ", m, ")"
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y           d_A_g_to_ls_z             sparse_index       full_index"
        do i=1,adjoint_mesh%panels(index)%d_A_g_to_ls(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_A_g_to_ls(m)%columns(i)%vector_values(:), &
            i, adjoint_mesh%panels(index)%d_A_g_to_ls(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%panels(index)%d_A_g_to_ls(m)%get_values(i) - d_A_g_to_ls_FD(m,i,:)
        end do
        
        write(*,'(A, I1, A)') "         d_A_g_to_ls panel 1 (row ", m, ") expanded "
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y           d_A_g_to_ls_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%panels(index)%d_A_g_to_ls(m)%get_values(i), residuals3(:,i)
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
            write(m_char,'(I1)') m
            failure_log(total_tests-passed_tests) = "d_A_g_to_ls (row "// trim(m_char) // ") test FAILED"
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1, A)') "d_A_g_to_ls (row ",m,") test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""

        
    ! end panel edge loop
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_A_ls_to_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_A_ls_to_g ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the global to ls transformation matrix for panel 1 WRT each design variable"
    write(*,*) ""
    
    ! allocate central difference variables
    allocate(A_ls_to_g_up(N_verts*3))
    allocate(A_ls_to_g_dn(N_verts*3))
    allocate(d_A_ls_to_g_FD(3,N_verts*3,3))

    


    ! do for each row of A_ls_to_g
    do m=1,3

        !!!!!!!!! Finite Difference  d_A_ls_to_g (row m) !!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_A_ls_to_g (row ", m, ")"


        ! for each x, y, z of A_ls_to_g (row m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    call test_mesh%panels(index)%calc_derived_geom()
                    
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    

                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_ls_to_g_up(j + (i-1)*N_verts) = test_mesh%panels(index)%A_ls_to_g(m,k)

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    call test_mesh%panels(index)%calc_derived_geom()
                    
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_ls_to_g_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%A_ls_to_g(m,k)
                
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    

                end do 
            end do 

            ! central difference 
            d_A_ls_to_g_FD(m,:,k) = (A_ls_to_g_up - A_ls_to_g_dn)/(2.*step)

        end do

        ! write results
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_ls_to_g_FD panel 1 (row ", m, ")"
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y            d_A_ls_to_g_z "
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x))') d_A_ls_to_g_FD(m,i, :)
        end do 


        !!!!!!!!!! ADJOINT d_A_ls_to_g (edge m)!!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  ADJOINT d_A_ls_to_g (row ", m, ")"
        write(*,*) ""


        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_ls_to_g panel 1 (row ", m, ")"
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y           d_A_ls_to_g_z             sparse_index       full_index"
        do i=1,adjoint_mesh%panels(index)%d_A_ls_to_g(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_A_ls_to_g(m)%columns(i)%vector_values(:), &
            i, adjoint_mesh%panels(index)%d_A_ls_to_g(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%panels(index)%d_A_ls_to_g(m)%get_values(i) - d_A_ls_to_g_FD(m,i,:)
        end do

        write(*,'(A, I1, A)') "         d_A_ls_to_g panel 1 (row ", m, ") expanded "
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y           d_A_ls_to_g_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%panels(index)%d_A_ls_to_g(m)%get_values(i), residuals3(:,i)
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
            write(m_char,'(I1)') m
            failure_log(total_tests-passed_tests) = "d_A_ls_to_g (row "// trim(m_char) // ") test FAILED"
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1, A)') "d_A_ls_to_g (row ",m,") test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""

        
    ! end panel edge loop
    end do




!!!!!d!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_vertices_ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_vertices_ls ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the panel 1 vertices in local scaled coordinates WRT each design variable"
    write(*,*) ""
    
    ! allocate central difference variables
    allocate(vertices_ls_up(N_verts*3))
    allocate(vertices_ls_dn(N_verts*3))
    allocate(d_vertices_ls_FD(N_verts*3))

    


    ! do for vertex
    do m=1,3

        write(*,*) "---------------------------------- TEST d_vertices_ls VERTEX", m, "  ----------------------------------"
        write(*,*) ""

        ! do for xi and eta
        do n = 1,2


            !!!!!!!!! Finite Difference  d_vertices_ls (vertex m, coordinate n) !!!!!!!!!
            write(*,*) ""
            write(*,*) "---------------------------------------------------------------------------------------------"
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_vertices_ls (vertex ", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_vertices_ls (vertex ", m, " eta coordinate)"
            end if
            

            ! ! for each x, y, z of n_hat_g (edge m) 
            ! do k=1,3
                ! do for each design variable
                do i=1,3
                    do j=1,N_verts
                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        
                        ! update panel geometry and calcs
                        deallocate(test_mesh%panels(index)%n_hat_g)
                        call test_mesh%panels(index)%calc_derived_geom()
                        
                        deallocate(test_mesh%panels(index)%vertices_ls)
                        deallocate(test_mesh%panels(index)%n_hat_ls)
                        deallocate(test_mesh%panels(index)%b)
                        deallocate(test_mesh%panels(index)%b_mir)  
                        deallocate(test_mesh%panels(index)%sqrt_b)
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                        
                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        vertices_ls_up(j + (i-1)*N_verts) = test_mesh%panels(index)%vertices_ls(n,m)

                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    
                        ! update panel geometry and calcs
                        deallocate(test_mesh%panels(index)%n_hat_g)
                        call test_mesh%panels(index)%calc_derived_geom()
                        
                        deallocate(test_mesh%panels(index)%vertices_ls)
                        deallocate(test_mesh%panels(index)%n_hat_ls)
                        deallocate(test_mesh%panels(index)%b)
                        deallocate(test_mesh%panels(index)%b_mir)  
                        deallocate(test_mesh%panels(index)%sqrt_b)
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                        
                        
                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        vertices_ls_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%vertices_ls(n,m)
                    
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        
                        

                    end do 
                end do 

                ! central difference 
                d_vertices_ls_FD(:) = (vertices_ls_up - vertices_ls_dn)/(2.*step)

            ! end do

            ! write results
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls_FD panel 1 (vertex ", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls_FD panel 1 (vertex ", m, " eta coordinate)"
            end if
            do i = 1, N_verts*3
                write(*, '(f14.10)') d_vertices_ls_FD(i)
            end do 
        


            !!!!!!!!!! ADJOINT d_vertices_ls (edge m)!!!!!!!!!!!!!
            write(*,*) ""
            write(*,*) "------------------------------------------------"
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  ADJOINT d_vertices_ls (vertex ", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  ADJOINT d_vertices_ls (vertex ", m, " eta coordinate)"
            end if
            
            write(*,*) ""

            
            ! write sparse matrix
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex ", m, " eta coordinate)"
            end if
            write(*,*) "  sparse value                  sparse_index       full_index"
            do i=1,adjoint_mesh%panels(index)%d_vertices_ls(n,m)%sparse_size
                write(*,'(f14.10, 20x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_vertices_ls(n,m)%elements(i)%value, &
                i, adjoint_mesh%panels(index)%d_vertices_ls(n,m)%elements(i)%full_index
            end do
            write(*,*) ""

            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = adjoint_mesh%panels(index)%d_vertices_ls(n,m)%get_value(i) - d_vertices_ls_FD(i)
            end do

            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex", m, " xi coordinate) expanded"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex ", m, " eta coordinate) expanded"
            end if
            write(*,*) "  adjoint value         residuals"
            do i = 1, N_verts*3
                write(*, '(f14.10,3x, f14.10)') adjoint_mesh%panels(index)%d_vertices_ls(n,m)%get_value(i), residuals(i)
            end do
            write(*,*) ""


            ! check if test failed
            do i=1,N_verts*3
                if (abs(residuals(i)) > 1.0e-8) then
                    test_failed = .true.
                    exit
                else 
                    test_failed = .false.
                end if
            end do
            if (test_failed) then
                total_tests = total_tests + 1
                write(m_char,'(I1)') m

                if (n==1) then
                    failure_log(total_tests-passed_tests) = "d_vertices_ls (vertex "// trim(m_char) // ", xi coordinate) &
                    test FAILED"
                else 
                    failure_log(total_tests-passed_tests) = "d_vertices_ls (vertex "// trim(m_char) // " eta coordinate) &
                    test FAILED"
                end if

                write(*,*) failure_log(total_tests-passed_tests)
            else
                if (n==1) then
                    write(*,'(A, I1, A)') "d_vertices_ls (vertex ",m," xi coordinate) test PASSED"
                else
                    write(*,'(A, I1, A)') "d_vertices_ls (vertex ",m," eta coordinate) test PASSED" 
                end if
               
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
            end if
            test_failed = .false.
            write(*,*) "" 
            write(*,*) ""
        
        ! end "n" loop    
        end do

    ! end "m"" loop
    end do





!!!!!d!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_n_hat_ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_n_hat_ls ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the edge vectors (ls) of panel 1 WRT each design variable"
    write(*,*) ""
    
    ! allocate central difference variables
    allocate(n_hat_ls_up(N_verts*3))
    allocate(n_hat_ls_dn(N_verts*3))
    allocate(d_n_hat_ls_FD(N_verts*3))

    


    ! do for vertex
    do m=1,3

        write(*,*) "---------------------------------- TEST d_n_hat_ls edge", m, "  ----------------------------------"
        write(*,*) ""

        ! do for xi and eta
        do n = 1,2


            !!!!!!!!! Finite Difference  d_n_hat_ls (edge m, coordinate n) !!!!!!!!!
            write(*,*) ""
            write(*,*) "---------------------------------------------------------------------------------------------"
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_n_hat_ls (edge ", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_n_hat_ls (edge ", m, " eta coordinate)"
            end if
           

            ! ! for each x, y, z of n_hat_ls (edge m) 
            ! do k=1,3
                ! do for each design variable
                do i=1,3
                    do j=1,N_verts
                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        
                        ! update panel geometry and calcs
                        deallocate(test_mesh%panels(index)%n_hat_g)
                        call test_mesh%panels(index)%calc_derived_geom()
                        
                        deallocate(test_mesh%panels(index)%vertices_ls)
                        deallocate(test_mesh%panels(index)%n_hat_ls)
                        deallocate(test_mesh%panels(index)%b)
                        deallocate(test_mesh%panels(index)%b_mir)  
                        deallocate(test_mesh%panels(index)%sqrt_b)
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    

                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        n_hat_ls_up(j + (i-1)*N_verts) = test_mesh%panels(index)%n_hat_ls(n,m)

                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    
                        ! update panel geometry and calcs
                        deallocate(test_mesh%panels(index)%n_hat_g)
                        call test_mesh%panels(index)%calc_derived_geom()
                        
                        deallocate(test_mesh%panels(index)%vertices_ls)
                        deallocate(test_mesh%panels(index)%n_hat_ls)
                        deallocate(test_mesh%panels(index)%b)
                        deallocate(test_mesh%panels(index)%b_mir)  
                        deallocate(test_mesh%panels(index)%sqrt_b)
                        call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    
                        
                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        n_hat_ls_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%n_hat_ls(n,m)
                    
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        
                        

                    end do 
                end do 

                ! central difference 
                d_n_hat_ls_FD(:) = (n_hat_ls_up - n_hat_ls_dn)/(2.*step)

            ! end do

            ! write results
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls_FD panel 1 (edge ", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls_FD panel 1 (edge ", m, " eta coordinate)"
            end if
            do i = 1, N_verts*3
                write(*, '(f14.10)') d_n_hat_ls_FD(i)
            end do 
        


            !!!!!!!!!! ADJOINT d_n_hat_ls (edge m)!!!!!!!!!!!!!
            write(*,*) ""
            write(*,*) "------------------------------------------------"
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  ADJOINT d_n_hat_ls (edge ", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  ADJOINT d_n_hat_ls (edge ", m, " eta coordinate)"
            end if
            
            write(*,*) ""

            
            ! write sparse matrix
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge ", m, " eta coordinate)"
            end if
            write(*,*) "  sparse value              sparse_index       full_index"
            do i=1,adjoint_mesh%panels(index)%d_n_hat_ls(n,m)%sparse_size
                write(*,'(f14.10, 20x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_n_hat_ls(n,m)%elements(i)%value, &
                i, adjoint_mesh%panels(index)%d_n_hat_ls(n,m)%elements(i)%full_index
            end do
            write(*,*) ""

            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = adjoint_mesh%panels(index)%d_n_hat_ls(n,m)%get_value(i) - d_n_hat_ls_FD(i)
            end do

            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge", m, " xi coordinate) expanded"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge ", m, " eta coordinate) expanded"
            end if
            write(*,*) "  adjoint value         residuals"
            do i = 1, N_verts*3
                write(*, '(f14.10,3x, f14.10)') adjoint_mesh%panels(index)%d_n_hat_ls(n,m)%get_value(i), residuals(i)
            end do
            write(*,*) ""


            ! check if test failed
            do i=1,N_verts*3
                if (abs(residuals(i)) > 1.0e-8) then
                    test_failed = .true.
                    exit
                else 
                    test_failed = .false.
                end if
            end do
            if (test_failed) then
                total_tests = total_tests + 1
                write(m_char,'(I1)') m

                if (n==1) then
                    failure_log(total_tests-passed_tests) = "d_n_hat_ls (edge "// trim(m_char) // ", xi coordinate) &
                    test FAILED"
                else 
                    failure_log(total_tests-passed_tests) = "d_n_hat_ls (edge "// trim(m_char) // " eta coordinate) &
                    test FAILED"
                end if

                write(*,*) failure_log(total_tests-passed_tests)
            else
                if (n==1) then
                    write(*,'(A, I1, A)') "d_n_hat_ls (edge ",m," xi coordinate) test PASSED"
                else
                    write(*,'(A, I1, A)') "d_n_hat_ls (edge ",m," eta coordinate) test PASSED" 
                end if
               
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
            end if
            test_failed = .false.
            write(*,*) "" 
            write(*,*) ""
        
        ! end "n" loop    
        end do

    ! end "m"" loopn_hat
    end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_T_mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_T_mu ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the global to ls transformation matrix for panel 1 WRT each design variable"
    write(*,*) ""
    
    ! allocate central difference variables
    allocate(T_mu_up(N_verts*3))
    allocate(T_mu_dn(N_verts*3))
    allocate(d_T_mu_FD(3,N_verts*3,3))

    


    ! do for each row of T_mu
    do m=1,3

        !!!!!!!!! Finite Difference  d_T_mu (row m) !!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_T_mu (row ", m, ")"


        ! for each x, y, z of A_ls_to_g (row m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%vertices)
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    deallocate(test_mesh%panels(index)%edge_is_discontinuous)
                    call test_mesh%panels(index)%init(test_mesh%vertices(1),test_mesh%vertices(2),test_mesh%vertices(3),index)
         
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    deallocate(test_mesh%panels(index)%i_vert_d)
                    deallocate(test_mesh%panels(index)%S_mu_inv)
                    deallocate(test_mesh%panels(index)%T_mu)
                    deallocate(test_mesh%panels(index)%i_panel_s)
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    call test_mesh%panels(index)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                    test_mesh%vertices,.false.)
                    

                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    T_mu_up(j + (i-1)*N_verts) = test_mesh%panels(index)%S_mu_inv(m,k)

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%vertices)
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    deallocate(test_mesh%panels(index)%edge_is_discontinuous)
                    call test_mesh%panels(index)%init(test_mesh%vertices(1),test_mesh%vertices(2),test_mesh%vertices(3),index)
                    
                    deallocate(test_mesh%panels(index)%vertices_ls)
                    deallocate(test_mesh%panels(index)%n_hat_ls)
                    deallocate(test_mesh%panels(index)%b)
                    deallocate(test_mesh%panels(index)%b_mir)  
                    deallocate(test_mesh%panels(index)%sqrt_b)
                    deallocate(test_mesh%panels(index)%i_vert_d)
                    deallocate(test_mesh%panels(index)%S_mu_inv)
                    deallocate(test_mesh%panels(index)%T_mu)
                    deallocate(test_mesh%panels(index)%i_panel_s)
                    call test_mesh%panels(index)%init_with_flow(freestream_flow, .false., 0)
                    call test_mesh%panels(index)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                    test_mesh%vertices,.false.)
                    
                    
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    T_mu_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%S_mu_inv(m,k)
                
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    

                end do 
            end do 

            ! central difference 
            d_T_mu_FD(m,:,k) = (T_mu_up - T_mu_dn)/(2.*step)

        end do

        ! write results
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_T_mu_FD panel 1 (row ", m, ")"
        write(*,*) "  d_T_mu_x           d_T_mu_y            d_T_mu_z "
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x))') d_T_mu_FD(m,i, :)
        end do 


        !!!!!!!!!! ADJOINT d_T_mu (edge m)!!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  ADJOINT d_T_mu (row ", m, ")"
        write(*,*) ""


        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_T_mu panel 1 (row ", m, ")"
        write(*,*) "  d_T_mu_x           d_T_mu_y           d_T_mu_z             sparse_index       full_index"
        do i=1,adjoint_mesh%panels(index)%d_T_mu(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_T_mu(m)%columns(i)%vector_values(:), &
            i, adjoint_mesh%panels(index)%d_T_mu(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%panels(index)%d_T_mu(m)%get_values(i) - d_T_mu_FD(m,i,:)
        end do

        write(*,'(A, I1, A)') "         d_T_mu panel 1 (row ", m, ") expanded "
        write(*,*) "  d_T_mu_x           d_T_mu_y           d_T_mu_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%panels(index)%d_T_mu(m)%get_values(i), residuals3(:,i)
        end do
        write(*,*) ""


        ! check if test failed
        do i=1,N_verts*3
            if (any(abs(residuals3(:,i)) > 1.0e-7)) then
                test_failed = .true.
                exit
            else 
                test_failed = .false.
            end if
        end do
        if (test_failed) then
            total_tests = total_tests + 1
            write(m_char,'(I1)') m
            failure_log(total_tests-passed_tests) = "d_T_mu (row "// trim(m_char) // ") test FAILED"
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1, A)') "d_T_mu (row ",m,") test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""

        
    ! end panel edge loop
    end do


    !!!!!!!!!!!!!! FLOW-DEPENDENT PANEL SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------FLOW-DEPENDENT PANEL SENSITIVITIES TEST RESULTS--------------"
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

end program panel_flow2