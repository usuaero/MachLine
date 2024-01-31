program panel_geom2

    ! tests various intermediate sensitivities 
    use adjoint_mod
    use base_geom_mod
    use panel_mod
    use flow_mod
    use surface_mesh_mod
    use json_mod
    use json_xtnsn_mod
    use panel_solver_mod
    use helpers_mod
    
    implicit none

    !!!!!!!!!!!!!!!!!!! STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!
    character(100) :: test_input, adjoint_input
    character(len=:),allocatable :: body_file, wake_file, control_point_file, points_file, points_output_file
    character(len=:),allocatable :: report_file, spanwise_axis
    character(len=:),allocatable :: formulation

    type(json_file) :: input_json
    type(json_value),pointer :: flow_settings, &
                                geom_settings, &
                                solver_settings, &
                                adjoint_flow_settings, &
                                adjoint_geom_settings, &
                                adjoint_solver_settings
                                
    type(surface_mesh) :: test_mesh, adjoint_mesh
    type(flow) :: freestream_flow
    integer :: start_count, end_count, i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, loc_up, loc_dn, centr_up, centr_dn, normal_up, normal_dn, &
    area_up, area_dn, d_area_FD, n_hat_g_up, n_hat_g_dn

    real,dimension(:,:),allocatable :: v, vertex_locs, d_loc_FD, d_centr_FD, d_n_g_FD, d_norm_d_g_FD, residuals3 

    real,dimension(:,:,:),allocatable :: d_n_hat_g_FD

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
    ! call input_json%get('post_processing', processing_settings, found)
    ! call input_json%get('output', output_settings, found)

    ! Initialize surface mesh
    call test_mesh%init(geom_settings)

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
    call input_json%load_file(filename=adjoint_input)
    call json_check()
    call input_json%get('flow', adjoint_flow_settings, found)
    call input_json%get('geometry', adjoint_geom_settings, found)
    call input_json%get('solver', adjoint_solver_settings, found)
    ! call input_json%get('post_processing', processing_settings, found)
    ! call input_json%get('output', output_settings, found)

    call adjoint_mesh%init(adjoint_geom_settings)

    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!
    index = 1

    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    step = 0.00001
    index = 1
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL GEOMETRY TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "-------------------------------------------- PANEL GEOMETRY TEST ----------------------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_loc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_loc -----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of vertex 1 WRT each design variable"
    write(*,*) ""
    
    ! following sensitivities are with respect to a perturbation in x y and z of vertex 1
    index = 1

    !!!!!!!!! Finite Difference  d_loc !!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  CENTRAL DIFFERENCE d_loc"
    
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

                ! put the x y or z component of the vertex of interest (index) in a list
                loc_up(j + (i-1)*N_verts) = test_mesh%vertices(index)%loc(k)

                ! perturb down the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                ! put the x y or z component of the vertex of interest (index) in a list
                loc_dn(j + (i-1)*N_verts) = test_mesh%vertices(index)%loc(k)
                
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
    
    !!!!!!!!!! ADJOINT d_loc!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  ADJOINT d_loc"
    write(*,*) ""
    
    !write sparse matrix
    write(*,*) ""
    write(*,*) "         d_loc point 1"
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z             sparse_index       full_index"

    do i=1,adjoint_mesh%vertices(index)%d_loc%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%vertices(index)%d_loc%columns(i)%vector_values, &
        i, adjoint_mesh%vertices(index)%d_loc%columns(i)%full_index
    end do
    write(*,*) ""


    ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = adjoint_mesh%vertices(index)%d_loc%get_values(i) - d_loc_FD(:,i)
    end do

    
    write(*,*) "         d_loc vertex 1 expanded "
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z                                 residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%vertices(index)%d_loc%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (any(abs(residuals3(:,i)) > 1.0e-10)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "d_loc test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "d_loc test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
        
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_centr !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_centr ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the CENTROID of panel 1 WRT each design variable"
    write(*,*) ""

    !!!!!!!!! Finite Difference  d_centr !!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  CENTRAL DIFFERENCE d_centr"

    ! perturb x1 up
    allocate(centr_up(N_verts*3))
    allocate(centr_dn(N_verts*3))
    allocate(d_centr_FD(3,N_verts*3))

    ! for each x, y, z of centr 1 
    do k=1,3
        ! do for each design variable
        do i=1,3
            do j=1,N_verts
                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                ! update panel centroid calculations
                call test_mesh%panels(index)%calc_centroid()

                ! put the x y or z component of the panel's perturbed centroid in a list
                centr_up(j + (i-1)*N_verts) = test_mesh%panels(index)%centr(k)

                ! perturb down the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                ! update panel centroid calculationsadjoint_mesh%panel    
                call test_mesh%panels(index)%calc_centroid()
                
                ! put the x y or z component of the panel's perturbed centroid in a list
                centr_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%centr(k)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

            end do 
        end do 
        
        ! central difference 
        d_centr_FD(k,:) = (centr_up - centr_dn)/(2.*step)
            
    end do

    ! write results
    write(*,*) ""
    write(*,*) "                d_centr_FD panel 1"
    write(*,*) "  d_centr_x         d_centr_y         d_centr_z "
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x))') d_centr_FD(:,i)
    end do 
    
    !!!!!!!!!! ADJOINT d_centr!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  ADJOINT d_centr"
    write(*,*) ""


    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_centr panel 1"
    write(*,*) "  d_centr_x         d_centr_y         d_centr_z           sparse_index       full_index"
    do i=1,adjoint_mesh%panels(index)%d_centr%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_centr%columns(i)%vector_values(:), &
        i, adjoint_mesh%panels(index)%d_centr%columns(i)%full_index
    end do
    write(*,*) ""

   ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = adjoint_mesh%panels(index)%d_centr%get_values(i) - d_centr_FD(:,i)
    end do
    

    write(*,*) "         d_centr panel 1 expanded "
    write(*,*) "  d_centr_x         d_centr_y         d_centr_z                               residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%panels(index)%d_centr%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (any(abs(residuals3(:,i)) > 1.0e-10)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "d_centr test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "d_centr test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_n_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_n_g ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the NORMAL VECTOR of panel 1 WRT each design variable"
    write(*,*) ""
    


    !!!!!!!!! Finite Difference  d_n_g !!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  CENTRAL DIFFERENCE d_n_g"

    ! perturb x1 up
    allocate(normal_up(N_verts*3))
    allocate(normal_dn(N_verts*3))
    allocate(d_n_g_FD(3,N_verts*3))

    ! for each x, y, z of normal 1 
    do k=1,3
        ! do for each design variable
        do i=1,3
            do j=1,N_verts
                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                ! update panel normal calculations
                call test_mesh%panels(index)%calc_normal()

                ! put the x y or z component of the panel's perturbed normal in a list
                normal_up(j + (i-1)*N_verts) = test_mesh%panels(index)%n_g(k)

                ! perturb down the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                ! update panel normal calculations
                call test_mesh%panels(index)%calc_normal()
                
                ! put the x y or z component of the panel's perturbed normal in a list
                normal_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%n_g(k)
            
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

            end do 
        end do 

        ! central difference 
        d_n_g_FD(k,:) = (normal_up - normal_dn)/(2.*step)

    end do

    ! write results
    write(*,*) ""
    write(*,*) "                d_n_g_FD panel 1"
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z "
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x))') d_n_g_FD(:,i)
    end do 


    !!!!!!!!!! ADJOINT d_n_g!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  ADJOINT d_n_g"
    write(*,*) ""
            
    

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_n_g panel 1"
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z             sparse_index       full_index"
    do i=1,adjoint_mesh%panels(index)%d_n_g%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_n_g%columns(i)%vector_values(:), &
        i, adjoint_mesh%panels(index)%d_n_g%columns(i)%full_index
    end do
    write(*,*) ""

    ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = adjoint_mesh%panels(index)%d_n_g%get_values(i) - d_n_g_FD(:,i)
    end do

    write(*,*) "         d_n_g panel 1 expanded "
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z                                residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%panels(index)%d_n_g%get_values(i), residuals3(:,i)
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
        failure_log(total_tests-passed_tests) = "d_n_g test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "d_n_g test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_area !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_area ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the area of panel 1 WRT each design variable"
    write(*,*) ""
    
    
    !!!!!!!!! Finite Difference  d_area !!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  CENTRAL DIFFERENCE d_area"


    ! perturb x1 up
    allocate(area_up(N_verts*3))
    allocate(area_dn(N_verts*3))
    allocate(d_area_FD(N_verts*3))

    ! do for each design variable 
    do i=1,3
        do j=1,N_verts
            ! perturb up the current design variable
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

            ! update panel area calculations
            call test_mesh%panels(index)%calc_area()

            ! put the x y or z component of the panel's perturbed area in a list
            area_up(j + (i-1)*N_verts) = test_mesh%panels(index)%A

            ! perturb down the current design variable
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

            ! update panel area calculations
            call test_mesh%panels(index)%calc_area()
            
            ! put the x y or z component of the panel's perturbed area in a list
            area_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%A
                    
            ! restore geometry
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

        end do 
    end do 

    ! central difference 
    d_area_FD = (area_up - area_dn)/(2.*step)

    ! write results
    write(*,*) ""
    write(*,*) "          d_area_FD panel 1"
    write(*,*) "  d_area "
    do i = 1, N_verts*3
        write(*, '(f14.10)') d_area_FD(i)
    end do 


    !!!!!!!!!! ADJOINT d_area!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,*) "  ADJOINT d_area"
    write(*,*) ""
            
    

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_area panel 1"
    write(*,*) "  d_area                  sparse_index       full_index"
    do i=1,adjoint_mesh%panels(index)%d_A%sparse_size
        write(*,'(f14.10, 20x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_A%elements(i)%value, &
        i, adjoint_mesh%panels(index)%d_A%elements(i)%full_index
    end do
    write(*,*) ""

    ! calculate residuals

    do i =1, N_verts*3
        residuals(i) = adjoint_mesh%panels(index)%d_A%get_value(i) - d_area_FD(i)
    end do

    write(*,*) "         d_area panel 1 expanded "
    write(*,*) "  d_area         residuals"
    do i = 1, N_verts*3
        write(*, '(f14.10,3x, f14.10)') adjoint_mesh%panels(index)%d_A%get_value(i), residuals(i)
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
        failure_log(total_tests-passed_tests) = "d_area test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "d_area test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""




!!!!!d!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_n_hat_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_n_hat_g ----------------------------------"
    write(*,*) ""
    write(*,*) "the sensitivity of the edge vectors (global) of panel 1 WRT each design variable"
    write(*,*) ""
    
    ! allocate central difference variables
    allocate(n_hat_g_up(N_verts*3))
    allocate(n_hat_g_dn(N_verts*3))
    allocate(d_n_hat_g_FD(3,N_verts*3,3))


    ! do for each edge
    do m=1,3

        !!!!!!!!! Finite Difference  d_n_hat_g (edge m) !!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_n_hat_g (edge ", m, ")"

        ! for each x, y, z of n_hat_g (edge m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    call test_mesh%panels(index)%calc_derived_geom()

                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    n_hat_g_up(j + (i-1)*N_verts) = test_mesh%panels(index)%n_hat_g(k,m)

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
                
                    ! update panel geometry and calcs
                    deallocate(test_mesh%panels(index)%n_hat_g)
                    call test_mesh%panels(index)%calc_derived_geom()
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    n_hat_g_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%n_hat_g(k,m)
                
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                    

                end do 
            end do 

            ! central difference 
            d_n_hat_g_FD(k,:,m) = (n_hat_g_up - n_hat_g_dn)/(2.*step)

        end do

        ! write results
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_n_hat_g_FD panel 1 (edge ", m, ")"
        write(*,*) "  d_n_hat_g_x           d_n_hat_g_y            d_n_hat_g_z "
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x))') d_n_hat_g_FD(:,i, m)
        end do 


        !!!!!!!!!! ADJOINT d_n_hat_g (edge m)!!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        write(*,'(A, I1, A)') "  ADJOINT d_n_hat_g (edge ", m, ")"
        write(*,*) ""


        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_n_hat_g panel 1 (edge ", m, ")"
        write(*,*) "  d_n_hat_g_x           d_n_hat_g_y           d_n_hat_g_z             sparse_index       full_index"
        do i=1,adjoint_mesh%panels(index)%d_n_hat_g(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_mesh%panels(index)%d_n_hat_g(m)%columns(i)%vector_values(:), &
            i, adjoint_mesh%panels(index)%d_n_hat_g(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%panels(index)%d_n_hat_g(m)%get_values(i) - d_n_hat_g_FD(:,i,m)
        end do

        write(*,'(A, I1, A)') "         d_n_hat_g panel 1 (edge ", m, ") expanded "
        write(*,*) "  d_n_hat_g_x           d_n_hat_g_y           d_n_hat_g_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_mesh%panels(index)%d_n_hat_g(m)%get_values(i), residuals3(:,i)
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
            failure_log(total_tests-passed_tests) = "d_n_hat_g (edge "// trim(m_char) // ") test FAILED"
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A, I1, A)') "d_n_hat_g (edge ",m,") test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""

        
    ! end panel edge loop
    end do



    !!!!!!!!!!!!!! GRADIENT TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------PANEL GEOMETRY SENSITIVITIES TEST RESULTS--------------"
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


end program panel_geom2