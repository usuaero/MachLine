program dirichlet_super_test1

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
    integer :: i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(:),allocatable :: residuals, X_beta, loc_up, loc_dn, centr_up, centr_dn, normal_up, normal_dn, &
    area_up, area_dn, d_area_FD, n_hat_g_up, n_hat_g_dn

    real,dimension(:,:),allocatable :: v, vertex_locs, d_loc_FD, d_centr_FD, d_n_g_FD, d_norm_d_g_FD, residuals3 

    real,dimension(:,:,:),allocatable :: d_n_hat_g_FD

    integer :: i,j,k,m,n,y,z, N_verts, N_panels, vert, index, vert_ind
    real :: step, error_allowed
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(500) :: failure_log
    character(len=10) :: m_char
    integer(8) :: start_count, end_count
    real(16) :: count_rate, time



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
    ! call input_json%get('post_processing', processing_settings, found)
    ! call input_json%get('output', output_settings, found)

    ! Initialize surface mesh
    call test_mesh%init(geom_settings)

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
    call input_json%load_file(filename=adjoint_input)
    call json_check()
    call input_json%get('flow', adjoint_flow_settings, found)
    call input_json%get('geometry', adjoint_geom_settings, found)
    call input_json%get('solver', adjoint_solver_settings, found)
    ! call input_json%get('post_processing', processing_settings, found)
    ! call input_json%get('output', output_settings, found)

    call adjoint_mesh%init(adjoint_geom_settings)

    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!

    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels

    allocate(loc_up(N_verts*3))
    allocate(loc_dn(N_verts*3))
    allocate(d_loc_FD(3,N_verts*3))
    allocate(centr_up(N_verts*3))
    allocate(centr_dn(N_verts*3))
    allocate(d_centr_FD(3,N_verts*3))
    allocate(normal_up(N_verts*3))
    allocate(normal_dn(N_verts*3))
    allocate(d_n_g_FD(3,N_verts*3))
    allocate(area_up(N_verts*3))
    allocate(area_dn(N_verts*3))
    allocate(d_area_FD(N_verts*3))
    allocate(n_hat_g_up(N_verts*3))
    allocate(n_hat_g_dn(N_verts*3))
    allocate(d_n_hat_g_FD(3,N_verts*3,3))
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    error_allowed = 1.0e-9
    step = 0.000001
    index = 1
    vert_ind = 1
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL GEOMETRY TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) "--------------------------------------------------------------------------"
    write(*,*) "                    Dirichlet Supersonic PANEL GEOMETRY TEST "
    write(*,*) "--------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""
    
    do z= 1,N_verts ! loop through vertices
        vert_ind = z
        write(*,'(A,I5)') "MESH VERTEX TEST ", z
        
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_loc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! for each x, y, z of each vertex
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts

                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                    ! get desired info
                    loc_up(j + (i-1)*N_verts) = test_mesh%vertices(vert_ind)%loc(k)

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    ! get desired info
                    loc_dn(j + (i-1)*N_verts) = test_mesh%vertices(vert_ind)%loc(k)
                    
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                end do 
            end do 
            
            ! central difference 
            d_loc_FD(k,:) = (loc_up - loc_dn)/(2.*step)
                
        end do
                


        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%vertices(vert_ind)%d_loc%get_values(i) - d_loc_FD(:,i)
        end do

        if (maxval(abs(residuals3(:,:)))>error_allowed) then
            write(*,*) ""
            write(*,*) "     FLAGGED VALUES :"
            do i = 1, N_verts*3
                if (any(abs(residuals3(:,i))>error_allowed)) then
                    write(*,*) "         Central Difference    d_loc_g      x, y, and z"
                    write(*, '(8x,3(f25.10, 4x))') d_loc_FD(:,i)
                    write(*,*) "        adjoint       d_loc_g      x, y, and z                  &
                    residuals"
                    write(*, '(8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') &
                    adjoint_mesh%vertices(vert_ind)%d_loc%get_values(i), residuals3(:,i)
                end if
            end do
        end if


        ! check if test failed
        do i=1,N_verts*3
            if (any(abs(residuals3(:,i)) > error_allowed)) then 
                do j = 1,3
                    if (abs(d_loc_FD(j,i))>100.0) then
                        if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (100.0>abs(d_loc_FD(j,i)) .and. abs(d_loc_FD(j,i)) >10.0) then
                        if (abs(residuals3(j,i)) > error_allowed*100.0) then
                            test_failed = .true.
                            exit
                        else
                            test_failed = .false.
                        end if
                    elseif (10.0>abs(d_loc_FD(j,i)) .and. abs(d_loc_FD(j,i)) >1.0) then
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
            ! if (any(abs(residuals3(:,i)) > error_allowed) .and. any(abs(d_loc_FD(:,i))<10.0) .and. &
            ! any(abs(residuals3(:,i)) > error_allowed*10.0)) then
            !     test_failed = .true.
            !     exit
            ! else 
            !     test_failed = .false.
            ! end if
        end do
        if (test_failed) then
            total_tests = total_tests + 1
            write(*,*)"                              d_loc_g vertex ", z," test FAILED"
            failure_log(total_tests-passed_tests) = "d_loc_g test FAILED"
        else
            ! write(*,*) "        CALC d_loc_g test PASSED"
            ! write(*,*) "" 
            ! write(*,*) ""
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.

    end do ! z loop end vertices check

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    FOR EACH PANEL       !!!!!!!!!!!!!!!!!!!!!
    do y = 1,N_panels
        index = y
        write(*,'(A,I5)') "PANEL TEST ", y

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_normal and d_centr !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! for each x y z coord of centroid or normal
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                    !!!!!!!!!!   update  !!!!!!!!!!!
                    
                    ! update panel normal
                    call test_mesh%panels(index)%calc_normal()

                    ! update panel centroid 
                    call test_mesh%panels(index)%calc_centroid()

                    
                    !!!!!!!!!!   end update  !!!!!!!!!!!
                    
                    ! get desired info
                    normal_up(j + (i-1)*N_verts) = test_mesh%panels(index)%n_g(k)
                    centr_up(j + (i-1)*N_verts) = test_mesh%panels(index)%centr(k)

                    ! perturb down the current design variable
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                    !!!!!!!!!!   update  !!!!!!!!!!!

                    ! update panel normal
                    call test_mesh%panels(index)%calc_normal()

                    ! update panel centroid 
                    call test_mesh%panels(index)%calc_centroid()

                    
                    !!!!!!!!!!   end update  !!!!!!!!!!!

                    ! get desired info
                    centr_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%centr(k)
                    normal_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%n_g(k)
                    
                    ! restore geometry
                    test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                end do 
            end do 
            
            ! central difference 
            d_centr_FD(k,:) = (centr_up - centr_dn)/(2.*step)
            d_n_g_FD(k,:) = (normal_up - normal_dn)/(2.*step)
                
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_n_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%panels(index)%d_n_g%get_values(i) - d_n_g_FD(:,i)
        end do

        
        if (maxval(abs(residuals3(:,:)))>error_allowed) then
            write(*,*) ""
            write(*,*) "     FLAGGED VALUES :"
            do i = 1, N_verts*3
                if (any(abs(residuals3(:,i))>error_allowed)) then
                    write(*,*) "         Central Difference    d_n_g      x, y, and z"
                    write(*, '(8x,3(f25.10, 4x))') d_n_g_FD(:,i)
                    write(*,*) "        adjoint       d_n_g      x, y, and z                  &
                    residuals"
                    write(*, '(8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') &
                    adjoint_mesh%panels(index)%d_n_g%get_values(i), residuals3(:,i)
                end if
            end do
        end if


        ! check if test failed
        do i=1,N_verts*3
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
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_mesh%panels(index)%d_centr%get_values(i) - d_centr_FD(:,i)
        end do

        
        if (maxval(abs(residuals3(:,:)))>error_allowed) then
            write(*,*) ""
            write(*,*) "     FLAGGED VALUES :"
            do i = 1, N_verts*3
                if (any(abs(residuals3(:,i))>error_allowed)) then
                    write(*,*) "         Central Difference    d_centr     x, y, and z"
                    write(*, '(8x,3(f25.10, 4x))') d_centr_FD(:,i)
                    write(*,*) "        adjoint       d_centr      x, y, and z                  &
                    residuals"
                    write(*, '(8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') &
                    adjoint_mesh%panels(index)%d_centr%get_values(i), residuals3(:,i)
                end if
            end do
        end if


        ! check if test failed
        do i=1,N_verts*3
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_area !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        ! do for each design variable 
        do i=1,3
            do j=1,N_verts
                ! perturb up the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

                !!!!!!!!!!   update  !!!!!!!!!!!
                ! update panel area calculations
                call test_mesh%panels(index)%calc_area()
                !!!!!!!!!!   end update  !!!!!!!!!!!

                ! get desired info
                area_up(j + (i-1)*N_verts) = test_mesh%panels(index)%A

                ! perturb down the current design variable
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step

                !!!!!!!!!!   update  !!!!!!!!!!!
                ! update panel area calculations
                call test_mesh%panels(index)%calc_area()
                !!!!!!!!!!   end update  !!!!!!!!!!!
                
                ! get desired info
                area_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%A
                        
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step

            end do 
        end do 

        ! central difference 
        d_area_FD = (area_up - area_dn)/(2.*step)

        
        
        ! calculate residuals
        do i =1, N_verts*3
            residuals(i) = adjoint_mesh%panels(index)%d_A%get_value(i) - d_area_FD(i)
        end do
        
        if (maxval(abs(residuals))>error_allowed) then
            write(*,*) ""
            write(*,*) "     FLAGGED VALUES :"
            write(*,*) "          d_area FD             adjoint d_area             residual"
            do i = 1, N_verts*3
                if (abs(residuals(i))>error_allowed) then
                    write(*, '(8x,(f25.10, 4x),3x, (f25.10, 4x),3x, (f25.10, 4x))') &
                    d_area_FD(i), adjoint_mesh%panels(index)%d_A%get_value(i), residuals(i)
                end if
            end do
        end if
        
        
        ! check if test failed
        do i=1,N_verts*3
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


    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_n_hat_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
        ! do for each edge
        do m=1,3

            ! for each x, y, z of n_hat_g (edge m) 
            do k=1,3
                ! do for each design variable
                do i=1,3
                    do j=1,N_verts
                        ! perturb up the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        
                        !!!!!!!!!!   update  !!!!!!!!!!!
                        deallocate(test_mesh%panels(index)%n_hat_g)
                        call test_mesh%panels(index)%calc_derived_geom()
                        !!!!!!!!!!   end update  !!!!!!!!!!!

                        ! get desired info
                        n_hat_g_up(j + (i-1)*N_verts) = test_mesh%panels(index)%n_hat_g(k,m)

                        ! perturb down the current design variable
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
                    
                        !!!!!!!!!!   update  !!!!!!!!!!!
                        deallocate(test_mesh%panels(index)%n_hat_g)
                        call test_mesh%panels(index)%calc_derived_geom()
                        !!!!!!!!!!   end update  !!!!!!!!!!!
                        
                        ! get desired info
                        n_hat_g_dn(j + (i-1)*N_verts) = test_mesh%panels(index)%n_hat_g(k,m)
                    
                        ! restore geometry
                        test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
                        

                    end do 
                end do 

                ! central difference 
                d_n_hat_g_FD(k,:,m) = (n_hat_g_up - n_hat_g_dn)/(2.*step)

            end do


            ! calculate residuals3
            do i =1, N_verts*3
                residuals3(:,i) = adjoint_mesh%panels(index)%d_n_hat_g(m)%get_values(i) - d_n_hat_g_FD(:,i,m)
            end do


            if (maxval(abs(residuals3(:,:)))>error_allowed) then
                write(*,*) ""
                write(*,*) "     FLAGGED VALUES :"
                do i = 1, N_verts*3
                    if (any(abs(residuals3(:,i))>error_allowed)) then
                        write(*,'(A,I5,A)') "         Central Difference    d_n_hat_g edge ",m,"      x, y, and z"
                        write(*, '(8x,3(f25.10, 4x))') d_n_hat_g_FD(:,i,m)
                        write(*,'(A,I5,A)') "        adjoint       d_n_hat_g edge ",m,"      x, y, and z                  &
                        residuals"
                        write(*, '(8x,3(f25.10, 4x),3x, 3(f25.10, 4x))') &
                        adjoint_mesh%panels(index)%d_n_hat_g(m)%get_values(i), residuals3(:,i)
                    end if
                end do
            end if
    
    
            ! check if test failed
            do i=1,N_verts*3
                if (any(abs(residuals3(:,i)) > error_allowed)) then 
                    do j = 1,3
                        if (abs(d_n_hat_g_FD(j,i,m))>1000.0) then
                            if (abs(residuals3(j,i)) > error_allowed*10000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (1000.0>abs(d_n_hat_g_FD(j,i,m)) .and. abs(d_n_hat_g_FD(j,i,m))>100.0) then
                            if (abs(residuals3(j,i)) > error_allowed*1000.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (100.0>abs(d_n_hat_g_FD(j,i,m)) .and. abs(d_n_hat_g_FD(j,i,m))>10.0) then
                            if (abs(residuals3(j,i)) > error_allowed*100.0) then
                                test_failed = .true.
                                exit
                            else
                                test_failed = .false.
                            end if
                        elseif (10.0>abs(d_n_hat_g_FD(j,i,m)) .and. abs(d_n_hat_g_FD(j,i,m))>1.0) then
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
                write(*,'(A,I5,A,I5,A)')"                              d_n_hat_g panel ",y," edge ",m," test FAILED"
                failure_log(total_tests-passed_tests) = "d_n_hat_g test FAILED"
            else
                ! write(*,*) "        CALC d_n_hat_g test PASSED"
                ! write(*,*) "" 
                ! write(*,*) ""
                passed_tests = passed_tests + 1
                total_tests = total_tests + 1
                
            end if
            test_failed = .false.

            

            
        ! end panel edge loop
        end do

    end do ! end panel loop

    !!!!!!!!!!!!!! PANEL GEOMETRY  RESULTS!!!!!!!!!!!!!
        write(*,*) "------------------------------------------------------------------------------"
        write(*,*) "                    Dirichlet Dirichlet SupersonicSupersonic PANEL GEOMETRY TEST RESULTS "
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
        write(*,'(A,f25.10, A)') " Total test time = ", time, " minutes"
        write(*,*) ""
        write(*,*) "----------------------"
        write(*,*) "Program Complete"
        write(*,*) "----------------------"


end program dirichlet_super_test1