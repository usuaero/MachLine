program calc_subsonic_geom_cp
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
    integer :: start_count, end_count, i_unit
    logical :: exists, found

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals, l1_up, l1_dn, d_l1_FD, l2_up, l2_dn, d_l2_FD,&
    a_up, a_dn, d_a_FD, g2_up, g2_dn, d_g2_FD, R1_up, R1_dn, d_R1_FD, R2_up, R2_dn, d_R2_FD,&
    dR_up, dR_dn, d_dR_FD
    real,dimension(:,:),allocatable ::  residuals3 

    integer :: i,j,k,m,n,p, N_verts, N_panels, vert, index
    real :: step
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
    test_geom = test_mesh%panels(1)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
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

    ! calc CALC BASIC GEOM geom sensitivity of relation between cp1 and panel1 
    adjoint_geom = adjoint_mesh%panels(1)%calc_subsonic_geom_adjoint(adjoint_mesh%cp(1)%loc,&
    adjoint_mesh%cp(1)%d_loc, adjoint_freestream_flow)
    !!!!!!!!!!!! END ADJOINT TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!


    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))
    
    step = 0.00001
    index = 1
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALC_SUBSONIC_GEOM (CONTROL POINTS) SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ CALC_SUBSONIC_GEOM (CONTROL POINTS) SENSITIVITIES TEST ---&
    ------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_a !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(l1_up(N_verts*3))
    allocate(l1_dn(N_verts*3))
    allocate(d_l1_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_l1 edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_l1 !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                l1_up(j + (i-1)*N_verts) = test_geom%l1(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                l1_dn(j + (i-1)*N_verts) = test_geom%l1(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_l1_FD = (l1_up - l1_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_l1"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_l1_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_l1 !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_l1 edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_l1              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_l1(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_l1(p)%elements(i)%value, &
            i, adjoint_geom%d_l1(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_l1(p)%get_value(i) - d_l1_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_l1 edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_l1                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_l1(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_l1 edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_l2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(l2_up(N_verts*3))
    allocate(l2_dn(N_verts*3))
    allocate(d_l2_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_l2 edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_l2 !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                l2_up(j + (i-1)*N_verts) = test_geom%l2(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                l2_dn(j + (i-1)*N_verts) = test_geom%l2(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_l2_FD = (l2_up - l2_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_l2"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_l2_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_l2 !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_l2 edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_l2              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_l2(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_l2(p)%elements(i)%value, &
            i, adjoint_geom%d_l2(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_l2(p)%get_value(i) - d_l2_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_l2 edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_l2                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_l2(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_l2 edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_a !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(a_up(N_verts*3))
    allocate(a_dn(N_verts*3))
    allocate(d_a_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_a edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_a !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                a_up(j + (i-1)*N_verts) = test_geom%a(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                a_dn(j + (i-1)*N_verts) = test_geom%a(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_a_FD = (a_up - a_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_a"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_a_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_a !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_a edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_a              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_a(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_a(p)%elements(i)%value, &
            i, adjoint_geom%d_a(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_a(p)%get_value(i) - d_a_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_a edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_a                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_a(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_a edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_g2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(g2_up(N_verts*3))
    allocate(g2_dn(N_verts*3))
    allocate(d_g2_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_g2 edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_g2 !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                g2_up(j + (i-1)*N_verts) = test_geom%g2(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                g2_dn(j + (i-1)*N_verts) = test_geom%g2(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_g2_FD = (g2_up - g2_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_g2"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_g2_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_g2 !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_g2 edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_g2              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_g2(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_g2(p)%elements(i)%value, &
            i, adjoint_geom%d_g2(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_g2(p)%get_value(i) - d_g2_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_g2 edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_g2                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_g2(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_g2 edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_R1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(R1_up(N_verts*3))
    allocate(R1_dn(N_verts*3))
    allocate(d_R1_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_R1 edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_R1 !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                R1_up(j + (i-1)*N_verts) = test_geom%R1(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                R1_dn(j + (i-1)*N_verts) = test_geom%R1(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_R1_FD = (R1_up - R1_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_R1"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_R1_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_R1 !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_R1 edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_R1              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_R1(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_R1(p)%elements(i)%value, &
            i, adjoint_geom%d_R1(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_R1(p)%get_value(i) - d_R1_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_R1 edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_R1                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_R1(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_R1 edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do
    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_R2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(R2_up(N_verts*3))
    allocate(R2_dn(N_verts*3))
    allocate(d_R2_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_R2 edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_R2 !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                R2_up(j + (i-1)*N_verts) = test_geom%R2(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                R2_dn(j + (i-1)*N_verts) = test_geom%R2(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_R2_FD = (R2_up - R2_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_R2"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_R2_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_R2 !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_R2 edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_R2              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_R2(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_R2(p)%elements(i)%value, &
            i, adjoint_geom%d_R2(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_R2(p)%get_value(i) - d_R2_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_R2 edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_R2                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_R2(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_R2 edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST CALC BASIC GEOM (panel 1, cp 1) d_dR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perturb x1 up
    allocate(dR_up(N_verts*3))
    allocate(dR_dn(N_verts*3))
    allocate(d_dR_FD(N_verts*3))
    
    do p=1,3

        write(*,*) "---------------------------------- TEST CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge ", p," -&
        --------------------------------"
        write(*,*) ""
        write(*,*) "the sensitivity of CALC BASIC GEOM (panel 1, cp 1) d_dR edge ", p, " WRT each design variable"
        write(*,*) ""

    
        !!!!!!!!! CENTRAL DIFFERENCE CALC BASIC GEOM (panel 1, cp 1) d_dR !!!!!!!!!
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        
    
        
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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)
    
                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
                
                ! put the x y or z component of the vertex of interest (index) in a list
                dR_up(j + (i-1)*N_verts) = test_geom%dR(p)

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
                test_geom = test_mesh%panels(index)%calc_subsonic_geom(test_mesh%cp(1)%loc,freestream_flow,.false.)

                !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!

                ! put the x y or z component of the vertex of interest (index) in a list
                dR_dn(j + (i-1)*N_verts) = test_geom%dR(p)
                
                ! restore geometry
                test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            end do 
        end do 
        
        ! central difference 
        d_dR_FD = (dR_up - dR_dn)/(2.*step)
                
        

        ! write results
        write(*,*) ""
    
        write(*,*) "--------------------------------------------------------------------------"
        write(*,'(A, I1)') "  CENTRAL DIFFERENCE CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge ",p
        write(*,*) "--------------------------------------------------------------------------"
        write(*,*) ""
        write(*,*) "  d_dR"
    
        do i = 1, N_verts*3
            write(*, '(f14.10, 4x)') d_dR_FD(i)
        end do 
        
        !!!!!!!!!! ADJOINT CALC BASIC GEOM (panel 1, cp 1) d_dR !!!!!!!!!!!!!
        write(*,*) ""
        write(*,*) "------------------------------------------------"
        write(*,'(A, I1)') "  ADJOINT  d_dR edge ",p
        write(*,*) "------------------------------------------------"
        write(*,*) ""
        
        !write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "  adjoint CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge ",p,"  (sparse)"
        write(*,*) ""
        write(*,*) "  d_dR              sparse_index       full_index"
        

        do i=1,adjoint_geom%d_dR(p)%sparse_size
            write(*,'((f14.10, 4x), 12x, I5, 12x, I5)') adjoint_geom%d_dR(p)%elements(i)%value, &
            i, adjoint_geom%d_dR(p)%elements(i)%full_index
        end do
        write(*,*) ""
        write(*,*) ""


        ! calculate residuals3
        do i =1, N_verts*3
            residuals(i) = adjoint_geom%d_dR(p)%get_value(i) - d_dR_FD(i)
        end do

        write(*,'(A, I1, A)') "  adjoint CALC SUNSONIC GEOM (panel 1, cp 1) d_dR edge ",p,",  expanded"
        write(*,*) ""
        write(*,*) "  d_dR                 residual"
        

        do i = 1, N_verts*3
            write(*, '((f14.10, 4x),3x, (f14.10, 4x))') adjoint_geom%d_dR(p)%get_value(i), residuals(i)
        end do
        write(*,*) ""
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
            if (p ==1) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge &
                1 test FAILED"
            elseif (p ==2) then
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge &
                2 test FAILED"
            else
                failure_log(total_tests-passed_tests) = "CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge &
                3 test FAILED"
            end if
            write(*,*) failure_log(total_tests-passed_tests)
        else
            write(*,'(A,I1,A)') "CALC SUBSONIC GEOM (panel 1, cp 1) d_dR edge ",p," test PASSED"
            passed_tests = passed_tests + 1
            total_tests = total_tests + 1
            
        end if
        test_failed = .false.
        write(*,*) "" 
        write(*,*) ""
    
    
    ! p loop
    end do


    !!!!!!!!!!!!!! CALC_SUBSONIC_GEOM (CONTROL POINTS) SENSITIVITIES RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------CALC_SUBSONIC_GEOM (CONTROL POINTS) SENSITIVITIES TEST RESULTS--------------"
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

end program calc_subsonic_geom_cp