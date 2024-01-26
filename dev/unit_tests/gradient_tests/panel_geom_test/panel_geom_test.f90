program gradient_test

    ! tests various intermediate sensitivities 
    use adjoint_mod
    use base_geom_mod
    use panel_mod
    use flow_mod
    
    implicit none

    type(json_file) :: input_json
    type(flow) :: freestream
    type(eval_point_geom) :: geom
    type(integrals) :: int
    real :: step
    real,dimension(3) :: values1, values2, values3, dummy_n, dummy_t, dummy_cross

    real,dimension(:),allocatable :: residuals, X_beta, loc_up, loc_dn, centr_up, centr_dn, normal_up, normal_dn, &
    area_up, area_dn, d_area_FD, d_g_up, d_g_dn, norm_d_g_up, norm_d_g_dn, t_hat_g_up, t_hat_g_dn, &
    n_hat_g_up, n_hat_g_dn, A_g_to_ls_up, A_g_to_ls_dn,  A_ls_to_g_up, A_ls_to_g_dn, vertices_ls_up, vertices_ls_dn, &
    d_vertices_ls_FD, n_hat_ls_up, n_hat_ls_dn, d_n_hat_ls_FD

    real,dimension(:,:),allocatable :: v, vertex_locs, d_loc_FD, d_centr_FD, d_n_g_FD, d_norm_d_g_FD, residuals3 

    real,dimension(:,:,:),allocatable :: d_d_g_FD, d_t_hat_g_FD, d_n_hat_g_FD, dt_cross_n_FD,t_cross_dn_FD, &
    d_A_g_to_ls_FD, d_A_ls_to_g_FD

    integer :: i,j,k,m,n, N_verts, N_panels, vert, index
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels    ! list of panels, this should be a mesh attribute
    character(len=:),allocatable :: spanwise_axis

    type(sparse_matrix),dimension(3) :: d_n_hat_g_adjoint
    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(20) :: failure_log
    character(len=10) :: m_char



    test_failed = .true. ! assume test failed, if the test condition is met, test passed
    ! NOTE: on the sparse vector test, I assume the test passes, if it fails a test condition, test fails
    passed_tests = 0
    total_tests = 0

    
    ! Initialize freestream
    freestream%supersonic = .false.
    freestream%v_inf = (/1., 0., 0./)
    freestream%c_hat_g = (/1., 0., 0./) ! freestream direction vector
    freestream%U = 1.                   ! freestream magnitude
    freestream%U_inv = 1.
    freestream%M_inf = 0.               ! incompressible, Mach = 0
    freestream%B = 1.                   ! supersonic compressibility scaling factor
    freestream%s = 1.                   ! flow type indicator (1 = subsonic)
    freestream%K = 4.*pi                ! Kappa factor (4 pi) for subsonic, (2 pi) for supersonic
    freestream%K_inv = 1./(4.*pi)       
    spanwise_axis = 'y+'
    call freestream%calc_metric_matrices()
    call freestream%calc_transforms(spanwise_axis)


    ! For this simple example, we don't need to collapse dublicates
    N_verts = 6
    N_panels = 8
    
    allocate(v(3,N_verts))
    ! Initialize vertices
    v(:,1) = (/ 0.0, 0.0,-0.2/)  ! top
    v(:,2) = (/ 1.0, 0.0, 0.0/)  ! forward 
    v(:,3) = (/ 0.0, 1.0, 0.0/)  ! right wingtip
    v(:,4) = (/-1.0, 0.0, 0.0/)  ! aft
    v(:,5) = (/ 0.0,-1.0, 0.0/)  ! left wingtip
    v(:,6) = (/ 0.0, 0.0, 0.2/)  ! bottom
    
    ! put vertex_locs in a list
    allocate(vertex_locs(3,N_verts))
    do i =1, N_verts
        vertex_locs(:,i) = v(:,i)
    end do

    ! populate vertices list
    allocate(vertices(N_verts))
    do i=1,N_verts
        call vertices(i)%init(vertex_locs(:,i), i)
    end do
    
    ! build design variable vector X_beta
    allocate(X_beta(N_verts*3))
    do i=1,3
        do j=1,N_verts
            X_beta(j + (i-1)*N_verts) = vertices(j)%loc(i)
        end do
    end do

    ! initialize 8 panels
    ! panel init calculates the following:
    !   normal vector
    !   area
    !   centroid
    !   radius
    !   n_hat_g (g edge vectors)
    allocate(panels(N_panels))
    
    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward

    
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

    ! sensitivity to vertex 1
    step = 0.0001
    
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
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

                ! put the x y or z component of the vertex of interest (index) in a list
                loc_up(j + (i-1)*N_verts) = vertices(index)%loc(k)

                ! perturb down the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                ! put the x y or z component of the vertex of interest (index) in a list
                loc_dn(j + (i-1)*N_verts) = vertices(index)%loc(k)
                
                ! restore geometry
                vertices(j)%loc(i) = vertices(j)%loc(i) + step
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
                
    ! call vertex init (do it for all vertices because they will be used later)
    do i =1,N_verts
        call vertices(i)%init_adjoint(N_verts)
    end do
    
    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_loc point 1"
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z             sparse_index       full_index"
    do i=1,vertices(index)%d_loc%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') vertices(index)%d_loc%columns(i)%vector_values(:), &
        i, vertices(index)%d_loc%columns(i)%full_index
    end do
    write(*,*) ""



    allocate(residuals3(3,N_verts*3))
    
    ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = vertices(index)%d_loc%get_values(i) - d_loc_FD(:,i)
    end do

    
    write(*,*) "         d_loc vertex 1 expanded "
    write(*,*) "  d_loc_x           d_loc_y           d_loc_z                                 residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') vertices(index)%d_loc%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (any(residuals3(:,i) > 1.0e-12)) then
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

    ! we want the sensitivity of the centroid of panel 1 WRT X(beta)
    index = 1

    ! sensitivity to vertex 1
    step = 0.0001

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
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

                ! update panel centroid calculations
                call panels(index)%calc_centroid()

                ! put the x y or z component of the panel's perturbed centroid in a list
                centr_up(j + (i-1)*N_verts) = panels(index)%centr(k)

                ! perturb down the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                ! update panel centroid calculations
                call panels(index)%calc_centroid()
                
                ! put the x y or z component of the panel's perturbed centroid in a list
                centr_dn(j + (i-1)*N_verts) = panels(index)%centr(k)
                
                ! restore geometry
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

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
  
    ! calculate d_centr
    
    ! calculate d_centr (do it for all panels because they will be used later)
    do i =1,N_panels
        call panels(i)%calc_d_centr()
    end do

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_centr panel 1"
    write(*,*) "  d_centr_x         d_centr_y         d_centr_z           sparse_index       full_index"
    do i=1,panels(index)%d_centr%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_centr%columns(i)%vector_values(:), &
        i, panels(index)%d_centr%columns(i)%full_index
    end do
    write(*,*) ""

   ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = panels(index)%d_centr%get_values(i) - d_centr_FD(:,i)
    end do
    

    write(*,*) "         d_centr panel 1 expanded "
    write(*,*) "  d_centr_x         d_centr_y         d_centr_z                               residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_centr%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (any(residuals3(:,i) > 1.0e-12)) then
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

    ! we want the sensitivity of the normal of panel 1 WRT X(beta)
    index = 1

    ! sensitivity to vertex 1
    step = 0.0001

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
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

                ! update panel normal calculations
                call panels(index)%calc_normal()

                ! put the x y or z component of the panel's perturbed normal in a list
                normal_up(j + (i-1)*N_verts) = panels(index)%n_g(k)

                ! perturb down the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                ! update panel normal calculations
                call panels(index)%calc_normal()
                
                ! put the x y or z component of the panel's perturbed normal in a list
                normal_dn(j + (i-1)*N_verts) = panels(index)%n_g(k)
            
                ! restore geometry
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

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
            
    ! calculate d_n_g (do it for all panels because they will be used later)
    do i =1,N_panels
        call panels(i)%calc_d_normal_and_d_area()
    end do

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_n_g panel 1"
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z             sparse_index       full_index"
    do i=1,panels(index)%d_n_g%sparse_num_cols
        write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_n_g%columns(i)%vector_values(:), &
        i, panels(index)%d_n_g%columns(i)%full_index
    end do
    write(*,*) ""

    ! calculate residuals3
    do i =1, N_verts*3
        residuals3(:,i) = panels(index)%d_n_g%get_values(i) - d_n_g_FD(:,i)
    end do

    write(*,*) "         d_n_g panel 1 expanded "
    write(*,*) "  d_n_g_x           d_n_g_y           d_n_g_z                                residuals"
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_n_g%get_values(i), residuals3(:,i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (any(residuals3(:,i) > 1.0e-8)) then
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

    ! we want the sensitivity of the area of panel 1 WRT X(beta)
    index = 1

    ! sensitivity to vertex 1
    step = 0.0001

    ! perturb x1 up
    allocate(area_up(N_verts*3))
    allocate(area_dn(N_verts*3))
    allocate(d_area_FD(N_verts*3))

    ! do for each design variable 
    do i=1,3
        do j=1,N_verts
            ! perturb up the current design variable
            vertices(j)%loc(i) = vertices(j)%loc(i) + step

            ! update panel area calculations
            call panels(index)%calc_area()

            ! put the x y or z component of the panel's perturbed area in a list
            area_up(j + (i-1)*N_verts) = panels(index)%A

            ! perturb down the current design variable
            vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

            ! update panel area calculations
            call panels(index)%calc_area()
            
            ! put the x y or z component of the panel's perturbed area in a list
            area_dn(j + (i-1)*N_verts) = panels(index)%A
                    
            ! restore geometry
            vertices(j)%loc(i) = vertices(j)%loc(i) + step

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
            
    ! calculate d_area for all panels (this was done in the previous d_normal test)
    

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_area panel 1"
    write(*,*) "  d_area                  sparse_index       full_index"
    do i=1,panels(index)%d_A%sparse_size
        write(*,'(f14.10, 20x, I5, 12x, I5)') panels(index)%d_A%elements(i)%value, &
        i, panels(index)%d_A%elements(i)%full_index
    end do
    write(*,*) ""

    ! calculate residuals

    allocate(residuals(N_verts*3))

    do i =1, N_verts*3
        residuals(i) = panels(index)%d_A%get_value(i) - d_area_FD(i)
    end do

    write(*,*) "         d_area panel 1 expanded "
    write(*,*) "  d_area         residuals"
    do i = 1, N_verts*3
        write(*, '(f14.10,3x, f14.10)') panels(index)%d_A%get_value(i), residuals(i)
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,N_verts*3
        if (residuals(i) > 1.0e-8) then
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

        ! we want the sensitivity of the edge vector (global) of panel 1 edge 1 WRT X(beta)
        index = 1

        ! sensitivity to vertex 1
        step = 0.0001

        ! for each x, y, z of n_hat_g (edge m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(panels)
                    allocate(panels(N_panels))
                    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                    
                    deallocate(panels(index)%n_hat_g)
                    call panels(index)%calc_g_edge_vectors()

                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    n_hat_g_up(j + (i-1)*N_verts) = panels(index)%n_hat_g(k,m)

                    ! perturb down the current design variable
                    vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                
                    ! update panel geometry and calcs
                    deallocate(panels)
                    allocate(panels(N_panels))
                    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                    
                    deallocate(panels(index)%n_hat_g)
                    call panels(index)%calc_g_edge_vectors()
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    n_hat_g_dn(j + (i-1)*N_verts) = panels(index)%n_hat_g(k,m)
                
                    ! restore geometry
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    

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


        ! Re init the body
        deallocate(panels)
        allocate(panels(N_panels))
        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
        
        
        ! calculate d_n_hat_g 
        call panels(index)%calc_derived_geom_adjoint()

        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_n_hat_g panel 1 (edge ", m, ")"
        write(*,*) "  d_n_hat_g_x           d_n_hat_g_y           d_n_hat_g_z             sparse_index       full_index"
        do i=1,panels(index)%d_n_hat_g(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_n_hat_g(m)%columns(i)%vector_values(:), &
            i, panels(index)%d_n_hat_g(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = panels(index)%d_n_hat_g(m)%get_values(i) - d_n_hat_g_FD(:,i,m)
        end do

        write(*,'(A, I1, A)') "         d_n_hat_g panel 1 (edge ", m, ") expanded "
        write(*,*) "  d_n_hat_g_x           d_n_hat_g_y           d_n_hat_g_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_n_hat_g(m)%get_values(i), residuals3(:,i)
        end do
        write(*,*) ""


        ! check if test failed
        do i=1,N_verts*3
            if (any(residuals3(:,i) > 1.0e-8)) then
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

    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_A_g_to_ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_A_g_to_ls ----------------------------------"
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

        ! we want the sensitivity of the A_g_to_ls panel 1 row m WRT X(beta)
        index = 1

        step = 0.00001

        ! for each x, y, z of A_g_to_ls (row m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(panels)
                    allocate(panels(N_panels))
                    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                    
                    call panels(index)%init_with_flow(freestream, .false., 0)

                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_g_to_ls_up(j + (i-1)*N_verts) = panels(index)%A_g_to_ls(m,k)

                    ! perturb down the current design variable
                    vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                
                    ! update panel geometry and calcs
                    deallocate(panels)
                    allocate(panels(N_panels))
                    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                    
                    call panels(index)%init_with_flow(freestream, .false., 0)
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_g_to_ls_dn(j + (i-1)*N_verts) = panels(index)%A_g_to_ls(m,k)
                
                    ! restore geometry
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    

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


        ! Re init the body
        deallocate(panels)
        allocate(panels(N_panels))
        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
        
        call panels(1)%init_with_flow(freestream, .false., 0)

        ! calculate d_A_g_to_ls 
        call panels(index)%init_adjoint(freestream)

        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_g_to_ls panel 1 (row ", m, ")"
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y           d_A_g_to_ls_z             sparse_index       full_index"
        do i=1,panels(index)%d_A_g_to_ls(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_A_g_to_ls(m)%columns(i)%vector_values(:), &
            i, panels(index)%d_A_g_to_ls(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = panels(index)%d_A_g_to_ls(m)%get_values(i) - d_A_g_to_ls_FD(m,i,:)
        end do

        write(*,'(A, I1, A)') "         d_A_g_to_ls panel 1 (row ", m, ") expanded "
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y           d_A_g_to_ls_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_A_g_to_ls(m)%get_values(i), residuals3(:,i)
        end do
        write(*,*) ""


        ! check if test failed
        do i=1,N_verts*3
            if (any(residuals3(:,i) > 1.0e-8)) then
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

        ! we want the sensitivity of the A_ls_to_g panel 1 row m WRT X(beta)
        index = 1

        ! sensitivity to vertex 1
        step = 0.00001

        ! for each x, y, z of A_ls_to_g (row m) 
        do k=1,3
            ! do for each design variable
            do i=1,3
                do j=1,N_verts
                    ! perturb up the current design variable
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    ! update panel geometry and calcs
                    deallocate(panels)
                    allocate(panels(N_panels))
                    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                    
                    call panels(index)%init_with_flow(freestream, .false., 0)

                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_ls_to_g_up(j + (i-1)*N_verts) = panels(index)%A_ls_to_g(m,k)

                    ! perturb down the current design variable
                    vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                
                    ! update panel geometry and calcs
                    deallocate(panels)
                    allocate(panels(N_panels))
                    call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                    call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                    call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                    call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                    call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                    call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                    call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                    call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                    
                    call panels(index)%init_with_flow(freestream, .false., 0)
                    
                    ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                    A_ls_to_g_dn(j + (i-1)*N_verts) = panels(index)%A_ls_to_g(m,k)
                
                    ! restore geometry
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    

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


        ! Re init the body
        deallocate(panels)
        allocate(panels(N_panels))
        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
        
        call panels(1)%init_with_flow(freestream, .false., 0)

        ! calculate d_A_ls_to_g 
        call panels(index)%init_adjoint(freestream)

        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_ls_to_g panel 1 (row ", m, ")"
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y           d_A_ls_to_g_z             sparse_index       full_index"
        do i=1,panels(index)%d_A_ls_to_g(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_A_ls_to_g(m)%columns(i)%vector_values(:), &
            i, panels(index)%d_A_ls_to_g(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = panels(index)%d_A_ls_to_g(m)%get_values(i) - d_A_ls_to_g_FD(m,i,:)
        end do

        write(*,'(A, I1, A)') "         d_A_ls_to_g panel 1 (row ", m, ") expanded "
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y           d_A_ls_to_g_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_A_ls_to_g(m)%get_values(i), residuals3(:,i)
        end do
        write(*,*) ""


        ! check if test failed
        do i=1,N_verts*3
            if (any(residuals3(:,i) > 1.0e-8)) then
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
            ! we want the sensitivity of the coordinates in local scaled of  panel 1 vertex 1 WRT X(beta)
            index = 1

            ! sensitivity to vertex 1
            step = 0.0001

            ! ! for each x, y, z of n_hat_g (edge m) 
            ! do k=1,3
                ! do for each design variable
                do i=1,3
                    do j=1,N_verts
                        ! perturb up the current design variable
                        vertices(j)%loc(i) = vertices(j)%loc(i) + step
                        
                        ! update panel geometry and calcs
                        deallocate(panels)
                        allocate(panels(N_panels))
                        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                        
                        call panels(index)%init_with_flow(freestream, .false., 0)

                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        vertices_ls_up(j + (i-1)*N_verts) = panels(index)%vertices_ls(n,m)

                        ! perturb down the current design variable
                        vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                    
                        ! update panel geometry and calcs
                        deallocate(panels)
                        allocate(panels(N_panels))
                        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                        
                        call panels(index)%init_with_flow(freestream, .false., 0)
                        
                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        vertices_ls_dn(j + (i-1)*N_verts) = panels(index)%vertices_ls(n,m)
                    
                        ! restore geometry
                        vertices(j)%loc(i) = vertices(j)%loc(i) + step
                        
                        

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


            ! Re init the body
            deallocate(panels)
            allocate(panels(N_panels))
            call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
            call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
            call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
            call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

            call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
            call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
            call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
            call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
            
            call panels(1)%init_with_flow(freestream, .false., 0)

            ! calculate d_vertices_ls 
            call panels(index)%init_adjoint(freestream)

            
            ! write sparse matrix
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex ", m, " eta coordinate)"
            end if
            write(*,*) "  sparse value                  sparse_index       full_index"
            do i=1,panels(index)%d_vertices_ls(n,m)%sparse_size
                write(*,'(f14.10, 20x, I5, 12x, I5)') panels(index)%d_vertices_ls(n,m)%elements(i)%value, &
                i, panels(index)%d_vertices_ls(n,m)%elements(i)%full_index
            end do
            write(*,*) ""

            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = panels(index)%d_vertices_ls(n,m)%get_value(i) - d_vertices_ls_FD(i)
            end do

            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex", m, " xi coordinate) expanded"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex ", m, " eta coordinate) expanded"
            end if
            write(*,*) "  adjoint value         residuals"
            do i = 1, N_verts*3
                write(*, '(f14.10,3x, f14.10)') panels(index)%d_vertices_ls(n,m)%get_value(i), residuals(i)
            end do
            write(*,*) ""


            ! check if test failed
            do i=1,N_verts*3
                if (residuals(i) > 1.0e-8) then
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
            ! we want the sensitivity of the edge vector of panel 1 edge 1 WRT X(beta)
            index = 1

            ! sensitivity to edge 1
            step = 0.0001

            ! ! for each x, y, z of n_hat_ls (edge m) 
            ! do k=1,3
                ! do for each design variable
                do i=1,3
                    do j=1,N_verts
                        ! perturb up the current design variable
                        vertices(j)%loc(i) = vertices(j)%loc(i) + step
                        
                        ! update panel geometry and calcs
                        deallocate(panels)
                        allocate(panels(N_panels))
                        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                        
                        call panels(index)%init_with_flow(freestream, .false., 0)

                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        n_hat_ls_up(j + (i-1)*N_verts) = panels(index)%n_hat_ls(n,m)

                        ! perturb down the current design variable
                        vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                    
                        ! update panel geometry and calcs
                        deallocate(panels)
                        allocate(panels(N_panels))
                        call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
                        call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
                        call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
                        call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

                        call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
                        call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
                        call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
                        call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                        
                        call panels(index)%init_with_flow(freestream, .false., 0)
                        
                        ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
                        n_hat_ls_dn(j + (i-1)*N_verts) = panels(index)%n_hat_ls(n,m)
                    
                        ! restore geometry
                        vertices(j)%loc(i) = vertices(j)%loc(i) + step
                        
                        

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


            ! Re init the body
            deallocate(panels)
            allocate(panels(N_panels))
            call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
            call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
            call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
            call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

            call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
            call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
            call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
            call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
            
            call panels(1)%init_with_flow(freestream, .false., 0)

            ! calculate d_vertices_ls 
            call panels(index)%init_adjoint(freestream)

            
            ! write sparse matrix
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge ", m, " eta coordinate)"
            end if
            write(*,*) "  sparse value              sparse_index       full_index"
            do i=1,panels(index)%d_n_hat_ls(n,m)%sparse_size
                write(*,'(f14.10, 20x, I5, 12x, I5)') panels(index)%d_n_hat_ls(n,m)%elements(i)%value, &
                i, panels(index)%d_n_hat_ls(n,m)%elements(i)%full_index
            end do
            write(*,*) ""

            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = panels(index)%d_n_hat_ls(n,m)%get_value(i) - d_n_hat_ls_FD(i)
            end do

            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge", m, " xi coordinate) expanded"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge ", m, " eta coordinate) expanded"
            end if
            write(*,*) "  adjoint value         residuals"
            do i = 1, N_verts*3
                write(*, '(f14.10,3x, f14.10)') panels(index)%d_n_hat_ls(n,m)%get_value(i), residuals(i)
            end do
            write(*,*) ""


            ! check if test failed
            do i=1,N_verts*3
                if (residuals(i) > 1.0e-8) then
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_d_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) "---------------------------------- TEST d_d_g ----------------------------------"
!     write(*,*) ""
!     write(*,*) "the sensitivity of the edge vectors of panel 1 WRT each design variable"
!     write(*,*) ""
    
!     ! allocate central difference variables
!     allocate(d_g_up(N_verts*3))
!     allocate(d_g_dn(N_verts*3))
!     allocate(d_d_g_FD(3,N_verts*3,3))

!     deallocate(panels(index)%n_hat_g)


!     ! do for each edge
!     do m=1,3

!         !!!!!!!!! Finite Difference  d_d_g (edge m) !!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_d_g (edge ", m, ")"

!         ! we want the sensitivity of the edge vector of panel 1 edge 1 WRT X(beta)
!         index = 1

!         ! sensitivity to vertex 1
!         step = 0.0001

!         ! for each x, y, z of n_hat_g (edge m) 
!         do k=1,3
!             ! do for each design variable
!             do i=1,3
!                 do j=1,N_verts
!                     ! perturb up the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step

!                     ! update panel edge outward normal unit vector calculations

!                     call panels(index)%calc_g_edge_vectors()

!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     d_g_up(j + (i-1)*N_verts) = panels(index)%d_g(k,m)

!                     ! perturb down the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                     ! update panel edge outward normal unit vector calculations
!                     call panels(index)%calc_g_edge_vectors()
                    
!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     d_g_dn(j + (i-1)*N_verts) = panels(index)%d_g(k,m)
                
!                     ! restore geometry
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                 end do 
!             end do 

!             ! central difference 
!             d_d_g_FD(k,:,m) = (d_g_up - d_g_dn)/(2.*step)

!         end do

!         ! write results
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          d_d_g_FD panel 1 (edge ", m, ")"
!         write(*,*) "  d_d_g_x           d_d_g_y            d_d_g_z "
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x))') d_d_g_FD(:,i, m)
!         end do 


!         !!!!!!!!!! ADJOINT d_d_g (edge m)!!!!!!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  ADJOINT d_d_g (edge ", m, ")"
!         write(*,*) ""

!         ! only need to calc this once:
!         if (m==1) then

!             ! calculate d_n_hat_g (do it for all panels because they will be used later)
!             do i =1,N_panels
!                 call panels(i)%calc_d_n_hat_g()
!             end do

!         end if

!         ! write sparse matrix
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          d_d_g panel 1 (edge ", m, ")"
!         write(*,*) "  d_d_g_x           d_d_g_y           d_d_g_z             sparse_index       full_index"
!         do i=1,panels(index)%d_d_g(m)%sparse_num_cols
!             write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_d_g(m)%columns(i)%vector_values(:), &
!             i, panels(index)%d_d_g(m)%columns(i)%full_index
!         end do
!         write(*,*) ""

!         ! calculate residuals3
!         do i =1, N_verts*3
!             residuals3(:,i) = panels(index)%d_d_g(m)%get_values(i) - d_d_g_FD(:,i,m)
!         end do

!         write(*,'(A, I1, A)') "         d_d_g panel 1 (edge ", m, ") expanded "
!         write(*,*) "  d_d_g_x           d_d_g_y           d_d_g_z                            residuals"
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_d_g(m)%get_values(i), residuals3(:,i)
!         end do
!         write(*,*) ""


!         ! check if test failed
!         do i=1,N_verts*3
!             if (any(residuals3(:,i) > 1.0e-8)) then
!                 test_failed = .true.
!                 exit
!             else 
!                 test_failed = .false.
!             end if
!         end do
!         if (test_failed) then
!             total_tests = total_tests + 1
!             write(m_char,'(I1)') m
!             failure_log(total_tests-passed_tests) = "d_d_g (edge "// trim(m_char) // ") test FAILED"
!             write(*,*) failure_log(total_tests-passed_tests)
!         else
!             write(*,'(A, I1, A)') "d_d_g (edge ",m,") test PASSED"
!             passed_tests = passed_tests + 1
!             total_tests = total_tests + 1
!         end if
!         test_failed = .false.
!         write(*,*) "" 
!         write(*,*) ""


!         !copy d_n_hat_g_adjoint data
!         call d_n_hat_g_adjoint(m)%init_from_sparse_matrix(panels(index)%d_n_hat_g(m))

!     ! end panel edge loop
!     end do

    
    

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_norm_d_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) "---------------------------------- TEST d_norm_d_g ----------------------------------"
!     write(*,*) ""
!     write(*,*) "the sensitivity of the edge vector normals of panel 1 WRT each design variable"
!     write(*,*) ""
    
!     ! allocate central difference variables
!     allocate(norm_d_g_up(N_verts*3))
!     allocate(norm_d_g_dn(N_verts*3))
!     allocate(d_norm_d_g_FD(N_verts*3, 3))

!     !deallocate(panels(index)%n_hat_g)
    


!     ! do for each edge
!     do m=1,3

!         !!!!!!!!! Finite Difference  d_norm_d_g (edge m) !!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_norm_d_g (edge ", m, ")"

!         ! we want the sensitivity of the edge vector normals of panel 1 edge m WRT X(beta)
!         index = 1

!         ! sensitivity to vertex 1
!         step = 0.0001

!         ! for each x, y, z of n_hat_g (edge m) 
!         !do k=1,3
!             ! do for each design variable
!             do i=1,3
!                 do j=1,N_verts
!                     ! perturb up the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step

!                     ! update panel edge outward normal unit vector calculations

!                     call panels(index)%calc_g_edge_vectors()

!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     norm_d_g_up(j + (i-1)*N_verts) = panels(index)%norm_d_g(m)

!                     ! perturb down the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                     ! update panel edge outward normal unit vector calculations
!                     call panels(index)%calc_g_edge_vectors()
                    
!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     norm_d_g_dn(j + (i-1)*N_verts) = panels(index)%norm_d_g(m)
                
!                     ! restore geometry
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                 end do 
!             end do 

!             ! central difference 
!             d_norm_d_g_FD(:,m) = (norm_d_g_up - norm_d_g_dn)/(2.*step)

!         !end do

!         ! write results
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          TEST d_norm_d_g_FD panel 1 (edge ", m, ")"
!         write(*,*) "  d_norm_d_g "
!         do i = 1, N_verts*3
!             write(*, '((f14.10, 4x))') d_norm_d_g_FD(i, m)
!         end do 


!         !!!!!!!!!! ADJOINT d_norm_d_g (edge m)!!!!!!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  ADJOINT d_norm_d_g (edge ", m, ")"
!         write(*,*) ""

!         ! ! only need to calc this once:
!         ! if (m==1) then

!         !     ! calculate d_norm_d_g (do it for all panels because they will be used later)
!         !     do i =1,N_panels
!         !         call panels(i)%calc_d_n_hat_g()
!         !     end do

!         ! end if

!         ! write sparse matrix
!         write(*,*) ""
!         write(*,'(A, I1, A)') "         d_norm_d_g panel 1 (edge ", m, ")"
!         write(*,*) "  d_norm_d_g         sparse_index       full_index"
!         do i=1,panels(index)%d_norm_d_g(m)%sparse_size
!             write(*,'(f14.10, 12x, I5, 12x, I5)') panels(index)%d_norm_d_g(m)%elements(i)%value, &
!             i, panels(index)%d_norm_d_g(m)%elements(i)%full_index
!         end do
!         write(*,*) ""

!         ! calculate residuals
!         do i =1, N_verts*3
!             residuals(i) = panels(index)%d_norm_d_g(m)%get_value(i) - d_norm_d_g_FD(i,m)
!         end do

!         write(*,'(A, I1, A)') "         d_norm_d_g panel 1 (edge ", m, ") expanded "
!         write(*,*) "  d_norm_d_g                            residuals"
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_norm_d_g(m)%get_value(i), residuals(i)
!         end do
!         write(*,*) ""


!         ! check if test failed
!         do i=1,N_verts*3
!             if (residuals(i) > 1.0e-8) then
!                 test_failed = .true.
!                 exit
!             else 
!                 test_failed = .false.
!             end if
!         end do
!         if (test_failed) then
!             total_tests = total_tests + 1
!             write(m_char,'(I1)') m
!             failure_log(total_tests-passed_tests) = "d_norm_d_g (edge "// trim(m_char) // ") test FAILED"
!             write(*,*) failure_log(total_tests-passed_tests)
!         else
!             write(*,'(A, I1, A)') "d_norm_d_g (edge ",m,") test PASSED"
!             passed_tests = passed_tests + 1
!             total_tests = total_tests + 1
!         end if
!         test_failed = .false.
!         write(*,*) "" 
!         write(*,*) ""

!     ! end panel edge loop
!     end do



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_t_hat_g !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) "---------------------------------- TEST d_t_hat_g ----------------------------------"
!     write(*,*) ""
!     write(*,*) "the sensitivity of the edge vectors of panel 1 WRT each design variable"
!     write(*,*) ""
    
!     ! allocate central difference variables
!     allocate(t_hat_g_up(N_verts*3))
!     allocate(t_hat_g_dn(N_verts*3))
!     allocate(d_t_hat_g_FD(3,N_verts*3,3))

!     !deallocate(panels(index)%n_hat_g)


!     ! do for each edge
!     do m=1,3

!         !!!!!!!!! Finite Difference  d_t_hat_g (edge m) !!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_t_hat_g (edge ", m, ")"

!         ! we want the sensitivity of the edge vector of panel 1 edge 1 WRT X(beta)
!         index = 1

!         ! sensitivity to vertex 1
!         step = 0.0001

!         ! for each x, y, z of n_hat_g (edge m) 
!         do k=1,3
!             ! do for each design variable
!             do i=1,3
!                 do j=1,N_verts
!                     ! perturb up the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step

!                     ! update panel edge outward normal unit vector calculations

!                     call panels(index)%calc_g_edge_vectors()

!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     t_hat_g_up(j + (i-1)*N_verts) = panels(index)%t_hat_g(k,m)

!                     ! perturb down the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                     ! update panel edge outward normal unit vector calculations
!                     call panels(index)%calc_g_edge_vectors()
                    
!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     t_hat_g_dn(j + (i-1)*N_verts) = panels(index)%t_hat_g(k,m)
                
!                     ! restore geometry
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                 end do 
!             end do 

!             ! central difference 
!             d_t_hat_g_FD(k,:,m) = (t_hat_g_up - t_hat_g_dn)/(2.*step)

!         end do

!         ! write results
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          d_t_hat_g_FD panel 1 (edge ", m, ")"
!         write(*,*) "  d_t_hat_g_x           d_t_hat_g_y            d_t_hat_g_z "
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x))') d_t_hat_g_FD(:,i, m)
!         end do 


!         !!!!!!!!!! ADJOINT d_n_hat_g (edge m)!!!!!!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  ADJOINT d_t_hat_g (edge ", m, ")"
!         write(*,*) ""

!         ! ! only need to calc this once:
!         ! if (m==1) then

!         !     ! calculate d_n_hat_g (do it for all panels because they will be used later)
!         !     do i =1,N_panels
!         !         call panels(i)%calc_d_n_hat_g()
!         !     end do

!         ! end if

!         ! write sparse matrix
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          d_t_hat_g panel 1 (edge ", m, ")"
!         write(*,*) "  d_t_hat_g_x           d_t_hat_g_y           d_t_hat_g_z             sparse_index       full_index"
!         do i=1,panels(index)%d_t_hat_g(m)%sparse_num_cols
!             write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_t_hat_g(m)%columns(i)%vector_values(:), &
!             i, panels(index)%d_t_hat_g(m)%columns(i)%full_index
!         end do
!         write(*,*) ""

!         ! calculate residuals3
!         do i =1, N_verts*3
!             residuals3(:,i) = panels(index)%d_t_hat_g(m)%get_values(i) - d_t_hat_g_FD(:,i,m)
!         end do

!         write(*,'(A, I1, A)') "         d_t_hat_g panel 1 (edge ", m, ") expanded "
!         write(*,*) "  d_t_hat_g_x           d_t_hat_g_y           d_t_hat_g_z                            residuals"
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_t_hat_g(m)%get_values(i), residuals3(:,i)
!         end do
!         write(*,*) ""


!         ! check if test failed
!         do i=1,N_verts*3
!             if (any(residuals3(:,i) > 1.0e-8)) then
!                 test_failed = .true.
!                 exit
!             else 
!                 test_failed = .false.
!             end if
!         end do
!         if (test_failed) then
!             total_tests = total_tests + 1
!             write(m_char,'(I1)') m
!             failure_log(total_tests-passed_tests) = "d_t_hat_g (edge "// trim(m_char) // ") test FAILED"
!             write(*,*) failure_log(total_tests-passed_tests)
!         else
!             write(*,'(A, I1, A)') "d_t_hat_g (edge ",m,") test PASSED"
!             passed_tests = passed_tests + 1
!             total_tests = total_tests + 1
!         end if
!         test_failed = .false.
!         write(*,*) "" 
!         write(*,*) ""

!     ! end panel edge loop
!     end do


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST dt_cross_n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) "---------------------------------- TEST dt_cross_n ----------------------------------"
!     write(*,*) ""
!     write(*,*) "the sensitivity of the edge vectors of panel 1 WRT each design variable"
!     write(*,*) ""
    
!     ! allocate central difference variables
    
!     allocate(dt_cross_n_FD(3,N_verts*3,3))

!     !deallocate(panels(index)%n_hat_g)


!     ! do for each edge
!     do m=1,3

!         !!!!!!!!! Finite Difference  dt_cross_n (edge m) !!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE dt_cross_n (edge ", m, ")"

!         ! we want the sensitivity of the edge vector of panel 1 edge 1 WRT X(beta)
!         index = 1

!         ! sensitivity to vertex 1
!         step = 0.0001

!         ! for each x, y, z of n_hat_g (edge m) 
!         do k=1,3
!             ! do for each design variable
!             do i=1,3
!                 do j=1,N_verts
!                     ! perturb up the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step

!                     ! update panel edge outward normal unit vector calculations

!                     call panels(index)%calc_g_edge_vectors()

!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     t_hat_g_up(j + (i-1)*N_verts) = panels(index)%t_hat_g(k,m)

!                     ! perturb down the current design variable
!                     vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                     ! update panel edge outward normal unit vector calculations
!                     call panels(index)%calc_g_edge_vectors()
                    
!                     ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                     t_hat_g_dn(j + (i-1)*N_verts) = panels(index)%t_hat_g(k,m)
                
!                     ! restore geometry
!                     vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
!                     ! panel_calc_g_edge_vectors allocates for n_hat_g each time, so deallocate old values
!                     deallocate(panels(index)%n_hat_g)

!                 end do 
!             end do 

!             ! central difference 
!             dt_cross_n_FD(k,:,m) = (t_hat_g_up - t_hat_g_dn)/(2.*step)

!         end do
        
       

!         ! write results
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          dt_cross_n_FD panel 1 (edge ", m, ")"
!         write(*,*) "  dt_cross_n_x           dt_cross_n_y            dt_cross_n_z "
!         do i = 1, N_verts*3
!             dt_cross_n_FD(:,i, m) = cross(dt_cross_n_FD(:,i, m), panels(index)%n_g)
!             write(*, '(3(f14.10, 4x))') dt_cross_n_FD(:,i, m)
!         end do 


!         !!!!!!!!!! ADJOINT d_n_hat_g (edge m)!!!!!!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  ADJOINT dt_cross_n (edge ", m, ")"
!         write(*,*) ""

!         ! ! only need to calc this once:
!         ! if (m==1) then

!         !     ! calculate d_n_hat_g (do it for all panels because they will be used later)
!         !     do i =1,N_panels
!         !         call panels(i)%calc_d_n_hat_g()
!         !     end do

!         ! end if

!         ! write sparse matrix
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          dt_cross_n panel 1 (edge ", m, ")"
!         write(*,*) "  dt_cross_n_x           dt_cross_n_y           dt_cross_n_z             sparse_index       full_index"
!         do i=1,panels(index)%dt_cross_n(m)%sparse_num_cols
!             write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%dt_cross_n(m)%columns(i)%vector_values(:), &
!             i, panels(index)%dt_cross_n(m)%columns(i)%full_index
!         end do
!         write(*,*) ""

!         ! calculate residuals3
!         do i =1, N_verts*3
!             residuals3(:,i) = panels(index)%dt_cross_n(m)%get_values(i) - dt_cross_n_FD(:,i, m)
!         end do

!         write(*,'(A, I1, A)') "         dt_cross_n panel 1 (edge ", m, ") expanded "
!         write(*,*) "  dt_cross_n_x           dt_cross_n_y           dt_cross_n_z                            residuals"
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%dt_cross_n(m)%get_values(i), residuals3(:,i)
!         end do
!         write(*,*) ""


!         ! check if test failed
!         do i=1,N_verts*3
!             if (any(residuals3(:,i) > 1.0e-8)) then
!                 test_failed = .true.
!                 exit
!             else 
!                 test_failed = .false.
!             end if
!         end do
!         if (test_failed) then
!             total_tests = total_tests + 1
!             write(m_char,'(I1)') m
!             failure_log(total_tests-passed_tests) = "dt_cross_n (edge "// trim(m_char) // ") test FAILED"
!             write(*,*) failure_log(total_tests-passed_tests)
!         else
!             write(*,'(A, I1, A)') "dt_cross_n (edge ",m,") test PASSED"
!             passed_tests = passed_tests + 1
!             total_tests = total_tests + 1
!         end if
!         test_failed = .false.
!         write(*,*) "" 
!         write(*,*) ""

!     ! end panel edge loop
!     end do



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST t_cross_dn !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) "---------------------------------- TEST t_cross_dn ----------------------------------"
!     write(*,*) ""
!     write(*,*) "the sensitivity of the edge vectors of panel 1 WRT each design variable"
!     write(*,*) ""
    
!     ! allocate central difference variables
    
!     allocate(t_cross_dn_FD(3,N_verts*3,3))
!     !deallocate(panels(index)%n_hat_g)


!     ! do for each edge
!     do m=1,3

!         !!!!!!!!! Finite Difference  t_cross_dn (edge m) !!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE t_cross_dn (edge ", m, ")"

!         ! we want the sensitivity of the edge vector of panel 1 edge 1 WRT X(beta)
!         index = 1

!         ! sensitivity to vertex 1
!         step = 0.0001

    
!         ! write results
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          t_cross_dn_FD panel 1 (edge ", m, ")"
!         write(*,*) "  t_cross_dn_x           t_cross_dn_y            t_cross_dn_z "
!         do i = 1, N_verts*3
!             t_cross_dn_FD(:,i, m) = cross(panels(index)%t_hat_g(:,m), d_n_g_FD(:,i))
!             write(*, '(3(f14.10, 4x))') t_cross_dn_FD(:,i, m)
!         end do 


!         !!!!!!!!!! ADJOINT d_t_cross_dn (edge m)!!!!!!!!!!!!!
!         write(*,*) ""
!         write(*,*) "------------------------------------------------"
!         write(*,*) ""
!         write(*,'(A, I1, A)') "  ADJOINT t_cross_dn (edge ", m, ")"
!         write(*,*) ""

!         ! ! only need to calc this once:
!         ! if (m==1) then

!         !     ! calculate d_n_hat_g (do it for all panels because they will be used later)
!         !     do i =1,N_panels
!         !         call panels(i)%calc_d_n_hat_g()
!         !     end do

!         ! end if

!         ! write sparse matrix
!         write(*,*) ""
!         write(*,'(A, I1, A)') "          t_cross_dn panel 1 (edge ", m, ")"
!         write(*,*) "  t_cross_dn_x           t_cross_dn_y           t_cross_dn_z             sparse_index       full_index"
!         do i=1,panels(index)%dt_cross_n(m)%sparse_num_cols
!             write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%t_cross_dn(m)%columns(i)%vector_values(:), &
!             i, panels(index)%t_cross_dn(m)%columns(i)%full_index
!         end do
!         write(*,*) ""

!         ! calculate residuals3
!         do i =1, N_verts*3
!             residuals3(:,i) = panels(index)%t_cross_dn(m)%get_values(i) - t_cross_dn_FD(:,i, m)
!         end do

!         write(*,'(A, I1, A)') "         t_cross_dn panel 1 (edge ", m, ") expanded "
!         write(*,*) "  t_cross_dn_x           t_cross_dn_y           t_cross_dn_z                            residuals"
!         do i = 1, N_verts*3
!             write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%t_cross_dn(m)%get_values(i), residuals3(:,i)
!         end do
!         write(*,*) ""


!         ! check if test failed
!         do i=1,N_verts*3
!             if (any(residuals3(:,i) > 1.0e-8)) then
!                 test_failed = .true.
!                 exit
!             else 
!                 test_failed = .false.
!             end if
!         end do
!         if (test_failed) then
!             total_tests = total_tests + 1
!             write(m_char,'(I1)') m
!             failure_log(total_tests-passed_tests) = "t_cross_dn (edge "// trim(m_char) // ") test FAILED"
!             write(*,*) failure_log(total_tests-passed_tests)
!         else
!             write(*,'(A, I1, A)') "t_cross_dn (edge ",m,") test PASSED"
!             passed_tests = passed_tests + 1
!             total_tests = total_tests + 1
!         end if
!         test_failed = .false.
!         write(*,*) "" 
!         write(*,*) ""

!     ! end panel edge loop
!     end do


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_v0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) "---------------------------------- TEST d_v0 ----------------------------------"
!     write(*,*) ""
!     write(*,*) "the sensitivity of the edge vectors of panel 1 WRT each design variable"
!     write(*,*) ""
    
!     ! allocate central difference variables
!     allocate(v0_up(N_verts*3))
!     allocate(v0_dn(N_verts*3))
!     allocate(d_v0_FD(3,N_verts*3))


!         !!!!!!!!! Finite Difference  d_v0 (row m) !!!!!!!!!
!     write(*,*) ""
!     write(*,*) "------------------------------------------------"
!     write(*,*) ""
!     write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_v0 (row ", m, ")"

!     ! we want the sensitivity of the v0 panel 1 row m WRT X(beta)
!     index = 1

!     ! sensitivity to vertex 1
!     step = 0.00001

!     ! for each x, y, z of v0 (row m) 
!     do k=1,3
!         ! do for each design variable
!         do i=1,3
!             do j=1,N_verts
!                 ! perturb up the current design variable
!                 vertices(j)%loc(i) = vertices(j)%loc(i) + step
                
!                 ! update panel edge outward normal unit vector calculations
!                 deallocate(panels)
!                 allocate(panels(N_panels))
!                 call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
!                 call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
!                 call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
!                 call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

!                 call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
!                 call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
!                 call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
!                 call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                
!                 call panels(index)%init_with_flow(freestream, .false., 0)

!                 ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                 v0_up(j + (i-1)*N_verts) = panels(index)%v0(k)

!                 ! perturb down the current design variable
!                 vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

            
!                 ! update panel edge outward normal unit vector calculations
!                 deallocate(panels)
!                 allocate(panels(N_panels))
!                 call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
!                 call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
!                 call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
!                 call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

!                 call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
!                 call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
!                 call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
!                 call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
                
!                 call panels(index)%init_with_flow(freestream, .false., 0)
                
!                 ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!                 v0_dn(j + (i-1)*N_verts) = panels(index)%v0(k)
            
!                 ! restore geometry
!                 vertices(j)%loc(i) = vertices(j)%loc(i) + step
                
                

!             end do 
!         end do 

!         ! central difference 
!         d_v0_FD(k,:) = (v0_up - v0_dn)/(2.*step)

!     end do

!     ! write results
!     write(*,*) ""
!     write(*,'(A, I1, A)') "          d_v0_FD panel 1 "
!     write(*,*) "  d_v0_x           d_v0_y            d_v0_z "
!     do i = 1, N_verts*3
!         write(*, '(3(f14.10, 4x))') d_v0_FD(:,i)
!     end do 


!     !!!!!!!!!! ADJOINT d_v0 (edge m)!!!!!!!!!!!!!
!     write(*,*) ""
!     write(*,*) "------------------------------------------------"
!     write(*,*) ""
!     write(*,'(A, I1, A)') "  ADJOINT d_v0 "
!     write(*,*) ""


!     ! Re init the body
!     deallocate(panels)
!     allocate(panels(N_panels))
!     call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
!     call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
!     call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
!     call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

!     call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
!     call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
!     call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
!     call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
    
!     call panels(1)%init_with_flow(freestream, .false., 0)

!     ! calculate d_A_g_to_ls 
!     call panels(index)%init_adjoint(freestream)

    
!     ! write sparse matrix
!     write(*,*) ""
!     write(*,'(A, I1, A)') "          d_v0 panel 1 "
!     write(*,*) "  d_v0_x           d_v0_y           d_v0_z             sparse_index       full_index"
!     do i=1,panels(index)%d_v0%sparse_num_cols
!         write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_v0%columns(i)%vector_values(:), &
!         i, panels(index)%d_v0%columns(i)%full_index
!     end do
!     write(*,*) ""

!     ! calculate residuals3
!     do i =1, N_verts*3
!         residuals3(:,i) = panels(index)%d_v0%get_values(i) - d_v0_FD(:,i)
!     end do

!     write(*,'(A, I1, A)') "         d_v0 panel 1  expanded "
!     write(*,*) "  d_v0_x           d_v0_y           d_v0_z                            residuals"
!     do i = 1, N_verts*3
!         write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_v0%get_values(i), residuals3(:,i)
!     end do
!     write(*,*) ""


!     ! check if test failed
!     do i=1,N_verts*3
!         if (any(residuals3(:,i) > 1.0e-8)) then
!             test_failed = .true.
!             exit
!         else 
!             test_failed = .false.
!         end if
!     end do
!     if (test_failed) then
!         total_tests = total_tests + 1
!         write(m_char,'(I1)') m
!         failure_log(total_tests-passed_tests) = "d_v0  test FAILED"
!         write(*,*) failure_log(total_tests-passed_tests)
!     else
!         write(*,'(A, I1, A)') "d_v0  test PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) "" 
!     write(*,*) ""



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_u0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(*,*) "---------------------------------- TEST d_u0 ----------------------------------"
! write(*,*) ""
! write(*,*) "the sensitivity of the edge vectors of panel 1 WRT each design variable"
! write(*,*) ""

! ! allocate central difference variables
! allocate(u0_up(N_verts*3))
! allocate(u0_dn(N_verts*3))
! allocate(d_u0_FD(3,N_verts*3))


! !!!!!!!!! Finite Difference  d_u0  !!!!!!!!!
! write(*,*) ""
! write(*,*) "------------------------------------------------"
! write(*,*) ""
! write(*,'(A, I1, A)') "  CENTRAL DIFFERENCE d_u0 "

! ! we want the sensitivity of the u0 panel 1 row m WRT X(beta)
! index = 1

! ! sensitivity to vertex 1
! step = 0.00001

! ! for each x, y, z of v0 (row m) 
! do k=1,3
!     ! do for each design variable
!     do i=1,3
!         do j=1,N_verts
!             ! perturb up the current design variable
!             vertices(j)%loc(i) = vertices(j)%loc(i) + step
            
!             ! update panel edge outward normal unit vector calculations
!             deallocate(panels)
!             allocate(panels(N_panels))
!             call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
!             call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
!             call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
!             call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

!             call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
!             call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
!             call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
!             call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
            
!             call panels(index)%init_with_flow(freestream, .false., 0)

!             ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!             u0_up(j + (i-1)*N_verts) = panels(index)%u0(k)

!             ! perturb down the current design variable
!             vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

        
!             ! update panel edge outward normal unit vector calculations
!             deallocate(panels)
!             allocate(panels(N_panels))
!             call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
!             call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
!             call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
!             call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

!             call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
!             call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
!             call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
!             call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward
            
!             call panels(index)%init_with_flow(freestream, .false., 0)
            
!             ! put the x y or z component of the panel's perturbed edge outward normal unit vector in a list
!             u0_dn(j + (i-1)*N_verts) = panels(index)%u0(k)
        
!             ! restore geometry
!             vertices(j)%loc(i) = vertices(j)%loc(i) + step
            
            

!         end do 
!     end do 

!     ! central difference 
!     d_u0_FD(k,:) = (u0_up - u0_dn)/(2.*step)

! end do

! ! write results
! write(*,*) ""
! write(*,'(A, I1, A)') "          d_u0_FD panel 1 "
! write(*,*) "  d_u0_x           d_u0_y            d_u0_z "
! do i = 1, N_verts*3
!     write(*, '(3(f14.10, 4x))') d_u0_FD(:,i)
! end do 


! !!!!!!!!!! ADJOINT d_u0 (edge m)!!!!!!!!!!!!!
! write(*,*) ""
! write(*,*) "------------------------------------------------"
! write(*,*) ""
! write(*,'(A, I1, A)') "  ADJOINT d_u0 "
! write(*,*) ""


! ! Re init the body
! deallocate(panels)
! allocate(panels(N_panels))
! call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
! call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
! call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
! call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

! call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
! call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
! call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
! call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward

! call panels(1)%init_with_flow(freestream, .false., 0)

! ! calculate d_A_g_to_ls 
! call panels(index)%init_adjoint(freestream)


! ! write sparse matrix
! write(*,*) ""
! write(*,'(A, I1, A)') "          d_u0 panel 1 "
! write(*,*) "  d_u0_x           d_u0_y           d_u0_z             sparse_index       full_index"
! do i=1,panels(index)%d_u0%sparse_num_cols
!     write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') panels(index)%d_u0%columns(i)%vector_values(:), &
!     i, panels(index)%d_u0%columns(i)%full_index
! end do
! write(*,*) ""

! ! calculate residuals3
! do i =1, N_verts*3
!     residuals3(:,i) = panels(index)%d_u0%get_values(i) - d_u0_FD(:,i)
! end do

! write(*,'(A, I1, A)') "         d_u0 panel 1  expanded "
! write(*,*) "  d_u0_x           d_u0_y           d_u0_z                            residuals"
! do i = 1, N_verts*3
!     write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') panels(index)%d_u0%get_values(i), residuals3(:,i)
! end do
! write(*,*) ""


! ! check if test failed
! do i=1,N_verts*3
!     if (any(residuals3(:,i) > 1.0e-8)) then
!         test_failed = .true.
!         exit
!     else 
!         test_failed = .false.
!     end if
! end do
! if (test_failed) then
!     total_tests = total_tests + 1
!     write(m_char,'(I1)') m
!     failure_log(total_tests-passed_tests) = "d_u0  test FAILED"
!     write(*,*) failure_log(total_tests-passed_tests)
! else
!     write(*,'(A, I1, A)') "d_u0  test PASSED"
!     passed_tests = passed_tests + 1
!     total_tests = total_tests + 1
! end if
! test_failed = .false.
! write(*,*) "" 
! write(*,*) ""

end program gradient_test