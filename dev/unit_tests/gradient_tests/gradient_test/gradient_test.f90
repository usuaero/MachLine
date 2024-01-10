program gradient_test

    ! tests various intermediate sensitivities 
    use adjoint_mod
    use base_geom_mod
    use panel_mod
    use flow_mod
    
    implicit none

    
    type(flow) :: freestream
    type(eval_point_geom) :: geom
    type(integrals) :: int
    real :: step
    real,dimension(:,:),allocatable :: v, vertex_locs, d_loc_FD, d_centr_FD
    real,dimension(3) :: values1, values2, values3
    real,dimension(:),allocatable :: X_beta, loc_up, loc_dn, centr_up, centr_dn
    integer :: i,j,k, N_verts, N_panels, vert
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels    ! list of panels, this should be a mesh attribute
    character(len=:),allocatable :: spanwise_axis

    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(20) :: failure_log


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

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST d_loc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "---------------------------------- TEST d_loc -----------------------------------"
    write(*,*) ""
    
    !!!!!!!!!! ADJOINT d_loc!!!!!!!!!!!!!
    write(*,*) "ADJOINT d_loc"
    write(*,*) ""
    ! init vertex attribute d_loc
    do i=1,N_verts
                
        ! call vertex init
        call vertices(i)%init_adjoint(N_verts)
        
    end do

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_loc point 1"
    write(*,*) "     d_x1        d_y1       d_z1          sparse_index      full_index"
    do i=1,vertices(1)%d_loc%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') vertices(1)%d_loc%columns(i)%vector_values(1), &
        vertices(1)%d_loc%columns(i)%vector_values(2),vertices(1)%d_loc%columns(i)%vector_values(3), &
        i, vertices(1)%d_loc%columns(i)%full_index
    end do
    write(*,*) ""

    vert = 1
    write(*,*) "         d_loc vertex 1 full "
    write(*,*) "    d_x1              d_y1            d_z1"
    do i = 1, N_verts*3
        write(*, '(3(f14.6, 2x))') vertices(vert)%d_loc%get_values(i)
    end do
    write(*,*) ""


    !!!!!!!!! Finite Difference  d_loc !!!!!!!!!
    write(*,*) ""
    write(*,*) "CENTRAL DIFFERENCE d_loc"

    ! sensitivity to vertex 1
    step = 0.0001

    ! perturb x1 up
    allocate(loc_up(N_verts*3))
    allocate(loc_dn(N_verts*3))
    allocate(d_loc_FD(N_verts*3,3))
    ! for each x, y, z of centr 1 
    do k=1,3

        do i=1,3
            do j=1,N_verts

                ! perturb up the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

                ! put the x y or z component of the vertex of interest (vert) in a list
                loc_up(j + (i-1)*N_verts) = vertices(vert)%loc(k)

                ! perturb down the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                ! put the x y or z component of the vertex of interest (vert) in a list
                loc_dn(j + (i-1)*N_verts) = vertices(vert)%loc(k)
                
                ! central difference 
                d_loc_FD(:,k) = (loc_up - loc_dn)/(2.*step)
            
                ! restore geometry
                vertices(j)%loc(i) = vertices(j)%loc(i) + step
            end do 
        end do 
    end do
    
    
    write(*,*) ""
    write(*,*) "                d_loc_FD vertex 1"
    write(*,*) "    d_loc_x              d_loc_y            d_loc_z"
    do i = 1, N_verts*3
        write(*, '(3(f14.6, 2x))') d_loc_FD(i,1), d_loc_FD(i,2), d_loc_FD(i,3)
    end do 


    ! check if test failed
    do i=1,vertices(vert)%d_loc%full_num_cols
        if (any(vertices(vert)%d_loc%get_values(i) - d_loc_FD(i,:) > 1.0e-12)) then
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
    
    !!!!!!!!!! ADJOINT d_centr!!!!!!!!!!!!!
    write(*,*) "ADJOINT d_centr"
    write(*,*) ""

    ! for each panel,
    do i=1,N_panels
                
        ! calculate d_centr
        call panels(i)%calc_d_centr(N_verts)
        
    end do

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         d_centr panel 1"
    write(*,*) "     d_centr_x        d_centr_y       d_centr_z          sparse_index      full_index"
    do i=1,panels(1)%d_centr%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') panels(1)%d_centr%columns(i)%vector_values(1), &
        panels(1)%d_centr%columns(i)%vector_values(2),panels(1)%d_centr%columns(i)%vector_values(3), &
        i, panels(1)%d_centr%columns(i)%full_index
    end do
    write(*,*) ""

    vert = 1
    write(*,*) "         d_centr panel 1 full "
    write(*,*) "    d_centr_x              d_centr_y            d_centr_z"
    do i = 1, N_verts*3
        write(*, '(3(f14.6, 2x))') panels(vert)%d_centr%get_values(i)
    end do
    write(*,*) ""


    !!!!!!!!! Finite Difference  d_centr !!!!!!!!!
    write(*,*) ""
    write(*,*) "CENTRAL DIFFERENCE d_centr"

    ! sensitivity to vertex 1
    step = 0.0001

    ! perturb x1 up
    allocate(centr_up(N_verts*3))
    allocate(centr_dn(N_verts*3))
    allocate(d_centr_FD(N_verts*3,3))

    ! for each x, y, z of centr 1 
    do k=1,3

        do i=1,3
            do j=1,N_verts
                ! perturb up the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

                ! update panel centroid calculations
                call panels(vert)%calc_centroid()

                ! put the x y or z component of the panel's perturbed centroid in a list
                centr_up(j + (i-1)*N_verts) = panels(vert)%centr(k)

                ! perturb down the current design variable
                vertices(j)%loc(i) = vertices(j)%loc(i) - 2.*step

                ! update panel centroid calculations
                call panels(vert)%calc_centroid()
                
                ! put the x y or z component of the panel's perturbed centroid in a list
                centr_dn(j + (i-1)*N_verts) = panels(vert)%centr(k)
                
                ! central difference 
                d_centr_FD(:,k) = (centr_up - centr_dn)/(2.*step)
            
                ! restore geometry
                vertices(j)%loc(i) = vertices(j)%loc(i) + step

            end do 
        end do 
    end do
    write(*,*) ""
    write(*,*) "                d_centr_FD panel 1"
    write(*,*) "    d_centr_x          d_centr_y          d_centr_z"
    do i = 1, N_verts*3
        write(*, '(3(f14.6, 2x))') d_centr_FD(i,1), d_centr_FD(i,2), d_centr_FD(i,3)
    end do 

    ! check if test failed
    do i=1,panels(vert)%d_centr%full_num_cols
        if (any(panels(vert)%d_centr%get_values(i) - d_centr_FD(i,:) > 1.0e-12)) then
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







    !!!!!!!!!!!!!! GRADIENT TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------GRADIENT TEST RESULTS--------------"
    write(*,*) ""
    write(*,'(I15,a14)') total_tests - passed_tests, " tests FAILED"
    write(*,*) ""
    write(*,'(I4,a9,I2,a14)') passed_tests, " out of ", total_tests, " tests PASSED"
    if (passed_tests < total_tests)then
        write(*,*) ""
        write(*,*) "Failure Log:"
        do i=1,total_tests-passed_tests
            write(*,*) failure_log(i)
        end do
    end if
    
    write(*,*) ""
    write(*,*) "Program Complete"
    write(*,*) ""





end program gradient_test