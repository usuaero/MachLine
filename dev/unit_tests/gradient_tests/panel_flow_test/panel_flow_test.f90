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

    real,dimension(:),allocatable :: residuals, X_beta, A_g_to_ls_up, A_g_to_ls_dn,  A_ls_to_g_up, &
    A_ls_to_g_dn, vertices_ls_up, vertices_ls_dn, d_vertices_ls_FD, n_hat_ls_up, n_hat_ls_dn, &
    d_n_hat_ls_FD, T_mu_up, T_mu_dn

    real,dimension(:,:),allocatable :: v, vertex_locs, residuals3 

    real,dimension(:,:,:),allocatable ::  d_A_g_to_ls_FD, d_A_ls_to_g_FD, d_T_mu_FD

    integer :: i,j,k,m,n, N_verts, N_panels, vert, index
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels    ! list of panels, this should be a mesh attribute
    character(len=:),allocatable :: spanwise_axis

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
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))


    ! run adjoint once

    ! call vertex init (do it for all vertices because they will be used later)
    do i =1,N_verts
        call vertices(i)%init_adjoint(N_verts)
    end do
    
    allocate(adjoint_panels(N_panels))
    
    call adjoint_panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
    call adjoint_panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
    call adjoint_panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
    call adjoint_panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward
    
    call adjoint_panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
    call adjoint_panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
    call adjoint_panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
    call adjoint_panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward

    ! following sensitivities are with respect to a perturbation in x y and z of vertex 1
    index = 1

    call adjoint_panels(index)%init_with_flow(freestream, .false., 0)
    
    call adjoint_panels(index)%init_adjoint()

    call adjoint_panels(index)%init_with_flow_adjoint(freestream)




!!!!!!!!!!!!!!!!!!!!!!!!!!FLOW-DEPENDENT PANEL SENSITIVITIES TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "---------------------------------- FLOW-DEPENDENT PANEL SENSITIVITIES TEST ------------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""

    
    
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
           
        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_g_to_ls panel 1 (row ", m, ")"
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y           d_A_g_to_ls_z             sparse_index       full_index"
        do i=1,adjoint_panels(index)%d_A_g_to_ls(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_panels(index)%d_A_g_to_ls(m)%columns(i)%vector_values(:), &
            i, adjoint_panels(index)%d_A_g_to_ls(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_panels(index)%d_A_g_to_ls(m)%get_values(i) - d_A_g_to_ls_FD(m,i,:)
        end do
        
        write(*,'(A, I1, A)') "         d_A_g_to_ls panel 1 (row ", m, ") expanded "
        write(*,*) "  d_A_g_to_ls_x           d_A_g_to_ls_y           d_A_g_to_ls_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_panels(index)%d_A_g_to_ls(m)%get_values(i), residuals3(:,i)
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


        
        ! write sparse matrix
        write(*,*) ""
        write(*,'(A, I1, A)') "          d_A_ls_to_g panel 1 (row ", m, ")"
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y           d_A_ls_to_g_z             sparse_index       full_index"
        do i=1,adjoint_panels(index)%d_A_ls_to_g(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_panels(index)%d_A_ls_to_g(m)%columns(i)%vector_values(:), &
            i, adjoint_panels(index)%d_A_ls_to_g(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_panels(index)%d_A_ls_to_g(m)%get_values(i) - d_A_ls_to_g_FD(m,i,:)
        end do

        write(*,'(A, I1, A)') "         d_A_ls_to_g panel 1 (row ", m, ") expanded "
        write(*,*) "  d_A_ls_to_g_x           d_A_ls_to_g_y           d_A_ls_to_g_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_panels(index)%d_A_ls_to_g(m)%get_values(i), residuals3(:,i)
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

            
            ! write sparse matrix
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex ", m, " eta coordinate)"
            end if
            write(*,*) "  sparse value                  sparse_index       full_index"
            do i=1,adjoint_panels(index)%d_vertices_ls(n,m)%sparse_size
                write(*,'(f14.10, 20x, I5, 12x, I5)') adjoint_panels(index)%d_vertices_ls(n,m)%elements(i)%value, &
                i, adjoint_panels(index)%d_vertices_ls(n,m)%elements(i)%full_index
            end do
            write(*,*) ""

            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = adjoint_panels(index)%d_vertices_ls(n,m)%get_value(i) - d_vertices_ls_FD(i)
            end do

            if (n==1) then
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex", m, " xi coordinate) expanded"
            else 
                write(*,'(A, I1, A)') "  d_vertices_ls panel 1 (vertex ", m, " eta coordinate) expanded"
            end if
            write(*,*) "  adjoint value         residuals"
            do i = 1, N_verts*3
                write(*, '(f14.10,3x, f14.10)') adjoint_panels(index)%d_vertices_ls(n,m)%get_value(i), residuals(i)
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

            
            ! write sparse matrix
            write(*,*) ""
            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge", m, " xi coordinate)"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge ", m, " eta coordinate)"
            end if
            write(*,*) "  sparse value              sparse_index       full_index"
            do i=1,adjoint_panels(index)%d_n_hat_ls(n,m)%sparse_size
                write(*,'(f14.10, 20x, I5, 12x, I5)') adjoint_panels(index)%d_n_hat_ls(n,m)%elements(i)%value, &
                i, adjoint_panels(index)%d_n_hat_ls(n,m)%elements(i)%full_index
            end do
            write(*,*) ""

            ! calculate residuals3
            do i =1, N_verts*3
                residuals(i) = adjoint_panels(index)%d_n_hat_ls(n,m)%get_value(i) - d_n_hat_ls_FD(i)
            end do

            if (n==1) then
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge", m, " xi coordinate) expanded"
            else 
                write(*,'(A, I1, A)') "  d_n_hat_ls panel 1 (edge ", m, " eta coordinate) expanded"
            end if
            write(*,*) "  adjoint value         residuals"
            do i = 1, N_verts*3
                write(*, '(f14.10,3x, f14.10)') adjoint_panels(index)%d_n_hat_ls(n,m)%get_value(i), residuals(i)
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

        ! we want the sensitivity of the T_mu panel 1 row m WRT X(beta)
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
                    T_mu_up(j + (i-1)*N_verts) = panels(index)%T_mu(m,k)

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
                    T_mu_dn(j + (i-1)*N_verts) = panels(index)%T_mu(m,k)
                
                    ! restore geometry
                    vertices(j)%loc(i) = vertices(j)%loc(i) + step
                    
                    

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
        do i=1,adjoint_panels(index)%d_T_mu(m)%sparse_num_cols
            write(*,'(3(f14.10, 4x), 12x, I5, 12x, I5)') adjoint_panels(index)%d_T_mu(m)%columns(i)%vector_values(:), &
            i, adjoint_panels(index)%d_T_mu(m)%columns(i)%full_index
        end do
        write(*,*) ""

        ! calculate residuals3
        do i =1, N_verts*3
            residuals3(:,i) = adjoint_panels(index)%d_T_mu(m)%get_values(i) - d_T_mu_FD(m,i,:)
        end do

        write(*,'(A, I1, A)') "         d_T_mu panel 1 (row ", m, ") expanded "
        write(*,*) "  d_T_mu_x           d_T_mu_y           d_T_mu_z                            residuals"
        do i = 1, N_verts*3
            write(*, '(3(f14.10, 4x),3x, 3(f14.10, 4x))') adjoint_panels(index)%d_T_mu(m)%get_values(i), residuals3(:,i)
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




end program gradient_test