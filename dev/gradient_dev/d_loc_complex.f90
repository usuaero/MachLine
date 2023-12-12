program d_loc_complex
    ! Finds dx/dX(beta), dy/dX(beta), dz/dX(beta) for each vertex in the mesh
    use base_geom_mod
    use panel_mod
    use flow_mod
    
    implicit none

    
    type(flow) :: freestream
    type(eval_point_geom) :: geom
    type(integrals) :: int
    complex :: step
    real,dimension(:,:),allocatable :: v, vertex_locs, d_loc_FD
    real,dimension(:),allocatable :: X_beta, X_beta_up, X_beta_dn
    integer :: i,j,k, N_verts, N_panels, vert
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels    ! list of panels, this should be a mesh attribute
    character(len=:),allocatable :: spanwise_axis

    
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

    
    !!!!!!!!!! Adjoint Method !!!!!!!!!
    ! allocate vertex attribute d_loc
    do i=1,N_verts
        allocate(vertices(i)%d_loc(N_verts*3,3), source=0.)
        ! use logic to put a 1 at the corresponding X(beta) spot 
        do j=1,3
            vertices(i)%d_loc(i + (j-1)*N_verts,j) = 1.
        end do
    end do

    vert = 2
    write(*,*) "vertex ", vert, " d_loc ="
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 2x))') vertices(vert)%d_loc(i,1), vertices(vert)%d_loc(i,2), vertices(vert)%d_loc(i,3)
    end do 
     !!!!!!!!!! End Adjoint Method !!!!!!!!!


    !!!!!!!!! Finite Difference !!!!!!!!!
    ! sensitivity to vertex 1, x1
    step = 0.0001
    ! perturb x1 up
    allocate(X_beta_up(N_verts*3))
    allocate(X_beta_dn(N_verts*3))
    allocate(d_loc_FD(N_verts*3,3))
    do k=1,3
        vertices(vert)%loc(k) = vertices(vert)%loc(k) + step
        ! get perturbation up
        do i=1,3
            do j=1,N_verts
                X_beta_up(j + (i-1)*N_verts) = vertices(j)%loc(i)
            end do
        end do
        
        ! get perturbation down
        vertices(vert)%loc(k) = vertices(vert)%loc(k) - 2.*step
        do i=1,3
            do j=1,N_verts
                X_beta_dn(j + (i-1)*N_verts) = vertices(j)%loc(i)
            end do
        end do
        
        ! restore geometry
        vertices(vert)%loc(k) = vertices(vert)%loc(k) + step
        

        d_loc_FD(:,k) = (X_beta_up - X_beta_dn)/(2.*step)
    end do 

    write(*,*) "vertex", vert, "d_loc_FD ="
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 2x))') d_loc_FD(i,1), d_loc_FD(i,2), d_loc_FD(i,3)
    end do 

    !!!!!!!!! End Finite Difference !!!!!!!!!




    ! Init with flow
    !call my_panel%init_with_flow(freestream, .false., 0)

end program d_loc_complex