program d_loc_complex
    ! Finds dx/dX(beta), dy/dX(beta), dz/dX(beta) for each vertex in the mesh
    use base_geom_mod
    !use panel_mod
    
    implicit none

    
    
    complex :: step
    complex,dimension(:,:),allocatable :: v, vertex_locs, d_loc_CX
    complex,dimension(:),allocatable :: X_beta, X_beta_CX
    integer :: i,j,k, N_verts, N_panels, vert
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    !type(panel),dimension(:),allocatable :: panels    ! list of panels, this should be a mesh attribute
 

    


    ! For this simple example, we don't need to collapse dublicates
    N_verts = 6
    N_panels = 8
    
    allocate(v(3,N_verts))
    ! Initialize vertices
    v(:,1) = (/ (0.0, 0.0), (0.0, 0.0),(-0.2, 0.0)/)  ! top
    v(:,2) = (/ (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)  ! forward 
    v(:,3) = (/ (0.0, 0.0), (1.0, 0.0), (0.0, 0.0)/)  ! right wingtip
    v(:,4) = (/(-1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)  ! aft
    v(:,5) = (/ (0.0, 0.0),(-1.0, 0.0), (0.0, 0.0)/)  ! left wingtip
    v(:,6) = (/ (0.0, 0.0), (0.0, 0.0), (0.2, 0.0)/)  ! bottom
    
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
    !allocate(panels(N_panels))
    
    !call panels(1)%init(vertices(1), vertices(2), vertices(3), 1, .false.) !    top, right, forward
    !call panels(2)%init(vertices(1), vertices(3), vertices(4), 2, .false.) !    top, right,     aft 
    !call panels(3)%init(vertices(1), vertices(4), vertices(5), 3, .false.) !    top,  left,     aft
    !call panels(4)%init(vertices(1), vertices(5), vertices(2), 4, .false.) !    top,  left, forward

    !call panels(5)%init(vertices(6), vertices(2), vertices(3), 5, .false.) ! bottom, right, forward
    !call panels(6)%init(vertices(6), vertices(3), vertices(4), 6, .false.) ! bottom, right,     aft
    !call panels(7)%init(vertices(6), vertices(4), vertices(5), 7, .false.) ! bottom,  left,     aft
    !call panels(8)%init(vertices(6), vertices(5), vertices(2), 8, .false.) ! bottom,  left, forward


    !!!!!!!!! Complex Step !!!!!!!!!

    ! sensitivity to vertex 1, x1
    step = (0.0, 1.0e-200)
    ! perturb x1 up
    allocate(X_beta_CX(N_verts*3))
    allocate(d_loc_CX(N_verts*3,3))
    do k=1,3
        vertices(vert)%loc(k) = vertices(vert)%loc(k) + step
        ! get perturbation up
        do i=1,3
            do j=1,N_verts
                X_beta_CX(j + (i-1)*N_verts) = vertices(j)%loc(i)
            end do
        end do
        
        ! restore geometry
        vertices(vert)%loc(k) = vertices(vert)%loc(k) - step
        

        d_loc_CX(:,k) = (X_beta_CX)/(step)
    end do 

    write(*,*) "vertex", vert, "d_loc_CX ="
    do i = 1, N_verts*3
        write(*, '(3(f14.10, 2x))') AIMAG(d_loc_CX(i,1)), AIMAG(d_loc_CX(i,2)), AIMAG(d_loc_CX(i,3))
    end do 

    !!!!!!!!! End Complex-Step !!!!!!!!!




    ! Init with flow
    !call my_panel%init_with_flow(freestream, .false., 0)

end program d_loc_complex