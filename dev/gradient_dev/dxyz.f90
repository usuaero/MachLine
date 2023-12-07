program dxyz

    use base_geom_mod
    use panel_mod
    use flow_mod
    
    implicit none

    type(panel) :: panel_1, panel_2, panel_3, panel_4, panel_5, panel_6
    type(flow) :: freestream
    type(eval_point_geom) :: geom
    type(integrals) :: int
    real,dimension(3) :: v1, v2, v3, v4, v5, v6
    real,dimension(3,6) :: points
    real,dimension(:),allocatable :: X_beta
    integer :: i, N_verts, N_panels
    real,dimension(:,:),allocatable :: vertex_locs
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),allocatable,dimension(:) :: panels    ! list of panels, this should be a mesh attribute
    character(len=:),allocatable :: spanwise_axis

    ! For this simple example, we don't need to collapse dublicates
    N_verts = 6
    N_panels = 8

    ! Initialize panels
    v1 = (/ 0.0, 0.0,-0.2/)  ! top
    v2 = (/ 1.0, 0.0, 0.0/)  ! forward 
    v3 = (/ 0.0, 1.0, 0.0/)  ! right wingtip
    v4 = (/-1.0, 0.0, 0.0/)  ! aft
    v5 = (/ 0.0,-1.0, 0.0/)  ! left wingtip
    v6 = (/ 0.0, 0.0, 0.2/)  ! bottom
    
    ! put vertex_locs in a list
    allocate(vertex_locs(3,N_verts))
    do i =1, N_verts
        vertex_locs(:,i)
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
            X_beta(j + (i-1)*N_verts) = vertex_locs(i,j)
        end do
    end do

    ! initialize 8 panels
    allocate(panels(N_panels))
    call panels(i)%init(vertices(i1), vertices(i2), vertices(i3), i)
    call panel_1%init(vert1, vert2, vert3, 1, .false.) !    top, right, forward
    call panel_2%init(vert1, vert3, vert4, 2, .false.) !    top, right,     aft 
    call panel_3%init(vert1, vert4, vert5, 3, .false.) !    top,  left,     aft
    call panel_4%init(vert1, vert5, vert2, 4, .false.) !    top,  left, forward

    call panel_5%init(vert6, vert2, vert3, 5, .false.) ! bottom, right, forward
    call panel_6%init(vert6, vert3, vert4, 6, .false.) ! bottom, right,     aft
    call panel_7%init(vert6, vert4, vert5, 7, .false.) ! bottom,  left,     aft
    call panel_8%init(vert6, vert5, vert2, 8, .false.) ! bottom,  left, forward





    
    ! Initialize points
    points(:,1) = (/0., 0., 1./)
    points(:,2) = (/0., 0., 2./)
    points(:,3) = (/1., 1., 1./)
    points(:,4) = (/2., 2., 2./)
    points(:,5) = (/0., 0., -1./)
    points(:,6) = (/0.5, 0.5, 0./)
    
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


    ! Init with flow
    call my_panel%init_with_flow(freestream, .false., 0)

    ! Loop through points
    do i=1,6

        ! Get geometry
        geom = my_panel%calc_subsonic_geom(points(:,i), freestream, .false.)

        ! Get integrals
        allocate(int%F111(3), source=0.)
        call my_panel%calc_subsonic_edge_integrals(geom, freestream, .false., int)
        call my_panel%calc_subsonic_panel_integrals(geom, freestream, .false., int)
        deallocate(int%F111)
        write(*,*) int%hH113/(4.*pi)
    end do
    
end program dxyz