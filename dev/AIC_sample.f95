program AIC_sample

    use base_geom_mod
    use panel_mod
    use flow_mod
    
    implicit none

    type(panel) :: my_panel
    type(vertex) :: vert1, vert2, vert3
    type(flow) :: freestream
    type(eval_point_geom) :: geom
    type(integrals) :: int
    real,dimension(3) :: v1, v2, v3
    real,dimension(3,6) :: points
    integer :: i
    character(len=:),allocatable :: spanwise_axis

    ! Initialize panel
    v1 = (/ 0.0, 0.0, 0.0/)
    v2 = (/ 1.0, 0.0, 0.0/)
    v3 = (/ 0.0, 1.0, 0.0/)
    call vert1%init(v1, 0, 1)
    call vert2%init(v2, 0, 1)
    call vert3%init(v3, 0, 1)
    call my_panel%init(vert1, vert2, vert3, 0)
    spanwise_axis = 'y+'

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
    freestream%c_hat_g = (/1., 0., 0./)
    freestream%U = 1.
    freestream%U_inv = 1.
    freestream%M_inf = 0.
    freestream%B = 1.
    freestream%s = 1.
    freestream%K = 4.*pi
    freestream%K_inv = 1./(4.*pi)
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
    
end program AIC_sample