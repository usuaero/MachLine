program test28

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
    use vtk_mod
    
    implicit none

    !!!!!!!!!!!!!!!!!!! STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!
    character(100) :: test_input, adjoint_input
    character(len=:),allocatable :: body_file, wake_file, control_point_file, points_file, &
    mirrored_body_file, points_output_file
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
    type(dod) :: test_dod_info, adjoint_dod_info
    type(integrals) :: test_int, adjoint_int
    integer :: i_unit
    logical :: exists, found 
    integer :: adjoint_solver_stat, test_solver_stat
    type(sparse_vector) :: zeros

    ! real,dimension(3) :: adjoint_P, test_P, test_v_d, test_v_s
    ! type(sparse_matrix) :: adjoint_d_P_term2
    ! type(sparse_matrix) :: adjoint_d_P
    ! type(sparse_matrix) :: adjoint_d_v_d_panel

    !!!!!!!!!!!!!!!!!!!!! END STUFF FROM MAIN !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! TESTING STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(:),allocatable :: residuals
    real,dimension(:,:),allocatable ::  residuals3, CF_up, CF_dn, d_CF_FD

    integer :: i,j,k,m,n,p, N_verts, N_panels, vert, index, cp_ind, row,col, stat
    real :: step, error_allowed, cp_offset
    type(vertex),dimension(:),allocatable :: vertices ! list of vertex types, this should be a mesh attribute
    type(panel),dimension(:),allocatable :: panels, adjoint_panels   ! list of panels, this should be a mesh attribute
    ! test stuff
    integer :: passed_tests, total_tests
    logical :: test_failed
    character(len=100),dimension(50) :: failure_log
    character(len=10) :: m_char
    integer(8) :: start_count, end_count
    real(16) :: count_rate, time

    type(vtk_out) :: body_vtk
    integer ::  N_cells
    real,dimension(:),allocatable :: panel_inclinations, orders, N_discont_edges, convex
    real,dimension(:,:),allocatable :: cents
    real,dimension(:,:),allocatable :: vertex_normals, d_CF_x, d_CF_y, d_CF_z, freestream_vector

    !!!!!!!!!!!!!!!!!!! END TESTING STUFF !!!!!!!!!!!!!!!!!!!!!11

    test_failed = .true. ! assume test failed, if the test condition is met, test passed
    ! NOTE: on the sparse vector test, I assume the test passes, if it fails a test condition, test fails
    passed_tests = 0
    total_tests = 0

    error_allowed = 1.0e-2
    
    step = 0.001
    index = 1
    cp_ind = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                             FROM MAIN

    !!!!!!!!!!!!!!! TEST INPUT (calc_adjoint = false) !!!!!!!!!!!!!!!!!!!!!!!
    ! Set up run
    call json_initialize()

    test_input = "dev/input_files/adjoint_inputs/central_difference_test.json"
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
    call test_solver%solve(test_mesh, test_solver_stat, formulation,freestream_flow)


    !!!!!!! write to vtk !!!!!!!!!!!!!!!!!!!
    
    call delete_file(body_file)
        
    ! Determine number of cells and verts to export
    N_cells = test_mesh%N_panels
    N_verts = test_mesh%N_verts
    
    ! Get panel inclinations, centroids, and distribution orders
    allocate(panel_inclinations(test_mesh%N_panels))
    allocate(orders(test_mesh%N_panels))
    allocate(N_discont_edges(test_mesh%N_panels))
    allocate(cents(3,test_mesh%N_panels))
    allocate(convex(test_mesh%N_verts), source=0.)
    do i=1,test_mesh%N_panels
        panel_inclinations(i) = test_mesh%panels(i)%r
        orders(i) = test_mesh%panels(i)%order
        N_discont_edges(i) = test_mesh%panels(i)%N_discont_edges
        cents(:,i) = test_mesh%panels(i)%centr
    end do
    do i=1,test_mesh%N_verts
        if (test_mesh%vertices(i)%convex) convex(i) = 1.
    end do

    ! Write geometry
    call body_vtk%begin(body_file)
    call body_vtk%write_points(test_mesh%vertices)
    call body_vtk%write_panels(test_mesh%panels, mirror=.false.)
    call body_vtk%write_cell_normals(test_mesh%panels)
    call body_vtk%write_cell_scalars(panel_inclinations, "inclination")
    call body_vtk%write_cell_scalars(orders, "distribution_order")
    call body_vtk%write_cell_scalars(N_discont_edges, "N_discontinuous_edges")
    call body_vtk%write_cell_vectors(cents, "centroid")


    ! Pressures
    if (allocated(test_mesh%C_p_inc)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_inc(1:N_cells), "C_p_inc")
    end if
    if (allocated(test_mesh%C_p_ise)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_ise(1:N_cells), "C_p_ise")
    end if
    if (allocated(test_mesh%C_p_2nd)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_2nd(1:N_cells), "C_p_2nd")
    end if
    if (allocated(test_mesh%C_p_lin)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_lin(1:N_cells), "C_p_lin")
    end if
    if (allocated(test_mesh%C_p_sln)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_sln(1:N_cells), "C_p_sln")
    end if

    ! Corrected pressures
    if (allocated(test_mesh%C_p_pg)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_pg(1:N_cells), "C_p_PG")
    end if
    if (allocated(test_mesh%C_p_kt)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_kt(1:N_cells), "C_p_KT")
    end if
    if (allocated(test_mesh%C_p_lai)) then
        call body_vtk%write_cell_scalars(test_mesh%C_p_lai(1:N_cells), "C_p_L")
    end if

    ! Source strengths
    call body_vtk%write_cell_scalars(test_mesh%sigma(1:test_mesh%N_panels), "sigma")

    ! Cell velocities and sources
    call body_vtk%write_cell_vectors(test_mesh%V_cells(:,1:N_cells), "v")
    call body_vtk%write_cell_vectors(test_mesh%V_cells_inner(:,1:N_cells), "v_inner")
    call body_vtk%write_cell_vectors(test_mesh%dC_f(:,1:N_cells), "dC_f")

    ! Surface potential values
    call body_vtk%write_point_scalars(test_mesh%mu(1:test_mesh%N_verts), "mu")
    call body_vtk%write_point_scalars(test_mesh%Phi_u(1:test_mesh%N_verts), "Phi_u")
        
    
    
    !!!!!!!!!!!!!!!!!!!!! END TEST MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    call system_clock(start_count, count_rate)
   
    
    
    N_verts = test_mesh%N_verts
    N_panels = test_mesh%N_panels
    
    
    allocate(residuals3(3,N_verts*3))
    allocate(residuals(N_verts*3))

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! d_CF_wrt_vars_test !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) ""
    write(*,*) ""
    write(*,*) "------------------------------ d_CF_sensitivities_Central_Difference ---&
    ------------------------------"
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""

    ! allocate data holders
    allocate(CF_up(3,N_verts*3))
    allocate(CF_dn(3,N_verts*3))
    allocate(d_CF_FD(3,N_verts*3))

        
    !!!!!!!!! CENTRAL DIFFERENCE (panel 1, cp 1) d_CF_wrt_vars_test column ", k, !!!!!!!!!
    write(*,*) ""
    write(*,*) "--------------------------------------------------------------------------------------"
    write(*,*) "                   d_CF_sensitivities_Central_Difference "
    write(*,*) "--------------------------------------------------------------------------------------"
    write(*,*) ""

    
    do i=1,3
        do j=1,N_verts

            ! perturb up the current design variable
            ! before perturbation:
            write(*,('(A,f16.10)')) "Before perturbation up: ", test_mesh%vertices(j)%loc(i)
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
            write(*,('(A,f16.10)')) "AFTER perturbation up: ", test_mesh%vertices(j)%loc(i)

            
            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
        
            ! update vertex normal
            do m =1,N_panels
                deallocate(test_mesh%panels(m)%n_hat_g)
                call test_mesh%panels(m)%calc_derived_geom()
            end do
            
            call test_mesh%calc_vertex_geometry()
            
            ! update with flow
            do m =1,N_panels
                deallocate(test_mesh%panels(m)%vertices_ls)
                deallocate(test_mesh%panels(m)%n_hat_ls)
                deallocate(test_mesh%panels(m)%b)
                deallocate(test_mesh%panels(m)%b_mir)  
                deallocate(test_mesh%panels(m)%sqrt_b)
                deallocate(test_mesh%panels(m)%i_vert_d)
                deallocate(test_mesh%panels(m)%S_mu_inv)
                deallocate(test_mesh%panels(m)%T_mu)
                ! deallocate(test_mesh%panels(m)%i_panel_s)
                call test_mesh%panels(m)%init_with_flow(freestream_flow, .false., 0)
                call test_mesh%panels(m)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                test_mesh%vertices,.false.)
            end do

            ! recalculates cp locations
            deallocate(test_solver%sigma_known)
            deallocate(test_mesh%cp)
            deallocate(test_solver%P)
            call test_solver%init(solver_settings, processing_settings, &
            test_mesh, freestream_flow, control_point_file)
            
            ! Check for errors
            if (test_solver_stat /= 0) return

            deallocate(test_solver%I_known, test_solver%BC)
            deallocate(test_mesh%sigma, test_mesh%mu)
            deallocate(test_mesh%Phi_u)
            deallocate(test_mesh%C_p_inc, test_mesh%dC_f)
            deallocate(test_mesh%V_cells_inner, test_mesh%V_cells)
    

            ! deallocate(test_solver%A, test_solver%b)
            
            call test_solver%solve(test_mesh, test_solver_stat, formulation,freestream_flow)
            write(*,*) "success1"
            
            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
            
            ! get the needed info
            CF_up(:,j + (i-1)*N_verts) = test_solver%C_F(:)
            
            
            ! perturb down the current design variable
            ! write(*,*) " perturb down"
            write(*,('(A,f16.10)')) "Before perturbation down: ", test_mesh%vertices(j)%loc(i)
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) - 2.*step
            write(*,('(A,f16.10)')) "AFTER perturbation down: ", test_mesh%vertices(j)%loc(i)
            !!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!
            ! update vertex normal
            do m =1,N_panels
                deallocate(test_mesh%panels(m)%n_hat_g)
                call test_mesh%panels(m)%calc_derived_geom()
            end do
            
            call test_mesh%calc_vertex_geometry()
            
            ! update with flow
            do m =1,N_panels
                deallocate(test_mesh%panels(m)%vertices_ls)
                deallocate(test_mesh%panels(m)%n_hat_ls)
                deallocate(test_mesh%panels(m)%b)
                deallocate(test_mesh%panels(m)%b_mir)  
                deallocate(test_mesh%panels(m)%sqrt_b)
                deallocate(test_mesh%panels(m)%i_vert_d)
                deallocate(test_mesh%panels(m)%S_mu_inv)
                deallocate(test_mesh%panels(m)%T_mu)
                ! deallocate(test_mesh%panels(m)%i_panel_s)
                call test_mesh%panels(m)%init_with_flow(freestream_flow, .false., 0)
                call test_mesh%panels(m)%set_distribution(test_mesh%initial_panel_order,test_mesh%panels,&
                test_mesh%vertices,.false.)
            end do

            ! recalculates cp locations
            deallocate(test_solver%sigma_known)
            deallocate(test_mesh%cp)
            deallocate(test_solver%P)
            call test_solver%init(solver_settings, processing_settings, &
            test_mesh, freestream_flow, control_point_file)
            
            ! Check for errors
            if (test_solver_stat /= 0) return

            deallocate(test_solver%I_known, test_solver%BC)
            deallocate(test_mesh%sigma, test_mesh%mu)
            deallocate(test_mesh%Phi_u)
            deallocate(test_mesh%C_p_inc, test_mesh%dC_f)
            deallocate(test_mesh%V_cells_inner, test_mesh%V_cells)
            ! deallocate(test_solver%A, test_solver%b)
            
            call test_solver%solve(test_mesh, test_solver_stat, formulation,freestream_flow)
            write(*,*) "success2"
            !!!!!!!!!!!! END UPDATE !!!!!!!!!!!!!!!
            
            ! get the needed info
            CF_dn(:, j + (i-1)*N_verts) = test_solver%C_F(:)

            ! restore geometry
            test_mesh%vertices(j)%loc(i) = test_mesh%vertices(j)%loc(i) + step
        
        end do 
    end do 

        
        ! central difference 
    d_CF_FD = (CF_up - CF_dn)/(2.*step)
            
    

    write(*,*) ""
        
    ! write results
    write(*,*) ""
    write(*,'(A)') "              CF sensitiviites Central Difference "
    write(*,*) "       d_CFx_FD                d_CFy_FD             d_CFz_FD "
    

    do i = 1, N_verts*3
        write(*, '(3(f20.10, 4x))') d_CF_FD(1,i), d_CF_FD(2,i), d_CF_FD(3,i)
    end do 


    !!!!!!!!!!!1 write sensitivities to vtk !!!!!!!!!!!1
    allocate(vertex_normals(3,N_verts))
            
    do i=1,N_verts
        vertex_normals(:,i) = test_mesh%vertices(i)%n_g
    end do

    ! organize sensitivitiy values
    allocate(d_CF_x(3,N_verts))
    allocate(d_CF_y(3,N_verts))
    allocate(d_CF_z(3,N_verts))
    allocate(freestream_vector(3,1))

    
    ! reshape arrays
    do i=1,3
        freestream_vector(i,:) = freestream_flow%v_inf(i)
        do j=1,N_verts
            
            d_CF_x(i,j) = d_CF_FD(1, j + (i-1)*N_verts)
            d_CF_y(i,j) = d_CF_FD(2, j + (i-1)*N_verts)
            d_CF_z(i,j) = d_CF_FD(3, j + (i-1)*N_verts)
            
        end do
    end do

    ! Write geometry
    ! call body_vtk%write_point_vectors(freestream_vector, "freestream")
    call body_vtk%write_point_vectors(vertex_normals, "vertex_outward_normal_vectors")
    call body_vtk%write_point_vectors(d_CF_x, "CFx_sensitivities_central_diff")
    call body_vtk%write_point_vectors(d_CF_y, "CFy_sensitivities_central_diff")
    call body_vtk%write_point_vectors(d_CF_z, "CFz_sensitivities_central_diff")

    call body_vtk%write_point_scalars(convex, "convex")

    ! Finalize
    call body_vtk%finish()

    !!!!! end writing to vtk
    

    !!!!!!!!!!!!!!  RESULTS!!!!!!!!!!!!!
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "                          d_CF Central Difference TEST RESULTS "
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""

    
    write(*,*) ""
    call system_clock(end_count)
    time = real(end_count - start_count)/(count_rate*60.0)
    write(*,'(A,f16.10, A)') " Total test time = ", time, " minutes"
    write(*,*) ""
    write(*,*) "----------------------"
    write(*,*) "Program Complete"
    write(*,*) "----------------------"


end program test28
