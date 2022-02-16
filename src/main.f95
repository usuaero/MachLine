program main

    use json_mod
    use json_xtnsn_mod
    use surface_mesh_mod
    use flow_mod
    use panel_solver_mod

    implicit none

    character(100) :: input_file
    character(len=:),allocatable :: body_file, wake_file, control_point_file, report_file, spanwise_axis

    type(json_file) :: input_json
    type(json_value), pointer :: flow_settings,&
                                 geometry_settings,&
                                 solver_settings,&
                                 processing_settings,&
                                 output_settings
    type(surface_mesh) :: body_mesh
    type(flow) :: freestream_flow
    type(panel_solver) :: linear_solver
    real :: start, end

    ! Initialize developer things
    eval_count = 0

    ! Start timer
    call cpu_time(start)

    ! Welcome message
    write(*,*) "           /"
    write(*,*) "          /"
    write(*,*) "         /"
    write(*,*) "        /          ____"
    write(*,*) "       /          /   /"
    write(*,*) "      /          /   /"
    write(*,*) "     /     MachLine (c) 2022 USU Aerolab"
    write(*,*) "    / _________/___/_______________"
    write(*,*) "   ( (__________________________"
    write(*,*) "    \          \   \"
    write(*,*) "     \          \   \"
    write(*,*) "      \          \   \"
    write(*,*) "       \          \___\"
    write(*,*) "        \"

    ! Set up run
    call json_initialize()

    ! Get input file from command line
    call getarg(1, input_file)

    ! Get input file from user
    if (input_file == '') then
        write(*,*) "Please specify an input file:"
        read(*,'(a)') input_file
        input_file = trim(input_file)
    end if

    ! Load settings from input file
    write(*,*) "Loading input file: ", input_file
    call input_json%load_file(filename=input_file)
    call json_check()
    call input_json%get('flow', flow_settings)
    call input_json%get('geometry', geometry_settings)
    call input_json%get('solver', solver_settings)
    call input_json%get('post_processing', processing_settings)
    call input_json%get('output', output_settings)

    ! Initialize surface mesh
    call body_mesh%init(geometry_settings)

    ! Initialize flow
    call json_xtnsn_get(geometry_settings, 'spanwise_axis', spanwise_axis, '+y')
    call freestream_flow%init(flow_settings, spanwise_axis)

    write(*,*)
    write(*,*) "Initializing"
    
    ! Get result files
    call json_xtnsn_get(output_settings, 'wake_file', wake_file, 'none')
    call json_xtnsn_get(output_settings, 'control_point_file', control_point_file, 'none')

    ! Perform flow-dependent initialization on the surface mesh
    call body_mesh%init_with_flow(freestream_flow, wake_file)

    ! Initialize panel solver
    call linear_solver%init(solver_settings, processing_settings, body_mesh, freestream_flow, control_point_file)

    write(*,*)
    write(*,*) "Running solver using ", linear_solver%formulation, " formulation"

    ! Run solver
    call json_xtnsn_get(output_settings, 'report_file', report_file, 'none')
    call linear_solver%solve(body_mesh, report_file)

    ! Output results
    write(*,*)
    write(*,*) "Writing results to file"
    call json_xtnsn_get(output_settings, 'body_file', body_file, 'none')

    call body_mesh%output_results(body_file, wake_file, control_point_file)

    ! Goodbye
    call cpu_time(end)
    write(*,*)
    write(*,'(a, f10.4, a)') "MachLine exited successfully. Execution time: ", end-start, " s"

end program main