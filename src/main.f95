program main

    use json_mod
    use json_xtnsn_mod
    use surface_mesh_mod
    use flow_mod
    use panel_solver_mod
    use helpers_mod
    use omp_lib

    implicit none

    character(100) :: input_file
    character(len=:),allocatable :: body_file, wake_file, control_point_file, points_file, points_output_file, wake_strip_file
    character(len=:),allocatable :: mirrored_body_file, mirrored_control_point_file
    character(len=:),allocatable :: report_file, spanwise_axis

    type(json_file) :: input_json
    type(json_value),pointer :: flow_settings, &
                                geom_settings, &
                                solver_settings, &
                                processing_settings, &
                                output_settings, &
                                report_json, p_parent
    type(surface_mesh) :: body_mesh
    type(flow) :: freestream_flow
    type(panel_solver) :: linear_solver
    integer :: start_count, end_count, i_unit
    real :: count_rate, runtime
    logical :: exists, found, exported

    ! Start timer
    call system_clock(start_count, count_rate)

    ! Set up run
    call json_initialize()

    ! Get input file from command line
    call getarg(1, input_file)
    input_file = trim(input_file)

    ! Get input file from user
    if (input_file == '') then
        write(*,*) "Please specify an input file:"
        read(*,'(a)') input_file
        input_file = trim(input_file)
    end if

    ! Check it exists
    inquire(file=input_file, exist=exists)
    if (.not. exists) then
        write(*,*) "!!! The file ", input_file, " does not exist. Quitting..."
        stop
    end if

    ! Load settings from input file
    call input_json%load_file(filename=input_file)
    call json_check()
    call input_json%get('flow', flow_settings, found)
    call input_json%get('geometry', geom_settings, found)
    call input_json%get('solver', solver_settings, found)
    call input_json%get('post_processing', processing_settings, found)
    call input_json%get('output', output_settings, found)

    ! Get checks and verbose toggles
    call json_xtnsn_get(solver_settings, 'run_checks', run_checks, .false.)
    call json_xtnsn_get(output_settings, 'verbose', verbose, .true.)

    ! Welcome message
    if (verbose) then
        write(*,*) "           /"
        write(*,*) "          /"
        write(*,*) "         /"
        write(*,*) "        /          ____"
        write(*,*) "       /          /   /"
        write(*,*) "      /     MachLine (c) 2022 USU Aerolab"
        write(*,*) "     /          /   /    v2.3"
        write(*,*) "    / _________/___/_______________"
        write(*,*) "   ( (__________________________"
        write(*,*) "    \          \   \"
        write(*,*) "     \          \   \"
        write(*,*) "      \          \   \"
        write(*,*) "       \          \___\"
        write(*,*) "        \"
        write(*,*)
        write(*,*) "Got input file: ", input_file
        write(*,*)
        write(*,*) "Reading and analyzing surface mesh"
    end if

    ! Initialize report JSON
    call json_value_create(report_json)
    call to_object(report_json, 'report') ! Meaningless, but necessary
    call json_value_create(p_parent)
    call to_object(p_parent, 'info')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, 'generated_by', 'MachLine (c) 2022 USU Aerolab')
    call json_value_add(p_parent, 'executed', fdate())
    nullify(p_parent)

    ! Initialize surface mesh
    call body_mesh%init(geom_settings)

    ! Initialize flow
    call json_xtnsn_get(geom_settings, 'spanwise_axis', spanwise_axis, '+y')
    call freestream_flow%init(flow_settings, spanwise_axis)

    if (verbose) then
        write(*,*)
        write(*,*) "Initializing based on flow properties"
    end if
    
    ! Get result files
    call json_xtnsn_get(output_settings, 'body_file', body_file, 'none')
    call json_xtnsn_get(output_settings, 'wake_file', wake_file, 'none')
    call json_xtnsn_get(output_settings, 'wake_strip_file', wake_strip_file, 'none')
    call json_xtnsn_get(output_settings, 'control_point_file', control_point_file, 'none')
    call json_xtnsn_get(output_settings, 'mirrored_body_file', mirrored_body_file, 'none')
    call json_xtnsn_get(output_settings, 'mirrored_control_point_file', mirrored_control_point_file, 'none')
    call json_xtnsn_get(output_settings, 'offbody_points.points_file', points_file, 'none')
    call json_xtnsn_get(output_settings, 'offbody_points.output_file', points_output_file, 'none')

    ! Perform flow-dependent initialization on the surface mesh
    call body_mesh%init_with_flow(freestream_flow, body_file, wake_file)

    ! Initialize panel solver
    call linear_solver%init(solver_settings, processing_settings, body_mesh, freestream_flow, control_point_file)

    if (verbose) then
        write(*,*)
        write(*,*) "Running solver using the ", linear_solver%formulation, " boundary-condition formulation"
    end if

    ! Run solver
    call json_xtnsn_get(output_settings, 'report_file', report_file, 'none')
    call linear_solver%solve(body_mesh)

    ! Update report
    call linear_solver%update_report(report_json, body_mesh)

    ! Output slice
    if (points_file /= 'none' .and. points_output_file /= 'none') then
        call linear_solver%export_off_body_points(points_file, points_output_file, body_mesh)
    end if

    ! Write input
    call json_value_create(p_parent)
    call to_object(p_parent, 'input')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, flow_settings) ! Somehow this writes all the settings...
    nullify(p_parent)

    if (verbose) then
        write(*,*)
        write(*,*) "Writing results to file"
    end if

    ! Output mesh results
    call body_mesh%output_results(body_file, wake_file, control_point_file, mirrored_body_file, mirrored_control_point_file)
    if (points_file /= 'none' .and. points_output_file /= 'none' .and. verbose) then
        write(*,'(a30 a)') "     Off-Body Points: ", points_output_file
    end if

    ! Wake strips
    if (wake_strip_file /= 'none') then
        call body_mesh%wake%write_strips(wake_strip_file, exported, body_mesh%mu)
    end if

    ! Figure out how long this took
    call system_clock(end_count)
    runtime = real(end_count - start_count)/count_rate
    call json_value_add(report_json, "total_runtime", runtime)

    ! Write report
    if (report_file /= 'none') then
        open(newunit=i_unit, file=report_file, status='REPLACE')
        call json_print(report_json, i_unit)
        close(i_unit)
        if (verbose) write(*,'(a30 a)') "    Report: ", report_file
    end if

    ! Destroy pointers
    call json_destroy(report_json)

    ! Goodbye
    if (verbose) write(*,*)
    write(*,'(a, f10.4, a)') " MachLine exited successfully. Execution time: ", runtime, " s"

end program main