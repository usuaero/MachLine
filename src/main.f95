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
    character(len=:),allocatable :: body_file, wake_file, control_point_file, slice_file
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
    real :: start, end
    logical :: exists
    integer :: i_unit

    ! Start timer
    call cpu_time(start)
    !$ start = omp_get_wtime()

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
    call input_json%get('flow', flow_settings)
    call input_json%get('geometry', geom_settings)
    call input_json%get('solver', solver_settings)
    call input_json%get('post_processing', processing_settings)
    call input_json%get('output', output_settings)

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
        write(*,*) "      /          /   /"
        write(*,*) "     /     MachLine (c) 2022 USU Aerolab"
        write(*,*) "    / _________/___/_______________"
        write(*,*) "   ( (__________________________"
        write(*,*) "    \          \   \"
        write(*,*) "     \          \   \"
        write(*,*) "      \          \   \"
        write(*,*) "       \          \___\"
        write(*,*) "        \"
        write(*,*) "Loading input file: ", input_file
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
        write(*,*) "Initializing"
    end if
    
    ! Get result files
    call json_xtnsn_get(output_settings, 'body_file', body_file, 'none')
    call json_xtnsn_get(output_settings, 'wake_file', wake_file, 'none')
    call json_xtnsn_get(output_settings, 'control_point_file', control_point_file, 'none')
    call json_xtnsn_get(output_settings, 'mirrored_body_file', mirrored_body_file, 'none')
    call json_xtnsn_get(output_settings, 'mirrored_control_point_file', mirrored_control_point_file, 'none')
    call json_xtnsn_get(output_settings, 'potential_slice.file', slice_file, 'none')

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
    if (slice_file /= 'none') then
        call linear_solver%export_potential_slice(slice_file, output_settings, body_mesh, freestream_flow)
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

    ! Write report
    if (report_file /= 'none') then
        open(newunit=i_unit, file=report_file, status='REPLACE')
        call json_print(report_json, i_unit)
        close(i_unit)
        if (verbose) write(*,'(a30 a)') "    Report: ", report_file
    end if

    ! Destroy pointers
    call json_destroy(report_json)

    ! Output mesh results
    call body_mesh%output_results(body_file, wake_file, control_point_file, mirrored_body_file, mirrored_control_point_file)
    if (slice_file /= 'none' .and. verbose) write(*,'(a30 a)') "               Slice: ", slice_file

    ! Goodbye
    call cpu_time(end)
    !$ end = omp_get_wtime()
    if (verbose) write(*,*)
    write(*,'(a, f10.4, a)') " MachLine exited successfully. Execution time: ", end-start, " s"

end program main