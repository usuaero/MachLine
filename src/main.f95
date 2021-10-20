program main

    use json_mod
    use json_xtnsn_mod
    use surface_mesh_mod
    use flow_mod
    use panel_solver_mod

    implicit none

    character(100) :: input_file
    character(len=:),allocatable :: body_file, wake_file, control_point_file

    type(json_file) :: input_json
    type(json_value), pointer :: flow_settings,&
                                 geometry_settings,&
                                 solver_settings,&
                                 output_settings,&
                                 surface_mesh_settings,&
                                 volume_mesh_settings
    type(surface_mesh) :: body_mesh
    type(flow) :: freestream_flow
    type(panel_solver) :: linear_solver

    ! Welcome message
    write(*,*) "           /"
    write(*,*) "          /"
    write(*,*) "         /"
    write(*,*) "        /          ____"
    write(*,*) "       /          /   /"
    write(*,*) "      /          /   /"
    write(*,*) "     /      MFTran (c) 2021 USU Aerolab"
    write(*,*) "    / _________/___/_______________"
    write(*,*) "   ( (__________________________"
    write(*,*) "    \          \   \"
    write(*,*) "     \          \   \"
    write(*,*) "      \          \   \"
    write(*,*) "       \          \___\"
    write(*,*) "        \"
    write(*,*)

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

    ! Let user know what input file is being used
    write(*,*)
    write(*,*) "MFTran called with input file ", input_file

    ! Load settings from input file
    write(*,*)
    write(*,*) "Loading input"
    call input_json%load_file(filename=input_file)
    call json_check()
    call input_json%get('flow', flow_settings)
    call input_json%get('geometry', geometry_settings)
    call input_json%get('solver', solver_settings)
    call input_json%get('output', output_settings)

    ! Initialize surface mesh
    call json_get(geometry_settings, 'surface_mesh', surface_mesh_settings)
    call body_mesh%init(surface_mesh_settings)

    ! Initialize flow
    call freestream_flow%init(flow_settings)

    write(*,*)
    write(*,*) "Initializing"

    ! Perform flow-dependent initialization on the surface mesh
    call body_mesh%init_with_flow(freestream_flow)

    write(*,*)
    write(*,*) "Running flow solvers"

    ! Initialize panel solver
    call linear_solver%init(solver_settings)

    ! Run solver
    call linear_solver%solve(body_mesh, freestream_flow)

    ! Output results
    write(*,*)
    write(*,*) "Post-processing"
    write(*,*)
    write(*,*) "    Writing results to file..."
    call json_get(output_settings, 'body_file', body_file)
    call json_get(output_settings, 'wake_file', wake_file)
    call json_get(output_settings, 'control_point_file', control_point_file)
    call body_mesh%output_results(body_file, wake_file, control_point_file)

    ! Goodbye
    write(*,*)
    write(*,*) "MFTran exited successfully."

end program main