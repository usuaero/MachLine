program main

    use math_mod
    use json_mod
    use json_xtnsn_mod
    use mesh_mod

    implicit none

    character(100) :: input_file
    character(len=:),allocatable :: output_file

    type(json_file) :: input_json
    type(json_value), pointer :: flow_settings,&
                                 geometry_settings,&
                                 solver_settings,&
                                 output_settings,&
                                 surface_mesh_settings,&
                                 volume_mesh_settings
    type(surface_mesh) :: body_mesh
    type(cart_volume_mesh) :: volume_mesh

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

    ! Get input file
    call getarg(1, input_file)
    write(*,*) "MFTran called with input file ", input_file

    ! Load settings from input file
    write(*,*)
    write(*,*) "Loading input..."
    call input_json%load_file(filename=input_file)
    call json_check()
    call input_json%get('flow', flow_settings)
    call input_json%get('geometry', geometry_settings)
    call input_json%get('solver', solver_settings)
    call input_json%get('output', output_settings)

    ! Initialize surface mesh
    call json_get(geometry_settings, 'surface_mesh', surface_mesh_settings)
    call body_mesh%initialize(surface_mesh_settings)

    ! Output results
    write(*,*)
    write(*,*) "Writing results to file..."
    call json_get(output_settings, 'file', output_file)
    call body_mesh%output_results(output_file)

    ! Goodbye
    write(*,*)
    write(*,*) "MFTran exited successfully."

end program main