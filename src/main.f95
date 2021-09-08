program main

    use math
    use json
    use json_xtnsn
    use mesh

    implicit none

    character(100) :: input_file

    type(json_file) :: input_json
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
    write(*,*)

    ! Load
    write(*,*) "Loading input..."
    call input_json%load_file(filename=input_file)
    call json_check()

    ! Initialize surface mesh
    call body_mesh%initialize(input_json)

end program main