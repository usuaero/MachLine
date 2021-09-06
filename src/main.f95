program main
    implicit none
    character(100) :: input_file

    ! Welcome message
    write(*,*) "       /"
    write(*,*) "      /"
    write(*,*) "     /"
    write(*,*) "    /          ____"
    write(*,*) "   /          /   /"
    write(*,*) "  /          /   /"
    write(*,*) " /      MFTran (c) 2021 USU Aerolab"
    write(*,*) "/ _________/___/________________"
    write(*,*) " (___________________________"
    write(*,*) "\          \   \"
    write(*,*) " \          \   \"
    write(*,*) "  \          \   \"
    write(*,*) "   \          \___\"
    write(*,*) "    \"
    write(*,*)

    ! Get input file
    call getarg(1, input_file)
    write(*,*) "MFTran called with input file ", input_file
    write(*,*)
    write(*,*) "Loading input..."
    
end program main