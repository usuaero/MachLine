module stl_mod

    use panel_mod
    use vertex_mod

    implicit none
    
contains

    subroutine load_surface_stl(mesh_file, N_verts, N_panels, vertices, panels)
        ! Loads a surface mesh from an stl file

        implicit none

        character(len=:),allocatable,intent(in) :: mesh_file
        integer,intent(out) :: N_verts, N_panels
        type(vertex),dimension(:),allocatable,intent(inout) :: vertices
        type(panel),dimension(:),allocatable,intent(inout) :: panels

        character(len=200) :: dummy_read
        real,dimension(3) :: vertex_loc
        integer :: N

        ! Open file
        open(1, file=mesh_file)

            ! Skip header
            read(1,*)

        close (1)
    
    end subroutine load_surface_stl
    
end module stl_mod