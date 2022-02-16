module tri_mod

    use panel_mod
    use vertex_mod
    use math_mod

    implicit none
    
contains


    subroutine load_surface_tri(mesh_file, N_verts, N_panels, vertices, panels)
        ! Loads a surface mesh from a tri file. Only a body.
        ! Needs to be updated to automatically delete duplicate vertices.

        implicit none

        character(len=:),allocatable,intent(in) :: mesh_file
        integer,intent(out) :: N_verts, N_panels
        type(vertex),dimension(:),allocatable,intent(out) :: vertices
        type(panel),dimension(:),allocatable,intent(out) :: panels

        integer :: i, i1, i2, i3, j
        real,dimension(3) :: vert
        character(len=200) :: dummy_read

        ! Open mesh file
        open(12, file=mesh_file)

            ! Get number of vertices and panels
            read(12,*) N_verts, N_panels

            ! Allocate vertex storate
            allocate(vertices(N_verts))

            ! Read in and initialize vertices
            do i=1,N_verts

                ! Get coords
                read(12,*) vert(1), vert(2), vert(3), dummy_read

                ! Initialize
                call vertices(i)%init(vert, i)

            end do

            ! Allocate panel storage
            allocate(panels(N_panels))

            ! Read in and initialize panels
            do i=1,N_panels

                ! Get vertex indices
                read(12,*) i1, i2, i3

                ! Initialize
                call panels(i)%init(vertices(i1), vertices(i2), vertices(i3), i)

                ! Add panel index to vertices
                call vertices(i1)%panels%append(i)
                call vertices(i2)%panels%append(i)
                call vertices(i3)%panels%append(i)

            end do

        close(12)

    end subroutine load_surface_tri
    
end module tri_mod