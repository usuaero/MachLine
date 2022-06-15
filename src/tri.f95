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

        integer :: i, i1, i2, i3, j, N_words
        real,dimension(3) :: vert
        character(len=200) :: dummy_read
        logical :: more_than_three, on_space

        ! Open mesh file
        open(12, file=mesh_file)

            ! Get number of vertices and panels
            read(12,*) N_verts, N_panels

            ! Allocate vertex storate
            allocate(vertices(N_verts))

            ! Determine number of elements per line
            read(12,'(a)') dummy_read
            more_than_three = .false.

            ! Loop through each character
            on_space = .true.
            N_words = 0
            do i=1,200

                ! Check if we're on a space
                if (dummy_read(i:i) == ' ') then

                    on_space = .true.

                ! If we're not on a space but we were, then we've moved onto a word
                else if (on_space) then

                    N_words = N_words + 1
                    on_space = .false.

                    ! Check number of words found
                    if (N_words > 3) then
                        more_than_three = .true.
                        exit
                    end if
                end if

            end do

            ! Go back to beginning
            rewind(12)
            read(12,'(a)') dummy_read

            ! Read in and initialize vertices
            do i=1,N_verts

                ! Get coords
                if (more_than_three) then
                    read(12,*) vert(1), vert(2), vert(3), dummy_read
                else
                    read(12,*) vert(1), vert(2), vert(3)
                end if

                ! Initialize
                call vertices(i)%init(vert, i, 1)

            end do

            ! Check for duplicate vertices
            do i=1,N_verts
                do j=i+1,N_verts
                    if (dist(vertices(i)%loc, vertices(j)%loc) < 1e-12) then
                        write(*,*) "!!! Detected duplicate vertices in ", mesh_file, " Solution quality may be reduced."
                        write(*,*) "!!! ", i, " and ", j, " are duplicate vertices."
                    end if
                end do
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