module stl_mod

    use panel_mod
    use vertex_mod
    use math_mod

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
        integer :: N, stat, i, j, i1, i2, i3, N_skipped, i_panel, i_vert
        integer,dimension(:,:),allocatable :: panel_vertex_indices
        real,dimension(:,:),allocatable :: vertex_coords
        logical,dimension(:),allocatable :: is_duplicate
        integer,dimension(:),allocatable :: duplicate_index

        ! Open mesh file
        open(1, file=mesh_file)

            ! Skip header
            read(1,*)

            ! Figure out how many panels are in the mesh
            N_panels = 0
            do

                ! Read line
                read(1,*,iostat=stat) dummy_read

                ! Check status
                if (stat == 0) then

                    ! If not at the end of the file, check for "facet", indicating a new panel
                    if (dummy_read == 'facet') then
                        N_panels = N_panels + 1
                    end if

                else
                    ! At end of file
                    exit
                end if

            end do

        close(1)

        ! Allocate storage
        N_verts = N_panels*3
        allocate(panel_vertex_indices(N_panels,3))
        allocate(vertex_coords(N_verts,3))

        ! Reopen file
        open(1, file=mesh_file)

            ! Skip header
            read(1,*)

            ! Read in vertex locations
            do i=1,N_panels

                ! Skip normal and outer loop
                read(1,*)
                read(1,*)

                ! Get vertex locations
                read(1,*) dummy_read, vertex_coords(i*3-2,1), vertex_coords(i*3-2,2), vertex_coords(i*3-2,3)
                read(1,*) dummy_read, vertex_coords(i*3-1,1), vertex_coords(i*3-1,2), vertex_coords(i*3-1,3)
                read(1,*) dummy_read, vertex_coords(i*3,1), vertex_coords(i*3,2), vertex_coords(i*3,3)

                ! Skip end loop and end facet
                read(1,*)
                read(1,*)

                ! Point panels to vertices
                panel_vertex_indices(i,1) = i*3-2
                panel_vertex_indices(i,2) = i*3-1
                panel_vertex_indices(i,3) = i*3

            end do

        close(1)

        ! Locate duplicate vertices
        allocate(is_duplicate(N_verts), source=.false.)
        allocate(duplicate_index(N_verts))
        N_skipped = 0

        ! Loop through each vertex and try to find any which are the same
        do i=1,N_verts

            ! Initialize duplicate indices
            duplicate_index(i) = i

            ! Check we don't already know this is a duplicate (this is always false for the first vertex)
            if (.not. is_duplicate(i)) then

                ! Loop through possible duplicates (don't need to check itself or any previous vertices)
                do j=i+1,N_verts

                    ! Check if the vertices are the same
                    if (dist(vertex_coords(i,:), vertex_coords(j,:)) < 1.e-12) then

                        ! Mark duplicate
                        is_duplicate(j) = .true.
                        duplicate_index(j) = i

                    end if

                end do

            else

                ! Update the number of vertices we've skipped
                N_skipped = N_skipped + 1

            end if

        end do

        ! Allocate vertex and panel arrays
        N_verts = N_verts - N_skipped
        allocate(vertices(N_verts))
        allocate(panels(N_panels))

        ! Initialize vertex objects, skipping duplicates, and panel objects
        j = 0

        ! Loop through panels
        do i_panel=1,N_panels

            ! Loop through vertices
            do i_vert=1,3

                i = 3*(i_panel-1) + i_vert

                ! If this vertex is not a duplicate, update its index and initialize its object
                if (.not. is_duplicate(i)) then

                    ! Update index in actual vertex array
                    j = j + 1

                    ! Initialize
                    call vertices(j)%init(vertex_coords(i,:), j)

                end if

                ! Point panel to non-duplicate vertex
                panel_vertex_indices(i_panel, i_vert) = duplicate_index(i)

            end do

            ! Get vertex indices
            i1 = panel_vertex_indices(i,1)
            i2 = panel_vertex_indices(i,2)
            i3 = panel_vertex_indices(i,3)

            if (i_panel==1) then
                write(*,*) i1
                write(*,*) i2
                write(*,*) i3
                write(*,*) vertices(i1)%loc
                write(*,*) vertices(i2)%loc
                write(*,*) vertices(i3)%loc
            end if

            ! Initialize
            call panels(i_panel)%init(vertices(i1), vertices(i2), vertices(i3), i1, i2, i3, i_panel)

            ! Add panel index to vertices
            call vertices(i1)%panels%append(i)
            call vertices(i2)%panels%append(i)
            call vertices(i3)%panels%append(i)

        end do

    end subroutine load_surface_stl
    
end module stl_mod