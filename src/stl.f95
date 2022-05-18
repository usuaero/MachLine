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
        integer :: N, stat, i, j, i1, i2, i3, N_duplicates, i_panel, i_vert, i_unique
        integer,dimension(:,:),allocatable :: panel_vertex_indices
        real,dimension(:,:),allocatable :: vertex_coords
        logical,dimension(:),allocatable :: is_duplicate
        integer,dimension(:),allocatable :: new_ind, duplicate_of

        ! Open mesh file
        open(12, file=mesh_file)

            ! Skip header
            read(12,*)

            ! Figure out how many panels are in the mesh
            N_panels = 0
            do

                ! Read line
                read(12,*,iostat=stat) dummy_read

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

        close(12)

        ! Allocate storage
        N_verts = N_panels*3
        allocate(panel_vertex_indices(N_panels,3))
        allocate(vertex_coords(N_verts,3))

        ! Reopen file
        open(12, file=mesh_file)

            ! Skip header
            read(12,*)

            ! Read in vertex locations
            do i=1,N_panels

                ! Skip normal and outer loop
                read(12,*)
                read(12,*)

                ! Get vertex locations
                read(12,*) dummy_read, vertex_coords(i*3-2,1), vertex_coords(i*3-2,2), vertex_coords(i*3-2,3)
                read(12,*) dummy_read, vertex_coords(i*3-1,1), vertex_coords(i*3-1,2), vertex_coords(i*3-1,3)
                read(12,*) dummy_read, vertex_coords(i*3,1), vertex_coords(i*3,2), vertex_coords(i*3,3)

                ! Skip end loop and end facet
                read(12,*)
                read(12,*)

                ! Point panels to vertices
                panel_vertex_indices(i,1) = i*3-2
                panel_vertex_indices(i,2) = i*3-1
                panel_vertex_indices(i,3) = i*3

            end do

        close(12)

        ! Locate duplicate vertices
        allocate(is_duplicate(N_verts), source=.false.)
        allocate(duplicate_of(N_verts))
        do i=1,N_verts
            duplicate_of(i) = i
        end do

        ! Loop through each vertex and try to find any which are the same
        !$OMP parallel do private(j) schedule(dynamic)
        do i=1,N_verts

            ! Check we don't already know this is a duplicate (this is always false for the first vertex)
            if (.not. is_duplicate(i)) then

                ! Loop through possible duplicates (don't need to check itself or any previous vertices)
                do j=i+1,N_verts

                    ! Check we don't already know this is a duplicate of another vertex
                    if (.not. is_duplicate(j)) then

                        ! Check if the vertices are the same
                        if (dist(vertex_coords(i,:), vertex_coords(j,:)) < 1.e-12) then

                            ! Mark duplicate
                            !$OMP critical
                            is_duplicate(j) = .true.
                            duplicate_of(j) = i
                            !$OMP end critical

                        end if
                    end if
                end do
            end if

        end do

        ! Collapse duplicates
        do i=1,N_verts

            ! Check if this one is a duplicate
            if (is_duplicate(i)) then

                ! Start at this vertex
                i_unique = i

                ! If it is not a duplicate of itself, move to the duplicate
                do while (.not. i_unique == duplicate_of(i_unique))
                    i_unique = duplicate_of(i_unique)
                end do

                ! Assign to new_ind
                duplicate_of(i) = i_unique

            end if

        end do

        ! Allocate new indices
        allocate(new_ind(N_verts))

        ! Loop through to determine the new indices and count up the number of duplicates
        N_duplicates = 0
        do i=1,N_verts

            ! Check if this vertex is a duplicate
            if (is_duplicate(i)) then

                ! The new index of a duplicate will be the new index of the unique vertex
                new_ind(i) = new_ind(duplicate_of(i))

                ! Update number of duplicates
                N_duplicates = N_duplicates + 1

            else

                ! If it is not a duplicate, then we simply have to shift its index down to account for vertices skipped
                new_ind(i) = i - N_duplicates

            end if
        end do

        ! Determine number of non-duplicate vertices
        N_verts = N_verts - N_duplicates

        ! Allocate vertex and panel object arrays
        allocate(vertices(N_verts))
        allocate(panels(N_panels))

        ! Initialize vertex objects, skipping duplicates
        do i=1,N_verts+N_duplicates

            ! If this vertex is not a duplicate, initialize its object
            if (.not. is_duplicate(i)) then

                ! Initialize
                j = new_ind(i)
                call vertices(j)%init(vertex_coords(i,:), j)

            end if

            ! Get index of panel and the index of the vertex in that panel
            i_panel = (i-1)/3+1
            i_vert = mod(i-1, 3)+1

            ! Point panel to non-duplicate vertex
            panel_vertex_indices(i_panel, i_vert) = new_ind(i)

        end do

        ! Initialize panel objects
        do i=1,N_panels

            ! Get vertex indices
            i1 = panel_vertex_indices(i,1)
            i2 = panel_vertex_indices(i,2)
            i3 = panel_vertex_indices(i,3)

            ! Initialize
            call panels(i)%init(vertices(i1), vertices(i2), vertices(i3), i)

            ! Add panel index to vertices
            call vertices(i1)%panels%append(i)
            call vertices(i2)%panels%append(i)
            call vertices(i3)%panels%append(i)

        end do

    end subroutine load_surface_stl
    

end module stl_mod