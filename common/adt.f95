! Describes and handles an alternating digital tree which stores a set of vertices in 3-space
module adt_mod

    use vertex_mod

    type alternating_digital_tree

        real,dimension(3) :: p_min, p_max
        type(adt_node),pointer :: root => null()
        integer :: N_nodes = 0

        contains

            procedure :: add => adt_add
            procedure :: recursive_add => adt_recursive_add
            procedure :: write_vtk => adt_write_vtk

    end type alternating_digital_tree


    type adt_node

        type(vertex),pointer :: point
        integer :: level, dir
        real,dimension(3) :: p_min, p_max
        real :: div
        type(adt_node),pointer :: left, right => null()

    end type adt_node

    
contains


    recursive subroutine adt_recursive_add(t, curr, point, level, p_min, p_max)

        class(alternating_digital_tree),intent(inout) :: t
        type(adt_node),pointer,intent(inout) :: curr
        real,dimension(3),intent(in) :: p_min, p_max
        type(vertex),pointer,intent(in) :: point
        real,dimension(3) :: p_mid
        integer,intent(in) :: level
        integer :: dir

        ! Check if current node already exists
        if (associated(curr)) then

            ! Send to left (lesser) node
            if (point%loc(curr%dir) < curr%div) then
                p_mid = p_max
                p_mid(curr%dir) = curr%div
                call t%recursive_add(curr%left, point, level+1, p_min, p_mid)

            ! Send to right (greater) node
            else
                p_mid = p_min
                p_mid(curr%dir) = curr%div
                call t%recursive_add(curr%right, point, level+1, p_mid, p_max)
            end if

        ! Add at this node
        else

            allocate(curr)
            curr%point => point
            curr%level = level
            curr%p_min = p_min
            curr%p_max = p_max

            ! Determine sorting direction
            curr%dir = modulo(level, 3)+1

            ! Get dividing location
            curr%div = 0.5*(p_min(curr%dir)+p_max(curr%dir))
        end if

    end subroutine adt_recursive_add


    subroutine adt_add(t, point)

        class(alternating_digital_tree),intent(inout) :: t
        type(vertex),intent(in),target :: point
        type(vertex),pointer :: ptr

        ! Create pointer to point passed in
        ptr => point

        ! Enter recursion
        call t%recursive_add(t%root, ptr, 0, t%p_min, t%p_max)

        ! Increment count
        t%N_nodes = t%N_nodes + 1

    end subroutine adt_add


    subroutine adt_write_vtk(t, filename)

        implicit none

        class(alternating_digital_tree),intent(in) :: t
        character(100),intent(in) :: filename

        ! Formats for writing out results
        100 format(f20.12, ' ', f20.12, ' ', f20.12) ! Vertices
        101 format(i1, ' ', i20, ' ', i20, ' ', i20) ! Panel indices

        ! Open file
        open(1, file=filename)

            ! Write header
            write(1,'(a)') "# vtk DataFile Version 3.0"
            write(1,'(a)') "MFTran ADT file. Generated by MFTran, USU AeroLab (c) 2021."
            write(1,'(a)') "ASCII"

            !! Write out vertices
            !N_vert = size(vertices)/3
            !write(1,'(a)') "DATASET POLYDATA"
            !write(1,'(a i20 a)') "POINTS", N_vert, " float"
            !do i=1,N_vert
            !    write(1,100) vertices(i,1), vertices(i,2), vertices(i,3)
            !end do

            !! Determine panel info size
            !panel_info_size = 0
            !N_panel = size(panels)
            !do i=1,N_panel
            !    panel_info_size = panel_info_size + panels(i)%N + 1
            !end do

            !! Write out panels
            !write(1,'(a i20 i20)') "POLYGONS", N_panel, panel_info_size
            !do i=1,N_panel
            !    write(1,101) panels(i)%N, panels(i)%i1, panels(i)%i2, panels(i)%i3
            !end do

        close(1)
    
    end subroutine adt_write_vtk

end module adt_mod