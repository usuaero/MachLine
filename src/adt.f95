! Describes and handles an alternating digital tree which stores a set of vertices
module adt

    type alternating_digital_tree

        integer :: N_dims
        type(adt_node),pointer :: root => null()

        contains

            procedure :: add => adt_add

    end type alternating_digital_tree


    type adt_node

        real,dimension(3) :: point
        integer :: level
        type(adt_node),pointer :: left, right => null()

    end type adt_node

    
    contains


        recursive subroutine adt_node_add(t, point, level, N_dims)

            type(adt_node),pointer,intent(inout) :: t
            real,dimension(3),intent(in) :: point
            integer,intent(in) :: level, N_dims
            integer :: dir

            ! Check if node already exists
            if (associated(t)) then

                ! Determine sorting
                dir = modulo(level, N_dims)

                ! Send to left (lesser) node
                if (point(dir) < t%point(dir)) then
                    call adt_node_add(t%left, point, level+1, N_dims)

                ! Send to right (greater) node
                else
                    call adt_node_add(t%right, point, level+1, N_dims)
                end if

            ! Add at this node
            else
                allocate(t)
                t%point = point
                t%level = level
            end if

        end subroutine adt_node_add


        subroutine adt_add(t, point)

            class(alternating_digital_tree),intent(inout) :: t
            real,dimension(3),intent(in) :: point

            ! Enter recursion
            call adt_node_add(t%root, point, 0, t%N_dims)

        end subroutine adt_add



end module adt