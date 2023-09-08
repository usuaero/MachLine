module filament_segment_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod 

!!!! mirrors panel module

    implicit none

    type :: filament_segment !!!! changed to filament_segment type

        type(vertex_pointer), dimension(:),allocatable :: vertices 
        integer :: index ! index of this filament in the filament mesh array
        integer :: N = 2 ! number of vertices

        contains
            ! initilizers 
            procedure :: init => filament_segment_init

            ! Getters
            procedure :: get_vertex_index => filament_segment_get_vertex_index
            procedure :: get_vertex_loc => filament_segment_get_vertex_loc

    end type filament_segment !!!! changed to filament_segment type 

contains

    subroutine filament_segment_init(this,v1,v2,index)
        implicit none

        class(filament_segment), intent(inout) :: this
        type(vertex),intent(inout),target :: v1, v2
        integer, intent(in) :: index

        ! allocate vertex array
        allocate(this%vertices(2))

        ! assign pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2

        ! store index
        this%index = index

        ! still need to calculated derived geometry (based on what is needed in solver)
    end subroutine filament_segment_init

    function filament_segment_get_vertex_loc(this, i) result(loc)

        implicit none
        class(filament_segment), intent(in) :: this
        integer,intent(in) :: i
        real,dimension(3) :: loc

        loc = this%vertices(i)%ptr%loc

    end function filament_segment_get_vertex_loc

    function filament_segment_get_vertex_index(this, i) result(index)

        implicit none

        class(filament_segment),intent(in) :: this
        integer,intent(in) :: i
        integer :: index

        index = this%vertices(i)%ptr%index

    end function filament_segment_get_vertex_index

end module filament_segment_mod