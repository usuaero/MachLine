module filament_segment_mod
    use panel_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod 
    use wake_strip_mod 

    ! use...
    ! use...

    !!!! call correct mods


!!!! mirrors panel module

    implicit none

    type :: filament_segment !!!! changed to filament_segment type

    type(vertex_pointer), dimension(:),allocatable :: verticies 
    integer :: index ! index of this filament in the filament mesh array
    

    end type filament_segment !!!! changed to filament_segment type 

contains

    subroutine filamenet_segment_init(this,v1,v2,index)
        implicit none

        class(filamenet_segment), intent(inout) :: this
        type(vertex),intent(inout),target :: v1, v2
        integer, intent(in) :: index

        ! allocate vertex array
        allocate(this%verticies(2))

        ! assign pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2

        ! store index
        this%index = index

        ! still need to calculated derived geometry (based on what is needed in solver)
    end subroutine filamenet_segment_init



end module filament_segment_mod