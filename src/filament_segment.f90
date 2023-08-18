module filament_segment_mod
    use panel_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod 
    use wake_strip_mod 

    ! use...
    ! use...

    use...
    use...
    !!!! call correct mods


!!!! mirrors panel module

    implicit none

    type :: filament_segment !!!! changed to filament_segment type

    integer :: filament !!!! can change later. Just needed to put something so it would compile 
    
    end type filament_segment !!!! changed to filament_segment type 


end module filament_segment_mod