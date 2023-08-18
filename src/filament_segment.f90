module filament_segment_mod
    use panel_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod 
    use wake_strip_mod 
<<<<<<< HEAD
    ! use...
    ! use...
=======
    use...
    use...
    !!!! call correct mods
>>>>>>> a2cd015b609efa6e84993236216ea6d5acc11766

!!!! mirrors panel module

    implicit none

<<<<<<< HEAD
    type :: filament_segment !!!! changed to filament_segment type

    integer :: filament !!!! can change later. Just needed to put something so it would compile 
    
    end type filament_segment !!!! changed to filament_segment type 
=======
    type integrals
    ...
    end type integrals


    type dod
        ! Container type for parameters of whether a panel lies in a point's domain of dependence

        logical :: in_dod = .true.
        logical,dimension(3) :: point_in_dod = .true.  !!!! checks if segment start point is in the DOD
                                                       !!!! maybe should check if buffer panel centerpoint is in DoD

    end type dod


    type filament

    integer ::
    ...


    contains
!!!! type bound procedures !!!!
        ! Initialization procedures
        ! Flow-dependent initialization procedures
        ! Mirror initialization
        ! Singularity distributions
        ! Integration
        ! Getters
        ! Checks
        ! Update Information
        ! DoD Checking
        ! Influence Calculations
        ...
        ...

    end type filament 
>>>>>>> a2cd015b609efa6e84993236216ea6d5acc11766


end module filament_segment_mod