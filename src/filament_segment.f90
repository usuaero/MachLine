module filament_segment_mod
    use panel_mod
    use base_geom_mod
    use helpers_mod
    use mesh_mod 
    use wake_strip_mod 
    use...
    use...
    !!!! call correct mods

!!!! mirrors panel module

    implicit none

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


end module filament_segment_mod