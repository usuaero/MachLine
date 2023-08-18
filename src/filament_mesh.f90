module filament_wake_mesh_mod  
    use json_mod
    use json_xtnsn_mod
    use linked_list_mod
    use helpers_mod
    use base_geom_mod
    use panel_mod
    use math_mod
    use flow_mod
    use vtk_mod
    use wake_strip_mod
    use filament_segment_mod
    !!!! use filament_mod
    use mesh_mod

    !!!! models wake_mesh

!****************************************
!*    filament_segment => panel         *
!*            filament => wake_strip    *
!*       filament_mesh => wake_mesh     *
!****************************************
    implicit none

    logical ::
  
    implicit none
!!!! flip definitions of filiment mesh
    !type, extends(mesh) ::  filament_mesh   !!!!(does the filament extend mesh or is it something new)
    type :: filament_wake_mesh

        type(filament_segment),allocatable,dimension(:) :: filaments
        integer :: N_segments = 0
        integer :: N_min_segments = 0 !!!! min == 1
        integer :: N_max_segments = 0 !!!! maybe set a max number of segments to limit computing time (will have to be verified experimentally)

        integer :: N_filaments !!!! shold we define this here or in filament?
        integer :: N_filaments_min
        integer :: N_filaments_max

        contains

            procedure :: init => wake_filament_mesh_init
            procedure :: init_filaments => wake_filament_mesh_init_filaments
            procedure :: write_filaments => wake_filament_mesh_write_filaments

    end type filament_wake_mesh


contains

    subroutine wake_mesh_init(this, body_edges, body_verts, freestream, asym_flow, mirror_plane, N_panels_streamwise, &
                                trefftz_dist, body_mirrored, initial_panel_order, N_body_panels) !!!!!! what else does it need
    ! Initializes the wake mesh

    !!!!.......continue

end module filament_wake_mesh_mod
