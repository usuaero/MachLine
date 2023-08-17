module filament_mod
    use helpers_mod
    use linked_list_mod
    use base_geom_mod
    use math_mod
    use flow_mod
    use linalg_mod

    !!!! models wake_strip

    implicit none

    type, extends(mesh) :: wake_filament

        integer :: i_top_parent_1, i_top_parent_2, i_bot_parent_1, i_bot_parent_2
        integer :: i_top_parent, i_bot_parent
        logical :: on_mirror_plane

        contains

            procedure :: init => wake_filament_seg_init !!!! how to differentiate between segment and mesh
            procedure :: init_vertices => wake_filament_seg_init_vertices !!!!
            procedure :: init_panels => wake_strip_init_panels !!!!
            procedure :: init_panel => wake_strip_init_panel !!!!

    end type wake_filament 
    
contains


    subroutine wake_filament_init(this, freestream, starting_edge, mirror_start, mirror_plane, &    !!!! N_filaments, N_segments, freestream_plane?, normals?
                                N_panels_streamwise, trefftz_dist, body_verts, wake_mirrored, initial_panel_order, N_body_panels) 
        !initializes wake filament based on the provided info

        implicit none

        !!!! do we need to change any of these?
        class(wake_strip),intent(inout) :: this
        type(flow),intent(in) :: freestream
        type(edge),intent(in) :: starting_edge
        logical,intent(in) :: mirror_start
        integer,intent(in) :: mirror_plane, N_panels_streamwise, initial_panel_order, N_body_panels
        real,intent(in) :: trefftz_dist
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        logical,intent(in) :: wake_mirrored

        real,dimension(3) :: start_1, start_2
        integer :: N_body_verts, i

        !!!!mirror stuff, do we need to change? (46 - 100)

    end subroutine


    subroutine wake_

end module filament_mod