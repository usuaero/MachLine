module filament_mod
    use helpers_mod
    use linked_list_mod
    use base_geom_mod
    use math_mod
    use flow_mod
    use linalg_mod
    use mesh_mod
    use filament_segment_mod
    !!!! call correct mods

    !!!! mirrors wake_strip

    implicit none

    type, extends(mesh) :: filament 

        integer :: i_top_parent_1, i_top_parent_2, i_bot_parent_1, i_bot_parent_2
        integer :: i_top_parent, i_bot_parent
        logical :: on_mirror_plane
        integer :: N_segments = 0
        type(filament_segment),allocatable,dimension(:) :: segments

        contains
            procedure :: init => filament_init !!!! added this. May take it off if we don't init the wake_filament here. 
            procedure :: init_segment => filament_init_segment !!!! how to differentiate between segment and mesh
            procedure :: init_vertices => filament_init_vertices  !!!! comment out all type bound procedure statements until they are used or MachLine won't compile. -SA
            procedure :: init_segments => filament_init_segments 
        !     procedure :: init_panel => wake_strip_init_panel !!!!

    end type filament 
    
contains


    subroutine filament_init(this, freestream, starting_edge, mirror_start, mirror_plane, &    !!!! N_filaments, N_segments, freestream_plane?, normals?
                                N_segments_streamwise, trefftz_dist, body_verts, wake_mirrored, initial_panel_order, N_body_panels) 
        !initializes wake filament based on the provided info

        implicit none

        !!!! do we need to change any of these?
        class(filament),intent(inout) :: this
        type(flow),intent(in) :: freestream
        type(edge),intent(in) :: starting_edge
        logical,intent(in) :: mirror_start
        integer,intent(in) :: mirror_plane, N_segments_streamwise, initial_panel_order, N_body_panels
        real,intent(in) :: trefftz_dist
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts
        logical,intent(in) :: wake_mirrored

        real,dimension(3) :: start_1, start_2, start_c
        integer :: N_body_verts, i

        !!!!mirror stuff, do we need to change? ( line (46 - 100) in wake_strip)
        ! Get number of vertices on the body
        N_body_verts = size(body_verts)

        ! Set mirroring
        this%on_mirror_plane = starting_edge%on_mirror_plane
        this%mirror_plane = mirror_plane
        this%mirrored = wake_mirrored .and. .not. this%on_mirror_plane

        ! Get starting locations and parent vertices
        if (mirror_start) then

            ! Get starting location
            ! Note order is flipped here to maintain convention
            start_1 = mirror_across_plane(body_verts(starting_edge%top_verts(2))%loc, mirror_plane)
            start_2 = mirror_across_plane(body_verts(starting_edge%top_verts(1))%loc, mirror_plane)

            ! Get parent vertices
            this%i_top_parent_1 = starting_edge%top_verts(2) + N_body_verts
            this%i_top_parent_2 = starting_edge%top_verts(1) + N_body_verts
            this%i_bot_parent_1 = starting_edge%bot_verts(2) + N_body_verts
            this%i_bot_parent_2 = starting_edge%bot_verts(1) + N_body_verts

            ! Get parent panels
            this%i_top_parent = starting_edge%panels(1) + N_body_panels
            this%i_top_parent = starting_edge%panels(2) + N_body_panels

        else

            ! Get starting location
            start_1 = body_verts(starting_edge%top_verts(1))%loc
            start_2 = body_verts(starting_edge%top_verts(2))%loc
            start_c = (start_1 + start_2) * 0.5

            ! Get parent vertices
            this%i_top_parent_1 = starting_edge%top_verts(1)
            this%i_top_parent_2 = starting_edge%top_verts(2)
            this%i_bot_parent_1 = starting_edge%bot_verts(1)
            this%i_bot_parent_2 = starting_edge%bot_verts(2)

            ! Get parent panels
            this%i_top_parent = starting_edge%panels(1)
            this%i_top_parent = starting_edge%panels(2)

        end if

        ! Initialize vertices
        call this%init_vertices(freestream, N_segments_streamwise, trefftz_dist, start_c, body_verts)

        ! Intialize segments
        call this%init_segments(N_segments_streamwise)


        ! initalize properties
        ! for panels this is init with flow and set distribution

        ! Initialize other segment properties
        ! do i=1,this%N_segments
        !     call this%segments(i)%init_with_flow(freestream, this%mirrored, mirror_plane)
        !     call this%segments(i)%set_distribution(1, this%segments, this%vertices, this%mirrored) ! With the current formulation, wake segments are always linear
        ! end do

    end subroutine

    subroutine filament_init_vertices(this, freestream, N_segments_streamwise, trefftz_dist, &
                start_c, body_verts)
        ! Initializes this wake strip's vertices based on the provided info

        implicit none

        class(filament),intent(inout) :: this
        type(flow),intent(in) :: freestream
        integer,intent(in) :: N_segments_streamwise
        real,intent(in) :: trefftz_dist
        real,dimension(3),intent(in) :: start_c
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts

        real,dimension(3) :: loc
        real :: d1, sep
        integer :: i, N_body_verts

        ! Allocate memory
        this%N_verts = N_segments_streamwise*2
        N_body_verts = size(body_verts)
        allocate(this%vertices(this%N_verts))

        ! Initialize starting vertex in wake filament
        call this%vertices(1)%init(start_c, 1)
        

        ! Calculate distances to Trefftz plane
        d1 = trefftz_dist - inner(start_c, freestream%c_hat_g)
        

        ! Calculate spacing between vertices
        sep = d1 / N_segments_streamwise
        

        ! Loop through following vertices
        do i=2,this%N_verts

            ! Calculate location of vertex
            loc = start_c + sep*(i-1)*freestream%c_hat_g
           
            ! Initialize vertex
            call this%vertices(i)%init(loc, i)

        end do

    end subroutine filament_init_vertices

    subroutine filament_init_segments(this, N_segments_streamwise)
        ! Initializes this wake strip's segments based on the provided info

        implicit none

        class(filament),intent(inout) :: this
        integer,intent(in) :: N_segments_streamwise

        real :: d1, d2
        integer :: i, i1, i2, advance, N_skipped

        ! Determine number of segments
        this%N_segments = N_segments_streamwise
        allocate(this%segments(this%N_segments))
        i1 = 1
        i2 = 2
        do i=1,this%N_segments

            ! initialize
            call this%init_segment(i,i1,i2)
            ! increment index
            i1 = i1 + 1
            i2 = i2 + 1
            
        end do

    end subroutine filament_init_segments

    subroutine filament_init_segment(this, i_segment, i1, i2)
        ! Initializes the specified panel

        implicit none
        
        class(filament),intent(inout) :: this
        integer,intent(in) :: i_segment, i1, i2
        integer, dimension(4) :: parents

        parents = [this%i_top_parent_1, this%i_bot_parent_1, this%i_top_parent_2, this%i_bot_parent_2]
    
    
        call this%segments(i_segment)%init(this%vertices(i1), this%vertices(i2), i_segment,parents)

    end subroutine filament_init_segment

end module filament_mod