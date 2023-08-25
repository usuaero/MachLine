module filament_mod
    use helpers_mod
    use linked_list_mod
    use base_geom_mod
    use math_mod
    use flow_mod
    use linalg_mod
    use mesh_mod
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

        real,dimension(3) :: start_1, start_2
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

            ! Get parent segments
            this%i_top_parent = starting_edge%segments(1) + N_body_panels
            this%i_top_parent = starting_edge%segments(2) + N_body_panels

        else

            ! Get starting location
            start_1 = body_verts(starting_edge%top_verts(1))%loc
            start_2 = body_verts(starting_edge%top_verts(2))%loc
            start_c = start_1 + start_2 / 2

            ! Get parent vertices
            this%i_top_parent_1 = starting_edge%top_verts(1)
            this%i_top_parent_2 = starting_edge%top_verts(2)
            this%i_bot_parent_1 = starting_edge%bot_verts(1)
            this%i_bot_parent_2 = starting_edge%bot_verts(2)

            ! Get parent segments
            this%i_top_parent = starting_edge%segments(1)
            this%i_top_parent = starting_edge%segments(2)

        end if

        ! Initialize vertices
        call this%init_vertices(freestream, N_segments_streamwise, trefftz_dist, start_1, start_2, body_verts)

        ! Intialize segments
        call this%init_segments(N_segments_streamwise)

        ! Initialize other panel properties
        do i=1,this%N_segments
            call this%segments(i)%init_with_flow(freestream, this%mirrored, mirror_plane)
            call this%segments(i)%set_distribution(1, this%segments, this%vertices, this%mirrored) ! With the current formulation, wake segments are always linear
        end do

    end subroutine

    subroutine filament_init_vertices(this, freestream, N_segments_streamwise, trefftz_dist, &
                start_1, start_2, body_verts)
        ! Initializes this wake strip's vertices based on the provided info

        implicit none

        class(filament),intent(inout) :: this
        type(flow),intent(in) :: freestream
        integer,intent(in) :: N_segments_streamwise
        real,intent(in) :: trefftz_dist
        real,dimension(3),intent(in) :: start_1, start_2
        type(vertex),dimension(:),allocatable,intent(in) :: body_verts

        real,dimension(3) :: loc
        real :: d1, d2, sep_1, sep_2
        integer :: i, N_body_verts

        ! Allocate memory
        this%N_verts = N_segments_streamwise*2
        N_body_verts = size(body_verts)
        allocate(this%vertices(this%N_verts))

        ! Initialize starting vertex in wake filament
        call this%vertices(1)%init(start_1, 1)
        

        ! Set parents
        this%vertices(1)%top_parent = this%i_top_parent_1
        this%vertices(1)%bot_parent = this%i_bot_parent_1
        this%vertices(2)%top_parent = this%i_top_parent_2
        this%vertices(2)%bot_parent = this%i_bot_parent_2

        ! Calculate distances to Trefftz plane
        d1 = trefftz_dist - inner(start_1, freestream%c_hat_g)
        d2 = trefftz_dist - inner(start_2, freestream%c_hat_g)

        ! Calculate spacing between vertices
        sep_1 = d1 / N_segments_streamwise
        sep_2 = d2 / N_segments_streamwise

        ! Loop through following vertices
        do i=2,this%N_verts

            ! Calculate location of vertices
            if (modulo(i, 2) == 0) then
                loc = start_2 + sep_2*(i-1)/2*freestream%c_hat_g
            else
                loc = start_1 + sep_1*i/2*freestream%c_hat_g
            end if

            ! Initialize vertices
            call this%vertices(i)%init(loc, i)

            ! Set parents
            if (modulo(i, 2) == 0) then
                this%vertices(i)%top_parent = this%i_top_parent_2
                this%vertices(i)%bot_parent = this%i_bot_parent_2
            else
                this%vertices(i)%top_parent = this%i_top_parent_1
                this%vertices(i)%bot_parent = this%i_bot_parent_1
            end if

        end do

    end subroutine filament_init_vertices

    subroutine filament_init_segments(this, N_segments_streamwise)
        ! Initializes this wake strip's segments based on the provided info

        implicit none

        class(wake_strip),intent(inout) :: this
        integer,intent(in) :: N_segments_streamwise

        real :: d1, d2
        integer :: i, j, i1, i2, advance, N_skipped

        ! Determine number of segments
        this%N_segments = N_segments_streamwise*2
        allocate(this%segments(this%N_segments))

        ! Create segments
        i1 = 1
        i2 = 2
        do i=1,this%N_segments
            !!!! commented this out until init_segment is completed. need more info -JH
            ! ! Initialize
            ! advance = 0

            ! ! See if one is at the end
            ! if (i1 == this%N_verts - 1) then
            !     advance = 2
            ! else if (i2 == this%N_verts) then
            !     advance = 1
            ! else

            !     ! Check lengths of hypotenuses
            !     d1 = dist(this%vertices(i1+2)%loc, this%vertices(i2)%loc)
            !     d2 = dist(this%vertices(i1)%loc, this%vertices(i2+2)%loc)

            !     ! Pick which one
            !     if (d1 < d2) then
            !         advance = 1
            !     else
            !         advance = 2
            !     end if
            ! end if

            ! ! Advance
            ! if (advance == 1) then

            !     ! Initialize
            !     call this%init_panel(i, i1, i1+2, i2, skipped_panels)
                
            !     ! Increment index
            !     i1 = i1 + 2
            ! else
                
            !     ! Initialize
            !     call this%init_panel(i, i1, i2+2, i2, skipped_panels)

            !     ! Increment index
            !     i2 = i2 + 2
            ! end if

        end do

    end subroutine filament_init_segments

    subroutine filament_init_segment(this, i_segment, i1, i2)
        ! Initializes the specified panel

        implicit none
        
        class(filament),intent(inout) :: filament
        integer,intent(in) :: i_segment, i1, i2, i3
    
        call this%segmetns(i_segment)%init(this%vertices(i1), this%vertices(i2), i_segment)

    end subroutine filament_init_segment
end module filament_mod