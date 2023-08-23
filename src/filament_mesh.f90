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
    !!!! call correct mods

    !!!! mirrors wake_mesh

!****************************************
!*    filament_segment => panel         *
!*            filament => wake_strip    *
!*       filament_mesh => wake_mesh     *
!****************************************
    implicit none

    ! logical ::  !!!! commenting out so it will compile. Need to put a variable here.  
  
    ! implicit none !!!! commenting out so it will compile 

    type, extends(mesh) ::  filament_mesh   !!!!(does the filament extend mesh or is it something new) (extends mesh - JH)

        type(filament),allocatable,dimension(:) :: filaments
        integer :: N_min_segments = 0 !!!! min == 1
        integer :: N_max_segments = 0 !!!! maybe set a max number of segments to limit computing time (will have to be verified experimentally)

        integer :: N_filaments 
        integer :: N_max_strip_verts = 0
        integer :: N_max_strip_panels = 0

        contains

            procedure :: init => filament_mesh_init 
            procedure :: init_filaments => filament_mesh_init_filaments             !!!! comment out all type bound procedure statements until they are used or MachLine won't compile. -SA
            ! procedure :: write_filaments => wake_filament_mesh_write_filaments

    end type filament_mesh


contains

            !!!! I am commenting wake_mesh_init out because if we don't declare the variables, it won't compile
    subroutine filament_mesh_init(this, body_edges, body_verts, freestream, asym_flow, mirror_plane, N_panels_streamwise, &
                                 trefftz_dist, body_mirrored, initial_panel_order, N_body_panels) !!!!!! what else does it need
    
    ! Initializes the wake mesh
        implicit none

        class(filament_wake_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_panels_streamwise, initial_panel_order, N_body_panels
        real,intent(in) :: trefftz_dist
        logical,intent(in) :: body_mirrored

        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Initializing wake with a Trefftz distance of ", trefftz_dist, "..."

        ! Set whether the wake will be mirrored
        this%mirrored = body_mirrored .and. .not. asym_flow
        this%mirror_plane = mirror_plane

        ! Initialize filaments
        call this%init_filaments(body_edges, body_verts, freestream, asym_flow, mirror_plane, N_panels_streamwise, &
                                trefftz_dist, initial_panel_order, N_body_panels)

        if (verbose) write(*,'(a i7 a i7 a i7 a)') "Done. Created ", this%N_verts, " wake vertices and ", &
                                                this%N_panels, " wake panels distributed between ", &
                                                this%N_filaments, " filaments."
    
    
    end subroutine filament_mesh_init 

    subroutine filament_mesh_init_filaments(body_edges, body_verts, freestream, asym_flow, mirror_plane, N_panels_streamwise, &
        trefftz_dist, initial_panel_order, N_body_panels)
        ! creates the filaments for the wake

        implicit none

        class(filament_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_panels_streamwise, initial_panel_order, N_body_panels
        real,intent(in) :: trefftz_dist

        integer :: i, i_filament, i_start_edge
        type(list) :: wake_shedding_edges

        ! Loop through edges to find which ones shed a wake and how many there are
        this%N_filaments = 0
        do i=1,size(body_edges)

            ! Check if it sheds a wake
            if (body_edges(i)%sheds_wake) then
                this%N_filaments = this%N_filaments + 1
                call wake_shedding_edges%append(i)

                ! If the flow is asymmetric, then wake-shedding edges not on the mirror plane will also be used to generate a wake strip
                if (asym_flow .and. .not. body_edges(i)%on_mirror_plane) then
                    this%N_filaments = this%N_filaments + 1
                end if
            end if

        end do

        ! Allocate strip storage
        allocate(this%filaments(this%N_filaments))
    
        ! Initialize filaments
        i = 0
        i_filament = 0
        do while (i_filament < this%N_filaments)

            ! Get starting edge index
            i = i + 1
            call wake_shedding_edges%get(i, i_start_edge)

            ! Initialize strip
            i_filament = i_filament + 1
            call this%filaments(i_filament)%init(freestream, body_edges(i_start_edge), .false., mirror_plane, &
                                           N_panels_streamwise, trefftz_dist, body_verts, this%mirrored, &
                                           initial_panel_order, N_body_panels)

            ! Check if we need to create a mirror strip
            if (asym_flow .and. .not. body_edges(i_start_edge)%on_mirror_plane) then
                i_filament = i_filament + 1
                call this%filaments(i_filament)%init(freestream, body_edges(i_start_edge), .true., mirror_plane, &
                                               N_panels_streamwise, trefftz_dist, body_verts, this%mirrored, &
                                               initial_panel_order, N_body_panels)
            end if
        end do

        !!!! this was in the wake mesh but I haven't found where it is used.  leaving it in so that if we need it we can change it -JH
        ! ! Find out the maximum number of vertices and panels between the filaments and totals
        ! this%N_verts = 0
        ! this%N_panels = 0
        ! do i=1,this%N_filaments

        !     ! Find maxima
        !     this%N_max_strip_panels = max(this%N_max_strip_panels, this%filaments(i)%N_panels)
        !     this%N_max_strip_verts = max(this%N_max_strip_verts, this%filaments(i)%N_verts)

        !     ! Sum totals
        !     this%N_verts = this%N_verts + this%filaments(i)%N_verts
        !     this%N_panels = this%N_panels + this%filaments(i)%N_panels
        ! end do

end module filament_wake_mesh_mod
