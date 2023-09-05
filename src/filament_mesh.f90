module filament_mesh_mod  
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
    use filament_mod
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
        integer :: N_segments

        contains

            procedure :: init => filament_mesh_init 
            procedure :: init_filaments => filament_mesh_init_filaments            
            procedure :: write_filaments => wake_filament_mesh_write_filaments

    end type filament_mesh


contains

    subroutine filament_mesh_init(this, body_edges, body_verts, freestream, asym_flow, mirror_plane, N_segments_streamwise, &
                                 trefftz_dist, body_mirrored, initial_panel_order, N_body_panels) 
    
    ! Initializes the wake mesh
        implicit none

        class(filament_mesh),intent(inout) :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_segments_streamwise, initial_panel_order, N_body_panels
        real,intent(in) :: trefftz_dist
        logical,intent(in) :: body_mirrored

        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Initializing wake with a Trefftz distance of ", trefftz_dist, "..."

        ! Set whether the wake will be mirrored
        this%mirrored = body_mirrored .and. .not. asym_flow
        this%mirror_plane = mirror_plane

        ! Initialize filaments
        call this%init_filaments(body_edges, body_verts, freestream, asym_flow, mirror_plane, N_segments_streamwise, &
                                trefftz_dist, initial_panel_order, N_body_panels)

        if (verbose) write(*,'(a i7 a i7 a i7 a)') "Done. Created ", this%N_verts, " wake vertices and ", &
                                                this%N_segments, " wake segments distributed between ", &
                                                this%N_filaments, " filaments."
    
    
    end subroutine filament_mesh_init 

    subroutine filament_mesh_init_filaments(body_edges, body_verts, freestream, asym_flow, mirror_plane, N_segments_streamwise, &
        trefftz_dist, initial_panel_order, N_body_panels)
        ! creates the filaments for the wake

        implicit none

        class(filament_mesh),intent(inout),target :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_segments_streamwise, initial_panel_order, N_body_panels
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
                                           N_segments_streamwise, trefftz_dist, body_verts, this%mirrored, &
                                           initial_panel_order, N_body_panels)

            ! Check if we need to create a mirror strip
            if (asym_flow .and. .not. body_edges(i_start_edge)%on_mirror_plane) then
                i_filament = i_filament + 1
                call this%filaments(i_filament)%init(freestream, body_edges(i_start_edge), .true., mirror_plane, &
                                               N_segments_streamwise, trefftz_dist, body_verts, this%mirrored, &
                                               initial_panel_order, N_body_panels)
            end if
        end do

        !!!! this was in the wake mesh but I haven't found where it is used.  leaving it in so that if we need it we can change it -JH
        ! ! Find out the maximum number of vertices and panels between the filaments and totals
        ! this%N_verts = 0
        ! this%N_segments = 0


        do i=1,this%N_filaments

        !     ! Find maxima
        !     this%N_max_strip_panels = max(this%N_max_strip_panels, this%filaments(i)%N_segments)
        !     this%N_max_strip_verts = max(this%N_max_strip_verts, this%filaments(i)%N_verts)

            ! Sum totals
            this%N_verts = this%N_verts + this%filaments(i)%N_verts
            this%N_segments = this%N_segments + this%filaments(i)%N_segments
        end do
    end subroutine filament_mesh_init_filaments

    !!!! need to change mu to circulation
    subroutine filament_mesh_write_filaments(this, wake_file_exported, mu)
        ! Writes the wake filaments out to file

        implicit none
        
        class(filament_mesh),intent(in) :: this
        character(len=:),allocatable,intent(in) :: wake_file
        logical,intent(out) :: exported
        real,dimension(:),allocatable,intent(in),optional :: mu

        type(vtk_out) :: wake_vtk  !!!! new name?
        integer :: i, j, k, l, m, n, shift
        real,dimension(:),allocatable :: parent_mu !!!!mu_on_wake ?   
        !real,dimension(:),allocatable :: circ_filament !!!!
        real,dimension(:,:),allocatable :: verts

        ! Clear old file
        call delete_file(wake_file)

        if (this%N_filaments > 0) then
        
            ! Get all vertices
            allocate(verts(2,this%N_verts)) !!!! 2 verts instead of 3
            i = 0
            do k=1,this%N_filaments
                do j=1,this%filaments(k)%N_verts
                    i = i + 1
                    verts(:,i) = this%filaments(k)%vertices(j)%loc
                end do
            end do

            ! Initialize and write out vertices
            call wake_vtk%begin(wake_file)
            call wake_vtk%write_points(verts)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! call wake_vtk%write_filament_segments(segments) !!!! removed for testing -jjh
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Write out segments
            shift = 0
            do k=1,this%N_filaments
                call wake_vtk%write_filament_segments(this%filaments(k)%segments, mirror=.false., &
                                           vertex_index_shift=shift, N_total_segmants=N_segments)
                shift = shift + this%strips(k)%N_verts
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if (present(mu)) then
!
            !    ! Calculate doublet strengths
            !    allocate(mu_on_wake(N_verts))
            !    i = 0
            !    do k=1,this%N_strips
            !        do j=1,this%strips(k)%N_verts
            !            i = i + 1
            !            mu_on_wake(i) = mu(this%strips(k)%vertices(j)%top_parent) - mu(this%strips(k)%vertices(j)%bot_parent)
            !        end do
            !    end do
!
            !    ! Write doublet strengths
            !    !call wake_vtk%write_point_scalars(mu_on_wake, "mu") !!!! do we need this???
!
            !    ! Calculate filament circulation strength
            !    allocate(circ_filament(N_filaments)) 
            !    l = 0
            !    do m=1,this%N_filaments
            !        !do n=1,this%filaments(m)%N_verts !!!! dose this need to loop through verticies ???
            !            l = l + 1
            !            circ_filament(i) = mu_on_wake(n + 1) - mu_on_wake(n)
            !        !end do 
            !    end do
!
                ! Write doublet strengths
                ! call wake_vtk%write_point_scalars(circ_filament, "Circulation") !!!! removed for to compile -jjh

                
            end if

             Finish up
            call wake_vtk%finish()
            exported = .true.

        else
            exported = .false.
        end if

    end subroutine filament_mesh_write_filaments


end module filament_mesh_mod
