! Types and subroutines for meshes
module mesh_mod

    use json_mod
    use json_xtnsn_mod
    use vtk_mod
    use vertex_mod
    use panel_mod
    use adt_mod
    use flow_mod
    use math_mod

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        character(len=:),allocatable :: mesh_file
        type(alternating_digital_tree) :: vertex_tree
        real :: kutta_angle, C_kutta_angle

        contains

            procedure :: init => surface_mesh_init
            procedure :: output_results => surface_mesh_output_results
            procedure :: locate_kutta_edges => surface_mesh_locate_kutta_edges

    end type surface_mesh


    type cart_volume_mesh

    end type cart_volume_mesh

    
contains


    subroutine surface_mesh_init(this, settings)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        character(len=:),allocatable :: extension
        real,dimension(3) :: p_min, p_max
        integer :: loc, i

        ! Get mesh file
        call json_get(settings, 'file', this%mesh_file)
        this%mesh_file = trim(this%mesh_file)
        write(*,*) "    Initializing surface mesh from file: ", this%mesh_file

        ! Determine the type of mesh file
        loc = index(this%mesh_file, '.')
        extension = this%mesh_file(loc:len(this%mesh_file))

        ! Load vtk
        if (extension .eq. '.vtk') then
            call load_surface_vtk(this%mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)
        end if

        ! Display mesh info
        write(*,*) "    Surface mesh has", this%N_verts, "vertices and", this%N_panels, "panels."

        ! Determine bounds of alternating digital tree
        p_min = this%vertices(1)%loc
        p_max = this%vertices(1)%loc
        do i=2,this%N_verts

            ! Check mins
            p_min(1) = min(this%vertices(i)%loc(1), p_min(1))
            p_min(2) = min(this%vertices(i)%loc(2), p_min(2))
            p_min(3) = min(this%vertices(i)%loc(3), p_min(3))

            ! Check maxs
            p_max(1) = max(this%vertices(i)%loc(1), p_max(1))
            p_max(2) = max(this%vertices(i)%loc(2), p_max(2))
            p_max(3) = max(this%vertices(i)%loc(3), p_max(3))

        end do
        
        ! Store
        this%vertex_tree%p_min = p_min
        this%vertex_tree%p_max = p_max

        ! Load vertices into alternating digital tree
        write(*,*)
        write(*,'(a)',advance='no') "     Loading vertices into ADT..."
        do i=1,this%N_verts
            call this%vertex_tree%add(this%vertices(i))
        end do
        write(*,*) "Done."

        ! Store other settings
        call json_get(settings, 'wake_shedding_angle', this%kutta_angle)
        this%C_kutta_angle = cos(this%kutta_angle*pi/180.0)
    
    end subroutine surface_mesh_init


    subroutine surface_mesh_output_results(t, output_file)

        implicit none

        class(surface_mesh),intent(inout) :: t
        character(len=:),allocatable,intent(in) :: output_file

        ! Write out data
        call write_surface_vtk(output_file, t%vertices, t%panels)
    
    end subroutine surface_mesh_output_results


    subroutine surface_mesh_locate_kutta_edges(this, freestream_flow)

        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        integer :: i, j, m, n, N_shared_verts, N_kutta_edges
        real,dimension(3) :: d
        logical :: abutting
        real :: distance

        write(*,*)
        write(*,'(a)',advance='no') "     Locating wake-shedding edges..."

        ! Loop through each pair of panels
        N_kutta_edges = 0
        do i=1,this%N_panels
            do j=i+1,this%N_panels

                ! Initialize for this panel pair
                n_shared_verts = 0
                abutting = .false.

                ! Check if the panels are abutting
                abutting_loop: do m=1,this%panels(i)%N
                    do n=1,this%panels(j)%N

                        ! Get distance between vertices
                        d = this%panels(i)%get_vertex_loc(m)-this%panels(j)%get_vertex_loc(n)
                        distance = norm(d)

                        ! Check distance
                        if (distance < 1e-10) then

                            ! Previously found a shared vertex
                            if (n_shared_verts == 1) then
                                abutting = .true.
                                exit abutting_loop

                            ! First shared vertex
                            else
                                n_shared_verts = 1
                            end if

                        end if

                    end do
                end do abutting_loop

                if (abutting) then

                    ! Check angle between panels
                    if (inner(this%panels(i)%normal, this%panels(j)%normal) < this%C_kutta_angle) then

                        ! Check angle with freestream
                        if (inner(this%panels(i)%normal, freestream_flow%V_inf) > 0.0 .or. &
                            inner(this%panels(j)%normal, freestream_flow%V_inf) > 0.0) then
                            N_kutta_edges = N_kutta_edges + 1
                        end if
                    end if
                end if

            end do
        end do

        write(*,*) "Done. Found", N_kutta_edges, "wake-shedding edges."

    end subroutine surface_mesh_locate_kutta_edges
    
end module mesh_mod