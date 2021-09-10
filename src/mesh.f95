! Types and subroutines for meshes
module mesh

    use json
    use json_xtnsn
    use vtk
    use geometry
    use adt

    implicit none


    type surface_mesh

        integer :: N_verts, N_panels
        type(vertex),allocatable,dimension(:) :: vertices
        type(panel_pointer),allocatable,dimension(:) :: panels
        character(len=:),allocatable :: mesh_file
        type(alternating_digital_tree) :: vertex_tree

        contains

            procedure :: initialize => surface_mesh_initialize
            procedure :: output_results => surface_mesh_output_results

    end type surface_mesh


    type cart_volume_mesh

    end type cart_volume_mesh

    
contains


    subroutine surface_mesh_initialize(this, settings)

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
        write(*,*) "    Loading vertices into ADT..."
        do i=1,this%N_verts
            call this%vertex_tree%add(this%vertices(i))
        end do
    
    end subroutine surface_mesh_initialize


    subroutine surface_mesh_output_results(t, output_file)

        implicit none

        class(surface_mesh),intent(inout) :: t
        character(len=:),allocatable,intent(in) :: output_file

        ! Write out data
        call write_surface_vtk(output_file, t%vertices, t%panels)
    
    end subroutine surface_mesh_output_results

    
end module mesh