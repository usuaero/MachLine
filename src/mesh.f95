module mesh

    use json
    use json_xtnsn
    use vtk
    use geometry

    implicit none


    type surface_mesh

        integer :: N_vert, N_panel
        real,allocatable,dimension(:,:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        character(len=:),allocatable :: mesh_file

        contains

            procedure :: initialize => surface_mesh_initialize
            procedure :: output_results => surface_mesh_output_results

    end type surface_mesh


    type cart_volume_mesh

    end type cart_volume_mesh

    
contains


    subroutine surface_mesh_initialize(t, settings)

        implicit none

        class(surface_mesh),intent(inout) :: t
        type(json_value),pointer,intent(in) :: settings
        character(len=:),allocatable :: extension
        integer :: loc

        ! Get mesh file
        call json_get(settings, 'file', t%mesh_file)
        t%mesh_file = trim(t%mesh_file)
        write(*,*) "    Initializing surface mesh from file:", t%mesh_file

        ! Determine the type of mesh file
        loc = index(t%mesh_file, '.')
        extension = t%mesh_file(loc:len(t%mesh_file))

        ! Load vtk
        if (extension .eq. '.vtk') then
            call load_surface_vtk(t%mesh_file, t%N_vert, t%N_panel, t%vertices, t%panels)
        end if

        ! Display mesh info
        write(*,*) "    Surface mesh has", t%N_vert, "vertices and", t%N_panel, "panels."
    
    end subroutine surface_mesh_initialize


    subroutine surface_mesh_output_results(t, output_file)

        implicit none

        class(surface_mesh),intent(inout) :: t
        character(len=:),allocatable,intent(in) :: output_file

        ! Write out data
        call write_surface_vtk(output_file, t%vertices, t%panels)
    
    end subroutine surface_mesh_output_results

    
end module mesh