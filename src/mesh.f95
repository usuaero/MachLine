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
        type(json_file) :: input_json

        contains

            procedure :: initialize => surface_mesh_initialize
            procedure :: output_results => surface_mesh_output_results

    end type surface_mesh


    type cart_volume_mesh

        type(json_file) :: input_json

    end type cart_volume_mesh

    
contains


    subroutine surface_mesh_initialize(t, input_json)

        implicit none

        class(surface_mesh),intent(out) :: t
        type(json_file),intent(in) :: input_json
        character(len=:),allocatable :: mesh_file, extension
        integer :: loc

        ! Get mesh file
        t%input_json = input_json
        call t%input_json%get('geometry.surface_mesh.file', mesh_file)
        mesh_file = trim(mesh_file)
        write(*,*) "    Initializing surface mesh from file:", mesh_file

        ! Determine the type of mesh file
        loc = index(mesh_file, '.')
        extension = mesh_file(loc:len(mesh_file))

        ! Load vtk
        if (extension .eq. '.vtk') then
            call load_vtk(mesh_file, t%N_vert, t%N_panel, t%vertices, t%panels)
        end if

        ! Display mesh info
        write(*,*) "    Surface mesh has", t%N_vert, "vertices and", t%N_panel, "panels."
    
    end subroutine surface_mesh_initialize


    subroutine surface_mesh_output_results(t)

        implicit none

        class(surface_mesh),intent(out) :: t
        character(len=:),allocatable :: output_file

        ! Get filename
        call t%input_json%get('output.file', output_file)

        ! Write out data
        call write_surface_vtk(output_file, t%vertices, t%panels)
    
    end subroutine surface_mesh_output_results

    
end module mesh