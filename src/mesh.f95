! Types and subroutines for meshes
module mesh

    use json
    use json_xtnsn
    use vtk
    use geometry
    use adt

    implicit none


    type surface_mesh

        integer :: N_vert, N_panel
        real,allocatable,dimension(:,:) :: vertices
        type(vertex),allocatable,dimension(:) :: vertex_objects
        type(panel),allocatable,dimension(:) :: panels
        character(len=:),allocatable :: mesh_file
        type(alternating_digital_tree) :: vertex_tree

        contains

            procedure :: initialize => surface_mesh_initialize
            procedure :: output_results => surface_mesh_output_results
            procedure :: calc_inv_panel_vertex_mapping => surface_mesh_calc_inverse_panel_vertex_mapping

    end type surface_mesh


    type cart_volume_mesh

    end type cart_volume_mesh

    
contains


    subroutine surface_mesh_initialize(t, settings)

        implicit none

        class(surface_mesh),intent(inout) :: t
        type(json_value),pointer,intent(in) :: settings
        character(len=:),allocatable :: extension
        integer :: loc, i

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

        ! Load vertices into object array
        allocate(t%vertex_objects(t%N_vert))
        do i=1,t%N_vert
            t%vertex_objects(i)%p = t%vertices(i,:)
        end do

        ! Display mesh info
        write(*,*) "    Surface mesh has", t%N_vert, "vertices and", t%N_panel, "panels."

        ! Set bounds of alternating digital tree
        t%vertex_tree%p_min = minval(t%vertices, 1)
        t%vertex_tree%p_max = maxval(t%vertices, 1)

        ! Load vertices into alternating digital tree
        write(*,*)
        write(*,*) "    Loading vertices into ADT..."
        do i=1,t%N_vert
            call t%vertex_tree%add(t%vertices(i,:))
        end do

        ! Run inverse vertex mapping
        call t%calc_inv_panel_vertex_mapping()
    
    end subroutine surface_mesh_initialize


    subroutine surface_mesh_calc_inverse_panel_vertex_mapping(t)

        implicit none

        class (surface_mesh),intent(inout) :: t
        integer :: i,j

        ! Loop through panels and add their indices to each vertex's list of neighboring panels
        do i=1,t%N_panel

            ! Vertex 1
            do j=1,20
                if (t%vertex_objects(t%panels(i)%i1)%neighboring_panels(j) == -1) then
                    t%vertex_objects(i)%neighboring_panels(j) = i
                    exit
                end if
            end do

            ! Vertex 2
            do j=1,20
                if (t%vertex_objects(t%panels(i)%i2)%neighboring_panels(j) == -1) then
                    t%vertex_objects(i)%neighboring_panels(j) = i
                    exit
                end if
            end do

            ! Vertex 3
            do j=1,20
                write(*,*) t%panels(i)%i3
                if (t%vertex_objects(t%panels(i)%i3)%neighboring_panels(j) == -1) then
                    t%vertex_objects(i)%neighboring_panels(j) = i
                    exit
                end if
            end do

        end do

        do i=1,t%N_vert
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*) t%vertex_objects(i)%neighboring_panels
        end do
    
    end subroutine surface_mesh_calc_inverse_panel_vertex_mapping


    subroutine surface_mesh_output_results(t, output_file)

        implicit none

        class(surface_mesh),intent(inout) :: t
        character(len=:),allocatable,intent(in) :: output_file

        ! Write out data
        call write_surface_vtk(output_file, t%vertices, t%panels)
    
    end subroutine surface_mesh_output_results

    
end module mesh