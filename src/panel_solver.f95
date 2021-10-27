module panel_solver_mod

    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use vertex_mod
    use surface_mesh_mod
    use flow_mod
    use math_mod

    implicit none


    type panel_solver


        character(len=:),allocatable :: type

        contains

            procedure :: init => panel_solver_init
            procedure :: solve => panel_solver_solve
            procedure :: set_source_strengths => panel_solver_set_source_strengths

    end type panel_solver


contains


    subroutine panel_solver_init(this, settings, body_mesh)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        type(surface_mesh),intent(inout) :: body_mesh
        real :: control_point_offset

        ! Get settings
        call json_xtnsn_get(settings, 'type', this%type, 'indirect')

        write(*,*)
        write(*,'(a)',advance='no') "     Placing control points..."

        ! Place control points
        call json_xtnsn_get(settings, 'control_point_offset', control_point_offset, 1e-5)
        call body_mesh%place_control_points(control_point_offset)
        write(*,*) "Done."


    end subroutine panel_solver_init


    subroutine panel_solver_solve(this, body_mesh, freestream_flow)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body_mesh
        type(flow),intent(in) :: freestream_flow

        integer :: i, j, k
        real,dimension(:),allocatable :: influence
        integer,dimension(:),allocatable :: vertex_indices
        real,dimension(:,:),allocatable :: A
        real,dimension(:),allocatable :: b

        write(*,*)
        write(*,'(a)') "     Running linear solver"

        if (this%type .eq. 'indirect') then ! Morino formulation

            ! Allocate linear system
            allocate(A(body_mesh%N_verts, body_mesh%N_verts), source=0.)
            allocate(b(body_mesh%N_verts), source=0.)

            ! Set source strengths
            call this%set_source_strengths(freestream_flow, body_mesh)

            ! Calculate source influences
            write(*,*)
            write(*,'(a)',advance='no') "     Calculating source influences..."
            do i=1,body_mesh%N_verts
                do j=1,body_mesh%N_panels

                    ! Get source influence
                    influence = body_mesh%panels(j)%get_source_potential(body_mesh%control_points(i,:), vertex_indices)

                    ! Add to RHS
                    if (source_order .eq. 0) then
                        b(i) = b(i) + influence(1)*body_mesh%sigma(j)
                    end if

                end do
            end do
            write(*,*) "Done."

            ! Calculate doublet influences
            write(*,*)
            write(*,'(a)',advance='no') "     Calculating doublet influences..."
            do i=1,body_mesh%N_verts
                do j=1,body_mesh%N_panels

                    ! Get source influence
                    influence = body_mesh%panels(j)%get_doublet_potential(body_mesh%control_points(i,:), vertex_indices)

                    ! Add to LHS
                    if (doublet_order .eq. 1) then
                        do k=1,size(vertex_indices)
                            A(i,vertex_indices(k)) = A(i,vertex_indices(k)) + influence(k)
                        end do
                    end if

                end do
            end do
            write(*,*) "Done."

            ! Solve
            write(*,*)
            write(*,'(a)',advance='no') "     Solving linear system..."
            allocate(body_mesh%mu(body_mesh%N_verts))
            call lu_solve(body_mesh%N_verts, A, b, body_mesh%mu)
            write(*,*) "Done."

        end if

    end subroutine panel_solver_solve


    subroutine panel_solver_set_source_strengths(this, freestream_flow, body_mesh)
        ! Determines the source strengths for the given situation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(flow),intent(in) :: freestream_flow
        type(surface_mesh),intent(inout) :: body_mesh

        integer :: i

        ! Constant source distribution
        if (source_order .eq. 0.) then

            ! Allocate source strength array
            allocate(body_mesh%sigma(body_mesh%N_panels))

            ! Loop through panels
            do i=1,body_mesh%N_panels
                body_mesh%sigma(i) = -inner(body_mesh%panels(i)%normal, freestream_flow%V_inf)
            end do

        end if


    end subroutine panel_solver_set_source_strengths


end module panel_solver_mod