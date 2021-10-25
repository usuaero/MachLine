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


        contains

            procedure :: init => panel_solver_init
            procedure :: solve => panel_solver_solve

    end type panel_solver


contains


    subroutine panel_solver_init(this, settings)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings


    end subroutine panel_solver_init


    subroutine panel_solver_solve(this, body_mesh, freestream_flow)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body_mesh
        type(flow),intent(in) :: freestream_flow
        integer :: i

        write(*,*)
        write(*,'(a)',advance='no') "     Running linear solver..."

        ! Set source strengths

        ! Calculate influence matrices

        ! Assemble linear system

        ! Solve

        write(*,*) "Done."

    end subroutine panel_solver_solve


end module panel_solver_mod