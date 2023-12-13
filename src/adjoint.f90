! module for calculating the gradient of CF, CM with respect to mesh points via the adjoint method
module adjoint_mod

    use linked_list_mod
    use math_mod
    use helpers_mod
    use base_geom_mod
    use mesh_mod
    use surface_mesh_mod

    implicit none

    type adjoint

    real,dimension(:),allocatable :: X_beta  ! list of x y z values of all mesh points (design variables)

    contains
            procedure :: init => adjoint_init
            procedure :: get_X_beta => adjoint_get_X_beta ! puts all x y z values of mesh points in a list

    end type adjoint

    subroutine adjoint_init(this, body)
        class(adjoint),intent(inout) :: this
        type(surface_mesh),intent(in) :: body

        call this%get_X_beta(body)

    end subroutine adjoint_init


    subroutine adjoint_get_X_beta(this, body)
        implicit none

        class(adjoint),intent(inout) :: this
        type(surface_mesh),intent(in) :: body
        integer :: i, j, N_verts

        N_verts = body%N_verts
        ! build design variable vector X_beta
        allocate(this%X_beta(N_verts*3))
        do i=1,3
            do j=1,N_verts
                this%X_beta(j + (i-1)*N_verts) = body%vertices(j)%loc(i)
            end do
        end do


    end subroutine adjoint_get_X_beta



end module adjoint_mod