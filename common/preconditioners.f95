module preconditioners_mod

    implicit none
    
contains


subroutine diagonal_preconditioner(N, A, b, A_p, b_p)
    ! Applies the preconditioning M^-1 A x = M^-1 b where
    ! M is the diagonal elements of A

    implicit none

    integer,intent(in) :: N
    real,dimension(N,N),intent(in) :: A
    real,dimension(N),intent(in) :: b
    real,dimension(:,:),allocatable,intent(out) :: A_p
    real,dimension(:),allocatable,intent(out) :: b_p

    real,dimension(N) :: A_ii_inv
    integer :: i, j

    ! Get preconditioning matrix
    do i=1,N
        A_ii_inv = 1./A(i,i)
    end do

    ! Allocate
    allocate(A_p(N,N))
    allocate(b_p(N))

    ! Apply preconditioning
    do j=1,N
        do i=1,N
            A_p(i,j) = A_ii_inv(i)*A(i,j)
        end do

        b_p(j) = b(j)*A_ii_inv(j)
    end do
    
end subroutine diagonal_preconditioner

    
end module preconditioners_mod