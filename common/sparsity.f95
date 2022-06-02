! Subroutines for rearranging matrices for desired sparsity
module sparsity_mod

    implicit none
    
contains


subroutine rearrange(N, A, P_row, P_col)
    ! Gives the indices which will rearrange A for better sparsity

    implicit none

    integer,intent(in) :: N
    real,dimension(N,N),intent(in) :: A
    integer,dimension(N),intent(out) :: P_row, P_col

    integer :: i

    ! Initialize indices
    do i=1,N
        P_row(i) = i
        P_col(i) = i
    end do
    
end subroutine rearrange

    
end module sparsity_mod