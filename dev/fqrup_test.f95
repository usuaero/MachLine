program fqrup_test

    use linalg_mod

    implicit none

    integer,parameter :: N = 4
    real,dimension(N,N) :: A, A_copy
    real,dimension(N) :: b, b_copy
    real,dimension(:),allocatable :: x
    integer :: i

    ! Set up A
    A(1,1) = 0.
    A(1,2) = 2.
    A(1,3) = 3.
    A(1,4) = 3.

    A(2,1) = 2.
    A(2,2) = -1.
    A(2,3) = 0.
    A(2,4) = 3.

    A(3,1) = 0.
    A(3,2) = 1.
    A(3,3) = 1.
    A(3,4) = 3.

    A(4,1) = 0.
    A(4,2) = 0.
    A(4,3) = -1.
    A(4,4) = -4.
    A_copy = A

    ! Set up b
    b(1) = 1.
    b(2) = 1.
    b(3) = 1.
    b(4) = 1.
    b_copy = b

    ! Solve
    call QR_fast_givens_solve_upper_pentagonal(N, A, b, x)

    write(*,*) "x"
    write(*,*) x

    write(*,*) "A tilde"
    do i=1,N
        write(*,*) A(i,:)
    end do

    write(*,*) "b tilde"
    write(*,*) b

    write(*,*) "Ax"
    write(*,*) matmul(A_copy, x)

    write(*,*) "b"
    write(*,*) b_copy

    stop
    
end program fqrup_test