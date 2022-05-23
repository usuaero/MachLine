! Linear algebra subroutines
module linalg_mod

    implicit none
    
contains


subroutine matinv(n, a, ai)
      ! This sobroutine inverts a matrix "a" and returns the inverse in "ai"
      ! n  - Input by user, an integer specifying the size of the matrix to be inverted.
      ! a  - Input by user, an n by n real array containing the matrix to be inverted.
      ! ai - Returned by subroutine, an n by n real array containing the inverted matrix.
      ! d  - Work array, an n by 2n real array used by the subroutine.
      ! io - Work array, a 1-dimensional integer array of length n used by the subroutine.

      ! NEVER INVERT A MARTIX EXPLICITLY!
      ! Unless you know what you're doing. Odds are you may not, so be careful.

      implicit none

      integer :: n,i,j,k,m,itmp
      real :: a(n,n),ai(n,n),tmp,r
      real,allocatable,dimension(:,:) :: d
      integer,allocatable,dimension(:) :: io

      allocate(d(n,2*n))
      allocate(io(n))

      d(:,:) = 0.0
      io(:) = 0

!     Fill in the "io" and "d" matrix.
!     ********************************
      do i=1,n
         io(i)=i
      end do
      do i=1,n
         do j=1,n
            d(i,j)=a(i,j)
            if(i.eq.j)then
               d(i,n+j)=1.
            else
               d(i,n+j)=0.
            endif
         end do
      end do

!     Scaling
!     *******
      do i=1,n
         m=1
         do k=2,n
            if(abs(d(i,k)).gt.abs(d(i,m))) m=k
         end do
         tmp=d(i,m)
         do k=1,2*n
            d(i,k)=d(i,k)/tmp
         end do
      end do

!     Lower Elimination
!     *****************
      do i=1,n-1
!        Pivoting
!        ********
         m=i
         do j=i+1,n
            if(abs(d(io(j),i)).gt.abs(d(io(m),i))) m=j
         end do
         itmp=io(m)
         io(m)=io(i)
         io(i)=itmp
!        Scale the Pivot element to unity
!        ********************************
         r=d(io(i),i)
         do k=1,2*n
            d(io(i),k)=d(io(i),k)/r
         end do
!        ********************************
         do j=i+1,n
            r=d(io(j),i)
            do k=1,2*n
               d(io(j),k)=d(io(j),k)-r*d(io(i),k)
            end do
         end do
      end do

!     Upper Elimination
!     *****************
      r=d(io(n),n)
      do k=1,2*n
         d(io(n),k)=d(io(n),k)/r
      end do
      do i=n-1,1,-1
         do j=i+1,n
            r=d(io(i),j)
            do k=1,2*n
               d(io(i),k)=d(io(i),k)-r*d(io(j),k)
            end do
         end do
      end do

!     Fill Out "ai" matrix
      do i=1,n
         do j=1,n
            ai(i,j)=d(io(i),n+j)
         end do
      end do

      ! Cleanup
      deallocate(d)
      deallocate(io)

end subroutine matinv


subroutine lu_solve(n, A, b, x)
  ! Solves a general [A]x=b on an nxn matrix
  ! This replaces A (in place) with its LU decomposition (permuted row-wise)

    implicit none

    integer,intent(in) :: n
    real,dimension(n,n),intent(inout) :: A
    real,dimension(n),intent(in) :: b
    real,dimension(:),allocatable,intent(out) :: x

    integer,allocatable,dimension(:) :: indx
    integer :: D, info

    allocate(indx(n))

    ! Compute decomposition
    call lu_decomp(A, n, indx, D, info)

    ! if the matrix is nonsingular, then backsolve to find X
    if (info == 1) then
        write(*,*) 'Subroutine lu_decomp() failed. The given matrix is singular (i.e. no unique solution). Quitting...'
        stop
    else
        call lu_back_sub(A, n, indx, b, x)
    end if

    ! Cleanup
    deallocate(indx)

end subroutine lu_solve


!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*    improved for F95 by Cory Goates, Logan, UT, USA  *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     *
!*******************************************************


subroutine lu_decomp(A, N, indx, D, code)
  ! Given an N x N matrix A, this routine replaces it by the LU
  ! decomposition of a rowwise permutation of itself. A and N  
  ! are input. indx is an output vector which records the row  
  ! permutation effected by the partial pivoting; D is output  
  ! as -1 or 1, depending on whether the number of row inter-  
  ! changes was even or odd, respectively. This routine is used
  ! in combination with lu_back_sub to solve linear equations or to 
  ! invert a matrix. Return code is 1 if matrix is singular.  

  implicit none

  real,dimension(N,N),intent(inout) :: A
  integer,intent(in) :: N
  integer,dimension(N),intent(out) :: indx
  integer,intent(out) :: code, D

  real,dimension(N) :: vv
  real,parameter :: tiny=1.5e-20
  integer :: i, j, k, imax
  real :: amax, dum, s

  ! Initialize
  D = 1
  code = 0
  imax = 0

  ! Loop over rows to get implicit scaling information
  do i=1,N

    ! Get largest element in this row
    amax=0.0
    do j=1,N
      if (abs(A(i,j)) > amax) then
        amax = abs(A(i,j))
      end if
    end do

    ! Check the largest element in this row is nonzero
    if (amax <= tiny) then
      code = 1 ! Singular matrix
    end if

    ! Store implicit scaling
    vv(i) = 1.0 / amax

  end do

  ! Check for singular matrix
  if (code == 1) return

  ! Loop over columns of Crout's method
  do j=1,N

    do i=1,j-1
      
      s = A(i,j)

      do k=1,i-1
        s = s - A(i,k)*A(k,j)
      end do

      A(i,j) = s

    end do

    ! Initialize search for largest pivot element
    amax = 0.0
    do i=j,N
  
      s = A(i,j)
      do k=1,j-1
        s = s - A(i,k)*A(k,j)
      end do
      A(i,j) = s
  
      ! Determine figure of merit for the pivot
      dum = vv(i)*abs(s)
      if (dum >= amax) then
        imax = i
        amax = dum
      end if
  
    end do

    ! Figure out if we need to interchange rows
    if (j /= imax) then
  
      ! Perform interchange
      do k=1,N
        dum = A(imax,k)
        A(imax,k) = A(j,k)
        A(j,k) = dum
      end do
  
      ! Update the sign of D since a row interchange has occurred
      D = -D
  
      ! Interchange the implicit scaling factor
      vv(imax) = vv(j)
  
    end if

    ! Store pivoting
    indx(j) = imax

    ! Divide by pivot element
    if (j /= N) then
      dum = 1.0 / A(j,j)
      do i=j+1,N
        A(i,j) = A(i,j)*dum
      end do
    end if

  end do

end subroutine lu_decomp


subroutine lu_back_sub(A, N, indx, b, x)
  ! Solves the set of N linear equations Ax = b.  Here A is     
  ! input, not as the matrix A but rather as its LU decomposition, 
  ! determined by the routine LUDCMP. indx is input as the permuta-
  ! tion vector returned by LUDCMP. b is input as the right-hand   
  ! side vector b. The solution vector is x. A, N, b and
  ! indx are not modified by this routine and can be used for suc- 
  ! cessive calls with different right-hand sides. This routine is 
  ! also efficient for plain matrix inversion.                     

  implicit none

  integer,intent(in) :: N
  real,dimension(N,N),intent(in) :: A
  real,dimension(N),intent(in) :: b
  integer,dimension(N),intent(in) :: indx
  real,dimension(:),allocatable,intent(out) :: x

  real :: sum
  integer :: ii,i,j,ll

  ! Initialize solution
  allocate(x, source=b)

  ! Set tracker to ignore leading zeros in b
  ii = 0

  ! Forward substitution
  do i=1,N

    ! Untangle pivoting
    ll = indx(i)
    sum = x(ll)
    x(ll) = x(i)

    ! If a nonzero element of b has already been encountered
    if (ii /= 0) then
      do J=ii,i-1
        sum = sum - A(i,J)*x(J)
      end do

    ! Check for first nonzero element of b
    else if(sum /= 0.0) then
      ii = i
    end if

    x(i) = sum

  end do

  ! Back substitution
  do i=N,1,-1
    sum = x(i)
    do j=i+1,N
      sum = sum - A(i,j)*x(j)
    end do
    x(i) = sum / A(i,i)
  end do

end subroutine lu_back_sub


subroutine quadratic_fit(pts, a, b, c)
  ! Fits a parabola through three specified   
  ! points and returns the coefficients a, b, c defining this parabola 
  ! according to the equation y = a * x**2 + b * x + c                 
  !            pts = list of three (x, y) points                       
  !   Outputs: a, b, c = quadratic coefficients                        

  implicit none

  real,dimension(3, 2),intent(in) :: pts
  real,intent(out) :: a, b, c

  integer :: i
  real,dimension(3,3) :: m
  real,dimension(:),allocatable :: coeff

  do i = 1, 3
    m(i,1) = pts(i,1)**2
    m(i,2) = pts(i,1)
    m(i,3) = 1.0
  end do

  call lu_solve(3, m, pts(:,2), coeff)

  a = coeff(1)
  b = coeff(2)
  c = coeff(3)

end subroutine quadratic_fit


subroutine lu_solve_replace(n, A, b, x)
  ! Solves a general [A]x=b on an nxn matrix without overwriting A

  implicit none

  integer,intent(in) :: n
  real,dimension(n,n),intent(in) :: A
  real,dimension(n),intent(in) :: b
  real,dimension(:),allocatable,intent(out) :: x

  real,dimension(n,n) :: A_copy

  ! Create copy of A
  A_copy = A

  ! Solve
  call lu_solve(n, A_copy, b, x)

end subroutine lu_solve_replace


subroutine decompose_blocks(N, A, N_blocks, block_size, N_last, ind_P, i_start_block, i_end_block)
  ! Decomposes the diagonal blocks of A using LU decomposition
  ! N is the size of the system
  ! A is the system matrix; block diagonals will be replaced with their LU decompositions
  ! block_size is the desired size of each block
  ! N_last is the number of elements in the last block
  ! ind_P is a matrix of the block permutation indices
  ! i_start_block is a vector of the start indices for each block
  ! i_end_block is a vector of the end indices for each block

  implicit none

  integer,intent(in) :: N
  real,dimension(N,N),intent(inout) :: A
  integer,intent(in) :: N_blocks, block_size
  integer,intent(out) :: N_last
  integer,dimension(:,:),allocatable,intent(out) :: ind_P
  integer,dimension(:),allocatable,intent(out) :: i_start_block, i_end_block

  integer :: i, info, D

  ! Allocate start and end indices
  allocate(i_start_block(N_blocks))
  allocate(i_end_block(N_blocks))

  ! Allocate permutation incides
  allocate(ind_P(block_size, N_blocks))

  ! Decompose blocks
  do i=1,N_blocks

    ! Last block
    if (i == N_blocks) then

      ! Determine start and end indices of this block
      i_start_block(i) = (i-1)*block_size + 1
      i_end_block(i) = N
      N_last = i_end_block(i)-i_start_block(i)+1

      ! Decompose
      call lu_decomp(A(i_start_block(i):N,i_start_block(i):N), N_last, ind_P(1:N_last,i), D, info)

      ! Check for singular submatrix
      if (info == 1) then
          write(*,*) 'Subroutine lu_decomp() failed on subblock ', i, '. Quitting...'
          stop
      end if

    ! Not last block
    else

      ! Determine start and end indices of this block
      i_start_block(i) = (i-1)*block_size + 1
      i_end_block(i) = i*block_size

      ! Decompose
      call lu_decomp(A(i_start_block(i):i_end_block(i),i_start_block(i):i_end_block(i)), block_size, ind_P(:,i), D, info)

      ! Check for singular submatrix
      if (info == 1) then
          write(*,*) 'Subroutine lu_decomp() failed on subblock ', i, '. Quitting...'
          stop
      end if

    end if

  end do

end subroutine decompose_blocks


subroutine block_sor(N, A, b, block_size, tol, rel, max_iterations, verbose, x)
  ! Iteratively solves the [A]x=b system using block Successive Overrelaxation
  ! N is the size of the system
  ! A is the system matrix; block diagonals will be replaced with their LU decompositions
  ! b is the RHS vector
  ! block_size is the desired size of each block
  ! tol is the convergence tolerance between iterations
  ! rel is a relaxation factor between 0 and 2
  ! x is the solution

  implicit none

  integer,intent(in) :: N, block_size, max_iterations
  real,dimension(N,N),intent(inout) :: A
  real,dimension(N),intent(in) :: b
  real,intent(in) :: tol, rel
  logical,intent(in) :: verbose
  real,dimension(:),allocatable,intent(out) :: x

  real :: err, dx
  real,dimension(N) :: x_new
  real,dimension(N,N) :: A_copy
  real,dimension(block_size) :: bi
  real,dimension(:),allocatable :: xi
  integer :: i, N_blocks, r, N_last, iteration
  integer,dimension(:),allocatable :: i_start_block, i_end_block
  integer,dimension(:,:),allocatable :: ind_P

  ! Check relaxation
  if (rel < 0. .or. rel > 2.) then
    write(*,*) "!!! Relaxation for SOR must be between 0 and 2 (exclusive). Quitting..."
    stop
  end if

  ! Give initial error estimate
  err = tol + 1.

  ! Make a copy of A for calculating the error
  A_copy = A

  ! Initialize solution vector
  allocate(x(N), source=0.)

  ! Calculate number of blocks
  N_blocks = N/block_size
  r = modulo(N, block_size)
  if (r > 0) then
    N_blocks = N_blocks + 1
  end if

  ! Decompose blocks
  if (verbose) then
    write(*,*)
    write(*,'(a)',advance='no') "         Decomposing blocks..."
  end if
  call decompose_blocks(N, A, N_blocks, block_size, N_last, ind_P, i_start_block, i_end_block)
  if (verbose) write(*,*) "Done."

  ! Progress
  if (verbose) then
    write(*,*)
    write(*,*) "        Running Block SOR..."
    write(*,*) "        Iteration      ||dx||            ||err||"
    write(*,*) "        ----------------------------------------"
  end if

  ! Iterate
  iteration = 0
  do while(err >= tol .and. iteration < max_iterations)

    iteration = iteration + 1

    ! Invert blocks
    do i=1,N_blocks

      ! Last block
      if (i == N_blocks) then

        ! Calculate new RHS vector
        bi(1:N_last) = b(i_start_block(i):N) - matmul(A(i_start_block(i):N,1:i_start_block(i)-1), x_new(1:i_start_block(i)-1))

        ! Solve
        call lu_back_sub(A(i_start_block(i):N,i_start_block(i):N), N_last, ind_P(1:N_last,i), bi, xi)

        ! Store
        x_new(i_start_block(i):N) = (1.-rel)*x_new(i_start_block(i):N) + rel*xi

      ! Not last block
      else

        ! Calculate new RHS vector
        bi = b(i_start_block(i):i_end_block(i))
        bi = bi - matmul(A(i_start_block(i):i_end_block(i),1:i_start_block(i)-1), x_new(1:i_start_block(i)-1))
        bi = bi - matmul(A(i_start_block(i):i_end_block(i),i_end_block(i)+1:N), x(i_end_block(i)+1:N))

        ! Solve
        call lu_back_sub(A(i_start_block(i):i_end_block(i),i_start_block(i):i_end_block(i)), block_size, ind_P(:,i), bi, xi)

        ! Store
        x_new(i_start_block(i):i_end_block(i)) = (1.-rel)*x_new(i_start_block(i):i_end_block(i)) + rel*xi

      end if
    end do

    ! Calculate error
    err = norm2(matmul(A_copy, x_new) - b)
    dx = norm2(x-x_new)

    ! Output progress
    if (modulo(iteration, 50) == 0) then
      if (verbose) write(*,'(i18, ES15.3, ES15.3)') iteration, dx, err
    end if

    ! Update
    x = x_new

  end do

  ! Final iteration count and error
  if (verbose) write(*,'(i18, ES15.3, ES15.3)') iteration, dx, err

end subroutine block_sor


subroutine block_sor_adaptive(N, A, b, block_size, tol, max_iterations, verbose, x)
  ! Iteratively solves the [A]x=b system using block Successive Overrelaxation
  ! N is the size of the system
  ! A is the system matrix; block diagonals will be replaced with their LU decompositions
  ! b is the RHS vector
  ! block_size is the desired size of each block
  ! tol is the convergence tolerance between iterations
  ! verbose
  ! x is the solution

  implicit none

  integer,intent(in) :: N, block_size, max_iterations
  real,dimension(N,N),intent(inout) :: A
  real,dimension(N),intent(in) :: b
  real,intent(in) :: tol
  logical,intent(in) :: verbose
  real,dimension(:),allocatable,intent(out) :: x

  real :: err, rel, dx
  real,dimension(N) :: x_new
  real,dimension(N,N) :: A_copy
  real,dimension(block_size) :: bi
  real,dimension(:),allocatable :: xi
  integer :: i, N_blocks, r, N_last, iteration, step, start, end
  integer,dimension(:),allocatable :: i_start_block, i_end_block
  integer,dimension(:,:),allocatable :: ind_P

  ! Give initial error estimate and relaxation
  err = tol + 1.
  rel = 0.1

  ! Initialize alterantion
  step = -1

  ! Make a copy of A for calculating the error
  A_copy = A

  ! Initialize solution vector
  allocate(x(N), source=0.)

  ! Calculate number of blocks
  N_blocks = N/block_size
  r = modulo(N, block_size)
  if (r > 0) then
    N_blocks = N_blocks + 1
  end if

  ! Decompose blocks
  if (verbose) then
    write(*,*)
    write(*,'(a)',advance='no') "         Decomposing blocks..."
  end if
  call decompose_blocks(N, A, N_blocks, block_size, N_last, ind_P, i_start_block, i_end_block)
  if (verbose) write(*,*) "Done."

  ! Progress
  if (verbose) then
    write(*,*)
    write(*,*) "        Running Adaptive Block SOR..."
    write(*,*) "        Iteration      ||dx||         Relaxation        ||err||"
    write(*,*) "        -------------------------------------------------------"
  end if

  ! Iterate
  iteration = 0
  do while(err >= tol .and. iteration < max_iterations)

    ! Update iteration number
    iteration = iteration + 1

    ! Alternate direction
    if (step == 1) then
      start = N_blocks
      end = 1
      step = -1
    else
      start = 1
      end = N_blocks
      step = 1
    end if

    ! Perform SOR on blocks
    do i=1,N_blocks

      ! Last block
      if (i == N_blocks) then

        ! Calculate new RHS vector
        if (step == 1) then
          bi(1:N_last) = b(i_start_block(i):N) - matmul(A(i_start_block(i):N,1:i_start_block(i)-1), x_new(1:i_start_block(i)-1))
        else
          bi(1:N_last) = b(i_start_block(i):N) - matmul(A(i_start_block(i):N,1:i_start_block(i)-1), x(1:i_start_block(i)-1))
        end if

        ! Solve
        call lu_back_sub(A(i_start_block(i):N,i_start_block(i):N), N_last, ind_P(1:N_last,i), bi, xi)

        ! Store
        x_new(i_start_block(i):N) = (1.-rel)*x(i_start_block(i):N) + rel*xi

      ! Not last block
      else

        ! Calculate new RHS vector
        bi = b(i_start_block(i):i_end_block(i))
        if (step == 1) then
          bi = bi - matmul(A(i_start_block(i):i_end_block(i),1:i_start_block(i)-1), x_new(1:i_start_block(i)-1))
          bi = bi - matmul(A(i_start_block(i):i_end_block(i),i_end_block(i)+1:N), x(i_end_block(i)+1:N))
        else
          bi = bi - matmul(A(i_start_block(i):i_end_block(i),1:i_start_block(i)-1), x(1:i_start_block(i)-1))
          bi = bi - matmul(A(i_start_block(i):i_end_block(i),i_end_block(i)+1:N), x_new(i_end_block(i)+1:N))
        end if

        ! Solve
        call lu_back_sub(A(i_start_block(i):i_end_block(i),i_start_block(i):i_end_block(i)), block_size, ind_P(:,i), bi, xi)

        ! Store
        x_new(i_start_block(i):i_end_block(i)) = (1.-rel)*x(i_start_block(i):i_end_block(i)) + rel*xi

      end if
    end do

    ! Calculate error
    err = norm2(matmul(A_copy, x_new) - b)
    dx = norm2(x - x_new)

    ! Output progress
    if (verbose .and. modulo(iteration, 50) == 0) then
      write(*,'(i18, ES15.3, ES15.3, ES15.3)') iteration, dx, rel, err
    end if

    ! Update
    x = x_new

    ! Aim for a step size that's of the same order of magnitude as the error
    rel = rel - 0.1*(log10(dx/err))

    ! Nowhere should we go above 2
    if (rel > 2.) rel = 2.

    ! Make sure the relaxation factor stays positive
    if (rel < 0.) rel = 0.01

  end do

  ! Final iteration count and error
  if (verbose) write(*,'(i18, ES15.3, ES15.3, ES15.3)') iteration, dx, rel, err

end subroutine block_sor_adaptive


subroutine iterative_solve(method, N, A, b, block_size, tol, rel, max_iterations, verbose, x)
  ! Iteratively solves the [A]x=b system using the specified block-iterative method
  ! Will alternate directions through the blocks on each iteration
  ! N is the size of the system
  ! A is the system matrix; block diagonals will be replaced with their LU decompositions
  ! b is the RHS vector
  ! block_size is the desired size of each block
  ! method is "BJAC", "BSOR", "ORBJ", or "ABSOR"
  ! tol is the convergence tolerance between iterations
  ! rel is a relaxation factor between 0 and 2; it is the initial factor for adaptive methods
  ! max_iterations
  ! verbose
  ! x is the solution

  ! BJAC method is standard block Jacobi iteration with a relaxation factor
  ! ORBJ method is block Jacobi with a relaxation factor calculated to minimize the error at each step; rel is not used for this method

  implicit none

  character(len=*),intent(in) :: method
  integer,intent(in) :: N, block_size, max_iterations
  real,dimension(N,N),intent(inout) :: A
  real,dimension(N),intent(in) :: b
  real,intent(inout) :: tol, rel
  logical,intent(in) :: verbose
  real,dimension(:),allocatable,intent(out) :: x

  real :: err, dx, vkTvk, bTvk, vkp1Tvkp1, vkp1Tvk, bTvkp1
  real,dimension(N) :: x_new, vk, vkp1
  real,dimension(N,N) :: A_copy
  real,dimension(block_size) :: bi
  real,dimension(:),allocatable :: xi
  integer :: i, N_blocks, r, N_last, iteration, step, start, end, unit
  integer,dimension(:),allocatable :: i_start_block, i_end_block
  integer,dimension(:,:),allocatable :: ind_P
  logical :: jacobi

  ! Give initial error estimate
  err = tol + 1.

  ! Initialize alterantion
  step = -1

  ! Check if we're using a Jacobi method
  jacobi = method == "BJAC" .or. method == "ORBJ"

  ! Initialize optimal relaxation determination for ORBJ method
  if (method == "ORBJ") then
    vk = 0.
    vkTvk = 0.
    bTvk = 0.
  end if

  ! Make a copy of A for calculating the error
  A_copy = A

  ! Initialize solution vector
  allocate(x(N), source=0.)

  ! Calculate number of blocks
  N_blocks = N/block_size
  r = modulo(N, block_size)
  if (r > 0) then
    N_blocks = N_blocks + 1
  end if

  ! Decompose blocks
  call decompose_blocks(N, A, N_blocks, block_size, N_last, ind_P, i_start_block, i_end_block)

  ! Write out header
  if (verbose) then
    open(newunit=unit, file='iterative_solver_prog.csv')
    write(unit,*) "method,", method
    write(unit,*) "iteration,||dx||,||err||,relaxation"
  end if

  ! Iterate
  iteration = 0
  do while(err >= tol .and. iteration < max_iterations)

    ! Update iteration number
    iteration = iteration + 1

    ! Alternate direction
    if (step == 1) then
      start = N_blocks
      end = 1
      step = -1
    else
      start = 1
      end = N_blocks
      step = 1
    end if

    ! Perform Jacobi iteration on blocks
    do i=1,N_blocks

      ! Last block
      if (i == N_blocks) then

        ! Calculate new RHS vector
        if (jacobi) then
          bi(1:N_last) = b(i_start_block(i):N) - matmul(A(i_start_block(i):N,1:i_start_block(i)-1), x(1:i_start_block(i)-1))

        else
          if (step == 1) then
            bi(1:N_last) = b(i_start_block(i):N) - matmul(A(i_start_block(i):N,1:i_start_block(i)-1), x_new(1:i_start_block(i)-1))
          else
            bi(1:N_last) = b(i_start_block(i):N) - matmul(A(i_start_block(i):N,1:i_start_block(i)-1), x(1:i_start_block(i)-1))
          end if
        end if

        ! Solve
        call lu_back_sub(A(i_start_block(i):N,i_start_block(i):N), N_last, ind_P(1:N_last,i), bi, xi)

        ! Store
        if (jacobi) then
          x_new(i_start_block(i):N) = xi
        else
          x_new(i_start_block(i):N) = (1.-rel)*x(i_start_block(i):N) + rel*xi
        end if

      ! Not last block
      else

        ! Calculate new RHS vector
        if (jacobi) then
          bi = b(i_start_block(i):i_end_block(i))
          bi = bi - matmul(A(i_start_block(i):i_end_block(i),1:i_start_block(i)-1), x(1:i_start_block(i)-1))
          bi = bi - matmul(A(i_start_block(i):i_end_block(i),i_end_block(i)+1:N), x(i_end_block(i)+1:N))

        else
          if (step == 1) then
            bi = bi - matmul(A(i_start_block(i):i_end_block(i),1:i_start_block(i)-1), x_new(1:i_start_block(i)-1))
            bi = bi - matmul(A(i_start_block(i):i_end_block(i),i_end_block(i)+1:N), x(i_end_block(i)+1:N))
          else
            bi = bi - matmul(A(i_start_block(i):i_end_block(i),1:i_start_block(i)-1), x(1:i_start_block(i)-1))
            bi = bi - matmul(A(i_start_block(i):i_end_block(i),i_end_block(i)+1:N), x_new(i_end_block(i)+1:N))
          end if
        end if

        ! Solve
        call lu_back_sub(A(i_start_block(i):i_end_block(i),i_start_block(i):i_end_block(i)), block_size, ind_P(:,i), bi, xi)

        ! Store
        if (jacobi) then
          x_new(i_start_block(i):i_end_block(i)) = xi
        else
          x_new(i_start_block(i):i_end_block(i)) = (1.-rel)*x(i_start_block(i):i_end_block(i)) + rel*xi
        end if

      end if
    end do

    ! Calculate optimal relaxation for the ORBJ method
    if (method == "ORBJ") then
      vkp1 = matmul(A_copy, x_new)
      bTvkp1 = dot_product(b, vkp1)
      vkp1Tvk = dot_product(vkp1, vk)
      vkp1Tvkp1 = dot_product(vkp1, vkp1)
      rel = (vkTvk - vkp1Tvk - bTvk + bTvkp1) / (vkTvk - vkp1Tvk + vkp1Tvkp1)
    end if

    ! Calculate new x for Jacobi methods
    if (jacobi) then
      x_new = (1.-rel)*x + rel*x_new
    end if

    ! Calculate error
    vk = matmul(A_copy, x_new)
    err = norm2(vk - b)
    dx = norm2(x - x_new)

    ! Update relaxation for ABSOR method
    if (method == "ABSOR") then

      ! Aim for a step size that's of the same order of magnitude as the error
      rel = rel - 0.1*(log10(dx/err))

      ! Nowhere should we go above 2
      if (rel > 2.) rel = 2.

      ! Make sure the relaxation factor stays positive
      if (rel < 0.) rel = 0.01
    end if

    ! Output progress
    if (verbose) write(unit,'(i6, a, ES10.3, a, ES10.3, a, ES10.3)') iteration, ',', dx, ',', err, ',', rel

    ! Prepare for next iteration
    x = x_new
    if (method == "ORBJ") then
      vkTvk = dot_product(vk, vk)
      bTvk = dot_product(b, vk)
    end if

  end do

  close(unit)

end subroutine iterative_solve


subroutine block_gauss_siedel(N, A, b, block_size, tol, x)
  ! Iteratively solves the [A]x=b system using block Gauss-Siedel
  ! THIS IS NOT GUARANTEED TO CONVERGE FOR MOST MATRICES.
  ! N is the size of the system
  ! A is the system matrix
  ! b is the RHS vector
  ! block_size is the desired size of each block
  ! tol is the convergence tolerance between iterations
  ! x is the solution

  implicit none

  integer,intent(in) :: N, block_size
  real,dimension(N,N),intent(in) :: A
  real,dimension(N),intent(in) :: b
  real,intent(in) :: tol
  real,dimension(:),allocatable,intent(out) :: x

  real :: err
  real,dimension(N) :: x_new
  real,dimension(block_size) :: bi
  real,dimension(:),allocatable :: xi
  integer :: i, N_blocks, r, i_start, i_end, N_last

  ! Give initial error estimate
  err = tol + 1.

  ! Initialize solution vector
  allocate(x(N), source=0.)

  ! Calculate number of blocks
  N_blocks = N/block_size
  r = modulo(N, block_size)
  if (r > 0) then
    N_blocks = N_blocks + 1
  end if

  ! Iterate
  do while(err >= tol)

    ! Invert blocks
    do i=1,N_blocks

      ! Last block
      if (i == N_blocks) then

        ! Determine start index of this block
        i_start = (i-1)*block_size + 1
        N_last = N-i_start+1

        ! Calculate new RHS vector
        bi(1:N_last) = b(i_start:N) - matmul(A(i_start:N,1:i_start-1), x_new(1:i_start-1))

        call lu_solve_replace(N_last, A(i_start:N,i_start:N), bi(1:N_last), xi)

        ! Store
        x_new(i_start:N) = xi

      ! Not last block
      else

        ! Determine start and end indices of this block
        i_start = (i-1)*block_size + 1
        i_end = i*block_size

        ! Calculate new RHS vector
        bi = b(i_start:i_end) - matmul(A(i_start:i_end,1:i_start-1), x_new(1:i_start-1))
        bi = bi - matmul(A(i_start:i_end,i_end+1:N), x(i_end+1:N))

        call lu_solve_replace(block_size, A(i_start:i_end,i_start:i_end), bi, xi)

        ! Store
        x_new(i_start:i_end) = xi

      end if
    end do

    ! Calculate error
    err = norm2(x-x_new)
    write(*,*) err

    ! Update
    x = x_new

  end do

end subroutine block_gauss_siedel

    
end module linalg_mod