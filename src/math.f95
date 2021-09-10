! Generic math subroutines
module math_mod

    implicit none
    REAL, parameter :: pi = 3.1415926535897932
    
contains


subroutine math_plane_normal(p1,p2,p3,ans)
    implicit none
    real,dimension(3) :: p1,p2,p3,a,b,ans

    a = p2 - p1
    b = p3 - p1
    call math_cross_product(a,b,ans)
end subroutine math_plane_normal


real function math_max(dim,vec)
    implicit none
    integer :: dim,i
    real :: vec(dim)
    math_max = vec(1)
    do i=1,dim
        if(vec(i) > math_max) math_max = vec(i)
    end do
end function math_max


real function math_length(dim,p1,p2)
    implicit none
    integer :: dim,i
    real,dimension(dim) :: p1,p2
    math_length = 0.0
    do i=1,dim
        math_length = math_length + (p2(i)-p1(i))**2
    end do
    math_length = sqrt(math_length)
end function math_length


subroutine math_reflect_point(A,B,C,D,P,ans)
    implicit none
    real :: A,B,C,D,P(3),ans(3),mult
    mult = 2.0*(A*P(1) + B*P(2) + C*P(3) + D)/(A**2 + B**2 + C**2)
    ans(1) = P(1) - mult*A
    ans(2) = P(2) - mult*B
    ans(3) = P(3) - mult*C
end subroutine math_reflect_point


real function math_mag(n,vec)
    implicit none
    integer :: n
    real :: vec(n)
    math_mag = sqrt(dot_product(vec,vec))
end function math_mag


subroutine math_cross_product(a,b,c)
    implicit none
    real :: a(3),b(3),c(3)
    c(1) = a(2)*b(3)-a(3)*b(2);
    c(2) = a(3)*b(1)-a(1)*b(3);
    c(3) = a(1)*b(2)-a(2)*b(1);
end subroutine math_cross_product


subroutine math_rot_x(vec,th)
    implicit none
    real :: vec(3),th,rm(3,3),ans(3)

    rm(1,1) = 1.0;
    rm(1,2) = 0.0;
    rm(1,3) = 0.0;
    rm(2,1) = 0.0;
    rm(2,2) = cos(th);
    rm(2,3) = -sin(th);
    rm(3,1) = 0.0;
    rm(3,2) = sin(th);
    rm(3,3) = cos(th);
    ans = matmul(rm,vec)
    vec = ans
end subroutine math_rot_x


subroutine math_rot_y(vec,th)
    implicit none
    real :: vec(3),th,rm(3,3),ans(3)

    rm(1,1) = cos(th);
    rm(1,2) = 0.0;
    rm(1,3) = sin(th);
    rm(2,1) = 0.0;
    rm(2,2) = 1.0;
    rm(2,3) = 0.0;
    rm(3,1) = -sin(th);
    rm(3,2) = 0.0;
    rm(3,3) = cos(th);
        ans = matmul(rm,vec)
        vec = ans
end subroutine math_rot_y


subroutine math_rot_z(vec,th)
    implicit none
    real :: vec(3),th,rm(3,3),ans(3)

    rm(1,1) = cos(th);
    rm(1,2) = -sin(th);
    rm(1,3) = 0.0;
    rm(2,1) = sin(th);
    rm(2,2) = cos(th);
    rm(2,3) = 0.0;
    rm(3,1) = 0.0;
    rm(3,2) = 0.0;
    rm(3,3) = 1.0;
        ans = matmul(rm,vec)
        vec = ans
end subroutine math_rot_z


subroutine math_matinv(n,a,ai)
      implicit none
!
! This sobroutine inverts a matrix "a" and returns the inverse in "ai"
! n  - Input by user, an integer specifying the size of the matrix to be inverted.
! a  - Input by user, an n by n real array containing the matrix to be inverted.
! ai - Returned by subroutine, an n by n real array containing the inverted matrix.
! d  - Work array, an n by 2n real array used by the subroutine.
! io - Work array, a 1-dimensional integer array of length n used by the subroutine.
!
      integer :: n,i,j,k,m,itmp
      real :: a(n,n),ai(n,n),tmp,r !,d(n,2*n)
!      integer :: io(n)
      real,allocatable,dimension(:,:) :: d
      integer,allocatable,dimension(:) :: io

      allocate(d(n,2*n))
      allocate(io(n))

      d(:,:) = 0.0
      io(:) = 0

!      write(6,*)'Inverting vortex panel matrix.  Please wait.'
!      itime1=mclock()
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
!     ********************
      do i=1,n
         do j=1,n
            ai(i,j)=d(io(i),n+j)
         end do
      end do
!      itime2=mclock()
!      rtime=real(itime2-itime1)/1000.
!      write(6,'(a,f7.2,a)')' Matrix inversion time =',rtime,' sec'
      deallocate(d)
      deallocate(io)
      return
      end subroutine math_matinv


SUBROUTINE math_AXB_LUD(n,A,B,X)
!Solves a general [A]*X=B on an nxn matrix
    IMPLICIT NONE
    INTEGER::n,D,info
    real,DIMENSION(n)::B,X
    real,DIMENSION(n,n)::A

    INTEGER,allocatable,DIMENSION(:) :: INDX

    allocate(INDX(n))

    CALL math_LUDCMP(A,n,INDX,D,info)
    if(info.eq.0) then
        CALL math_LUBKSB(A,n,INDX,B)
    end if
    if(info.eq.1) then
        write(*,*) ' The system matrix is singular. No solution.'
    end if
    X = B
    deallocate(INDX)
    RETURN
END SUBROUTINE math_AXB_LUD

!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     *
!*******************************************************
!MODULE LU

!CONTAINS

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
Subroutine math_LUDCMP(A,N,INDX,D,CODE)
  implicit none
  integer, PARAMETER :: NMAX=100
  REAL, parameter :: TINY=1.5D-16
  real  AMAX,DUM, SUM, A(N,N)!,VV(N)
  real,allocatable,dimension(:) :: VV
  INTEGER N, CODE, D, INDX(N)
  integer :: I,J,K,IMAX

  allocate(VV(N))

  D=1; CODE=0; IMAX = 0

  DO I=1,N
    AMAX=0.0
    DO J=1,N
      IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
    END DO ! j loop
    IF(AMAX.LT.TINY) THEN
      CODE = 1
      RETURN
    END IF
    VV(I) = 1.0 / AMAX
  END DO ! i loop

  DO J=1,N
    DO I=1,J-1
      SUM = A(I,J)
      DO K=1,I-1
        SUM = SUM - A(I,K)*A(K,J)
      END DO ! k loop
      A(I,J) = SUM
    END DO ! i loop
    AMAX = 0.0
    DO I=J,N
      SUM = A(I,J)
      DO K=1,J-1
        SUM = SUM - A(I,K)*A(K,J)
      END DO ! k loop
      A(I,J) = SUM
      DUM = VV(I)*ABS(SUM)
      IF(DUM.GE.AMAX) THEN
        IMAX = I
        AMAX = DUM
      END IF
    END DO ! i loop

    IF(J.NE.IMAX) THEN
      DO K=1,N
        DUM = A(IMAX,K)
        A(IMAX,K) = A(J,K)
        A(J,K) = DUM
      END DO ! k loop
      D = -D
      VV(IMAX) = VV(J)
    END IF

    INDX(J) = IMAX
    IF(ABS(A(J,J)) < TINY) A(J,J) = TINY

    IF(J.NE.N) THEN
      DUM = 1.0 / A(J,J)
      DO I=J+1,N
        A(I,J) = A(I,J)*DUM
      END DO ! i loop
    END IF
  END DO ! j loop

  deallocate(VV)
  RETURN
END subroutine math_LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine math_LUBKSB(A,N,INDX,B)
 implicit none
 integer :: N
 real  SUM, A(N,N),B(N)
 INTEGER INDX(N)
 integer :: II,I,J,LL

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.0) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine math_LUBKSB


!C
!C--------------------------------------------------------------------C
!C The following subroutine computes the LU decomposition for a       C
!C diagonally dominant matrix (no pivoting is done).                  C
!C   Inputs:  n = number of equations/unknowns                        C
!C            a = nxn coefficient matrix                              C
!C   Outputs: a = nxn matrix containing the LU matrices               C
!C                                                                    C
!C Deryl Snyder, 10-16-98                                             C
!C--------------------------------------------------------------------C
      subroutine math_snyder_ludcmp(a,n)
      implicit none
      integer :: n,i,j,k
      real a(n,n),z

      do k=1,n-1
        do i=k+1,n
          z=a(i,k)/a(k,k)                     !compute gauss factor
          a(i,k)=z                            !store gauss factor in matrix
          do j=k+1,n
            a(i,j)=a(i,j)-z*a(k,j)            !apply row operation
          enddo
        enddo
      enddo
      return
      end subroutine math_snyder_ludcmp
!C--------------------------------------------------------------------C
!C The following subroutine solves for the unknowns (x) given the LU  C
!C matrix and the right hand side.                                    C
!C   Inputs:  n = number of equations/unknowns                        C
!C            a = nxn matrix containing the L and U values            C
!C            b = n vector containing right hand side values          C
!C   Outputs: x = n vector containing solution                        C
!C                                                                    C
!C Deryl Snyder, 10-16-98                                             C
!C--------------------------------------------------------------------C
      subroutine math_snyder_lusolv(a,b,x,n)
      implicit none
      integer :: n,i,j,k
      real a(n,n),b(n),x(n)
      do i=1,n
         x(i)=b(i)
      end do
      do k=1,n-1                                 !do forward substitution
        do i=k+1,n
          x(i)=x(i)-a(i,k)*x(k)
        enddo
      enddo
      do i=n,1,-1                                !do back subsitution
        do j=i+1,n
          x(i)=x(i)-a(i,j)*x(j)
        enddo
        x(i)=x(i)/a(i,i)
      enddo
      return
      end subroutine math_snyder_lusolv


!C--------------------------------------------------------------------C
!C The following subroutine fits a parabola through three specified   C
!C points and returns the coefficients a, b, c defining this parabola C
!C according to the equation y = a * x**2 + b * x + c                 C
!C            pts = list of three (x, y) points                       C
!C   Outputs: a, b, c = quadratic coefficients                        C
!C                                                                    C
!C--------------------------------------------------------------------C
      subroutine quadratic_fit(pts, a, b, c)
      implicit none
      real, dimension(3, 2), intent(in) :: pts
      real, intent(out) :: a, b, c

      integer :: i
      real, dimension(3, 3) :: m, m_inv
      real, dimension(3) :: v, coeff

      do i = 1, 3
        m(i, 1) = pts(i, 1)**2
        m(i, 2) = pts(i, 1)
        m(i, 3) = 1.0
      end do
      call math_matinv(3, m, m_inv)

      v(:) = pts(:, 2)
      coeff = matmul(m_inv, v)

      a = coeff(1)
      b = coeff(2)
      c = coeff(3)

      end subroutine quadratic_fit


end module math_mod
