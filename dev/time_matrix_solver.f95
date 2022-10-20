program time_matrix_solver

    use linalg_mod

    implicit none

    character(100) :: solver, A_matrix_file, b_vector_file, dummy_read
    character(len=:),allocatable :: output_file
    integer :: N, unit, i, stat
    real,dimension(:,:),allocatable :: A
    real,dimension(:),allocatable :: b, x
    real :: count_rate, rel
    integer :: start_count, end_count

    ! Initialize
    rel = 0.5
    !allocate(output_file(4))
    output_file = 'none'

    ! Get arguments
    call getarg(1, solver)
    call getarg(2, A_matrix_file)
    call getarg(3, b_vector_file)

    ! Figure out the size of the system
    N = 0
    open(newunit=unit, file=b_vector_file)
        do
            ! Read line
            read(unit,*,iostat=stat) dummy_read

            ! Check status
            if (stat /= 0) then
                ! At end of file
                exit
            else
                N = N + 1
            end if

        end do
    close(unit)

    ! Allocate space for A and b
    allocate(A(N,N))
    allocate(b(N))
    
    ! Read in A
    open(newunit=unit, file=A_matrix_file)
    do i=1,N
        read(unit,*) A(i,:)
    end do
    close(unit)

    ! Read in b
    open(newunit=unit, file=b_vector_file)
    do i=1,N
        read(unit,*) b(i)
    end do
    close(unit)

    ! Solve
    select case (solver)

    ! LU decomposition
    case ('LU')
        call system_clock(start_count, count_rate)
        call lu_solve(N, A, b, x)
        call system_clock(end_count)

    ! QR via Givens rotations for upper-pentagonal
    case ('QRUP')
        call system_clock(start_count, count_rate)
        call QR_row_givens_solve_UP(N, A, b, x)
        call system_clock(end_count)

    ! QR via fast Givens rotations for upper-pentagonal
    case ('FQRUP')
        call system_clock(start_count, count_rate)
        call QR_fast_givens_solve_upper_pentagonal(N, A, b, x)
        call system_clock(end_count)

    ! GMRES
    case ('GMRES')
        call system_clock(start_count, count_rate)
        call GMRES(N, A, b, 1.e-12, 1000, output_file, i, x)
        call system_clock(end_count)

    ! Block successive over-relaxation
    case ('BSOR')
        call system_clock(start_count, count_rate)
        call block_sor_solve(N, A, b, 400, 1.e-12, rel, 1000, output_file, i, x)
        call system_clock(end_count)
    
    ! Block Jacobi
    case ('BJAC')
        call system_clock(start_count, count_rate)
        call block_jacobi_solve(N, A, b, 400, 1.e-12, rel, 1000, output_file, i, x)
        call system_clock(end_count)

    ! Improper specification
    case default
        write(*,*) "!!! ", solver, " is not a valid option. Quitting..."
        stop
    end select

    ! Report time
    write(*,*) "Solution time: ", end_time-start_time, " s"
    
end program time_matrix_solver