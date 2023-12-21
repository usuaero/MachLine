program sparse_test
    ! tests sparse vector and sparse matrix operations
    
    use adjoint_mod
    
    implicit none

    type(sparse_vector):: sparse_a, sparse_b 
    type(sparse_matrix) :: sparse_matrix_a, sparse_matrix_b
    real,dimension(:),allocatable :: vector_a, expanded_a, residual
    integer :: i, full_size_a, add_index, add_shift_index
    real :: add_value

    write(*,*) ""
    write(*,*) ""
    write(*,*) "-----------------------------SPARSE TEST------------------------------"
    write(*,*) ""
    write(*,*) ""


    !!!!!!!!!!!! INITIALIZE VECTOR !!!!!!!!!!!!!!
    write(*,*) "--------INITIALIZE A REGULAR VECTOR---------"
    vector_a = (/ 1.0, 0.0, 11.0, 0.0, 87.0, -1.0, 0.0, 0.002121, 4.0, 0.0, 3.14 /)
    ! allocate(vector_a(100), source=0.0)
    ! do i= 1,9
    !     vector_a(i*10) = real(i) 
    ! end do

    ! write results
    write(*,*) "full vector_a"
    write(*, '(f14.10, 2x)') vector_a
    write(*,*) ""
    write(*,*) ""


    !!!!!!!!!! TEST COMPRESS VECTOR !!!!!!!!!!
    write(*,*) "--------TEST COMPRESS VECTOR---------"
    ! compress vector_a into sparse_a
    call sparse_a%compress(vector_a)
    
    ! write results
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) " "
    write(*,*) ""


    !!!!!!!!!! TEST EXPAND VECTOR  !!!!!!!!!!
    write(*,*) "--------TEST EXPAND VECTOR---------"
    ! expand sparse_a
    expanded_a = sparse_a%expand()
    residual = vector_a - expanded_a
    
    ! write results
    write(*,*) "expanded_a          residual"
    do i=1,sparse_a%full_size
        write(*, '(2(f14.10, 2x))') expanded_a(i), residual(i)
    end do
    write(*,*) ""
    write(*,*) ""

    !!!!!!!!!! TEST INCREASE_SIZE VECTOR  !!!!!!!!!!
    write(*,*) "--------TEST INCREASE SIZE VECTOR---------"
    write(*,*) "sparse_size = ", sparse_a%sparse_size
    write(*,*) ""

    ! write sparse_a
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    
    write(*,*) ""
    write(*,*) "increase size....."
    write(*,*) ""
    
    call sparse_a%increase_size()
    
    write(*,*) "sparse_size = ", sparse_a%sparse_size
    write(*,*) ""

    ! write sparse_a again
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""
    write(*,*) ""

    ! return sparse_a to original size
    deallocate(sparse_a%elements)
    sparse_a%sparse_size = sparse_a%sparse_size -1
    call sparse_a%compress(vector_a)

    !!!!!!!!!! TEST ADD_ELEMENT VECTOR  !!!!!!!!!!
    write(*,*) "--------TEST ADD ELEMENT VECTOR---------"
    write(*,*) ""

    ! write sparse_a increased size
    write(*,*) "original sparse_a;"
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do

    write(*,*) ""
    write(*,*) "add nonzero value of 5.1 at index", 4
    write(*,*) ""
    write(*,*) "adding element...."
    write(*,*) ""
    add_value = 5.1
    add_index = 4
    add_shift_index = 3
    call sparse_a%add_element(add_value, add_index, add_shift_index)
    ! write sparse_a 
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    



end program sparse_test