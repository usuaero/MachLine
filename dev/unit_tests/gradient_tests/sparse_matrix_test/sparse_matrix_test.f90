program sparse_matrix_test
    ! tests sparse vector and sparse matrix operations
    
    use adjoint_mod
    
    implicit none

    type(sparse_vector) :: sparse_a, sparse_b, sparse_c 
    type(sparse_matrix) :: sparse_matrix_a, sparse_matrix_b
    real,dimension(:),allocatable :: vector_a, vector_b, vector_c
    real,dimension(3) :: vec, values, add_values
    real,dimension(:,:),allocatable :: scaled_vecs, full_matrix_a, full_matrix_b, expanded_matrix_b, residual
    integer :: i, passed_tests, total_tests, nonzeros, old_size, add_full_index, add_shift_index
    logical :: test_failed

    write(*,*) ""
    write(*,*) ""
    write(*,*) "-----------------------------SPARSE MATRIX TEST------------------------------"
    write(*,*) ""
    write(*,*) ""
    test_failed = .true. ! assume test failed, if the test condition is met, test passed
    ! NOTE: on the sparse vector test, I assume the test passes, if it fails a test condition, test fails
    passed_tests = 0
    total_tests = 0
    


!!!!!!!!!!!! INITIALIZE SPARSE VECTORS !!!!!!!!!!!!!!
    write(*,*) "--------INITIALIZE SPARSE VECTORS---------"
    write(*,*) ""

    vector_a = (/ 0.0, 1.0, 11.0, 0.0, 87.0, -1.0, 0.0, 0.002121, 0.0, 0.0, 3.14 /)
    vector_b = (/ 0.0, 6.0,  0.0, 0.0, -2.0, -1.0, 0.0, 0.0,      0.0, 0.0, 3.14 /)
    vector_c = (/ 0.0, 1.0, -7.0, 0.0,  0.0, 17.0, 0.0, 1.5,      0.0, 0.0, 3.14 /)


    ! write vectors
    write(*,*) "   vector_a      vector_b      vector_c"
    do i=1,11
        write(*, '(3(f12.5, 2x))') vector_a(i), vector_b(i), vector_c(i)
    end do
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!! TEST INIT FROM SPARSE VECTORS !!!!!!!!!!
    write(*,*) "-------------TEST INIT FROM SPARSE VECTORS--------------"
    write(*,*) ""
    write(*,*) "sparse vectors:"
    write(*,*) ""

     ! initialize and compress vector_a into sparse_a
    call sparse_a%init(vector_a)
    call sparse_b%init(vector_b)
    call sparse_c%init(vector_c)
    
    ! write sparse a
    write(*,*) "sparse_a                 sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f12.5, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""

    ! write sparse b
    write(*,*) "sparse_b                 sparse_index      full_index"
    do i=1,sparse_b%sparse_size
        write(*,'(f12.5, 12x, I5, 12x, I5)') sparse_b%elements(i)%value, i, sparse_b%elements(i)%full_index
    end do
    write(*,*) ""

    ! write sparse c
    write(*,*) "sparse_c                 sparse_index      full_index"
    do i=1,sparse_c%sparse_size
        write(*,'(f12.5, 12x, I5, 12x, I5)') sparse_c%elements(i)%value, i, sparse_c%elements(i)%full_index
    end do
    write(*,*) ""

    ! make a full matrix
    allocate(full_matrix_a(3,11))
    do i=1,11 
        full_matrix_a(:,i) = (/vector_a(i), vector_b(i), vector_c(i) /)
    end do 

    ! write full matrix
    write(*,*) ""
    write(*,*) "             full matrix a   "
    write(*,*) ""
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! init from sparse vectors
    call sparse_matrix_a%init_from_sparse_vectors(sparse_a, sparse_b, sparse_c)
    write(*,*) "sparse matrix a size: ", size(sparse_matrix_a%columns)

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! check if test failed
    if (any((full_matrix_a(:,2) /= sparse_matrix_a%columns(1)%vector_values)) &
    .or. any((full_matrix_a(:,3) /= sparse_matrix_a%columns(2)%vector_values)) &
    .or. any((full_matrix_a(:,5) /= sparse_matrix_a%columns(3)%vector_values)) &
    .or. any((full_matrix_a(:,6) /= sparse_matrix_a%columns(4)%vector_values)) &
    .or. any((full_matrix_a(:,8) /= sparse_matrix_a%columns(5)%vector_values)) &
    .or. any((full_matrix_a(:,11) /= sparse_matrix_a%columns(6)%vector_values))) then    
    
        test_failed = .true.
    else 
        test_failed = .false.
        
    end if


    if (test_failed) then
        write(*,*) "init from sparse vectors FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "init from sparse vectors PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) "" 
     


!!!!!!!!!! TEST INIT FROM FULL MATRIX !!!!!!!!!!
    write(*,*) "-------------TEST INIT FROM FULL MATRIX--------------"
    write(*,*) ""
    
    !reset values
    write(*,*) "deallocate sparse matrix a, so it can be initialized from a full matrix"
    deallocate(sparse_matrix_a%columns)
    write(*,*) "deallocation successful"

    ! write full matrix
    write(*,*) ""
    write(*,*) "             full matrix a   "
    write(*,*) ""
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! init from full matrix
    call sparse_matrix_a%init_from_full_matrix(full_matrix_a)
    write(*,*) "sparse matrix a size: ", sparse_matrix_a%sparse_num_cols

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! check if test failed
    if (any((full_matrix_a(:,2) /= sparse_matrix_a%columns(1)%vector_values)) &
    .or. any((full_matrix_a(:,3) /= sparse_matrix_a%columns(2)%vector_values)) &
    .or. any((full_matrix_a(:,5) /= sparse_matrix_a%columns(3)%vector_values)) &
    .or. any((full_matrix_a(:,6) /= sparse_matrix_a%columns(4)%vector_values)) &
    .or. any((full_matrix_a(:,8) /= sparse_matrix_a%columns(5)%vector_values)) &
    .or. any((full_matrix_a(:,11) /= sparse_matrix_a%columns(6)%vector_values))) then    
    
        test_failed = .true.
    else 
        test_failed = .false.
        
    end if


    if (test_failed) then
        write(*,*) "init from sparse vectors FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "init from sparse vectors PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) "" 
    


!!!!!!!!!! TEST INCREASE_SIZE (SPARSE MATRIX)  !!!!!!!!!!
    write(*,*) "-------------TEST INCREASE_SIZE (SPARSE MATRIX)--------------"
    write(*,*) ""
    old_size = sparse_matrix_a%sparse_num_cols
    write(*,*) "sparse num cols = ", old_size
    write(*,*) ""

    ! write sparse matrix a again
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I10)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    write(*,*) ""
    write(*,*) "increase size....."
    write(*,*) ""
    
    call sparse_matrix_a%increase_size()
    
    write(*,*) "sparse num cols = ", sparse_matrix_a%sparse_num_cols
    write(*,*) ""

    ! write sparse matrix a again
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! check 
    if (sparse_matrix_a%sparse_num_cols /= (old_size + 1))then
        test_failed = .true.
    else
        test_failed = .false.
    end if

    write(*,*) ""
    write(*,*) "return to original size"
    write(*,*) ""

    ! reset values
    deallocate(sparse_matrix_a%columns)

    !write(*,*) "deallocation successful" 
    call sparse_matrix_a%init_from_sparse_vectors(sparse_a, sparse_b, sparse_c)

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! result of test
    if (test_failed) then
        write(*,*) "increase size (sparse matrix) test FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "increase size (sparse matrix) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""


!!!!!!!!!! TEST ADD ELEMENT (SPARSE MATRIX)   !!!!!!!!!!
    write(*,*) "-------------TEST ADD ELEMENT (SPARSE MATRIX) --------------"
    write(*,*) ""

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    write(*,*) ""
    write(*,*) "add a sparse element with values (1.0, 2.0, 3.0) at full index", 4
    write(*,*) ""
    write(*,*) "adding element...."
    write(*,*) ""
    add_values = (/1.0,2.0,3.0/)
    add_full_index = 4
    add_shift_index = 3
    call sparse_matrix_a%add_element(add_values, add_full_index, add_shift_index)
    
   ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! check if test failed
    full_matrix_a(:,4) = (/1.0, 2.0, 3.0/)
    do i=1,sparse_matrix_a%full_num_cols
        if (any(abs(sparse_matrix_a%get_values(i) - full_matrix_a(:,i)) > 1.0e-12))  then
            test_failed = .true.
            exit
        else
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        write(*,*) "add element (sparse matrix) test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "add element (sparse matrix) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""

    

!!!!!!!!!! TEST GET_VALUES  !!!!!!!!!!
    write(*,*) "-------------TEST GET_VALUES--------------"
    write(*,*) ""

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    do i=1,11    
        write(*,'(a,I3)') "Get element values of full index",i 
        values = sparse_matrix_a%get_values(i)
        if (any(abs(sparse_matrix_a%get_values(i)-full_matrix_a(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if 
        write(*,'(a, 3(f12.5))') "values = ", values
        write(*,*) ""
    end do
    
    ! check if test failed
    if (test_failed) then
        write(*,*) "get element values test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "get element values test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""



!!!!!!!!!! TEST SET_VALUES (SPARSE MATRIX)  !!!!!!!!!!
    write(*,*) "-------------TEST SET_VALUES (SPARSE MATRIX)--------------"
    write(*,*) ""

    ! make b matrices
    allocate(full_matrix_b(3,11))
    full_matrix_b(:,:) = full_matrix_a(:,:)
    call sparse_matrix_b%init_from_sparse_vectors(sparse_a, sparse_b, sparse_c)

    ! set values
    write(*,*) "Set values of full index 2 to (-9.9, -5.2, 1.0)"
    vec = (/ -9.9, -5.2, 1.0/)
    full_matrix_b(:,2) = vec 
    call sparse_matrix_b%set_values(vec, 2)
    
    write(*,*) "Set values of full index 3 to (0.0, 11.11, 2.0)" 
    vec = (/0.0, 11.11, 2.0/)
    full_matrix_b(:,3) = vec 
    call sparse_matrix_b%set_values(vec, 3)

    write(*,*) "Set values of full index 4 to (-2.22, 3.14, 6.3)" 
    vec = (/-2.22, 3.14, 6.3/)
    full_matrix_b(:,4) = vec 
    call sparse_matrix_b%set_values(vec, 4)

    write(*,*) "Set values of full index 10 to (1.0, 1.0, 1.0)" 
    vec = (/1.0, 1.0, 1.0/)
    full_matrix_b(:,10) = vec
    call sparse_matrix_b%set_values(vec, 10)

    write(*,*) "Set values of full index 6,7,8,9,11 to 0.0" 
    vec = (/ 0.0, 0.0, 0.0 /)
    full_matrix_b(:,6) = vec
    call sparse_matrix_b%set_values(vec, 6)
    full_matrix_b(:,7) = vec
    call sparse_matrix_b%set_values(vec, 7)
    full_matrix_b(:,8) = vec
    call sparse_matrix_b%set_values(vec, 8)
    full_matrix_b(:,9) = vec
    call sparse_matrix_b%set_values(vec, 9)
    full_matrix_b(:,11) = vec
    call sparse_matrix_b%set_values(vec, 11)
    write(*,*) ""
    
    ! write full matrix
    write(*,*) "resulting full matrix b:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_b(:,i)
    end do
    write(*,*) ""
    
    write(*,*) ""
    write(*,*) "resulting sparse matrix:"
    
    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(1), &
        sparse_matrix_b%columns(i)%vector_values(2),sparse_matrix_b%columns(i)%vector_values(3), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_b%full_num_cols
        if (any(abs(sparse_matrix_b%get_values(i) - full_matrix_b(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        write(*,*) "set values (sparse matrix) test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "set values (sparse matrix) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""
    

    !!!!!!!!!! TEST GET_VALUES  !!!!!!!!!!
    write(*,*) "-------------TEST GET_VALUES again, matrix b--------------"
    write(*,*) ""

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(1), &
        sparse_matrix_b%columns(i)%vector_values(2),sparse_matrix_b%columns(i)%vector_values(3), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    do i=1,11    
        write(*,'(a,I3)') "Get element values of full index",i 
        values = sparse_matrix_b%get_values(i)
        if (any(abs(sparse_matrix_b%get_values(i)-full_matrix_b(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if 
        write(*,*) "values = ", values
        write(*,*) ""
    end do
    
    ! check if test failed
    if (test_failed) then
        write(*,*) "get element values test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "get element values test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""

    
    
!!!!!!!!!! TEST COMPRESS (SPARSE MATRIX)  !!!!!!!!!!
    write(*,*) "-------------TEST COMPRESS (SPARSE MATRIX)--------------"
    write(*,*) ""
    write(*,*) "compressing sparse matrix:"
    write(*,*) ""

    call sparse_matrix_b%compress()

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(1), &
        sparse_matrix_b%columns(i)%vector_values(2),sparse_matrix_b%columns(i)%vector_values(3), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""


    ! check if test failed
    do i=1,11
        if (any(abs(sparse_matrix_b%get_values(i) - full_matrix_b(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else
            test_failed = .false.
        end if
    end do


    if (sparse_matrix_b%sparse_num_cols /= 5) then
        test_failed = .true.
    else 
        test_failed = .false.
    end if


    if (test_failed) then
        write(*,*) "compress (sparse matrix) test FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "compress (sparse matrix) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""
    

 !!!!!!!!!! TEST GET_VALUES  !!!!!!!!!!
    write(*,*) "-------------TEST GET_VALUES after compress, matrix b--------------"
    write(*,*) ""

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(1), &
        sparse_matrix_b%columns(i)%vector_values(2),sparse_matrix_b%columns(i)%vector_values(3), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    do i=1,11    
        write(*,'(a,I3)') "Get element values of full index",i 
        values = sparse_matrix_b%get_values(i)
        if (any(abs(sparse_matrix_b%get_values(i)-full_matrix_b(:,i)) > 1.0e-12)) then
            test_failed = .true.
            write(*,*) "failed b value= ", sparse_matrix_b%get_values(i)
            write(*,*) "failed b value= ", sparse_matrix_b%get_values(i)
            write(*,*) "failed b value= ", sparse_matrix_b%get_values(i)
            write(*,*) "failed b value= ", sparse_matrix_b%get_values(i)
            write(*,*) "failed b value= ", sparse_matrix_b%get_values(i)
            exit
        else 
            test_failed = .false.
        end if 
        write(*,*) "values = ", values
        write(*,*) ""
    end do
    
    
    ! check if test failed
    if (test_failed) then
        write(*,*) "get element values test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "get element values test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""
   
    
    !!!!!!!!!! TEST EXPAND SPARSE MATRIX  !!!!!!!!!!
    write(*,*) "-------------TEST EXPAND SPARSE MATRIX--------------"
    write(*,*) ""
    
    ! expand sparse_matrix_b
    allocate(expanded_matrix_b(3,sparse_matrix_b%full_num_cols))
    allocate(residual(3,sparse_matrix_b%full_num_cols))
    expanded_matrix_b = sparse_matrix_b%expand()
    residual = full_matrix_b - expanded_matrix_b
    
    ! write results
    write(*,*) "           expanded_matrix_b                        residual"
    do i=1,sparse_matrix_b%full_num_cols
        write(*, '(3(f10.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_b(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_b%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        write(*,*) "expand sparse matrix FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "expand sparse matrix PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""

!!!!!!!!!! TEST GET_VALUES  !!!!!!!!!!
    write(*,*) "-------------TEST GET_VALUES after expand, matrix b--------------"
    write(*,*) ""

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(1), &
        sparse_matrix_b%columns(i)%vector_values(2),sparse_matrix_b%columns(i)%vector_values(3), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    do i=1,11    
        write(*,'(a,I3)') "Get element values of full index",i 
        values = sparse_matrix_b%get_values(i)
        if (any(abs(sparse_matrix_b%get_values(i)-full_matrix_b(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if 
        write(*,*) "values = ", values
        write(*,*) ""
    end do
    
    ! check if test failed
    if (test_failed) then
        write(*,*) "get element values test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "get element values test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) ""





!!!!!!!!!! TEST SPARSE ADD (MATRIX)  !!!!!!!!!!
    write(*,*) "-------------TEST SPARSE ADD (MATRIX)--------------"
    write(*,*) ""
    
    ! update full_matrix_a by expanding the most recent sparse matrix a
    write(*,*) "update full matrix a to reflect changes made to sparse matrix a"
    write(*,*) ""
    full_matrix_a = sparse_matrix_a%expand()

    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! write full matrix a
    write(*,*) " full matrix b:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_b(:,i)
    end do
    write(*,*) ""

    ! add full matrix b to full matrix c
    write(*,*) "add full matrix b to full matrix a"
    write(*,*) ""
    do i=1,11 
        full_matrix_a(:,i) = full_matrix_a(:,i) + full_matrix_b(:,i)
    end do 

    ! display new full matrix a
    write(*,*) " new full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! add full matrix b to full matrix c
    
    write(*,*) ""
    write(*,*) "display sparse matrices before adding"
    write(*,*) ""

    ! write sparse matrix a
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! write sparse matrix b
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(1), &
        sparse_matrix_b%columns(i)%vector_values(2),sparse_matrix_b%columns(i)%vector_values(3), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! add sparse matrix b to sparse matrix a
    write(*,*) ""
    write(*,*) "adding sparse_b to sparse_a...."
    write(*,*) ""
    write(*,*) "sparse b column 11 get before call:", sparse_matrix_b%get_values(11)
    call sparse_matrix_a%sparse_add(sparse_matrix_b)

    write(*,*) "sparse b column 11 get after after call:", sparse_matrix_b%get_values(11)

    write(*,*) "     new sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(1), &
        sparse_matrix_a%columns(i)%vector_values(2),sparse_matrix_a%columns(i)%vector_values(3), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    !write(*,*) "sparse b column 6 values:", sparse_matrix_b%columns(6)%vector_values
    write(*,*) "sparse b column 6 get :", sparse_matrix_b%get_values(6)
    write(*,*) "sparse b column 11 get :", sparse_matrix_b%get_values(11)

    ! ! expand sparse_a
    ! expanded_ab = sparse_ab%expand()
    ! residual = vector_ab - expanded_ab
    
    ! ! write results
    ! write(*,*) "residual = vector_ab - expanded_ab"
    ! write(*,*) ""
    ! write(*,*) "expanded_ab          residual"
    ! do i=1,sparse_ab%full_size
    !     write(*, '(2(f14.10, 2x))') expanded_ab(i), residual(i)
    ! end do
    ! write(*,*) ""

    ! ! check if test failed
    ! do i=1,sparse_a%full_size
    !     if (abs(residual(i)) > 1.0e-12) then
    !         test_failed = .true.
    !     end if
    ! end do
    ! if (test_failed) then
    !     write(*,*) "sparse add (matrix) FAILED"
    !     total_tests = total_tests + 1
    ! else
    !     write(*,*) "sparse add (matrix) PASSED"
    !     passed_tests = passed_tests + 1
    !     total_tests = total_tests + 1
    ! end if
    ! test_failed = .false.
    ! write(*,*) ""
    ! write(*,*) ""

! !!!!!!!!!! TEST SPARSE SUBTRACT (MATRIX, a - b and b - a)  !!!!!!!!!!
!     write(*,*) "-------------TEST SPARSE SUBTRACT (MATRIX, a - b and b - a)--------------"
!     write(*,*) ""
!     write(*,*) "vectors to be subtracted:"
!     vector_a = (/ 1.0, 0.0, 11.0, 0.0, 87.0, 0.0, 0.0, 0.0, 4.0, 0.0, 3.14 /)
!     vector_b = (/ 1.0, 0.0, 1.0, 0.0, 1.0, 7.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
!     vector_a_minus_b = (/ 0.0, 0.0, 10.0, 0.0, 86.0, -7.0, 0.0, 0.0, 3.0, 0.0, 2.14 /)
!     vector_b_minus_a = (/ 0.0, 0.0, -10.0, 0.0, -86.0, 7.0, 0.0, 0.0, -3.0, 0.0, -2.14 /)
    
!     ! write  full vectors
!     write(*,*) ""
!     write(*,*) " vector_a     -     vector_b   =   vector_a_minus_b"
!     do i=1,11
!         write(*, '(3(f14.10, 2x))') vector_a(i), vector_b(i), vector_a_minus_b(i)
!     end do
!     write(*,*) ""
    
!     ! check if test failed
!     do i=1,sparse_a%full_size
!         if (abs(residual(i)) > 1.0e-12) then
!             test_failed = .true.
!         end if
!     end do
!     if (test_failed) then
!         write(*,*) "sparse subtract (matrix) a - b FAILED"
!         total_tests = total_tests + 1
!     else
!         write(*,*) "sparse subtract (matrix) a - b PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) ""
!     write(*,*) ""
    
    
!     ! subtract and display 
!     write(*,*) "subtracting sparse_a minus sparse_b...."
!     write(*,*) ""
!     sparse_a_minus_b = sparse_a%sparse_subtract(sparse_b)
!     write(*,*) "sparse_a_minus_b value    sparse_index      full_index "
!     do i=1,sparse_a_minus_b%sparse_size
!         write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a_minus_b%elements(i)%value, i, sparse_a_minus_b%elements(i)%full_index
!     end do
!     write(*,*) ""

!     ! expand 
!     expanded_a_minus_b = sparse_a_minus_b%expand()
!     residual = vector_a_minus_b - expanded_a_minus_b
    
!     ! write results
!     write(*,*) "residual = vector_a_minus_b - expanded_a_minus_b"
!     write(*,*) ""
!     write(*,*) "expanded_a_minus_b          residual"
!     do i=1,sparse_a_minus_b%full_size
!         write(*, '(2(f14.10, 2x))') expanded_a_minus_b(i), residual(i)
!     end do
!     write(*,*) ""
    
!     ! write  full vectors
!     write(*,*) "-------now do vector b minus vector a-------"
!     write(*,*) ""
!     write(*,*) " vector_b     -     vector_a   =   vector_b_minus_a"
!     do i=1,11
!         write(*, '(3(f14.10, 2x))') vector_b(i), vector_a(i), vector_b_minus_a(i)
!     end do
!     write(*,*) "" 
    
    
!     ! subtract and display 
!     write(*,*) "subtracting sparse_b minus sparse_a...."
!     write(*,*) ""
!     sparse_b_minus_a = sparse_b%sparse_subtract(sparse_a)
!     write(*,*) "sparse_a_minus_b value    sparse_index      full_index "
!     do i=1,sparse_b_minus_a%sparse_size
!         write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_b_minus_a%elements(i)%value, i, sparse_b_minus_a%elements(i)%full_index
!     end do
!     write(*,*) ""

!     ! expand 
!     expanded_b_minus_a = sparse_b_minus_a%expand()
!     residual = vector_b_minus_a - expanded_b_minus_a
    
!     ! write results
!     write(*,*) "residual = vector_b_minus_a - expanded_b_minus_a"
!     write(*,*) ""
!     write(*,*) "expanded_b_minus_a          residual"
!     do i=1,sparse_b_minus_a%full_size
!         write(*, '(2(f14.10, 2x))') expanded_b_minus_a(i), residual(i)
!     end do
!     write(*,*) ""
    
!     ! check if test failed
!     do i=1,sparse_a%full_size
!         if (abs(residual(i)) > 1.0e-12) then
!             test_failed = .true.
!         end if
!     end do
!     if (test_failed) then
!         write(*,*) "sparse subtract (matrix) b - a FAILED"
!         total_tests = total_tests + 1
!     else
!         write(*,*) "sparse subtract (matrix) b - a PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) ""
!     write(*,*) ""
    


! !!!!!!!!!! TEST BROADCAST VECTOR CROSS ELEMENT   !!!!!!!!!!

!     write(*,*) "-------------TEST BROADCAST VECTOR CROSS ELEMENT--------------"
!     write(*,*) ""
!     vec = (/1.0,-2.0,3.0/)
!     write(*,'(A, 3(f10.5, 2x))') "new vector to be scaled: ", vec
!     write(*,*) ""
!     write(*,*) "unscaled vectors"
!     do i=1,sparse_a%full_size
!         write(*, '(3(f14.5, 2x))') vec
!     end do
!     write(*,*) ""
    
!     ! scale and write results
!     sparse_matrix_a = sparse_a%broadcast_element_times_vector(vec)
    
!     write(*,*) "scale vec by each sparse element...."
!     write(*,*) ""
!     write(*,*) "resulting sparse_matrix:"
!     write(*,*) ""
!     write(*,*) "                   scaled vector                    sparse_index       &
!               full_index "
!     write(*,*) ""
!     do i=1,sparse_matrix_a%sparse_num_cols
!         write(*, '(3(f14.5, 2x), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values, i, &
!         sparse_matrix_a%columns(i)%full_index  
!     end do
!     write(*,*) ""
    
!     write(*,*) "scaled vectors"
!     write(*,*) ""
!     do i=1,sparse_matrix_a%full_num_cols
!         write(*, '(3(f14.5, 2x))') sparse_matrix_a%get_values(i)
!     end do
!     write(*,*) ""

!     ! check if test failed
!     if (sparse_matrix_a%sparse_num_cols /= 5) then
!         test_failed = .true.
!     end if
!     do i=1, sparse_a%full_size
!         if (any(abs(sparse_matrix_a%get_values(i) - (vector_a(i)*vec)) > 1.0e-12)) then
!             test_failed = .true.
!         end if
!     end do   
    
!     if (test_failed) then
!         write(*,*) "BROADCAST VECTOR CROSS ELEMENT FAILED"
!         total_tests = total_tests + 1
!     else
!         write(*,*) "BROADCAST VECTOR CROSS ELEMENT PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) ""
!     write(*,*) ""


! !!!!!!!!!! TEST BROADCAST ELEMENT CROSS VECTOR  !!!!!!!!!!
!     write(*,*) "-------------TEST BROADCAST ELEMENT CROSS VECTOR--------------"
!     write(*,*) ""
   

!     ! check if test failed

!     if (test_failed) then
!         write(*,*) "broadcast element cross vector test  FAILED"
!         total_tests = total_tests + 1
!     else
!         write(*,*) "broadcast element cross vector test PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) "" 
!     write(*,*) ""


! !!!!!!!!!! TEST BROADCAST VECTOR DOT ELEMENT  !!!!!!!!!!
!     write(*,*) "-------------TEST BROADCAST VECTOR DOT ELEMENT--------------"
!     write(*,*) ""
   

!     ! check if test failed

!     if (test_failed) then
!         write(*,*) "broadcast vector dot element test FAILED"
!         total_tests = total_tests + 1
!     else
!         write(*,*) "broadcast vector dot element test PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) "" 
!     write(*,*) ""


!     !!!!!!!!!! TEST BROADCAST ELEMENT TIMES SCALAR  !!!!!!!!!!
!     write(*,*) "-------------TEST BROADCAST ELEMENT TIMES SCALAR--------------"
!     write(*,*) ""
   

!     ! check if test failed

!     if (test_failed) then
!         write(*,*) "broadcast element times scalar test FAILED"
!         total_tests = total_tests + 1
!     else
!         write(*,*) "broadcast element times scalar test PASSED"
!         passed_tests = passed_tests + 1
!         total_tests = total_tests + 1
!     end if
!     test_failed = .false.
!     write(*,*) "" 
!     write(*,*) ""




!!!!!!!!!!!!!! SPARSE MATRIX TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------SPARSE MATRIX TEST RESULTS--------------"
    write(*,*) ""
    write(*,'(I15,a14)') total_tests - passed_tests, " tests FAILED"
    write(*,*) ""
    write(*,'(I4,a9,I2,a14)') passed_tests, " out of ", total_tests, " tests PASSED"

    



end program sparse_matrix_test