program sparse_matrix_test
    ! tests sparse vector and sparse matrix operations
    
    use adjoint_mod
    
    implicit none

    type(sparse_vector) :: sparse_a, sparse_b, sparse_c, sparse_vector_result 
    type(sparse_matrix) :: sparse_matrix_a, sparse_matrix_b, sparse_matrix_result, sparse_matrix_a_copy

    real :: scalar
    real,dimension(3) :: vec, values, add_values, check
    real,dimension(:),allocatable :: vector_a, vector_b, vector_c, vector_result, expanded_vector_result, &
    residual_vector
    real,dimension(:,:),allocatable :: scaled_vecs, full_matrix_a, full_matrix_b, expanded_matrix_a,&
    expanded_matrix_b, residual, full_matrix_result, expanded_matrix_result
    integer :: i, passed_tests, total_tests, nonzeros, old_size, add_full_index, add_shift_index
    logical :: test_failed
    character(len=100),dimension(20) :: failure_log

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
    call sparse_a%init_from_full_vector(vector_a)
    call sparse_b%init_from_full_vector(vector_b)
    call sparse_c%init_from_full_vector(vector_c)
    
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! check if test failed
    if (any((full_matrix_a(:,2) - sparse_matrix_a%columns(1)%vector_values > 1.0e-12)) &
    .or. any((full_matrix_a(:,3) - sparse_matrix_a%columns(2)%vector_values > 1.0e-12)) &
    .or. any((full_matrix_a(:,5) - sparse_matrix_a%columns(3)%vector_values > 1.0e-12)) &
    .or. any((full_matrix_a(:,6) - sparse_matrix_a%columns(4)%vector_values > 1.0e-12)) &
    .or. any((full_matrix_a(:,8) - sparse_matrix_a%columns(5)%vector_values > 1.0e-12)) &
    .or. any((full_matrix_a(:,11) - sparse_matrix_a%columns(6)%vector_values > 1.0e-12))) then    
    
        test_failed = .true.
    else 
        test_failed = .false.
        
    end if


    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "init from sparse vectors FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
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
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "init from sparse vectors FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "init from sparse vectors PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .true.
    write(*,*) "" 
    write(*,*) "" 



!!!!!!!!!! TEST INIT FROM SPARSE MATRIX !!!!!!!!!!
    write(*,*) "-------------TEST INIT FROM SPARSE MATRIX--------------"
    write(*,*) ""
    
    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a                       sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! init from sparse matrix
    write(*,*) "init a new sparse matrix from sparse_matrix_a "
    call sparse_matrix_a_copy%init_from_sparse_matrix(sparse_matrix_a)

    ! write sparse matrix
    write(*,*) ""
    write(*,*) "         sparse_matrix a COPY                       sparse_index      full_index"
    do i=1,sparse_matrix_a_copy%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a_copy%columns(i)%vector_values(:), &
        i, sparse_matrix_a_copy%columns(i)%full_index
    end do
    write(*,*) ""
    

    ! check if test failed
    do i = 1,sparse_matrix_a_copy%sparse_num_cols
        
        if ((any(sparse_matrix_a%columns(i)%vector_values /= sparse_matrix_a_copy%columns(i)%vector_values)) &
        .or. (sparse_matrix_a%columns(i)%full_index /= sparse_matrix_a_copy%columns(i)%full_index)) then    

            test_failed = .true.
        else 
            test_failed = .false.
            
        end if

    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "init from sparse matrix FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "init from sparse matrix PASSED"
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
        write(*,'(3(f12.5), 12x, I5, 12x, I10)') sparse_matrix_a%columns(i)%vector_values(:), &
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! result of test
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "increase size (sparse matrix) test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
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
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "add element (sparse matrix) test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
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
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "get element values test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(:), &
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
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "set values (sparse matrix) test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "set values (sparse matrix) test PASSED"
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(:), &
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
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "compress (sparse matrix) test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "compress (sparse matrix) test PASSED"
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
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "expand sparse matrix FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "expand sparse matrix PASSED"
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

    ! write full matrix b
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
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! write sparse matrix b
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(:), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! add sparse matrix b to sparse matrix a
    write(*,*) ""
    write(*,*) "add:  sparse_matrix_a + sparse_matrix_b...."
    write(*,*) ""

    call sparse_matrix_a%sparse_add(sparse_matrix_b)

    !write(*,*) "sparse b column 11 get after after call:", sparse_matrix_b%get_values(11)

    write(*,*) "     new sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! expand sparse_a
    expanded_matrix_a = sparse_matrix_a%expand()
    residual = full_matrix_a - expanded_matrix_a
    
    ! write results
    write(*,*) "residual =  new full matrix_a - new expanded_matrix_a"
    write(*,*) ""
    write(*,*) "           expanded_matrix_a                        residual"
    do i=1,sparse_matrix_b%full_num_cols
        write(*, '(3(f10.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_a(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_a%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "sparse add (matrix) FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "sparse add (matrix) PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""

    

!!!!!!!!!! TEST SPARSE SUBTRACT (MATRIX, a - b and b - a)  !!!!!!!!!!
    write(*,*) "-------------TEST SPARSE SUBTRACT (MATRIX) (a - b) --------------"
    write(*,*) ""
    
    write(*,*) ""
    write(*,*) "return sparse_matrix_a and full_matrix_a to original values (init from sparse vectors)"
    write(*,*) ""

    ! reset values
    deallocate(sparse_matrix_a%columns)

    !write(*,*) "deallocation successful" 
    call sparse_matrix_a%init_from_sparse_vectors(sparse_a, sparse_b, sparse_c)

    ! update full_matrix_a
    full_matrix_a = sparse_matrix_a%expand()

    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! write full matrix b
    write(*,*) " full matrix b:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_b(:,i)
    end do
    write(*,*) ""

    ! add full matrix b to full matrix c
    write(*,*) "subtract: full matrix a - full matrix b"
    write(*,*) ""
    do i=1,11 
        full_matrix_a(:,i) = full_matrix_a(:,i) - full_matrix_b(:,i)
    end do 

    ! display new full matrix a
    write(*,*) " new full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! add full matrix b to full matrix c
    
    write(*,*) ""
    write(*,*) "display sparse matrices before subtracting"
    write(*,*) ""

    ! write sparse matrix a
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! write sparse matrix b
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(:), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    ! subrac: sparse matrix a - sparse matrix b
    write(*,*) ""
    write(*,*) "subtract sparse_matrix_a - sparse_matrix_b...."
    write(*,*) ""

    call sparse_matrix_a%sparse_subtract(sparse_matrix_b)

    ! write new sparse a results
    write(*,*) "     new sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""

    ! expand sparse_a
    expanded_matrix_a = sparse_matrix_a%expand()
    residual = full_matrix_a - expanded_matrix_a
    
    ! write results
    write(*,*) "residual =  new full matrix_a - new expanded_matrix_a"
    write(*,*) ""
    write(*,*) "           expanded_matrix_a                        residual"
    do i=1,sparse_matrix_a%full_num_cols
        write(*, '(3(f10.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_a(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_a%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "sparse subtract (matrix) (a - b) FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "sparse subtract (matrix) (a - b) PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""

    write(*,*) "-------------TEST SPARSE SUBTRACT (MATRIX) (b - a) --------------"
    write(*,*) ""
    
    write(*,*) ""
    write(*,*) "return sparse_matrix_a and full_matrix_a to original values (init from sparse vectors)"
    write(*,*) ""

    ! reset values
    deallocate(sparse_matrix_a%columns)

    !write(*,*) "deallocation successful" 
    call sparse_matrix_a%init_from_sparse_vectors(sparse_a, sparse_b, sparse_c)

    ! update full_matrix_a
    full_matrix_a = sparse_matrix_a%expand()

    ! write full matrix b
    write(*,*) " full matrix b:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_b(:,i)
    end do
    write(*,*) ""

    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""

    ! subtract: full matrix b - full matrix a
    write(*,*) "subtract: full matrix b - full matrix a"
    write(*,*) ""
    do i=1,11 
        full_matrix_b(:,i) = full_matrix_b(:,i) - full_matrix_a(:,i)
    end do 

    ! display new full matrix b
    write(*,*) " new full matrix b:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_b(:,i)
    end do
    write(*,*) ""

    
    write(*,*) ""
    write(*,*) "display sparse matrices before subtracting"
    write(*,*) ""

    ! write sparse matrix b
    write(*,*) ""
    write(*,*) "         sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(:), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    ! write sparse matrix a
    write(*,*) ""
    write(*,*) "         sparse_matrix a                          sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    ! subrac: sparse matrix b - sparse matrix a
    write(*,*) ""
    write(*,*) "subtract sparse_matrix_b - sparse_matrix_a...."
    write(*,*) ""

    call sparse_matrix_b%sparse_subtract(sparse_matrix_a)

    write(*,*) "     new sparse_matrix b                          sparse_index      full_index"
    do i=1,sparse_matrix_b%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_b%columns(i)%vector_values(:), &
        i, sparse_matrix_b%columns(i)%full_index
    end do
    write(*,*) ""

    ! expand sparse_b
    expanded_matrix_b = sparse_matrix_b%expand()
    residual = full_matrix_b - expanded_matrix_b
    
    ! write results
    write(*,*) "residual =  new full matrix_b - new expanded_matrix_b"
    write(*,*) ""
    write(*,*) "           expanded_matrix_b                        residual"
    do i=1,sparse_matrix_b%full_num_cols
        write(*, '(3(f10.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_b(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_a%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "sparse subtract (matrix) (a - b) FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "sparse subtract (matrix) (a - b) PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""

    


!!!!!!!!!! TEST BROADCAST VECTOR-CROSS-ELEMENT   !!!!!!!!!!

    write(*,*) "-------------TEST BROADCAST VECTOR-CROSS-ELEMENT--------------"
    write(*,*) ""

    
    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""


    vec = (/1.0,-2.0,3.0/)
    write(*,'(A, 3(f10.5, 2x))') "vector: ", vec
    write(*,*) ""

    ! perform cross product manually
    write(*,*) "broadcast vector-cross-element in full matrix"
    write(*,*) ""
    allocate(full_matrix_result(3,11))
    do i=1,11
        full_matrix_result(:,i) = cross(vec(:),full_matrix_a(:,i)) 
    end do

    ! write full matrix result
    write(*,*) " full matrix result:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_result(:,i)
    end do
    write(*,*) ""

    
    ! perform cross product
    write(*,*) "broadcast vector crossed with each sparse element...."
    
    sparse_matrix_result = sparse_matrix_a%broadcast_vector_cross_element(vec)
    
    write(*,*) ""
    write(*,*) "resulting sparse_matrix:"
    write(*,*) ""
    write(*,*) "     sparse_matrix_result                      sparse_index      full_index"
    do i=1,sparse_matrix_result%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_result%columns(i)%vector_values(:), &
        i, sparse_matrix_result%columns(i)%full_index
    end do
    write(*,*) ""
    
    
    ! expand sparse_b
    expanded_matrix_result = sparse_matrix_result%expand()
    residual = full_matrix_result - expanded_matrix_result
    
    ! write results
    write(*,*) "residual =  full matrix_result - expanded_matrix_result"
    write(*,*) ""
    write(*,*) "           expanded_matrix_result                          residual"
    do i=1,sparse_matrix_result%full_num_cols
        write(*, '(3(f10.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_result(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_result%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "broadcast vector-cross-element FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "broadcast vector-cross-element PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!! TEST BROADCAST ELEMENT-CROSS-VECTOR   !!!!!!!!!!

    write(*,*) "-------------TEST BROADCAST ELEMENT-CROSS-VECTOR--------------"
    write(*,*) ""

    
    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""


    vec = (/1.0,-2.0,3.0/)
    write(*,'(A, 3(f10.5, 2x))') "vector: ", vec
    write(*,*) ""

    ! perform cross product manually
    write(*,*) "broadcast element-cross-vector in full matrix"
    write(*,*) ""
    
    do i=1,11
        full_matrix_result(:,i) = cross(full_matrix_a(:,i), vec(:)) 
    end do

    ! write full matrix result
    write(*,*) " full matrix result:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_result(:,i)
    end do
    write(*,*) ""

    
    ! perform cross product
    write(*,*) "broadcast each sparse element crossed with given vector...."
    
    sparse_matrix_result = sparse_matrix_a%broadcast_element_cross_vector(vec)
    
    write(*,*) ""
    write(*,*) "resulting sparse_matrix:"
    write(*,*) ""
    write(*,*) "     sparse_matrix_result                      sparse_index      full_index"
    do i=1,sparse_matrix_result%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_result%columns(i)%vector_values(:), &
        i, sparse_matrix_result%columns(i)%full_index
    end do
    write(*,*) ""
    
    
    ! expand sparse_b
    expanded_matrix_result = sparse_matrix_result%expand()
    residual = full_matrix_result - expanded_matrix_result
    
    ! write results
    write(*,*) "residual =  full matrix_result - expanded_matrix_result"
    write(*,*) ""
    write(*,*) "           expanded_matrix_result                          residual"
    do i=1,sparse_matrix_result%full_num_cols
        write(*, '(3(f12.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_result(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_result%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "broadcast element-cross-vector FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "broadcast element-cross-vector PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!! TEST BROADCAST VECTOR DOT ELEMENT  !!!!!!!!!!
    write(*,*) "-------------TEST BROADCAST VECTOR-DOT-ELEMENT--------------"
    write(*,*) ""

    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""


    vec = (/1.0,-2.0,3.0/)
    write(*,'(A, 3(f10.5, 2x))') "vector: ", vec
    write(*,*) ""

    ! perform cross product manually
    write(*,*) "broadcast vector-dot-element in full matrix"
    write(*,*) ""
    allocate(vector_result(11))
    do i=1,11
        vector_result(i) = inner(full_matrix_a(:,i), vec(:)) 
    end do

    ! write full matrix result
    write(*,*) " resulting full VECTOR (broadcasting dot product results in an array of scalars):"
    do i=1,11
        write(*, '(f12.5)') vector_result(i)
    end do
    write(*,*) ""

    
    ! perform cross product
    write(*,*) "broadcast vector dotted with each sparse element ...."
    
    sparse_vector_result = sparse_matrix_a%broadcast_vector_dot_element(vec)
    
    write(*,*) ""
    write(*,*) "resulting sparse_vector:"
    write(*,*) ""
    write(*,*) "  sparse_vector_result          sparse_index      full_index"
    do i=1,sparse_vector_result%sparse_size
        write(*,'(f12.5, 18x, I5, 12x, I5)') sparse_vector_result%elements(i)%value, &
        i, sparse_vector_result%elements(i)%full_index
    end do
    write(*,*) ""
    
    
    ! expand sparse_b
    expanded_vector_result = sparse_vector_result%expand()
    allocate(residual_vector(sparse_vector_result%full_size))
    residual_vector = vector_result - expanded_vector_result
    
    ! write results
    write(*,*) "residual =  full vector_result - expanded_matrix_result"
    write(*,*) ""
    write(*,*) "expanded_vector_result         residual"
    do i=1,sparse_vector_result%full_size
        write(*, '(f12.6 ,18x, f10.6)') expanded_vector_result(i), residual_vector(i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_vector_result%full_size
        if (abs(residual_vector(i)) > 1.0e-12) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "broadcast vector-dot-element test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "broadcast vector-dot-element test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""


!!!!!!!!!! TEST BROADCAST ELEMENT-TIMES-SCALAR  !!!!!!!!!!
    write(*,*) "-------------TEST BROADCAST ELEMENT-TIMES-SCALAR--------------"
    write(*,*) ""
   

    ! write full matrix a
    write(*,*) " full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""


    scalar = 7.812
    write(*,*) "scalar: ", scalar
    write(*,*) ""

    ! perform scalar multiplication manually
    write(*,*) "broadcast scalar multiplication in full matrix"
    write(*,*) ""
    
    do i=1,11
        full_matrix_a(:,i) = full_matrix_a(:,i)*scalar 
    end do

    ! write full matrix a
    write(*,*) " new full matrix a:"
    do i=1,11
        write(*, '(3(f12.5, 2x))') full_matrix_a(:,i)
    end do
    write(*,*) ""
    
    ! perform cross product
    write(*,*) "broadcast each sparse element times scalar ...."
    
    call sparse_matrix_a%broadcast_element_times_scalar(scalar)
    
    write(*,*) ""
    write(*,*) "resulting sparse_matrix:"
    write(*,*) ""
    write(*,*) "     sparse_matrix_a                      sparse_index      full_index"
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*,'(3(f12.5), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values(:), &
        i, sparse_matrix_a%columns(i)%full_index
    end do
    write(*,*) ""
    
    
    ! expand sparse_a
    expanded_matrix_a = sparse_matrix_a%expand()
    residual = full_matrix_a - expanded_matrix_a
    
    ! write results
    write(*,*) "residual =  new full matrix_a - new expanded_matrix_a"
    write(*,*) ""
    write(*,*) "        new expanded_matrix_a                        residual"
    do i=1,sparse_matrix_a%full_num_cols
        write(*, '(3(f10.6, 2x),3x, 3(f10.6, 2x))') expanded_matrix_a(:,i), residual(:,i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_matrix_a%full_num_cols
        if (any(abs(residual(:,i)) > 1.0e-12)) then
            test_failed = .true.
            exit
        else 
            test_failed = .false.
        end if
    end do
    if (test_failed) then
        total_tests = total_tests + 1
        failure_log(total_tests-passed_tests) = "broadcast element-times-scalar test FAILED"
        write(*,*) failure_log(total_tests-passed_tests)
    else
        write(*,*) "broadcast element-times-scalar test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""




!!!!!!!!!!!!!! SPARSE MATRIX TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------SPARSE MATRIX TEST RESULTS--------------"
    write(*,*) ""
    write(*,'(I15,a14)') total_tests - passed_tests, " tests FAILED"
    write(*,*) ""
    write(*,'(I4,a9,I2,a14)') passed_tests, " out of ", total_tests, " tests PASSED"
    
    if (passed_tests < total_tests)then
        write(*,*) ""
        write(*,*) "Failure Log:"
        do i=1,total_tests-passed_tests
            write(*,*) failure_log(i)
        end do
    end if
    
    write(*,*) ""
    write(*,*) "Program Complete"
    write(*,*) ""
    



end program sparse_matrix_test