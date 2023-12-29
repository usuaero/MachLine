program sparse_vector_test
    ! tests sparse vector and sparse matrix operations
    
    use adjoint_mod
    
    implicit none

    type(sparse_vector) :: sparse_a, sparse_b, sparse_ab, sparse_a_minus_b, sparse_b_minus_a
    type(sparse_matrix) :: sparse_matrix_a, sparse_matrix_b
    real,dimension(:),allocatable :: vector_a, vector_b, vector_ab, expanded_a, expanded_ab, expanded_a_minus_b, &
    expanded_b_minus_a, vector_a_minus_b, vector_b_minus_a, residual
    real,dimension(3) :: vec
    real,dimension(:,:),allocatable :: scaled_vecs
    integer :: i, full_size_a, add_index, add_shift_index, passed_tests, total_tests
    real :: add_value, value
    logical :: test_failed

    write(*,*) ""
    write(*,*) ""
    write(*,*) "-----------------------------SPARSE VECTOR TEST------------------------------"
    write(*,*) ""
    write(*,*) ""
    test_failed = .false.
    passed_tests = 0
    total_tests = 0
    


!!!!!!!!!!!! INITIALIZE VECTOR !!!!!!!!!!!!!!
    write(*,*) "--------INITIALIZE A REGULAR VECTOR---------"
    write(*,*) ""

    vector_a = (/ 1.0, 0.0, 11.0, 0.0, 87.0, -1.0, 0.0, 0.002121, 4.0, 0.0, 3.14 /)
    ! allocate(vector_a(100), source=0.0)
    ! do i= 1,9
    !     vector_a(i*10) = real(i) 
    ! end do

    ! write full vector
    write(*,*) "full vector_a"
    write(*, '(f14.10, 2x)') vector_a
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!! TEST INIT SPARSE VECTOR !!!!!!!!!!
    write(*,*) "-------------TEST INIT SPARSE VECTOR--------------"
    write(*,*) ""
    write(*,*) "initialize a sparse vector"
    write(*,*) ""

    ! initialize and compress vector_a into sparse_a
    call sparse_a%init(vector_a)
    
    ! write results
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_a%full_size
        if (abs(sparse_a%get_value(i) - vector_a(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "init sparse vector FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "init sparse vector PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) "" 


!!!!!!!!!! TEST GET_VALUE (SPARSE VECTOR)  !!!!!!!!!!
    write(*,*) "-------------TEST GET_VALUE (SPARSE VECTOR)--------------"
    write(*,*) ""
    write(*,*) "Get value of full index 1" 
    value = sparse_a%get_value(1)
    if (abs(sparse_a%get_value(1)-vector_a(1)) > 1.0e-12) then
        test_failed = .true.
    end if 
    write(*,*) "value = ", value
    write(*,*) ""

    write(*,*) "Get value of full index 2" 
    value = sparse_a%get_value(2)
    if (abs(sparse_a%get_value(2)-vector_a(2)) > 1.0e-12) then
        test_failed = .true.
    end if
    write(*,*) "value = ", value
    write(*,*) ""
    
    write(*,*) "Get value of full index 3" 
    value = sparse_a%get_value(3)
    if (abs(sparse_a%get_value(3)-vector_a(3)) > 1.0e-12) then
        test_failed = .true.
    end if
    write(*,*) "value = ", value
    write(*,*) ""

    write(*,*) "Get value of full index 8" 
    value = sparse_a%get_value(8)
    if (abs(sparse_a%get_value(8)-vector_a(8)) > 1.0e-12) then
        test_failed = .true.
    end if
    write(*,*) "value = ", value
    write(*,*) ""

    write(*,*) "Get value of full index 9 " 
    value = sparse_a%get_value(9)
    if (abs(sparse_a%get_value(9)-vector_a(9)) > 1.0e-12) then
        test_failed = .true.
    end if
    write(*,*) "value = ", value
    write(*,*) ""
    
    ! check if test failed
    if (test_failed) then
        write(*,*) "get value (sparse vector) test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "get value (sparse vector) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) "" 


!!!!!!!!!! TEST EXPAND SPARSE VECTOR  !!!!!!!!!!
    write(*,*) "-------------TEST EXPAND SPARSE VECTOR--------------"
    ! expand sparse_a
    expanded_a = sparse_a%expand()
    residual = vector_a - expanded_a
    
    ! write results
    write(*,*) "expanded_a          residual"
    do i=1,sparse_a%full_size
        write(*, '(2(f14.10, 2x))') expanded_a(i), residual(i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_a%full_size
        if (abs(residual(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "expand sparse vector FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "expand sparse vector PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""


!!!!!!!!!! TEST INCREASE_SIZE (SPARSE VECTOR)  !!!!!!!!!!
    write(*,*) "-------------TEST INCREASE_SIZE (SPARSE VECTOR)--------------"
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
    write(*,*) "return to original size"
    write(*,*) ""

    deallocate(sparse_a%elements)
    call sparse_a%init(vector_a)

    ! write sparse_a again
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""
    
    ! result of test
    if (test_failed) then
        write(*,*) "increase size (sparse vector) test FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "increase size (sparse vector) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""


!!!!!!!!!! TEST ADD_ELEMENT (SPARSE VECTOR)  !!!!!!!!!!
    write(*,*) "-------------TEST ADD ELEMENT (SPARSE VECTOR)--------------"
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
    write(*,*) ""

    ! check if test failed
    vector_a(4) = 5.1
    do i=1,sparse_a%full_size
        if (abs(sparse_a%get_value(i) - vector_a(i)) > 1.0e-12)  then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "add element (sparse vector) test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "add element (sparse vector) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""



!!!!!!!!!! TEST SET_VALUE (SPARSE VECTOR)  !!!!!!!!!!
    write(*,*) "-------------TEST SET_VALUE (SPARSE VECTOR)--------------"
    write(*,*) ""

    write(*,*) "Set value of full index 2 to 3.14" 
    call sparse_a%set_value(3.14, 2)
    write(*,*) "Set value of full index 3 to -1.234" 
    call sparse_a%set_value(-1.234, 3)
    write(*,*) "Set value of full index 4 to 0.0" 
    call sparse_a%set_value(0.0, 4)
    write(*,*) "Set value of full index 6,7,8,9,10,11 to 0.0" 
    call sparse_a%set_value(0.0, 6)
    call sparse_a%set_value(0.0, 7)
    call sparse_a%set_value(0.0, 8)
    call sparse_a%set_value(0.0, 9)
    call sparse_a%set_value(0.0, 10)
    call sparse_a%set_value(0.0, 11)
    write(*,*) ""
    
    write(*,*) "resulting vector:"
    write(*,*) ""
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""

    ! check if test failed
    vector_a = (/ 1.0, 3.14, -1.234, 0.0, 87.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    do i=1,sparse_a%full_size
        if (abs(sparse_a%get_value(i) - vector_a(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "set value (sparse vector) test  FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "set value (sparse vector) test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""



!!!!!!!!!! TEST COMPRESS SPARSE VECTOR  !!!!!!!!!!
    write(*,*) "-------------TEST COMPRESS SPARSE VECTOR--------------"
    write(*,*) ""
    write(*,*) "compressing sparse vector:"
    write(*,*) ""
    
    call sparse_a%compress()
    
    write(*,*) "sparse_a value    sparse_index      full_index"
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""
    
     ! check if test failed
    do i=1,3
        if (abs(sparse_a%elements(i)%value - vector_a(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if ((sparse_a%sparse_size /= 4) .or. (sparse_a%elements(4)%full_index /= 5)) then
        test_failed = .true.
    end if
    if (test_failed) then
        write(*,*) "compress sparse vector test FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "compress sparse vector test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) "" 
    write(*,*) ""
    deallocate(sparse_a%elements)


!!!!!!!!!! TEST SPARSE ADD (VECTORS)  !!!!!!!!!!
    write(*,*) "-------------TEST SPARSE ADD (VECTORS)--------------"
    write(*,*) ""
    write(*,*) "vectors to be added:"
    vector_a = (/ 1.0, 0.0, 11.0, 0.0, 87.0, 0.0, 0.0, 0.0, 4.0, 0.0, 3.14 /)
    vector_b = (/ 1.0, 0.0, 1.0, 0.0, 1.0, 7.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    vector_ab = (/ 2.0, 0.0, 12.0, 0.0, 88.0, 7.0, 0.0, 0.0, 5.0, 0.0, 4.14 /)
    
    ! write  full vectors
    write(*,*) ""
    write(*,*) "full vector_a  +  full vector_b   =   vector_ab"
    do i=1,11
        write(*, '(3(f14.10, 2x))') vector_a(i), vector_b(i), vector_ab(i)
    end do
    write(*,*) ""
    
    ! compress vectors
    write(*,*) "compressing..."
    write(*,*) ""
    call sparse_a%init(vector_a)
    call sparse_b%init(vector_b) 
    
    ! write compressed vectors
    write(*,*) "sparse_a value    sparse_index      full_index "
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a%elements(i)%value, i, sparse_a%elements(i)%full_index
    end do
    write(*,*) ""
    write(*,*) "sparse_b value    sparse_index      full_index "
    do i=1,sparse_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_b%elements(i)%value, i, sparse_b%elements(i)%full_index
    end do
    write(*,*) ""
    
    ! add and display 
    write(*,*) "adding sparse_b to sparse_a...."
    write(*,*) ""

    sparse_ab = sparse_a%sparse_add(sparse_b) !sparse_b%sparse_add(sparse_a)
    write(*,*) "sparse_ab value    sparse_index      full_index "
    do i=1,sparse_ab%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_ab%elements(i)%value, i, sparse_ab%elements(i)%full_index
    end do
    write(*,*) ""

    ! expand sparse_a
    expanded_ab = sparse_ab%expand()
    residual = vector_ab - expanded_ab
    
    ! write results
    write(*,*) "residual = vector_ab - expanded_ab"
    write(*,*) ""
    write(*,*) "expanded_ab          residual"
    do i=1,sparse_ab%full_size
        write(*, '(2(f14.10, 2x))') expanded_ab(i), residual(i)
    end do
    write(*,*) ""

    ! check if test failed
    do i=1,sparse_a%full_size
        if (abs(residual(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "sparse add (vectors) FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "sparse add (vectors) PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""

!!!!!!!!!! TEST SPARSE SUBTRACT (VECTORS) a-b and b-a  !!!!!!!!!!
    write(*,*) "-------------TEST SPARSE SUBTRACT (VECTORS) a-b and b-a--------------"
    write(*,*) ""
    write(*,*) "vectors to be subtracted:"
    vector_a = (/ 1.0, 0.0, 11.0, 0.0, 87.0, 0.0, 0.0, 0.0, 4.0, 0.0, 3.14 /)
    vector_b = (/ 1.0, 0.0, 1.0, 0.0, 1.0, 7.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    vector_a_minus_b = (/ 0.0, 0.0, 10.0, 0.0, 86.0, -7.0, 0.0, 0.0, 3.0, 0.0, 2.14 /)
    vector_b_minus_a = (/ 0.0, 0.0, -10.0, 0.0, -86.0, 7.0, 0.0, 0.0, -3.0, 0.0, -2.14 /)
    
    ! write  full vectors
    write(*,*) ""
    write(*,*) " vector_a     -     vector_b   =   vector_a_minus_b"
    do i=1,11
        write(*, '(3(f14.10, 2x))') vector_a(i), vector_b(i), vector_a_minus_b(i)
    end do
    write(*,*) ""
    
    ! check if test failed
    do i=1,sparse_a%full_size
        if (abs(residual(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "sparse subtract (vector) a - b FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "sparse subtract (vector) a - b PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""
    
    
    ! subtract and display 
    write(*,*) "subtracting sparse_a minus sparse_b...."
    write(*,*) ""
    sparse_a_minus_b = sparse_a%sparse_subtract(sparse_b)
    write(*,*) "sparse_a_minus_b value    sparse_index      full_index "
    do i=1,sparse_a_minus_b%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_a_minus_b%elements(i)%value, i, sparse_a_minus_b%elements(i)%full_index
    end do
    write(*,*) ""

    ! expand 
    expanded_a_minus_b = sparse_a_minus_b%expand()
    residual = vector_a_minus_b - expanded_a_minus_b
    
    ! write results
    write(*,*) "residual = vector_a_minus_b - expanded_a_minus_b"
    write(*,*) ""
    write(*,*) "expanded_a_minus_b          residual"
    do i=1,sparse_a_minus_b%full_size
        write(*, '(2(f14.10, 2x))') expanded_a_minus_b(i), residual(i)
    end do
    write(*,*) ""
    
    ! write  full vectors
    write(*,*) "-------now do vector b minus vector a-------"
    write(*,*) ""
    write(*,*) " vector_b     -     vector_a   =   vector_b_minus_a"
    do i=1,11
        write(*, '(3(f14.10, 2x))') vector_b(i), vector_a(i), vector_b_minus_a(i)
    end do
    write(*,*) "" 
    
    
    ! subtract and display 
    write(*,*) "subtracting sparse_b minus sparse_a...."
    write(*,*) ""
    sparse_b_minus_a = sparse_b%sparse_subtract(sparse_a)
    write(*,*) "sparse_a_minus_b value    sparse_index      full_index "
    do i=1,sparse_b_minus_a%sparse_size
        write(*,'(f14.10, 12x, I5, 12x, I5)') sparse_b_minus_a%elements(i)%value, i, sparse_b_minus_a%elements(i)%full_index
    end do
    write(*,*) ""

    ! expand 
    expanded_b_minus_a = sparse_b_minus_a%expand()
    residual = vector_b_minus_a - expanded_b_minus_a
    
    ! write results
    write(*,*) "residual = vector_b_minus_a - expanded_b_minus_a"
    write(*,*) ""
    write(*,*) "expanded_b_minus_a          residual"
    do i=1,sparse_b_minus_a%full_size
        write(*, '(2(f14.10, 2x))') expanded_b_minus_a(i), residual(i)
    end do
    write(*,*) ""
    
    ! check if test failed
    do i=1,sparse_a%full_size
        if (abs(residual(i)) > 1.0e-12) then
            test_failed = .true.
        end if
    end do
    if (test_failed) then
        write(*,*) "sparse subtract (vector) b - a FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "sparse subtract (vector) b - a PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""
    


!!!!!!!!!! TEST BROADCAST ELEMENT TIMES VECTOR    !!!!!!!!!!

    write(*,*) "-------------TEST BROADCAST ELEMENT TIMES VECTOR --------------"
    write(*,*) ""
    vec = (/1.0,-2.0,3.0/)
    write(*,'(A, 3(f10.5, 2x))') "new vector to be scaled: ", vec
    write(*,*) ""
    write(*,*) "unscaled vectors"
    do i=1,sparse_a%full_size
        write(*, '(3(f14.5, 2x))') vec
    end do
    write(*,*) ""
    
    ! scale and write results
    sparse_matrix_a = sparse_a%broadcast_element_times_vector(vec)
    
    write(*,*) "scale vec by each sparse element...."
    write(*,*) ""
    write(*,*) "resulting sparse_matrix:"
    write(*,*) ""
    write(*,*) "                   scaled vector                    sparse_index       &
              full_index "
    write(*,*) ""
    do i=1,sparse_matrix_a%sparse_num_cols
        write(*, '(3(f14.5, 2x), 12x, I5, 12x, I5)') sparse_matrix_a%columns(i)%vector_values, i, &
        sparse_matrix_a%columns(i)%full_index  
    end do
    write(*,*) ""
    
    write(*,*) "scaled vectors"
    write(*,*) ""
    do i=1,sparse_matrix_a%full_num_cols
        write(*, '(3(f14.5, 2x))') sparse_matrix_a%get_values(i)
    end do
    write(*,*) ""

    ! check if test failed
    if (sparse_matrix_a%sparse_num_cols /= 5) then
        test_failed = .true.
    end if
    do i=1, sparse_a%full_size
        if (any(abs(sparse_matrix_a%get_values(i) - (vector_a(i)*vec)) > 1.0e-12)) then
            test_failed = .true.
        end if
    end do   
    
    if (test_failed) then
        write(*,*) "broadcast element times vector test FAILED"
        total_tests = total_tests + 1
    else
        write(*,*) "broadcast element times vector test PASSED"
        passed_tests = passed_tests + 1
        total_tests = total_tests + 1
    end if
    test_failed = .false.
    write(*,*) ""
    write(*,*) ""


!!!!!!!!!!!!!! SPARSE VECTOR TEST RESULTS!!!!!!!!!!!!!
    write(*,*) "-------------SPARSE VECTOR TEST RESULTS--------------"
    write(*,*) ""
    write(*,*) total_tests - passed_tests, " tests FAILED"
    write(*,*) ""
    write(*,*) passed_tests, " out of ", total_tests, " tests PASSED"

    



end program sparse_vector_test