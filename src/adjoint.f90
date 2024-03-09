! module for calculating the gradient of CF, CM with respect to mesh points via the adjoint method
module adjoint_mod

    !use linked_list_mod
    use math_mod
    use helpers_mod
    

    implicit none

    type sparse_vector_element
    
        real :: value   ! non-zero value of an element in a sparse vector
        integer :: full_index     ! index of value in the full matrix
        
    end type sparse_vector_element


    type sparse_matrix_element
    
        real,dimension(3) :: vector_values ! values is used for sparse 3 matrices             
        integer :: full_index     ! index of value in the full matrix
    
    end type sparse_matrix_element


    type sparse_vector
        ! used to save memory when we have many sparse vectors (most of the vector elements are 0) 
        ! "full vector" has all the zero and nonzero values
        ! "sparse vector" has only the nonzero values, with their corresponding positions in the full_vector

        integer :: sparse_size,full_size    ! sparse_size is the number of non-zero numbers, full_size is full vector
        type(sparse_vector_element),dimension(:),allocatable :: elements        ! non zero values in sparse vector

        contains
            procedure :: init => sparse_vector_init   
            procedure :: init_from_full_vector => sparse_vector_init_from_full_vector
            procedure :: init_from_sparse_vector => sparse_vector_init_from_sparse_vector
            
            procedure :: increase_size => sparse_vector_increase_size
            procedure :: add_element => sparse_vector_add_element

            procedure :: get_value => sparse_vector_get_value
            procedure :: set_value => sparse_vector_set_value

            procedure :: compress => sparse_vector_compress
            procedure :: expand => sparse_vector_expand

            procedure :: sparse_add => sparse_vector_sparse_add
            procedure :: sparse_subtract => sparse_vector_sparse_subtract
            
            procedure :: broadcast_element_times_vector => sparse_vector_broadcast_element_times_vector
            procedure :: broadcast_element_times_scalar => sparse_vector_broadcast_element_times_scalar

    end type sparse_vector
            
            
    type sparse_matrix
        ! used to save memory when we have matrices made of sparse vectors
        
        integer :: sparse_num_cols,full_num_cols    ! sparse_num_cols is the number of columns with a nonzero element
        type(sparse_matrix_element),dimension(:),allocatable :: columns       ! non zero values in sparse vector
        
            
        contains
            procedure :: init => sparse_matrix_init 
            procedure :: init_from_sparse_vectors => sparse_matrix_init_from_sparse_vectors
            procedure :: init_from_full_matrix => sparse_matrix_init_from_full_matrix
            procedure :: init_from_sparse_matrix => sparse_matrix_init_from_sparse_matrix

            procedure :: increase_size => sparse_matrix_increase_size
            procedure :: add_element => sparse_matrix_add_element

            procedure :: get_values => sparse_matrix_get_values
            procedure :: set_values => sparse_matrix_set_values

            procedure :: compress => sparse_matrix_compress
            procedure :: expand => sparse_matrix_expand
            
            procedure :: sparse_add => sparse_matrix_sparse_add
            
            procedure :: sparse_subtract => sparse_matrix_sparse_subtract
            
            
            procedure :: broadcast_vector_cross_element => sparse_matrix_broadcast_vector_cross_element
            procedure :: broadcast_element_cross_vector => sparse_matrix_broadcast_element_cross_vector
            
            procedure :: broadcast_vector_dot_element => sparse_matrix_broadcast_vector_dot_element
            procedure :: broadcast_element_times_scalar => sparse_matrix_broadcast_element_times_scalar
            
            procedure :: broadcast_matmul_3x3_times_element => sparse_matrix_broadcast_matmul_3x3_times_element
            procedure :: broadcast_matmul_element_times_3x3 => sparse_matrix_broadcast_matmul_element_times_3x3
            
            procedure :: split_into_sparse_vectors => sparse_matrix_split_into_sparse_vectors
            
            
            
            
    end type sparse_matrix
            
            
    type sparse_3D
            
            integer :: num_rows    ! sparse_num_cols is the number of columns with a nonzero element
            type(sparse_matrix),dimension(:), allocatable :: rows ! 1 row of a 3x3xN is a sparse_matrix 
            
            
            contains
                procedure :: init_from_sparse_matrices => sparse_3D_init_from_sparse_matrices
                procedure :: sparse_add_3 => sparse_3D_sparse_add_3
                procedure :: transpose_3 => sparse_3D_transpose_3

                procedure :: broadcast_matmul_3row_times_3x3 => sparse_3D_broadcast_matmul_3row_times_3x3
                procedure :: broadcast_matmul_3x3_times_3row => sparse_3D_broadcast_matmul_3x3_times_3row

                procedure :: broadcast_matmul_3row_times_3x1 => sparse_3D_broadcast_matmul_3row_times_3x1
                procedure :: broadcast_matmul_1x3_times_3row => sparse_3D_broadcast_matmul_1x3_times_3row

                procedure :: convert_to_sparse_vector_3x3 => sparse_3D_convert_to_sparse_vector_3x3
                

    end type sparse_3D


contains


!!!!!!!!!!!!!!!!!!!!! SPARSE VECTOR TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine sparse_vector_init(this,N)
        ! initializes a sparse vector with one element =0. if a derivative WRT X(beta), N = N_verts*3
        
        implicit none

        class(sparse_vector), intent(inout) :: this
        integer,intent(in) :: N

        ! specify number of design variables N (if WRT to X(beta), N = N_verts*3)
        this%full_size = N

        ! allocate and populate an initial sparse vector element
        allocate(this%elements(1))
        this%elements(1)%value = 0.0
        this%elements(1)%full_index = 1
        this%sparse_size = 1

    end subroutine sparse_vector_init



    subroutine sparse_vector_init_from_full_vector(this, full_vector)
        ! takes a real full vector that is sparse (has a lot of 0.0's) and compresses it to a sparse vector type

        implicit none 

        class(sparse_vector),intent(inout) :: this
        real,dimension(:),intent(inout) :: full_vector
        integer :: i,count,full_size, sparse_iter
        integer,dimension(:),allocatable :: indices 

        ! get size of vector
        full_size = size(full_vector)
        allocate(indices(full_size))

        ! count how many nonzero elements there are
        count = 0
        do i=1,full_size
            if (abs(full_vector(i)) > 1.0e-16 ) then
                count = count + 1
                indices(count) = i
            end if
        end do
        
        
        ! now that we have the number of non zero numbers, we can allocate the space for the sparse_vector
        this%sparse_size = count
        this%full_size = full_size
        allocate(this%elements(count))

        ! populate the sparse_vector elements
        do i=1,count
            this%elements(i)%value = full_vector(indices(i))
            this%elements(i)%full_index = indices(i)
        end do
        
        ! to save memory, you can deallocate the original array after converting to a sparse vector

    end subroutine sparse_vector_init_from_full_vector


    subroutine sparse_vector_init_from_sparse_vector(this, sparse_input)
        ! inits a sparse vector as a copy of the input sparse vector

        implicit none

        class(sparse_vector),intent(inout) :: this
        type(sparse_vector) :: sparse_input

        integer :: i

        ! copy info over
        this%sparse_size = sparse_input%sparse_size
        this%full_size = sparse_input%full_size
        
        ! allocate the same number of sparse vector elements
        allocate(this%elements(sparse_input%sparse_size))
        
        do i=1,this%sparse_size
            
            ! copy each element
            this%elements(i) = sparse_input%elements(i)

        end do    

    end subroutine sparse_vector_init_from_sparse_vector


    subroutine sparse_vector_increase_size(this)
        ! increases the size of the sparse_matrix by 1 sparse element

        implicit none
        class(sparse_vector),intent(inout) :: this

        type(sparse_vector) :: temp_vector
        integer :: new_size

        new_size = size(this%elements) + 1

        ! allocate a temporary array
        allocate(temp_vector%elements(new_size))

        ! copy array info
        temp_vector%elements(1:new_size-1) = this%elements(:)

        ! move allocation of temporary array and overwrite this
        call move_alloc(temp_vector%elements, this%elements) 
        
        ! update the sparse number of columns
        this%sparse_size = this%sparse_size + 1

    end subroutine sparse_vector_increase_size


    subroutine sparse_vector_add_element(this, value, full_index, shift_index)
        ! adds a sparse element to the current sparse vector
        ! could be made more efficient by checking if the shift index is closer to the beginning or 
        ! the end. if closer to the beginning, cshift, then pull back the first few elements until shift_index
        
        implicit none

        class(sparse_vector),intent(inout) :: this
        integer, intent(in) :: full_index, shift_index
        real, intent(in) :: value
        
        integer :: i
        
        ! increase the size of the array by 1
        call this%increase_size()

        ! starting from the last index of the old array (new_size-1), shift the element up one index.
        !do this up to and including the given shift index
        do i=this%sparse_size,shift_index,-1
            this%elements(i)%value = this%elements(i-1)%value
            this%elements(i)%full_index = this%elements(i-1)%full_index
        end do

        ! insert the new element into its proper position
        this%elements(shift_index)%value = value
        this%elements(shift_index)%full_index = full_index

           

    end subroutine sparse_vector_add_element


    function sparse_vector_get_value(this, full_index) result(value)
        ! given a full vector index, this function returns the full vector value at that index

        implicit none

        class(sparse_vector),intent(inout) :: this
        integer :: full_index, i

        real :: value 

        
        value = 1.0e-16
              
        ! if the sparse vector element has the same full index as the given full index value is set to
        ! that element's value, if not, the value stays zero 
        do i=1,this%sparse_size
            if (this%elements(i)%full_index == full_index) then
                value = this%elements(i)%value
                exit                
            end if  
        end do
       
    end function sparse_vector_get_value


    subroutine sparse_vector_set_value(this, value, full_index) 
        ! given a full_vector index and a value, this adds or updates the corresponding value in the sparse vector

        implicit none

        class(sparse_vector),intent(inout) :: this
        integer :: full_index
        real :: value
        
        integer :: i, shift_index
        logical :: new_sparse_element_needed
        
        ! initialize logical
        new_sparse_element_needed = .FALSE.

        ! check to see if a new sparse element is needed or if it can be updated
        do i=1,this%sparse_size
            if (this%elements(i)%full_index == full_index) then
                this%elements(i)%value = value  
                exit  
            else if (this%elements(i)%full_index > full_index) then
                new_sparse_element_needed = .TRUE.
                shift_index = i
                exit

            ! if given full index is neither equal to, nor less than current full, and if we are on the
            ! last sparse element, then add an element 
            else if (i == this%sparse_size) then
                new_sparse_element_needed = .TRUE.
                shift_index = i + 1
                exit
            end if
        end do

        if (new_sparse_element_needed) then
            call this%add_element(value, full_index, shift_index)
        end if
       
    end subroutine sparse_vector_set_value


    subroutine sparse_vector_compress(this)
        ! compresses 0's in a sparse vector type 
        
        implicit none 
        
        class(sparse_vector),intent(inout) :: this
        
        type(sparse_vector) :: temp_vector
        integer :: i,count
        integer, dimension(this%sparse_size) :: indices
        
        ! count how many nonzero elements there are
        count = 0
        do i=1, this%sparse_size
            if (abs(this%elements(i)%value) > 1.0e-16) then
                count = count + 1
                indices(count) = i
            end if
        end do

        ! check to see if its worth it to compress. 0.8 is arbitrary, adjust as necessary
        ! 0.8 means only 80 percent of the sparse elements are actually nonzero. If count/sparse_size
        ! is less than 0.8, it is worth it to compress, if its something like 0.95, don't bother
        if (real(count)/real(this%sparse_size) < 0.8) then
             
            ! now that we have the number of non zero elements, allocate temp array
            temp_vector%sparse_size = count
            allocate(temp_vector%elements(count))

            ! store the nonzero sparse elemnets in a temp array
            do i=1,count
                temp_vector%elements(i) = this%elements(indices(i))
            end do

            ! move allocation of temporary array and overwrite 'this'
            call move_alloc(temp_vector%elements, this%elements) 
            
            ! update the sparse number of columns
            this%sparse_size = count
        else

        end if

    end subroutine sparse_vector_compress

    
    function sparse_vector_expand(this) result(full_vector)
        ! takes a sparse vector and expands it, returning a real full vector

        implicit none

        class(sparse_vector),intent(inout) :: this

        real,dimension(:),allocatable :: full_vector
        integer :: i

        allocate(full_vector(this%full_size), source=0.)
        
        ! put nonzero values in their corresponding full_index location
        do i=1,this%sparse_size
            full_vector(this%elements(i)%full_index) = this%elements(i)%value
        end do

    end function sparse_vector_expand


    subroutine sparse_vector_sparse_add(this, sparse_input) 
        ! subroutine to add a sparse vector to this
        ! this + sparse_input = this

        implicit none

        class(sparse_vector),intent(inout) :: this
        type(sparse_vector) :: sparse_input

        integer :: i
        real :: this_i, sparse_input_i, added
        
        ! make sure the input vector has the same full size as this
        if (this%full_size /= sparse_input%full_size) then
            write(*,*) "Error!!! sparse_vector_sparse_add requires input to have the same full_size. Quitting..."
            stop
        end if

        ! loop through full index
        do i=1, this%full_size
    
            ! get vector values at full index i
            sparse_input_i = sparse_input%get_value(i)
        
            ! if sparse_input_i is populated, add them
            if (abs(sparse_input_i) > 1.0e-16) then
                                
                this_i = this%get_value(i)
                added = this_i + sparse_input_i
                call this%set_value(added, i)

            end if
            
        end do 
        
    end subroutine sparse_vector_sparse_add


    subroutine sparse_vector_sparse_subtract(this, sparse_input)
        ! subroutine to subtract a sparse vector from this
        ! this - sparse_input = this

        implicit none

        class(sparse_vector),intent(inout) :: this
        type(sparse_vector) :: sparse_input

        integer :: i
        real :: this_i, sparse_input_i, subtracted

        ! make sure the input vector has the same full size as this
        if (this%full_size /= sparse_input%full_size) then
            write(*,*) "Error!!! sparse_vector_sparse_subtract requires input to have the same full_size. Quitting..."
            stop
        end if
        
       
        ! loop through full index
        do i=1, this%full_size
    
            ! get vector values at full index i
            sparse_input_i = sparse_input%get_value(i)
        
            ! if sparse_input_i is populated, subtract them
            if (abs(sparse_input_i) > 1.0e-16) then
                                
                this_i = this%get_value(i)
                subtracted = this_i - sparse_input_i
                call this%set_value(subtracted, i)
                
            end if
            
        end do 

    end subroutine sparse_vector_sparse_subtract 


    function sparse_vector_broadcast_element_times_vector(this, vec) result(result_matrix)
        ! multiply a vector by each scalar in the sparse_vector
        ! returns a sparse_matrix

        implicit none

        class(sparse_vector),intent(inout) :: this
        real,dimension(3), intent(in) :: vec

        type(sparse_matrix) :: result_matrix
        integer :: i

        ! allocate sparse_matrix elements
        allocate(result_matrix%columns(this%sparse_size))

        ! perform scalar multiplication of vec by each sparse vector scalar
        do i=1,this%sparse_size
            result_matrix%columns(i)%vector_values = vec*this%elements(i)%value
            result_matrix%columns(i)%full_index = this%elements(i)%full_index
        end do 

        ! give sparse_matrix sizes
        result_matrix%sparse_num_cols = this%sparse_size
        result_matrix%full_num_cols = this%full_size
        
        
    end function sparse_vector_broadcast_element_times_vector 


    subroutine sparse_vector_broadcast_element_times_scalar(this, scalar) 
        ! multiply a scalar by each scalar in the sparse_vector

        implicit none

        class(sparse_vector),intent(inout) :: this
        real,intent(in) :: scalar

        integer :: i


        ! perform scalar multiplication of vec by each sparse vector scalar
        do i=1,this%sparse_size
            this%elements(i)%value = scalar*this%elements(i)%value
        end do 
        
    end subroutine sparse_vector_broadcast_element_times_scalar


!!!!!!!!!!!!!!!!!!! END SPARSE VECTOR TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

!!!!!!!!!!!!!!!!!!! SPARSE MATRIX TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine sparse_matrix_init(this,N)
        ! initializes a sparse matrix with one element full of zeros. if a derivative WRT X(beta), N = N_verts*3
        
        implicit none

        class(sparse_matrix), intent(inout) :: this
        integer,intent(in) :: N

        real,dimension(3) :: zeros

        ! specify number of design variables N (if WRT to X(beta), N = N_verts*3)
        this%full_num_cols = N

        ! allocate and populate an initial sparse matrix element
        allocate(this%columns(1))
        zeros = (/0.0, 0.0, 0.0/)
        this%columns(1)%vector_values = zeros
        this%columns(1)%full_index = 1
        this%sparse_num_cols = 1

    end subroutine sparse_matrix_init


    subroutine sparse_matrix_init_from_sparse_vectors(this, sparse_v1, sparse_v2, sparse_v3)
        ! combines 3 sparse vectors into a sparse matrix
        implicit none

        class(sparse_matrix),intent(inout) :: this
        type(sparse_vector),intent(inout) :: sparse_v1,sparse_v2,sparse_v3
        
        integer :: i, count
        real,dimension(3) :: values

        ! check to see if the given sparse vectors are the same full size
        if (sparse_v1%full_size == sparse_v2%full_size .and. sparse_v1%full_size == sparse_v3%full_size) then
            
            this%full_num_cols = sparse_v1%full_size
        else
            write(*,*) "!!! sparse_matrix_init_from_sparse_vectors requires all sparse vector inputs &
            to have same full_size. Quitting..."
            stop
        end if
        count = 0
        do i=1,this%full_num_cols
            ! check to see if the given sparse vectors have a nonzero value at full index i 
            ! if at least one sparse vector has a nonzero value at i, allocate a and populate 3 vector
            if (abs(sparse_v1%get_value(i)) > 1.0e-16 .or. abs(sparse_v2%get_value(i)) > 1.0e-16 &
            .or. abs(sparse_v3%get_value(i)) > 1.0e-16) then           
                count = count + 1
                
                ! update array of values (at least 1 will be nonzero)
                values = (/sparse_v1%get_value(i), sparse_v2%get_value(i), sparse_v3%get_value(i) /)
            
                ! the first element needs to be initialized 
                if (count == 1) then
                    ! allocate the first sparse matrix element
                    allocate(this%columns(1))

                    ! populate the first sparse matrix element
                    this%columns(1)%vector_values = values
                    this%columns(1)%full_index = i
                    this%sparse_num_cols = 1

                else
                    ! add an element (a column)
                    call this%add_element(values,i,count)

                end if

            end if
        end do
    
    end subroutine sparse_matrix_init_from_sparse_vectors


    subroutine sparse_matrix_init_from_full_matrix(this, full_matrix) 
        ! takes a full matrix that is sparse (has a lot of elements with 0.0s) and converts it to a sparse matrix

        implicit none 

        class(sparse_matrix),intent(inout) :: this
        real,dimension(:,:),intent(inout) :: full_matrix
        integer :: i, full_size,count
        integer, dimension(:),allocatable :: indices

        ! get size of vector
        full_size = size(full_matrix,2)
        allocate(indices(full_size))

        ! count how many nonzero elements there are
        count = 0
        do i=1,full_size
            if (any(abs(full_matrix(:,i)) > 1.0e-16)) then
                count = count + 1
                indices(count) = i
            end if
        end do
        
        ! now that we have the number of non zero numbers, we can allocate the space for the sparse_vector
        this%sparse_num_cols = count
        this%full_num_cols = full_size
        allocate(this%columns(count))

        ! populate the sparse_vector elements
        do i=1,count
            ! populate sparse columns with the nonzero full matrix columns
            this%columns(i)%vector_values = full_matrix(:,indices(i))
            this%columns(i)%full_index = indices(i)
        end do
        
        ! to save memory, you can deallocate the original array after compressing it

    end subroutine sparse_matrix_init_from_full_matrix


    subroutine sparse_matrix_init_from_sparse_matrix(this, sparse_input)
        ! inits a sparse matrix as a copy of the input sparse matrix

        implicit none

        class(sparse_matrix),intent(inout) :: this
        type(sparse_matrix) :: sparse_input

        integer :: i

        ! copy info over
        this%sparse_num_cols = sparse_input%sparse_num_cols
        this%full_num_cols = sparse_input%full_num_cols
        
        ! allocate the same number of sparse matrix elements
        allocate(this%columns(sparse_input%sparse_num_cols))
        
        do i=1,this%sparse_num_cols
            
            ! copy each element
            this%columns(i) = sparse_input%columns(i)

        end do    

    end subroutine sparse_matrix_init_from_sparse_matrix

    
    subroutine sparse_matrix_increase_size(this)
        ! increases the size of the sparse_matrix by 1 sparse element

        implicit none
        class(sparse_matrix),intent(inout) :: this
        
        type(sparse_matrix) :: temp_matrix
        integer :: new_size
        
        new_size = size(this%columns) + 1

        ! allocate a temporary array
        allocate(temp_matrix%columns(new_size))
        
        ! copy array info
        temp_matrix%columns(1:new_size-1) = this%columns(:)
        
        ! move allocation of temporary array and overwrite this
        call move_alloc(temp_matrix%columns, this%columns)

        ! update the sparse number of columns
        this%sparse_num_cols = this%sparse_num_cols + 1
        

    end subroutine sparse_matrix_increase_size


    subroutine sparse_matrix_add_element(this, values,  full_index, shift_index)
        ! adds a sparse element to the current sparse matrix
        ! could be made more efficient by checking if the shift index is closer to the beginning or 
        ! the end. if closer to the beginning, cshift, then pull back the first few elements until shift_index
        
        implicit none

        class(sparse_matrix),intent(inout) :: this
        integer, intent(in) :: full_index, shift_index
        real,dimension(:),intent(in) :: values
        
        integer :: i
        
        ! original information is preserved, but an extra space is allocated
        call this%increase_size()
        
        ! shift necessary elements up in the increased size matrix
        do i=this%sparse_num_cols,shift_index,-1
            this%columns(i)%vector_values(:) = this%columns(i-1)%vector_values(:)
            this%columns(i)%full_index = this%columns(i-1)%full_index
        end do

        ! place the new element in the correct spot
        this%columns(shift_index)%vector_values(:) = values(:)
        this%columns(shift_index)%full_index = full_index

    end subroutine sparse_matrix_add_element


    function sparse_matrix_get_values(this, full_index) result(values)
        ! given a full matrix column index, this function returns the 3 values at that index

        implicit none

        class(sparse_matrix),intent(inout) :: this
        integer :: full_index, i

        real,dimension(3) :: values 

        values = (/1.0e-16, 1.0e-16, 1.0e-16/)

        ! if the sparse matrix element has the same full index as the given full index, 
        ! return that element's value, if not, return the default zeros        
        do i=1,this%sparse_num_cols
            if (this%columns(i)%full_index == full_index) then
                values(:) = this%columns(i)%vector_values(:)
                exit                
            end if  
        end do
        

    end function sparse_matrix_get_values
  

    subroutine sparse_matrix_set_values(this, values, full_index) 
        ! given a full_vector index and a value, this adds or updates the corresponding value in the sparse vector

        implicit none

        class(sparse_matrix),intent(inout) :: this
        integer, intent(in):: full_index
        real,dimension(3),intent(in) :: values
        
        integer :: i, shift_index
        logical :: new_sparse_element_needed

        ! initialize logical
        new_sparse_element_needed = .FALSE.

        do i=1,this%sparse_num_cols
                
            ! check to see the given full_index matches an existing sparse element full index
            if (this%columns(i)%full_index == full_index) then
                this%columns(i)%vector_values(:) = values(:)  
                exit  
                
            ! check to see if the given full index is less than the current element's full index (we skipped it)
            else if (this%columns(i)%full_index > full_index) then
                new_sparse_element_needed = .TRUE.
                shift_index = i
                exit
            
            ! if given full index is neither equal to, nor less than current full, and if we are on the
            ! last sparse column, then add an element 
            else if (i == this%sparse_num_cols) then
                new_sparse_element_needed = .TRUE.
                shift_index = i + 1
                exit

            end if
            
        end do

        if (new_sparse_element_needed) then
            call this%add_element(values, full_index, shift_index)
        end if
       
    end subroutine sparse_matrix_set_values


    subroutine sparse_matrix_compress(this) 
        ! compresses all zero elements in a sparse matrix type 
        
        implicit none 
        
        class(sparse_matrix),intent(inout) :: this
        
        type(sparse_matrix) :: temp_matrix
        integer :: i,count
        integer, dimension(:), allocatable :: indices
        
        allocate(indices(this%sparse_num_cols))

        ! count and store indices of NONZERO elements
        count = 0
        do i=1,this%sparse_num_cols
            if (any(abs(this%columns(i)%vector_values) > 1.0e-16)) then
                count = count + 1
                indices(count) = i
            end if
        end do
       

        ! check to see if its worth it to compress. 0.8 is arbitrary, adjust as necessary
        ! 0.8 means only 80 percent of the sparse elements are actually nonzero. If count/sparse_num_cols
        ! is less than 0.8, it is worth it to compress, if its something like 0.95, don't bother
        if (real(count)/real(this%sparse_num_cols) < 0.8) then
             
            ! now that we have the number of NON ZERO elements, allocate temp sparse matrix
            temp_matrix%sparse_num_cols = count
            allocate(temp_matrix%columns(count))

            ! store the nonzero sparse elemnets in a temp array
            do i=1,count
                temp_matrix%columns(i) = this%columns(indices(i))
            end do

            ! move allocation of temporary array and overwrite 'this' (does it overwrite?)
            call move_alloc(temp_matrix%columns, this%columns) 
            
            ! update the sparse number of columns
            this%sparse_num_cols = count
        else

        end if
    end subroutine sparse_matrix_compress


    function sparse_matrix_expand(this, tall) result(full_matrix)
        ! takes a sparse matrix and expands it, returning a real full matrix

        implicit none

        class(sparse_matrix),intent(inout) :: this
        logical,intent(in) :: tall

        real,dimension(:,:),allocatable :: full_matrix
        integer :: i

        allocate(full_matrix(3,this%full_num_cols), source=0.)
        
        ! put nonzero values in their corresponding full_index location (column major)
        do i=1,this%sparse_num_cols
            full_matrix(:,this%columns(i)%full_index) = this%columns(i)%vector_values
        end do

        ! if tall is the desired output, make N rows by 3 column (used in matrix split to vectors)
        if (tall) then
            full_matrix = transpose(full_matrix)
        end if 
    end function sparse_matrix_expand


    subroutine sparse_matrix_sparse_add(this, sparse_input) 
        ! subroutine to add a sparse matrix to this
        ! this + sparse_input = this

        implicit none

        class(sparse_matrix),intent(inout) :: this
        type(sparse_matrix) :: sparse_input

        integer :: i
        real,dimension(3) :: this_i, sparse_input_i, added

        ! make sure the input matrix has the same full size as this
        if (this%full_num_cols /= sparse_input%full_num_cols) then
            write(*,*) "Error!!! sparse_matrix_add requires input to have the same full_size. Quitting..."
            stop
        end if

       
        ! loop through full index
        do i=1, this%full_num_cols
    
            ! get vector values at full index i
            sparse_input_i = sparse_input%get_values(i)
        
            ! if sparse_input_i is populated, add them
            if (any(abs(sparse_input_i) > 1.0e-16)) then
                                
                this_i = this%get_values(i)
                added = this_i + sparse_input_i
                call this%set_values(added, i)

            end if
            
        end do 
        
    end subroutine sparse_matrix_sparse_add



    subroutine sparse_matrix_sparse_subtract(this, sparse_input)
        ! subroutine to subtract a sparse matrix from this
        ! this - sparse_input = this

        implicit none

        class(sparse_matrix),intent(inout) :: this
        type(sparse_matrix) :: sparse_input

        integer :: i
        real,dimension(3) :: this_i, sparse_input_i, subtracted

        ! make sure the input matrix has the same full size as this
        if (this%full_num_cols /= sparse_input%full_num_cols) then
            write(*,*) "Error!!! sparse_matrix_subtract requires input to have the same full_size. Quitting..."
            stop
        end if
        
       
        ! loop through full index
        do i=1, this%full_num_cols
    
            ! get vector values at full index i
            sparse_input_i = sparse_input%get_values(i)
        
            ! if sparse_input_i is populated, subtract them
            if (any(abs(sparse_input_i) > 1.0e-16)) then
                                
                this_i = this%get_values(i)
                subtracted = this_i - sparse_input_i
                call this%set_values(subtracted, i)
                
            end if
            
        end do 

    end subroutine sparse_matrix_sparse_subtract


    
    
    function sparse_matrix_broadcast_vector_cross_element(this, vec) result(result_matrix)
        ! cross product of a 3-vector and a "list" of vectors (derivative of a vector)
        ! returns a "list" of 3-vectors ( a sparse matrix 3)
        
        implicit none
        
        class(sparse_matrix),intent(inout) :: this
        real,dimension(3),intent(in) :: vec 
        
        type(sparse_matrix) :: result_matrix
        real,dimension(3) :: temp_vec
        
        integer :: i

        ! set full and sparse column numbers/size
        result_matrix%full_num_cols = this%full_num_cols
        result_matrix%sparse_num_cols = this%sparse_num_cols
        
        ! allocate result matrix elements
        allocate(result_matrix%columns(size(this%columns)))
        
        ! go through each sparse column index
        do i=1,this%sparse_num_cols
            ! do cross product
            temp_vec = this%columns(i)%vector_values(:)
            temp_vec = cross(vec,temp_vec)
            
            ! update the results matrix
            result_matrix%columns(i)%vector_values(:) = temp_vec(:)
            result_matrix%columns(i)%full_index = this%columns(i)%full_index
            
            
        end do       
        
        
    end function sparse_matrix_broadcast_vector_cross_element
    
    
    function sparse_matrix_broadcast_element_cross_vector(this, vec) result(result_matrix)
        ! cross product of a "list" of vectors (derivative of a vector) and a 3-vector
        ! returns a "list" of 3-vectors ( a sparse matrix 3)
        
        implicit none
        
        class(sparse_matrix),intent(inout) :: this
        real,dimension(3),intent(in) :: vec 
        
        type(sparse_matrix) :: result_matrix
        real,dimension(3) :: temp_vec
        
        integer :: i
        
        ! set full and sparse column numbers/size
        result_matrix%full_num_cols = this%full_num_cols
        result_matrix%sparse_num_cols = this%sparse_num_cols
        
        ! allocate result matrix elements
        allocate(result_matrix%columns(size(this%columns)))
        
        ! go through each sparse column index
        do i=1,this%sparse_num_cols
            ! do cross product
            temp_vec = this%columns(i)%vector_values(:)
            temp_vec = cross(temp_vec, vec)
            
            ! update the results matrix 
            result_matrix%columns(i)%vector_values(:) = temp_vec(:)
            result_matrix%columns(i)%full_index = this%columns(i)%full_index
            
        end do       
        
    end function sparse_matrix_broadcast_element_cross_vector
    
    
    function sparse_matrix_broadcast_vector_dot_element(this, vec) result(result_vector)
        ! dot product of a 3-vector and a "list" of vectors (derivative of a vector)
        ! returns a sparse_vector
        
        implicit none
        
        class(sparse_matrix),intent(inout) :: this
        real,dimension(3),intent(in) :: vec 
        
        type(sparse_vector) :: result_vector
        real :: temp_val
        real,dimension(3) :: temp_vec
        
        integer :: i
        
        ! set full and sparse column numbers/size
        result_vector%full_size = this%full_num_cols
        result_vector%sparse_size = this%sparse_num_cols
        ! allocate result matrix elements
        allocate(result_vector%elements(result_vector%sparse_size))
        
        ! go through each sparse column index
        do i=1,this%sparse_num_cols
            ! do dot product
            temp_vec = this%columns(i)%vector_values(:)
            temp_val = inner(vec,temp_vec)
            
            ! update results vector
            result_vector%elements(i)%value = temp_val
            result_vector%elements(i)%full_index = this%columns(i)%full_index
            
        end do       
        
    end function sparse_matrix_broadcast_vector_dot_element
    
    
    subroutine sparse_matrix_broadcast_element_times_scalar(this,scalar)
        ! multiplies a single scalar value by each of the sparse matrix elements
        ! (this could also be a function)
        
        implicit none
        class(sparse_matrix),intent(inout) :: this
        real,intent(in) :: scalar
        
        integer :: i
        
        
        ! mutiply the scalar
        do i=1,this%sparse_num_cols
            this%columns(i)%vector_values = this%columns(i)%vector_values * scalar
        end do
        
    end subroutine sparse_matrix_broadcast_element_times_scalar
    

    function sparse_matrix_broadcast_matmul_3x3_times_element(this, matrix3) result(result_matrix)
        ! matmul a 3x3 matrix with a "list" of vectors (sparse matrix)
        ! returns a "list" of 3-vectors ( a sparse matrix )
        
        implicit none
        
        class(sparse_matrix),intent(inout) :: this
        real,dimension(3,3),intent(in) :: matrix3 
        
        type(sparse_matrix) :: result_matrix
        real,dimension(3) :: temp_vec
        
        integer :: i
        
        ! set full and sparse column numbers/size
        result_matrix%full_num_cols = this%full_num_cols
        result_matrix%sparse_num_cols = this%sparse_num_cols
        
        ! allocate result matrix elements
        allocate(result_matrix%columns(size(this%columns)))
        
        ! go through each sparse column index
        do i=1,this%sparse_num_cols
            ! do matmul
            temp_vec = this%columns(i)%vector_values(:)
            temp_vec = matmul(matrix3,  temp_vec)
            
            ! update the results matrix 
            result_matrix%columns(i)%vector_values(:) = temp_vec(:)
            result_matrix%columns(i)%full_index = this%columns(i)%full_index
            
        end do       
        
    end function sparse_matrix_broadcast_matmul_3x3_times_element
    
    
    function sparse_matrix_broadcast_matmul_element_times_3x3(this, matrix3) result(result_matrix)
        ! matmul a  "list" of vectors (sparse matrix) with a 3x3 matrix
        ! returns a "list" of 3-vectors ( a sparse matrix )
        
        implicit none
        
        class(sparse_matrix),intent(inout) :: this
        real,dimension(3,3),intent(in) :: matrix3 
        
        type(sparse_matrix) :: result_matrix
        real,dimension(3) :: temp_vec
        
        integer :: i
        
        ! set full and sparse column numbers/size
        result_matrix%full_num_cols = this%full_num_cols
        result_matrix%sparse_num_cols = this%sparse_num_cols
        
        ! allocate result matrix elements
        allocate(result_matrix%columns(size(this%columns)))
        
        ! go through each sparse column index
        do i=1,this%sparse_num_cols
            ! do matmul
            temp_vec = this%columns(i)%vector_values(:)
            ! temp_vec = matmul(matrix3,  temp_vec)
            temp_vec = matmul(temp_vec,matrix3)
            
            ! update the results matrix 
            result_matrix%columns(i)%vector_values(:) = temp_vec(:)
            result_matrix%columns(i)%full_index = this%columns(i)%full_index
            
        end do       
        
    end function sparse_matrix_broadcast_matmul_element_times_3x3
    
    
    
    
    function sparse_matrix_split_into_sparse_vectors(this) result(sparse_vecs)
        ! takes a sparse matrix and converts it to a sparse_vector dimension 3
        
        implicit none
        
        class(sparse_matrix),intent(inout) :: this
        
        integer :: i
        real,dimension(:,:), allocatable :: full_matrix
        type(sparse_vector),dimension(3) :: sparse_vecs
        logical :: tall 
        
        ! tall means we want the matrix expanded as a Nx3 rather than 3xN
        tall = .true.
        full_matrix = this%expand(tall)
        
        do i=1,3
            call sparse_vecs(i)%init_from_full_vector(full_matrix(:,i))
        end do
        
    end function sparse_matrix_split_into_sparse_vectors
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!! START SPARSE 3D TYPE BOUND PROCEDURES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sparse_3D_init_from_sparse_matrices(this, sparse_matrices)
        ! takes a sparse_matrix dimension(x) and makes a sparse 3D object

        implicit none

        class(sparse_3D), intent(inout) :: this
        type(sparse_matrix), dimension(:), intent(in) :: sparse_matrices

        integer :: i

        this%num_rows = size(sparse_matrices)
        allocate(this%rows(this%num_rows))

        do i=1,this%num_rows
            call this%rows(i)%init_from_sparse_matrix(sparse_matrices(i))
        end do

    end subroutine sparse_3D_init_from_sparse_matrices

    subroutine sparse_3D_sparse_add_3(this, sparse_3D_3rows)
        ! adds a sparse_matrix dimension 3
        
        implicit none
        
        class(sparse_3D),intent(inout) :: this
        type(sparse_3D), intent(in) :: sparse_3D_3rows 
        
        integer :: i
        
        do i=1,3
            call this%rows(i)%sparse_add(sparse_3D_3rows%rows(i))
        end do
        
    end subroutine sparse_3D_sparse_add_3
    
    
    function sparse_3D_transpose_3(this) result(transposed)
        ! takes a sparse matrix and converts it to a sparse_vector dimension 3
        
        implicit none
        
        class(sparse_3D),intent(inout) :: this
        
        integer :: i
        type(sparse_vector),dimension(3) :: row1, row2, row3
        type(sparse_3D) :: transposed
        
        allocate(transposed%rows(3))
        
        row1 = this%rows(1)%split_into_sparse_vectors()
        row2 = this%rows(2)%split_into_sparse_vectors()
        row3 = this%rows(3)%split_into_sparse_vectors()
        
        do i=1,3
            call transposed%rows(i)%init_from_sparse_vectors(row1(i), row2(i), row3(i))
        end do
        
    end function sparse_3D_transpose_3

    
    function sparse_3D_broadcast_matmul_3row_times_3x3(this, matrix3) result(result_rows)
        ! matmuls a sparse_matrix dim(3) (rows) by a 3x3 matrix
    
        implicit none
    
        class(sparse_3D),intent(inout) :: this
        real,dimension(3,3),intent(in) :: matrix3 
    
        integer :: i
        type(sparse_3D):: result_rows

        allocate(result_rows%rows(3))
    
        do i=1,3
            result_rows%rows(i) = this%rows(i)%broadcast_matmul_element_times_3x3(matrix3)
        end do
    
    end function sparse_3D_broadcast_matmul_3row_times_3x3
    
    
    function sparse_3D_broadcast_matmul_3x3_times_3row(this, matrix3) result(result_rows)
        ! matmuls a 3x3 matrix times a sparse_matrix dim(3) (rows)
    
        implicit none
    
        class(sparse_3D),intent(inout) :: this
        real,dimension(3,3),intent(in) :: matrix3 
    
        integer :: i
        type(sparse_3D) :: this_cols, temp_cols, result_rows

        allocate(this_cols%rows(3))
        allocate(temp_cols%rows(3))
        allocate(result_rows%rows(3))

        this_cols = this%transpose_3()
    
        do i=1,3                           
            temp_cols%rows(i) = this_cols%rows(i)%broadcast_matmul_3x3_times_element(matrix3)
            ! cols%rows above sounds confusing, but it means this row is used like a column in the matrix multiplication
        end do

        result_rows = temp_cols%transpose_3()
    
    end function sparse_3D_broadcast_matmul_3x3_times_3row


    function sparse_3D_broadcast_matmul_3row_times_3x1(this, vector3) result(sparse_mat)
        ! matmuls a sparse_matrix dim(3) (rows) by a 3x1 matrix
    
        implicit none
    
        class(sparse_3D),intent(inout) :: this
        real,dimension(3,3),intent(in) :: vector3 
    
        integer :: i
        type(sparse_vector), dimension(3) :: temp
        type(sparse_matrix):: sparse_mat
    
        do i=1,3
            temp(i) = this%rows(i)%broadcast_vector_dot_element(vector3)
        end do

        call sparse_mat%init_from_sparse_vectors(temp(1), temp(2), temp(3))
    
    end function sparse_3D_broadcast_matmul_3row_times_3x1
    
    
    function sparse_3D_broadcast_matmul_1x3_times_3row(this, vector3) result(sparse_mat)
        ! matmuls a 1x3 matrix times a sparse_matrix dim(3) (rows)
    
        implicit none
    
        class(sparse_3D),intent(inout) :: this
        real,dimension(3,3),intent(in) :: vector3 

        integer :: i
        type(sparse_3D) :: this_cols
        type(sparse_vector), dimension(3) :: temp
        type(sparse_matrix):: sparse_mat

        allocate(this_cols%rows(3))

        this_cols = this%transpose_3()
    
        do i=1,3                           
            tmep(i) = this_cols%rows(i)%broadcast_vector_dot_element(vector3)
            ! cols%rows above sounds confusing, but it means this row is used like a column in the matrix multiplication
        end do

        call sparse_mat%init_from_sparse_vectors(temp(1), temp(2), temp(3)) 
    
    end function sparse_3D_broadcast_matmul_1x3_times_3row


    function sparse_3D_convert_to_sparse_vector_3x3(this) result(sparse_vec_3x3)
        ! converts a sparse 3D (3row) to a sparse vector 3x3

        implicit none

        class(sparse_3D),intent(inout) :: this
    
        integer :: i
        type(sparse_vector),dimension(3,3) :: sparse_vec_3x3

        do i=1,3
            sparse_vec_3x3(i,:) = this%rows(i)%split_into_sparse_vectors()
        end do

    end function sparse_3D_convert_to_sparse_vector_3x3
    
    
    
    
    ! subroutine adjoint_init(this)
    !     ! initializes adjoint type
    
    !     implicit none
    
    !     class(adjoint),intent(inout) :: this
    !     !type(surface_mesh),intent(in) :: body
    
    !     integer :: i, j
    
    !     !this%size = body%N_verts*3
    
    !     !call this%get_X_beta(body)
    
    ! end subroutine adjoint_init
    
    
    ! subroutine adjoint_get_X_beta(this)
    !     ! builds a vector of design variables
    
    !     implicit none
    
    !     class(adjoint),intent(inout) :: this
    !     !type(surface_mesh),intent(in) :: body
    
    !     integer :: i, j
    
    !     ! build design variable vector X_beta
    !     !allocate(this%X_beta(this%size))

        
    !     ! ! place all x values in X_beta, followed by all y, then z values 
    !     ! ! X_beta will look like this: x1, x2, x3, ...xn, y1, y2, y3, ...yn, z1, z2, z3, ...zn
    !     ! do i=1,3 
    !     !     do j=1,body%N_verts
    !     !         this%X_beta(j + (i-1)*body%N_verts) = body%vertices(j)%loc(i)
    !     !     end do
    !     ! end do
        
        

    ! end subroutine adjoint_get_X_beta



end module adjoint_mod