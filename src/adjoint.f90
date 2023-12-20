! module for calculating the gradient of CF, CM with respect to mesh points via the adjoint method
module adjoint_mod

    !use linked_list_mod
    use math_mod
    use helpers_mod
    use base_geom_mod
    use mesh_mod
    use surface_mesh_mod

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
            procedure :: compress => sparse_vector_compress
            procedure :: expand => sparse_vector_expand
            procedure :: increase_size => sparse_vector_increase_size
            procedure :: add_element => sparse_vector_add_element
            procedure :: get_value => sparse_vector_get_value
            procedure :: set_value => sparse_vector_set_value

            procedure :: multiply_single_vector => sparse_vector_multiply_single_vector
            
            !procedure :: sparse_add => sparse_vector_sparse_add
            !procedure :: sparse_subtract => sparse_vector_sparse_subtract

            !procedure :: fill_vector => sparse_vector_fill_vector
            
            
    end type sparse_vector
            
            
    type sparse_matrix
        ! used to save memory when we have matrices made of sparse vectors
        
        integer :: sparse_num_cols,full_num_cols    ! sparse_num_cols is the number of columns with a nonzero element
        type(sparse_matrix_element),dimension(:),allocatable :: columns       ! non zero values in sparse vector
        
            
        contains
            procedure :: init => sparse_matrix_init
            procedure :: count_nonzero_matrix_elements => sparse_matrix_count_nonzero_matrix_elements
            procedure :: compress => sparse_matrix_compress
            !procedure :: expand => sparse_matrix_expand
            procedure :: increase_size => sparse_matrix_increase_size
            procedure :: add_element => sparse_matrix_add_element
            procedure :: get_values => sparse_matrix_get_values
            procedure :: set_values => sparse_matrix_set_values
            
            procedure :: vec_cross_matrix => sparse_matrix_vec_cross_matrix
            procedure :: matrix_cross_vec => sparse_matrix_matrix_cross_vec

            procedure :: dot_vec_and_matrix => sparse_matrix_dot_vec_and_matrix

            procedure :: multiply_scalar => sparse_matrix_multiply_scalar

            procedure :: sparse_add => sparse_matrix_sparse_add
            procedure :: sparse_subtract => sparse_matrix_sparse_subtract

        

    end type sparse_matrix


    type adjoint

        integer :: size    ! number of design varibles, the length of X_beta
        ! we may never use X_beta
        real,dimension(:),allocatable :: X_beta  ! list of x y z values of all mesh points (design variables)
        type(sparse_vector),dimension(:),allocatable :: d_xyz
    

        contains
            procedure :: adjoint => adjoint_init
            procedure :: get_X_beta => adjoint_get_X_beta ! puts all x y z values of mesh points in a list

    end type adjoint


contains


    subroutine sparse_vector_init(this, full_size)
        ! initializes a sparse_vector type, to be used when creating a sparse vector of known length and values

        implicit none

        class(sparse_vector),intent(inout) :: this
        integer :: full_size

        this%full_size = full_size

    end subroutine sparse_vector_init
        

    subroutine sparse_vector_compress(this, full_vector, full_size)
        ! takes a full vector that is sparse (has a lot of 0.0's) and compresses it to a sparse vector type

        implicit none 

        class(sparse_vector),intent(inout) :: this
        integer,intent(in) :: full_size
        real,dimension(:),intent(in) :: full_vector
        integer :: i,count

        
        ! count how many nonzero elements there are
        count = count_nonzero_vector_values(full_vector)
        
        ! now that we have the number of non zero numbers, we can allocate the space for the sparse_vector
        this%sparse_size = count
        this%full_size = full_size
        allocate(this%elements(count))

        ! populate the sparse_vector elements
        do i=1,full_size
            if (full_vector(i) /= 0.0) then
                this%elements(i)%value = full_vector(i)
                this%elements(i)%full_index = i
            end if
        end do
        
        ! to save memory, you can deallocate the original array after compressing it

    end subroutine sparse_vector_compress

    function sparse_vector_expand(this) result(full_vector)
        ! takes a full vector that is sparse (has a lot of 0.0's) and compresses it to a sparse vector type

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


    subroutine sparse_vector_increase_size(this)
        ! increases the size of the sparse_matrix by 1 sparse element

        implicit none
        class(sparse_vector),intent(inout) :: this

        type(sparse_vector),dimension(:),allocatable :: temp_vector
        integer :: new_size

        new_size = size(this) + 1
        ! allocate a temporary array
        allocate(temp_vector(new_size))

        ! copy array info
        temp_vector%elements(1:new_size-1) = this%elements(:)

        ! move allocation of temporary array and overwrite this
        call move_alloc(temp_vector, this) 
        
        ! update the sparse number of columns
        this%sparse_size = this%sparse_size + 1

    end subroutine sparse_vector_increase_size


    subroutine sparse_vector_add_element(this, full_index, value, shift_index)
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
            this%elements(i+1)%value = this%elements(i)%value
            this%elements(i+1)%full_index = this%elements(i)%full_index
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

        ! if the sparse vector element has the same full index as the given full index value is set to
        ! that element's value, if not, the value stays zero 
        do i=1,this%sparse_size
            if (this%elements(i)%full_index == full_index) then
                value = this%elements(i)%value
                exit
            else if (this%elements(i)%full_index > full_index) then
                value = 0.0
                exit
            end if  
        end do
       
    end function sparse_vector_get_value


    subroutine sparse_vector_set_value(this, full_index, value) 
        ! given a full_vector index and a value, this adds or updates the corresponding value in the sparse vector

        implicit none

        class(sparse_vector),intent(inout) :: this
        integer, intent(in):: full_index
        real,intent(inout) :: value
        
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
            end if
        end do

        if (new_sparse_element_needed) then
            call this%add_element(full_index,value, shift_index)
        end if
       
    end subroutine sparse_vector_set_value


    function sparse_vector_multiply_single_vector(this, vec) result(result_matrix)
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
        
        
    end function sparse_vector_multiply_single_vector


    subroutine sparse_matrix_init(this, sparse_v1, sparse_v2, sparse_v3)
        ! combines 3 sparse vectors into a sparse matrix
        implicit none

        class(sparse_matrix),intent(inout) :: this
        type(sparse_vector),intent(in) :: sparse_v1,sparse_v2,sparse_v3
        
        integer :: i, column_count = 0

        ! check to see if the given sparse vectors are the same full size
        if (sparse_v1%full_size == sparse_v2%full_size .and. sparse_v1%full_size == sparse_v3%full_size) then
            
            this%full_num_cols = sparse_v1%full_size
        else
            write(*,*) "!!! sparse_matrix_init requires all sparse vector inputs to have same full_size. Quitting..."
            stop
        end if
        
        do i=1,this%full_num_cols
            ! check to see if the given sparse vectors have a nonzero value at full index i 
            ! if at least one sparse vector has a nonzero value at i, allocate a and populate 3 vector
            if (sparse_v1%elements(i)%full_index == i .or. sparse_v2%elements(i)%full_index == i &
                .or. sparse_v3%elements(i)%full_index) then           
                column_count = column_count+1
                
                ! allocate space for one more element
                call this%increase_size()

                ! associate this element with the full index i
                this%columns(i)%full_index = i

                ! populate the vector_values  NOTE: this could be a bit more memory efficient if the 0.0 vector_values were collapse
                this%columns(i)%vector_values(1) = sparse_v1%elements(i)%value
                this%columns(i)%vector_values(2) = sparse_v2%elements(i)%value
                this%columns(i)%vector_values(3) = sparse_v3%elements(i)%value
                
            end if
        end do
    
    end subroutine sparse_matrix_init


    function sparse_matrix_count_nonzero_matrix_elements(this) result(count) 
        ! counts the number of nonzero vectors (columns) in a matrix

        implicit none

        class(sparse_matrix),intent(in) :: this

        integer :: i, count

        count = 0
        do i=1,this%sparse_num_cols
            if (abs(this%columns(i)%vector_values(1)) > 1.0e-12 .and. &
                abs(this%columns(i)%vector_values(2)) > 1.0e-12 .and. &
                abs(this%columns(i)%vector_values(3)) > 1.0e-12 ) then
                count = count + 1
            end if
        end do

    end function sparse_matrix_count_nonzero_matrix_elements


    function sparse_matrix_compress(this) result(result_matrix)
        ! takes a full vector that is sparse (has a lot of 0.0's) and compresses it to a sparse vector type

        implicit none 

        class(sparse_matrix),intent(in) :: this
        type(sparse_matrix) :: result_matrix
        integer :: i,count

        
        ! count how many nonzero elements there are
        ! count how many nonzero elements there are
        count = this%count_nonzero_matrix_elements()
            
        ! now that we have the number of non zero numbers, we can allocate the space for the sparse_vector
        result_matrix%sparse_num_cols = count
        result_matrix%full_num_cols = this%full_num_cols
        allocate(result_matrix%columns(count))

        ! populate the nonzero sparse_vector elements
        do i=1,this%full_num_cols
            if (abs(this%columns(:)%vector_values(1)) > 1.0e-12 .and. &
                abs(this%columns(:)%vector_values(2)) > 1.0e-12 .and. &
                abs(this%columns(:)%vector_values(3)) > 1.0e-12) then
                result_matrix%columns(i)%vector_values = this%columns(i)%vector_values
                result_matrix%columns(i)%full_index = i
            end if
        end do
        
        ! to save memory, maybe deallocate the original array after compressing it

    end function sparse_matrix_compress


    subroutine sparse_matrix_set_values(this, full_index, values) 
        ! given a full_vector index and a value, this adds or updates the corresponding value in the sparse vector

        implicit none

        class(sparse_matrix),intent(inout) :: this
        integer, intent(in):: full_index
        real,dimension(3),intent(inout) :: values
        
        integer :: i
        integer :: shift_index
        logical :: new_sparse_element_needed = .FALSE.

        do i=1,this%sparse_num_cols
            if (this%columns(i)%full_index == full_index) then
                this%columns(i)%vector_values(:) = values(:)  
                exit  
            else if (this%columns(i)%full_index > full_index) then
                new_sparse_element_needed = .TRUE.
                shift_index = i
                exit
            end if
        end do

        if (new_sparse_element_needed) then
            call this%add_element(full_index, values, shift_index)
        end if
       
    end subroutine sparse_matrix_set_values


    function sparse_matrix_get_values(this, full_index) result(values)
        ! given a full matrix column index, this function returns the 3 values at that index

        implicit none

        class(sparse_matrix),intent(inout) :: this
        integer :: full_index, i

        real,dimension(3) :: values 

        ! if the sparse matrix element has the same full index as the given full index, 
        ! returnthat element's value, if not, return the default zeros 
        do i=1,this%sparse_num_cols
            if (this%columns(i)%full_index == full_index) then
                values(:) = this%columns(i)%vector_values(:)
                exit
            else if (this%columns(i)%full_index > full_index) then
                values(:) = (/0.0, 0.0, 0.0/)
                exit
            end if  
        end do
       
    end function sparse_matrix_get_values


    subroutine sparse_matrix_add_element(this, full_index, values, shift_index)
        ! adds a sparse element to the current sparse matrix
        ! could be made more efficient by checking if the shift index is closer to the beginning or 
        ! the end. if closer to the beginning, cshift, then pull back the first few elements until shift_index
        
        implicit none

        class(sparse_matrix),intent(inout) :: this
        integer, intent(in) :: full_index, shift_index
        real,dimension(:),intent(in) :: values
        
        integer :: i,new_size

        new_size = this%sparse_num_cols + 1

        ! original information is preserved, but an extra space is allocated
        call this%increase_size()

        ! shift necessary elements up in the increased size matrix
        do i=this%sparse_num_cols,shift_index,-1
            this%columns(i+1)%vector_values(:) = this%columns(i)%vector_values(:)
            this%columns(i+1)%full_index = this%columns(i)%full_index
        end do

        ! place the new element in the correct spot
        this%columns(shift_index)%vector_values(:) = values(:)
        this%columns(shift_index)%full_index = full_index

    end subroutine sparse_matrix_add_element


    subroutine sparse_matrix_increase_size(this)
        ! increases the size of the sparse_matrix by 1 sparse element

        implicit none
        class(sparse_matrix),intent(inout) :: this

        type(sparse_matrix) :: temp_matrix
        integer :: new_size

        new_size = size(this) + 1
        ! allocate a temporary array
        allocate(temp_matrix%columns(new_size))

        ! copy array info
        temp_matrix%columns(1:new_size-1) = this%columns(:)

        ! move allocation of temporary array and overwrite this
        call move_alloc(temp_matrix, this)

        ! update the sparse number of columns
        this%sparse_num_cols = this%sparse_num_cols + 1
        

    end subroutine sparse_matrix_increase_size


    function sparse_matrix_vec_cross_matrix(this, vec) result(result_matrix)
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


    end function sparse_matrix_vec_cross_matrix


    function sparse_matrix_matrix_cross_vec(this, vec) result(result_matrix)
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

    end function sparse_matrix_matrix_cross_vec


    function sparse_matrix_dot_vec_and_matrix(this, vec) result(result_vector)
        ! dot product of a 3-vector and a "list" of vectors (derivative of a vector)
        ! returns a sparse_vector

        implicit none

        class(sparse_matrix),intent(inout) :: this
        real,intent(in) :: vec 

        type(sparse_vector),dimension(:),allocatable :: result_vector
        real :: temp_val
        real,dimension(3) :: temp_vec

        integer :: i

        ! set full and sparse column numbers/size
        result_vector%full_size = this%full_num_cols
        result_vector%sparse_size = this%sparse_num_cols

        ! allocate result matrix elements
        allocate(result_vector%elements(size(this%columns)))

        ! go through each sparse column index
        do i=1,this%sparse_num_cols
            ! do cross product
            temp_vec = this%columns(i)%vector_values(:)
            temp_val = dot(vec,temp_val)

            ! update results vector
            result_vector%elements(i)%value = temp_val
            result_vector%elements(i)%full_index = this%columns(i)%full_index

        end do       

    end function sparse_matrix_dot_vec_and_matrix


    subroutine sparse_matrix_multiply_scalar(this,scalar)
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

    end subroutine sparse_matrix_multiply_scalar


    subroutine sparse_matrix_sparse_add(this, matrix_a) 
        ! subroutine to add a sparse matrix from this
        ! this + matrix_a = this

        implicit none

        class(sparse_matrix),intent(inout) :: this
        type(sparse_matrix),intent(inout) :: matrix_a

        integer :: i

        ! make sure the input matrix has the same full size as this
        if (this%full_num_cols /= matrix_a%full_num_cols) then
            write(*,*) "Error!!! sparse_matrix_add requires input to have the same full_size. Quitting..."
            stop
        end if
        
        ! add matrix_a to 'this'
        do i=1, matrix_a%full_num_cols
            if (this%get_values(i) = (/0.0, 0.0, 0.0/) )
            
            this%columns(i)%vector_values =  this%columns(i)%vector_values + matrix_a%get_values(i)
        end do 
        
        

    end subroutine sparse_matrix_sparse_add


    subroutine sparse_matrix_sparse_subtract(this, matrix_a)
        ! function to add two sparse matrices (adding sets of vectors) returns a sparse matrix
        ! matrix_a - matrix b = result matrix

        implicit none
        
        class(sparse_matrix),intent(inout) :: this

        type(sparse_matrix),intent(inout) :: matrix_a
        type(sparse_matrix),intent(inout) :: matrix_b

        integer :: i

        ! allocate a full size result matrix
        allocate(this%columns(matrix_a%full_num_cols))
        
        ! add matrix_a and matrix_b vector values
        do i=1, matrix_a%full_num_cols
            this(i)%vector_values = matrix_a%get_values(i) - matrix_b%get_values(i)
        end do 
        
        ! collapse the 'this' matrix after this

    end subroutine sparse_matrix_sparse_subtract


    subroutine adjoint_init(this, body)
        ! initializes adjoint type

        implicit none

        class(adjoint),intent(inout) :: this
        type(surface_mesh),intent(in) :: body

        integer :: i, j

        this%size = body%N_verts*3

        call this%get_X_beta(body)

    end subroutine adjoint_init


    subroutine adjoint_get_X_beta(this, body)
        ! builds a vector of design variables
        
        implicit none

        class(adjoint),intent(inout) :: this
        type(surface_mesh),intent(in) :: body
    
        integer :: i, j
        
        ! build design variable vector X_beta
        allocate(this%X_beta(this%size))

        ! place all x values in X_beta, followed by all y, then z values 
        ! X_beta will look like this: x1, x2, x3, ...xn, y1, y2, y3, ...yn, z1, z2, z3, ...zn
        do i=1,3 
            do j=1,body%N_verts
                this%X_beta(j + (i-1)*body%N_verts) = body%vertices(j)%loc(i)
            end do
        end do


    end subroutine adjoint_get_X_beta



end module adjoint_mod