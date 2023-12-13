! module for calculating the gradient of CF, CM with respect to mesh points via the adjoint method
module adjoint_mod

    use linked_list_mod
    use math_mod
    use helpers_mod
    use base_geom_mod
    use mesh_mod
    use surface_mesh_mod

    implicit none

    type sparse_element
    

        real :: value             ! non-zero value of an element in a sparse matrix
        integer :: full_index     ! index of value in the full matrix

    end type sparse_element


    type sparse_vector
        ! used to save memory when we have many sparse vectors (most of the vector elements are 0) 
        ! "full vector" has all the zero and nonzero values
        ! "sparse vector" has only the nonzero values, with their corresponding positions in the full_vector

        integer :: sparse_size,full_size    ! sparse_size is the number of non-zero numbers, full_size is full vector
        type(sparse_element),dimension(:),allocatable :: element        ! non zero values in sparse vector

        contains
            procedure :: init => sparce_vector_init
            procedure :: get_value => sparce_vector_get_value
            procedure :: set_value => sparce_vector_set_value
            
            ! these should maybe go in the math module
            procedure :: sparse_add => sparce_vector_sparse_add
            procedure :: sparse_subtract => sparce_vector_sparse_subtract
            procedure :: sparse_inner => sparse_vector_sparse_inner
            procedure :: sparse_cross => sparse_vector_sparse_cross

            procedure :: fill_vector => sparse_vector_fill_vector


    end type sparse_vector




    type adjoint

        integer :: size    ! number of design varibles, the length of X_beta
        ! we may never use X_beta
        real,dimension(:),allocatable :: X_beta  ! list of x y z values of all mesh points (design variables)
    

        contains
            procedure :: init => adjoint_init
            procedure :: get_X_beta => adjoint_get_X_beta ! puts all x y z values of mesh points in a list

    end type adjoint


contains


    subroutine sparse_vector_init(this, full_size)
        ! initializes a sparse_vector

        class(sparse_vector),intent(inout) :: this
        integer :: full_size

        this%full_size = full_size

    end subroutine sparse_vector_init


    function sparse_vector_get_value(this, full_index) result(value)
        ! given a full vector index, this function returns the full vector value at that index

        implicit none

        class(sparse_vector),intent(inout) :: this
        integer :: full_index
        real :: value = 0.0

        ! if the sparse vector element has the same full index as the given full index value is set to
        ! that element's value, if not, the value stays zero 
        do i=1,this%sparse_size
            if (this%elements(i)%full_index == full_index) then
                value = this%elements(i)%value
                exit
            end if  
        end do
       
    end function sparse_vector_get_value


    subroutine adjoint_init(this, body)
        ! initializes adjoint type

        implicit none

        class(adjoint),intent(inout) :: this
        type(surface_mesh),intent(in) :: body

        integer :: i, j, N_verts

        N_verts = body%N_verts
        this%size = N_verts*3

        call this%get_X_beta(body)

    end subroutine adjoint_init


    subroutine adjoint_get_X_beta(this, body)
        ! builds a vector of design variables
        
        implicit none

        class(adjoint),intent(inout) :: this
        type(surface_mesh),intent(in) :: body
        
        N_verts = body%N_verts

        ! build design variable vector X_beta
        allocate(this%X_beta(this%size))


        ! place all x values in X_beta, followed by all y, then z values 
        ! X_beta will look like this: x1, x2, x3, ...xn, y1, y2, y3, ...yn, z1, z2, z3, ...zn
        do i=1,3 
            do j=1,N_verts
                this%X_beta(j + (i-1)*N_verts) = body%vertices(j)%loc(i)
            end do
        end do


    end subroutine adjoint_get_X_beta



end module adjoint_mod