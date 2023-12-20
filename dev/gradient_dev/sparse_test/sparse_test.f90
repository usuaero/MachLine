program sparse_test
    ! tests sparse vector and sparse matrix operations
    
    use adjoint_mod
    
    implicit none

    type(sparse_vector),dimension(:),allocatable :: sparse_vector_a, sparse_vector_b 
    type(sparse_matrix),dimension(:),allocatable :: sparse_matrix_a, sparse_matrix_b
    real,dimension(:),allocatable :: vector_a, vector_b, vector_x, vector_y, vector_z

    ! initialize full vector
    vector_a = (/ 1.0, 0.0, 2.0, 0.0, 2.0, -1.0, 0.0, 0.1, 2.0, 0.0, 1.0 /)
    write(*, '(f14.10, 2x)') vector_a


end program sparse_test