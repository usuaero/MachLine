! module for calculating the gradient of CF, CM with respect to mesh points via the adjoint method
module adjoint_mod

    use linked_list_mod
    use math_mod
    use helpers_mod
    use base_geom_mod

    implicit none

    contains

            procedure :: calc_d_loc => adjoint_calc




end module adjoint_mod