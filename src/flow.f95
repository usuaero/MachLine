module flow_mod

    use json_mod
    use json_xtnsn_mod
    use math_mod

    implicit none
    
    type flow

        real,dimension(:),allocatable :: v_inf ! Freestream velocity
        real :: M_inf ! Freestream Mach number
        real :: gamma ! Ratio of specific heats
        real :: U ! Freestream velocity magnitude
        real,dimension(3) :: c0 ! Compressibility direciton (freestream direction)

        contains

            procedure :: init => flow_init

    end type flow


contains


    subroutine flow_init(this, settings)

        implicit none

        class(flow),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        logical :: found

        ! Get flow params
        call json_get(settings, 'V_inf', this%v_inf)
        call json_xtnsn_get(settings, 'M_inf', this%M_inf, 0.0)
        call json_xtnsn_get(settings, 'gamma', this%gamma, 1.4)

        ! Derived quantities
        this%U = norm(this%v_inf)
        this%c0 = this%v_inf/this%U

    end subroutine flow_init


end module flow_mod