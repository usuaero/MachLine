module flow_mod

    use json_mod
    use json_xtnsn_mod
    use math_mod

    implicit none
    
    type flow

        real,dimension(:),allocatable :: V_inf ! Freestream velocity
        real :: M_inf ! Freestream Mach number
        real :: gamma ! Ratio of specific heats
        real :: V_inf_mag
        real,dimension(3) :: u_inf

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
        call json_get(settings, 'V_inf', this%V_inf)
        call json_get(settings, 'M_inf', this%M_inf, found)
        call json_get(settings, 'gamma', this%gamma, found)

        ! Set default M (incompressible)
        if (.not. found) then
            this%M_inf = 0.0
        end if

        ! Set default gamma
        if (.not. found) then
            this%gamma = 1.4
        end if

        ! Derived quantities
        this%V_inf_mag = norm(this%V_inf)
        this%u_inf = this%V_inf/this%V_inf_mag

    end subroutine flow_init


end module flow_mod