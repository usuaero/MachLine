module flow_mod

    use json_mod
    use json_xtnsn_mod
    use math_mod

    implicit none
    
    type flow

        real,dimension(:),allocatable :: v_inf ! Freestream velocity
        real :: M_inf ! Freestream Mach number
        real :: gamma ! Ratio of specific heats
        real :: U, U_inv ! Freestream velocity magnitude
        real :: B ! Compressibility scale factor
        real :: s ! Sign of 1-M^2; determines character of governing PDE (hyperbolic (s=-1) vs elliptic(s=1))
        real :: c ! Freestream speed of sound
        real :: mu ! Mach angle
        real,dimension(3) :: c0 ! Compressibility axis (assumed in TriPan to be aligned with the freestream direction)
        logical,dimension(3) :: sym_about ! Whether the flow condition is symmetric about any plane
        real,dimension(3,3) :: psi ! Dual metric matrix, expressed in global coords
        logical :: supersonic
        real,dimension(3,3) :: A_g_to_c, A_c_to_s, A_g_to_s ! Coordinate transformation matrices

        contains

            procedure :: init => flow_init

    end type flow


contains


    subroutine flow_init(this, settings)

        implicit none

        class(flow),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings

        logical :: found
        integer :: i

        ! Get flow params
        call json_get(settings, 'freestream_velocity', this%v_inf, found)
        if (.not. found) then
            write(*,*) "Freestream velocity was not specified. Quitting..."
            stop
        end if
        call json_xtnsn_get(settings, 'freestream_mach_number', this%M_inf, 0.0)
        call json_xtnsn_get(settings, 'gamma', this%gamma, 1.4)

        ! Derived quantities
        this%U = norm(this%v_inf)
        this%U_inv = 1./this%U
        this%c0 = -this%v_inf*this%U_inv

        ! Determine condition
        if (this%M_inf == 1.) then
            write(*,*) "A freestream Mach number of 1.0 is not allowed in TriPan. Quitting..."
            stop
        end if
        this%supersonic = this%M_inf > 1.0

        ! Calculate B and s
        if (this%supersonic) then
            this%B = sqrt(this%M_inf**2-1.)
            this%s = -1.
        else
            this%B = sqrt(1.-this%M_inf**2)
            this%s = 1.
        end if

        ! Calculate freestream speed of sound
        this%c = this%M_inf*this%U

        ! Calculate Mach angle
        this%mu = asin(1.0/this%M_inf)

        ! Assemble dual metric matrix
        do i=1,3
            this%psi(i,i) = 1.
        end do
        this%psi = this%psi-this%M_inf**2*outer(this%c0, this%c0)

        ! Check symmetry
        this%sym_about = this%v_inf == 0.

    end subroutine flow_init


end module flow_mod