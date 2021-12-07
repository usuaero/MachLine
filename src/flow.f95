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
        real :: c ! Freestream speed of sound
        real,dimension(3) :: c0 ! Compressibility direciton (freestream direction)
        logical,dimension(3) :: sym_about ! Whether the flow condition is symmetric about any plane
        real,dimension(3,3) :: psi ! Dual metric matrix, expressed in global coords

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
        this%c0 = this%v_inf*this%U_inv

        ! Calculate B
        if (this%M_inf == 1.) then
            write(*,*) "A freestream Mach number of 1.0 is not allowed in TriPan. Quitting..."
            stop
        else if (this%M_inf < 1.) then ! Subsonic
            this%B = sqrt(1.-this%M_inf**2)
        else ! Supersonic
            this%B = sqrt(this%M_inf**2-1.)
        end if

        ! Calculate freestream speed of sound
        this%c = this%M_inf*this%U

        ! Assemble dual metric matrix
        do i=1,3
            this%psi(i,i) = 1.
            this%psi(:,i) = this%psi(:,i) - this%M_inf**2*this%c0(:)*this%c0(i)
        end do

        ! Check symmetry
        this%sym_about = this%v_inf == 0.

    end subroutine flow_init


end module flow_mod