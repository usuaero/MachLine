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
        real :: mu, C_mu ! Mach angle
        real,dimension(3) :: c0 ! Compressibility axis (assumed in TriPan to be aligned with the freestream direction)
        logical,dimension(3) :: sym_about ! Whether the flow condition is symmetric about any plane
        real,dimension(3,3) :: B_mat_g, B_mat_c ! Dual metric matrix
        real,dimension(3,3) :: C_mat_g, C_mat_c ! Metric matrix
        logical :: supersonic, incompressible
        real,dimension(3,3) :: A_g_to_c, A_c_to_s, A_g_to_s ! Coordinate transformation matrices

        contains

            procedure :: init => flow_init
            procedure :: C_g_inner => flow_C_g_inner
            procedure :: point_in_dod => flow_point_in_dod

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
        this%incompressible = this%M_inf == 0.

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
        if (this%supersonic) then
            this%mu = asin(1.0/this%M_inf)
            this%C_mu = cos(this%mu)
        end if

        ! Assemble dual metric matrix
        ! Global
        do i=1,3
            this%B_mat_g(i,i) = 1.
        end do
        this%B_mat_g = this%B_mat_g-this%M_inf**2*outer(this%c0, this%c0)

        ! Compressible
        this%B_mat_c = 0.
        this%B_mat_c(1,1) = this%s*this%B**2
        this%B_mat_c(2,2) = 1.
        this%B_mat_c(3,3) = 1.
        
        ! Assemble metric matrix
        ! Global
        do i=1,3
            this%C_mat_g(i,i) = this%s*this%B**2
        end do
        this%C_mat_g = this%C_mat_g+this%M_inf**2*outer(this%c0, this%c0)

        ! Compressible
        this%C_mat_c = 0.
        this%C_mat_c(1,1) = 1.
        this%C_mat_c(2,2) = this%s*this%B**2
        this%C_mat_c(3,3) = this%s*this%B**2

        ! Check symmetry
        this%sym_about = this%v_inf == 0.

        ! Calculate transform from global to compressible coordinates

        ! Calculate transform from compressible to scaled coordinates


    end subroutine flow_init


    function flow_C_g_inner(this, a, b) result(c)
        ! Calculates the inner product a*C_g*b = c

        implicit none

        class(flow),intent(in) :: this
        real,dimension(3),intent(in) :: a, b
        real :: c

        c = inner(a, matmul(this%C_mat_g, b))

    end function flow_C_g_inner


    function flow_point_in_dod(this, Q, P) result(in_dod)
        ! Calculates whether the point Q lies in the dod of point P

        implicit none

        class(flow),intent(in) :: this
        real,dimension(3),intent(in) :: Q, P
        logical :: in_dod

        real,dimension(3) :: d

        ! Calculate displacement
        d = P-Q

        ! Check upstream
        if (inner(d, this%c0) >=0.) then

            ! Check in dod
            if (this%C_g_inner(d, d) >= 0.) then

                in_dod = .true.
                return

            end if

        end if

        in_dod = .false.
        
    end function flow_point_in_dod


end module flow_mod