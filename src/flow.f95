module flow_mod

    use json_mod
    use json_xtnsn_mod
    use math_mod
    use linalg_mod
    use helpers_mod

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
        real :: K, K_inv ! Kappa factor
        real,dimension(3) :: c_hat_g ! Compressibility axis (assumed in TriPan to be aligned with the freestream direction)
        logical,dimension(3) :: sym_about ! Whether the flow condition is symmetric about any plane
        real,dimension(3,3) :: B_mat_g, B_mat_c, B_mat_g_inv ! Dual metric matrix
        real,dimension(3,3) :: C_mat_g, C_mat_c ! Metric matrix
        logical :: supersonic, incompressible
        real,dimension(3,3) :: A_g_to_c, A_c_to_s, A_g_to_s ! Coordinate transformation matrices

        contains

            procedure :: init => flow_init
            procedure :: calc_metric_matrices => flow_calc_metric_matrices
            procedure :: calc_transforms => flow_calc_transforms
            procedure :: B_g_inner => flow_B_g_inner
            procedure :: C_g_inner => flow_C_g_inner
            procedure :: point_in_dod => flow_point_in_dod

    end type flow


contains


    subroutine flow_init(this, settings, spanwise_axis)

        implicit none

        class(flow),intent(inout) :: this
        type(json_value),pointer,intent(in) :: settings
        character(len=:),allocatable,intent(in) :: spanwise_axis

        logical :: found
        integer :: i

        ! Get flow params
        call json_get(settings, 'freestream_velocity', this%v_inf, found)
        if (.not. found) then
            write(*,*) "!!! Freestream velocity was not specified. Quitting..."
            stop
        end if
        call json_xtnsn_get(settings, 'freestream_mach_number', this%M_inf, 0.)
        call json_xtnsn_get(settings, 'gamma', this%gamma, 1.4)

        ! Check for positive Mach number (idiot-proofing)
        if (this%M_inf < 0.) then
            write(*,*) "!!! Invalid freestream Mach number selected. Cannot be a negative number. Quitting..."
            stop
        end if

        ! Check symmetry
        this%sym_about = this%v_inf == 0.

        ! Derived quantities
        this%U = norm2(this%v_inf)
        this%U_inv = 1./this%U
        this%c_hat_g = this%v_inf*this%U_inv

        ! Determine condition
        if (this%M_inf == 1.) then
            write(*,*) "!!! A freestream Mach number of 1.0 is not allowed in MachLine. Quitting..."
            stop
        end if
        this%supersonic = this%M_inf > 1.0
        this%incompressible = this%M_inf == 0.

        ! Calculate B and s
        if (this%supersonic) then
            this%B = sqrt(this%M_inf**2 - 1.)
            this%s = -1.
            this%K = 2.*pi
        else
            this%B = sqrt(1. - this%M_inf**2)
            this%s = 1.
            this%K = 4.*pi
        end if
        this%K_inv = 1./this%K

        ! Calculate freestream speed of sound
        this%c = this%M_inf*this%U

        ! Calculate Mach angle (if supersonic; it is meaningless otherwise)
        if (this%supersonic) then
            this%mu = asin(1.0/this%M_inf)
            this%C_mu = cos(this%mu)
        end if

        ! Calculate relevant matrices
        call this%calc_metric_matrices()
        call this%calc_transforms(spanwise_axis)

    end subroutine flow_init


    subroutine flow_calc_metric_matrices(this)

        implicit none

        class(flow),intent(inout) :: this

        integer :: i

        ! Assemble dual metric matrix
        ! Global (E&M Eq. (E.3.9))
        do i=1,3
            this%B_mat_g(i,i) = 1.
        end do
        this%B_mat_g = this%B_mat_g - this%M_inf**2*outer(this%c_hat_g, this%c_hat_g)

        ! Invert (used for source-free formulation)
        call matinv(3, this%B_mat_g, this%B_mat_g_inv)

        ! Compressible (E&M Eq. (E.3.8))
        this%B_mat_c = 0.
        this%B_mat_c(1,1) = this%s*this%B**2
        this%B_mat_c(2,2) = 1.
        this%B_mat_c(3,3) = 1.
        
        ! Assemble metric matrix
        ! Global (E&M Eq. (E.3.9))
        do i=1,3
            this%C_mat_g(i,i) = 1.-this%M_inf**2
        end do
        this%C_mat_g = this%C_mat_g + this%M_inf**2*outer(this%c_hat_g, this%c_hat_g)

        ! Compressible (E&M Eq. (E.3.8))
        this%C_mat_c = 0.
        this%C_mat_c(1,1) = 1.
        this%C_mat_c(2,2) = this%s*this%B**2
        this%C_mat_c(3,3) = this%s*this%B**2
    
    end subroutine flow_calc_metric_matrices


    subroutine flow_calc_transforms(this, spanwise_axis)

        implicit none

        class(flow),intent(inout) :: this
        character(len=:),allocatable,intent(in) :: spanwise_axis

        real,dimension(3) :: j_g, c_hat_c, c_hat_s

        ! Set positive spanwise axis
        j_g = 0.
        select case (spanwise_axis)

        case ('+x')
            j_g(1) = 1.

        case ('-x')
            j_g(1) = -1.

        case ('+y')
            j_g(2) = 1.

        case ('-y')
            j_g(2) = -1.

        case ('+z')
            j_g(3) = 1.

        case ('-z')
            j_g(3) = -1.

        case default
            j_g(2) = 1.

        end select

        ! Calculate transform from global to compressible coordinates
        this%A_g_to_c = 0.
        this%A_g_to_c(1,:) = this%c_hat_g
        this%A_g_to_c(3,:) = cross(this%c_hat_g, j_g)
        this%A_g_to_c(3,:) = this%A_g_to_c(3,:)/norm2(this%A_g_to_c(3,:))
        this%A_g_to_c(2,:) = cross(this%A_g_to_c(3,:), this%c_hat_g)

        ! Calculate transform from compressible to scaled coordinates
        this%A_c_to_s = 0.
        this%A_c_to_s(1,1) = 1.
        this%A_c_to_s(2,2) = this%B
        this%A_c_to_s(3,3) = this%B

        ! Check calculation
        if (run_checks) then
            c_hat_c = matmul(this%A_g_to_c, this%c_hat_g)
            if (abs(c_hat_c(1)-1.)>1e-12 .or. abs(c_hat_c(2))>1e-12 .or. abs(c_hat_c(3))>1e-12) then
                write(*,*) "!!! Transformation to the compressible coordinate system failed. Quitting..."
                stop
            end if
        end if

        ! Calculate transform from global to scaled coordinates
        this%A_g_to_s = matmul(this%A_c_to_s, this%A_g_to_c)

    end subroutine flow_calc_transforms


    function flow_B_g_inner(this, a, b) result(c)
        ! Calculates the inner product a*B_g*b = c

        implicit none

        class(flow),intent(in) :: this
        real,dimension(3),intent(in) :: a, b
        real :: c

        c = inner(a, matmul(this%B_mat_g, b))

    end function flow_B_g_inner


    function flow_C_g_inner(this, a, b) result(c)
        ! Calculates the inner product a*C_g*b = c

        implicit none

        class(flow),intent(in) :: this
        real,dimension(3),intent(in) :: a, b
        real :: c

        c = inner(a, matmul(this%C_mat_g, b))

    end function flow_C_g_inner


    function flow_point_in_dod(this, Q, P) result(in_dod)
        ! Calculates whether the point Q lies in the DoD of point P

        implicit none

        class(flow),intent(in) :: this
        real,dimension(3),intent(in) :: Q, P
        logical :: in_dod

        real,dimension(3) :: d

        in_dod = .false.

        ! Calculate displacement
        d = P-Q

        ! Check upstream
        if (inner(d, this%c_hat_g) >= 0.) then ! E&M Eq. (J.3.1)

            ! Check in dod
            if (this%C_g_inner(d, d) >= 0.) then ! E&M Eq. (J.3.2)

                in_dod = .true.

            end if

        end if
        
    end function flow_point_in_dod


end module flow_mod