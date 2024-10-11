module filament_segment_mod
    use base_geom_mod
    use flow_mod 
    use helpers_mod
    use mesh_mod 
    use math_mod 

!!!! mirrors panel module

    implicit none


    type filament_segment_integrals 
        ! Container type for the fundamental integrals used to calculate influence coefficients
        
        integer :: r, s, rs ! Parameters that will be most convenient to keep here (flow type indicators)
        real :: u, v, w ! unrelaxed influence integrals in the compressible frame. !!!! may need them to be allocatable


    end type filament_segment_integrals


    type filament_dod
    ! container for whether a filament lies in a point's domain of dependence
    
        logical :: first_in_dod = .false. !!!! the nomenclature can certainly change here
        logical :: both_in_dod = .false. !!!! the nomenclature can certainly change here

    end type filament_dod


    type filament_segment 
        ! A section of a filament that contains two vertices. 
        type(vertex_pointer), dimension(:),allocatable :: vertices 
        integer :: index ! index of this filament in the filament mesh array
        integer :: N = 2 ! number of vertices
        real,dimension(3,3) :: A_g_to_c, A_c_to_g ! Coordinate transformation matrices
        real,dimension(3,3) :: A_g_to_c_mir, A_c_to_g_mir
        integer, dimension(4) :: parents ! stores index for parent vertices.  Order is top1 bot1 top2 bot2
        integer :: M_dim



        contains
            procedure :: init => filament_segment_init 
            procedure :: calc_integrals => filament_segment_calc_integrals
            ! Getters
            procedure :: get_vertex_index => filament_segment_get_vertex_index
            !velocities
            procedure :: calc_velocity_influences => filament_segment_calc_velocity_influences
            procedure :: assemble_v_d_M_space => filament_segment_assemble_v_d_M_space 
            procedure :: calc_influence_integrals => filament_segment_calc_influence_integrals 
            procedure :: get_vertex_loc => filament_segment_get_vertex_loc
            procedure :: get_doublet_strengths => filament_segment_get_doublet_strengths
            procedure :: calc_velocities => filament_segment_calc_velocities
            ! DOD stuff
            procedure :: check_dod => filament_segment_check_dod
    end type filament_segment 


contains


    subroutine filament_segment_init(this,v1,v2,index,input_parents)
        ! initializes important filament segment parameters 

        implicit none

        class(filament_segment), intent(inout) :: this
        type(vertex),intent(inout),target :: v1, v2
        integer, intent(in) :: index
        integer,dimension(4),intent(in) :: input_parents

        ! set number of sides 
        this%N = 2

        ! allocate vertex array
        allocate(this%vertices(this%N))

        ! assign pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        
        !store parent vertex indicies
        this%parents = input_parents

        ! store index
        this%index = index

        ! still need to calculated derived geometry (based on what is needed in solver)
    end subroutine filament_segment_init


    ! likely need filament_segment_init_with_flow?


    function filament_segment_get_vertex_loc(this, i) result(loc)

        implicit none
        class(filament_segment), intent(in) :: this
        integer,intent(in) :: i
        real,dimension(3) :: loc

        loc = this%vertices(i)%ptr%loc

    end function filament_segment_get_vertex_loc


    function filament_segment_get_vertex_index(this, i) result(index)
        ! gets a vertex index

        implicit none

        class(filament_segment),intent(in) :: this
        integer,intent(in) :: i
        integer :: index
        index = this%vertices(i)%ptr%index

    end function filament_segment_get_vertex_index


    ! subroutine filament_segment_calc_g_to_c_transform(this, freestream)
    !     ! Calculates the necessary transformations to move from global to compressible coordinates (Eq. (E.0.1) in Epton and Magnus)
        
    ! end subroutine filament_segment_calc_g_to_c_transform


    subroutine filament_segment_calc_velocity_influences(this,eval_point, freestream, mirror_filament, v_d_M_space)
        ! calculates the velocity influence of a vortex filament 
        
        implicit none 

        class(filament_segment),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point 
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_filament
        real, dimension(:,:),allocatable,intent(out) :: v_d_M_space

        type(filament_dod) :: dod_info !!!! might need to make a filament_dod type in this mod (filament_segment_mod) 
        type(filament_segment_integrals) :: int 
        real,dimension(:,:),allocatable :: v_d_mu_space
        real :: x2, y2 !!!! not sure if this is needed
        integer :: i
        
        dod_info = this%check_dod(eval_point, freestream, mirror_filament)
        
        ! Get integrals 
        int = this%calc_integrals(eval_point, freestream, mirror_filament, dod_info)
        v_d_M_space =  this%assemble_v_d_M_space(int, freestream, mirror_filament)
        
    end subroutine filament_segment_calc_velocity_influences


    function filament_segment_check_dod(this, eval_point, freestream, mirror_filament) result(dod_info)
        ! Determines how (if) this panel lies within the domain of dependence of the evaluation point
        
        implicit none 

        class(filament_segment),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream 
        logical,intent(in),optional :: mirror_filament 
        type(filament_dod) :: dod_info

        logical :: mirrored, point1_inside,point2_inside

        if (present(mirror_filament)) then
            mirrored = mirror_filament
        else
            mirrored = .false.
        end if  


        if (freestream%supersonic) then
            
            ! check if vertices are in DOD
            point1_inside = freestream%point_in_dod(this%get_vertex_loc(1),eval_point)
            point2_inside = freestream%point_in_dod(this%get_vertex_loc(2),eval_point)
            
            if (point1_inside .and. point2_inside) then
                dod_info%both_in_dod = .true.
            else if (point1_inside .or. point2_inside) then
                dod_info%first_in_dod = .true.
                dod_info%both_in_dod = .false.
            else
                dod_info%both_in_dod = .false.
                dod_info%first_in_dod = .false.
            end if
        ! Will always be entirely inside for subsonic flow 
        else
            dod_info%both_in_dod = .true.
        end if
        
    end function filament_segment_check_dod
         

    function filament_segment_calc_integrals(this, eval_point, freestream, mirror_filament, dod_info) result(int)
        ! Calculates the integrals necessary for the given influence
        
        implicit none 
        
        class(filament_segment),intent(in) :: this 
        real,dimension(3),intent(in) :: eval_point
        ! character(len=*),intent(in) :: influence_type !!!! don't think this will apply - SA 
        type(flow), intent(in) :: freestream
        logical,intent(in) :: mirror_filament

        type(filament_segment_integrals) :: int
        type(filament_dod),intent(in) :: dod_info

        int%s = freestream%s 
        int%rs = int%r*int%s 

        ! Calculate necessary integrals based on the flow condition
        call this%calc_influence_integrals(eval_point, dod_info, freestream, mirror_filament, int) !!!! I think the int input in this subroutine call is what will output from this function
        !!!! should dod_info carry whether the dod is supersonic or subsonic? -SA 

 
    end function filament_segment_calc_integrals 


    subroutine filament_segment_calc_influence_integrals(this, eval_point, dod_info, freestream, mirror_filament, int)
        ! calculates the supersonic and subsonic vortex filament influence integrals

        implicit none 

        class(filament_segment),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(filament_dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream 
        logical,intent(in) :: mirror_filament
        type(filament_segment_integrals),intent(inout) :: int
        real :: unrel_v_numer_one, unrel_v_numer_two
        real :: unrel_v_denom_one, unrel_v_denom_two
        real :: unrel_w_numer_one, unrel_w_numer_two
        real :: unrel_w_denom_one, unrel_w_denom_two
        real :: xf, xi, yf, yi, zf, zi
        real :: xo, yo, zo
        real :: bsq  
        real,dimension(3) :: loc_1, loc_2,r,V_local,V_global
        integer :: k 
        real :: fac = 1
        real :: c_a, c_b, s_a, s_b, a,b,L
        real :: C1, C2, C3
        real :: x0, y0, z0
        real :: y, z
        real,dimension(3,3) :: A_c_to_f
        real :: term1, term2, denominator, v_numerator, w_numerator, tx, ty, tz,t, R_bf, R_bi,f_c
        real :: tol = 1e-4
        f_c = 0.2
        bsq = 1 - (freestream%M_inf**2)

        if (dod_info%both_in_dod .or. dod_info%first_in_dod) then 
            loc_1 = this%get_vertex_loc(1) !!!! put in appropriate input and get 
            loc_2 = this%get_vertex_loc(2) !!!! put in appropriate input and get 

            loc_1 = matmul(freestream%A_g_to_c, loc_1) !!!! Did we do this right? 
            loc_2 = matmul(freestream%A_g_to_c, loc_2)
            r = matmul(freestream%A_g_to_c, eval_point)
            
            xf = loc_2(1)
            xi = loc_1(1)
            
            yf = loc_2(2)
            yi = loc_1(2)
            

            zf = loc_2(3)
            zi = loc_1(3)

            x0 = r(1)
            y0 = r(2)
            z0 = r(3)

            ! check to see if cp is on the line
            if ((xf - xi)>tol) tx = (x0 - xi) / (xf - xi) 
            if ((yf - yi)>tol) ty = (y0 - yi) / (yf - yi) 
            if ((zf - zi) >tol) tz = (z0 - zi) / (zf - zi) 

            ! Check if tx, ty, and tz are approximately equal
            if (abs(tx - ty) < tol .and. abs(ty - tz) < tol) then
                t = tx  ! All t values are approximately the same
                if (t >= -tol .and. t <= 1.0 + tol) then
                    write(*,*) "Point is on the line"
                    int%u = 0
                    int%v = 0
                    int%w = 0
                    return
                endif
            endif

            ! Check if t is within the segment range [0, 1] with tolerance
            

        ! compute trig functions
            L = ((xf-xi)**2 + (yf-yi)**2 + (zf-zi)**2)**0.5
            a = asin((zf-zi)/L)
            b = asin((yf-yi)/L)
            
            c_a = cos(a)
            c_b = cos(b)
            s_a = sin(a)
            s_b = sin(b)
            ! transform to filament frame
            ! A_c_to_f = [[c_a * c_b, -s_a, c_a * s_b],[s_a * c_b, c_a, s_a * s_b],[-s_b, 0., c_b]]
            A_c_to_f(1,1) = c_a * c_b
            A_c_to_f(1,2) = -s_a
            A_c_to_f(1,3) = c_a * s_b
            A_c_to_f(2,1) = s_a * c_b
            A_c_to_f(2,2) = c_a
            A_c_to_f(2,3) = s_a * s_b
            A_c_to_f(3,1) = -s_b
            A_c_to_f(3,2) = 0.
            A_c_to_f(3,3) = c_b

            A_c_to_f = transpose(A_c_to_f)
            
            ! write(*,*) "before", loc_1, loc_2, r
            loc_1 = matmul(A_c_to_f, loc_1)
            loc_2 = matmul(A_c_to_f, loc_2)
            r = matmul(A_c_to_f, r)
            ! write(*,*) "after", loc_1, loc_2, r

            xf = loc_2(1)
            xi = loc_1(1)
            
            yf = loc_2(2)
            yi = loc_1(2)
            

            zf = loc_2(3)
            zi = loc_1(3)

            y = yi
            z = zi
            x0 = r(1)
            y0 = r(2)
            z0 = r(3)

        
            

        
            ! compute C coefs
            C1 = bsq * s_b**2 + c_a**2 * c_b**2 + bsq * s_a**2 * c_b**2
            C2 = ((bsq - 1) * 2 * s_a * c_a * c_b * (z - z0) - 2 * bsq * c_b * s_a**2 * s_b * (y - y0) + 2 * bsq * s_b * c_b &
            * (y - y0) - 2 * c_a**2 * c_b * s_b * (y - y0))
            C3 = (bsq * s_a**2 * s_b**2 + bsq * c_b**2 + c_a**2 * s_b**2) * (y-y0)**2 + (1 - bsq) * 2 * c_a * s_a * s_b * (y - y0) &
            * (z - z0) + (bsq * c_a**2 + s_a**2) * (z - z0)**2

            R_bf = C1 * (xf - x0)**2 + C2 * (xf - x0) + C3
            R_bi = C1 * (xi - x0)**2 + C2 * (xi - x0) + C3

            if (R_bf < f_c .and. R_bi < f_c) then
                int%u = 0
                int%v = 0
                int%w = 0
                return
            end if

            ! get terms for influence coefficients
            
            ! v local
            v_numerator =  bsq * (z - z0)
            denominator = freestream%K/2 * (4 * C1 * C3 - C2**2) ! k is 2 k pi   Needed pi k)
            
            
            ! w local
            w_numerator = - bsq * (y - y0)
            

            if (dod_info%both_in_dod .and. (R_bf > f_c .and. R_bi > f_c)) then 
                term1 = (2 * C1 * (xf - x0) + C2) / ((C1 * (xf - x0)**2 + C2 * (xf - x0) + C3)**0.5)
                term2 = (2 * C1 * (xi - x0) + C2) / ((C1 * (xi - x0)**2 + C2 * (xi - x0) + C3)**0.5)
                V_local(1) = 0
                V_local(2) = (v_numerator / denominator) * (term1 - term2)
                V_local(3) =  (w_numerator / denominator) * (term1 - term2)
            else if (dod_info%first_in_dod .and. (R_bi > f_c .and. R_bf < f_c)) then 
                term2 = (2 * C1 * (xi - x0) + C2) / ((C1 * (xi - x0)**2 + C2 * (xi - x0) + C3)**0.5)
                V_local(1) = 0
                V_local(2) = (v_numerator / denominator) * (- term2)
                V_local(3) = (w_numerator / denominator) * (- term2)
            else if (R_bi < f_c .and. R_bf > f_c) then
                ! write(*,*) "This should not happen"
                term1 = (2 * C1 * (xf - x0) + C2) / ((C1 * (xf - x0)**2 + C2 * (xf - x0) + C3)**0.5)
                V_local(1) = 0
                V_local(2) = (v_numerator / denominator) * (term1)
                V_local(3) = (w_numerator / denominator) * (term1)
            else
                V_local(1) = 0
                V_local(2) = 0
                V_local(3) = 0
            end if  

            V_global = matmul(transpose(A_c_to_f), V_local)
            
            int%u = V_global(1)
            int%v = V_global(2)
            int%w = V_global(3)
        else 
            int%u = 0
            int%v = 0
            int%w = 0
        end if

        ! ! Miranda straight filament

        ! xf = loc_2(1)
        ! xi = loc_1(1)
        
        ! yf = loc_2(2)
        ! yi = loc_1(2)
        

        ! zf = loc_2(3)
        ! zi = loc_1(3)
        ! xo = r(1)
        ! yo = r(2)
        ! zo = r(3)
        ! ! write out the integrals here based on the different dod stuff. 
        ! if (-bsq*((yo-yf)**2+(zo-zf)**2)>fac*(xo-xf)**2) then
        !     unrel_v_numer_one = 0
        !     unrel_v_denom_one = 1
        !     unrel_w_numer_one = 0
        !     unrel_w_denom_one = 1
        ! else
        !     unrel_v_numer_one = (zo-zf)*(xo-xf)
        !     unrel_v_denom_one = ((yo-yf)**2+(zo-zf)**2)*((xo-xf)**2+bsq*((yo-yf)**2+(zo-zf)**2))**0.5
        !     unrel_w_numer_one = (yo-yf)*(xo-xf)
        !     unrel_w_denom_one = ((yo-yf)**2+(zo-zf)**2)*((xo-xf)**2+bsq*((yo-yf)**2+(zo-zf)**2))**0.5
        ! end if

        ! if (-bsq*((yo-yi)**2+(zo-zi)**2)>fac*(xo-xi)**2) then
        !     unrel_v_numer_two = 0
        !     unrel_v_denom_two = 1
        !     unrel_w_numer_two = 0
        !     unrel_w_denom_two = 1
        ! else
        !     unrel_v_numer_two = (zo-zi)*(xo-xi)
        !     unrel_v_denom_two = ((yo-yi)**2+(zo-zi)**2)*((xo-xi)**2+bsq*((yo-yi)**2+(zo-zi)**2))**0.5
        !     unrel_w_numer_two = (yo-yi)*(xo-xi)
        !     unrel_w_denom_two = ((yo-yi)**2+(zo-zi)**2)*((xo-xi)**2+bsq*((yo-yi)**2+(zo-zi)**2))**0.5
        ! end if

        ! if (dod_info%both_in_dod) then 
        !     int%u = 0
        !     int%v = (freestream%k_inv)*((unrel_v_numer_one/unrel_v_denom_one)-(unrel_v_numer_two/unrel_v_denom_two)) !!!! declare all of these variables and stuff - SA
        !     int%w = (-freestream%k_inv)*((unrel_w_numer_one/unrel_w_denom_one)-(unrel_w_numer_two/unrel_w_denom_two)) !!!! declare all of these variables and stuff - SA
        ! else if (dod_info%first_in_dod) then 
        !     int%u = 0
        !     int%v = (-freestream%k_inv)*((unrel_v_numer_two/unrel_v_denom_two)) !!!! declare all of these variables and stuff - SA
        !     int%w = (freestream%k_inv)*((unrel_w_numer_two/unrel_w_denom_two)) !!!! declare all of these variables and stuff - SA
        ! else 
        !     int%u = 0
        !     int%v = 0
        !     int%w = 0

        ! end if 




        

    end subroutine filament_segment_calc_influence_integrals


    function filament_segment_assemble_v_d_M_space(this, int, freestream, mirror_filament) result(v_d_M_space) !!!! need to create this space for wakes
        ! Assembles the doublet-induced velocity influence coefficient matrix from the previously-calculated influence integrals
        implicit none

        class(filament_segment),intent(in) :: this
        type(filament_segment_integrals),intent(in) :: int
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_filament

        real,dimension(:,:),allocatable :: v_d_M_space

        real,dimension(:,:),allocatable :: v_d_mu_space

        ! Allocate space
        allocate(v_d_M_space(3,4), source=0.)

        v_d_M_space(1,1) = int%u
        v_d_M_space(1,2) = -int%u
        v_d_M_space(1,3) = -int%u
        v_d_M_space(1,4) = int%u

        v_d_M_space(2,1) = int%v
        v_d_M_space(2,2) = -int%v
        v_d_M_space(2,3) = -int%v
        v_d_M_space(2,4) = int%v

        v_d_M_space(3,1) = int%w
        v_d_M_space(3,2) = -int%w
        v_d_M_space(3,3) = -int%w
        v_d_M_space(3,4) = int%w
            

        ! ! Wake bottom influence is opposite the top influence
        ! v_d_M_space(:,this%M_dim+1:this%M_dim*2) = -v_d_M_space(:,1:this%M_dim)

        ! Transform to global coordinates
        ! if (mirror_filament) then !!!! need to get mirroring going 
        !     v_d_M_space = matmul(transpose(freestream%A_g_to_c_mir), v_d_M_space) !!!! need to go from compressible to global
        ! else
        v_d_M_space = matmul(transpose(freestream%A_g_to_c), v_d_M_space) !!!! need to go from compressible to global
        ! end if
        
    end function filament_segment_assemble_v_d_M_space


    function filament_segment_get_doublet_strengths(this, mu, mirror, asym_flow) result(mu_strengths)
        ! Returns the relevant doublet strengths for this panel

        implicit none
        
        class(filament_segment),intent(in) :: this
        real,dimension(:),allocatable,intent(in) :: mu
        logical,intent(in) :: mirror
        logical,intent(in) :: asym_flow

        real,dimension(:),allocatable :: mu_strengths

        integer :: shift, i, i_top1, i_bot1, i_top2, i_bot2

        ! Determine shift
        if (mirror) then
            shift = size(mu)/2
        else
            shift = 0
        end if

        ! Allocate
        allocate(mu_strengths(this%N)) !!!! should M_dim be 3 or 2? -SA I just switched it to N -JJH

        ! Get doublet strengths based on parents
        
        ! do i=1,4
            i_top1 = this%parents(1) + shift
            i_bot1 = this%parents(2) + shift
            i_top2 = this%parents(3) + shift
            i_bot2 = this%parents(4) + shift
            if (i_top1 > size(mu)) i_top1 = i_top1 - size(mu)
            if (i_bot1 > size(mu)) i_bot1 = i_bot1 - size(mu)
            if (i_top2 > size(mu)) i_top2 = i_top2 - size(mu)
            if (i_bot2 > size(mu)) i_bot2 = i_bot2 - size(mu)
            mu_strengths = [mu(i_top1),mu(i_bot1),mu(i_top2),mu(i_bot2)]
            ! mu_strengths(i) = (mu(i_top1) - mu(i_bot1)) - (mu(i_top2) - mu(i_bot2))
        ! end do
    
    end function filament_segment_get_doublet_strengths


    subroutine filament_segment_calc_velocities(this, P, freestream, mirror_filament, sigma, mu, asym_flow, v_d)
        ! Calculates the velocity induced at the given point

        implicit none

        class(filament_segment),intent(in) :: this
        real,dimension(3),intent(in) :: P
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_filament, asym_flow
        real,dimension(:),allocatable,intent(in) :: sigma, mu
        real,dimension(3),intent(out) :: v_d

        real,dimension(:,:),allocatable :: doublet_inf
        real,dimension(:),allocatable :: doublet_strengths
    

        ! Get influences
        call this%calc_velocity_influences(P, freestream, mirror_filament, doublet_inf)

        ! Get strengths
        doublet_strengths = this%get_doublet_strengths(mu, mirror_filament, asym_flow)

        ! Apply strengths to calculate potentials
        ! v_s = matmul(source_inf, source_strengths) !!!! commenting out since we don't have sources

        ! if (this%in_wake) then !!!! commetning out since we are always in wake here
        ! v_d = matmul((doublet_inf(:,1:this%M_dim)+doublet_inf(:,this%M_dim+1:)), doublet_strengths)
        ! else
        v_d = matmul(doublet_inf, doublet_strengths)
        ! end if

    end subroutine filament_segment_calc_velocities



end module filament_segment_mod