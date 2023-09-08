module filament_segment_mod
    use base_geom_mod
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
    
        logical :: first_in_dod = .true. !!!! the nomenclature can certainly change here
        logical,dimension(2) :: both_in_dod = .true. !!!! the nomenclature can certainly change here
        logical,dimension(2) :: neither_in_dod = .true. !!!! the nomenclature can certainly change here

    end type filament_dod


    type filament_segment 
        ! A section of a filament that contains two vertices. 
        type(vertex_pointer), dimension(:),allocatable :: vertices 
        integer :: index ! index of this filament in the filament mesh array
        integer :: N = 2 ! number of vertices
        integer :: 

        contains
            procedure :: init => filament_segment_init 
            procedure :: check_dod => filament_check_dod
            procedure :: calc_integrals => filament_segment_calc_integrals
            ! Getters
            procedure :: get_vertex_index => filament_segment_get_vertex_index
            !velocities
            procedure :: calc_velocity_influences => filament_segment_calc_velocity_influences
            procedure :: assemble_v_d_M_space => filament_segment_assemble_v_d_M_space 
            procedure :: calc_influence_integrals => filament_segment_calc_influence_integrals 
            procedure :: get_vertex_loc => filament_segment_get_vertex_loc

    end type filament_segment 


contains


    subroutine filament_segment_init(this,v1,v2,index)
        ! initializes important filament segment parameters 

        implicit none

        class(filament_segment), intent(inout) :: this
        type(vertex),intent(inout),target :: v1, v2
        integer, intent(in) :: index

        ! set number of sides 
        this%N = 2

        ! allocate vertex array
        allocate(this%vertices(this%N))

        ! assign pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2

        ! store index
        this%index = index

        ! still need to calculated derived geometry (based on what is needed in solver)
    end subroutine filament_segment_init

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


    subroutine filament_segment_calc_velocity_influences(this, freestream, mirror_filament, v_d_M_space)
        ! calculates the velocity influence of a vortex filament 
        
        implicit none 

        class(filament_segment),intent(in) :: this
        real,dimension(3),intent(in) :: F !!!! not sure what P is in panel_calc_velocity_influences
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_filament
        real, dimension(:,:),allocatable,intent(out) :: v_d_M_space

        type(filament_dod) :: dod_info !!!! might need to make a filament_dod type in this mod (filament_segment_mod)
        type(eval_point_geom) :: geom !!!! go to base_geom and make this type? Or come up with another way. 
        type(filament_segment_integrals) :: int 
        real,dimension(:,:),allocatable :: v_d_mu_space
        real :: x2, y2 !!!! not sure if this is needed
        integer :: i
        
        dod_info = this%check_dod(F, freestream, mirror_filament)
        if(dod_info%in_dod) then 
            if(freestream%supersonic) then 
                geom = this%calc_supersonic_geom(F, freestream, mirror_filament, dod_info)
            else
                geom = this%calc_subsonic_geom(F, freestream, mirror_filament)    
            end if 
             ! Get integrals 
            int = this%calc_integrals(geom, 'velocity', freestream, mirror_filament, dod_info)
            v_d_M_space =  this%assemble_v_d_M_space(int, geom, freestream, mirror_filament)
        end if 

    end subroutine filament_segment_calc_velocity_influences


    function filament_segment_check_dod(this, eval_point, freestream, mirror_filament) result(dod_info)
        ! Determines how (if) this panel lies within the domain of dependence of the evaluation point
        
        implicit none 

        class(filament_segment),intent(in) :: this
        real,dimension(3),intent(in) :: eval_point
        type(flow),intent(in) :: freestream 
        logical,intent(in),optional :: mirror_filament 

        type(filament_dod) :: dod_info

        integer :: i, i_next, inside_outside
        logical :: mirrored, in_filament

        if (present(mirror_filament)) then
            mirrored = mirror_filament
        else
            mirrored = .false.
        end if  

        ! check if the filament segment is entirely inside or ouside the DoD
        inside_outside = this%entirely_inside_outside_dod(eval_point, freestream, mirrored) !!!! make this function

        if (inside_outside == 1) then
            dod_info%in_dod = .true.
            dod_info%edges_in_dod = .true.

        else if (inside_outside == -1) then
            dod_info%in_dod = .false.
            dod_info%edges_in_dod = .false.

        ! Neither guaranteed, so we go further
        else
            
            ! Read in vertex information
            do i=1,this%N
                if (mirrored) then
                    these_verts_in_dod(i) = freestream%point_in_dod( &
                                                                    mirror_across_plane(this%get_vertex_loc(i), this%mirror_plane),&
                                                                    eval_point)
                else
                     these_verts_in_dod(i) = freestream%point_in_dod(this%get_vertex_loc(i), eval_point)
                end if
            end do
        
                   ! If all the vertices are in, then the panel is totally in and we can be done
            if (all(these_verts_in_dod)) then
                dod_info%in_dod = .true.
                dod_info%edges_in_dod = .true.
                
            end if 
        end if !!!! use Josh's stuff here - SA 
    end do 
        

    end function filament_segment_check_dod


    function filament_segment_calc_integrals(this, geom, influence_type, freestream, mirror_filament, dod_info) result(int)
        ! Calculates the integrals necessary for the given influence
        
        implicit none 
        
        class(filamenet_segment),intent(in) :: this 
        type(eval_point_geom),intent(in) :: geom
        ! character(len=*),intent(in) :: influence_type !!!! don't think this will apply - SA 
        type(flow), intent(in) :: freestream
        logical,intent(in) :: mirror_filament

        type(filament_segment_integrals) :: int
        type(filament_dod),intent(in) :: dod_info

        int%s = freestream%s 
        int%rs = int%r*int%s 

        ! Calculate necessary integrals based on the flow condition
        call this%calc_influence_integrals(geom, dod_info, freestream, mirror_filament, int) !!!! I think the int input in this subroutine call is what will output from this function
        !!!! should dod_info carry whether the dod is supersonic or subsonic? -SA 

 
    end function filament_segment_calc_integrals 


    subroutine filament_segment_calc_influence_integrals(this, geom, dod_info, freestream, mirror_filament, int)
        ! calculates the supersonic and subsonic vortex filament influence integrals

        implicit none 

        class(filament_segment),intent(in) :: this
        type(eval_point_geom),intent(in) :: geom 
        type(dod),intent(in) :: dod_info
        type(flow),intent(in) :: freestream 
        logical,intent(in) :: mirror_filament
        type(filament_segment_integrals),intent(inout) :: int
        real :: unrel_v_numer_one, unrel_v_numer_two
        real :: unrel_v_denom_one, unrel_v_denom_two
        real :: unrel_w_numer_one, unrel_w_numer_two
        real :: unrel_w_denom_one, unrel_w_denom_two
        real :: xf, xi, yf, yi, zf, zi
        real :: xo, yo, zo 
        real,dimension(3) :: loc_1, loc_2

        loc_1 = this%get_vertex_loc(1) !!!! put in appropriate input and get 
        loc_2 = this%get_vertex_loc(2) !!!! put in appropriate input and get 

        xf = loc_2(1)
        xi = loc_1(1)
        
        yf = loc_2(2)
        yi = loc_1(2)

        zf = loc_2(3)
        zi = loc_1(3)

        xo = geom(1)
        yo = geom(2)
        zo = geom(3)

        ! write out the integrals here based on the different dod stuff. 
        unrel_v_numer_one = (zo-zf)*(xo-xf)
        unrel_v_denom_one = ((yo-yf)**2+(zo-zf)**2)*(((xo-xf)**2+bsq*((yo-yf)**2+(zo-zf)**2))**(1/2))
        unrel_v_numer_two = (zo-zi)*(xo-xi)
        unrel_v_denom_two = ((yo-yi)**2+(zo-zi)**2)*(((xo-xi)**2+bsq*((yo-yi)**2+(zo-zi)**2))**(1/2))
        
        unrel_w_numer_one = (yo-yf)*(xo-xf)
        unrel_w_denom_one = ((yo-yf)**2+(zo-zf)**2)*(((xo-xf)**2+bsq*((yo-yf)**2+(zo-zf)**2))**(1/2))
        unrel_w_numer_two = (yo-yi)*(xo-xi)
        unrel_w_denom_two = ((yo-yi)**2+(zo-zi)**2)*(((xo-xi)**2+bsq*((yo-yi)**2+(zo-zi)**2))**(1/2))


       !if (this%relaxed) !!!! put in logic like this later when we employ relaxation. -SA  
        if (dod_info%both_in_dod) then 
            int%u = 0
            int%v = (1/(2*pi*k))*((unrel_v_numer_one/unrel_v_denom_one)-(unrel_v_numer_two/unrel_v_denom_two)) !!!! declare all of these variables and stuff - SA
            int%w = (-1/(2*pi*k))*((unrel_w_numer_one/unrel_w_denom_one)-(unrel_w_numer_two/unrel_w_denom_two)) !!!! declare all of these variables and stuff - SA
        else if (dod_info%first_in_dod) then 
            int%u = 0
            int%v = (-1/(2*pi*k))*((unrel_v_numer_two/unrel_v_denom_two)) !!!! declare all of these variables and stuff - SA
            int%w = (1/(2*pi*k))*((unrel_w_numer_two/unrel_w_denom_two)) !!!! declare all of these variables and stuff - SA
        else if (dod_info%neither_in_dod) then
            int%u = 0
            int%v = 0
            int%w = 0
        else 
            int%u = 0
            int%v = 0
            int%w = 0
        end if 
    
    end subroutine filament_segment_calc_influence_integrals


    function filament_segment_assemble_v_d_M_space(this, int, geom, freestream, mirror_panel) result(v_d_M_space) !!!! need to create this space for wakes
        ! Assembles the doublet-induced velocity influence coefficient matrix from the previously-calculated influence integrals

        implicit none

        class(filamenet_segment),intent(in) :: this
        type(filament_segment_integrals),intent(in) :: int
        type(eval_point_geom),intent(in) :: geom
        type(flow),intent(in) :: freestream
        logical,intent(in) :: mirror_filament

        real,dimension(:,:),allocatable :: v_d_M_space

        real,dimension(:,:),allocatable :: v_d_mu_space

        ! Allocate space
        allocate(v_d_mu_space(3,this%mu_dim), source=0.)
        if (this%in_wake) then
            allocate(v_d_M_space(3,2*this%M_dim), source=0.)
        else
            allocate(v_d_M_space(3,this%M_dim), source=0.)
        end if

        v_d_mu_space(1,1) = 0
        v_d_mu_space(1,2) = int%hH113
        v_d_mu_space(1,3) = 0

        v_d_mu_space(2,1) = 0
        v_d_mu_space(2,2) = 0
        v_d_mu_space(2,3) = int%hH113

        v_d_mu_space(3,1) = 0
        v_d_mu_space(3,2) = int%H213
        v_d_mu_space(3,3) = int%H123

        if (this%order == 2) then

            ! Quadratic terms
            v_d_mu_space(1,4) = 3.*int%r*(0.5*int%H215*(geom%P_ls(1)**2)*geom%h + int%hH315*geom%P_ls(1) + 0.5*int%H415*geom%h)
            v_d_mu_space(1,5) = 3.*int%r*(int%H215*geom%h*geom%P_ls(1)*geom%P_ls(2) + int%hH315*geom%P_ls(2) &
                                + int%H225*geom%h*geom%P_ls(1) + int%H325*geom%h)
            v_d_mu_space(1,6) = 3.*int%r*(0.5*int%H215*(geom%P_ls(2)**2)*geom%h + int%H225*geom%h*geom%P_ls(2) &
                                + 0.5*int%H235*geom%h)

            v_d_mu_space(2,4) = 3.*int%s*(0.5*int%H125*(geom%P_ls(1)**2)*geom%h + int%H225*geom%h*geom%P_ls(1) &
                                + 0.5*int%H325*geom%h)
            v_d_mu_space(2,5) = 3.*int%s*(int%H125*geom%h*geom%P_ls(1)*geom%P_ls(2) + int%hH135*geom%P_ls(1) &
                                + int%H225*geom%h*geom%P_ls(2) + int%H235*geom%h)
            v_d_mu_space(2,6) = 3.*int%s*(0.5*int%H125*(geom%P_ls(2)**2)*geom%h + int%hH135*geom%P_ls(2) + 0.5*int%H145*geom%h)

            v_d_mu_space(3,4) = 0.5*(geom%P_ls(1)**2)*int%H113_3rsh2H115 + geom%P_ls(1)*(int%H213 - 3.*int%rs*geom%h2*int%H215) &
                                + 0.5*(int%H313 - 3.*int%rs*geom%h*int%hH315)
            v_d_mu_space(3,5) = geom%P_ls(1)*geom%P_ls(2)*(int%H113_3rsh2H115) &
                                + geom%P_ls(2)*(int%H213 - 3.*int%rs*geom%h2*int%H215) &
                                + geom%P_ls(1)*(int%H123 - 3.*int%rs*geom%h2*int%H125) + int%H223 - 3.*int%rs*geom%h2*int%H225
            v_d_mu_space(3,6) = 0.5*(geom%P_ls(2)**2)*int%H113_3rsh2H115 + geom%P_ls(2)*(int%H123 - 3.*int%rs*geom%h2*int%H125) &
                                + 0.5*(int%H133 - 3.*int%rs*geom%h*int%hH135)

        end if 
            
        ! Convert to strength influences (Davis Eq. (4.41))
        if (mirror_panel) then
            v_d_M_space(:,1:this%M_dim) = int%s*freestream%K_inv*matmul(v_d_mu_space, this%T_mu_mir)
        else
            v_d_M_space(:,1:this%M_dim) = int%s*freestream%K_inv*matmul(v_d_mu_space, this%T_mu)
        end if

        ! Wake bottom influence is opposite the top influence
        if (this%in_wake) then
            v_d_M_space(:,this%M_dim+1:this%M_dim*2) = -v_d_M_space(:,1:this%M_dim)
        end if

        ! Transform to global coordinates
        if (mirror_panel) then
            v_d_M_space = matmul(transpose(this%A_g_to_ls_mir), v_d_M_space)
        else
            v_d_M_space = matmul(transpose(this%A_g_to_ls), v_d_M_space)
        end if
        
    end function filament_segment_assemble_v_d_M_space


end module filament_segment_mod