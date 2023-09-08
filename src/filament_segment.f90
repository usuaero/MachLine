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


    subroutine filament_segment_calc_velocity_influences(this, F, freestream, mirror_filament, v_d_M_space)
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

        call this%get_vertex_location() !!!! put in appropriate input and get 
        
        ! write out the integrals here based on the different dod stuff. 

       !if (this%relaxed) !!!! put in logic like this later when we employ relaxation. -SA  
        if (dod_info%both_in_dod) then 
            int%u = 0
            int%v = (1/(2*pi*k))*((((zo-zf)*(xo-xf))/(()()))) !!!! declare all of these variables and stuff - SA
        
        end if 
    
    end subroutine filament_segment_calc_influence_integrals




end module filament_segment_mod