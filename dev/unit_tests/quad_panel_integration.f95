program main

    implicit none

    integer :: i, j, k, k_next
    real,dimension(3) :: d_eta, d_xi
    real,dimension(:,:,:,:),allocatable :: II
    real,dimension(:,:),allocatable :: vertices_ls, C
    real,dimension(:),allocatable :: C_P_params
    integer :: Ni, Nj
    real :: A, C_P_avg

    ! Set up panel
    allocate(vertices_ls(2,3), source=0.)
    vertices_ls(1,2) = 1.
    vertices_ls(1,3) = 1.
    vertices_ls(2,3) = 1.
    A = 0.5

    ! Set up distribution
    allocate(C_P_params(6))
    C_P_params(1) = 2.
    C_P_params(2) = 1.
    C_P_params(3) = -5.
    C_P_params(4) = 2.
    C_P_params(5) = 3.
    C_P_params(6) = 1.

    ! Calculate C integrals

    ! Set how many C integrals we need
    Ni = 3
    Nj = 3

    ! Calculate offsets for each edge
    do k=1,3

        ! Get endpoint index
        k_next = modulo(k, 3) + 1

        ! Get offset
        d_xi(k) = vertices_ls(1,k_next) - vertices_ls(1,k)
        d_eta(k) = vertices_ls(2,k_next) - vertices_ls(2,k)

    end do

    ! Initialize H integrals (remember II(i,j,0,:) = H(i,j,:))
    allocate(II(0:Ni+2,0:Nj,0:Ni+1,3))
    do i=0,Ni+2
        II(i,0,0,:) = 1./(i+1.)
    end do

    ! Calculate other H integrals using eta recursion
    do j=1,Nj
        do i=0,Ni-j
            II(i,j,0,:) = vertices_ls(2,:)*II(i,j-1,0,:) + d_eta*II(i+1,j-1,0,:)
        end do 
    end do

    ! Calculate I integrals using xi recursion
    do j=0,Nj
        do k=1,Ni-j
            do i=Ni-j,k,-1
                II(i,j,k,:) = vertices_ls(1,:)*II(i-1,j,k-1,:) + d_xi*II(i,j,k-1,:)
            end do
        end do
    end do

    ! Calculate C integrals (remember G(i,j,:) = II(i,j,i,:))
    allocate(C(0:Ni,0:Nj))
    do i=0,Ni
        do j=0,Nj
            C(i,j) = sum(d_eta*II(i+1,j,i+1,:))/(i+1)
        end do
    end do

    ! Integrate
    C_P_avg = C(0,0)*C_P_params(1) + C(1,0)*C_P_params(2) + C(0,1)*C_P_params(3) + &
              0.5*C(2,0)*C_P_params(4) + C(1,1)*C_P_params(5) + 0.5*C(0,2)*C_P_params(6)

    write(*,*) C_P_avg

    !function panel_get_moment_about_centroid(this, mu, sigma, mirrored, N_body_panels, N_body_verts, asym_flow, &
    !                                               freestream, inner_flow, mirror_plane, rule, M_corr) result(C_M)
    !    ! Calculates the moment coefficient (units of length) about the panel centroid

    !    implicit none

    !    class(panel),intent(in) :: this
    !    real,dimension(:),allocatable,intent(in) :: mu, sigma
    !    logical,intent(in) :: mirrored, asym_flow
    !    integer,intent(in) :: N_body_panels, N_body_verts
    !    type(flow),intent(in) :: freestream
    !    real,dimension(3),intent(in) :: inner_flow
    !    integer,intent(in) :: mirror_plane
    !    character(len=*),intent(in) :: rule
    !    real,intent(in),optional :: M_corr

    !    real,dimension(3) :: C_M

    !    real,dimension(6) :: C_P_params

    !    ! Constant pressure
    !    if (this%order == 1) then

    !        C_M = 0.

    !    ! Quadratic pressure
    !    else

    !        ! Get pressure distribution parameters
    !        if (present(M_corr)) then
    !            C_P_params = this%get_quadratic_pressure_params(mu, sigma, mirrored, N_body_panels, N_body_verts, asym_flow, &
    !                                                            freestream, inner_flow, mirror_plane, rule, M_corr)
    !        else
    !            C_P_params = this%get_quadratic_pressure_params(mu, sigma, mirrored, N_body_panels, N_body_verts, asym_flow, &
    !                                                            freestream, inner_flow, mirror_plane, rule)
    !        end if

    !        ! Integrate
    !        if (mirrored) then

    !            ! Integrate
    !            C_M(1) = this%C_mir(1,0)*C_P_params(1) + this%C_mir(2,0)*C_P_params(2) + this%C_mir(1,1)*C_P_params(3) + &
    !                     0.5*this%C_mir(3,0)*C_P_params(4) + this%C_mir(2,1)*C_P_params(5) + 0.5*this%C_mir(1,2)*C_P_params(6)
    !            C_M(2) = this%C_mir(0,1)*C_P_params(1) + this%C_mir(1,1)*C_P_params(2) + this%C_mir(0,2)*C_P_params(3) + &
    !                     0.5*this%C_mir(2,1)*C_P_params(4) + this%C_mir(1,2)*C_P_params(5) + 0.5*this%C_mir(0,3)*C_P_params(6)
    !            C_M(3) = 0.

    !            ! Apply area factor and transform
    !            C_M = -this%J_mir*cross(this%n_g_mir, matmul(this%A_ls_to_g_mir, C_M))

    !        else

    !            ! Integrate
    !            C_M(1) = this%C(1,0)*C_P_params(1) + this%C(2,0)*C_P_params(2) + this%C(1,1)*C_P_params(3) + &
    !                     0.5*this%C(3,0)*C_P_params(4) + this%C(2,1)*C_P_params(5) + 0.5*this%C(1,2)*C_P_params(6)
    !            C_M(2) = this%C(0,1)*C_P_params(1) + this%C(1,1)*C_P_params(2) + this%C(0,2)*C_P_params(3) + &
    !                     0.5*this%C(2,1)*C_P_params(4) + this%C(1,2)*C_P_params(5) + 0.5*this%C(0,3)*C_P_params(6)
    !            C_M(3) = 0.

    !            ! Apply area factor and transform
    !            C_M = -this%J*cross(this%n_g, matmul(this%A_ls_to_g, C_M))

    !        end if

    !    end if
    !    
    !end function panel_get_moment_about_centroid

end program main