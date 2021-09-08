module vtk

    use geometry

    implicit none
    
contains

    subroutine load_vtk(mesh_file, N_vert, N_panel, vertices, panels)

        implicit none

        character(len=:),allocatable,intent(in) :: mesh_file
        integer,intent(out) :: N_vert, N_panel
        character(len=:),allocatable :: dummy_read
        real,dimension(:,:),allocatable,intent(out) :: vertices
        type(panel),dimension(:),allocatable,intent(out) :: panels
        integer :: i, N, i1, i2, i3

        ! Open file
        open(1, file=mesh_file)

            ! Determine number of vertices
            read(1,*) ! Header
            read(1,*) ! Header
            read(1,*) ! Header
            read(1,*) ! Header
            read(1,*) dummy_read, N_vert, dummy_read

            ! Allocate vertex array
            allocate(vertices(N_vert,3))

            ! Store vertices
            do i=1,N_vert
                read(1,*) vertices(i,1), vertices(i,2), vertices(i,3)
            end do

            ! Determine number of panels
            read(1,*) dummy_read, N_panel, dummy_read

            ! Allocate panel array
            allocate(panels(N_panel))

            ! Initialize panels
            do i=1,N_panel
                read(1,*) N, i1, i2, i3
                panels(i) = panel(N, i1, i2, i3, vertices(i1,:), vertices(i2,:), vertices(i3,:))
            end do

        close(1)
    
    end subroutine load_vtk

    
end module vtk