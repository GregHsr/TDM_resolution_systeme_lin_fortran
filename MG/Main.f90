!-----------------------------------------------------------------------------!
!  Multigrid method for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!-----------------------------------------------------------------------------!
program Poisson_2D
  use Multigrid_2D_Dirichlet_BC
  implicit none
  integer::i,j,k,m,nx,ny,exactsolver,Level_num,I_cycle,nn
  character(len=7) :: mode
  character(len=1) :: Level,version
  real*8,dimension(:,:),allocatable ::u,f,ue,e
  real*8,dimension(:),allocatable ::x,y
  real*8 ::dx,dy,tol,rms,x0,xL,y0,yL,start,finish
  character*80 :: file_name,taille,file_cpu

  !Domain
  x0 =-1.0d0 !left
  xL = 1.0d0 !right

  y0 =-1.0d0 !bottom
  yL = 1.0d0 !up

 
!reading the data
open(unit=10,file='input_mg.in')
     read(10,*) nx
     write(*,*) 'nx=',nx
     read(10,*) ny
     write(*,*) 'ny=',ny
     read(10,*) tol
     write(*,*) 'tol',tol
     read(10,*) exactsolver
     write(*,*) 'exactsolver=',exactsolver
close(10)
write(taille,'(I4.4)') nx

       file_name = 'Residu_Multigrille_'//taille//'.dat'
       file_name=trim(file_name)//'.dat'
       file_cpu='Temps_CPU_Multigrille_'//taille
       file_cpu=trim(file_cpu)//'.dat'       


  !---------------------------------------------------------------------------!
  ! Exact Solver
  !     1. Gauss Seidel Solver
  !     2. LU decomposition
  !     3. CG Solver 
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!

  

  Level_num = 0
   
  k = 0
  nn = nx
  do while (nn > 4 )
     nn = nn / 2
     Level_num = Level_num +1
  enddo
  
  write(*,*) 'Level_num =', Level_num
  !do k=8,11

    ! Multigrid Level
    !do Level_num=6,7

    !number of points
    !nx = 2**k !number of grid points in x
    !ny = nx   !number of grid points in y

    !grid spacing (spatial)
    dx = (xL-x0)/dfloat(nx)
    dy = (yL-y0)/dfloat(ny)

    !spatial coordinates
    allocate(x(0:nx))
    do i=0,nx
    x(i) = x0 + dfloat(i)*dx
    end do

    allocate(y(0:ny))
    do j=0,ny
    y(j) = y0 + dfloat(j)*dy
    end do

    allocate(u(0:nx,0:ny))
    allocate(f(0:nx,0:ny))
    allocate(e(0:nx,0:ny))
    allocate(ue(0:nx,0:ny))

    !---------------------------------------------!
    ! Exact solution
    !---------------------------------------------!
    do j=0,ny
    do i=0,nx
    !f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
    !ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
    f(i,j) = ((y(j)**2-1.)**2)*(12*x(i)**2-4.)+ ((x(i)**2-1.)**2)*(12*y(j)**2-4.)
    ue(i,j)= ((x(i)*x(i)-1.0d0)**2)*((y(j)*y(j)-1.0d0)**2)
   end do
    end do


    !Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do

    !Boundary conditions has to satisfy exact solution
    do i=0,nx
    u(i,0)  = ue(i,0)
    u(i,ny) = ue(i,ny)
    end do

    do j=0,ny
    u(0,j)  = ue(0,j)
    u(nx,j) = ue(nx,j)
    end do

    !---------------------------------------------------------------------------!
    ! Solver:
    !---------------------------------------------------------------------------!

  
    ! Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do
    call cpu_time(start)
    call MG_Vcycle(Nx,Ny,dx,dy,F,U,Level_num,tol,exactsolver,I_cycle,file_name)
    call cpu_time(finish)

    !----------------------!
    ! Error analysis:
    !----------------------!
    do j=0,nx
    do i=0,ny
    e(i,j) = dabs(u(i,j)-ue(i,j))
    end do
    end do

    !L-2 Norm:
    call L2norm(nx,ny,e,rms)

    !write(*,'(A7,i8,i8,E18.7,E15.7)') mode,nx,I_cycle,finish-start,rms
    !write(66,'(A7,i8,i8,E18.12,E15.10)') mode,nx,I_cycle,finish-start,rms

    open(unit=20,file=file_cpu)
    write(20,*) 'temps cpu=', finish-start
    close(20)

  close(66)


end program
