module Multigrid_2D_Dirichlet_BC
  use Exact_Solver
  implicit none
contains
  !---------------------------------------------------------------------------!
  ! Multigrid V cycle scheme
  !---------------------------------------------------------------------------!
  subroutine MG_Vcycle(Nx,Ny,dx,dy,RHS,U,Level_num,tol,exact_solver,I_cycle,file_name)
    implicit none
    integer,intent(in) :: Nx,Ny,Level_num,exact_solver
    real*8 ,intent(in) :: dx,dy,tol
    real*8,dimension(0:Nx,0:Ny), intent(in)     :: RHS
    real*8,dimension(0:Nx,0:Ny), intent(inout)  :: U
    integer,dimension(Level_num) :: Level_Nx,Level_Ny
    real*8 ,dimension(Level_num) :: Level_dx,Level_dy
    real*8  :: rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,j,k,Level,iter
    integer, intent(inout) :: I_cycle
    character*80 , intent(in) :: file_name

    ! Defined all data space
    type space
      real*8,allocatable :: data(:,:)
    end type space

    type (space) UL(Level_num),R(Level_num),F(Level_num),P(Level_num)

    allocate(UL(1)%data(0:Nx,0:Ny),R(1)%data(0:Nx,0:Ny),F(1)%data(0:Nx,0:Ny))
    allocate(P(1)%data(0:Nx,0:Ny))

    Level_Nx(1) = Nx
    Level_Ny(1) = Ny
    Level_dx(1) = dx
    Level_dy(1) = dy

    do i=2,Level_num
      k=2**(i-1)
      Level_Nx(i) = Nx / k
      Level_Ny(i) = Ny / k
      Level_dx(i) = dx * dble(k)
      Level_dy(i) = dy * dble(k)
      allocate(UL(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( R(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( F(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( P(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
    enddo

    Max_Iter   = 100000     ! Allowed maximum number of outer iteration
    Relax_Iter = 2          ! Number of relaxation for restriction in V-cycle

    !Check the coarsest grid
    if (Level_Nx(Level_num).le.3) then
    write(*,*) Level_num," level is high for ",Nx," grid "
    stop
    end if

    do j=0,Level_Ny(1)
    do i=0,Level_Nx(1)
      F(1)%data(i,j)  = RHS(i,j)
      UL(1)%data(i,j) = U(i,j)
      R(1)%data(i,j)  = 0.0d0
    enddo
    enddo

    !Compute initial resitual:
    call Residual(Nx,Ny,dx,dy,F(1)%data,UL(1)%data,R(1)%data)
    call L2norm(Nx,Ny,R(1)%data,rms0)

    !open(40,file='residual_All_MG4V3.plt')
    !write(40,*) 'variables ="Cycle","Iteration","rms","Residual"'

    iter=0
    open  (unit = 66, file=file_name) 
    write(66,*) 'Cycle     ', ' Residual'
    do I_cycle=1,Max_Iter

      call Residual(Level_Nx(1),Level_Ny(1),Level_dx(1),Level_dy(1),F(1)%data,UL(1)%data,R(1)%data)

        ! Check for convergence on finest grid
      call L2norm(Level_Nx(1),Level_Ny(1),R(1)%data,rms)
      write(66,*) I_cycle,rms/rms0

    !---------------------------------------------------------------------------!
      do Level=1,Level_num-1

        ! Relax
        do i=1,Relax_Iter
          call Relax(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data)
        end do

        ! Compute residual
        call Residual(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data,R(Level)%data)

        ! Check for convergence on finest grid
        call L2norm(Level_Nx(Level),Level_Ny(Level),R(Level)%data,rms)
        !write(40,*) I_cycle,iter,rms,rms/rms0

        if (rms/rms0.le.tol .and. Level .eq. 1) goto 10

        ! Restriction
        call Restriction(Level_Nx(Level),Level_Ny(Level),Level_Nx(Level+1),Level_Ny(Level+1),R(Level)%data,F(Level+1)%data)

        do j=0,Level_Ny(Level+1)
        do i=0,Level_Nx(Level+1)
          UL(Level+1)%data(i,j) = 0.0d0
        end do
        end do

      end do
      !---------------------------------------------------------------------------!
      ! Compute residual on coarsest grid
      call Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
      call L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rmsc)
      !iter=iter+1
      !write(*,*) I_cycle,iter,rmsc,rmsc/rms0

      ! Solve exact solution on coarsest grid
      if (exact_solver .eq. 1) then
        do while (rms/rmsc .gt. tol)
          call Relax(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
          call Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
          call L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
        end do
      else if (exact_solver .eq. 2) then
        call LU_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
      else
        call CG_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
      endif

    !---------------------------------------------------------------------------!

      do Level=Level_num-1,1,-1

        ! Prolongation
        call Prolongation(Level_Nx(Level+1),Level_Ny(Level+1),Level_Nx(Level),Level_Ny(Level),UL(Level+1)%data,P(Level)%data)

        ! Correct
        do j=1,Level_Ny(Level)-1
        do i=1,Level_Nx(Level)-1
          UL(Level)%data(i,j) = UL(Level)%data(i,j) + P(Level)%data(i,j)
        end do
        end do

        ! Relax
        do i=1,Relax_Iter
          call Relax(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data)
        end do

      end do


    end do  ! Outer iteration loop

    10 continue


    do j=0,Ny
    do i=0,Nx
      U(i,j) = UL(1)%data(i,j)
    end do
    end do

    close(66)

    do i=1,Level_num
      deallocate(UL(i)%data,R(i)%data,F(i)%data,P(i)%data)
    enddo

    return
  end subroutine
  
  !---------------------------------------------------------------------------!
  ! Relaxation formula for Poisson equation
  ! Uses GS relaxation
  !---------------------------------------------------------------------------!
  subroutine Relax(Nx,Ny,dx,dy,F,U)
    implicit none
    integer ,intent(in) :: Nx,Ny
    real*8  ,intent(in) :: dx,dy
    real*8  ,dimension(0:Nx,0:Ny),intent(in)    :: F
    real*8  ,dimension(0:Nx,0:Ny),intent(inout) :: U

    real*8  :: a
    integer :: i,j

    a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)

    do j=1,Ny-1
    do i=1,Nx-1
      U(i,j) = (1.0d0/a)*(F(i,j) &
              &  - (U(i+1,j)+U(i-1,j))/(dx*dx) &
              &  - (U(i,j+1)+U(i,j-1))/(dy*dy) )
    end do
    end do

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Residual formula for Poisson equation
  !---------------------------------------------------------------------------!
  subroutine Residual(Nx,Ny,dx,dy,F,U,R)
    implicit none
    integer,intent(in) :: Nx,Ny
    real*8 ,intent(in) :: dx,dy
    real*8, dimension(0:Nx,0:Ny),intent(in)  :: U,F
    real*8, dimension(0:Nx,0:Ny),intent(out) :: R

    integer :: i,j

    do j=1,Ny-1
    do i=1,Nx-1
      R(i,j) = F(i,j) - (U(i+1,j) - 2.0d0*U(i,j) + U(i-1,j))/(dx*dx) &
                   &  - (U(i,j+1) - 2.0d0*U(i,j) + U(i,j-1))/(dy*dy)
    end do
    end do

    !Boundary conditions for residuals
    do i=0,Nx
      R(i,0)  = 0.0d0
      R(i,Ny) = 0.0d0
    end do

    do j=0,Ny
      R(0,j)  = 0.0d0
      R(Nx,j) = 0.0d0
    end do

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Restriction operators
  !---------------------------------------------------------------------------!
  subroutine Restriction(Nxf,Nyf,Nxh,Nyh,R,F)
    implicit none
    integer ,intent(in) :: Nxf,Nyf,Nxh,Nyh
    real*8, dimension(0:Nxf,0:Nyf),intent(in)  :: R       !on higher grid
    real*8, dimension(0:Nxh,0:Nyh),intent(out) :: F       !on lower grid
    integer :: i,j


    do j=1,Nyh-1
    do i=1,Nxh-1
      F(i,j) = 1.0d0/16.0d0 * (4.0d0*R(2*i,2*j) &
            &      + 2.0d0 * (R(2*i+1,2*j)+R(2*i-1,2*j)+R(2*i,2*j+1)+R(2*i,2*j-1))        &
            &      + 1.0d0 * (R(2*i+1,2*j+1)+R(2*i-1,2*j-1)+R(2*i-1,2*j+1)+R(2*i+1,2*j-1)))
    end do
    end do

    !Boundary conditions
    do i=0,Nxh
      F(i,0)   = R(2*i,0)
      F(i,Nyh) = R(2*i,Nyf)
    end do

    do j=0,Nyh
      F(0,j)   = R(0,2*j)
      F(Nxh,j) = R(Nxf,2*j)
    end do

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Prolongation operator
  !---------------------------------------------------------------------------!
  subroutine Prolongation(Nxh,Nyh,Nxf,Nyf,U,P)
    implicit none
    integer ,intent(in) :: Nxf,Nyf,Nxh,Nyh
    real*8, dimension(0:Nxh,0:Nyh) ,intent(in)  :: U       !on lower grid
    real*8, dimension(0:Nxf,0:Nyf) ,intent(out) :: P       !on higher grid
    integer :: i,j


    do j=0,Nyh-1
    do i=0,Nxh-1
      P(2*i,2*j)     = U(i,j)
      P(2*i+1,2*j)   = 1.0d0/2.0d0*(U(i,j)+U(i+1,j))
      P(2*i,2*j+1)   = 1.0d0/2.0d0*(U(i,j)+U(i,j+1))
      P(2*i+1,2*j+1) = 1.0d0/4.0d0*(U(i,j)+U(i,j+1)+U(i+1,j)+U(i+1,j+1))
    end do
    end do

    !Boundary conditions
    do j=0,Nyh
      P(Nxf,2*j)     = U(Nxh,j)
    end do
    do i=0,Nxh
      P(2*i,Nyf)     = U(i,Nyh)
    end do

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Compute L2-norm for an array
  !---------------------------------------------------------------------------!
  subroutine L2norm(Nx,Ny,R,RMS)
    implicit none
    integer ,intent(in) :: Nx,Ny
    real*8, dimension(0:Nx,0:Ny),intent(in) :: R
    real*8  ,intent(out):: RMS
    integer :: i,j

    RMS=0.0d0

    do j=1,Ny-1
    do i=1,Nx-1
      RMS = RMS + R(i,j)*R(i,j)
    end do
    end do

    RMS= dsqrt(RMS/dfloat((Nx-1)*(Ny-1)))

    return
  end subroutine

end Module
