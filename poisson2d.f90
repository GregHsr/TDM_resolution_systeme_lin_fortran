
!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Iterative Methods for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!     
!     Point Jacobi (PJ)
!     Gauss-Seidel (GS)
!     Successive over relaxation (SOR)
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 22, 2015
!-----------------------------------------------------------------------------!

program poisson2d
implicit none
integer::i,j,k,nitr,nn,nx,ny,ichoice
real*8 ::dx,dy,x0,xL,y0,yL,omega,tol,err,start,finish,pi
real*8 ,allocatable ::f(:,:),u(:,:),ue(:,:),x(:),y(:),uu(:,:),uue(:,:)
real*8 ,allocatable :: xEnsight(:),yEnsight(:)


 integer :: mdim=3
 real*8 ,allocatable :: coef(:,:) 
 real*8 ,allocatable :: b(:),p(:),r(:),r2(:)
 real*8 ,allocatable :: q(:),s(:),xx(:)
 real*8 ,allocatable :: l(:,:)
 integer, allocatable:: jcoef(:)
 integer, dimension(1:5):: ldiag


character*80 :: file_name,taille,file_cpu

pi =  4.d0*atan(1.d0)
write(*,*) 'pi =',pi
!Domain
x0 = -1.0d0 !left
xL = 1.0d0 !right

y0 = -1.0d0 !bottom
yL = 1.0d0 !up

!reading the data
open(unit=10,file='input.in')
     read(10,*) nx
     read(10,*) ny
     read(10,*) ichoice
     read(10,*) tol
close(10)
write(taille,'(I4.4)') nx

start=0.d0 ; finish = 0.d0

if    (ichoice == 1 ) then
       file_name = 'Residu_Jacobi_'//taille
       file_name=trim(file_name)//'.dat'
       file_cpu='Temps_CPU_Jacobi_'//taille
       file_cpu=trim(file_cpu)//'.dat'
elseif(ichoice == 2 ) then
       file_name = 'Residu_Gauss_Seidel_'//taille//'.dat'
       file_name=trim(file_name)//'.dat'
       file_cpu='Temps_CPU_Gauss_Seidel_'//taille
       file_cpu=trim(file_cpu)//'.dat'
elseif(ichoice == 3 ) then
       file_name = 'Residu_SOR_'//taille//'.dat'
       file_name=trim(file_name)//'.dat'
       file_cpu='Temps_CPU_SOR_'//taille
       file_cpu=trim(file_cpu)//'.dat'       
elseif(ichoice == 4 ) then
       file_name = 'Residu_Gradient_conjugue_'//taille//'.dat'
       file_name=trim(file_name)//'.dat'
       file_cpu='Temps_CPU_Gradient_conjugue_'//taille
       file_cpu=trim(file_cpu)//'.dat'
elseif(ichoice == 5 ) then
       file_name = 'Residu_Gradient_conjugue_precondionne_'//taille//'.dat'
       file_name=trim(file_name)//'.dat'
       file_cpu='Temps_CPU_Gradient_conjugue_precondionne_'//taille
       file_cpu=trim(file_cpu)//'.dat'       
endif

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

!given source term and the exact solution on unit square domain
allocate(f(0:nx,0:ny))
allocate(ue(0:nx,0:ny))

f  = 0.d0
ue = 0.d0


do j=0,ny
do i=0,nx
 !   f(i,j) =-2.0d0*(1.0d0-x(i)*x(i)-y(j)*y(j))
 !   ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
     f(i,j) =-2.*pi*pi*(sin(pi*x(i))*sin(pi*y(j)))
     ue(i,j)= sin(pi*x(i))*sin(pi*y(j))
end do
end do


!do j=1,ny-1
!    do i=1,nx -1
!        f(i,j) = 0.d0
!        ue(i,j)= 1.d0
!    end do
!end do
!  i = 1   ; j = 1; f(i,j) = -2./dx/dx ; i = 1     ; j = ny-1; f(i,j) = -2./dx/dx 
!  i= nx-1 ; j = 1; f(i,j) = -2./dx/dx ; i = nx -1 ; j = ny-1; f(i,j) = -2./dx/dx

!  do j=2, ny-2
!     f(1,j)    = -1./dx/dx
!     f(nx-1,j) = -1./dx/dx
!  enddo

!  do i=2, nx-2
!    f(i,1)    = -1./dx/dx
!    f(i,ny-1) = -1./dx/dx
! enddo


open(19,file=file_name)
write(19,*) '#"iteration"   ,"Residual",  "Error"'

!Iterative schemes to solve Poisson Equation:
allocate(u(0:nx,0:ny))

!initial guess which satisfies boundary condtions
!Homogeneous Drichlet B.C.
!we are not updating boundary points, they are zero
do j=0,ny
do i=0,nx
u(i,j) = 0.0d0
end do
end do

call CPU_TIME(start)

if    (ichoice == 1 ) then
       !Point Jacobi (PJ)

   call PJ(nx,ny,dx,dy,f,ue,u,tol,nitr,x,y)
   
elseif(ichoice == 2 ) then
       !Gauss-Seidel (GS)

   call GS(nx,ny,dx,dy,f,ue,u,tol,nitr,x,y) 
    
elseif(ichoice == 3 ) then
       !Successive over relaxation (SOR)

   !omega = 2./(1. + sin(pi/nx))
   omega = 1.
   write(*,*) ' omega Opt = ' , omega
   call SOR(nx,ny,dx,dy,f,ue,u,omega,tol,nitr,x,y)
   
elseif(ichoice == 4 ) then
       !Conjugate Gradient (CG)

   allocate(coef((nx-1)*(ny-1),mdim)) ; allocate(jcoef(mdim)) ; allocate (xx((nx-1)*(ny-1)))
   allocate(r((nx-1)*(ny-1))) ; allocate(b((nx-1)*(ny-1))) ; allocate(p((nx-1)*(ny-1)))
   call matgen(coef,jcoef,nx-1,ny-1,mdim,dx,dy)
   
    
   k=1
   do j=1,ny-1
      do i=1,nx-1
         b(k)=-f(i,j)
         xx(k)=u(i,j)
         p(k)=0.0d0
         r(k)=0.0d0
         k=k+1
      enddo
   enddo
   
   nn=(nx-1)*(ny-1)
   
   call CG(coef,jcoef,b,xx,nn,mdim,p,r,tol,nitr,ue,nx,ny,x,y)
   
   k=1
   do j=1,ny-1
      do i=1,nx-1
         u(i,j)=xx(k)
         k=k+1
      enddo
   enddo
 
elseif(ichoice == 5 ) then
   !Preconditionned Conjugate Gradient
   allocate(coef((nx-1)*(ny-1),mdim)) ; allocate(l((nx-1)*(ny-1),5))
   allocate(jcoef(mdim)) ; allocate (xx((nx-1)*(ny-1)))
   allocate(r ((nx-1)*(ny-1))) ; allocate(b((nx-1)*(ny-1))) ; allocate(p((nx-1)*(ny-1)))
   allocate(r2((nx-1)*(ny-1))) ; allocate(q((nx-1)*(ny-1))) ; allocate(s((nx-1)*(ny-1)))

   call matgen(coef,jcoef,nx-1,ny-1,mdim,dx,dy)
   
    
   k=1
   do j=1,ny-1
      do i=1,nx-1
         b(k) =-f(i,j)
         xx(k)=u(i,j)
         p(k) =0.0d0
         r(k) =0.0d0
         r2(k)=0.0d0
         q(k) =0.0d0
         s(k) =0.d0
         k=k+1
      enddo
   enddo
 
   
   nn=(nx-1)*(ny-1)
   
   call ICCG2(coef,jcoef,l,ldiag,b,xx,nn,mdim,p,r,r2,q,s,tol,nitr,ue,nx,ny)
        
   k=1
   do j=1,ny-1
      do i=1,nx-1
         u(i,j)=xx(k)
         k=k+1
      enddo
   enddo
      
!elseif(ichoice == 6 ) then
!call multigrid(nx,ny,dx,dy,ue,u,tol,nitr)    !geometric multigrid
endif

close(19)

allocate (uu(1:nx-1,1:ny-1)) ; allocate (uue(1:nx-1,1:ny-1))

do j =1 , ny-1
    do i = 1 , nx-1
    uu(i,j) = u(i,j)
    uue(i,j) = ue(i,j)
    enddo
enddo

do j = 1, ny-1
 !   write(*,'(5(Xe12.5))') (uu(i,j), i=1,nx-1)
enddo
!write(*,*) '****************************************'
do j = 1, ny-1
!   write(*,'(5(Xe12.5))') (uue(i,j), i=1,nx-1)
enddo

allocate (xEnsight(1:nx-1)) ; allocate(yEnsight(1:ny-1))
xEnsight(1:nx-1) = x(1:nx-1) ; yEnsight(1:ny-1) = y(1:ny-1)
call write_result_ensight(xEnsight,yEnsight,u(1:nx-1,1:ny-1),ue(1:nx-1,1:ny-1),nx-1,ny-1,1)

call CPU_TIME(finish)

open(unit=20,file=file_cpu)
write(20,*) 'temps cpu=', finish-start
close(20)






end

!-----------------------------------------------------------------------------!
!Point Jacobi (PJ)
!-----------------------------------------------------------------------------!
subroutine PJ(nx,ny,dx,dy,f,ue,u,tol,nitr,x,y)  
implicit none
integer::nx,ny,i,j,nitr
real*8 ::dx,dy,tol,err,err0,a,err_f,err_f0
real*8 ::f(0:nx,0:ny),ue(0:nx,0:ny),u(0:nx,0:ny),v(0:nx,0:ny),e(0:nx,0:ny)
real*8 :: x(0:nx)
real*8 :: y(0:ny)

a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)

    nitr = 0
    err0=0.d0
    err_f0=0.d0
 
    !compute initial residual
 	do j=1,ny-1
	do i=1,nx-1
	err0 = err0+f(i,j)*f(i,j)
	end do
	end do   
    !max value of the error
    err0 = dsqrt(err0)

    !compute initial error
 	do j=1,ny-1
        do i=1,nx-1
        err_f0 = err_f0+ue(i,j)*ue(i,j)
        end do
    end do   
    !max value of the error
    err_f0 = dsqrt(err_f0)

    write(19,*) nitr, err0/err0, err_f0/err_f0
    write(*,*)  nitr, err0/err0, err_f0/err_f0
    

    err=err0
    err_f=err_f0
    
do while(err.ge.tol)

	nitr = nitr + 1
  
    !known iteration
	do j=0,ny
	do i=0,nx
	v(i,j) = u(i,j)
	end do
	end do

    !update
	do j=1,ny-1
	do i=1,nx-1
	u(i,j) = (1.0d0/a)*(f(i,j) &
                   - (v(i+1,j)+v(i-1,j))/(dx*dx) &
                   - (v(i,j+1)+v(i,j-1))/(dy*dy) )
	end do
	end do

    err=0.d0
    !compute residual
	do j=1,ny-1
	    do i=1,nx-1
	    err = err+(f(i,j)-(u(i+1,j)+u(i-1,j)-4.*u(i,j)+u(i,j+1)+u(i,j-1))/dx/dx)**2
 	    end do
	end do   
    !max value of the residual
    err = dsqrt(err)/err0

    err_f = 0.d0
  
    !compute error
    do j=1,ny-1
	    do i=1,nx-1	   
        err_f = err_f + (u(i,j)-ue(i,j))**2
	    end do
	end do  

    err_f=dsqrt(err_f)/err_f0

 
	!write error
    write(19,*) nitr, err,err_f
    write(*,*)  nitr, err,err_f
         
end do

end 


!-----------------------------------------------------------------------------!
!Gauss-Seidel (GS)
!-----------------------------------------------------------------------------!
subroutine GS(nx,ny,dx,dy,f,ue,u,tol,nitr,x,y)   
implicit none
integer::nx,ny,i,j,nitr
real*8 ::dx,dy,tol,err,err0,a,err_f,err_f0,usa
real*8 ::f(0:nx,0:ny),ue(0:nx,0:ny),u(0:nx,0:ny),e(0:nx,0:ny)
real*8 :: x(0:nx)
real*8 :: y(0:ny)


a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)
usa = 1.d0/a

nitr = 0
err0=0.d0
err_f0=0.d0
 
!compute initial residual
do j=1,ny-1
    do i=1,nx-1
    err0 = err0+f(i,j)*f(i,j)
    end do
end do   
!max value of the residual
err0 = dsqrt(err0)

!compute initial error
 do j=1,ny-1
    do i=1,nx-1
    err_f0 = err_f0+ue(i,j)*ue(i,j)
    end do
end do   
!max value of the error
err_f0 = dsqrt(err_f0)

write(19,*) nitr, err0/err0, err_f0/err_f0
write(*,*)  nitr, err0/err0, err_f0/err_f0


err = err0   
err_f=err_f0
 
do while(err.ge.tol)

	nitr = nitr + 1
  
	do j=1,ny-1
	do i=1,nx-1
	u(i,j) = usa*(f(i,j) &
                   - (u(i+1,j)+u(i-1,j))/(dx*dx) &
                   - (u(i,j+1)+u(i,j-1))/(dy*dy) )
	end do
	end do
    
    err=0.d0
    !compute residual
 	do j=1,ny-1
	    do i=1,nx-1
	    err=err+(f(i,j)-(u(i+1,j)+u(i-1,j)-4.*u(i,j)+u(i,j+1)+u(i,j-1))/dx/dx)**2
	  	end do
	end do   
     err = dsqrt(err)/err0

    err_f = 0.d0 
  !compute error
    do j=1,ny-1
        do i=1,nx-1
        err_f = err_f + (ue(i,j)-u(i,j))**2
        end do
    end do   
    err_f =err_f/err_f0

   
	!write error
    write(19,*) nitr, err, err_f 
    write(*,*)  nitr, err, err_f
         
end do

end 

!-----------------------------------------------------------------------------!
!Successive over relaxation (SOR)
!-----------------------------------------------------------------------------!
subroutine SOR(nx,ny,dx,dy,f,ue,u,omega,tol,nitr,x,y) 
implicit none
integer::nx,ny,i,j,nitr
real*8 ::dx,dy,tol,err,err0,a,usa, err_f,err_f0,omega,up
real*8 ::f(0:nx,0:ny),ue(0:nx,0:ny),u(0:nx,0:ny),e(0:nx,0:ny),v(0:nx,0:ny)
real*8 :: x(0:nx)
real*8 :: y(0:ny)


a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)
usa = 1./a

    nitr = 0
    err0=0.d0
    err_f0=0.d0
  
    !compute initial residual
 	do j=1,ny-1
	    do i=1,nx-1
	    err0 = err0+f(i,j)*f(i,j)
	    end do
	end do   
    !max value of the error
    err0 = dsqrt(err0)

    !compute initial error
 	do j=1,ny-1
        do i=1,nx-1
        err_f0 = err_f0+ue(i,j)*ue(i,j)
        end do
    end do   
    !max value of the error
    err_f0 = dsqrt(err_f0)

    write(19,*) nitr, err0/err0, err_f0/err_f0
    write(*,*)  nitr, err0/err0, err_f0/err_f0
    
    v = 0.d0
    err=err0    
    err_f=err_f0
 
do while(err.ge.tol)

	nitr = nitr + 1
  

    v = u
    !update
	do j=1,ny-1
	do i=1,nx-1
	u(i,j) = (1.0d0/a)*(f(i,j) &
                   - (u(i+1,j)+u(i-1,j))/(dx*dx) &
                   - (u(i,j+1)+u(i,j-1))/(dy*dy) )
    
	u(i,j) = v(i,j) + omega*(u(i,j)-v(i,j))
	end do
	end do
  
    err=0.d0  
    !compute residual
 	do j=1,ny-1
	    do i=1,nx-1
	    err=err+ (f(i,j)-(u(i+1,j)+u(i-1,j)-4.*u(i,j)+u(i,j+1)+u(i,j-1))/dx/dx)**2
	    end do
	end do   
    err = dsqrt(err)/err0

    err_f=0.d0
    !compute the error
 	do j=1,ny-1
	    do i=1,nx-1
	    err_f=err_f+ (u(i,j)-ue(i,j))**2
	    end do
	end do   
    err_f = dsqrt(err_f)/err_f0

    
	!write error
    write(19,*) nitr, err, err_f
    write(*,*)  nitr, err, err_f
         
end do

end 






