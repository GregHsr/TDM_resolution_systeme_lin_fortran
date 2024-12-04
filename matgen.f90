	subroutine matgen(coef,jcoef,nnx,nny,mdim,dx,dy) 
	implicit none

	integer :: nnx,nny,i,j,k,mdim
	integer ,dimension(1:mdim) :: jcoef 


	real*8 ,dimension(1:nnx*nny,1:mdim)::coef 
	real*8 :: dx,dy,aa,bb

	k=1 
	do j=1,nny
	   do i=1,nnx
              coef(k,1) =2./dx/dx+2./dy/dy 
              coef(k,2) =-1./dx/dx 
              coef(k,3) =-1./dy/dy 
	      k=k+1
	   enddo
	enddo


	k=1
	do j=1,nny
	   do i=1,nnx
	      if(i.eq.nnx)then
	         coef(k,2)=0.0
	      endif
	      if(j.eq.nny)then
	         coef(k,3)=0.0
	      endif
	      k=k+1
	   enddo
	enddo

	jcoef(1)=0
	jcoef(2)=1
	jcoef(3)=nnx
	
	return
	end subroutine matgen
