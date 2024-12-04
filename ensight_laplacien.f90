SUBROUTINE write_result_ensight(xEnsight,yEnsight,u,ue,nx,ny,nz)


   implicit none

   real*8, dimension(nx,ny), intent(in) :: u,ue
   real*8, dimension(1,nx,ny,nz) :: scalar
   integer, parameter :: iskip = 1
   integer, parameter :: nbVar = 3
   integer :: nbImages,ntini
   character(len=80) :: caseName,sol_e_Name,sol_c_Name,diff_sol_Name
   character(len=80), dimension(nbVar) :: varName
   integer, dimension(nbVar) :: varType
   real*8,  dimension(nx) :: xEnsight
   real*8,  dimension(ny) :: yEnsight
   real*8,  dimension(1 ) :: zEnsight
   integer :: i,j,imin,imax,jmin,jmax,kmin,kmax,istep_rep
   integer :: nx,ny,nz,istep,isto,nstep
   
   write(caseName,'(A)') 'sortieEnsight'
   write(varName(1),'(A)') 'sol_calc'
   varType(1) = 1
   write(varName(2),'(A)') 'sol_exac'
   varType(2) = 1
   write(varName(3),'(A)') 'diff_sol'
   varType(3) = 1
    
   
   istep = 1
   
   zEnsight=0.
   
   kmin = 1
   kmax = 1 
   jmin = 1
   jmax = ny
   imin = 1
   imax = nx



   write(sol_c_Name,'(A)') 'sol_calc'
   write(sol_e_Name,'(A)') 'sol_exac'
   write(diff_sol_Name,'(A)') 'diff_sol'
   
   call EnsightCase(nbVar,varName,caseName,VarType,1,1,1)
   call WriteRectEnsightGeo(imin,imax,jmin,jmax,kmin,kmax,xEnsight,yEnsight,zEnsight,caseName)
  
   scalar = 0.d0
   do j=1,ny
      do i=1,nx   
      scalar(1,i,j,1)=u(i,j)       
      enddo
   enddo
   call WriteEnsightVar(1,nx,ny,nz,scalar,sol_c_Name ,imin,imax,jmin,jmax,jmin,kmax)

   scalar = 0.d0
   do j=1,ny
      do i=1,nx   
      scalar(1,i,j,1)=ue(i,j)       
      enddo
   enddo
    call WriteEnsightVar(1,nx,ny,nz,scalar,sol_e_Name,imin,imax,jmin,jmax,kmin,kmax)
    
    scalar = 0.d0
    do j=1,ny
      do i=1,nx   
      scalar(1,i,j,1)=ue(i,j)-u(i,j)       
      enddo
   enddo
    call WriteEnsightVar(1,nx,ny,nz,scalar,diff_sol_Name,imin,imax,jmin,jmax,kmin,kmax)




END SUBROUTINE write_result_ensight 

!
!    ******************************************************************************
!                      subrout  WriteRectEnsightGeo
!    ******************************************************************************
!
!    writes mesh data in Ensight's ascii format for rectilinear geometry
!
!    n1,n2,n3....: number of nodes in the x1,x2,x3 direction
!    x1,x2,x3....: coordinates 
!
      subroutine WriteRectEnsightGeo(imin,imax,jmin,jmax,kmin,kmax,x1,x2,x3,FileName) 
      implicit none

      INTEGER,INTENT(IN)::imin,imax,jmin,jmax,kmin,kmax
      REAL*8,DIMENSION(imax),INTENT(IN)::x1
      REAL*8,DIMENSION(jmax),INTENT(IN)::x2
      REAL*8,DIMENSION(1   ),INTENT(IN)::x3
      CHARACTER(LEN=80)     ,INTENT(IN)::FileName

	character(LEN=80)::binary_form
	character(LEN=80)::file_description1,file_description2
	character(LEN=80)::node_id,element_id
	character(LEN=80)::part,description_part,blocck

      integer::FileUnit,i,j,k,npart,isize,jsize,ksize
      integer::reclength

      FileUnit = 40

      binary_form      ='C Binary'
      file_description1='Ensight Model Geometry File Created by '
      file_description2='WriteRectEnsightGeo Routine'
      node_id          ='node id off'
      element_id       ='element id off'
      part             ='part'
      npart            =1
      description_part ='Sortie Ensight'
      blocck            ='block rectilinear'
      isize=(imax-imin+1)
      jsize=(jmax-jmin+1)
      ksize=(kmax-kmin+1)

      reclength=80*8+4*(4+isize+jsize+ksize)

        open (unit=FileUnit,file=trim(FileName)//'.geo',form='UNFORMATTED',access="direct",recl=reclength,status='replace')
        write(unit=FileUnit,rec=1) binary_form &
                                 ,file_description1 &
                                 ,file_description2 &
                                 ,node_id &
                                 ,element_id &
                                 ,part,npart &
                                 ,description_part &
                                 ,blocck &
                                 ,isize,jsize,ksize &
                                 ,(real(sngl(x1(i)),4),i=imin,imax) &
                                 ,(real(sngl(x2(j)),4),j=jmin,jmax) &
                                 ,(real(sngl(x3(k)),4),k=kmin,kmax) 
      
      close(FileUnit)

      end subroutine
!
!    ******************************************************************************
!
!    WriteEnsightSca writes result data in Ensight's format
!
!    m1,m2,m3....: size of the variable in the x1,x2,x3 direction
!    ndv.........: number of dimension of the variable (1=>scalar   3=>vector) 
!    var.........: data to be written
!    Varname.....: word used to build filenames
!    imin,imax...: range of writting data in the x1 direction
!    jmin,jmax...: range of writting data in the x2 direction
!    kmin,kmax...: range of writting data in the x3 direction
!
!    ******************************************************************************
      subroutine WriteEnsightVar(ndv,nx,ny,nz,var,VarName,imin,imax,jmin,jmax,kmin,kmax)

      implicit none

      INTEGER     ,INTENT(IN)::ndv,nx,ny,nz
      REAL*8,DIMENSION(ndv,nx,ny,nz),INTENT(IN)::var
      CHARACTER*80,INTENT(IN)::Varname

      character(len=80):: VarFileName
      character(len=80):: part,blocck
      integer::FileUnit,i,j,k,ii,npart,m,descripteur,tailleDoublePrecision,tailleCharacter,tailleInteger,tailleReel
      integer :: type_carre,type_cube,imin,imax,jmin,jmax,kmin,kmax,reclength
	 
 
      FileUnit = 41
      part ='part'
      npart=1
      blocck='block rectilinear'


      if (ndv.eq.1)VarFileName = trim(Varname)//'.scl'
      if (ndv.eq.3)VarFileName = trim(Varname)//'.vec'
       
      reclength=80*3+4*(1+(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)*ndv)

      open (unit=FileUnit,file=VarFileName,form='UNFORMATTED',access="direct",recl=reclength)



       write(unit=FileUnit,rec=1) VarFileName &
                                 ,part,npart,blocck &
                                 ,((((real(sngl(var(m,i,j,k)),4) &
                                 ,i=imin,imax) &
                                 ,j=jmin,jmax) &
                                 ,k=kmin,kmax) &
                                 ,m=1,ndv)
    
      close(FileUnit)

      end  subroutine
!
!    ******************************************************************************
!   EnsightCase helps to write a Ensight's case file
! 
!    VarName.....: Name of the variable
!    GeoName.....: Name of the geometrie
!    VarType.....: 1 => Scalar       3 => Vector
!    ntini.......: filename start number
!    nstop.......: filename end number
!    nprint......: filename increment
!
!    nfile.......: number of result files (time steps)
! 
      subroutine EnsightCase(nbVar,VarName,GeoName,VarType,istep,nstep,isto)
      implicit none

      INTEGER,INTENT(IN)::nbVar,nstep,isto,istep
      INTEGER, DIMENSION(nbVar),INTENT(IN)::VarType
      CHARACTER(LEN=80),DIMENSION(nbVar),INTENT(IN)::Varname
      CHARACTER(LEN=80),INTENT(IN)::GeoName
      integer::FileUnit,i,j,nfile

      nfile=1

      FileUnit = 40
      open(FileUnit,file=trim(GeoName)//'.case',status='replace')

      write(FileUnit,10) trim(GeoName)//'.geo'
   10 format('FORMAT'            ,/ ,'type: ensight gold',//,'GEOMETRY'          ,/ ,'model:	',A         ,//,'VARIABLE')

      do i=1,nbVar
        if(VarType(i).eq.1) write(FileUnit,15)trim(Varname(i)),trim(Varname(i))//'.scl'
        if(VarType(i).eq.2) write(FileUnit,25)trim(Varname(i)),trim(Varname(i))//'.vec'
      enddo
      
      write(FileUnit,45) nfile,isto,isto
      write(FileUnit,'(i10)') (j*isto,j=1,nfile)
      
      close(FileUnit)

   15 format('scalar per node: ',A,'   ', A)
   25 format('vector per node: ',A,'   ', A)

   45 format(/,'TIME            '      ,/,'time set: 1     '      ,/,'number of steps:'      ,i6 ,/, &
   'filename start number:',i10/,'filename increment:'   ,i10/,'time values: ')

      end subroutine
