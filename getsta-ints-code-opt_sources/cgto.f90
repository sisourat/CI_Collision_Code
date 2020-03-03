module cgto 
use tools
use newtypes
use general
implicit none

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_cgto(ACenter)

  type(center), intent(inout) :: ACenter

  character(lenmax) :: h, finput
  integer, parameter :: iunit=10

  integer :: found
  integer :: i, j, k, li, nli, jj, kk

  finput = ACenter%fbasis
  write(*,*)" # Read Basis set in", finput

  open(unit=iunit,file=finput)

   do i = 1, nYlmMax
     call mysearch(YlmBlock(i),iunit,found)

     if(found==0) then
       read(iunit,*)h,ACenter%block(i)%ngto,ACenter%block(i)%ncgto
 
       li = i - 1
       nli = (li+1)*(li+2)/2     !---the number of GTO type in Cartesian---
       ACenter%nsizeBlocks = ACenter%nsizeBlocks + nli*ACenter%block(i)%ncgto
       Acenter%nsizePerBlock(i) = nli*ACenter%block(i)%ncgto

       ACenter%nsizegtoBlocks =ACenter%nsizegtoBlocks+ nli*ACenter%block(i)%ngto
       Acenter%nsizegtoPerBlock(i) = nli*ACenter%block(i)%ngto


       allocate(ACenter%block(i)%gto(ACenter%block(i)%ngto),ACenter%block(i)%cgto(ACenter%block(i)%ncgto,ACenter%block(i)%ngto))
       write(*,*)ACenter%block(i)%ngto,YlmBlock(i),'found'

       do j = 1, ACenter%block(i)%ngto
         read(iunit,*)ACenter%block(i)%gto(j),(ACenter%block(i)%cgto(k,j),k=1,ACenter%block(i)%ncgto)
       enddo

!--------------------allocate it-------------------

       allocate( ACenter%block(i)%ncgtoreal(ACenter%block(i)%ngto) )

       allocate( ACenter%block(i)%nincgto(ACenter%block(i)%ncgto, ACenter%block(i)%ngto) )


       

    

       ACenter%block(i)%ncgtoreal = 0

       do jj = 1, ACenter%block(i)%ngto
         do kk = 1, ACenter%block(i)%ncgto

             if( ACenter%block(i)%cgto(kk,jj) .ne. 0.0) then


                 ACenter%block(i)%ncgtoreal(jj) =  ACenter%block(i)%ncgtoreal(jj)   +  1  

             
                 ACenter%block(i)%nincgto( ACenter%block(i)%ncgtoreal(jj),jj) = kk 


 !                 write(*,*)  ACenter%block(i)%ncgtoreal(jj),jj 
!
           !       write(*,*)  ACenter%block(i)%nincgto(ACenter%block(i)%ncgtoreal(jj) ,jj) 

             endif

         enddo
       enddo

!-------------------------------------------------



 

     else

       write(*,*)YlmBlock(i),' not found'
       ACenter%block(i)%ngto=0
       ACenter%block(i)%ncgto=0

     endif

   enddo

  close(iunit)

  end subroutine read_cgto

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine normgto(lmax,alp,ngto)

  integer, intent(in) :: lmax
  double precision, intent(in) :: alp
  double precision, dimension(0:lmax,0:lmax,0:lmax), intent(inout) :: ngto

  integer :: l, m, n

!!$OMP THREADPRIVATE(l, m, n)

!N = sqrt( alpi**(li+0.5)*2**(2*li+0.5)/(pi**0.5*fact2(2*li-1))  )
  do l = 0, lmax
    do m = 0, lmax
      do n = 0, lmax
        ngto(n,m,l) = dsqrt( alp**(l+m+n+1.5d0)*2d0**(2*(l+m+n)+1.5)/(pi**1.5*fact2(2*l-1)*fact2(2*m-1)*fact2(2*n-1))  )
      enddo
    enddo
  enddo

  end subroutine normgto

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine genPolgto(li,nli,ind)

  integer, intent(in) :: li, nli
  integer, dimension(3,nli), intent(inout) :: ind

  integer :: i, j, k 
!!$OMP THREADPRIVATE(i, j, k )


!------The ordering of GTOs now is consistent with GAMESS-US, Aug 7, 2019, by Junwen----- 

  if((li .eq. 0) .or.(li .eq. 1) )then
    ! S and P orbital

   k = 0
   do i = li, 0, -1
     do j = li - i, 0, -1
       k = k + 1
       ind(1,k) = i
       ind(2,k) = j
       ind(3,k) = li-i-j
!       write(*,*)k,ind(:,k),"DDDD"
     enddo
   enddo


  else if(li .eq. 2 )then
    ! D orbital

    ind=reshape( (/2,0,0,&               !XX
                   0,2,0,&               !YY
                   0,0,2,&               !ZZ
                   1,1,0,&               !XY
                   1,0,1,&               !XZ
                   0,1,1/),(/3,nli/))    !YZ

  else if(li .eq. 3)then 
    ! F orbital

    ind=reshape( (/3,0,0,&               !XXX
                   0,3,0,&               !YYY
                   0,0,3,&               !ZZZ
                   2,1,0,&               !XXY
                   2,0,1,&               !XXZ
                   1,2,0,&               !XYY
                   0,2,1,&               !YYZ
                   1,0,2,&               !XZZ
                   0,1,2,&               !YZZ
                   1,1,1/),(/3,nli/))    !XYZ

  else if(li .eq. 4)then 

    ! G orbital

    ind=reshape( (/4,0,0,&              !XXXX              
                   0,4,0,&              !YYYY
                   0,0,4,&              !ZZZZ
                   3,1,0,&              !XXXY
                   3,0,1,&              !XXXZ
                   1,3,0,&              !XYYY
                   0,3,1,&              !YYYZ
                   1,0,3,&              !XZZZ
                   0,1,3,&              !YZZZ
                   2,2,0,&              !XXYY
                   2,0,2,&              !XXZZ
                   0,2,2,&              !YYZZ
                   2,1,1,&              !XXYZ
                   1,2,1,&              !XYYZ
                   1,1,2/),(/3,nli/))   !XYZZ
  else

  write(*,*) "STOP: The ordering of H orbital is not implenmeted, check genPolgto in cgto.f90"
  stop 
  endif


 


  end subroutine genPolgto

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cgto


