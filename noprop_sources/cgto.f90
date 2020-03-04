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
  integer :: i, j, k, li, nli

  finput = ACenter%fbasis
  write(*,*)" # Read Basis set in", finput

  open(unit=iunit,file=finput)

   ACenter%nsizeBlocks=0
   do i = 1, nYlmMax
     call mysearch(YlmBlock(i),iunit,found)

     if(found==0) then
       read(iunit,*)h,ACenter%block(i)%ngto,ACenter%block(i)%ncgto
 
       li = i - 1
       nli = (li+1)*(li+2)/2
       ACenter%nsizeBlocks = ACenter%nsizeBlocks + nli*ACenter%block(i)%ncgto
       Acenter%nsizePerBlock(i) = nli*ACenter%block(i)%ncgto

       allocate(ACenter%block(i)%gto(ACenter%block(i)%ngto),ACenter%block(i)%cgto(ACenter%block(i)%ncgto,ACenter%block(i)%ngto))
       write(*,*)ACenter%block(i)%ngto,YlmBlock(i),'found'

       do j = 1, ACenter%block(i)%ngto
         read(iunit,*)ACenter%block(i)%gto(j),(ACenter%block(i)%cgto(k,j),k=1,ACenter%block(i)%ncgto)
       enddo 

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

!   k = 0
!   do i = li, 0, -1
!     do j = li - i, 0, -1
!       k = k + 1
!       ind(1,k) = i
!       ind(2,k) = j
!       ind(3,k) = li-i-j
!       write(*,*)k,ind(:,k)
!     enddo
!   enddo

! to follow GUS ordering
!      DATA BFNAM1/'  S ','  X ','  Y ','  Z ',
!     *            ' XX ',' YY ',' ZZ ',' XY ',' XZ ',' YZ ',
!     *            ' XXX',' YYY',' ZZZ',' XXY',' XXZ',
!     *            ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',
!     *            'XXXX','YYYY','ZZZZ','XXXY','XXXZ',
!     *            'YYYX','YYYZ','ZZZX','ZZZY','XXYY',
!     *            'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
!      DATA BFNAM2/' XXXXX',' YYYYY',' ZZZZZ',' XXXXY',' XXXXZ',
!     *            ' YYYYX',' YYYYZ',' ZZZZX',' ZZZZY',' XXXYY',
!     *            ' XXXZZ',' YYYXX',' YYYZZ',' ZZZXX',' ZZZYY',
!     *            ' XXXYZ',' YYYXZ',' ZZZXY',' XXYYZ',' XXZZY',
!     *            ' YYZZX',
!     *            '    X6','    Y6','    Z6','   X5Y','   X5Z',
!     *            '   Y5X','   Y5Z','   Z5X','   Z5Y','  X4Y2',
!     *            '  X4Z2','  Y4X2','  Y4Z2','  Z4X2','  Z4Y2',
!     *            '  X4YZ','  Y4XZ','  Z4XY','  X3Y3','  X3Z3',
!     *            '  Y3Z3',' X3Y2Z',' X3Z2Y',' Y3X2Z',' Y3Z2X',

  k = 0 
  if(li==0) then
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 0
   ind(3,k) = 0
  elseif(li==1) then
!px
   k = k + 1
   ind(1,k) = 1
   ind(2,k) = 0
   ind(3,k) = 0
!py
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 1
   ind(3,k) = 0
!pz
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 0
   ind(3,k) = 1
  elseif(li==2) then
!dxx
   k = k + 1
   ind(1,k) = 2
   ind(2,k) = 0
   ind(3,k) = 0
!dyy
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 2
   ind(3,k) = 0
!dzz
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 0
   ind(3,k) = 2
!dxy
   k = k + 1
   ind(1,k) = 1
   ind(2,k) = 1
   ind(3,k) = 0
!dxz
   k = k + 1
   ind(1,k) = 1
   ind(2,k) = 0
   ind(3,k) = 1
!dyz
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 1
   ind(3,k) = 1
  elseif(li==3) then
!fxxx
   k = k + 1
   ind(1,k) = 3
   ind(2,k) = 0
   ind(3,k) = 0
!fyyy
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 3
   ind(3,k) = 0
!fzzz
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 0
   ind(3,k) = 3
!fxxy
   k = k + 1
   ind(1,k) = 2
   ind(2,k) = 1
   ind(3,k) = 0
!fxxz
   k = k + 1
   ind(1,k) = 2
   ind(2,k) = 0
   ind(3,k) = 1
!fyyx
   k = k + 1
   ind(1,k) = 1
   ind(2,k) = 2
   ind(3,k) = 0
!fyyz
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 2
   ind(3,k) = 1
!fzzx
   k = k + 1
   ind(1,k) = 1
   ind(2,k) = 0
   ind(3,k) = 2
!fzzy
   k = k + 1
   ind(1,k) = 0
   ind(2,k) = 1
   ind(3,k) = 2
!fxyz
   k = k + 1
   ind(1,k) = 1
   ind(2,k) = 1
   ind(3,k) = 1
  else
   write(*,*)'Cartesian GTO with l=',li,'not implemented, I stop'
  endif

  end subroutine genPolgto

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cgto


