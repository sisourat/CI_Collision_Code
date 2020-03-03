module phis
!!reads info from gamess outpt
!! only 1e and 2e matrix elements in the MO basis are needed
implicit none

double precision, dimension(:,:), allocatable :: h1eMO, kinMO, potMO
double complex, dimension(:,:), allocatable :: vp1eMO
double precision, dimension(:), allocatable :: int2e
integer :: nmos
integer*8 :: nint2e

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine freephis
if(allocated(h1eMO)) deallocate(h1eMO)
if(allocated(kinMO)) deallocate(kinMO)
if(allocated(potMO)) deallocate(potMO)
if(allocated(int2e)) deallocate(int2e)
end subroutine freephis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gusphis
implicit none

integer :: i, j

 open(unit=11,file='h1eMO.txt')
 read(11,*)nmos
 allocate(h1eMO(nmos,nmos),kinMO(nmos,nmos),potMO(nmos,nmos))
     do i=1,nmos
        do j=1,i
            read(11,*)h1eMO(i,j)
            h1eMO(j,i) = h1eMO(i,j)
        end do
     end do
 close(11)

 open(unit=11,file='kinMO.txt')
 read(11,*)nmos
     do i=1,nmos
        do j=1,i
            read(11,*)kinMO(i,j)
            kinMO(j,i) = kinMO(i,j)
        end do
     end do
 close(11)

 open(unit=11,file='potMO.txt')
 read(11,*)nmos
     do i=1,nmos
        do j=1,i
            read(11,*)potMO(i,j)
            potMO(j,i) = potMO(i,j)
        end do
     end do
 close(11)

write(*,*)"1e integrals read"

call read2e

write(*,*)"2e integrals read"

end subroutine gusphis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine read2e
implicit none

integer*8 :: intindex, ind
integer*8 ::  nt
double precision :: intt

integer(8) :: i,j,k,l,iloop

nt = nmos+1
nint2e = nt*(nt+1)*(nt*nt+nt+2)/8
 write(*,*) nint2e,"EEEEEE"

!write(*,*)nint2e,intindex(nmos,nmos,nmos,nmos)
!stop

allocate(int2e(nint2e))
!int2e=0.d0
iloop = 0
open(unit=10,file='moint.txt')
!open(unit=10,file='moint.dat',form='unformatted')

do while(iloop==0)
 read(10,*,end=100)ind,intt
             ! read(10,*,end=100)i,j,k,l,intt
             ! read(10,end=100)i,j,k,l,intt
             ! ind  = intindex(i,j,k,l) 
 int2e(ind) = intt
             ! write(*,*)i,j,k,l,ind, intt,nmos, nint2e
enddo

100 close(10)
end subroutine read2e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module phis
