program scaling
use phis
!use input
use SlaterDeterminant
use diagonalizer
implicit none

integer :: ndim
double precision, dimension(:,:), allocatable :: cimat, cimat1e, cimat2e, cikin, cipot, mattmp, matkin, matpot, matscal 
double precision, dimension(:), allocatable :: eigci, tempa

type csf
integer :: ndets
double precision, dimension(:), allocatable :: coeffs
type(Sdeterminant), dimension(:), allocatable :: dets 
end type csf

integer :: ncsfs, nalpha, nbeta
type(csf), dimension(:), allocatable :: csfs
double precision :: mat, mat1e, mat2e, kin, pot, ovl

type(diag) :: mydiag
character(8) :: cityp
character(1) :: ca, cb

double precision :: prttol
integer :: nsta, nmaincsfs

integer :: nruns, irun
character(50), dimension(:), allocatable :: fileigvec, filinput

integer :: nscale
double precision :: scalfact, scalmin, scalmax

double precision :: t1, t2, t3

integer :: i, j, i1, j1, k

mydiag%typ = 'exact'

 read(*,*)nruns
allocate(filinput(nruns),fileigvec(nruns))
do i = 1, nruns
 read(*,*)filinput(i)
enddo

  call cpu_time(t1)
  call gusphis
  call cpu_time(t2)
  write(*,*)"time read int",t2-t1

do irun = 1, nruns
 open(unit=21,file=filinput(irun))

 read(21,*)nsta, prttol
!  write(*,*)nsta, prttol
 read(21,*)fileigvec(irun)
!  write(*,*)fileigvec(irun)
 open(unit=22,file=fileigvec(irun))

! reads CSF
 read(21,*)nalpha,nbeta
 read(21,*)ncsfs !,nmaincsfs
 allocate(csfs(ncsfs))
 do i = 1, ncsfs
   read(21,*)csfs(i)%ndets
   allocate(csfs(i)%dets(csfs(i)%ndets),csfs(i)%coeffs(csfs(i)%ndets))
   do j = 1, csfs(i)%ndets
     csfs(i)%dets(j)%nalpha = nalpha
     csfs(i)%dets(j)%nbeta = nbeta
     allocate(csfs(i)%dets(j)%alpha(nalpha))
     allocate(csfs(i)%dets(j)%beta(nbeta))
     read(21,*)csfs(i)%coeffs(j),ca,(csfs(i)%dets(j)%alpha(k),k=1,nalpha),cb,(csfs(i)%dets(j)%beta(k),k=1,nbeta)
   enddo
 enddo
 read(21,*)nscale, scalmin, scalmax

!! creates CI matrix

 call cpu_time(t1)
 ndim = ncsfs
 allocate(cimat(ndim,ndim),cimat1e(ndim,ndim),cimat2e(ndim,ndim),eigci(ndim),tempa(ndim),cikin(ndim,ndim),cipot(ndim,ndim))
 allocate(mattmp(ndim,ndim),matkin(ndim,ndim),matpot(ndim,ndim),matscal(ndim,ndim))

 cimat(:,:) = 0d0
 cikin(:,:) = 0d0
 cipot(:,:) = 0d0

 cimat1e(:,:) = 0d0
 cimat2e(:,:) = 0d0

 do i = 1, ndim
   do j = 1, ndim 
     
    do i1 = 1, csfs(i)%ndets
      do j1 = 1, csfs(j)%ndets
         call mat2dets(csfs(i)%dets(i1),csfs(j)%dets(j1),mat1e,mat2e,kin,pot,ovl)

         cimat(j,i) = cimat(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*(mat1e+mat2e)
         cimat1e(j,i) = cimat1e(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*mat1e
         cimat2e(j,i) = cimat2e(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*mat2e

         cikin(j,i) = cikin(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*kin
         cipot(j,i) = cipot(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*pot

!         write(*,'(4(i3,1X),3(f20.15,1X))')j,i,j1,i1,csfs(i)%coeffs(i1),csfs(j)%coeffs(j1),mat
!         cimat(j,i) = i*j*1d0
      enddo
    enddo
!         write(*,'(2(i3,1X),2(f20.15,1X))')j,i,(cimat1e(j,i)+cimat2e(j,i))
 
   enddo
enddo
 call cpu_time(t2)
 write(*,*)"time build CI mat",t2-t1

!! scales and diagonalizes the CI matrix
  call cpu_time(t1)
  call eig(mydiag,ndim,cimat,eigci)
  call cpu_time(t2)
  write(*,*)"time for diag CI mat",t2-t1
  
!  mattmp = matmul(cimat1e+cimat2e,cimat)
!  matscal = matmul(transpose(cimat),mattmp)
!  do k = 1, ndim
!    write(*,*)eigci(i),matscal(i,i)
!  enddo
!  write(*,*)matmul(transpose(cimat),mattmp)


  call cpu_time(t1)
 do k = 1, nscale
  scalfact = scalmin+(k-1)*(scalmax-scalmin)/nscale
  matscal = scalfact**2*cikin(:,:) + scalfact*cipot(:,:)
  mattmp = matmul(matscal,cimat)
  matscal = matmul(transpose(cimat),mattmp)
!  cimat(:,:) = scalfact**2*cikin(:,:) + scalfact*cipot(:,:)
!  call eig(mydiag,ndim,cimat,eigci)
  write(22,'(10000(f15.8,1X))')scalfact,(matscal(j,j),j=1,ndim)
 enddo
  call cpu_time(t2)
  write(*,*)"time for scal",t2-t1

 close(21)
 close(22)

 deallocate(csfs,cimat,cimat1e,cimat2e,eigci,tempa,cikin,cipot,mattmp,matkin,matpot,matscal)
enddo

call freephis
deallocate(filinput,fileigvec)

end program scaling
