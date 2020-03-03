program myCI
use phis
!use input
use SlaterDeterminant
use diagonalizer
use diagolib
implicit none

integer :: ndim
double precision, dimension(:,:), allocatable :: cimat, cimat1e, cimat2e, ciovl
double precision, dimension(:), allocatable :: eigci, tempa

type csf
integer :: ndets
double precision, dimension(:), allocatable :: coeffs
type(Sdeterminant), dimension(:), allocatable :: dets 
end type csf

integer :: ncsfs, nalpha, nbeta
type(csf), dimension(:), allocatable :: csfs
double precision :: mat, mat1e, mat2e, cikin, cipot, ovl

type(diag) :: mydiag
character(8) :: cityp
character(1) :: ca, cb

double precision :: prttol
integer :: nsta, nmaincsfs

integer :: nruns, irun
character(50), dimension(:), allocatable :: fileigvec, filinput

logical :: file_exists

integer :: i, j, i1, j1, k

mydiag%typ = 'exact'

 read(*,*)nruns
allocate(filinput(nruns),fileigvec(nruns))
do i = 1, nruns
 read(*,*)filinput(i)
 INQUIRE(FILE=filinput(i), EXIST=file_exists)
 if(file_exists .eqv. 'F') then
   write(*,*)filinput(i),"does not exist, I stop"
   stop
 endif
enddo

call gusphis

do irun = 1, nruns

 open(unit=21,file=filinput(irun))

 read(21,*)nsta, prttol
 read(21,*)fileigvec(irun)
 write(*,*)fileigvec(irun)

 open(unit=22,file=fileigvec(irun))

! reads CSF
 read(21,*)nalpha,nbeta
 read(21,*)ncsfs !,nmaincsfs
 write(22,*)nalpha,nbeta
 write(22,*)ncsfs !,nmaincsfs
 allocate(csfs(ncsfs))
 do i = 1, ncsfs
   read(21,*)csfs(i)%ndets
   write(22,*)csfs(i)%ndets
!   write(*,*)i
   allocate(csfs(i)%dets(csfs(i)%ndets),csfs(i)%coeffs(csfs(i)%ndets))
   do j = 1, csfs(i)%ndets
     csfs(i)%dets(j)%nalpha = nalpha
     csfs(i)%dets(j)%nbeta = nbeta
     allocate(csfs(i)%dets(j)%alpha(nalpha))
     allocate(csfs(i)%dets(j)%beta(nbeta))
     read(21,*)csfs(i)%coeffs(j),ca,(csfs(i)%dets(j)%alpha(k),k=1,nalpha),cb,(csfs(i)%dets(j)%beta(k),k=1,nbeta)
     write(22,*)csfs(i)%coeffs(j)
     write(22,'(a,500(i4,1X))')ca,(csfs(i)%dets(j)%alpha(k),k=1,nalpha)
     write(22,'(a,500(i4,1X))')cb,(csfs(i)%dets(j)%beta(k),k=1,nbeta)
   enddo
 enddo

!! creates CI matrix

 ndim = ncsfs
 allocate(cimat(ndim,ndim),cimat1e(ndim,ndim),cimat2e(ndim,ndim),eigci(ndim),tempa(ndim),ciovl(ndim,ndim))

 cimat(:,:) = 0d0
 cimat1e(:,:) = 0d0
 cimat2e(:,:) = 0d0
 ciovl(:,:)= 0d0
 ovl=0d0
 do i = 1, ndim
   do j = 1, ndim 
     
    do i1 = 1, csfs(i)%ndets
      do j1 = 1, csfs(j)%ndets
         call mat2dets(csfs(i)%dets(i1),csfs(j)%dets(j1),mat1e,mat2e,cikin,cipot,ovl)
         cimat(j,i) = cimat(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*(mat1e+mat2e)
         cimat1e(j,i) = cimat1e(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*mat1e
         cimat2e(j,i) = cimat2e(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*mat2e
         ciovl(j,i) = ciovl(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*ovl


!         write(*,'(4(i3,1X),3(f20.15,1X))')j,i,j1,i1,csfs(i)%coeffs(i1),csfs(j)%coeffs(j1),mat
!         write(*,*)mat1e,mat2e,cikin+cipot
!         cimat(j,i) = i*j*1d0
      enddo
    enddo
!         write(*,'(2(i3,1X),2(f20.15,1X))')j,i,(cimat1e(j,i)+cimat2e(j,i))
 
   enddo
enddo

!do i=1,ndim
!write(100,'(10f12.7)') (ciovl(1:10,i))
!enddo
!write(*,*) ndim
!write(*,*) "*********************"
!do j=1,ndim
!write(*,'(10f12.7)') ( cimat(j,:))
!enddo

!write(*,*) "*********************"
!! diagonalizes the CI matrix
!ciovl=0d0
!do i=1,28
!  read(100,*) ciovl(i,:)
!enddo

 write(*,*)"DIAG"


!write(*,*) ciovl
!write(*,"(28(f10.5))") cimat

!!--The Diagonalizer used in GetSta-2e code------
    call eiglib(ndim, ciovl, cimat, eigci, cimat )
 


!!--The original Diagonalizer------
 !  call eig(mydiag,ndim,cimat,eigci)

    
 write(*,*)"NSTATES=",ndim

  do i = 1, nsta
    write(22,*)
    write(22,*)'E=',eigci(i),i
    write(*,*)'E=',eigci(i)
    do j = 1, ndim
!     if(abs(cimat(j,i))<prttol) cycle
      write(22,'(f20.16,1X,a,1X,i5,1X,a,500(i4))')cimat(j,i),"  ",j, "  ",(csfs(j)%dets(1)%alpha(k),k=1,nalpha),(csfs(j)%dets(1)%beta(k),k=1,nbeta)
    enddo
     write(22,*)
     tempa(:) = matmul(cimat1e(:,:),cimat(:,i))
     write(*,*)" ONE ELECTRON ENERGY = ", dot_product(cimat(:,i),tempa)
     tempa(:) = matmul(cimat2e(:,:),cimat(:,i))
     write(*,*)" TWO ELECTRON ENERGY = ", dot_product(cimat(:,i),tempa)
  enddo
 close(21)
 close(22)

 deallocate(csfs,cimat,cimat1e,cimat2e,ciovl,eigci,tempa)
enddo

call freephis
deallocate(filinput,fileigvec)

end program myCI
