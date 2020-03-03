program matColl
use phis
use SlaterDeterminant 
use omp_lib
implicit none

integer :: ndim
double precision, dimension(:,:), allocatable :: istavec, fstavec, mattmp, rmat, imat
double complex, dimension(:,:), allocatable :: matsta, ovlsta, matcsfs, ovlcsfs
double precision, dimension(:), allocatable :: istaeig, fstaeig

type csf
integer :: ndets
double precision, dimension(:), allocatable :: coeffs
type(Sdeterminant), dimension(:), allocatable :: dets 
end type csf

integer :: nicsfs, nialpha, nibeta
integer :: nfcsfs, nfalpha, nfbeta
type(csf), dimension(:), allocatable :: icsfs, fcsfs

character(2) :: ca, cb

character(50) :: fista, ffsta, filout, filename, filcoup, filovl
double precision :: prttol, mat, cikin, cipot, mat1e, ovl
double complex :: vp1e, c_ovl
integer :: nista, nfsta, nmainfsta

logical :: file_exists

integer :: nruns, irun, n1, n2, n3, iout

integer :: i, j, i1, j1, k

!for Collision matrix elements
double precision :: vproj, bproj, ta, tb
double complex, dimension(:,:,:), allocatable ::  mcoup, movl


integer :: ntime, ntotmo, itime
double precision, dimension(:), allocatable :: tgrid, esta

! for LAPACK S-1M
double complex, dimension(:), allocatable ::  WORK
integer :: LWORK, INFO
integer, dimension(:), allocatable :: ihpsi
double complex, dimension(:,:), allocatable ::  matdyn

  ta=OMP_GET_WTIME()
! read Collision matrix elements from previous Coll collision code
  read(*,*)filcoup,filovl
  open(unit=10,file=filcoup,form='unformatted')
    read(10)vproj,bproj
    write(*,*)vproj,bproj
    read(10)ntime,ntotmo
    write(*,*)ntime,ntotmo
    allocate(mcoup(ntime,ntotmo,ntotmo),tgrid(ntime),esta(ntotmo))
    allocate(h1eMO(ntotmo,ntotmo),vp1eMO(ntotmo,ntotmo))

    do i = 1, ntotmo
     read(10)esta(i)
    enddo

    do i = 1, ntime
      read(10)tgrid(i)
      do j = 1, ntotmo
        read(10)mcoup(i,j,:)
!!        write(2000+j,'(2(i4,1X),1000(f15.8,1X))')i,j,mcoup(i,j,:)
      enddo
    enddo
  close(10)

!  mcoup(:,1,5) = - mcoup(:,1,5)
!  mcoup(:,5,1) = - mcoup(:,5,1)

  open(unit=10,file=filovl,form='unformatted')
    read(10)vproj,bproj
    write(*,*)vproj,bproj 
    read(10)ntime,ntotmo
    write(*,*)ntime,ntotmo
    allocate(movl(ntime,ntotmo,ntotmo))

    do i = 1, ntotmo
     read(10)esta(i)
    enddo

    do i = 1, ntime
      read(10)tgrid(i)
      do j = 1, ntotmo
        read(10)movl(i,j,:)
      enddo
    enddo
  close(10)

! read CI stuffs
   read(*,*)nista, nfsta, fista, ffsta, filout
   if(nista .ne. nfsta) then
    write(*,*)'nista should be equal nfsta, I stop'
    stop
   endif
   iout=index(filout,'_ ')
   INQUIRE(FILE=fista, EXIST=file_exists)
   if(not(file_exists)) then
    write(*,*)fista,"does not exist, I stop"
    stop
   endif
   INQUIRE(FILE=ffsta, EXIST=file_exists)
   if(not(file_exists)) then
    write(*,*)ffsta,"does not exist, I stop"
    stop
   endif

  open(unit=22,file=fista)
! reads CSF
   read(22,*)nialpha,nibeta
   read(22,*)nicsfs
   allocate(icsfs(nicsfs),istavec(nicsfs,nista),istaeig(nista))

  do i = 1, nicsfs
    read(22,*)icsfs(i)%ndets
!    write(*,*)i,icsfs(i)%ndets
    allocate(icsfs(i)%dets(icsfs(i)%ndets),icsfs(i)%coeffs(icsfs(i)%ndets))
    do j = 1, icsfs(i)%ndets
      icsfs(i)%dets(j)%nalpha = nialpha
      icsfs(i)%dets(j)%nbeta = nibeta
      allocate(icsfs(i)%dets(j)%alpha(nialpha))
      allocate(icsfs(i)%dets(j)%beta(nibeta))
      read(22,*)icsfs(i)%coeffs(j)
      read(22,*)ca,(icsfs(i)%dets(j)%alpha(k),k=1,nialpha)
      read(22,*)cb,(icsfs(i)%dets(j)%beta(k),k=1,nibeta)
    enddo
  enddo

   do i = 1, nista
     read(22,*)ca,istaeig(i)
     do j = 1, nicsfs
       read(22,*)istavec(j,i)
!       write(*,*) istavec(j,i)

     enddo
!       write(*,*)sum(istavec(:,i)**2)
   enddo
  close(22)


  open(unit=23,file=ffsta)
! reads CSF
  read(23,*)nfalpha,nfbeta
  read(23,*)nfcsfs!, nmainfsta
  allocate(fcsfs(nfcsfs),fstavec(nfcsfs,nfsta),fstaeig(nfsta))
  do i = 1, nfcsfs
    read(23,*)fcsfs(i)%ndets
    allocate(fcsfs(i)%dets(fcsfs(i)%ndets),fcsfs(i)%coeffs(fcsfs(i)%ndets))
!    write(*,*)i,fcsfs(i)%ndets!fcsfs(i)%coeffs(j),ca,(fcsfs(i)%dets(j)%alpha(k),k=1,nfalpha),cb,(fcsfs(i)%dets(j)%beta(k),k=1,nfbeta)
    do j = 1, fcsfs(i)%ndets
      fcsfs(i)%dets(j)%nalpha = nfalpha
      fcsfs(i)%dets(j)%nbeta = nfbeta
      allocate(fcsfs(i)%dets(j)%alpha(nfalpha))
      allocate(fcsfs(i)%dets(j)%beta(nfbeta))
      read(23,*)fcsfs(i)%coeffs(j),ca,(fcsfs(i)%dets(j)%alpha(k),k=1,nfalpha),cb,(fcsfs(i)%dets(j)%beta(k),k=1,nfbeta)
!write(*,*)fcsfs(i)%coeffs(j),ca,(fcsfs(i)%dets(j)%alpha(k),k=1,nfalpha),cb,(fcsfs(i)%dets(j)%beta(k),k=1,nfbeta)
    enddo
  enddo

    fstavec(:,:) = 0d0
   do i = 1, nfsta
     read(23,*)ca,fstaeig(i)
     write(*,*)ca,fstaeig(i)

     do j = 1, nfcsfs
       read(23,*)fstavec(j,i)
     enddo
   enddo

  close(23)

! print the S-1M matrix coupling and energies for Prop
  open(unit=10,file="matCI_b",form='unformatted')
   write(10)vproj,bproj
   write(10)ntime,nista
   do i = 1, nista
     write(10)istaeig(i)
     write(*,*)istaeig(i)
   enddo
  allocate(matcsfs(nfcsfs,nicsfs),ovlcsfs(nfcsfs,nicsfs))
  LWORK = 3*nista
  allocate(ihpsi(nista),WORK(LWORK))
  allocate(matsta(nfsta,nista),ovlsta(nfsta,nista),matdyn(nfsta,nista))
  allocate(mattmp(nfcsfs,nista), rmat(nfsta,nista), imat(nfsta, nista))

!write(*,"(7f10.5)") mcoup(1,:,:)

 do itime = 1, ntime

  write(10)tgrid(itime)
!  write(*,*)
!  write(*,*)tgrid(itime)
  matcsfs(:,:) = 0d0
  ovlcsfs(:,:) = 0d0
  vp1eMO(:,:)=mcoup(itime,:,:)

   do i = 1, nicsfs
    do j = 1, nfcsfs
     do i1 = 1, icsfs(i)%ndets
       do j1 = 1, fcsfs(j)%ndets
          call collmat2dets(icsfs(i)%dets(i1),fcsfs(j)%dets(j1),vp1e,c_ovl)
          matcsfs(j,i) = matcsfs(j,i) + icsfs(i)%coeffs(i1)*fcsfs(j)%coeffs(j1)*vp1e
            ovlcsfs(j,i) = ovlcsfs(j,i) + icsfs(i)%coeffs(i1)*fcsfs(j)%coeffs(j1)*c_ovl
       enddo
     enddo
!!       if(itime==1) then
!!              write(1600,'(28f10.7)') abs( matcsfs(:,:))
!!       endif
    enddo
   enddo

!     if(itime==1) then
!        write(1000,'(28f12.7)') abs( matcsfs(:,:))
!     endif
 
  matsta(:,:) = 0d0
  ovlsta(:,:) = 0d0
  
!!!--------------------------------------------------------------------------
!!!--replace loops by matri-xmultiplication-----Junwen, nov.15,2019----------
!!!--------------------------------------------------------------------------

  rmat=0d0
  imat=0d0
  mattmp=0d0

!  matsta=matmul(transpose(fstavec),matmul(matcsfs,istavec))
  CALL DGEMM('N','N',nfcsfs,nista,nicsfs,1.d0,real(matcsfs),nfcsfs,istavec,nicsfs,0.d0,mattmp,nfcsfs)
  CALL DGEMM('N','N',nfsta,nista,nfcsfs,1.d0,transpose(fstavec),nfsta,mattmp,nfcsfs,0.d0,rmat,nfsta)
 

  mattmp=0d0
  CALL DGEMM('N','N',nfcsfs,nista,nicsfs,1.d0,aimag(matcsfs),nfcsfs,istavec,nicsfs,0.d0,mattmp,nfcsfs)
  CALL DGEMM('N','N',nfsta,nista,nfcsfs,1.d0,transpose(fstavec),nfsta,mattmp,nfcsfs,0.d0,imat,nfsta)

  matsta = dcmplx(rmat,imat)

! do i = 1, nista
!   do j = 1, nfsta
!    
!    do i1 = 1, nicsfs
!     do j1 = 1, nfcsfs

!         matsta(j,i) = matsta(j,i) + fstavec(j1,j)*istavec(i1,i)*matcsfs(j1,i1)
!           ovlsta(j,i) = ovlsta(j,i) + fstavec(j1,j)*istavec(i1,i)*ovlcsfs(j1,i1)

!     enddo
!   enddo

!  enddo
! enddo

!    if(itime==1) then
!       write(1500,'(28f10.7)') real( matsta(:,:))
!    endif


!!  do j = 1, nista
!!    write(200+j,'(2(i4,1X),10000(f15.8,1X))')i,j,matsta(j,:)
!!    write(300+j,'(1000(f15.8,1X))')tgrid(itime),ovlsta(j,:)
!!  enddo

!! S should actually be IDENTITY
! computes S-1
!  write(*,*)ovlsta(:,:)
!!  matdyn(:,:) = ovlsta(:,:) 
!!  CALL ZGETRF( nista, nista, matdyn, nista, ihpsi, INFO )
!!  CALL ZGETRI( nista, matdyn, nista, ihpsi, WORK, LWORK, INFO )

!!  if(INFO/=0) then
!!    write(*,*)"Error when computing S^-1"
!!    stop
!!  endif
!!  matsta(:,:) = matmul(matdyn,matsta(:,:))

  matdyn(:,:) = 0d0
  do j = 1, nista
!!    write(100+j,'(2(i4,1X),10000(f15.8,1X))')i,j,matsta(j,:)
    write(10)matsta(j,:)
  enddo
 enddo !itime

 close(10)

 deallocate(ihpsi,WORK)
 deallocate(rmat,imat,mattmp)
 deallocate(icsfs,istavec,istaeig)
 deallocate(fcsfs,fstavec,fstaeig)
 deallocate(matcsfs,ovlcsfs)
 deallocate(matsta,ovlsta,matdyn)

 tb=OMP_GET_WTIME()
 write(*,*)
 write(*,*)  "##############################"
 write(*,*)  "wall time for computing mats: ", tb-ta
 write(*,*)  "##############################"



call freephis

end program matColl
