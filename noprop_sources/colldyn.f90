module colldyn
!$ use OMP_LIB
use colllib
use splineinterp
use linearinterp
use matrices
implicit none

double precision, dimension(:,:), allocatable :: matintrp
double precision, dimension(:,:,:), allocatable :: mat2intrp
double complex, dimension(:), allocatable :: psi
double complex, dimension(:,:), allocatable ::  mat

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initdyn

double complex, dimension(:), allocatable ::  WORK
integer :: LWORK, INFO
integer, dimension(:), allocatable :: ihpsi

integer :: i

LWORK = 3*ntotsta
allocate(ihpsi(ntotsta),WORK(LWORK),mat(ntotsta,ntotsta))

! computes S-1
do i = 1, zgrid%na
  mat(:,:) = movl(i,:,:)
  CALL ZGETRF( ntotsta, ntotsta, mat, ntotsta, ihpsi, INFO )
  CALL ZGETRI( ntotsta, mat, ntotsta, ihpsi, WORK, LWORK, INFO )

!  write(*,*)i,matmul(mat,movl(i,:,:))

  if(INFO/=0) then
    write(*,*)"Error when computing S^-1"
    stop
  endif
  mcoup(i,:,:) = matmul(mat,mcoup(i,:,:))

enddo

  mat(:,:) = 0d0
deallocate(ihpsi,WORK,mat)

end subroutine initdyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hpsit(time,psiin,psiout,psidim,lhpsi,ihpsi,rhpsi,chpsi)
implicit none

logical    lhpsi(*)
integer    psidim, ihpsi(*)
real*8     rhpsi(*), time
complex*16 psiin(psidim), psiout(psidim), chpsi(*)

double precision :: rmat, cmat

integer :: i, j, ista, jsta

! --- INTERPOLATION OF MCOUP MATRIX

!$OMP PARALLEL DO COLLAPSE(2)  PRIVATE(i,j,rmat,cmat)
! SCHEDULE(DYNAMIC) 
do i = 1, ntsta
 do j = 1, ntsta
   call interp(tgrid%a,mcoup(:,j,i),tgrid%na,time,rmat,cmat)
   mat(j,i) = dcmplx(rmat,cmat)*exp(-dcmplx(0d0,1d0)*(testa(i)-testa(j))*time)
 enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)  PRIVATE(i,j,rmat,cmat,jsta,ista)
! SCHEDULE(DYNAMIC) 
do i = 1, npsta
 do j = 1, npsta
   ista = i +ntsta
   jsta = j +ntsta
   call interp(tgrid%a,mcoup(:,jsta,ista),tgrid%na,time,rmat,cmat)
   mat(jsta,ista) = dcmplx(rmat,cmat)*exp(-dcmplx(0d0,1d0)*(pesta(i)-pesta(j))*time)
 enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)  PRIVATE(i,j,rmat,cmat,jsta,ista)
! SCHEDULE(DYNAMIC) 
do i = 1, npsta
 do j = 1, ntsta
   ista = i +ntsta
   jsta = j 
   call interp(tgrid%a,mcoup(:,jsta,ista),tgrid%na,time,rmat,cmat)
   mat(jsta,ista) = dcmplx(rmat,cmat)*exp(-dcmplx(0d0,1d0)*(pesta(i)+0.5d0*vproj**2-testa(j))*time)
 enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)  PRIVATE(i,j,rmat,cmat,jsta,ista)
! SCHEDULE(DYNAMIC) 
do i = 1, ntsta
 do j = 1, npsta
   ista = i 
   jsta = j + ntsta
   call interp(tgrid%a,mcoup(:,jsta,ista),tgrid%na,time,rmat,cmat)
   mat(jsta,ista) = dcmplx(rmat,cmat)*exp(dcmplx(0d0,1d0)*(pesta(j)+0.5d0*vproj**2-testa(i))*time)
 enddo
enddo
!$OMP END PARALLEL DO

! --- COMPUTE ACTION OF HAMILTONIAN ON WAVEFUNCTION ---
 psiout = -dcmplx(0d0,1d0)*matmul(mat,psiin)

! call pmatvec(mat,psiin,psiout,psidim)
! write(60,'(60(f15.6,1X))')time,mat(1,1)*dcmplx(0d0,-1d0)/vproj
! write(156,'(60(f15.6,1X))')time,(cdabs(psiin(i))**2,i=1,ntotsta)


end subroutine hpsit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dyn

double complex, dimension(:), allocatable :: psit, dtpsit
double precision :: machprec
integer :: psidim

double precision :: time

double precision :: IntPeriod,AbsTime,LargeStep,TolError,ActualLargeStep,NextLargeStep, Step
integer :: IntOrder,SmallSteps,ErrorCode
external  AbsBSError,PolyExtrapol

double complex, dimension(:,:), allocatable :: AuxPsi !(psidim,intorder+2)
double precision, dimension(:), allocatable :: RData
double complex, dimension(:), allocatable :: CData
integer, dimension(:), allocatable :: IData
logical, dimension(:), allocatable :: LData

logical :: RestartABM
integer :: Steps, RepeatedSteps
double precision :: InitStep
External   AbsABMError

integer :: i

 psidim = ntotsta
 IntOrder = 6
 TolError = 1d-06

 allocate(Psit(ntotsta),dtPsit(ntotsta))
 allocate(RData(ntotsta),CData(ntotsta),IData(ntotsta),LData(ntotsta))
 allocate(AuxPsi(psidim,IntOrder+2))

! psi is declared in module general and must be given in main with the proper initial conditions
 Psit(:) = psi(:)

 RestartABM = .true.
 InitStep = abs(tgrid%a(2)-tgrid%a(1))*10
 Abstime = tgrid%a(2)
 IntPeriod = (tgrid%a(tgrid%na)-tgrid%a(5))

 call hpsit(Abstime,Psit,DtPsit,PsiDim,LData,IData,RData,CData)
 call  ABM (Psit,DtPsit,PsiDim,IntPeriod,AbsTime,IntOrder,              &
                           InitStep,TolError,RestartABM,Steps,          &
                           RepeatedSteps,ErrorCode,AuxPsi,              &
                           hpsit,AbsABMError,CData,RData,               &
                           IData,LData)

 Abstime = Abstime + IntPeriod

 if(ErrorCode /= 0) then
   write(*,*)'Error in ABM, I stop', ErrorCode
   stop
 endif

 psi(:) = Psit(:)

 deallocate(mat)
 deallocate(RData,CData,IData,LData)
 deallocate(AuxPsi)
 deallocate(Psit,dtPsit)

end subroutine dyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine pmatvec(mat,vecin,vecout,n)
!
! multiply n*n matrices
      integer i, j, k, n
      double complex  mat(n,n),vecin(n),vecout(n)

       vecout(:)=0d0
      do 30 i=1,n
        do 15 k=1,n
         vecout(i)=vecout(i)+mat(i,k)*vecin(k)
15      enddo
         vecout(i) = -dcmplx(0d0,1d0)*vecout(i)
30    enddo
      return
      end subroutine

end module colldyn
