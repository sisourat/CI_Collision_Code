module propdyn
use general, only : imag
use splineinterp
use linearinterp
implicit none

double complex, dimension(:,:), allocatable :: matintrp
double precision, dimension(:,:,:), allocatable :: rmat2intrp, cmat2intrp
double complex, dimension(:), allocatable :: psi
double complex, dimension(:,:), allocatable ::  mat
double complex, dimension(:,:,:), allocatable ::  mcoup

integer :: ntime, ntotsta
double precision, dimension(:), allocatable :: tgrid
double precision, dimension(:), allocatable :: esta

contains

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
integer :: klo, khi, k

! --- INTERPOLATION OF MCOUP MATRIX
      if(tgrid(ksave).lt. time .and. tgrid(ksave+1).gt. time) then
       klo=ksave
       khi=ksave+1
      else
       klo=1
       khi=ntime
1      if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(tgrid(k).gt.time)then
           khi=k
         else
           klo=k
         endif
       goto 1
       endif
       ksave = klo
      endif

do i = 1, ntotsta
 do j = 1, ntotsta
!   call interp(klo,tgrid,mcoup(:,j,i),ntime,time,rmat,cmat)
   call splint(klo,khi,tgrid,real(mcoup(:,j,i)),aimag(mcoup(:,j,i)),rmat2intrp(:,j,i),cmat2intrp(:,j,i),ntime,time,rmat,cmat)
   matintrp(j,i) = dcmplx(rmat,cmat)*exp(-imag*(esta(i)-esta(j))*time)
 enddo
enddo

! --- COMPUTE ACTION OF HAMILTONIAN ON WAVEFUNCTION ---
 psiout = -imag*matmul(matintrp,psiin)

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

integer :: i, j

 psidim = ntotsta
 IntOrder = 6
 TolError = 1d-09

 allocate(Psit(ntotsta),dtPsit(ntotsta))
 allocate(RData(ntotsta),CData(ntotsta),IData(ntotsta),LData(ntotsta))
 allocate(AuxPsi(psidim,IntOrder+2))

! for spline interpolation
ksave = 1
do i = 1, ntotsta
 do j = 1, ntotsta
   call spline(tgrid,real(mcoup(:,j,i)),ntime,0d0,0d0,rmat2intrp(:,j,i))
   call spline(tgrid,aimag(mcoup(:,j,i)),ntime,0d0,0d0,cmat2intrp(:,j,i))
 enddo
enddo

! psi is declared in module general and must be given in main with the proper initial conditions
 Psit(:) = psi(:)

 RestartABM = .true.
 InitStep = abs(tgrid(2)-tgrid(1))
 Abstime = tgrid(2)
 IntPeriod = (tgrid(ntime)-tgrid(5))

 call hpsit(Abstime,Psit,DtPsit,PsiDim,LData,IData,RData,CData)
 call  ABM (Psit,DtPsit,PsiDim,IntPeriod,AbsTime,IntOrder,              &
                           InitStep,TolError,RestartABM,Steps,          &
                           RepeatedSteps,ErrorCode,AuxPsi,              &
                           hpsit,AbsABMError,CData,RData,               &
                           IData,LData)

 !call rk4(hpsit, Abstime, Psit, 0.005d0, PsiDim, tgrid(ntime-3), lData,iData,rData,cData) 

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

subroutine rk4 (f,t,x,h,n,tmax,lhpsi,ihpsi,rhpsi,chpsi)   
        integer, intent (in) :: n
        double complex, dimension(n) :: f1, f2, f3, f4, ta
        double precision :: t, tmax, h
        double complex, dimension(n) , intent(inout) :: x

        logical    lhpsi(*)
        integer    psidim, ihpsi(*)
        real*8     rhpsi(*)
        complex*16 chpsi(*)
        
        external f

        do while(t<tmax)
           call f(t, x, f1, n, lhpsi,ihpsi,rhpsi,chpsi)       
           call f(t + h/2.0, x + 0.5*h*f1, f2, n, lhpsi,ihpsi,rhpsi,chpsi)   
           call f(t + h/2.0, x + 0.5*h*f2, f3, n, lhpsi,ihpsi,rhpsi,chpsi)   
           call f(t + h, x + h*f3, f4, n, lhpsi,ihpsi,rhpsi,chpsi)
           x = x + h*(f1 + 2*f2  +2*f3 + f4)/6.0
           t = t + h
        end do

end subroutine rk4

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
         vecout(i) = -imag*vecout(i)
30    enddo
      return
      end subroutine

end module propdyn
