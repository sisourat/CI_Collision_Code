module diagonalizer
use recipies
implicit none

 type diag
  character(4) :: typ
  integer :: maxit
 end type diag

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine eig(mydiag,ndim,mat,eigv)

 type(diag), intent(in) :: mydiag
 integer, intent(in) :: ndim
 double precision, dimension(ndim,ndim), intent(inout) :: mat
 double precision, dimension(ndim), intent(out) :: eigv

! for LAPACK diagonalizer
 double precision, dimension(:,:), allocatable :: U, t1, t2
 double precision, dimension(:), allocatable :: t3, d, e
 integer :: INFO, errflag

! for Lanczos
 integer :: order, maxorder
 integer :: psidim, iworkdim, cworkdim, rworkdim, neigs
 double precision, dimension(:), allocatable :: eigval
 double precision, dimension(:), allocatable :: eigintens, eigerr
 double precision, dimension(:), allocatable :: rwork
 double complex, dimension(:), allocatable :: cwork, psi
 integer, dimension(:), allocatable :: iwork

 logical :: calceigvec
 logical, dimension(:), allocatable :: lhpsi
 integer, dimension(:), allocatable :: ihpsi
 double precision, dimension(:), allocatable :: rhpsi
 double complex, dimension(:), allocatable :: chpsi
 double complex :: sclrprod
 external sclrprod, lanczout, eigvecout, restartout, hpsi

 double precision :: machprec
 double precision :: deriv, a, norme

 integer :: i, j

 if(mydiag%typ=='exac') then

  allocate(U(ndim,ndim))
  allocate(t1(ndim,ndim),t2(ndim,ndim),t3(4*ndim))
  allocate(lhpsi(ndim),ihpsi(ndim),rhpsi(ndim),chpsi(ndim))

  CALL DGEEV( 'N', 'V', ndim, mat, ndim, eigv, t1, t2, ndim, U, ndim, t3, 4*ndim, INFO )
  call eigsrt2(eigv,U,ndim,ndim)

  mat(:,:) = U(:,:)

  deallocate(U,t1,t2,t3)
  deallocate(lhpsi,rhpsi,ihpsi,chpsi)

 elseif(mydiag%typ=='lanc') then
   order = 0
   psidim = ndim
   maxorder = mydiag%maxit

   iworkdim = maxorder
   rworkdim = maxorder**2+4*maxorder-1
   cworkdim = 2*psidim

   allocate(psi(psidim),iwork(iworkdim),rwork(rworkdim),cwork(cworkdim))
   allocate(eigerr(maxorder),eigintens(maxorder),eigval(maxorder))
   allocate(lhpsi(ndim),ihpsi(ndim),rhpsi(ndim),chpsi(ndim))

     calceigvec=.false.
     call machin
     machprec = mprec

    do i = 1, ndim
      psi(i) = 1d0
    enddo
      psi(:)=psi(:)/sqrt(dble(sclrprod(psi,psi,psidim)))

   call lanczos (order, psi, psidim, maxorder, machprec,             &
                         calceigvec, iworkdim, rworkdim, cworkdim,   &
                         neigs, eigval, eigintens, eigerr, errflag,  &
                         sclrprod, lanczout, eigvecout, restartout,  &
                         hpsi, lhpsi, ihpsi, rhpsi, chpsi, iwork,    &
                         rwork, cwork)
   eigv(1:maxorder) = eigv(1:maxorder)

   deallocate(psi,iwork,rwork,cwork)
   deallocate(lhpsi,rhpsi,ihpsi,chpsi)

  endif

 end subroutine eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module diagonalizer
