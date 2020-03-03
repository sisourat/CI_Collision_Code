!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hpsi (psiin,psiout,psidim,lhpsi,ihpsi,rhpsi,chpsi)
      implicit none

      logical    lhpsi(*)
      integer    psidim, ihpsi(*)
      real*8     rhpsi(*), time
      complex*16 psiin(1:psidim), psiout(1:psidim), chpsi(*)

      integer :: i
! --- COMPUTE ACTION OF HAMILTONIAN ON WAVEFUNCTION ---

      double complex :: sclrprod
      external sclrprod

      do i = 1, psidim
         psiout(i) = 0d0
      enddo
!         psiout = matmul(mat,psiin)

      end subroutine hpsi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hampsi (psiin,psiout,psidim,lhpsi,ihpsi,rhpsi,chpsi)
      implicit none

      integer i
      complex*16 imag
      parameter (imag = (0.0d0,1.0d0))

      logical    lhpsi(*)
      integer    psidim, ihpsi(*)
      real*8     rhpsi(*), time
      complex*16 psiin(1:psidim), psiout(1:psidim), chpsi(*)

C --- COMPUTE ACTION OF HAMILTONIAN ON WAVEFUNCTION ---

C just a test
      double complex, dimension(256,256) :: mat
      common/matrix/mat

      double complex :: sclrprod
      external sclrprod

      do i = 1, psidim
         psiout(i) = 0d0
      enddo
         psiout = matmul(mat,psiin)

      end subroutine hampsi

C **********************************************************************
C *                                                                    *
C *                           LANCZOS (lanczos.f)                      *
C *                                                                    *
C * Library module containing a Lanczos algorithm for diagonalising    *
C * real symmetric or complex Hermitian matrices.                      *
C *                                                                    *
C * References:                                                        *
C *    [1] B. N. Parlett: "The Symmetric Eigenvalue Problem."          *
C *        Prentice-Hall, 1980.                                        *
C *    [2] J. K. Cullum and R. A. Willoughby: "Lanczos Algorithms for  *
C *        Large Symmetric Eigenvalue Computations. Vol. I: Theory."   *
C *        Birkhaeuser, 1985.                                          *
C *                                                                    *
C * Michael Beck, May 1999                                             *
C *                                                                    *
C **********************************************************************


C **********************************************************************
C *                                                                    *
C *                          SUBROUTINE LANCZOS                        *
C *                                                                    *
C * Purpose:                                                           *
C *    Master routine for diagonalising a real symmetric or complex    *
C *    Hermitian matrix.                                               *
C *                                                                    *
C * Input parameters:                                                  *
C *    order:       Number of Lanczos iterations made so far. Set      *
C *                 "order" to 0 for the initial call. (Allows one to  *
C *                 restart the Lanczos routine.)                      *
C *    psi:         Initial vector of length "psidim". If "order" > 0  *
C *                 this array must contain the first Lanczos vector.  *
C *    psidim:      Length of "psi".                                   *
C *    maxorder:    Number of Lanczos iterations to be made.           *
C *    machprec:    Machine precision (smallest x with 1.0+x > 1.0).   *
C *    calceigvec:  If true, the eigenvectors will be computed and     *
C *                 written to file by routine "eigvecout".            *
C *    iworkdim:    Length of "iwork" (must be >= maxorder).           *
C *    rworkdim:    Length of "rwork" (must be >= 6*maxorder-1 if      *
C *                 calceigvec = false, maxorder^2+4*maxorder-1 else). *
C *    cworkdim:    Length of "cwork" (must be >= 2*psidim).           *
C *    rwork:       For "order" > 0 the elements 1,...,order of this   *
C *                 real*8 array must contain the on-diagonal terms    *
C *                 "alpha(1:order)", and the elements order+1,...,    *
C *                 2*order the off-diagonal terms "beta(2:order+1)".  *
C *                 If "order" = 0 the contents of "rwork" is ignored. *
C *    cwork:       For "order" > 0 the first "psidim" elements of     *
C *                 this complex*16 array must contain the second      *
C *                 Lanczos vector. If "order" = 0 the contents of     *
C *                 "cwork" is ignored.                                *
C *                                                                    *
C * Output parameters:                                                 *
C *    psi:         This array is overwritten.                         *
C *    neigs:       Number of accurate eigenvalues found.              *
C *    eigval:      Vector containing the eigenvalues.                 *
C *    eigintens:   Vector containing the corresponding "intensities", *
C *                 i.e. overlaps between "psi" and the eigenvectors.  *
C *    eigerr:      Vector containing the error of the eigenvalues.    *
C *    errflag:     Error/warning flag having the following meaning:   *
C *                  0: everything o.k.                                *
C *                  1: Not enough integer work array                  *
C *                  2: Not enough real work array                     *
C *                  3: Not enough complex work array                  *
C *                  4: Illegal order specified                        *
C *                 11: Diagonalisation of tridiagonal matrix failed   *
C *                                                                    *
C * Other parameters:                                                  *
C *    lhpsi:       Logical data passed to "hpsi" and output routines. *
C *    ihpsi:       Integer data passed to "hpsi" and output routines. *
C *    rhpsi:       Real*8 data passed to "hpsi" and output routines.  *
C *    chpsi:       Complex*16 data passed to "hpsi" and output        *
C *                 routines.                                          *
C *    iwork:       Integer work array of length "iworkdim".           *
C *    rwork:       Real*8 work array of length "rworkdim".            *
C *    cwork:       Complex*16 work array of length "cworkdim".        *
C *                                                                    *
C * External routines:                                                 *
C *    sclrprod:    Complex function computing the scalar product of   *
C *                 two vectors. Called:                               *
C *                 z = sclrprod(bra-vector,ket-vector,vector_length)  *
C *    lanczout:    Routine that writes the current Lanczos vector to  *
C *                 file. Called:                                      *
C *                 call lanczout(current_order,lanczos_vector,        *
C *                               vector_length)                       *
C *    eigvecout:   Routine that writes the current eigenvector to     *
C *                 file. Called:                                      *
C *                 call eigvecout(number_of_eigenvector,eigenvector,  *
C *                                eigenvector_length,lhpsi,ihpsi,     *
C *                                rhpsi,chpsi)                        *
C *    restartout:  Routine that writes two Lanczos vectors and the    *
C *                 on- and off-diagonal elements alpha(1,...,order)   *
C *                 and beta(2,...,order+1) to file. (This data is     *
C *                 needed to restart the Lanczos routine.) Called:    *
C *                 call restartout(lanczos_vector1,lanczos_vector2,   *
C *                                 vector_length,alpha,beta,order,    *
C *                                 lhpsi,ihpsi,rhpsi,chpsi)           *
C *    hpsi:        Routine that multiplies the matrix to be           *
C *                 diagonalised with a vector. Called:                *
C *                 call hpsi(input_vector,resulting_vector,           *
C *                           vector_length,lhpsi,ihpsi,rhpsi,chpsi)   *
C *                                                                    *
C **********************************************************************

      subroutine lanczos (order, psi, psidim, maxorder, machprec,
     +                    calceigvec, iworkdim, rworkdim, cworkdim,
     +                    neigs, eigval, eigintens, eigerr, errflag,
     +                    sclrprod, lanczout, eigvecout, restartout,
     +                    hpsi, lhpsi, ihpsi, rhpsi, chpsi, iwork,
     +                    rwork, cwork)

      implicit none

      logical    calceigvec, lhpsi(*)
      integer    order, psidim, maxorder, iworkdim, rworkdim, cworkdim,
     +           neigs, errflag, iwork(1:iworkdim), ihpsi(*)
      real*8     machprec, eigval(1:maxorder),
     +           eigerr(1:maxorder), eigintens(1:maxorder),
     +           rwork(1:rworkdim), rhpsi(*)
      complex*16 psi(1:psidim), cwork(1:cworkdim), chpsi(*), sclrprod
      external   sclrprod, lanczout, eigvecout, restartout, hpsi

      integer neigvec, iwrk1, iwrk2, rwrk1, rwrk2, rwrk3, rwrk4, rwrk5,
     +        rwrk6, cwrk1, cwrk2, i

C --- CHECK THE ORDER ---

      if ((order .lt. 0) .or. (order .ge. maxorder)) then
         errflag = 4
         return
      endif

C --- COMPUTE POINTERS TO DIFFERENT PARTS OF WORK ARRAYS ---

      iwrk1       = 1
      iwrk2       = iwrk1+maxorder
      rwrk1       = 1
      rwrk2       = rwrk1+maxorder
      rwrk3       = rwrk2+maxorder
      rwrk4       = rwrk3+maxorder-1
      if (.not. calceigvec) then
         rwrk5 = rwrk4+2*maxorder
      else
         rwrk5 = rwrk4+maxorder**2
      endif
      rwrk6       = rwrk5+maxorder
      cwrk1       = 1
      cwrk2       = cwrk1+2*psidim

C --- CHECK THE LENGTH OF THE WORK SPACE ---

      if (iworkdim .lt. iwrk2-iwrk1) then
         errflag = 1
         return
      endif
      if (rworkdim .lt. rwrk6-rwrk1) then
         errflag = 2
         return
      endif
      if (cworkdim .lt. cwrk2-cwrk1) then
         errflag = 3
         return
      endif

C --- IN A RESTART RUN MOVE BETA TO CORRECT POSITION ---

      if (order .gt. 0) then
         do i = order,1,-1
            rwork(rwrk2-1+i) = rwork(rwrk1-1+order+i)
         enddo
      endif

C --- DEFINE HOW MANY COMPONENTS OF AN EIGENVECTOR ARE NEEDED ---

      if (calceigvec) then
         neigvec = maxorder
      else
         neigvec = 2
      endif

C --- CALL LANCZOS ROUTINE ---
      
      call lancz_slave(order,psi,psidim,maxorder,machprec,neigvec,neigs,
     +                 eigval,eigintens,eigerr,errflag,sclrprod,
     +                 lanczout,eigvecout,restartout,hpsi,lhpsi,ihpsi,
     +                 rhpsi,chpsi,iwork(iwrk1),rwork(rwrk1),
     +                 rwork(rwrk2),rwork(rwrk3),rwork(rwrk4),
     +                 rwork(rwrk5),cwork(cwrk1))

      return
      end

C **********************************************************************
C *                                                                    *
C *                       SUBROUTINE LANCZOS_SLAVE                     *
C *                                                                    *
C * Purpose:                                                           *
C *    Slave routine called by master Lanczos routine.                 *
C *                                                                    *
C * Parameters:                                                        *
C *    For most parameters see "lanczos". Others are given below.      *
C *                                                                    *
C * Input parameters:                                                  *
C *    neigvec: Number of components of eigenvectors to be computed.   *
C *                                                                    *
C * Output parameters:                                                 *
C *    eigstat: "Status" of an eigenvalue. See comment in code for a   *
C *             definition.                                            *
C *    alpha:   On-diagonal elements of tri-diagonal matrix.           *
C *    beta:    Off-diagonal elements of tri-diagonal matrix.          *
C *    eigval2: Eigenvalues of "reduced" tridiagonal matrix; needed    *
C *             for detecting spurious eigenvalues.                    *
C *    eigvec:  Matrix the columns of which contain the eigenvectors.  *
C *             Note that if no eigenvectors are desired, "eigvec"     *
C *             still contains the first and last component of each    *
C *             eigenvector.                                           *
C *    lancz:   Stores two Lanczos vectors (one is held in "psi").     *
C *                                                                    *
C * Other parameters:                                                  *
C *    work:    Real work array of length "maxorder".                  *
C *                                                                    *
C **********************************************************************

      subroutine lancz_slave (order, psi, psidim, maxorder, machprec,
     +                        neigvec, neigs, eigval, eigintens,
     +                        eigerr, errflag, sclrprod, lanczout,
     +                        eigvecout, restartout, hpsi, lhpsi, ihpsi,
     +                        rhpsi, chpsi, eigstat, alpha, beta,
     +                        eigval2, eigvec, work, lancz)

      implicit none

      logical    lhpsi(*)
      integer    order, psidim, maxorder, neigvec, eigstat(1:maxorder),
     +           neigs, errflag, ihpsi(*),wrorder
      real*8     machprec, alpha(1:maxorder), beta(2:maxorder+1),
     +           eigval(1:maxorder), eigval2(1:maxorder-1),
     +           eigerr(1:maxorder), eigintens(1:maxorder),
     +           eigvec(1:neigvec,1:maxorder), work(1:maxorder),
     +           rhpsi(*)
      complex*16 psi(1:psidim), lancz(1:psidim,1:2), sclrprod, chpsi(*)
      external   sclrprod, lanczout, eigvecout, restartout, hpsi

      integer p2, p3, swap, neigvec2, order2, i, j
      real*8  norm, invnorm, matsum, eps, eigdiff, error, dummy(1:1)

C --- INITIALISE VARIABLES ---

      errflag = 0
      p2      = 2
      p3      = 1
      if(maxorder .le. 1000) then 
         wrorder = max(maxorder/50,2)
      elseif(maxorder .le. 5000) then
         wrorder = maxorder/100
      else
         wrorder = maxorder/200
      endif
      if(wrorder .gt. 200) wrorder=100*((wrorder+33)/100) 
      if(wrorder .gt. 20)  wrorder=10*((wrorder+3)/10) 

C --- COMPUTE NORM OF INITIAL VECTOR ---

      if (order .eq. 0) then
         norm = sqrt(dble(sclrprod(psi,psi,psidim)))
      else
         norm = beta(order+1)
      endif

C --- LANCZOS ITERATION ---

C     This version of the Lanczos scheme is numerically the most stable,
C     see Chap. 2.3 of Ref. [2].

 100  continue
      order = order+1

C --- COMPUTE NORMALISED LANCZOS VECTOR ---

      invnorm = 1.0d0/norm
      do i = 1,psidim
         lancz(i,p2) = invnorm*psi(i)
      enddo

C --- WRITE LANCZOS VECTOR TO FILE IF DESIRED ---

      call lanczout(order,lancz(1,p2),psidim)

C --- EVALUATE ACTION OF HAMILTONIAN ON LANCZOS VECTOR ---

      call hpsi(lancz(1,p2),psi,psidim,lhpsi,ihpsi,rhpsi,chpsi)

C --- DETERMINE ON-DIAGONAL ELEMENT OF LANCZOS MATRIX ---

      alpha(order) = dble(sclrprod(lancz(1,p2),psi,psidim))

C --- GRAM-SCHMIDT-ORTHONORMALISE ON LAST LANCZOS VECTOR ---

      if (order .gt. 1) then
         do i = 1,psidim
            psi(i) = psi(i)-beta(order)*lancz(i,p3)
         enddo
      endif

C --- GRAM-SCHMIDT-ORTHONORMALISE ON LAST BUT ONE LANCZOS VECTOR ---

      do i = 1,psidim
         psi(i) = psi(i)-alpha(order)*lancz(i,p2)
      enddo

C --- DETERMINE OFF-DIAGONAL ELEMENT OF LANCZOS MATRIX ---

      norm = sqrt(dble(sclrprod(psi,psi,psidim)))
      beta(order+1) = norm

C --- SAVE CURRENT LANCZOS VECTOR ---

      swap = p2
      p2   = p3
      p3   = swap

C --- WRITE TO FILE DATA NEEDED FOR A RESTART ---

      call restartout(psi,lancz(1,p3),psidim,alpha,beta,order,lhpsi,
     +                ihpsi,rhpsi,chpsi)

C --- WRITE NUMBER OF ITERATION TO TO OUTPUT FILE

CNICO
C      if(mod(order,wrorder).eq.0. or. order.eq.maxorder)
C     +   call lancz_writer(order,wrorder)

C --- ITERATE UNTIL MAXIMUM ORDER REACHED ---

      if (order .lt. maxorder) goto 100

C --- CALCULATE EIGENVALUES AND -VECTORS OF TRIDIAGONAL MATRIX ---
       write(*,*)"NICO",order
      call eigtridiag(order,alpha,beta,machprec,neigvec,eigval,eigvec,
     +                errflag,work)
      if (errflag .ne. 0) then
         errflag = errflag+10
         return
      endif

C --- DETERMINE PRECISION PARAMETER ---

C     See Chap. 3.4 of Ref. [2].

      matsum = abs(alpha(1))
      do i = 2,order
         matsum = matsum+abs(alpha(i))+beta(i)
      enddo
      eps = dble(order)*matsum*machprec

C --- UNLESS STATED OTHERWISE ALL EIGENVALUES ARE GOOD ONES ---

C     The status "eigstat" of an eigenvalue is defined as:
C     eigstat =-1: inaccurate eigenvalue
C     eigstat = 0: single accurate eigenvalue
C     eigstat = i: multiple eigenvalue being identical to eigenvalue i
C     Note that multiple eigenvalues are always accurate (see Chap. 4.5
C     of Ref. [2]).

      do i = 1,order
         eigstat(i) = 0
      enddo

C --- DETECT AND LABEL MULTIPLE EIGENVALUES ---

      do i = 1,order
         if (eigstat(i) .eq. 0) then
            do j = i+1,order
               if (abs(eigval(i)-eigval(j)) .lt. eps) then
                  eigstat(i) = i
                  eigstat(j) = i
               endif
            enddo
         endif
      enddo

C --- CALCULATE EIGENVALUES OF REDUCED TRIDIAGONAL MATRIX ---

C     Eigenvalues of the full tridiagonal matrix which are also
C     eigenvalues of the same matrix, but with the first row and column
C     removed, are spurious (see Chap. 4.5 of Ref. [2]).

      neigvec2 = 0
      order2   = order-1
      call eigtridiag(order2,alpha(2),beta(3),machprec,neigvec2,eigval2,
     +                dummy,errflag,work)
      if (errflag .ne. 0) then
         errflag = errflag+10
         return
      endif

C --- DETECT AND LABEL SPURIOUS EIGENVALUES ---

      do i = 1,order
         if (eigstat(i) .eq. 0) then
            do j = 1,order-1
               if (abs(eigval(i)-eigval2(j)) .lt. eps) eigstat(i) = -1
            enddo
         endif
      enddo

C --- ESTIMATE ERROR OF EACH GOOD EIGENVALUE ---

      do i = 1,order
         if (eigstat(i) .ge. 0) then

C       --- DETERMINE INTENSITY ---

C           For multiple eigenvalues the intensities have to be added.

            if ((eigstat(i) .eq. 0) .or. (eigstat(i) .eq. i)) then
               eigintens(i) = eigvec(1,i)**2
            else
               eigintens(eigstat(i)) = eigintens(eigstat(i))
     +                                 +eigvec(1,i)**2
            endif

C       --- ESTIMATE ERROR ---

C           See Chap. 13.2 of Ref. [1] for the error formula.

            if (eigstat(i) .eq. 0) then
               eigdiff = 1.0d0/eps
               do j = 1,order
                  if (j .ne. i) then
                     eigdiff = min(eigdiff,abs(eigval(j)-eigval(i)))
                  endif
               enddo
               error     = beta(order+1)*abs(eigvec(neigvec,i))
               eigerr(i) = max(min(error,(error**2)/eigdiff),
     +                         eps/dble(order))
            else
               eigerr(i) = eps/dble(order)
            endif
         endif
      enddo

C --- REMOVE MULTIPLE AND SPURIOUS EIGENVALUES ---

      neigs = 0
      do i = 1,order
         if ((eigstat(i) .eq. 0) .or. (eigstat(i) .eq. i)) then
            neigs = neigs+1
            eigval(neigs)    = eigval(i)
            eigerr(neigs)    = eigerr(i)
            eigintens(neigs) = eigintens(i)
         endif
      enddo

C --- OUTPUT "GOOD" EIGENVECTORS ---

      if (neigvec .eq. maxorder) then
         neigs = 0
         do i = 1,order
            if ((eigstat(i) .eq. 0) .or. (eigstat(i) .eq. i)) then
               neigs = neigs+1
               call eigvecout(neigs,eigvec(1,i),order,lhpsi,ihpsi,rhpsi,
     +                        chpsi)
            endif
         enddo
      endif

      return
      end

C **********************************************************************
C *                                                                    *
C *                        SUBROUTINE LANCZOUT                         *
C *                                                                    *
C * Writes the current Lanczos vector to file. (In the current form    *
C * this is a dummy routine doing absolutely nothing; it is included   *
C * here for formal reasons but can (if desired) be extended such that *
C * it indeed outputs the Lanczos vector.                              *
C *                                                                    *
C * Input parameters:                                                  *
C *   order:    Current order of the Lanczos iteration.                *
C *   lanczvec: Current Lanczos vector.                                *
C *   dim:      Length of Lanczos vector.                              *
C *                                                                    *
C * Output parameters:                                                 *
C *   none                                                             *
C *                                                                    *
C * Other parameters:                                                  *
C *    lhpsi, ihpsi, rhpsi, chpsi: see above                           *
C *                                                                    *
C **********************************************************************

      subroutine lanczout (order, lanczvec, dim, lhpsi, ihpsi,
     +                         rhpsi, chpsi)

      implicit none

      logical    lhpsi(*)
      integer    order, dim, ihpsi(*)
      real*8     rhpsi(*)
      complex*16 lanczvec(1:dim), chpsi(*)

C --- WRITE LANCZOS VECTOR TO FILE ---

C Fill in your personal output commands here if desired.

      return
      end

! **********************************************************************
! *                                                                    *
! *                         FUNCTION SCLRPROD                    *
! *                                                                    *
! * Computes the Euklidian scalar product of two vectors.              *
! *                                                                    *
! * Input parameters:                                                  *
! *   psi1: "Bra" vector.                                              *
! *   psi2: "Ket" vector.                                              *
! *   dim:  Length of vectors.                                         *
! *                                                                    *
! * Output parameters:                                                 *
! *   none                                                             *
! *                                                                    *
! **********************************************************************

      complex*16 function sclrprod (psi1, psi2, dim)

      implicit none

      integer    dim
      complex*16 psi1(1:dim), psi2(1:dim)

      integer    i
      complex*16 s

! --- COMPUTE SCALAR PRODUCT ---

      s = dconjg(psi1(1))*psi2(1)
      do i = 2,dim
         s = s+dconjg(psi1(i))*psi2(i)
      enddo
      sclrprod = s

      return
      end function

C **********************************************************************
C *                                                                    *
C *                       SUBROUTINE EIGVECOUT                         *
C *                                                                    *
C * Writes the current eigenvector to file. (In the current form       *
C * this is a dummy routine doing absolutely nothing; it is included   *
C * here for formal reasons but can (if desired) be extended such that *
C * it indeed outputs the eigenvectors.                                *
C *                                                                    *
C * Input parameters:                                                  *
C *   number: Number of the current eigenvector.                       *
C *   eigvec: Current eigenvector.                                     *
C *   dim:    Length of eigenvector.                                   *
C *                                                                    *
C * Output parameters:                                                 *
C *   none                                                             *
C *                                                                    *
C * Other parameters:                                                  *
C *    lhpsi, ihpsi, rhpsi, chpsi: see above                           *
C *                                                                    *
C **********************************************************************

      subroutine eigvecout (number, eigvec, dim, lhpsi, ihpsi,
     +                            rhpsi, chpsi)

      implicit none

      logical    lhpsi(*)
      integer    number, dim, ihpsi(*)
      real*8     eigvec(1:dim), rhpsi(*)
      complex*16 chpsi(*)

      integer :: i, j

C --- WRITE EIGENVECTOR TO FILE ---

C Fill in your personal output commands here if desired.
      open(unit=10,file='lanceigvecout.txt',status='unknown',
     . position='append')
       write(10,*)number,(eigvec(j),j=1,dim)
      close(10)

      return
      end

C **********************************************************************
C *                                                                    *
C *                      SUBROUTINE RESTARTOUT                         *
C *                                                                    *
C * Writes the two current Lanczos vectors as well as the on- and      *
C * off-diagonal elements of the tridiagonal Lanczos matrix to file.   *
C * This data is required to restart the Lanczos routine. (In its      *
C * current form this is a dummy routine doing absolutely nothing; it  *
C * is included here for formal reasons but can (if desired) be        *
C * extended such that it indeed outputs the data.                     *
C *                                                                    *
C * Input parameters:                                                  *
C *   lanczvec1: first Lanczos vector.                                 *
C *   lanczvec2: second Lanczos vector.                                *
C *   lanczdim:  length of Lanczos vectors.                            *
C *   alpha:     on-diagonal elements of tridiagonal Lanczos matrix.   *
C *   beta:      off-diagonal elements of tridiagonal Lanczos matrix.  *
C *   order:     length of "alpha" and "beta" vectors.                 *
C *                                                                    *
C * Output parameters:                                                 *
C *   none                                                             *
C *                                                                    *
C * Other parameters:                                                  *
C *    lhpsi, ihpsi, rhpsi, chpsi: see above                           *
C *                                                                    *
C **********************************************************************

      subroutine restartout (lanczvec1, lanczvec2, lanczdim,
     +                             alpha, beta, order, lhpsi, ihpsi,
     +                             rhpsi, chpsi)

      implicit none

      logical    lhpsi(*)
      integer    lanczdim, order, ihpsi(*)
      real*8     alpha(1:order), beta(1:order), rhpsi(*)
      complex*16 lanczvec1(1:lanczdim), lanczvec2(1:lanczdim), chpsi(*)

C --- WRITE RESTART DATA TO FILE ---

C Fill in your personal output commands here if desired.

      return
      end

C **********************************************************************
C *                                                                    *
C *                        SUBROUTINE EIGTRIDIAG                       *
C *                                                                    *
C * Computes the eigenvalues and the first and last component of the   *
C * eigenvectors of a real symmetric tridiagonal matrix. This is a     *
C * modified version of the "tql2" routine.                            *
C *                                                                    *
C * Input parameters:                                                  *
C *   order:    Order (i.e. size) of the matrix.                       *
C *   ondiag:   Vector with the on-diagonal elements of the matrix.    *
C *   offdiag:  Vector with the off-diagonal elements of the matrix.   *
C *   machprec: Machine precision.                                     *
C *   neigvec:  Number of components of eigenvectors to be computed.   *
C *             Can be "0" (no eigenvectors), "2" (only first and last *
C *             component) or "order" (all components).                *
C *                                                                    *
C * Output parameters:                                                 *
C *   eigval:   Vector containing the eigenvalues of the matrix.       *
C *   eigvec:   Matrix with desired components of the eigenvectors.    *
C *   errflag:  Error flag having the meaning:                         *
C *             0: everything o.k.                                     *
C *             1: diagonalisation failed.                             *
C *             2: illegal value for "neigvec"                         *
C *                                                                    *
C * Other parameters:                                                  *
C *   work:     Real work array of length >= "order".                  *
C *                                                                    *
C **********************************************************************

      subroutine eigtridiag (order, ondiag, offdiag, machprec, neigvec,
     +                       eigval, eigvec, errflag, work)

      implicit none

      integer order, neigvec, errflag
      real*8  machprec, ondiag(1:order), offdiag(1:order-1),
     +        eigval(1:order), eigvec(1:neigvec,1:order),
     +        work(1:order)

      integer i, j, k, l, m, ii, l1, mml
      real*8  b, c, f, g, h, p, r, s

C --- INITIALISE VARIABLES ---

      errflag = 0
      do i = 1,order
         eigval(i) = ondiag(i)
      enddo
      do i = 1,order-1
         work(i) = offdiag(i)
      enddo
      work(order) = 0.0d0
      do i = 1,order
         do j = 1,neigvec
            eigvec(j,i) = 0.0d0
         enddo
      enddo
      if (neigvec .eq. 2) then
         eigvec(1,1)     = 1.0d0
         eigvec(2,order) = 1.0d0
      else
         do i = 1,neigvec
            eigvec(i,i) = 1.0d0
         enddo
      endif
      b = 0.0d0
      f = 0.0d0

C --- RETURN IMMEDIATELY IF MATRIX IS ONE-DIMENSIONAL ---

      if (order .eq. 1) return

C --- CHECK NUMBER OF EIGENVECTOR COMPONENTS ---

      if ((neigvec .ne. 0) .and. (neigvec .ne. 2) .and.
     +    (neigvec .ne. order)) then
         errflag = 2
         return
      endif

C --- COMPUTE EIGENVALUES ---

      do l = 1,order
         j = 0
         h = machprec*(abs(eigval(l))+abs(work(l)))
         if (b .lt. h) b = h

C    --- LOOK FOR SMALL SUB-DIAGONAL ELEMENT ---

         do m = l,order
            if (abs(work(m)) .le. b) goto 120
         enddo
 120     continue
         if (m .eq. l) goto 220
 130     continue
         if (j .eq. 30) then
            errflag = 1
            return
         endif
         j = j+1

C    --- FORM SHIFT ---

         l1 = l+1
         g = eigval(l)
         p = (eigval(l1)-g)/(2.0d0*work(l))
         r = sqrt(p*p+1.0d0)
         eigval(l) = work(l)/(p+sign(r,p))
         h = g-eigval(l)
         do i = l1,order
            eigval(i) = eigval(i)-h
         enddo
         f = f+h

C    --- QL TRANSFORMATION ---
         
         p = eigval(m)
         c = 1.0d0
         s = 0.0d0
         mml = m-l

C    --- FOR I = M-1 STEP -1 UNTIL L DO ---

         do ii = 1,mml
            i = m-ii
            g = c*work(i)
            h = c*p
            if (abs(p) .lt. abs(work(i))) goto 150
            c = work(i)/p
            r = sqrt(c*c+1.0d0)
            work(i+1) = s*p*r
            s = c/r
            c = 1.0d0/r
            goto 160
 150        continue
            c = p/work(i)
            r = sqrt(c*c+1.0d0)
            work(i+1) = s*work(i)*r
            s = 1.0d0/r
            c = c*s
 160        continue
            p = c*eigval(i)-s*g
            eigval(i+1) = h+s*(c*g+s*eigval(i))

C       --- FORM EIGENVECTOR ---

            do k = 1,neigvec
               h = eigvec(k,i+1)
               eigvec(k,i+1) = s*eigvec(k,i)+c*h
               eigvec(k,i) = c*eigvec(k,i)-s*h
            enddo
         enddo
         work(l)   = s*p
         eigval(l) = c*p
         if (abs(work(l)) .gt. b) goto 130
 220     continue
         eigval(l) = eigval(l)+f
      enddo

C --- ORDER EIGENVALUES AND EIGENVECTORS ---

      do ii = 2,order
         i = ii-1
         k = i
         p = eigval(i)
         do j = ii,order
            if (eigval(j) .lt. p) then
               k = j
               p = eigval(j)
            endif
         enddo
         if (k .ne. i) then
            eigval(k) = eigval(i)
            eigval(i) = p
            do j = 1,neigvec
               p = eigvec(j,i)
               eigvec(j,i) = eigvec(j,k)
               eigvec(j,k) = p
            enddo
         endif
      enddo

      return
      end
      
C **********************************************************************
C *                                                                    *
C *                       SUBROUTINE LANCZ_ERRMSG                      *
C *                                                                    *
C * Generates for a given error number returned by "lanczos" a corres- *
C * ponding error or warning message.                                  *
C *                                                                    *
C * Input parameters:                                                  *
C *   errflag: Error flag returned by "lanczos".                       *
C *                                                                    *
C * Output parameters:                                                 *
C *   errmsg: Error message.                                           *
C *                                                                    *
C **********************************************************************

      subroutine lancz_errmsg (errflag, errmsg)

      implicit none

      integer       errflag
      character*(*) errmsg

C --- GENERATE ERROR MESSAGE ---

      if (errflag .eq. 1) then
         errmsg = 'Not enough integer work array.'
      elseif (errflag .eq. 2) then
         errmsg = 'Not enough real work array.'
      elseif (errflag .eq. 3) then
         errmsg = 'Not enough complex work array.'
      elseif (errflag .eq. 4) then
         errmsg = 'Illegal order specified.'
      elseif (errflag .eq. 11) then
         errmsg = 'Diagonalisation of tridiagonal matrix failed.'
      elseif (errflag .eq. 12) then
         errmsg = 'Illegal number of eigenvector components.'
      else
         errmsg = 'Unknown error occurred.'
      endif

      return
      end
