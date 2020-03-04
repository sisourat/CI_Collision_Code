module ComputeMat
use newtypes
use overlap_mod
use nuclear_mod
use fdn_mod
use cgto, only : normgto, genPolgto
use centerlib
use matrices
implicit none
 

  contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ComputeStaMat(geom)

  integer, dimension(*) :: geom

  integer :: i, j, iy, jy, ist, jst, ifin, jfin, li, lj, nli, nlj
  integer :: igto, jgto, icgto, jcgto, im, jm, kicgto, kjcgto
  integer :: lmax, tmax, ipot, jpot, kpot

  double precision :: xi, yi, zi, xj, yj, zj, xp, yp, zp, alpi, alpj, Gam, rab2, kfactor
  double precision :: vx, vy, vz, cx, cy, cz, bet, alp

  double complex, dimension(:,:,:), allocatable :: ovmat, potmat, pottmp
  double complex, dimension(:,:), allocatable :: ovGto, potGto, kinGto
  double precision, dimension(:,:,:), allocatable ::  fkx, fky, fkz
  double precision, dimension(:,:,:), allocatable :: normGtoI, normGtoJ

  do i = 1, nCenters
   do iy = 1, nYlmMax

    if(Centers(i)%block(iy)%ncgto==0) cycle

    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, nCenters
      do jy = 1, nYlmMax

        if(Centers(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

        do igto = 1, Centers(i)%block(iy)%ngto
          do jgto = 1, Centers(j)%block(jy)%ngto

             xi = Centers(i)%x(geom(i)); yi = Centers(i)%y(geom(i)); zi = Centers(i)%z(geom(i))
             xj = Centers(j)%x(geom(j)); yj = Centers(j)%y(geom(j)); zj = Centers(j)%z(geom(j))

             alpi = Centers(i)%block(iy)%gto(igto); alpj = Centers(j)%block(jy)%gto(jgto)
             Gam =  alpi + alpj

             rab2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
             kfactor = exp(-alpi*alpj*rab2/Gam)
             xp  =  ( alpi*xi + alpj*xj ) / Gam ; yp  =  ( alpi*yi + alpj*yj ) / Gam ; zp  =  ( alpi*zi + alpj*zj ) / Gam
             lmax = li+lj+2 ; tmax = 3*lmax ! +2 is for kinetic energy

! ovGto contains overlap matrix elements needed for all possible x**l*y**m*z**n*exp(-Gam*r**2)  
! potGto contains nuclear attraction matrix elements 

             allocate(ovmat(0:lmax,0:lmax,0:lmax))
             allocate(potmat(0:lmax,0:lmax,0:lmax),pottmp(0:lmax,0:lmax,0:lmax))

             allocate(normGtoI(0:lmax,0:lmax,0:lmax),normGtoJ(0:lmax,0:lmax,0:lmax))
             call normgto(lmax,alpi,normGtoI)
             call normgto(lmax,alpj,normGtoJ)

             call overlap_driver(lmax,Gam,0d0,0d0,0d0,ovmat) 
             ovmat(:,:,:) = kfactor*ovmat(:,:,:)

! computes nuclear integrals over all centers
             potmat(:,:,:) = 0d0

             do ipot = 1, nCenters
! cx, cy, cz are shifted by -xp, -yp, -zp respectively
! rp is considered as the center of referential to compute the integrals
               cx = -Centers(ipot)%x(geom(ipot)) + xp ; cy = -Centers(ipot)%y(geom(ipot)) + yp ; cz = -Centers(ipot)%z(geom(ipot)) + zp

               do jpot = 1, Centers(ipot)%pot%n
                 bet = Centers(ipot)%pot%alp(jpot)

! no need to compute nuclear up to li+lj+2, li+lj is enough; save some cpu time
                 call nuclear_driver(lmax-2,0,Gam,0d0,0d0,0d0,bet,cx,cy,cz,tmax,pottmp(0:lmax-2,0:lmax-2,0:lmax-2))
                 potmat(:,:,:) = potmat(:,:,:) - Centers(ipot)%pot%c(jpot)*pottmp(:,:,:)

!                WRITE(*,'(i3,100(f20.12,1X))')ipot,cx,cy,cz,Gam
!                 call printmat(lmax,lmax,lmax,pottmp)

               enddo ! jpot
             enddo ! ipot

             potmat(:,:,:) = kfactor*potmat(:,:,:)

             allocate(fkx(0:lmax,0:lmax,0:lmax),fky(0:lmax,0:lmax,0:lmax),fkz(0:lmax,0:lmax,0:lmax))

! fkz, fky and fkz contain the coeff. from Gaussian product rule 
             call fdn(xi-xp,xj-xp,lmax,fkx) 
             call fdn(yi-yp,yj-yp,lmax,fky) 
             call fdn(zi-zp,zj-zp,lmax,fkz) 

! now the matrices are computed according to the Gaussian product rule for all
! (nli+1)*(nli+2)/2 and (nlj+1)*(nlj+2)/2 cartesian gaussian functions

             allocate(ovGto(nlj,nli),potGto(nlj,nli),kinGto(nlj,nli)) 

             call sumfk(lmax,li,lj,nli,nlj,fkx,fky,fkz,normGtoI,normGtoJ,ovmat,ovGto)
             call sumfk(lmax,li,lj,nli,nlj,fkx,fky,fkz,normGtoI,normGtoJ,potmat,potGto)

! kinetic matrix is obtained from the overlap (Delta_r = Delta_rj)
             call kinetic(alpj,lmax,li,lj,nli,nlj,fkx,fky,fkz,normGtoI,normGtoJ,ovmat,kinGto)

! finally the matrices are computed in the CGto basis
             kicgto = 0; kjcgto = 0
             do icgto = 1, Centers(i)%block(iy)%ncgto
               do jcgto = 1, Centers(j)%block(jy)%ncgto

                do im = 1, nli
                 do jm = 1, nlj

                  kicgto = sum(Centers(1:i-1)%nsizeBlocks) + sum(Centers(i)%nsizePerBlock(1:iy-1)) + (icgto-1)*(nli) + im
                  kjcgto = sum(Centers(1:j-1)%nsizeBlocks) + sum(Centers(j)%nsizePerBlock(1:jy-1)) + (jcgto-1)*(nlj) + jm

                  tCGTO(kjcgto,kicgto) = tCGTO(kjcgto,kicgto) + Centers(i)%block(iy)%cgto(icgto,igto)*Centers(j)%block(jy)%cgto(jcgto,jgto)*ovGto(jm,im)
                  PotCGTO(kjcgto,kicgto) = PotCGTO(kjcgto,kicgto) + Centers(i)%block(iy)%cgto(icgto,igto)*Centers(j)%block(jy)%cgto(jcgto,jgto)*potGto(jm,im)
                  KinCGTO(kjcgto,kicgto) = KinCGTO(kjcgto,kicgto) + Centers(i)%block(iy)%cgto(icgto,igto)*Centers(j)%block(jy)%cgto(jcgto,jgto)*kinGto(jm,im)

                 enddo
                enddo

               enddo
             enddo


             deallocate(ovmat,ovGto,fkx,fky,fkz,normGtoI,normGtoJ,potmat,pottmp,potGto,kinGto)

          enddo
        enddo

      enddo
     enddo

    enddo
  enddo

! renormalize the CGTO
       do kicgto = 1, nsizeCGto  
         do kjcgto = 1, nsizeCGto  
             OvCGTO(kjcgto,kicgto) = tCGTO(kjcgto,kicgto)/sqrt(tCGTO(kicgto,kicgto)*tCGTO(kjcgto,kjcgto))
             PotCGTO(kjcgto,kicgto) = PotCGTO(kjcgto,kicgto)/sqrt(tCGTO(kicgto,kicgto)*tCGTO(kjcgto,kjcgto))
             KinCGTO(kjcgto,kicgto) = KinCGTO(kjcgto,kicgto)/sqrt(tCGTO(kicgto,kicgto)*tCGTO(kjcgto,kjcgto))
         enddo
       enddo

  end subroutine ComputeStaMat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sumfk(lmax, li, lj, nli, nlj, fkx, fky, fkz, normgtoi, normgtoj, mat, mgto)
  integer, intent(in) :: lmax, li, lj, nli, nlj
  double complex, dimension(0:lmax,0:lmax,0:lmax) :: mat
  double precision, dimension(0:lmax,0:lmax,0:lmax) :: fkx, fky, fkz, normgtoi, normgtoj
  double complex, dimension(nlj, nli) :: mgto
  
  integer, dimension(3,nli) :: indi
  integer, dimension(3,nlj) :: indj

  integer :: ii, ji, ij, jj, kc, lc
  integer :: ix, iy, iz, jx, jy, jz
  integer :: kx, ky, kz

!! order i : x,y,z |  xx, xy, xz, yy, yz, zz | xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz | ....
  
  mgto(:,:) = 0d0
! creates all possible (li+1)*(li+2)/2  polynomials for a given angular momentum li
  call genPolgto(li,nli,indi)
! creates all possible (lj+1)*(lj+2)/2  polynomials for a given angular momentum lj
  call genPolgto(lj,nlj,indj)

  do kc = 1, nli
    ix = indi(1,kc)
    iy = indi(2,kc)
    iz = indi(3,kc)
    do lc = 1, nlj
      jx = indj(1,lc)
      jy = indj(2,lc)
      jz = indj(3,lc)

!       WRITE(*,*)ix,iy,iz,"----",jx,jy,jz
       do kx = 0, ix + jx 
         do ky = 0, iy + jy 
           do kz = 0, iz + jz 
             mgto(lc,kc) = mgto(lc,kc) + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
!             WRITE(*,'(3i3,100(f22.16,1X))')kx,ky,kz,fkx(kx,jx,ix),fky(ky,jy,iy),fkz(kz,jz,iz),mat(kz,ky,kx),normgtoi(ix,iy,iz),normgtoj(jx,jy,jz)
           enddo
         enddo
       enddo
!              WRITE(*,*)lc,kc,mgto(lc,kc)

    enddo
  enddo

  end subroutine sumfk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kinetic(bet,lmax, li, lj, nli, nlj, fkx, fky, fkz, normgtoi, normgtoj, mat, mgto)
  integer, intent(in) :: lmax, li, lj, nli, nlj
  double complex, dimension(0:lmax,0:lmax,0:lmax) :: mat
  double precision, dimension(0:lmax,0:lmax,0:lmax) :: fkx, fky, fkz, normgtoi, normgtoj
  double complex, dimension(nlj, nli) :: mgto
  double precision :: bet

  integer, dimension(3,nli) :: indi
  integer, dimension(3,nlj) :: indj

  double complex :: mx, my, mz, mxyz

  integer :: ii, ji, ij, jj, kc, lc
  integer :: ix, iy, iz, jx, jy, jz
  integer :: kx, ky, kz

  mgto(:,:) = 0d0

! creates all possible (li+1)*(li+2)/2  polynomials for a given angular momentum li
  call genPolgto(li,nli,indi)
! creates all possible (lj+1)*(lj+2)/2  polynomials for a given angular momentum lj
  call genPolgto(lj,nlj,indj)

  do kc = 1, nli
    ix = indi(1,kc)
    iy = indi(2,kc)
    iz = indi(3,kc)
    do lc = 1, nlj
      jx = indj(1,lc)
      jy = indj(2,lc)
      jz = indj(3,lc)

! first (common) term in d2/dx2, d2/dy2 and d2/dz2

              mxyz = 0d0
              do kx = 0, ix + jx
                do ky = 0, iy + jy
                  do kz = 0, iz + jz
                      mxyz = mxyz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) - 2d0*bet*mxyz*( (2d0*jx+1) + (2d0*jy+1) + (2d0*jz+1) )

! two other terms in d2/dx2
            if(jx>1) then
              mx = 0d0
              do kx = 0, ix + jx - 2
                do ky = 0, iy + jy
                  do kz = 0, iz + jz
                      mx = mx + fkx(kx,jx-2,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) + jx*(jx-1)*mx
            endif

              mx = 0d0
              do kx = 0, ix + jx + 2
                do ky = 0, iy + jy
                  do kz = 0, iz + jz
                      mx = mx + fkx(kx,jx+2,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) + 4d0*bet**2*mx

!two other terms in d2/dy2
            if(jy>1) then
              my = 0d0
              do kx = 0, ix + jx 
                do ky = 0, iy + jy - 2
                  do kz = 0, iz + jz
                      my = my + fkx(kx,jx,ix)*fky(ky,jy-2,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) + jy*(jy-1)*my
            endif

              my = 0d0
              do kx = 0, ix + jx 
                do ky = 0, iy + jy + 2
                  do kz = 0, iz + jz
                      my = my + fkx(kx,jx,ix)*fky(ky,jy+2,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) + 4d0*bet**2*my

!two other terms in d2/dz2
            if(jz>1) then
              mz = 0d0
              do kx = 0, ix + jx
                do ky = 0, iy + jy 
                  do kz = 0, iz + jz - 2
                      mz = mz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz-2,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) + jz*(jz-1)*mz
            endif

              mz = 0d0
              do kx = 0, ix + jx
                do ky = 0, iy + jy 
                  do kz = 0, iz + jz + 2
                      mz = mz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz+2,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              mgto(lc,kc) = mgto(lc,kc) + 4d0*bet**2*mz

!            WRITE(*,*)lc,kc,mgto(lc,kc)
    enddo
  enddo

  mgto(:,:) = -0.5d0*mgto(:,:)

  end subroutine kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printmat(l,m,n,mat)
  integer, intent(in) :: l, m, n
  double complex, dimension(0:n,0:m,0:l), intent(in) :: mat
  integer :: iprint, jprint, kprint 

    do iprint = 0, l
      do jprint = 0, m 
       do kprint = 0, n 
           write(*,'(3i3,100(f20.12,1X))')iprint, jprint, kprint, mat(kprint,jprint,iprint)
       enddo
     enddo
   enddo

  end subroutine printmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ComputeMat
