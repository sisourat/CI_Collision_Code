module collint
!$ use OMP_LIB
use newtypes
use overlap_mod
use nuclear_mod
use fdn_mod
use cgto, only : normgto, genPolgto
use centerlib
use matrices
use colllib
implicit none
 

  contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine computegtomat(nli,li,xi,yi,zi,alpi,nlj,lj,xj,yj,zj,alpj,vx,vy,vz,npot,px,py,pz,pot,ovGto,potGto,kiniGto,kinjGto,kin)

! vx,vy,vz are supposed to be related to j

  character(len=1), intent(in) :: kin
  double precision, intent(in) :: xi, yi, zi, xj, yj, zj, alpi, alpj, vx, vy, vz
  integer, intent(in) :: nli, li, nlj, lj, npot
  double precision, dimension(npot), intent(in) :: px, py, pz
  type(potential), dimension(npot), intent(in) :: pot

  double complex, dimension(nlj,nli), intent(out) :: ovGto, potGto, kiniGto, kinjGto

  double complex, dimension(:,:,:), allocatable :: ovmat, potmat, pottmp
  double precision, dimension(:,:,:), allocatable ::  fkx, fky, fkz
  double precision, dimension(:,:,:), allocatable :: normGtoI, normGtoJ

  double precision :: xp, yp, zp
  double precision :: cx, cy, cz, bet, alp
  double precision :: Gam, rab2
  double complex   :: kfactor

  integer :: lmax, tmax, ipot, jpot, kpot

    Gam =  alpi + alpj
    rab2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
    xp  =  ( alpi*xi + alpj*xj ) / Gam ; yp  =  ( alpi*yi + alpj*yj ) / Gam ; zp  =  ( alpi*zi + alpj*zj ) / Gam
    kfactor = exp(-alpi*alpj*rab2/Gam)*exp(dcmplx(0d0,vx*xp+vy*yp+vz*zp))
    lmax = li+lj+2 ; tmax = 3*lmax ! +2 is for kinetic energy

! ovGto contains overlap matrix elements needed for all possible x**l*y**m*z**n*exp(-Gam*r**2)  
! potGto contains nuclear attraction matrix elements 

    allocate(ovmat(0:lmax,0:lmax,0:lmax))
    allocate(potmat(0:lmax,0:lmax,0:lmax),pottmp(0:lmax,0:lmax,0:lmax))

    allocate(normGtoI(0:lmax,0:lmax,0:lmax),normGtoJ(0:lmax,0:lmax,0:lmax))
    call normgto(lmax,alpi,normGtoI)
    call normgto(lmax,alpj,normGtoJ)

    call overlap_driver(lmax,Gam,vx,vy,vz,ovmat) 

! computes nuclear integrals over all centers
    potmat(:,:,:) = 0d0

    do ipot = 1, npot
! cx, cy, cz are shifted by -xp, -yp, -zp respectively
! rp is considered as the center of referential to compute the integrals
      cx = -px(ipot) + xp ; cy = -py(ipot) + yp ; cz = -pz(ipot) + zp
!nico China 2019 change ci in -ci to be consistent with fdn
      cx = -cx; cy = -cy; cz = -cz

      do jpot = 1, pot(ipot)%n
        bet = pot(ipot)%alp(jpot)
        pottmp(:,:,:) = 0d0

! no need to compute nuclear up to li+lj+2, li+lj is enough; save some cpu time
        call nuclear_driver(lmax-2,0,Gam,vx,vy,vz,bet,cx,cy,cz,tmax,pottmp(0:lmax-2,0:lmax-2,0:lmax-2))
        potmat(:,:,:) = potmat(:,:,:) - pot(ipot)%c(jpot)*pottmp(:,:,:)
      enddo ! jpot
    enddo ! ipot

    allocate(fkx(0:lmax,0:lmax,0:lmax),fky(0:lmax,0:lmax,0:lmax),fkz(0:lmax,0:lmax,0:lmax))

! fkz, fky and fkz contain the coeff. from Gaussian product rule 
    call fdn(xi-xp,xj-xp,lmax,fkx) 
    call fdn(yi-yp,yj-yp,lmax,fky) 
    call fdn(zi-zp,zj-zp,lmax,fkz) 

! now the matrices are computed according to the Gaussian product rule for all
! (nli+1)*(nli+2)/2 and (nlj+1)*(nlj+2)/2 cartesian gaussian functions

    call sumfk(lmax,li,lj,nli,nlj,fkx,fky,fkz,normGtoI,normGtoJ,ovmat,ovGto)
    call sumfk(lmax,li,lj,nli,nlj,fkx,fky,fkz,normGtoI,normGtoJ,potmat,potGto)

    ovGto(:,:) =  dconjg(ovGto(:,:))*kfactor
    potGto(:,:) = dconjg(potGto(:,:))*kfactor

    if(kin=='t' .or. kin=='T') then
       kiniGto(:,:) = 0d0
       kinjGto(:,:) = 0d0

      call kinetic(alpi,alpj,lmax,li,lj,nli,nlj,fkx,fky,fkz,normGtoI,normGtoJ,ovmat,kiniGto,kinjGto)
      kiniGto(:,:) = dconjg(kfactor)*kiniGto(:,:)
      kinjGto(:,:) = kfactor*dconjg(kinjGto(:,:))

    endif

    deallocate(ovmat,fkx,fky,fkz,normGtoI,normGtoJ,potmat,pottmp)

  end subroutine computegtomat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ComputeCollMat(izgrid)

  integer, intent(in) :: izgrid
  integer :: i, j, iy, jy, ist, jst, ifin, jfin, li, lj, nli, nlj
  integer :: igto, jgto, icgto, jcgto, im, jm, kicgto, kjcgto

  double precision :: vx, vy, vz

  integer :: npot
  double precision, dimension(:), allocatable :: px, py, pz
  type(potential), dimension(:), allocatable :: pot
  double complex, dimension(:,:), allocatable :: ovGto, potGto, kinTGto, kinPGto

  vx = v%vx
  vy = v%vy
  vz = v%vz

! T-T block

  ttcgtoovl(:,:) = 0d0
  ttcgtocoup(:,:) = 0d0

  npot = npCenters
  allocate(px(npot),py(npot),pz(npot),pot(npot))
  do i = 1, npot
   px(i) = pCenters(i)%x(1)
   py(i) = pCenters(i)%y(1)
   pz(i) = zgrid%a(izgrid) !pCenters(i)%z(1)
   pot(i) = pCenters(i)%pot
  enddo
!  write(*,'(2(i3),3(F15.6,1X))')OMP_GET_THREAD_NUM(),izgrid,vx,vy,vz

  do i = 1, ntCenters
   do iy = 1, nYlmMax
    if(tCenters(i)%block(iy)%ncgto==0) cycle
    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, ntCenters
      do jy = 1, nYlmMax
        if(tCenters(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

        allocate(ovGto(nlj,nli),potGto(nlj,nli),kinTGto(nlj,nli)) 
! kinTGto not used for the moment

        do igto = 1, tCenters(i)%block(iy)%ngto
          do jgto = 1, tCenters(j)%block(jy)%ngto

           ovGto(:,:) = 0d0; potGto(:,:) = 0d0 
           call computegtomat(nli,li,tCenters(i)%x(1),tCenters(i)%y(1),tCenters(i)%z(1),tCenters(i)%block(iy)%gto(igto), &
                              nlj,lj,tCenters(j)%x(1),tCenters(j)%y(1),tCenters(j)%z(1),tCenters(j)%block(jy)%gto(jgto), &
                              0d0,0d0,0d0,npot,px,py,pz,pot,ovGto,potGto,kinTGto,kinTGto,'f')


!!           write(*,'(1000(f15.6,1X))')tCenters(i)%x(1),tCenters(i)%y(1),tCenters(i)%z(1),tCenters(i)%block(iy)%gto(igto)
!!           write(*,'(1000(f15.6,1X))')tCenters(j)%x(1),tCenters(j)%y(1),tCenters(j)%z(1),tCenters(j)%block(jy)%gto(jgto)
!!           write(*,'(1000(f15.6,1X))')px,py,pz
!!          do im = 1, nli
!!             do jm = 1, nlj
!!                write(*,'(2(i4,1X),1000(f15.6,1X))')jm,im,ovGto(jm,im),potGto(jm,im)
!!             enddo
!!           enddo
!!           write(*,*)
!!           write(*,*)
!!           write(*,*)

! finally the matrices are computed in the CGto basis

             kicgto = 0; kjcgto = 0
             do icgto = 1, tCenters(i)%block(iy)%ncgto
               if(tCenters(i)%block(iy)%cgto(icgto,igto)==0d0) cycle !nico090118
               do jcgto = 1, tCenters(j)%block(jy)%ncgto
                 if(tCenters(j)%block(jy)%cgto(jcgto,jgto)==0d0) cycle !nico090118

                do im = 1, nli
                 do jm = 1, nlj
                  kicgto = sum(tCenters(1:i-1)%nsizeBlocks) + sum(tCenters(i)%nsizePerBlock(1:iy-1)) + (icgto-1)*(nli) + im
                  kjcgto = sum(tCenters(1:j-1)%nsizeBlocks) + sum(tCenters(j)%nsizePerBlock(1:jy-1)) + (jcgto-1)*(nlj) + jm
                  ttcgtoovl(kjcgto,kicgto) = ttcgtoovl(kjcgto,kicgto) + tCenters(i)%block(iy)%cgto(icgto,igto)*tCenters(j)%block(jy)%cgto(jcgto,jgto)*ovGto(jm,im)
                  ttcgtocoup(kjcgto,kicgto) = ttcgtocoup(kjcgto,kicgto) + tCenters(i)%block(iy)%cgto(icgto,igto)*tCenters(j)%block(jy)%cgto(jcgto,jgto)*potGto(jm,im)
                 enddo
                enddo

               enddo
             enddo

          enddo ! jgto
        enddo   ! igto

       deallocate(ovGto,potGto,kinTGto)

      enddo ! jYlm
     enddo  ! jCenters

    enddo   ! iYlm
  enddo     ! iCenters

  deallocate(px,py,pz,pot)

! renormalize the CGTO

!!$OMP PARALLEL DO COLLAPSE(2)  
!! SCHEDULE(DYNAMIC) 
       do kicgto = 1, ntcgto  
         do kjcgto = 1, ntcgto  
             mcgtoovl(izgrid,kjcgto,kicgto) = ttcgtoovl(kjcgto,kicgto)/sqrt(ttcgtoovl(kicgto,kicgto)*ttcgtoovl(kjcgto,kjcgto))
             mcgtocoup(izgrid,kjcgto,kicgto) = ttcgtocoup(kjcgto,kicgto)/sqrt(ttcgtoovl(kicgto,kicgto)*ttcgtoovl(kjcgto,kjcgto))
         enddo
       enddo
!!$OMP END PARALLEL DO

! P-P block

  ppcgtoovl(:,:) = 0d0
  ppcgtocoup(:,:) = 0d0

  npot = ntCenters
  allocate(px(npot),py(npot),pz(npot),pot(npot))
  do i = 1, npot
   px(i) = tCenters(i)%x(1)
   py(i) = tCenters(i)%y(1)
   pz(i) = tCenters(i)%z(1)
   pot(i) = tCenters(i)%pot
  enddo

  do i = 1, npCenters
   do iy = 1, nYlmMax
    if(pCenters(i)%block(iy)%ncgto==0) cycle
    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, npCenters
      do jy = 1, nYlmMax
        if(pCenters(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

        allocate(ovGto(nlj,nli),potGto(nlj,nli),kinPGto(nlj,nli)) 

        do igto = 1, pCenters(i)%block(iy)%ngto
          do jgto = 1, pCenters(j)%block(jy)%ngto

           ovGto(:,:) = 0d0; potGto(:,:) = 0d0 
           call computegtomat(nli,li,pCenters(i)%x(1),pCenters(i)%y(1),zgrid%a(izgrid),pCenters(i)%block(iy)%gto(igto), &
                              nlj,lj,pCenters(j)%x(1),pCenters(j)%y(1),zgrid%a(izgrid),pCenters(j)%block(jy)%gto(jgto), &
                              0d0,0d0,0d0,npot,px,py,pz,pot,ovGto,potGto,kinPGto,kinPGto,'f')

! finally the matrices are computed in the CGto basis
             kicgto = 0; kjcgto = 0
             do icgto = 1, pCenters(i)%block(iy)%ncgto
              if(pCenters(i)%block(iy)%cgto(icgto,igto)==0d0) cycle !nico090118
               do jcgto = 1, pCenters(j)%block(jy)%ncgto
                if(pCenters(j)%block(jy)%cgto(jcgto,jgto)==0d0) cycle !nico090118

                do im = 1, nli
                 do jm = 1, nlj
                  kicgto = sum(pCenters(1:i-1)%nsizeBlocks) + sum(pCenters(i)%nsizePerBlock(1:iy-1)) + (icgto-1)*(nli) + im
                  kjcgto = sum(pCenters(1:j-1)%nsizeBlocks) + sum(pCenters(j)%nsizePerBlock(1:jy-1)) + (jcgto-1)*(nlj) + jm
                  ppcgtoovl(kjcgto,kicgto) = ppcgtoovl(kjcgto,kicgto) + pCenters(i)%block(iy)%cgto(icgto,igto)*pCenters(j)%block(jy)%cgto(jcgto,jgto)*ovGto(jm,im)
                  ppcgtocoup(kjcgto,kicgto) = ppcgtocoup(kjcgto,kicgto) + pCenters(i)%block(iy)%cgto(icgto,igto)*pCenters(j)%block(jy)%cgto(jcgto,jgto)*potGto(jm,im)
                 enddo
                enddo

               enddo
             enddo

          enddo ! jgto
        enddo   ! igto

       deallocate(ovGto,potGto,kinPGto)

      enddo ! jYlm
     enddo  ! jCenters

    enddo   ! iYlm
  enddo     ! iCenters

  deallocate(px,py,pz,pot)

! renormalize the CGTO
!!$OMP PARALLEL DO COLLAPSE(2)  
!! SCHEDULE(DYNAMIC) 
       do kicgto = 1, npcgto  
         do kjcgto = 1, npcgto  
             mcgtoovl(izgrid,kjcgto+ntcgto,kicgto+ntcgto) = ppcgtoovl(kjcgto,kicgto)/sqrt(ppcgtoovl(kicgto,kicgto)*ppcgtoovl(kjcgto,kjcgto))
             mcgtocoup(izgrid,kjcgto+ntcgto,kicgto+ntcgto) = ppcgtocoup(kjcgto,kicgto)/sqrt(ppcgtoovl(kicgto,kicgto)*ppcgtoovl(kjcgto,kjcgto))
         enddo
       enddo
!!$OMP END PARALLEL DO

! T-P/P-T block
!  WRITE(*,*)"T-P/P-T block"

  tpcgtoovl(:,:) = 0d0
  tpcgtocoup(:,:) = 0d0
  ptcgtoovl(:,:) = 0d0
  ptcgtocoup(:,:) = 0d0

  npot = ntCenters+npCenters
  allocate(px(npot),py(npot),pz(npot),pot(npot))
  do i = 1, ntCenters
   px(i) = tCenters(i)%x(1)
   py(i) = tCenters(i)%y(1)
   pz(i) = tCenters(i)%z(1)
   pot(i) = tCenters(i)%pot
  enddo
  do i = 1, npCenters
   j = i + ntCenters
   px(j) = pCenters(i)%x(1)
   py(j) = pCenters(i)%y(1)
   pz(j) = zgrid%a(izgrid) !pCenters(i)%z(1)
   pot(j) = pCenters(i)%pot
  enddo

  do i = 1, ntCenters
   do iy = 1, nYlmMax
    if(tCenters(i)%block(iy)%ncgto==0) cycle
    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, npCenters
      do jy = 1, nYlmMax
        if(pCenters(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

        allocate(ovGto(nlj,nli),potGto(nlj,nli),kinTGto(nlj,nli),kinPGto(nlj,nli)) 

        do igto = 1, tCenters(i)%block(iy)%ngto
          do jgto = 1, pCenters(j)%block(jy)%ngto

           ovGto(:,:) = 0d0; potGto(:,:) = 0d0 
           call computegtomat(nli,li,tCenters(i)%x(1),tCenters(i)%y(1),tCenters(i)%z(1),tCenters(i)%block(iy)%gto(igto), &
                              nlj,lj,pCenters(j)%x(1),pCenters(j)%y(1),zgrid%a(izgrid),pCenters(j)%block(jy)%gto(jgto), &
                              vx,vy,vz,npot,px,py,pz,pot,ovGto,potGto,kinTGto,kinPGto,'t')

! finally the matrices are computed in the CGto basis
             kicgto = 0; kjcgto = 0
             do icgto = 1, tCenters(i)%block(iy)%ncgto
               if(tCenters(i)%block(iy)%cgto(icgto,igto)==0d0) cycle !nico090118
               do jcgto = 1, pCenters(j)%block(jy)%ncgto
                if(pCenters(j)%block(jy)%cgto(jcgto,jgto)==0d0) cycle !nico090118

                do im = 1, nli
                 do jm = 1, nlj
                  kicgto = sum(tCenters(1:i-1)%nsizeBlocks) + sum(tCenters(i)%nsizePerBlock(1:iy-1)) + (icgto-1)*(nli) + im
                  kjcgto = sum(pCenters(1:j-1)%nsizeBlocks) + sum(pCenters(j)%nsizePerBlock(1:jy-1)) + (jcgto-1)*(nlj) + jm

                  tpcgtocoup(kjcgto,kicgto) = tpcgtocoup(kjcgto,kicgto) + tCenters(i)%block(iy)%cgto(icgto,igto)*pCenters(j)%block(jy)%cgto(jcgto,jgto)*(potGto(jm,im)+kinPGto(jm,im))
                  ptcgtocoup(kjcgto,kicgto) = ptcgtocoup(kjcgto,kicgto) + tCenters(i)%block(iy)%cgto(icgto,igto)*pCenters(j)%block(jy)%cgto(jcgto,jgto)*(dconjg(potGto(jm,im))+kinTGto(jm,im))
                  ptcgtoovl(kjcgto,kicgto) = ptcgtoovl(kjcgto,kicgto) + tCenters(i)%block(iy)%cgto(icgto,igto)*pCenters(j)%block(jy)%cgto(jcgto,jgto)*dconjg(ovGto(jm,im))

                 enddo
                enddo

               enddo
             enddo

          enddo ! jgto
        enddo   ! igto

       deallocate(ovGto,potGto,kinTGto,kinPGto)

      enddo ! jYlm
     enddo  ! jCenters

    enddo   ! iYlm
  enddo     ! iCenters

  deallocate(px,py,pz,pot)

! renormalize the CGTO
!!$OMP PARALLEL DO COLLAPSE(2)  
!! SCHEDULE(DYNAMIC) 
       do kicgto = 1, ntcgto  
         do kjcgto = 1, npcgto  
             mcgtoovl(izgrid,kjcgto+ntcgto,kicgto) = ptcgtoovl(kjcgto,kicgto)/sqrt(ttcgtoovl(kicgto,kicgto)*ppcgtoovl(kjcgto,kjcgto))
             mcgtocoup(izgrid,kjcgto+ntcgto,kicgto) = ptcgtocoup(kjcgto,kicgto)/sqrt(ttcgtoovl(kicgto,kicgto)*ppcgtoovl(kjcgto,kjcgto))
 
             mcgtoovl(izgrid,kicgto,kjcgto+ntcgto) = dconjg(mcgtoovl(izgrid,kjcgto+ntcgto,kicgto))
             mcgtocoup(izgrid,kicgto,kjcgto+ntcgto) = tpcgtocoup(kjcgto,kicgto)/sqrt(ttcgtoovl(kicgto,kicgto)*ppcgtoovl(kjcgto,kjcgto))

         enddo
       enddo
!!$OMP END PARALLEL DO
  end subroutine ComputeCollMat

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

  subroutine kinetic(alpi, alpj, lmax, li, lj, nli, nlj, fkx, fky, fkz, normgtoi, normgtoj, mat, migto, mjgto)

  integer, intent(in) :: lmax, li, lj, nli, nlj
  double complex, dimension(0:lmax,0:lmax,0:lmax) :: mat
  double precision, dimension(0:lmax,0:lmax,0:lmax) :: fkx, fky, fkz, normgtoi, normgtoj
  double complex, dimension(nlj, nli), intent(inout) :: migto, mjgto
  double precision, intent(in) :: alpi, alpj

  integer, dimension(3,nli) :: indi
  integer, dimension(3,nlj) :: indj

  double complex :: mix, miy, miz, mxyz
  double complex :: mjx, mjy, mjz

  integer :: ii, ji, ij, jj, kc, lc
  integer :: ix, iy, iz, jx, jy, jz
  integer :: kx, ky, kz

  migto(:,:) = 0d0
  mjgto(:,:) = 0d0

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
!              write(*,*)"1",mxyz
              migto(lc,kc) = migto(lc,kc) - 2d0*alpi*mxyz*( 2d0*(ix+iy+iz) + 3 )
              mjgto(lc,kc) = mjgto(lc,kc) - 2d0*alpj*mxyz*( 2d0*(jx+jy+jz) + 3 )

!two other terms in d2/dx2
              mix = 0d0; mjx = 0d0
              do kx = 0, ix + jx - 2
                do ky = 0, iy + jy
                  do kz = 0, iz + jz

                    if(jx>1) then
                      mjx = mjx + fkx(kx,jx-2,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                    endif
                    if(ix>1) then
                      mix = mix + fkx(kx,jx,ix-2)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                    endif

                  enddo
                enddo
              enddo

              migto(lc,kc) = migto(lc,kc) + ix*(ix-1)*mix
              mjgto(lc,kc) = mjgto(lc,kc) + jx*(jx-1)*mjx

              mix = 0d0; mjx = 0d0
              do kx = 0, ix + jx + 2
                do ky = 0, iy + jy
                  do kz = 0, iz + jz
                      mjx = mjx + fkx(kx,jx+2,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                      mix = mix + fkx(kx,jx,ix+2)*fky(ky,jy,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              migto(lc,kc) = migto(lc,kc) + 4d0*alpi**2*mix
              mjgto(lc,kc) = mjgto(lc,kc) + 4d0*alpj**2*mjx

!two other terms in d2/dy2
              miy = 0d0; mjy = 0d0
              do kx = 0, ix + jx 
                do ky = 0, iy + jy - 2
                  do kz = 0, iz + jz

                    if(jy>1) then
                      mjy = mjy + fkx(kx,jx,ix)*fky(ky,jy-2,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                    endif
                    if(iy>1) then
                      miy = miy + fkx(kx,jx,ix)*fky(ky,jy,iy-2)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                    endif

                  enddo
                enddo
              enddo

              migto(lc,kc) = migto(lc,kc) + iy*(iy-1)*miy
              mjgto(lc,kc) = mjgto(lc,kc) + jy*(jy-1)*mjy

              miy = 0d0; mjy = 0d0
              do kx = 0, ix + jx 
                do ky = 0, iy + jy + 2
                  do kz = 0, iz + jz
                      mjy = mjy + fkx(kx,jx,ix)*fky(ky,jy+2,iy)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                      miy = miy + fkx(kx,jx,ix)*fky(ky,jy,iy+2)*fkz(kz,jz,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              migto(lc,kc) = migto(lc,kc) + 4d0*alpi**2*miy
              mjgto(lc,kc) = mjgto(lc,kc) + 4d0*alpj**2*mjy

!two other terms in d2/dz2
              miz = 0d0; mjz = 0d0
              do kx = 0, ix + jx 
                do ky = 0, iy + jy 
                  do kz = 0, iz + jz - 2

                    if(jz>1) then
                      mjz = mjz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz-2,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                    endif
                    if(iz>1) then
                      miz = miz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz-2)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                    endif

                  enddo
                enddo
              enddo

              migto(lc,kc) = migto(lc,kc) + iz*(iz-1)*miz
              mjgto(lc,kc) = mjgto(lc,kc) + jz*(jz-1)*mjz

              miz = 0d0; mjz = 0d0
              do kx = 0, ix + jx 
                do ky = 0, iy + jy 
                  do kz = 0, iz + jz + 2
                      mjz = mjz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz+2,iz)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                      miz = miz + fkx(kx,jx,ix)*fky(ky,jy,iy)*fkz(kz,jz,iz+2)*mat(kz,ky,kx)*normgtoi(ix,iy,iz)*normgtoj(jx,jy,jz)
                  enddo
                enddo
              enddo
              migto(lc,kc) = migto(lc,kc) + 4d0*alpi**2*miz
              mjgto(lc,kc) = mjgto(lc,kc) + 4d0*alpj**2*mjz
!            WRITE(*,*)lc,kc,mgto(lc,kc)
    enddo
  enddo

  migto(:,:) = -0.5d0*migto(:,:)
  mjgto(:,:) = -0.5d0*mjgto(:,:)

  end subroutine kinetic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inttrans

double complex, dimension(:,:), allocatable :: mattmp
integer :: i, j, k, ii, jj

movl(:,:,:) = 0d0
mcoup(:,:,:) = 0d0

! Transforms T-T block
!$OMP PARALLEL PRIVATE(mattmp)

allocate(mattmp(ntcgto,ntsta))
!$OMP DO 
! SCHEDULE(DYNAMIC) 
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mcgtoovl(i,1:ntcgto,1:ntcgto),teigvec)
  movl(i,1:ntsta,1:ntsta) = matmul(transpose(teigvec),mattmp)

  mattmp(:,:) = matmul(mcgtocoup(i,1:ntcgto,1:ntcgto),teigvec)
  mcoup(i,1:ntsta,1:ntsta) = matmul(transpose(teigvec),mattmp)
enddo
!$OMP END DO
deallocate(mattmp)

! Transforms P-P block
allocate(mattmp(npcgto,npsta))
!$OMP DO 
! SCHEDULE(DYNAMIC) 
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mcgtoovl(i,ntcgto+1:ntotcgto,ntcgto+1:ntotcgto),peigvec)
  movl(i,ntsta+1:ntotsta,ntsta+1:ntotsta) = matmul(transpose(peigvec),mattmp)

  mattmp(:,:) = matmul(mcgtocoup(i,ntcgto+1:ntotcgto,ntcgto+1:ntotcgto),peigvec)
  mcoup(i,ntsta+1:ntotsta,ntsta+1:ntotsta) = matmul(transpose(peigvec),mattmp)
enddo
!$OMP END DO
deallocate(mattmp)

! Transforms T-P block
allocate(mattmp(ntcgto,npsta))
!$OMP DO 
! SCHEDULE(DYNAMIC) 
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mcgtoovl(i,1:ntcgto,ntcgto+1:ntotcgto),peigvec)
  movl(i,1:ntsta,ntsta+1:ntotsta) = matmul(transpose(teigvec),mattmp)

  mattmp(:,:) = matmul(mcgtocoup(i,1:ntcgto,ntcgto+1:ntotcgto),peigvec)
  mcoup(i,1:ntsta,ntsta+1:ntotsta) = matmul(transpose(teigvec),mattmp)

  do k = ntsta+1, ntotsta
    do j = 1, ntsta
       mcoup(i,j,k) = mcoup(i,j,k) - pesta(k-ntsta)*movl(i,j,k)
    enddo
  enddo

enddo
!$OMP END DO
deallocate(mattmp)

! Transforms P-T block
allocate(mattmp(npcgto,ntsta))
!$OMP DO 
! SCHEDULE(DYNAMIC) 
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mcgtoovl(i,ntcgto+1:ntotcgto,1:ntcgto),teigvec)
  movl(i,ntsta+1:ntotsta,1:ntsta) = matmul(transpose(peigvec),mattmp)

  mattmp(:,:) = matmul(mcgtocoup(i,ntcgto+1:ntotcgto,1:ntcgto),teigvec)
  mcoup(i,ntsta+1:ntotsta,1:ntsta) = matmul(transpose(peigvec),mattmp)

  do k = 1, ntsta
    do j = ntsta+1, ntotsta
       mcoup(i,j,k) = mcoup(i,j,k) - testa(k)*movl(i,j,k)
    enddo
  enddo

enddo
!$OMP END DO
deallocate(mattmp)

!$OMP END PARALLEL 

end subroutine inttrans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module collint
