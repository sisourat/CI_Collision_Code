program get2estamoints
use newtypes
use setup
use inputlib
use matrices
use centerlib
use cgto
use ComputeMat
use MatRep
use ComputeMat_tt
use general, only : lenmax
use diagolib
!use diagolib_2e

implicit none

  character(len=lenmax) :: input, option, suff, finputxml
  integer :: i, j, k, l, iy, jy, igto, jgto, im, jm, li, lj, nli, nlj, ic, jc, kc, lc

  integer :: ii,jj,kk,ll,aa,bb, imo, jmo, kmo, lmo, itrans, jtrans
!  integer*8 :: intindex, intde

  double precision, dimension(:,:), allocatable :: h1eMO,kinMO, potMO
  double precision, dimension(:,:), allocatable :: eigsta, eigsta_2
  double complex, dimension(:,:),  allocatable :: SrepCGto_2, mointtmp 

!  double precision, dimension(:,:,:,:), allocatable :: MOint
  double precision, dimension(:,:), allocatable :: avkin, avpot
  double precision, dimension(:), allocatable :: esta,temp_eigtp

!  double precision, dimension(4,4) :: hhhh
  double precision, dimension(:,:), allocatable :: eigsta_2e, eigsta_tp

  double precision, dimension(:,:), allocatable :: avkin_2e, avpot_2e, avrep_2e
  double precision, dimension(:), allocatable :: esta_2e,esta_tp
  double precision :: t1, t2, t11,t22
  integer, dimension(:,:), allocatable :: indi, indj
  integer, dimension(:), allocatable :: geom

  t1=OMP_get_wtime()

  call init_constants

! set up the calculations : directories, input, allocation, etc...
  call getarg(1,input)
  call getarg(2,option)

  suff = "_sta "
  call setdir(input,suff,option)

  finputxml='input.xml'
  call inputsta(finputxml)

  do i = 1, nCenters
    call read_cgto(Centers(i))
  enddo

! determine the total size of the matrices
  nsizeCGto =sum(Centers(:)%nsizeBlocks)
  allocate(OvCGto(nsizeCGto,nsizeCGto),potCGto(nsizeCGto,nsizeCGto),KinCGto(nsizeCGto,nsizeCGto),tCGto(nsizeCGto,nsizeCGto))

  allocate(Ov2eCGTO(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))
  allocate(t2eCGTO(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))
  allocate(pot2eCGTO(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))
  allocate(kin2eCGTO(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))
  allocate(repCGTO(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))

  allocate(SrepCGTO(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))
                

  allocate(t2eCGTO_2d(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))
  allocate(Ov2eCGTO_2d(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))
  allocate(pot2eCGTO_2d(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))
  allocate(kin2eCGTO_2d(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))
  allocate(repCGTO_2d(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))


  OvCGto(:,:) =  0d0
  PotCGto(:,:) = 0d0
  KinCGto(:,:) = 0d0
  tCGTo(:,:) =   0d0
  Ov2eCGto_2d(:,:) =  0d0
  kin2eCGto_2d(:,:) = 0d0
  pot2eCGto_2d(:,:) = 0d0
  repCGto_2d(:,:) =  0d0



  SrepCGTO(:,:,:,:)=0.d0
  Ov2eCGto(:,:,:,:) = 0d0
  Pot2eCGto(:,:,:,:) = 0d0
  Kin2eCGto(:,:,:,:) = 0d0
  t2eCGTo(:,:,:,:) = 0d0
  repCGTO(:,:,:,:)=0d0

! build the overlap, potential, kinetic matrices in the CGto basis
  allocate(geom(nCenters))

  geom(:) = 1

!---------compute-matrices-for-T-------------
  call ComputeStaMat(geom)

!-----two-electrons-repulsions-integrals-----

  call StaMatRep(geom)

! write(*,"(4f10.5)")real( SrepCGTO)

! solve HC = ESC
  allocate(esta(nsizeCGto),eigsta(nsizeCGto,nsizeCGto))
  allocate(avkin(nsizeCGto,nsizeCGto),avpot(nsizeCGto,nsizeCGto))

!------1e-orbitals--------------------------------------------------
  call eig(nsizeCgto,real(ovCGto),real(potCGTO(:,:))+real(kinCGTO(:,:)),esta,eigsta)



!-------compute-h1eAO.txt, h1eMO.txt, mocoef.txt, moint.txt, overlapAO.txt-------------------

! allocate(MOint(nsizeCGto,nsizeCGto,nsizeCGto,nsizeCGto))
 allocate(mointtmp(nsizeCGto**2,nsizeCGto**2), eigsta_2(nsizeCGto**2,nsizeCGto**2), SrepCGto_2(nsizeCGto**2,nsizeCGto**2))

 mointtmp=0d0
 eigsta_2=0d0
 SrepCGto_2=0d0


!-----If using eigenvectors from GAMESS, they should be read here. by Junwen------------
!   eigsta=0.d0
!   call system("cp ../eigsta .")
!   open(unit=2241,file='eigsta',status='old')

!   read(2241,*)        
!      read(2241,*)   
!   do i=1,nsizeCGto   
 !     read(2241,*) 
!     do j=1,nsizeCGto
!         read (2241,*)aa,bb,   eigsta(i,j)  
!         write(1212,*)  eigsta(i,j)     
!     enddo
!   enddo
!   close(2241)

!-----Transfer from 4d to 2d------------

 itrans=0
 do imo=1, nsizeCGto
   do jmo=1, nsizeCGto
      itrans=itrans+1; jtrans=0
      do kmo=1, nsizeCGto
        do lmo=1, nsizeCGto
           jtrans=jtrans+1

           eigsta_2(jtrans,itrans)   =   eigsta(imo,kmo) * eigsta(jmo,lmo)
           SrepCGto_2(jtrans,itrans) =   SrepCGto(imo,kmo,jmo,lmo)
         
        enddo
      enddo
   enddo
 enddo

!------here is the moints, if needed, lapack could be used here--by Junwen,Oct.15, 2019----

 mointtmp = matmul( transpose(eigsta_2), matmul( SrepCGto_2, eigsta_2)  )

!-----To save moints in 4d--------------------------------------------------------------
 open(unit=2100,file='moint.txt')

 itrans=0
 do imo=1, nsizeCGto
   do jmo=1, nsizeCGto
      itrans=itrans+1; jtrans=0
      do kmo=1, nsizeCGto
        do lmo=1, nsizeCGto
           jtrans=jtrans+1

             if((real( mointtmp(jtrans,itrans)) .ge. 1d-8) .or. (real(mointtmp(jtrans,itrans)) .le.-1d-8)   ) then

                intde=0d0
                call intindexnum(imo,jmo,kmo,lmo,intde)

                write(2100,"(i16,1x,f20.15,1x,i4,1x,i4,1x,i4,1x,i4,1x)") intde, real(mointtmp(jtrans,itrans)), imo,jmo,kmo,lmo

             endif

        enddo
      enddo
   enddo
 enddo
 close(2100)

 allocate(h1eMO(nsizeCGto,nsizeCGto), kinMO(nsizeCGto,nsizeCGto), potMO(nsizeCGto,nsizeCGto))

 
  MOint=0d0
  h1eMO=0d0
  kinMO=0d0
  potMO=0d0



 
  h1eMO(:,:)=matmul( eigsta,matmul( (real(potCGTO(:,:))+real(kinCGTO(:,:))) ,transpose(eigsta)) )
  kinMO(:,:)=matmul( eigsta,matmul(real(kinCGTO(:,:)),transpose(eigsta)) )
  potMO(:,:)=matmul( eigsta,matmul(real(potCGTO(:,:)),transpose(eigsta)) )


  open(unit=1030,file='h1eAO.txt')
  open(unit=1031,file='h1eMO.txt')
  open(unit=1032,file='kinMO.txt')
  open(unit=1033,file='potMO.txt')
  open(unit=1034,file='mocoef.txt')
  open(unit=1035,file='overlapAO.txt')



  do i=1,nsizeCGto
    do j=i,nsizeCGto
        write(1030,'(2i8,4X, f20.15)') j, i,  real(potCGTO(j,i))+ real(kinCGTO(j,i))
        write(1035,'(2i8,4X, f20.15)') j, i,  real(ovCGto(j,i)) 
!        write(*,*) j,i, real(ovCGto(j,i))
    enddo
  enddo

        write(1031,'(i4)')  sum(Centers(:)%nsizeBlocks)   
        write(1032,'(i4)')  sum(Centers(:)%nsizeBlocks)   
        write(1033,'(i4)')  sum(Centers(:)%nsizeBlocks)   
 
! write(*,*)  sum(Centers(:)%nsizeBlocks)

  do i=1,nsizeCGto
    do j=1,i
        write(1031,'(f20.15)')   h1eMO(j,i)       
        write(1032,'(f20.15)')   kinMO(j,i)    
        write(1033,'(f20.15)')   potMO(j,i)    
    enddo
  enddo

  write(1034,'(2i8)') nsizeCGto, nsizeCGto
  do i=1,nsizeCGto
    do j=1,nsizeCGto
        write(1034,'(2i8,4X, f20.15)') i, j,  eigsta(i,j)       
    enddo
  enddo

  

  close(1030)
  close(1031)
  close(1032)
  close(1033)
  close(1034)
  close(1035)


  t2=OMP_get_wtime()

  write(*,*) "########################################### "
  write(*,*) "We have got 2e-integrals! You could kill the code now!"
  write(*,*) "########################################### "
  write(*,*) "wall time: ", t2-t1
 


!---------compute-matrices-for-TT, moint.txt is computed here-------------

  call ComputeStaMat_tt(geom,eigsta)

!  write(*,'(49f10.5)') real(t2eCGTO)

  deallocate(geom)


!----------------------------------------------------------------------------------------


! solve HC = ESC for TT

  allocate(esta_2e(nsizeCGto*(nsizeCGto+1)/2),eigsta_2e(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))
  allocate(avkin_2e(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2),avpot_2e(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2))
  allocate(avrep_2e(nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2) )


  eigsta_2e(:,:) =real(repCGTO_2d)+real(pot2eCGTO_2d)+real(kin2eCGTO_2d)


  call eig(nsizeCGto*(nsizeCGto+1)/2,real(ov2eCGto_2d),eigsta_2e ,esta_2e,eigsta_2e)



  avkin(:,:) = matmul(eigsta,matmul(real(kinCGTO(:,:)),transpose(eigsta)))
  avpot(:,:) = matmul(eigsta,matmul(real(potCGTO(:,:)),transpose(eigsta)))
  avkin_2e(:,:) = matmul(eigsta_2e,matmul(real(kin2eCGTO_2d(:,:)),transpose(eigsta_2e)))
  avpot_2e(:,:) = matmul(eigsta_2e,matmul(real(pot2eCGTO_2d(:,:)),transpose(eigsta_2e)))

  avrep_2e(:,:) = matmul(eigsta_2e,matmul(real(repCGTO_2d(:,:)),transpose(eigsta_2e)))





!---print 2e-eigenvectors and eigenvalues-------------------------------


  open(unit=20,file='eigenvalues-t')
  open(unit=21,file='eigenvalues-tt')


   write(*,'(/,a,/)') "------T(1e) states------------"
   do i = 1, nsizeCGto
      write(*,'(a,i4,5f16.9)')"        State ",i,esta(i), avkin(i,i), avpot(i,i), avkin(i,i)+avpot(i,i), (2d0*avkin(i,i)+avpot(i,i))
      write(20,'(i3,5f20.15)')i,esta(i), avkin(i,i), avpot(i,i), avkin(i,i)+avpot(i,i), (2d0*avkin(i,i)+avpot(i,i))
   enddo
 
   write(*,'(/,a,/)') "------TT(2e) states------------"
   do i=1,nsizeCGto*(nsizeCGto+1)/2
      write(*,'(a,i8,f24.15)')"        State ",i,esta_2e(i)
!!, avkin_2e(i,i), avpot_2e(i,i),avrep_2e(i,i), avkin_2e(i,i)+avpot_2e(i,i)+avrep_2e(i,i)&
!!                                            &, (2d0*avkin_2e(i,i)+avpot_2e(i,i)+avrep_2e(i,i))
      write(21,'(i8,6f20.15)') i,esta_2e(i), avkin_2e(i,i), avpot_2e(i,i),avrep_2e(i,i), avkin_2e(i,i)+avpot_2e(i,i)+avrep_2e(i,i)&
                                             &, (2d0*avkin_2e(i,i)+avpot_2e(i,i)+avrep_2e(i,i))
   enddo

  close(20)
  close(21)


  open(unit=30,file='eigenvectors-t')
  open(unit=31,file='eigenvectors-tt')

  open(unit=40,file='PotMat-t')
  open(unit=41,file='PotMat-tt')

  write(30,'(a)')"# for each eigenvalue : ci, icenter, icgto, ix, iy, iz " 
  write(31,'(a)')"# for each eigenvalue : ci, icenter, icgto, ix, iy, iz, jcenter, jcgto, jx, jy, jz " 

  write(30,*)nsizeCGto, nsizeCGto
  write(31,*)nsizeCGto*(nsizeCGto+1)/2,nsizeCGto*(nsizeCGto+1)/2

  do l = 1, nsizeCGto
    k = 0
    write(30,*)l,esta(l)

    write(40,'(5000(f20.15,1X))')(avpot(j,l),j=1,nsizeCGto)

    do i = 1, nCenters
     do iy = 1, nYlmMax

      if(Centers(i)%block(iy)%ncgto==0) cycle
      li = iy - 1
      nli = (li+1)*(li+2)/2

      allocate(indi(3,nli))
      call genPolgto(li,nli,indi)
          do igto = 1, Centers(i)%block(iy)%ncgto
             do im = 1, nli
                   k = k + 1
                   write(30,'(f33.15,1X,5(i4,1X),f22.15)')eigsta(l,k), i, igto, indi(1,im), indi(2,im),indi(3,im),avpot(l,k)
             enddo
          enddo
      deallocate(indi)
     
     enddo
    enddo
   enddo

   close(30)
   close(40)
!-----correct the mark of eigenvectors for TT--by Junwen-Mar.1st,2019---------------
   do l = 1, nsizeCGto*(nsizeCGto+1)/2


     k = 0
     write(31,*)l,esta_2e(l)
     write(41,'(5000(f20.15,1X))') (avrep_2e(j,l)+avpot_2e(j,l),j=1,nsizeCGto*(nsizeCGto+1)/2)

    do i = 1, nCenters
     do iy = 1, nYlmMax

      if(Centers(i)%block(iy)%ncgto==0) cycle
      li = iy - 1
      nli = (li+1)*(li+2)/2

      do igto = 1, Centers(i)%block(iy)%ncgto
       do im = 1, nli

         do j = i, nCenters
          do jy = 1, nYlmMax

          if((i.eq.j) .and.(jy.lt.iy) ) cycle
 
          if(Centers(j)%block(jy)%ncgto==0) cycle
             lj = jy - 1
             nlj = (lj+1)*(lj+2)/2


             allocate(indi(3,nli))
             allocate(indj(3,nlj))

             call genPolgto(lj,nlj,indj)
             call genPolgto(li,nli,indi)

             do jgto = 1, Centers(j)%block(jy)%ncgto

             if((i.eq.j).and.(iy.eq.jy) .and.(jgto.lt.igto) ) cycle 

              do jm = 1, nlj
 
              if((i.eq.j).and.(iy.eq.jy) .and.(igto.eq.jgto) .and.(jm.lt.im) ) cycle 
       
                   k = k + 1              
                   write(31,'(f33.15,1X,5(i4,1X),5(i4,1X),f22.15)')eigsta_2e(l,k), i, igto, indi(1,im), indi(2,im),indi(3,im),&
                                                          j, jgto, indj(1,jm), indj(2,jm),indj(3,jm),avpot_2e(l,k)+avrep_2e(l,k)
              enddo
             enddo

             deallocate(indi)
             deallocate(indj)

           enddo
          enddo

   
         enddo
        enddo
      enddo
     enddo

    enddo

   close(31)
   close(41)

  t2=OMP_get_wtime()

  write(*,*) 
  write(*,*) "Wall time: ", t2-t1


  deallocate(mointtmp, eigsta_2, SrepCGto_2)


  deallocate(esta,eigsta)
  deallocate(OvCGto,PotCGto,KinCGto,tCGto)
  deallocate(avkin,avpot)

  deallocate(esta_2e,eigsta_2e,h1eMO,kinMO,potMO)
  deallocate(avkin_2e,avpot_2e, avrep_2e)
  deallocate(Ov2eCGto,Pot2eCGto,Kin2eCGto,t2eCGto,repCGto,SrepCGto)
  deallocate(Ov2eCGto_2d,Pot2eCGto_2d,Kin2eCGto_2d,repCGto_2d)

  call freecenters()

end program get2estamoints
