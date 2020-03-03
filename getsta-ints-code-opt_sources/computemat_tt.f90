!--------------22.Feb.2017--by J.W.GAO.--follow-Nico's-one-electron-code--------------
module ComputeMat_TT
use newtypes
use repulsion_mod
use fdn_mod
use cgto, only : normgto, genPolgto
use centerlib
use matrices
use omp_lib
implicit none
 

  contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ComputeStaMat_TT(geom,eigsta)
    
  integer, dimension(*) :: geom
  double precision, dimension(nsizeCGto,nsizeCGto), intent(inout) ::  eigsta
  integer :: i, j, k, l, i1, j1, k1, l1, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll
  integer :: k2ij, k2kl, k2cij, k2ckl
  integer :: igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto
  integer :: itrans, jtrans, ktrans, ltrans
  integer :: imo, jmo, kmo, lmo
  integer*8 :: intde
  integer :: lmax, lmax1, lmax2, tmax

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

        do k = 1, nCenters
         do ky = 1, nYlmMax

          if(Centers(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nCenters
            do ly = 1, nYlmMax

             if(Centers(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2

! !$OMP PARALLEL PRIVATE(igto,jgto,kgto,lgto,kicgto,kjcgto,kkcgto,klcgto,itrans,jtrans,ktrans,ltrans,im,jm,km,lm,imo,jmo,kmo,lmo)


!!$OMP DO  COLLAPSE(2)  SCHEDULE(DYNAMIC)
 !$OMP DO  SCHEDULE(DYNAMIC)
 
                do igto = 1, Centers(i)%block(iy)%ngto
                 do jgto = 1, Centers(j)%block(jy)%ngto


                  do kgto = 1, Centers(k)%block(ky)%ngto
                   do lgto = 1, Centers(l)%block(ly)%ngto
                      


                      kicgto = 0; kjcgto = 0; kkcgto = 0; klcgto = 0

                      do icgto = 1, Centers(i)%block(iy)%ncgtoreal(igto)                           
                       do jcgto = 1, Centers(j)%block(jy)%ncgtoreal(jgto)
 

                        do kcgto = 1, Centers(k)%block(ky)%ncgtoreal(kgto)
                         do lcgto = 1, Centers(l)%block(ly)%ncgtoreal(lgto)  !---to remove the cgto with coefficient 0.0---



                           itrans =  Centers(i)%block(iy)%nincgto(icgto,igto)
                           jtrans =  Centers(j)%block(jy)%nincgto(jcgto,jgto) 
                           ktrans =  Centers(k)%block(ky)%nincgto(kcgto,kgto) 
                           ltrans =  Centers(l)%block(ly)%nincgto(lcgto,lgto) !---to remove the cgto with coefficient 0.0---

 

                            do im = 1, nli
                             do jm = 1, nlj

                              do km = 1, nlk
                               do lm = 1, nll


                                 kicgto = sum(Centers(1:i-1)%nsizeBlocks) + sum(Centers(i)%nsizePerBlock(1:iy-1)) + (itrans-1)*(nli) + im
                                 kjcgto = sum(Centers(1:j-1)%nsizeBlocks) + sum(Centers(j)%nsizePerBlock(1:jy-1)) + (jtrans-1)*(nlj) + jm
                                 kkcgto = sum(Centers(1:k-1)%nsizeBlocks) + sum(Centers(k)%nsizePerBlock(1:ky-1)) + (ktrans-1)*(nlk) + km
                                 klcgto = sum(Centers(1:l-1)%nsizeBlocks) + sum(Centers(l)%nsizePerBlock(1:ly-1)) + (ltrans-1)*(nll) + lm

!--------------here for two electrons, only singlet is considered------------------------------------------
!--------------repulsion-matrices---------------------------------------------------------------------------

                             
             
                                repCGTO(klcgto,kkcgto,kjcgto,kicgto) =  repCGTO(klcgto,kkcgto,kjcgto,kicgto) &
                                   & +( SrepCGto(klcgto,kkcgto,kjcgto,kicgto) + SrepCGto(kkcgto,klcgto,kjcgto,kicgto)&
                                   & +  SrepCGto(klcgto,kkcgto,kicgto,kjcgto) + SrepCGto(kkcgto,klcgto,kicgto,kjcgto)  ) &
                                   & *  Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & *  Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)

                                                                    
!-------------overlap-matrices----------------------------------------------------------------------------

                                t2eCGTO(klcgto,kkcgto,kjcgto,kicgto) =  t2eCGTO(klcgto,kkcgto,kjcgto,kicgto)&
                                   & +( (  OvCGTO(kkcgto,kicgto) * OvCGTO(klcgto,kjcgto) &
                                   & + OvCGTO(klcgto,kicgto) * OvCGTO(kkcgto,kjcgto)&
                                   & + OvCGTO(kkcgto,kjcgto) * OvCGTO(klcgto,kicgto)&
                                   & + OvCGTO(klcgto,kjcgto) * OvCGTO(kkcgto,kicgto)   )&
                                   & * Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & * Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto) )
                                

!-------------kinetic-matrices----------------------------------------------------------------------------
                                kin2eCGTO(klcgto,kkcgto,kjcgto,kicgto) =  kin2eCGTO(klcgto,kkcgto,kjcgto,kicgto)&
                                   & + (  KinCGTO(kkcgto,kicgto)* OvCGTO(klcgto,kjcgto) &
                                   & +    KinCGTO(klcgto,kicgto) * OvCGTO(kkcgto,kjcgto)&
                                   & +    KinCGTO(kkcgto,kjcgto) * OvCGTO(klcgto,kicgto)&
                                   & +    KinCGTO(klcgto,kjcgto) * OvCGTO(kkcgto,kicgto) &
                                   & +  OvCGTO(kkcgto,kicgto)* KinCGTO(klcgto,kjcgto) &
                                   & +  OvCGTO(klcgto,kicgto) * KinCGTO(kkcgto,kjcgto)&
                                   & +  OvCGTO(kkcgto,kjcgto) * KinCGTO(klcgto,kicgto)&
                                   & +  OvCGTO(klcgto,kjcgto) * KinCGTO(kkcgto,kicgto)  )&
                                   & * Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & * Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)
                                 
!-------------potential-matrices--------------------------------------------------------------------------
                                pot2eCGTO(klcgto,kkcgto,kjcgto,kicgto) =  pot2eCGTO(klcgto,kkcgto,kjcgto,kicgto)&
                                   & + (  PotCGTO(kkcgto,kicgto)* OvCGTO(klcgto,kjcgto) &
                                   & + PotCGTO(klcgto,kicgto) * OvCGTO(kkcgto,kjcgto)&
                                   & + PotCGTO(kkcgto,kjcgto) * OvCGTO(klcgto,kicgto)&
                                   & + PotCGTO(klcgto,kjcgto) * OvCGTO(kkcgto,kicgto) &
                                   & + OvCGTO(kkcgto,kicgto)  * PotCGTO(klcgto,kjcgto) &
                                   & + OvCGTO(klcgto,kicgto) * PotCGTO(kkcgto,kjcgto)&
                                   & + OvCGTO(kkcgto,kjcgto) * PotCGTO(klcgto,kicgto)&
                                   & + OvCGTO(klcgto,kjcgto) * PotCGTO(kkcgto,kicgto)   )&
                                   & * Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & * Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)





                               enddo
                              enddo
                             enddo
                            enddo


                         enddo
                        enddo
                       enddo
                      enddo


                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto

  !$OMP END DO 



            enddo  !ly
           enddo   !l
 
         enddo     !ky
        enddo      !k
 
      enddo        !jy
     enddo         !j

   enddo           !iy
  enddo            !i

!---------------transfer from 4-d to 2-d---------------------------------------

       k2ij=0

       do kicgto = 1, nsizeCGto  
         do kjcgto = kicgto, nsizeCGto
            
            k2ij=k2ij+1; k2kl=0

           do kkcgto = 1, nsizeCGto  
             do klcgto = kkcgto, nsizeCGto 

             k2kl=k2kl+1
 
  

             t2eCGTO_2d(k2kl,k2ij) = t2eCGTO(klcgto,kkcgto,kjcgto,kicgto)
             kin2eCGTO_2d(k2kl,k2ij) = kin2eCGTO(klcgto,kkcgto,kjcgto,kicgto)
             pot2eCGTO_2d(k2kl,k2ij) = pot2eCGTO(klcgto,kkcgto,kjcgto,kicgto)
             repCGTO_2d(k2kl,k2ij) = repCGTO(klcgto,kkcgto,kjcgto,kicgto)




             enddo
           enddo
         enddo
       enddo

!-------------------renormalize the CGTO----------------------------------------


        do kicgto = 1,(nsizeCGto+1)*nsizeCGto/2 
          do kjcgto = 1, (nsizeCGto+1)*nsizeCGto/2 

             Ov2eCGTO_2d(kjcgto,kicgto)  = t2eCGTO_2d(kjcgto,kicgto)  /sqrt(real(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto)))
             pot2eCGTO_2d(kjcgto,kicgto) = pot2eCGTO_2d(kjcgto,kicgto)/sqrt(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto))
             kin2eCGTO_2d(kjcgto,kicgto) = kin2eCGTO_2d(kjcgto,kicgto)/sqrt(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto))
             repCGTO_2d(kjcgto,kicgto)  = repCGTO_2d(kjcgto,kicgto)/sqrt(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto))

         enddo
       enddo



  end subroutine ComputeStaMat_TT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ComputeStaMat_TT_TRI(geom)
  
  integer, dimension(*) :: geom
  integer :: i, j, k, l, i1, j1, k1, l1, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll
  integer :: k2ij, k2kl, k2cij, k2ckl
  integer :: igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto
  integer :: itrans, jtrans, ktrans, ltrans
  
  integer :: lmax, lmax1, lmax2, tmax


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

        do k = 1, nCenters
         do ky = 1, nYlmMax

          if(Centers(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nCenters
            do ly = 1, nYlmMax

             if(Centers(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2


 
                do igto = 1, Centers(i)%block(iy)%ngto
                 do jgto = 1, Centers(j)%block(jy)%ngto


                  do kgto = 1, Centers(k)%block(ky)%ngto
                   do lgto = 1, Centers(l)%block(ly)%ngto
                      



                      kicgto = 0; kjcgto = 0; kkcgto = 0; klcgto = 0

                      do icgto = 1, Centers(i)%block(iy)%ncgtoreal(igto)                           
                       do jcgto = 1, Centers(j)%block(jy)%ncgtoreal(jgto)
 

                        do kcgto = 1, Centers(k)%block(ky)%ncgtoreal(kgto)
                         do lcgto = 1, Centers(l)%block(ly)%ncgtoreal(lgto)  !---to remove the cgto with coefficient 0.0---



                           itrans = Centers(i)%block(iy)%nincgto(icgto,igto)
                           jtrans = Centers(j)%block(jy)%nincgto(jcgto,jgto) 
                           ktrans = Centers(k)%block(ky)%nincgto(kcgto,kgto) 
                           ltrans = Centers(l)%block(ly)%nincgto(lcgto,lgto) !---to remove the cgto with coefficient 0.0---
 

                            do im = 1, nli
                             do jm = 1, nlj

                              do km = 1, nlk
                               do lm = 1, nll


                                kicgto = sum(Centers(1:i-1)%nsizeBlocks) + sum(Centers(i)%nsizePerBlock(1:iy-1)) + (itrans-1)*(nli) + im
                                kjcgto = sum(Centers(1:j-1)%nsizeBlocks) + sum(Centers(j)%nsizePerBlock(1:jy-1)) + (jtrans-1)*(nlj) + jm
                                kkcgto = sum(Centers(1:k-1)%nsizeBlocks) + sum(Centers(k)%nsizePerBlock(1:ky-1)) + (ktrans-1)*(nlk) + km
                                klcgto = sum(Centers(1:l-1)%nsizeBlocks) + sum(Centers(l)%nsizePerBlock(1:ly-1)) + (ltrans-1)*(nll) + lm
   

        

!--------------here for two electrons, only singlet is considered------------------------------------------
!--------------repulsion-matrices---------------------------------------------------------------------------

                             
             
                                repCGTO(klcgto,kkcgto,kjcgto,kicgto) =  repCGTO(klcgto,kkcgto,kjcgto,kicgto) &
                                   & +( SrepCGto(klcgto,kkcgto,kjcgto,kicgto) - SrepCGto(kkcgto,klcgto,kjcgto,kicgto)&
                                   & -  SrepCGto(klcgto,kkcgto,kicgto,kjcgto) + SrepCGto(kkcgto,klcgto,kicgto,kjcgto)  ) &
                                   & *  Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & *  Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)

                                                                    
!-------------overlap-matrices----------------------------------------------------------------------------

                                t2eCGTO(klcgto,kkcgto,kjcgto,kicgto) =  t2eCGTO(klcgto,kkcgto,kjcgto,kicgto)&
                                   & +( (  OvCGTO(kkcgto,kicgto) * OvCGTO(klcgto,kjcgto) &
                                   & - OvCGTO(klcgto,kicgto) * OvCGTO(kkcgto,kjcgto)&
                                   & - OvCGTO(kkcgto,kjcgto) * OvCGTO(klcgto,kicgto)&
                                   & + OvCGTO(klcgto,kjcgto) * OvCGTO(kkcgto,kicgto)   )&
                                   & * Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & * Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto) )
                                

!-------------kinetic-matrices----------------------------------------------------------------------------
                                kin2eCGTO(klcgto,kkcgto,kjcgto,kicgto) =  kin2eCGTO(klcgto,kkcgto,kjcgto,kicgto)&
                                   & + (  KinCGTO(kkcgto,kicgto)* OvCGTO(klcgto,kjcgto) &
                                   & -    KinCGTO(klcgto,kicgto) * OvCGTO(kkcgto,kjcgto)&
                                   & -    KinCGTO(kkcgto,kjcgto) * OvCGTO(klcgto,kicgto)&
                                   & +    KinCGTO(klcgto,kjcgto) * OvCGTO(kkcgto,kicgto) &
                                   & +  OvCGTO(kkcgto,kicgto)* KinCGTO(klcgto,kjcgto) &
                                   & -  OvCGTO(klcgto,kicgto) * KinCGTO(kkcgto,kjcgto)&
                                   & -  OvCGTO(kkcgto,kjcgto) * KinCGTO(klcgto,kicgto)&
                                   & +  OvCGTO(klcgto,kjcgto) * KinCGTO(kkcgto,kicgto)  )&
                                   & * Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & * Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)
                                 
!-------------potential-matrices--------------------------------------------------------------------------
                                pot2eCGTO(klcgto,kkcgto,kjcgto,kicgto) =  pot2eCGTO(klcgto,kkcgto,kjcgto,kicgto)&
                                   & + (  PotCGTO(kkcgto,kicgto)* OvCGTO(klcgto,kjcgto) &
                                   & - PotCGTO(klcgto,kicgto) * OvCGTO(kkcgto,kjcgto)&
                                   & - PotCGTO(kkcgto,kjcgto) * OvCGTO(klcgto,kicgto)&
                                   & + PotCGTO(klcgto,kjcgto) * OvCGTO(kkcgto,kicgto) &
                                   & + OvCGTO(kkcgto,kicgto)  * PotCGTO(klcgto,kjcgto) &
                                   & - OvCGTO(klcgto,kicgto) * PotCGTO(kkcgto,kjcgto)&
                                   & - OvCGTO(kkcgto,kjcgto) * PotCGTO(klcgto,kicgto)&
                                   & + OvCGTO(klcgto,kjcgto) * PotCGTO(kkcgto,kicgto)   )&
                                   & * Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & * Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)


                               enddo
                              enddo
                             enddo
                            enddo


                         enddo
                        enddo
                       enddo
                      enddo


                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto




            enddo  !ly
           enddo   !l
 
         enddo     !ky
        enddo      !k
 
      enddo        !jy
     enddo         !j

   enddo           !iy
  enddo            !i




!---------------complement the symetric matrix---------------------------

! do i1=1,kicgto
!   do j1=1,kjcgto
!     do k1=1,kkcgto
!       do l1=1,klcgto

!         if(repCGTO(l1,k1,j1,i1).eq.0.d0)then

!           repCGTO(l1,k1,j1,i1) = repCGTO(j1,i1,l1,k1) 

!         endif

!         if(t2eCGTO(l1,k1,j1,i1).eq.0.d0)then

!           t2eCGTO(l1,k1,j1,i1) = t2eCGTO(j1,i1,l1,k1) 

!         endif

!         if(kin2eCGTO(l1,k1,j1,i1).eq.0.d0)then

!           kin2eCGTO(l1,k1,j1,i1) = kin2eCGTO(j1,i1,l1,k1) 

!         endif

!         if(pot2eCGTO(l1,k1,j1,i1).eq.0.d0)then

!           pot2eCGTO(l1,k1,j1,i1) = pot2eCGTO(j1,i1,l1,k1) 

!         endif



!       enddo      
!     enddo         
!   enddo           
! enddo  



!---------------transfer from 4-d to 2-d---------------------------------------

       k2ij=0

       do kicgto = 1, nsizeCGto  
         do kjcgto = kicgto+1, nsizeCGto
            
            k2ij=k2ij+1; k2kl=0

           do kkcgto = 1, nsizeCGto  
             do klcgto = kkcgto+1, nsizeCGto 

             k2kl=k2kl+1
 
  

             t2eCGTO_2d(k2kl,k2ij) = t2eCGTO(klcgto,kkcgto,kjcgto,kicgto)
             kin2eCGTO_2d(k2kl,k2ij) = kin2eCGTO(klcgto,kkcgto,kjcgto,kicgto)
             pot2eCGTO_2d(k2kl,k2ij) = pot2eCGTO(klcgto,kkcgto,kjcgto,kicgto)
             repCGTO_2d(k2kl,k2ij) = repCGTO(klcgto,kkcgto,kjcgto,kicgto)




             enddo
           enddo
         enddo
       enddo

!-------------------renormalize the CGTO----------------------------------------


        do kicgto = 1,(nsizeCGto-1)*nsizeCGto/2 
          do kjcgto = 1, (nsizeCGto-1)*nsizeCGto/2 

             Ov2eCGTO_2d(kjcgto,kicgto)  = t2eCGTO_2d(kjcgto,kicgto)  /sqrt(real(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto)))
             pot2eCGTO_2d(kjcgto,kicgto) = pot2eCGTO_2d(kjcgto,kicgto)/sqrt(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto))
             kin2eCGTO_2d(kjcgto,kicgto) = kin2eCGTO_2d(kjcgto,kicgto)/sqrt(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto))
             repCGTO_2d(kjcgto,kicgto)  = repCGTO_2d(kjcgto,kicgto)/sqrt(t2eCGTO_2d(kicgto,kicgto)*t2eCGTO_2d(kjcgto,kjcgto))

         enddo
       enddo

  end subroutine ComputeStaMat_TT_TRI


subroutine intindexnum(i1,j1,k1,l1,inde)
implicit none
integer, intent(in) :: i1,j1,k1,l1
integer :: i, j, k, l, it, jt, kt, lt
integer*8 :: ij, kl,inde

i=i1; j=j1; k=k1; l=l1

if(i<j) then
it = i
i = j
j = it
endif

if(k<l) then
kt = k
k = l
l = kt
endif

ij = i*(i+1)/2+j
kl = k*(k+1)/2+l

if(ij<kl) then
it = kl
kl = ij
ij = it
endif

inde = ij*(ij+1)/2+kl

end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ComputeMat_TT
