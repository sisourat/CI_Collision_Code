program Coll
use general
use colllib
use setup
use centerlib
use inputlib
use cgto
use matrices
use colldyn
use collint
implicit none

  character(len=lenmax) :: input, option, suff, finputxml, finput
  character(len=lenmax) :: fout
  character(len=lenmax) :: char_ib

  double precision :: spac
 
  integer :: ib, is, iunit, izgrid
  integer :: i, j, k

  double precision :: start, finish
  double precision :: ostart,oend,rang

  call init_constants

!$OMP PARALLEL PRIVATE(rang)
!$ rang = OMP_GET_THREAD_NUM ()
!$OMP END PARALLEL

! set up the calculations : directories, input, allocation, etc...
  call getarg(1,input)
  call getarg(2,option)

  suff = "_coll "
  call setdir(input,suff,option)
  
  finputxml='input.xml'
  call inputcoll(finputxml)
  call cpfiles()

! read target eigenvalues/eigenvectors
  write(*,*)" # Target data"
  finputxml='tinput.xml'
  call inputsta(finputxml)
  ntCenters = nCenters

  allocate(tCenters(ntCenters))

  tCenters = Centers
  deallocate(Centers)

  do i = 1, ntCenters
    call read_cgto(tCenters(i))
  enddo
  
  ntcgto = sum(tCenters(:)%nsizeBlocks)
  allocate(testa(ntsta),teigvec(ntcgto,ntsta))

  finput='teigenvectors'
  call get_eig(ntsta,ntcgto,testa,teigvec,tstai,tstaf,finput)

! read projectile eigenvalues/eigenvectors
  write(*,*)" # Projectile data"
  finputxml='pinput.xml'
  call inputsta(finputxml)
  npCenters = nCenters

  allocate(pCenters(npCenters))
  pCenters = Centers
  deallocate(Centers)

  do i = 1, npCenters
    call read_cgto(pCenters(i))
  enddo

  npcgto = sum(pCenters(:)%nsizeBlocks)
  allocate(pesta(npsta),peigvec(npcgto,npsta))

  finput='peigenvectors'
  call get_eig(npsta,npcgto,pesta,peigvec,pstai,pstaf,finput)

!! print eigeval/vec in workdir 

  j = 0
  open(unit=20,file='sta')
  write(20,*)ntsta,npsta
  do i = 1, ntsta
   j = j + 1
   write(20,'(a,i3,f16.9,1X,a)')"        State ",j,testa(i),"t"
  enddo
  do i = 1, npsta
    j = j + 1
    write(20,'(a,i3,f16.9,1X,a)')"        State ",j,pesta(i),"p"
  enddo
  close(20)

! here comes the real stuffs

  vproj = v%vz
  if(v%vx/=0d0 .or. v%vy/=0d0) then
   write(*,*)"Vproj should be along z axis"
   stop 
  endif

  tgrid%na = zgrid%na
  allocate(zgrid%a(zgrid%na),tgrid%a(tgrid%na))
  if(zgrid%typ=='linear' .or. zgrid%typ=='Linear' .or. zgrid%typ=='LINEAR') then

    do i = 1, zgrid%na
      zgrid%a(i) = zgrid%amin + (i-1)*(zgrid%amax-zgrid%amin)/zgrid%na
      tgrid%a(i) = zgrid%a(i)/vproj
    enddo

  elseif(zgrid%typ=='exp' .or. zgrid%typ=='Exp' .or. zgrid%typ=='EXP') then

    spac=0d0 
    do i = 1, zgrid%na/2
     spac=spac+1.1**i
    enddo
     spac=(zgrid%amax-zgrid%amin)/(2d0*spac)

    zgrid%a(zgrid%na/2)=0d0
     j = 0
    do i = zgrid%na/2+1, zgrid%na-1
      j = j + 1
      zgrid%a(i) = zgrid%a(i-1)+1.1**j*spac
      zgrid%a(zgrid%na-i) = -zgrid%a(i)
    enddo
      zgrid%a(zgrid%na) = zgrid%amax

    do i = 1, zgrid%na
      tgrid%a(i) = zgrid%a(i)/vproj
    enddo

   else
     write(*,*)"No other z grids implemented"
     stop
   endif

  nbproj = b%n
  allocate(bproj(nbproj))
  do i = 1, nbproj
   if(b%typ=='linear' .or. b%typ=='Linear' .or. b%typ=='LINEAR') then
    bproj(i) = b%bmin + (i-1)*(b%bmax-b%bmin)/b%n
   else
     write(*,*)"No other b grids implemented"
     stop
   endif
  enddo

  ntotsta = ntsta + npsta
  ntotcgto = ntcgto + npcgto

!$OMP PARALLEL 
  allocate(ttcgtocoup(ntcgto,ntcgto),ttcgtoovl(ntcgto,ntcgto))
  allocate(tpcgtocoup(npcgto,ntcgto),tpcgtoovl(npcgto,ntcgto))
  allocate(ptcgtocoup(npcgto,ntcgto),ptcgtoovl(npcgto,ntcgto))
  allocate(ppcgtocoup(npcgto,npcgto),ppcgtoovl(npcgto,npcgto))
!$OMP END PARALLEL

  allocate(mcoup(zgrid%na,ntotsta,ntotsta),movl(zgrid%na,ntotsta,ntotsta))
  allocate(mcgtocoup(zgrid%na,ntotcgto,ntotcgto),mcgtoovl(zgrid%na,ntotcgto,ntotcgto))
  allocate(matintrp(2*ntotsta,2*ntotsta))
  allocate(psi(ntotsta))

  write(*,*)"# Coll is running"
   
  do is = 1, ninitsta
   iunit = 10 + is
   fout = 'prob'//achar(is+48)
   open(iunit,file=fout)
   write(iunit,*)nbproj
  enddo

  do ib = 1, nbproj
 
    pCenters(1)%x(1) = bproj(ib) 
    pCenters(1)%y(1) = 0d0

    call cpu_time(start)
    ostart = omp_get_wtime()
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(mcoup,movl,mcgtocoup,mcgtoovl,zgrid,ntotsta,ntsta,npsta,ntcgto,npcgto,v,vproj,pCenters,tCenters) 
! SCHEDULE(DYNAMIC) 
   do izgrid = 1, zgrid%na
!    pCenters(1)%z(1) = zgrid%a(izgrid)  pCenters(1)%z(1) is replaced by  zgrid%a(izgrid) everywhere in collint for parall. purpose
    call ComputeCollMat(izgrid)
   enddo
!$OMP END PARALLEL DO
!!       do j = 1, ntcgto  
!!          write(4000+j,'(2(i4,1X),10000(f15.8,1X))')1,j,mcgtoovl(1,j,:)
!!       enddo
   call cpu_time(finish)
   oend = omp_get_wtime()
   WRITE(*,*)"#MEC time",finish-start,oend-ostart

!  do i=1,ntcgto
!    do j=1,ntcgto
!       write(*,*)j,i,mcgtoovl(1,j,i)
!     enddo
!   enddo
!
   call cpu_time(start)
   ostart = omp_get_wtime()
   call inttrans
   oend = omp_get_wtime()
   WRITE(*,*)"#INT Trans time",finish-start,oend-ostart

   write(char_ib,*)ib

   open(unit=10,file="matcoupb"//trim(adjustl(char_ib)),form='unformatted')
    write(10)vproj,bproj(ib)
    write(10)zgrid%na,ntotsta
    do i = 1, ntsta
      write(10)testa(i)
    enddo
    do i = 1, npsta
      write(10)pesta(i)+0.5d0*vproj**2
    enddo
    do i = 1, zgrid%na
      write(10)tgrid%a(i)
      do j = 1, ntotsta
       write(10)mcoup(i,j,:)
!!       write(2000+j,'(2(i4,1X),10000(f15.8,1X))')i,j,mcoup(i,j,:)
      enddo
    enddo
   close(10)

   open(unit=10,file="matovlb"//trim(adjustl(char_ib)),form='unformatted')
    write(10)vproj,bproj(ib)
    write(10)zgrid%na,ntotsta
    do i = 1, ntsta
      write(10)testa(i)
    enddo
    do i = 1, npsta
      write(10)pesta(i)+0.5d0*vproj**2
    enddo
    do i = 1, zgrid%na
      write(10)tgrid%a(i)
      do j = 1, ntotsta
       write(10)movl(i,j,:)
!       write(*,*)'s',i,j,movl(i,j,:)
!!       write(3000+j,'(2(i4,1X),10000(f15.8,1X))')i,j,movl(i,j,:)
      enddo
    enddo
   close(10)

   call cpu_time(start)
   ostart = omp_get_wtime()
   call initdyn
   oend = omp_get_wtime()
   WRITE(*,*)"#Init Dyn time",finish-start,oend-ostart

   open(unit=10,file="b"//trim(adjustl(char_ib)),form='unformatted')
    write(10)vproj,bproj(ib)
    write(10)zgrid%na,ntotsta
    do i = 1, ntsta
      write(10)testa(i)
    enddo
    do i = 1, npsta
      write(10)pesta(i)+0.5d0*vproj**2
    enddo
    do i = 1, zgrid%na
      write(10)tgrid%a(i)
      do j = 1, ntotsta
       write(10)mcoup(i,j,:)
      enddo
    enddo
   close(10)
 
! NO DYNAMICS, IT WILL BE DONE SEPARATELY
!   do is = 1, ninitsta
!    iunit = 10 + is
!    psi(:) = 0d0
!    psi(initsta(is)) = 1d0
!   
!    call cpu_time(start)
!    ostart = omp_get_wtime()
!    call dyn
!    call cpu_time(finish)
!    oend = omp_get_wtime()
!    WRITE(*,*)"#DYN time",finish-start,oend-ostart
!   
!    write(*,'(i3,500(f12.6,1X))')is, bproj(ib),(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2)
!    write(iunit,'(500(f12.6,1X))')bproj(ib), (cdabs(psi(i))**2,i=1,ntotsta),sum(cdabs(psi(:))**2)
!   enddo
!  
  enddo ! ib
  
!  do is = 1, ninitsta
!   iunit = 10 + is
!   close(iunit)
!  enddo

  deallocate(testa,teigvec,pesta,peigvec,bproj)
  deallocate(mcoup,movl,mcgtocoup,mcgtoovl,matintrp,psi)
  deallocate(ttcgtocoup,ttcgtoovl,tpcgtocoup,tpcgtoovl,ptcgtocoup,ptcgtoovl,ppcgtocoup,ppcgtoovl)
  call freeCenters()
  call freeColl()

end program Coll
