module colllib
use general, only : lenmax
use newtypes, only : impactparameter, velocity, defgrid
use flib_dom
implicit none

 type(impactparameter) :: b
 type(velocity) :: v

 integer :: ninitsta
 integer, dimension(:), allocatable :: initsta
 
 character(len=lenmax) :: tinput, pinput

 integer :: ntsta, npsta, ntgto, npgto, ntcgto, npcgto
 integer :: tstai, tstaf, pstai, pstaf

 double precision, dimension(:), allocatable :: testa, pesta
 double precision, dimension(:,:), allocatable :: teigvec, peigvec

 integer :: nbproj
 double precision :: vproj
 double precision, dimension(:), allocatable :: bproj

 type(defgrid) :: zgrid, tgrid

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inputcoll(finput)
  character(len=lenmax) :: k, finput

  type(fnode), pointer     :: myInput, Node, InitNode, TargetNode,ProjNode, GridNode
  type(fnodeList), pointer :: List, ImpactVel, ImpactParam, InitList, TargetList, ProjList, GridList

  integer :: i

  myInput => parsefile(finput,verbose=.false.)

! Get the Collision parameters

  ImpactParam => getElementsByTagName(myInput,"ImpactParam")
  ImpactVel => getElementsByTagName(myInput,"ImpactVel")
  GridList => getElementsByTagName(myInput,"Zgrid")
  InitList => getElementsByTagName(myInput,"InitState")

! print the input info

  Node => item(ImpactParam,0)
  b%typ = getAttribute(Node,"type")
  k = getAttribute(Node,"bmin"); read (k,*) b%bmin
  k = getAttribute(Node,"bmax"); read (k,*) b%bmax
  k = getAttribute(Node,"nb"); read (k,*) b%n
  write(*,'(a,a10,2(f12.6,1X),i5)')'b grid ',  b%typ, b%bmin, b%bmax, b%n

  Node => item(ImpactVel,0)
  k = getAttribute(Node,"vx"); read (k,*) v%vx
  k = getAttribute(Node,"vy"); read (k,*) v%vy
  k = getAttribute(Node,"vz"); read (k,*) v%vz
  write(*,'(a,3(f12.6,1X))') 'vx, vy, vz', v

  GridNode => item(GridList,0)
  zgrid%typ = getAttribute(GridNode,"type")
  k = getAttribute(GridNode,"zmin"); read (k,*) zgrid%amin
  k = getAttribute(GridNode,"zmax"); read (k,*) zgrid%amax
  k = getAttribute(GridNode,"nzgrid"); read (k,*) zgrid%na
  write(*,'(a,2(f12.6,1X),i5)') 'zgrid', zgrid%amin, zgrid%amax, zgrid%na

  ninitsta = getLength(InitList)
  write(*,*)ninitsta
  allocate(initsta(ninitsta))
 
  do i = 0, ninitsta - 1
     InitNode => item(InitList,i)
     k = getAttribute(InitNode,"state"); read (k,*) initsta(i+1)
     write(*,*)"Init state", i, initsta(i+1)
  enddo

! Target info
   TargetList => getElementsByTagName(myInput,"Target")
   TargetNode => item(TargetList,0)
   tinput = getAttribute(TargetNode,"file")
   k = getAttribute(TargetNode,"stamin"); read (k,*)tstai
   k = getAttribute(TargetNode,"stamax"); read (k,*)tstaf
   ntsta = tstaf - tstai + 1
   if(tstaf<=0) ntsta=0

! Proj info
   ProjList => getElementsByTagName(myInput,"Projectile")
   ProjNode => item(ProjList,0)
   pinput = getAttribute(ProjNode,"file")
   k = getAttribute(ProjNode,"stamin"); read (k,*)pstai
   k = getAttribute(ProjNode,"stamax"); read (k,*)pstaf
   npsta = pstaf - pstai + 1
   if(pstaf<=0) npsta=0

  end subroutine inputcoll

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpfiles()
! cp target and projectile eigenvalues, eigenvectors in workdir

  character(lenmax) ::  finput
  character(lenmax) :: workdir

  logical :: dir_e, file_e
  integer :: ilen, jlen

! target input xml
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/input.xml'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/tinput.xml')

! target eigenvalues
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/eigenvalues'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/teigenvalues')

! target eigenvectors
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/eigenvectors'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/teigenvectors')

! projectile input xml
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/input.xml'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/pinput.xml')

! projectile eigenvalues
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/eigenvalues'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/peigenvalues')

! projectile eigenvectors
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/eigenvectors'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/peigenvectors')

  end subroutine cpfiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_eig(nsta,ncgto,esta,eigvec,ni,nf,finput)

  integer, intent(in) :: nsta, ncgto
  double precision, dimension(nsta), intent(inout) :: esta
  double precision, dimension(ncgto,nsta), intent(inout) :: eigvec
  integer, intent(in) :: ni, nf

  character(len=lenmax), intent(in) :: finput
  character(len=lenmax) :: a

  double precision :: etmp

  integer :: nv, ne
  integer :: iunit

  integer :: i, j, k, ista

  iunit = 10
  write(*,*)"# List of states"

  open(unit=iunit,file=finput)
   read(iunit,*)a
   read(iunit,*)nv, ne

   if(ne/=ncgto) then
     write(*,*)"error in get_eig input"
     write(*,*)"wrong numbers of ncgto in ", finput
     stop
   endif

   if(nv<nsta) then
     write(*,*)"error in get_eig input"
     write(*,*)"wrong numbers of nstates in ", finput
     stop
   endif
   
  
   k = 0
   do i = 1, nv
     read(iunit,*)ista, etmp
     if(ista>=ni .and. ista<=nf) then
       k = k + 1
       write(*,*)k,etmp
       esta(k) = etmp
       do j = 1, ne
         read(iunit,*)eigvec(j,k) 
       enddo
     else
       do j = 1, ne
         read(iunit,*)
       enddo
     endif
   enddo
   
  close(iunit)

end subroutine get_eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine freeColl()
   if(allocated(initsta)) deallocate(initsta)
  end subroutine freeColl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module colllib

