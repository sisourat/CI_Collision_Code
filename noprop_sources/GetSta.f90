program getsta
use newtypes
use setup
use inputlib
use matrices
use centerlib
use cgto
use ComputeMat
use general, only : lenmax
use diagolib
implicit none

  character(len=lenmax) :: input, option, suff, finputxml
  integer :: i, j, k, l, iy, igto, im, li, nli, jlen

  double precision, dimension(:,:), allocatable :: eigsta
  double precision, dimension(:), allocatable :: esta
  integer, dimension(:,:), allocatable :: indi
  integer, dimension(:), allocatable :: geom

  call init_constants

! set up the calculations : directories, input, allocation, etc...
  call getarg(1,input)
  call getarg(2,option)

  suff = "_sta "
  call setdir(input,suff,option)

  finputxml='input.xml'
  call inputsta(finputxml)

  do i = 1, nCenters
    jlen=index(Centers(i)%fbasis,' ')
    call system('cp '//Centers(i)%fbasis(1:jlen-1)//' .')
    call read_cgto(Centers(i))
  enddo

! determine the total size of the matrices
  nsizeCGto = sum(Centers(:)%nsizeBlocks)
  allocate(OvCGto(nsizeCGto,nsizeCGto),potCGto(nsizeCGto,nsizeCGto),KinCGto(nsizeCGto,nsizeCGto),tCGto(nsizeCGto,nsizeCGto))

  OvCGto(:,:) = 0d0
  PotCGto(:,:) = 0d0
  KinCGto(:,:) = 0d0
  tCGTo(:,:) = 0d0

! build the overlap, potential, kinetic matrices in the CGto basis
  allocate(geom(nCenters))

  geom(:) = 1
  call ComputeStaMat(geom)

  deallocate(geom)

! solve HC = ESC
  allocate(esta(nsizeCGto),eigsta(nsizeCGto,nsizeCGto))
  call eig(nsizeCgto,real(ovCGto),real(potCGTO(:,:))+real(kinCGTO(:,:)),esta,eigsta)

! check overlap
!  eigsta = transpose(eigsta)
!  WRITE(*,*)matmul(transpose(eigsta),matmul(ovCGto,eigsta))
!  STOP

! print eigenvalues and eigenvectors

  open(unit=20,file='eigenvalues')
      write(20,*)nsizeCGto
   do i = 1, nsizeCGto
      write(*,'(a,i3,f16.9)')"        State ",i,esta(i)
      write(20,*)i,esta(i)
   enddo
  close(20)

  open(unit=30,file='eigenvectors')
  write(30,*)"# for each eigenvalue : ci, icenter, icgto, ix, iy, iz " 
  write(30,*)nsizeCGto, nsizeCGto
  do l = 1, nsizeCGto
    k = 0
    write(30,*)l,esta(l)

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
                   write(30,'(f32.15,1X,5(i4,1X))')eigsta(l,k), i, igto, indi(1,im), indi(2,im),indi(3,im)
             enddo
          enddo
      deallocate(indi)
     
     enddo
    enddo
   enddo

  deallocate(esta,eigsta)
  deallocate(OvCGto,PotCGto,KinCGto,tCGto)

  call freecenters()

end program getsta
