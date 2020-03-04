program getpec
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
  integer :: i, j, k, l, iy, igto, im, li, nli

  double precision, dimension(:,:), allocatable :: eigsta
  double precision, dimension(:), allocatable :: esta
  integer, dimension(:,:), allocatable :: indi
  integer, dimension(:), allocatable :: geom

  call init_constants

! set up the calculations : directories, input, allocation, etc...
  call getarg(1,input)
  call getarg(2,option)

  suff = "_pec "
  call setdir(input,suff,option)

  finputxml='input.xml'
  call inputsta(finputxml)

  do i = 1, nCenters
    call read_cgto(Centers(i))
  enddo

! determine the total size of the matrices
  nsizeCGto = sum(Centers(:)%nsizeBlocks)
  allocate(OvCGto(nsizeCGto,nsizeCGto),potCGto(nsizeCGto,nsizeCGto),KinCGto(nsizeCGto,nsizeCGto),tCGto(nsizeCGto,nsizeCGto))
  allocate(esta(nsizeCGto),eigsta(nsizeCGto,nsizeCGto))


! build the overlap, potential, kinetic matrices in the CGto basis
  allocate(geom(nCenters))
  geom(:) = 1

   open(unit=30,file='pec') 

! moves only the last center
  do j = 1, Centers(nCenters)%npos
  
    geom(nCenters) = j

    OvCGto(:,:) = 0d0
    PotCGto(:,:) = 0d0
    KinCGto(:,:) = 0d0
    tCGTo(:,:) = 0d0

    call ComputeStaMat(geom)

! solve HC = ESC
    call eig(nsizeCgto,real(ovCGto),real(potCGTO(:,:))+real(kinCGTO(:,:)),esta,eigsta)
    write(30,'(3(f12.6,1X),5X,5000(f12.6))')Centers(nCenters)%x(j),Centers(nCenters)%y(j),Centers(nCenters)%z(j), esta(:)

  enddo

  close(30)

  write(*,*)"# GetPec ended normally"

  deallocate(geom)
  deallocate(esta,eigsta)
  deallocate(OvCGto,PotCGto,KinCGto,tCGto)

  call freecenters()

end program getpec
