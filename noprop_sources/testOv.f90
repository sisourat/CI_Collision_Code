program DriverOv
use overlap_mod
implicit none

integer, parameter :: lmax = 1
double complex, dimension(0:lmax,0:lmax,0:lmax) :: ovmat
integer :: l, m, n
double precision :: alp, vx, vy, vz

alp=0.1d0
vx=0d0; vy=0d0; vz=0.2d0

call overlap_driver(lmax,alp,vx,vy,vz,ovmat)

do l = 0, lmax
  do m = 0, lmax
    do n = 0, lmax
      write(*,*)n,l,m,ovmat(n,l,m)
    enddo
  enddo
enddo

end program DriverOv
