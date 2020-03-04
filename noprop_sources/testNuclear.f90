program DriverNuclear
use f11mod
use nuclear_mod
implicit none


integer, parameter :: lmax = 4
double complex, dimension(0:lmax,0:lmax,0:lmax) :: potmat

integer :: l, m, n
double precision :: Gam, bet, rc, rc2, alp
double precision :: vx, vy, vz, cx, cy, cz

integer :: tmax
double complex :: nuclear

integer :: i, j

! compute int dx dy dz x**l y**m z**n * Exp(-Gam*r**2) * Exp(I*v.r) *
! Exp(-Bet*Rc**2)/Rc
Gam = 0.06d0
vx=0d0; vy=0d0; vz=0.2d0

cx=0.6d0; cy=0d0; cz=0.5d0
Bet=0d0

tmax=3*lmax

call nuclear_driver(lmax,0,Gam,vx,vy,vz,Bet,cx,cy,cz,tmax,potmat)

do l = 0, lmax
  do m = 0, lmax
    do n = 0, lmax
      write(*,'(3(i3,1X),3X,2(f30.16,1X))')n,l,m,potmat(n,l,m)
    enddo
  enddo
enddo

end program DriverNuclear
