implicit none

integer, parameter :: ntime = 1000, nsta = 2
integer :: i
double precision :: time
double complex :: coup, diag

open(unit=10,file='testcoup.bin',form='unformatted')
write(10)ntime, nsta
write(10)-2.0d0
write(10)-1.0d0

do i = 1, ntime
 time = (i-1)*0.2d0
 diag = 0d0
 coup = 0.04d0*sin(time)

 write(10)time
 write(10) diag, coup
 write(10) coup, diag

! write(*,*)time,0.01d0*sin(time)
! write(*,*)-2.0, 0.0, 0.01d0*sin(time) 
! write(*,*)-1.0, 0.01d0*sin(time), 0.0 
enddo

end
