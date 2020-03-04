module general
implicit none

integer, parameter :: nYlmMax = 7
character(14), dimension(nYlmMax), parameter :: YlmBlock = ['! s functions&', '! p functions&', '! d functions&', '! f functions&', '! g functions&', '! h functions&', '! i functions&' ] 
integer, parameter :: lenmax = 200

double precision, dimension(-1:4*nYlmMax+3) :: fact, fact2
double precision, dimension(-1:2*nYlmMax,0:2*nYlmMax) ::  cij

double precision, parameter :: pi = dacos(-1d0)

 contains

 subroutine init_constants

 integer :: i, j

   fact(:)=1d0
  do i=1,2*nYlmMax
   fact(i)=fact(i-1)*i
  enddo

    fact2(:)=1d0
  do i=2,2*nYlmMax
    fact2(i)=i    
    do j=i-2,2,-2
     fact2(i)=fact2(i)*j
    enddo
  end do

  do i = 0, 2*nYlmMax
    do j = 0, i! 2*nYlmMax
          cij(i,j) = fact(i)/(fact(j)*fact(i-j))
    enddo
  enddo

 end subroutine init_constants

end module general
