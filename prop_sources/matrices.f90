module matrices

  integer :: nsizeCGto
  double complex, dimension(:,:), allocatable :: ovCGto, potCGto, kinCGto, tCGto

  integer :: ntotsta, ntotcgto
  double complex, dimension(:,:), allocatable :: ttcgtocoup, ttcgtoovl, tpcgtocoup, tpcgtoovl
  double complex, dimension(:,:), allocatable :: ppcgtocoup, ppcgtoovl, ptcgtocoup, ptcgtoovl
  double complex, dimension(:,:,:), allocatable :: mcgtocoup, mcgtoovl 
  double complex, dimension(:,:,:), allocatable :: mcoup, movl

end module matrices
