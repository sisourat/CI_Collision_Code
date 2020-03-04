module matrices

  integer :: nsizeCGto
  double complex, dimension(:,:), allocatable :: ovCGto, potCGto, kinCGto, tCGto
!$OMP THREADPRIVATE(ovCGto, potCGto, kinCGto, tCGto)  

  integer :: ntotsta, ntotcgto
  double complex, dimension(:,:), allocatable :: ttcgtocoup, ttcgtoovl, tpcgtocoup, tpcgtoovl
!$OMP THREADPRIVATE(ttcgtocoup,ttcgtoovl,tpcgtocoup,tpcgtoovl)    
  double complex, dimension(:,:), allocatable :: ppcgtocoup, ppcgtoovl, ptcgtocoup, ptcgtoovl
!$OMP THREADPRIVATE(ppcgtocoup, ppcgtoovl, ptcgtocoup, ptcgtoovl)  
  double complex, dimension(:,:,:), allocatable :: mcgtocoup, mcgtoovl 
  double complex, dimension(:,:,:), allocatable :: mcoup, movl

!!!$OMP THREADPRIVATE(ttcgtocoup,ttcgtoovl)

end module matrices
