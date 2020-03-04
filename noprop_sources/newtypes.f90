!! created by Nico 23.07.2014

module newtypes
use general, only : nYlmMax, lenmax

!! new types for collision code

type impactparameter
 character(len=lenmax) :: typ
 double precision ::  bmin, bmax
 integer :: n
end type impactparameter

type velocity
 double precision :: vx, vy, vz
end type velocity

type defgrid
 character(len=lenmax) :: typ
 double precision :: amin, amax, da
 integer :: na
 double precision, dimension(:), allocatable :: a
end type defgrid

!! new types for getsta code

type cgtoblock
 integer :: ngto, ncgto
 double precision, dimension(:), allocatable :: gto
 double precision, dimension(:,:), allocatable :: cgto
end type cgtoblock

type potential
 integer :: n
 double precision, dimension(:), allocatable :: c, alp
end type potential

type center
 integer :: npos
 double precision, dimension(:), allocatable :: x, y, z
 type(potential) :: pot
 type(cgtoblock), dimension(nYlmMax) :: block
 integer :: nsizeBlocks
 integer, dimension(nYlmMax) :: nsizePerBlock
 character(len=lenmax) :: fbasis
end type center

end module newtypes
