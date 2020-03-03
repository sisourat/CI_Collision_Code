
integer*8 function intindex(i1,j1,k1,l1)
implicit none
integer, intent(in) :: i1,j1,k1,l1
integer :: i, j, k, l, it, jt, kt, lt
integer*8 :: ij, kl

i=i1; j=j1; k=k1; l=l1

if(i<j) then
it = i
i = j
j = it
endif

if(k<l) then
kt = k
k = l
l = kt
endif

ij = i*(i+1)/2+j
kl = k*(k+1)/2+l

if(ij<kl) then
it = kl
kl = ij
ij = it
endif

intindex = ij*(ij+1)/2+kl

end function intindex

