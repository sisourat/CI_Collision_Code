program sigma
use general, only : pi
use newtypes
use gnufor2
implicit none

integer :: ntsta, npsta, nsta, initsta, nmax
double precision, dimension(:), allocatable :: esta
character(1), dimension(:), allocatable :: typsta

integer :: nbproj
double precision, dimension(:), allocatable :: bproj, sig, pcap, pexc, pion
double precision, dimension(:,:), allocatable :: prob

double precision :: scapt, sexc, stion, spion, sion
double precision, dimension(10) :: signs, signp, signd, signf, signg, signh
double precision :: enl, signl
integer :: ins, inp, ind, inf, ing, inh
integer :: istate

character(40) :: input, suff, option, a
character(1) :: again, ista

integer :: ilen, jlen
logical :: dir_e, file_e
integer :: i, j, k, opt

! set up the calculations : directories, input, allocation, etc...
call getarg(1,input)
call getarg(2,ista)

ilen=index(input,' ')
inquire( directory="./"//input(1:ilen-1)//"/.", exist=dir_e )
if ( dir_e .eqv. .false. ) then
 write(*,*) input(1:ilen-1), " does not exist"
 stop
endif

call chdir(input(:ilen-1))

open(unit=10,file='sta')
 read(10,*)ntsta,npsta
 nsta = ntsta + npsta
 allocate(esta(nsta),typsta(nsta),sig(nsta))
 do i = 1, nsta
   read(10,*)a,j,esta(i),typsta(i)
 enddo
close(10)

a='prob'//ista
open(unit=10,file=a)
 read(10,*)nbproj
 allocate(bproj(nbproj),prob(nsta,nbproj),pcap(nbproj),pexc(nbproj),pion(nbproj))
 do i = 1, nbproj
   read(10,*)bproj(i),(prob(j,i),j=1,nsta)
 enddo
close(10)

!! computes the cross sections
 initsta = 1

 scapt = 0d0
 sexc = 0d0
 stion = 0d0
 spion = 0d0

 pcap(:) = 0d0
 pexc(:) = 0d0
 pion(:) = 0d0

 sig(:) = 0d0
do i = 1, nsta
   sig(i) = 0.5d0*bproj(1)*bproj(1)*prob(i,1)
 do j = 1, nbproj-1
   sig(i) =  sig(i) + 0.5d0*(bproj(j+1)-bproj(j))*(bproj(j)*prob(i,j)+bproj(j+1)*prob(i,j+1))
 enddo
  sig(i) =  sig(i)*2d0*pi
  if(typsta(i) == "p" .and. esta(i) < 0d0) then
   scapt = scapt + sig(i)
   do j = 1, nbproj
    pcap(j) = pcap(j) + 2d0*pi*bproj(j)*prob(i,j)
   enddo 
  elseif(typsta(i) == "p" .and. esta(i) > 0d0) then
   spion = spion + sig(i)
   do j = 1, nbproj
    pion(j) = pion(j) + 2d0*pi*bproj(j)*prob(i,j)
   enddo 
  elseif(typsta(i) == "t" .and. esta(i) > 0d0) then
   stion = stion + sig(i)
   do j = 1, nbproj
    pion(j) = pion(j) + 2d0*pi*bproj(j)*prob(i,j)
   enddo 
  elseif(typsta(i) == "t" .and. esta(i) < 0d0 .and. i/= initsta) then
   sexc = sexc + sig(i)
   do j = 1, nbproj
    pexc(j) = pexc(j) + 2d0*pi*bproj(j)*prob(i,j)
   enddo 
  endif
enddo
sion = stion + spion

!write(*,'(3(f15.9,1X))')scapt/3.57d0,sexc/3.57d0,sion/3.57d0
!stop !NICO
i=1
j=1
k = 1
ins=1; inp=1; ind=1; inf=1; ing=1; inh=1
enl = esta(i)
signl = sig(i)
do while (j<16) 
 i=i+1
 if(esta(i)==enl) then
   k = k + 1
   signl = signl + sig(i)
!   write(*,*)i,j,esta(i),enl
 else
!   write(*,*)'nl',j, k
   if(k==1) then
     signs(ins) = signl/3.57d0
     ins=ins+1
   elseif(k==3) then
     signp(inp) = signl/3.57d0
     inp=inp+1
   elseif(k==5) then
     signd(ind) = signl/3.57d0
     ind=ind+1
   elseif(k==7) then
     signf(inf) = signl/3.57d0
     inf=inf+1
   elseif(k==9) then
     signg(ing) = signl/3.57d0
     ing=ing+1
   elseif(k==11) then
     signh(inh) = signl/3.57d0
     inh=inh+1
   endif
   j=j+1
   k = 1
   enl = esta(i)
   signl = sig(i)
 endif 
enddo
!write(*,*)sum(sig(6:14))/3.57d0
!stop
!write(*,'(300(f15.9,1X))')(sig(i)/3.57d0,i=1,89)
!write(*,'(300(f15.9,1X))')(sig(i)/3.57d0,i=ntsta+1,ntsta+89)
!write(*,'(300(i3,1X))')ins-1,inp-1,ind-1,inf-1,ing-1,inh-1
write(*,'(300(f15.9,1X))')(signs(i),i=1,ins-1), (signp(i),i=1,inp-1), (signd(i),i=1,ind-1), (signf(i),i=1,inf-1), (signg(i),i=1,ing-1), (signh(i),i=1,inh-1)
i=ntsta+1
j=1
k = 1
ins=1; inp=1; ind=1; inf=1; ing=1; inh=1
enl = esta(i)
signl = sig(i)
do while (j<16)
 i=i+1
 if(esta(i)==enl) then
   k = k + 1
   signl = signl + sig(i)
!   write(*,*)i,j,esta(i),enl
 else
!   write(*,*)'nl',j, k
   if(k==1) then
     signs(ins) = signl/3.57d0
     ins=ins+1
   elseif(k==3) then
     signp(inp) = signl/3.57d0
     inp=inp+1
   elseif(k==5) then
     signd(ind) = signl/3.57d0
     ind=ind+1
   elseif(k==7) then
     signf(inf) = signl/3.57d0
     inf=inf+1
   elseif(k==9) then
     signg(ing) = signl/3.57d0
     ing=ing+1
   elseif(k==11) then
     signh(inh) = signl/3.57d0
     inh=inh+1
   endif
   j=j+1
   k = 1
   enl = esta(i)
   signl = sig(i)
 endif
enddo
!write(*,'(300(f15.9,1X))')(sig(i)/3.57d0,i=1,89)
!write(*,'(300(f15.9,1X))')(sig(i)/3.57d0,i=ntsta+1,ntsta+89)
!write(*,'(300(i3,1X))')ins-1,inp-1,ind-1,inf-1,ing-1,inh-1
write(*,'(300(f15.9,1X))')(signs(i),i=1,ins-1), (signp(i),i=1,inp-1), (signd(i),i=1,ind-1), (signf(i),i=1,inf-1), (signg(i),i=1,ing-1), (signh(i),i=1,inh-1) 

do j = 1, nbproj
  write(101,'(10(f20.15,1X))')bproj(j),pcap(j),pexc(j),pion(j)
enddo
!stop
opt = 0
do while(opt/=-1)
  write(*,*)"Menu"
  write(*,*)"-1 exit "
  write(*,*)"0 print the energies"
  write(*,*)"1 print the total cross sections"
  write(*,*)"2 print the partial (shell-by-shell) cross sections"
  write(*,*)"3 print the partial (state-by-state) cross sections"
  write(*,*)"4 plot reduced prob. of the state"
  read(*,*)opt

  if(opt==0) then
     do i = 1, nsta
       write(*,*)i,esta(i),typsta(i)
     enddo
  elseif(opt==1) then
      write(*,*)"Cross sections in 10^-16 cm^2"
      write(*,*)"Capture", scapt/3.57d0
      write(*,*)"Excitation", sexc/3.57d0
      write(*,*)"Ionisation (t,p,tot)", stion/3.57d0, spion/3.57d0, sion/3.57d0
      do j = 1, nbproj
       write(101,'(10(f20.15,1X))')bproj(j),pcap(j),pexc(j),pion(j)
      enddo
  elseif(opt==2) then
      write(*,*)"Cross sections in 10^-16 cm^2"
      nmax = ntsta**0.5
      write(*,*)
      write(*,*)"On target"
      j=1
      do i = 1, nmax
       j=j+(i-1)**2
       write(*,'(i6,1X,2(f20.15,1X))')i,sum(sig(j:j+i**2-1))/3.57d0
      enddo
      nmax = npsta**0.5
      write(*,*)
      write(*,*)"On projectile"
       j=ntsta+1
      do i = 1, nmax
       j=j+(i-1)**2
!       write(*,*)j,j+i**2-1
       write(*,'(2(f10.5,1X))')sum(sig(j:j+i**2-1))/3.57d0
      enddo
      write(*,*)
  elseif(opt==3) then
      write(*,*)"Cross sections in 10^-16 cm^2"
      do i = 1, nsta
        write(*,'(i6,1X,a,2(f20.15,1X))')i,typsta(i),esta(i),sig(i)/3.57d0
      enddo
  elseif(opt==4) then
10  write(*,*)"State ?"
    read(*,*)istate
    write(*,*),'********************************************************************'
    call plot(bproj,2d0*pi*bproj*prob(istate,:),pause=0.0,persist='no')
    write(*,*)"Plotting another orbital? (type n to return to main menu)"
    read(*,*)again
    if(again/='n' .and. again/='N') goto 10
  endif

enddo


! set off the calculations :: free the memory
deallocate(esta,typsta,bproj,prob)
write(*,*)"# last plot in data_file.txt"
write(*,*)"# sigma ended normally"

end program sigma
