program Prop
use propdyn
use general, only : lenmax
implicit none

  character(len=lenmax) :: input, finput
  character(len=lenmax) :: fout

  double precision :: spac, vproj, bproj

  double precision :: t1, t2, tmec, tdyn, ttot
  integer :: i, j, k, ista

  call getarg(1,finput)
  call getarg(2,input)

  read(input,*)ista
  write(*,*)ista

  open(unit=10,file=finput,form='unformatted')
    read(10)vproj,bproj
    write(*,*)vproj,bproj
    read(10)ntime,ntotsta
    write(*,*)ntime,ntotsta
    allocate(mcoup(ntime,ntotsta,ntotsta),tgrid(ntime),esta(ntotsta),mat(ntotsta,ntotsta))

    do i = 1, ntotsta
     read(10)esta(i)
    enddo

    do i = 1, ntime
      read(10)tgrid(i)
!      write(*,*)tgrid(i)
      do j = 1, ntotsta
        read(10)mcoup(i,j,:)
!        mcoup(i,j,:)=-mcoup(i,j,:) 
!        write(*,*)i,j,mcoup(i,j,:)
!!         write(3000+j,'(2(i4,1X),10000(f15.8,1X))')i,j,mcoup(i,j,:)
      enddo
    enddo
  close(10)

!  mcoup(:,2,1)=-mcoup(:,1,2)
!  mcoup(:,1,2)=-mcoup(:,1,2)
!  mcoup(:,2,2)=1.15*mcoup(:,2,2)

  allocate(rmat2intrp(ntime,ntotsta,ntotsta),cmat2intrp(ntime,ntotsta,ntotsta))
  allocate(matintrp(ntotsta,ntotsta))
  allocate(psi(ntotsta))
  
  psi(:) = 0d0
  psi(ista) = 1d0
   
  call cpu_time(t1)
  call dyn
  call cpu_time(t2)
  tdyn = t2-t1
  write(*,*)'DYN takes',tdyn
  write(*,'(5000(f12.6,1X))')bproj,(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2),vproj
  write(5000,'(5000(f12.6,1X))')bproj,(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2),vproj

  deallocate(mcoup,matintrp,psi,rmat2intrp,cmat2intrp)

end program Prop
