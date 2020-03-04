!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inttrans

double complex, dimension(:,:), allocatable :: mattmp
integer :: i, j, k, ii, jj

movl(:,:,:) = 0d0
mcoup(:,:,:) = 0d0

! Transforms T-T block
allocate(mattmp(ntgto,ntsta))
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mgtoovl(i,1:ntgto,1:ntgto),teigvec)
  movl(i,1:ntsta,1:ntsta) = matmul(transpose(teigvec),mattmp)

  mattmp(:,:) = matmul(mgtocoup(i,1:ntgto,1:ntgto),teigvec)
  mcoup(i,1:ntsta,1:ntsta) = matmul(transpose(teigvec),mattmp)
enddo
deallocate(mattmp)

! Transforms P-P block
allocate(mattmp(npgto,npsta))
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mgtoovl(i,ntgto+1:ntotgto,ntgto+1:ntotgto),peigvec)
  movl(i,ntsta+1:ntotsta,ntsta+1:ntotsta) = matmul(transpose(peigvec),mattmp)

  mattmp(:,:) = matmul(mgtocoup(i,ntgto+1:ntotgto,ntgto+1:ntotgto),peigvec)
  mcoup(i,ntsta+1:ntotsta,ntsta+1:ntotsta) = matmul(transpose(peigvec),mattmp)
enddo
deallocate(mattmp)

allocate(mattmp(ntgto,npsta))
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mgtoovl(i,1:ntgto,ntgto+1:ntotgto),peigvec)
  movl(i,1:ntsta,ntsta+1:ntotsta) = matmul(transpose(teigvec),mattmp)

  mattmp(:,:) = matmul(mgtocoup(i,1:ntgto,ntgto+1:ntotgto),peigvec)
  mcoup(i,1:ntsta,ntsta+1:ntotsta) = matmul(transpose(teigvec),mattmp)

  do k = ntsta+1, ntotsta
    do j = 1, ntsta
       mcoup(i,j,k) = mcoup(i,j,k) - pesta(k-ntsta)*movl(i,j,k)
    enddo
  enddo

enddo
deallocate(mattmp)

! Transforms P-T block
allocate(mattmp(npgto,ntsta))
do i = 1, zgrid%na
  mattmp(:,:) = matmul(mgtoovl(i,ntgto+1:ntotgto,1:ntgto),teigvec)
  movl(i,ntsta+1:ntotsta,1:ntsta) = matmul(transpose(peigvec),mattmp)

  mattmp(:,:) = matmul(mgtocoup(i,ntgto+1:ntotgto,1:ntgto),teigvec)
  mcoup(i,ntsta+1:ntotsta,1:ntsta) = matmul(transpose(peigvec),mattmp)

  do k = 1, ntsta
    do j = ntsta+1, ntotsta
       mcoup(i,j,k) = mcoup(i,j,k) - testa(k)*movl(i,j,k)
    enddo
  enddo

enddo
deallocate(mattmp)

!do i = 1, zgrid%na
!  write(50,'(10(f15.6,1X))')zgrid%a(i),mcoup(i,3,1),mcoup(i,1,3)
!  write(51,'(10(f15.6,1X))')zgrid%a(i),movl(i,3,1),movl(i,1,3)
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine inttrans
