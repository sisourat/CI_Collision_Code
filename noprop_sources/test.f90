
   integer :: ib
   character(len=60) :: char_ib


   ib=1
   write(char_ib,*)ib
   write(*,*)ib,char_ib
   write(*,*)"matcoupb"//trim(adjustl(char_ib))
    end
