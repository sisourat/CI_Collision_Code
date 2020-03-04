      IMPLICIT real (A-h,o-z)
      real *8 ALPHAS(1000), C(1000)
      CHARACTER*1 ylm
      real *8 alps,alpm,tinys

       NBS=0
       ALPS=0.0
       TINYS=0.0

       WRITE(*,*)"NB, ALPMAX, ALPMIN, YLM ?"
       READ(*,*)NBS,ALPS,alpm,ylm
       TINYS=alpm/ALPS
       write(*,*)TINYS

       EPS= 1.0D0/DBLE(NBS-1)
       EPS=TINYS**(EPS)

       ALPHAS(1)=ALPS
 
       DO J = 2,NBS
          ALPHAS(J) = ALPHAS(J-1)*(EPS)
       ENDDO

          write(*,'(a,a,a)')"! ",ylm," functions"
          write(*,*)"H", NBS, NBS
       DO J = 1,NBS
          C(:)=0d0
          C(J)=1d0 
          write(*,'(500(f15.8,1X))')ALPHAS(J),(C(I),I=1,NBS)
       ENDDO

       END
