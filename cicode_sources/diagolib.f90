module diagolib
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!double precision, private :: pythag

 contains

 subroutine eiglib(ndiag,ov,mat,eigval,eigvec)

 integer :: ndiag
 double precision, dimension(ndiag,ndiag), intent(in) :: ov, mat
 double precision, dimension(ndiag,ndiag), intent(inout) ::  eigvec
 double precision, dimension(ndiag,ndiag) :: mat1, mat2, mat3
 double precision, dimension(ndiag), intent(inout) :: eigval
 double precision, dimension(ndiag) :: d, e

! for LAPACK diagonalizer
 double precision, dimension(ndiag) ::  W
 double precision, dimension(3*ndiag) ::  WORK
 integer :: LWORK 

 double precision, dimension(:,:), allocatable :: U, t1, t2
 double precision, dimension(:), allocatable :: t3
 integer :: INFO, errflag
 logical, dimension(:), allocatable :: lhpsi
 integer, dimension(:), allocatable :: ihpsi
 double precision, dimension(:), allocatable :: rhpsi
 double complex, dimension(:), allocatable :: chpsi

!write(*,'(28f10.5)') ov


 LWORK = 3*ndiag 
 allocate(U(ndiag,ndiag))
 allocate(t1(ndiag,ndiag),t2(ndiag,ndiag),t3(8*ndiag))
 allocate(lhpsi(ndiag),ihpsi(ndiag),rhpsi(ndiag),chpsi(ndiag))




 mat1(:,:) = ov(:,:)
!write(*,*) ov
! computes S-1
! CALL DGETRF( ndiag, ndiag, mat1, ndiag, ihpsi, INFO )
! CALL DGETRI( ndiag, mat1, ndiag, ihpsi, WORK, LWORK, INFO )

! diagonalizes S-1H
! mat2(:,:) = matmul(mat1,mat)
! CALL DGEEV( 'N', 'V', ndiag, mat2, ndiag, d, e, t2, ndiag, U, ndiag, t3, 8*ndiag, INFO )
! write(*,*)"from S-1 H diag", d(:)
! stop

 mat1(:,:) = ov(:,:)
 mat2(:,:) = 0d0
 t2(:,:) = 0d0
 t3(:) = 0d0
 d(:) = 0d0
 e(:) = 0d0

! call tred2(mat1,ndiag,ndiag,d,e)
! call tqli2(d,e,ndiag,ndiag,mat1,0)

! U(:,:) = 0d0
! CALL DSYEV( 'V', 'L', ndiag, mat1, ndiag, W, WORK, LWORK, INFO )

! CALL DGEEV( 'N', 'V', ndiag, mat1, ndiag, d, e, t2, ndiag, U, ndiag, t3, 8*ndiag, INFO )

 call jacobi(mat1,ndiag,ndiag,d,U,INFO)
 mat1(:,:) = U(:,:)

 do i=1,ndiag
   do j=1,ndiag
     mat2(i,j)=mat1(j,i)
   enddo
 enddo

 do i=1,ndiag
   do j=1,ndiag
     mat1(j,i)=mat1(j,i)/dsqrt(d(i))
   enddo
 enddo

 call mulmat(mat1,mat2,mat3,ndiag,ndiag)
 call mulmat(mat3,mat,mat1,ndiag,ndiag)
 call mulmat(mat1,mat3,mat2,ndiag,ndiag)

! call tred2(mat2,ndiag,ndiag,d,e)
! call tqli2(d,e,ndiag,ndiag,mat2,0)
! eigval(:) = d(:)

 call jacobi(mat2,ndiag,ndiag,eigval,U,INFO)
 mat2(:,:) = U(:,:)

! CALL DSYEV( 'V', 'U', ndiag, mat1, ndiag, W, WORK, LWORK, INFO )
! eigval(:) = W(:)
 !CALL DGEEV( 'N', 'V', ndiag, mat2, ndiag, d, t1, t2, ndiag, U, ndiag, t3, 4*ndiag, INFO )
! eigval(:) = d(:)
! mat2(:,:) = U(:,:)

!  -> last transform (see Eric Cormier PhD)
 call mulmat(mat3,mat2,mat1,ndiag,ndiag)
 call eigsrt2(eigval,mat1,ndiag,ndiag)

!!--There is a transpose differece in GetSta-2e----
 eigvec(:,:) = (mat1(:,:)) 

! write(*,*)
! write(*,*)"from canonical ortho.",eigval
! write(*,*)
! write(*,*)eigvec(:,:)
! write(*,*)
! stop

!!write(*,*) "###################"
!!write(*,"(10f12.7)") eigvec(1:10,:)
  deallocate(U,t1,t2,t3)
  deallocate(lhpsi,rhpsi,ihpsi,chpsi)

 end subroutine eiglib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mulmat(a,b,c,n,np)
!
! multiply n*n matrices
!
      implicit real*8 (a-h,p-z)
      implicit integer (i-n)
      dimension a(np,np),b(np,np),c(np,np)
      do 30 i=1,n
       do 20 j=1,n
        c(i,j)=0
        do 15 k=1,n
         c(i,j)=c(i,j)+a(i,k)*b(k,j)
15      enddo
20     enddo
30    enddo
      return
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.
      FUNCTION pythag(a,b)
!      double precision :: pythag
      REAL*8 a,b
      REAL*8 absa,absb
      absa=dabs(a)
      absb=dabs(b)
      if(absa.gt.absb)then
        pythag=absa*dsqrt(1+(absb/absa)**2)
      else
        if(absb.eq.0)then
          pythag=0
        else
          pythag=absb*dsqrt(1+(absa/absb)**2)
        endif
      endif
      return
      END FUNCTION pythag
!  (C) Copr. 1986-92 Numerical Recipes Softwa

      SUBROUTINE tred2(a,n,np,d,e)
      INTEGER n,np
      REAL*8 a(np,np),d(np),e(np)
      INTEGER i,j,k,l
      REAL*8 f,g,h,hh,scale
      do 18 i=n,2,-1
        l=i-1
        h=0.
        scale=0.
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+abs(a(i,k))
11        continue
          if(scale.eq.0.)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(sqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.
            do 15 j=1,l
!     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
!     Omit following line if finding only eigenvalues.
      d(1)=0.
      e(1)=0.
      do 24 i=1,n
!!     Delete lines from here ...
        l=i-1
        if(d(i).ne.0.)then
          do 22 j=1,l
            g=0.
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
!     ... to here when finding only eigenvalues.
       d(i)=a(i,i)
!     Also delete lines from here ...
        a(i,i)=1.
        do 23 j=1,l
         a(i,j)=0.
          a(j,i)=0.
23      continue
!     ... to here when finding only eigenvalues.
24    continue
      return
      END SUBROUTINE tred2
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

      SUBROUTINE tqli2(d,e,n,np,z,ivap)
      INTEGER n,np
      REAL*8 d(np),e(np),z(np,np)
!U    USES pythag
      INTEGER i,iter,k,l,m
      REAL*8 b,c,dd,f,g,p,r,s
!      double precision :: pythag
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.40000) then
            print*, 'too many iterations in tqli',iter
            stop
          endif
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=pythag(g,1d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.)then
              d(i+1)=d(i+1)-p
              e(m)=0.
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
!     Omit lines from here ...
         if(ivap.ne.1) then
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
         endif
!     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.
          goto 1
        endif
15    continue
      return
      END SUBROUTINE tqli2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE eigsrt2(d,v,n,np)
      INTEGER n,np
      REAL*8 d(np),v(np,np)
      INTEGER i,j,k
      REAL*8 p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END SUBROUTINE eigsrt2
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE eigsrt3(d,v,w,n,np)
      INTEGER n,np
      REAL*8 d(np),v(np,np), w(np,np)
      INTEGER i,j,k
      REAL*8 p,q
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
            q=w(j,i)
            w(j,i)=w(j,k)
            w(j,k)=q
12        continue
        endif
13    continue
      return
      END SUBROUTINE eigsrt3
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=40000)
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+ &
     g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      write(*,*)'too many iterations in jacobi'
      d(:) = 1000d0
      return
      END subroutine jacobi
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

end module diagolib
