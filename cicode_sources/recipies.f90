module recipies

   double precision :: mprec

   contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE tqli(d,e,n,np,z)
      INTEGER n,np
      REAL*8 d(np),e(np),z(np,np)
!U    USES pythag
      INTEGER i,iter,k,l,m
      REAL*8 b,c,dd,f,g,p,r,s,pythag
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
          if(iter.eq.30)pause 'too many iterations in tqli'
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
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
!     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.
          goto 1
        endif
15    continue
      return
      END SUBROUTINE tqli
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine machin
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      real*8 half,one,two
      parameter ( half=0.5d0 , one=1.0d0 , two=2.0d0 )
      real*8 fouru,twou,u,halfu,t
!      common / small / fouru,twou,u
      halfu = half
 1    t = one+halfu
      if  ( t.le.one )  then
            u     = two*halfu
            twou  = two*u
            fouru = two*twou
      else
            halfu = half*halfu
            go to 1
      endif
      mprec = u
      return
      end subroutine machin

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

      SUBROUTINE eigsrt(d,n,np)
      INTEGER n,np
      REAL*8 d(np)
      INTEGER i,j,k
      REAL*8 p
      do 13 i=1,n-1
      k=i
      p=d(i)
       do 11 j=i+1,n
      if(d(j).ge.p)then
          k=j
       p=d(j)
      endif
11    continue
      if(k.ne.i)then
       d(k)=d(i)
       d(i)=p
      endif
13    continue
      return
      END SUBROUTINE eigsrt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE tridagi(a,b,c,r,u,n)
         INTEGER n,NMAX
         double complex a(n),b(n),c(n),r(n)
         double complex u(n)
         PARAMETER (NMAX=50000)
         INTEGER j
         double complex bet,gam(NMAX)

         if(b(1).eq.0.)then
         write(*,*) 'tridagi: rewrite equations'
         stop
         endif
         bet=b(1)
         u(1)=r(1)/bet
          do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.)then
           write(*,*) 'tridagi failed'
           stop
          endif
           u(j)=(r(j)-a(j)*u(j-1))/bet
        11 continue
           do 12 j=n-1,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
        12 continue
          return
        END subroutine tridagi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module recipies
