      subroutine mmtx(a,b,c,kt,jt,mt)
      real a(kt,jt),b(jt,mt),c(kt,mt)
      integer kt,jt,mt
      integer k,j,m

C**** real matrix a(kt,jt) * real matrix b(jt,mt) 
C****    = real matrix c(kt,mt)

c---- When using, make sure that a long vector length is kt
c---- and the short (or variable eg 1 to LW/2) length is mt
c---- to obtain best vectorization results.

      do m=1,mt
       do k=1,kt
        c(k,m)=0.0
       enddo
      enddo

c---- Note : The following will be replaced on the NEC
c---- with an inhouse simple matrix multiply upon compilation
      do m=1,mt
       do j=1,jt
        do k=1,kt
         c(k,m)=c(k,m)+a(k,j)*b(j,m)
        enddo
       enddo
      enddo

      return
      end
