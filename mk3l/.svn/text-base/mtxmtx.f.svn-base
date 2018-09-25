c $Log: mtxmtx.f,v $
c Revision 1.4  1997/12/19 02:03:14  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.3  1992/12/09  14:44:00  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c Revision 1.2  91/03/13  12:59:04  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:42  ldr
c Initial release V3-0
c 
      subroutine mtxmtx(a,b,c,x)
      include 'PARAMS.f'
      real a(nl,nl),b(nl,nl),c(nl,nl)
      real sumu(nl)
C**** REAL MATRIX A* REAL MATRIX B* REAL CONSTANT X= REAL MATRIX C
      do 20 m=1,nl
      do k=1,nl
        sumu(k)=0.0
      enddo
      do j=1,nl
      do k=1,nl
        sumu(k)=sumu(k)+a(k,j)*b(j,m)
      enddo
      enddo
      do k=1,nl
        c(k,m)=sumu(k)*x
      enddo
   20 continue
      return
      end
