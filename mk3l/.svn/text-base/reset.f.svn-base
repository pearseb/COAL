c $Log: reset.f,v $
c Revision 1.2  1992/12/09 14:44:24  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c Revision 1.1  91/02/22  16:38:04  ldr
c Initial release V3-0
c 
      subroutine reset(x1,x2,x3,x4,nl)
      real x1(nl), x2(nl), x3(nl), x4(nl)
      do 100 i=1,nl
          avg=.25*(x1(i)+x2(i)+x3(i)+x4(i))
          a1=.5*(x2(i)-x4(i))
          b1=.5*(x1(i)-x3(i))
          b2=.25*((x1(i)+x3(i))-(x2(i)+x4(i)))
          x1(i)=avg
          x2(i)=a1
          x3(i)=b1
          x4(i)=b2
 100  continue
      return
      end
