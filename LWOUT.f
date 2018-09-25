c $Log: LWOUT.f,v $
c Revision 1.3  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/07/27  14:58:10  ldr
c Merge contents of common blocks clrtemp and lwoutclr into lwout.
c
c Revision 1.1  92/04/15  11:13:37  mrd
c Initial revision
c 
c 
c     common block lwout contains the quantities outputted by the 
c     longwave radiation code to the external module: 

      real heatra(imax,l) ! Heating rate at data levels (K/day) 
      real grnflx(imax)   ! Net longwave flux at the ground (CGS units) 
      real topflx(imax)   ! Net longwave flux at the top    (CGS units) 
      real grnflxclr(imax)
      real exctsclr(imax,l)
      real ctso3clr(imax,l)
c 
      common / lwout / heatra, grnflx, topflx,
     &                 grnflxclr, exctsclr, ctso3clr
