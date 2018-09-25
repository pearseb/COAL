c $Log: TRACEBLK.f,v $
c Revision 1.8  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.7  1994/08/08  14:17:14  ldr
c Implement memory-saving trick: NX for tracers and NC for qcloud.
c
c Revision 1.6  94/08/04  16:53:55  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.5  93/09/24  10:33:31  ldr
c Corrected memory saving trick.
c 
c Revision 1.4  93/08/10  12:15:42  ldr
c Re-introduce (modified) memory saving trick.
c 
c Revision 1.3  93/06/28  15:18:46  ldr
c Remove memory saving trick - caused compiler problems.
c 
c Revision 1.2  93/03/15  15:35:12  ldr
c Memory saving trick.
c 
c Revision 1.1  93/01/26  16:15:43  ldr
c Initial revision
c 
      integer ntrace, nx, lonx, lat2x
      parameter (ntrace=3)
      parameter (nx=0) !Memory saving trick - set NX=1 to use tracer.
      parameter (lonx=nx*(lon-1)+1, lat2x=nx*(lat2-1)+1) ! = lon, lat if nx=1
      real con, conm, source2, source3, conmnth
      integer nslfirst
      common/traceblk/con(-1:lonx+2,-1:lat2x+2,nl,ntrace),
     &                conm(-1:lonx+2,-1:lat2x+2,nl,ntrace),
     &                source2(lonx,lat2x),
     &                source3(lonx,lat2x),
     &                conmnth(lonx,lat2x,nl,ntrace),nslfirst
