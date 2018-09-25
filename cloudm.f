c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: cloudm.f,v $
c Revision 1.8  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.7  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.6  1998/12/10  00:55:36  ldr
c HBG changes to V5-1-21
c
c Revision 1.5  1997/12/17  23:22:46  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.4  1996/10/24  01:02:32  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.2  91/03/13  12:57:08  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:04  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT
c     Input:   from common/giant3 in GIANT3.f
c                  cl  - total cloud       clm - mid level cloud
c                  clh - high level cloud  cll - low level cloud
c
c              from arguments
c                  ns - hemisphere index   lg - latitude index
c
c     Output:  from common/chmap in CHMAP.f
c                  character arrays for plotting:
c                  chhic  - high cloud   chmidc - mid cloud
c                  chlowc - low cloud    chtotc - total cloud
c     
c 
      subroutine cloudm(lg)

      implicit none

!$OMP THREADPRIVATE ( /GIANT3/ )

C     TO GENERATE DIGITAL CLOUD MAPS

C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'GIANT3.f'

C Global data blocks
      include 'CHMAP.f'

C Local work arrays and variables
      integer icl
      integer mg

C Local data, functions etc
      character*1 cont(13)
      data cont/'0','1','2','3','4','5','6','7','8','9','t','t','t'/

C Start code : ----------------------------------------------------------

      do 579 mg=1,ln2
      if(cll(mg).gt.0.0)then
        icl=int(cll(mg)*10.0)+1
        chlowc(mg,lg)=cont(icl)
      end if
      if(clm(mg).gt.0.0)then
        icl=int(clm(mg)*10.0)+1
        chmidc(mg,lg)=cont(icl)
      end if
      if(clh(mg).gt.0.0)then
        icl=int(clh(mg)*10.0)+1
        chhic(mg,lg)=cont(icl)
      end if
      if(cl(mg).gt.0.0)then
        icl=int(cl(mg)*10.0)+1
        chtotc(mg,lg)=cont(icl)
      end if
  579 continue
      return
      end
