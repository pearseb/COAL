c $Log: icestat.f,v $
c Revision 1.6  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.5  1996/10/24 01:02:56  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1996/06/13  02:06:54  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.3  1993/11/03  11:44:17  ldr
c Changes to V4-4-24l from HBG for coupled model
c
c Revision 1.2  93/10/05  13:06:27  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.1  93/07/06  16:25:15  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from arguments
c                  ind
c
c     In/Out:  from common/aticstr in this subroutine
c                  ttaux, ttauy - ice stress accumulator arrays
c
c 
      subroutine icestat(ind)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ind

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      real ttaux,ttauy
      common/aticstr/ttaux(ln2,lat),ttauy(ln2,lat)

C Local work arrays and variables
      integer lg
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(ind.eq.1)then
c removed
      end if

      if((ind.eq.2).or.(ind.eq.3))then
c.... Set these ice stress accumulator arrays to zero at start of 
c.... model run, and reset to zero after the ice step.
        do 480 lg=1,lat
        do 480 mg=1,ln2
        ttaux(mg,lg)=0.0
  480   ttauy(mg,lg)=0.0
      end if

      return
      end
