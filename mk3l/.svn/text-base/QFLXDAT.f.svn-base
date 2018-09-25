c $Log: QFLXDAT.f,v $
c Revision 1.10  2001/02/22 05:34:36  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.9  2000/11/14 03:11:34  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.8  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.7  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.6  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.5  1994/03/30  12:33:39  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.4  94/01/17  10:38:15  ldr
c Changes from Hal to make qfluxes interpolate and be implicit. Also move
c read to datard so it's done every month.
c 
c Revision 1.3  94/01/04  17:18:28  ldr
c A little fix so that new qflux runs are initialized properly.
c 
c Revision 1.2  93/03/16  15:43:11  ldr
c Memory saver: Alter gm1cur (qfluxes) to cover one month only.
c 
c     Include file for the mixed layer ocean with qflux correction
c     Note the different order of the dimensions in itdiff. This is consistent
c     with other map type arrays used for output.

      real gm1cur(ln2,lat)
      real ochfa (ln2,lat)
      real gm2cur(ln2,lat)
      real gm3cur(ln2,lat)
      integer itdiff(ln2,lat)
      logical newqfluxrun
      common /qflxdat/ gm1cur, ochfa, itdiff, gm2cur, gm3cur
     & , newqfluxrun
