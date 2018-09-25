c $Log: BGRND.f,v $
c Revision 1.9  2001/02/22 05:34:36  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.8  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.7  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.6  1997/12/19  02:03:08  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.5  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1993/08/19  15:06:56  ldr
c Minor cosmetic changes.
c
c Revision 1.3  93/08/03  10:12:44  ldr
c Rearrange common  blocks to avoid misaligned data on SGI.
c 
c Revision 1.2  93/06/25  12:16:09  ldr
c Replace 4 with nbelow to generalize to N levels (HBG).
c 
c Revision 1.1  92/12/10  10:27:08  ldr
c Initial revision
c 
      real    tmonth(ln2,lat,nl),
     &        umonth(ln2,lat,nl),
     &        vmonth(ln2,lat,nl),
     &        qmonth(ln2,lat,nl),
     &        pmonth(ln2,lat),
     &        rmonth(ln2,lat,nl),
     &        gmonth(ln2,lat,nl)
      integer kdynm,
     &        ibeta(ln2,lat,nbelow)
      common/bgrnd/tmonth, umonth, vmonth, qmonth, pmonth,
     &             rmonth, gmonth,
     &             kdynm, ibeta

