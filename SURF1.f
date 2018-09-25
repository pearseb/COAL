c $Log: SURF1.f,v $
c Revision 1.8  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.7  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.6  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.5  1997/12/23  00:34:48  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.4  1997/12/19  02:03:13  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.3  1996/12/05  03:42:32  ldr
c Fix to snow albedo from EAK.
c
c Revision 1.2  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1993/08/31  17:21:30  ldr
c Initial revision
c
      real tgg   (ln2,lat)
      real tgf   (ln2,lat)
      real mc    (ln2,lat)
      real osnowd(ln2,lat)
      real snage (ln2,lat)
      real ssdnn (ln2,lat)
      real gflux (ln2,lat)
      real sgflux(ln2,lat)

      real wb   (ln2,ms,lat)
      real wbice(ln2,ms,lat) 
      real tggsl(ln2,ms,lat)
      real tggsn(ln2,3,lat)
      real smass(ln2,3,lat)
      real ssdn3(ln2,3,lat)

      integer isflag(ln2,lat)
      common /surf1/ tgg, tgf, mc, osnowd, snage
     & , tggsl, tggsn, wb, wbice, smass, ssdn3, ssdnn, gflux, sgflux
     & , isflag

