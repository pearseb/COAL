c $Log: RADAV.f,v $
c Revision 1.12  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.11  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.10  1995/08/08  02:02:12  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.9  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.8  1994/08/04  10:34:05  ldr
c Split zonal mean cloud water diagnostic into separate liquid and frozen
c components.
c
c Revision 1.7  94/05/26  14:08:11  ldr
c Add new diagnostic for zonal mean cloud cover over sea.
c 
c Revision 1.6  94/03/30  10:17:26  ldr
c Add new zqlpath diagnostic (zonal mean liquid water path).
c 
c Revision 1.5  93/11/03  11:29:18  ldr
c Extra diagnostic for qcloud scheme.
c 
c Revision 1.4  93/10/19  10:31:24  ldr
c Add extra diagnostics for qcloud scheme.
c 
c Revision 1.3  92/12/09  14:42:32  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.2  91/05/15  11:12:50  ldr
c Modified for coupled model.
c 
c Revision 1.1  91/03/14  09:49:51  ldr
c Initial revisionb
c 
      integer nradav
      parameter (nradav=lat*2*(37+10*nl))
      real albz   (lat,2)
      real asbalz (lat,2)
      real atbalz (lat,2)
      real clcz   (lat,2)
      real clhz   (lat,2)
      real cllz   (lat,2)
      real clmz   (lat,2)
      real clz    (lat,2)
      real colwvz (lat,2)
      real evz    (lat,2)
      real flxicz (lat,2)
      real flxmlz (lat,2)
      real flxsez (lat,2)
      real flxmjz (lat,2)
      real rebalz (lat,2)
      real gwz    (lat,2)
      real hfz    (lat,2)
      real qlpz   (lat,2)
      real rgz    (lat,2)
      real rncz   (lat,2)
      real rnsz   (lat,2)
      real rnz    (lat,2)
      real rtclrz (lat,2)
      real rtz    (lat,2)
      real sblz   (lat,2)
      real sblz_l (lat,2)
      real sblz_i (lat,2)
      real sblz_s (lat,2)
      real sfbalz (lat,2)
      real sgz    (lat,2)
      real sinz   (lat,2)
      real souz   (lat,2)
      real souclrz(lat,2)
      real taz    (lat,2)
      real tslz   (lat,2)
      real tssz   (lat,2)
      real tsz    (lat,2)
      real cfracz (lat,2,nl)
      real htz    (lat,2,nl)
      real qaccrz (lat,2,nl)
      real qautoz (lat,2,nl)
      real qcollz (lat,2,nl)
      real qevapz (lat,2,nl)
      real qfz    (lat,2,nl)
      real qlz    (lat,2,nl)
      real qsublz (lat,2,nl)
      real rhz    (lat,2,nl)

      common/radav/albz, asbalz, atbalz, clcz, clhz, cllz, clmz
     &            ,clz, colwvz, evz
     &            ,flxicz, flxmlz, flxsez, flxmjz, rebalz
     &            ,gwz, hfz, qlpz, rgz, rncz, rnsz, rnz, rtclrz
     &            ,rtz, sblz, sblz_l, sblz_i, sblz_s
     &            ,sfbalz, sgz, sinz, souz, souclrz
     &            ,taz, tslz, tssz, tsz
     &            ,cfracz, htz, qaccrz, qautoz, qcollz, qevapz
     &            ,qfz, qlz, qsublz, rhz

      real clzx(nradav)
      equivalence (clzx(1),albz(1,1))
