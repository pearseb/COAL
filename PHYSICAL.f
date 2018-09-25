c $Log: PHYSICAL.f,v $
c Revision 1.7  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.6  1997/12/19  02:03:10  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.5  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1993/12/17  15:31:25  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  92/12/09  14:42:30  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.2  92/04/16  16:29:07  ldr
c Added new field alsf2 and replaced alsf by alsf1.
c 
c Revision 1.1  91/03/14  09:49:49  ldr
c Initial revision
c 
c If adding any new fields before ron(), ngrid should be increased too!

      real z4(ln2)
      real he(ln2)
      real tstar(ln2)
      real alsf2(ln2)
      real tstar1(ln2)
      real tstar2(ln2)
      real alsf1(ln2)
      real wg(ln2)
      real tb3(ln2)
      real wg2(ln2)
      real tb2(ln2)
      real snowd(ln2)
      real siced(ln2)
      real ron(ln2,nl)
      real son(ln2,nl)

      common /physical/ z4, he, tstar, alsf2, tstar1, tstar2, alsf1,
     &                 wg, tb3, wg2, tb2, snowd, siced, ron, son
      real griddata(ln2,ngrid)
      equivalence (griddata(1,1),z4(1))
