c $Log: ZMEAN.f,v $
c Revision 1.6  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.5  1998/12/10  00:55:51  ldr
c HBG changes to V5-1-21
c
c Revision 1.4  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.3  1993/10/19  10:31:26  ldr
c Add extra diagnostics for qcloud scheme.
c
c Revision 1.2  92/12/09  14:42:45  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  92/04/23  11:06:29  ldr
c Initial revision
c 
      real ztsk(lat,2)
      real zrgk(lat,2)
      real zsgk(lat,2)
      real zrtk(lat,2)
      real zsink(lat,2)
      real zsouk(lat,2)
      real zclfork(lat,2)
      real zalbk(lat,2)
      real zclhk(lat,2)
      real zclmk(lat,2)
      real zclck(lat,2)
      real zcllk(lat,2)
      real zclk(lat,2)
      real zevk(lat,2)
      real zrnk(lat,2)
      real zhfk(lat,2)
      real zsndk(lat,2)
      real zsidk(lat,2)
      real zrunofk(lat,2)
      real zsbalk(lat,2)
      real z_tmink(lat,2)
      real z_tmaxk(lat,2)
      real zfrmk(lat,2)
      real zflxik(lat,2)
      real zfrack(lat,2)
      real ztslk(lat,2)
      real zgwk(lat,2)
      real edifzk(lat,2,nl)
      real tzonk(lat,2,nl)
      real ezonk(lat,2,nl)
      real zpslk(lat,2)
      common/zmean/ztsk,zrgk,zrtk,zsgk,zsink,zsouk,
     &     zclfork,zalbk,zclhk,zclmk,zclck,
     &     zcllk,zclk,zevk,zrnk,zhfk,zsndk,
     &     zsidk,zrunofk,zsbalk,z_tmink,z_tmaxk,
     &     zfrmk,zflxik,zfrack,ztslk,zgwk,
     &     edifzk,tzonk,ezonk,zpslk

