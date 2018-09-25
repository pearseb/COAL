c $Log: RVARSG.f,v $
c Revision 1.6  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.5  1998/12/10  00:55:29  ldr
c HBG changes to V5-1-21
c
c Revision 1.4  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.3  1993/12/17  15:31:30  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  92/12/09  14:42:37  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  91/03/14  09:49:53  ldr
c Initial revision
c 
      real pg(ln2)
      real ttg(ln2,nl)
      real tscrn(ln2)
      real als(ln2)
      common /rvarsg/ pg, ttg, tscrn, als
