c $Log: WORKRID.f,v $
c Revision 1.3  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.2  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1993/10/06  11:25:39  ldr
c Initial revision
c
c Revision 1.2  92/12/09  14:42:19  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  91/03/14  09:49:43  ldr
c Initial revision
c 
      real ipr(lw,mw,nl),ipi(lw,mw,nl),ixr(lw,mw,nl),ixi(lw,mw,nl),
     &     itr(lw,mw,nl),iti(lw,mw,nl),elr(lw,mw,nl),eli(lw,mw,nl)
      common /workrid/ ipr, ipi, ixr, ixi, itr, iti, elr, eli

