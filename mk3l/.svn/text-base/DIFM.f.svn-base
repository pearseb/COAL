c $Log: DIFM.f,v $
c Revision 1.5  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1993/12/17  15:31:20  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  93/10/05  13:04:41  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.2  92/12/09  14:42:08  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  91/03/14  09:49:40  ldr
c Initial revision
c 
      real frur(lw1,mw,nl),
     &     frui(lw1,mw,nl),
     &     frvr(lw1,mw,nl),
     &     frvi(lw1,mw,nl)
      common /difm/ frur, frui, frvr, frvi
      real dpsir(lw,mw,nl),
     &     dpsii(lw,mw,nl),
     &     dxhir(lw,mw,nl),
     &     dxhii(lw,mw,nl)
       equivalence (dpsir,frur),(dpsii,frui),(dxhir,frvr),(dxhii,frvi)
