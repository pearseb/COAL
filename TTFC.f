c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Removed declaration of redundant variable TFAC.
c SJP 2007/06/01
c
c $Log: TTFC.f,v $
c Revision 1.4  1997/12/19 02:03:12  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.3  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/10/05  13:05:20  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.1  92/12/11  11:10:06  ldr
c Initial revision
c 
      real dzlev1
      integer itfs
      common/ttfc/itfs,dzlev1
