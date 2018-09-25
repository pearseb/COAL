c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere, including a change of name from COUPLE1.f to A2O.f.
c SJP 2007/12/20
c
c $Log: COUPLE1.f,v $
c Revision 1.7  2000/06/20 02:08:34  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.6  1998/12/10 00:55:52  ldr
c HBG changes to V5-1-21
c
c Revision 1.5  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1993/11/29  11:38:19  ldr
c Changes to V4-4-32l from HBG for coupled model
c
c Revision 1.3  93/11/03  11:44:01  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.2  93/02/03  12:44:13  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.1  92/12/10  10:26:05  ldr
c Initial revision
c 
      real otx, oty, osalf, osurh, osrad,oswnd, osice
      common /a2o/ otx(imt-2, jmt-2), oty(imt-2, jmt-2),
     &             osalf(imt-2, jmt-2), osurh(imt-2, jmt-2)
     &           ,osrad(imt-2,jmt-2), oswnd(imt-2,jmt-2)
     &	          ,osice(imt-2,jmt-2)
