c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c $Log: ATM2OC.f,v $
c Revision 1.4  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.3  1993/10/05  13:04:40  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.2  93/08/10  15:24:58  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.1  92/12/10  10:26:28  ldr
c Initial revision
c 
      integer nmth, ittx
      common/atm2oc/nmth,ittx
