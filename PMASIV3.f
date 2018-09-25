c $Log: PMASIV3.f,v $
c Revision 1.5  2001/02/22 05:34:46  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.4  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.3  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/07/22  09:22:03  mrd
c Removed saving of purely diagnostic radiative quanitities. Output
c diagnostics for these are now accumulated only on radiation steps.
c
c Revision 1.1  93/01/26  15:52:51  ldr
c Initial revision
c 
      real prgsav(ln2,lat)
      real psgsav(ln2,lat)
      real psgamp(ln2,lat)
      real phtksav(ln2,nl,lat)
      common /pmasiv3/ prgsav, psgsav, psgamp, phtksav
