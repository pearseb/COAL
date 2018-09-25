c $Log: RMGRID.f,v $
c Revision 1.5  2001/02/22 05:34:36  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.4  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.3  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/12/17  15:31:28  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  91/03/14  09:49:52  ldr
c Initial revision
c 
      real rmg (ln2,nl,lat)
      real rmmg(ln2,nl,lat)
      real rgt (ln2,nl,lat)

      common /rmgrid/ rmg, rmmg, rgt
