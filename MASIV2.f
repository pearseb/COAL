c $Log: MASIV2.f,v $
c Revision 1.4  2001/02/22 05:34:44  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.3  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.2  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1992/12/10  10:24:52  ldr
c Initial revision
c
      real radstm
      common/masiv2/radstm(ln2,14:max_radstm,lat)
