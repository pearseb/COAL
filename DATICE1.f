c $Log: DATICE1.f,v $
c Revision 1.6  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.5  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.4  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.3  1993/11/03  11:44:02  ldr
c Changes to V4-4-24l from HBG for coupled model
c
c Revision 1.2  93/02/03  12:44:31  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.1  92/12/10  10:26:56  ldr
c Initial revision
c 
      real rhocw1, antk, arck, diz, hoice(ln2,lat)
      common/datice/rhocw1,antk,arck,diz,hoice
