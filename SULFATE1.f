c $Log: SULFATE1.f,v $
c Revision 1.6  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.5  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.4  1997/11/27  05:35:00  ldr
c Flexible treatment of direct and indirect aerosol effects.
c
c Revision 1.3  1997/11/24  23:25:26  ldr
c Indirect sulfates changes to give version used for GRL paper
c
c Revision 1.2  1997/10/03  05:45:47  ldr
c Changes for sulfates from LDR
c
c Revision 1.1  1997/06/11  02:51:54  ldr
c Initial revision
c
      real so4dir   !Sulfate used to calculate direct effect (in g/m**2)
      real so4rad   !Sulfate used for indirect effect in radiation (in g/m**2)
      real so4rain  !Sulfate used for indirect effect in rainfall (in g/m**2)
      common/sulfate1/so4dir(ln2,lat),so4rad(ln2,lat),
     &                so4rain(ln2,lat)

