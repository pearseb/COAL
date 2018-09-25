c $Log: wrapu.f,v $
c Revision 1.2  1996/03/21 03:19:13  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.1  1992/06/16  11:55:07  ldr
c Initial revision
c
      subroutine wrapu (au)

c Sets longitudinal wraparound points for u-grid array au

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension au(0:plon+1,plat+1)

      do 10 j=latua,latub
        au(0,j) = au(plon,j)
        au(plon+1,j) = au(1,j)
  10  continue

      return
      end
