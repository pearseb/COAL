c $Log: ltou.f,v $
c Revision 1.4  1996/10/24 01:03:01  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1996/03/21  03:18:55  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.2  1993/12/17  15:33:02  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  92/06/16  11:55:00  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/dicegrid in COMDICE.f
c                  latua, latub - latitudinal indices bounding the extent
c                                 of ice (u-grid)
c                  slat, slatu - sine latitude scales
c
c     Output:  from arguments
c                  au - u-grid array 
c
      subroutine ltou (al, au)

c Interpolates lsx-grid array al to u-grid array au, and sets au's
c wraparound points. Note that the interpolation is area weighted,
c by the use of the sine latitude scales slat,slatu.

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension  al(plon,plat),  au(0:plon+1, plat+1)
      dimension  ah(0:plon+1,0:plat+1)

      call ltoh (al, ah)

      do 10 j=latua,latub
        x1 = slat(j)  - slatu(j)
        x2 = slatu(j) - slat(j-1)
        w1 = 0.5*x1/(x1+x2)
        w2 = 0.5*x2/(x1+x2)
        do 12 i=1,plon+1
          au(i,j) = w1*ah(i,j)   + w1*ah(i-1,j)
     &            + w2*ah(i,j-1) + w2*ah(i-1,j-1)
   12   continue
   10 continue

      call wrapu (au)

      return
      end
