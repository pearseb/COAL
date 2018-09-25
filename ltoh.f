c $Log: ltoh.f,v $
c Revision 1.3  1996/10/24 01:03:00  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.2  1996/03/21  03:18:54  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.1  1992/06/16  11:54:58  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/dicegrid in COMDICE.f
c                  latha, lathb - latitudinal indices bounding the
c                                 extent of ice (h-grid)
c
c              from arguments
c                  al - lsx grid(h grid with no wrap around)
c
c     In/Out:  from arguments
c                  ah - h grid array
c
      subroutine ltoh (al, ah)

c Copies lsx-grid (= h-grid with no wraparound) array al to h-grid
c array ah, sets ah's extra latitude points, and sets ah's
c longitudinal wraparound points.

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension al(plon,plat), ah(0:plon+1,0:plat+1)


c Copy to ah, including extra points at latha-1 and lathb+1 if not
c next to poles.

      do 5 j = max(latha-1,1), min (lathb+1,plat)
       do 6 i=1,plon
         ah(i,j) = al(i,j)
    6   continue
    5 continue


c If latha is next to south pole, set extra points ah(i,0).

      if (latha.eq.1) then
        do 10 i=1,plon
          ah(i,0) = al(i,1)
   10   continue
      endif

c If lathb is next to north pole, set extra points ah(i,plat+1).

      if (lathb.eq.plat) then
        do 15 i=1,plon
          ah(i,plat+1) = al(i,plat)
   15   continue
      endif


c Set longitudinal wraparound

      do 20 j=latha-1,lathb+1 
        ah(0,j) = ah(plon,j)
        ah(plon+1,j) = ah(1,j)
   20 continue

      return
      end
