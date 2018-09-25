c $Log: icebound.f,v $
c Revision 1.4  1996/10/24 01:02:53  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1996/03/21  03:18:49  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.2  1993/12/17  15:32:48  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  92/06/16  11:54:49  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from arguments
c                  fice - ice concentration, on h-grid
c                  ihem = 1 for southern hemisphere, 2 for northern hemisphere 
c
c     In/Out:  from common/dicegrid in COMDICE.f
c                  latha, lathb - h-grid)latitudinal indices bounding
c                  latua, latub - u-grid)the extent of ice 
c
      subroutine icebound (fice, ihem)


c For the current hemisphere, finds the latitudinal indices bounding
c the extent of ice. These latitudinal indices are latha,lathb for the
c h-grid, and latua,latub for the u-grid (in common dicegrid). If ice
c doesn't extend to the poles, the southernmost ice point is at latha+1,
c the northernmost ice point is at lathb-1, latua is the southern
c boundary of the ice (=latha+1), and latub is the northernmost boundary
c of the ice (=lathb). If ice does extend to the poles, then latha=1
c (89S) or lathb=plat(89N), and latua=1 (90S) or latub=plat+1 (90N).

c (Note that extra h-grid points are set in ltoh at latha-1 and lathb+1,
c which for the polar case are at j=0 (91S) or j=plat+1 (91N). So in all
c cases latua,latub are "inside" the region of defined h-grid points.)

c fice = ice concentration, on h-grid (supplied)
c ihem = 1 for southern hemisphere, 2 for northern hemisphere (supp)

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension fice(plon,plat)
      dimension ifany (plat)


      if (ihem.eq.1) then
        nla = 1
        nlb = plat/2
      else
        nla = plat/2 + 1
        nlb = plat
      endif


c First set ifany, flag array for whether ice exists at each latitude 

      do 100 j=nla,nlb
        ifany(j) = 0
        do 102 i=1,plon
          if (fice(i,j).gt.0.) then
            ifany(j) = 1
            goto 100
          endif
  102   continue
  100 continue


c Search northward from southernmost latitude, and set latha,latua.
c If no ice in this hemisphere, set latha=0 and return.

      do 200 j=nla,nlb
        if (ifany(j).eq.1) then
          latha = max (1, j-1)
          latua = latha + 1
          goto 202
        endif
  200 continue
      latha = 0
      return
  202 if (latha.eq.1) latua = 1


c Search southward from northernmost latitude, and set lathb,latub

      do 300 j=nlb,nla,-1
        if (ifany(j).eq.1) then
          lathb = min (plat, j+1)
          latub = lathb
          goto 302
        endif
  300 continue
  302 if (lathb.eq.plat) latub = plat+1

      if ( (ihem.eq.1 .and. lathb.eq.plat/2+1) .or.
     &     (ihem.eq.2 .and. latha.eq.plat/2-1) ) then
        write(ioterm,900) ihem
  900   format(' *** Warning. icebound 900. ice at equator. ihem=',i3)
      endif



      return
      end
