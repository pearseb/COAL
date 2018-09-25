c $Log: polefilt.f,v $
c Revision 1.5  1996/10/24 01:03:08  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1996/03/21  03:18:58  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.3  1993/12/17  15:33:25  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  92/12/10  09:56:34  ldr
c Replace np's with nlp's for compatibility with ogcm.
c 
c Revision 1.1  92/06/16  11:55:02  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/dicegrid in COMDICE.f
c                  alatu - latitudes on u-grid
c                  clonu - cos longitudes on u-grid (alonu)
c                  slonu - sin longitudes on u-grid (alonu)
c
c              from arguments
c                  icode - see below
c
c     In/Out:  from arguments
c                  u1 - eastward  stress/current on u-grid (modified)
c                  v1 - northward stress/current on u-
c 
      subroutine polefilt (u1, v1, latslab1, icode)

c Spatially filters wind stresses or ocean currents near North Pole
c towards pure slab motion. 

c u        = eastward  stress/current on u-grid (modified)
c v        = northward stress/current on u-grid (modified)
c latslab1 = lat (deg) northward of which to apply filter (inclusive)
c icode    = 1 for seaice wind stresses, 
c            2 for open-ocean wind stresses, 
c            3 for ocean currents,
c            4 for ice velocities.

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension u1(0:plon+1,plat+1),  v1(0:plon+1,plat+1)

c     umaxpd = max zonal-mean stress/speed per deg from NLP
c     wslab1 = filtered vs unfiltered weight at latslab1, (=1 at NLP)

      dimension umaxpd(4), utmp(plon), vtmp(plon)
      
      data umaxpd /.02, .02, .01, .01/
      data wslab1 /0.5/


c Set latitudinal boundaries for the filter. If outside current
c hemispheric range of ice, do nothing.

c     jslab1 =max (latua, nint (plat*(latslab1+90.)/180. + 1.))
c     jslab2 = min (latub, plat+1)
      jslab1=plat
      jslab2=plat+1
      if (jslab1.gt.jslab2) return


c Calculate mean slab motion. Weight all pts equally, not by true area,
c since loosely speaking the more polar values are more "significant".

      uslab = 0.
      vslab = 0.

      do 100 j=jslab1,jslab2
        do 102 i=1,plon
          utmp(i) = -u1(i,j)*slonu(i) - v1(i,j)*clonu(i)
          vtmp(i) =  u1(i,j)*clonu(i) - v1(i,j)*slonu(i)
  102   continue
        do 104 i=1,plon                                        !masscomp
          uslab = uslab + utmp(i)                              !masscomp
          vslab = vslab + vtmp(i)                              !masscomp
  104   continue                                               !masscomp
cray    uslab = uslab + ssum (plon, utmp, 1)
cray    vslab = vslab + ssum (plon, vtmp, 1)
  100 continue

      uslab = uslab / (plon*(jslab2-jslab1+1))
      vslab = vslab / (plon*(jslab2-jslab1+1))


c Calculate zonal mean eastward and northward values uav,vav, and
c constrain them by umax. Set the final field, weighted with the
c unfiltered field by w, which varies linearly with latitude to 1 at NP.
c The filtered field is pure slab motion + uav (rot) + vav(conv).
c Note if latslab1 = 90, umax = uav = vav = 0.

      do 200 j=jslab1,jslab2

        uavn = 0.
        vavn = 0.
        do 202 i=1,plon
          uavn = uavn + u1(i,j)
          vavn = vavn + v1(i,j)
  202   continue
        uavn = uavn / plon
        vavn = vavn / plon

        umax = umaxpd(icode) * (alatu(plat+1)-alatu(j)) * (180./pi)
        uavn = min ( umax, max(-umax,uavn) )
        vavn = min ( umax, max(-umax,vavn) )

        if (jslab1.lt.jslab2) then
          w = (1.*(j-jslab1) + wslab1*(jslab2-j)) / (jslab2-jslab1)
        else
          w = 1.
        endif

        do 204 i=1,plon
          u1(i,j)= w * (-uslab*slonu(i)   + vslab*clonu(i) + uavn)
     &         + (1.-w)*u1(i,j) 
          v1(i,j) = w * (-uslab*clonu(i) - vslab*slonu(i) + vavn)
     &          + (1.-w)*v1(i,j)
  204   continue

  200 continue


      call wrapu(u1)
      call wrapu(v1)

      return
      end
