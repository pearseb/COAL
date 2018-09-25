c Tidying up the source code for the sea ice model, particularly: fixing the
c bug in the calculation of the latitudes; and removing the buggy and 
c experimental code which sought to impose "limited slip at low resolution".
c SJP 2009/04/25
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: icefree.f,v $
c Revision 1.14  1998/12/10 00:55:48  ldr
c HBG changes to V5-1-21
c
c Revision 1.13  1997/12/19  02:03:16  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.12  1996/10/24  01:02:55  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.11  1996/06/13  02:06:52  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.10  1996/03/21  03:18:51  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.9  1995/08/31  04:30:44  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.8  1995/07/26  07:28:38  ldr
c Merge Mickles speedups to ice scheme (V4-5-30mic) into V4-7-3l.
c
c Revision 1.7  1994/08/08  17:11:33  ldr
c Move wrapu call from icefree to dynice.
c
c Revision 1.6  94/08/04  16:55:45  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.5.1.1  1995/07/26  07:24:08  ldr
c Speedups from Mickles replacing character mask with integers.
c
c Revision 1.5  94/03/30  12:34:33  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.4  93/12/17  15:32:53  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.3  93/11/03  11:44:16  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.2  93/02/03  12:44:48  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.1  92/06/16  11:54:55  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/dicegrid in COMDICE.f
c                  deltx - east-west widths of h-grid boxes
c                  deltyu - nl-south widths of u-grid boxes
c                  latua,latub - latitudinal indices bounding the extent
c                                of ice (u-grid)
c                  slatu - sine latitude scale
c
c              from common/dicemask in COMDICE.f
c                  umask - u-grid land-ocean mask 
c
c              from common/diceturn in COMDICE.f
c                  coswturn - water stress turning angle
c                  sinwturn - water stress turning angle
c
c              from arguments
c                  fuwindis - eastward wind stress on u-grid 
c                  fvwindis -  northward wind stress on u-grid ned)
c                  hiceu - total ice thickness on u-grid 
c                  uocn - eastward  surface ocean velocity on u-grid 
c                  vocn - northward surface ocean velocity on u-grid  
c
c     In/Out:  from common/dicework in COMDICE.f
c                  workh - h grid array
c
c              from common/ocwind in COMDICE.f
c                  u - eastward  velocity on u-grid 
c                  v - northward velocity on u-grid 
c 
c
      subroutine icefree(hiceu, fuwindis, fvwindis, uocn, vocn, pres)

c Solves for free-drift velocities u and v. 
c Note that total stress between agcm and ice is not conserved due
c to earlier transfer from h-grid f[u,v]windi to u-grid f[u,v]windis.
c (Could fix by weighting with fice before transfer and divide by ficeu
c after, but that would complicate polar filtering and time-smoothing.)

c u         = eastward  velocity on u-grid (returned)
c v        = northward velocity on u-grid (returned)
c hiceu    = total ice thickness on u-grid (supplied)
c fuwindis = eastward  wind stress on u-grid (supplied)
c fvwindis = northward wind stress on u-grid (supplied)
c uocn     = eastward  surface ocean velocity on u-grid (supplied) 
c vocn     = northward surface ocean velocity on u-grid (supplied) 
c pres     = internal ice "pressure" on h-grid (supplied)

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension
     &  hiceu(0:plon+1,plat+1),       fuwindis(0:plon+1,plat+1),
     &  fvwindis(0:plon+1,plat+1),    uocn(0:plon+1,plat+1),
     &  vocn(0:plon+1,plat+1),        pres(plon,plat)

c Transfer internal ice pressure to a wrapped h-grid array, since its 
c gradients on the u-grid are computed below

      call ltoh (pres, workh)

c Solve free-drift equations

      do 100 j=latua,latub
        a = csubw*coswturn(j)

        do 102 i=1,plon

c            Compute gradients of internal ice pressure

          prese = 0.5* ( (workh(i,j)  -workh(i-1,j)  )/deltx(j)
     &                  +(workh(i,j-1)-workh(i-1,j-1))/deltx(j-1) )
          presn = 0.5* ( (workh(i,j)  -workh(i,j-1)  )/deltyu(j)
     &                  +(workh(i-1,j)-workh(i-1,j-1))/deltyu(j) )

c            Compute gradients of dynamic topography assuming geostrophy

          slope =  2.*omega*slatu(j)*vocn(i,j) / grav
          slopn = -2.*omega*slatu(j)*uocn(i,j) / grav

c            add up all forcing terms (those indep of u and v) 

          uforce = fuwindis(i,j)*cosaturn(j)-sinaturn(j)*fvwindis(i,j)
     &           - prese
     &           - grav*rhoi*hiceu(i,j)*slope
     &           + csubw*coswturn(j)*uocn(i,j)
     &           - csubw*sinwturn(j)*vocn(i,j)

          vforce = fvwindis(i,j)*cosaturn(j)+sinaturn(j)*fuwindis(i,j)
     &           - presn
     &           - grav*rhoi*hiceu(i,j)*slopn
     &           + csubw*coswturn(j)*vocn(i,j)
     &           + csubw*sinwturn(j)*uocn(i,j)

c            solve the 2 linearized equations for u and v

          b = rhoi*hiceu(i,j)*2.*omega*slatu(j) + csubw*sinwturn(j)
          u(i,j) = umask(i,j) * (vforce*b + uforce*a) / (a*a+b*b)
          v(i,j) = umask(i,j) * (vforce*a - uforce*b) / (a*a+b*b)

  102   continue
  100 continue

c Set wraparound u,v points

      call wrapu (u)
      call wrapu (v)

      return
      end
