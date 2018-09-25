c Tidying up the source code for the sea ice model, particularly: fixing the
c bug in the calculation of the latitudes; and removing the buggy and 
c experimental code which sought to impose "limited slip at low resolution".
c SJP 2009/04/25
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: dynice.f,v $
c Revision 1.16  1998/12/10 00:55:34  ldr
c HBG changes to V5-1-21
c
c Revision 1.15  1997/12/19  02:03:11  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.14  1996/03/21  03:18:35  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.13  1994/08/08  17:11:05  ldr
c Move wrapu call from icefree to dynice.
c
c Revision 1.12  94/08/04  16:55:00  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.11  94/03/30  12:58:37  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.10  94/03/22  15:32:36  ldr
c Removed spurious common block lsmi.
c 
c Revision 1.9.1.1  94/03/30  12:34:13  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.9  93/12/17  15:32:08  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.8  93/11/03  11:44:12  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.7  93/08/10  15:26:14  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.6  93/07/27  17:12:30  ldr
c Changes to ensure coherent approach to pure restarts, monthly
c interpolation factors.
c 
      subroutine dynice(fice, hice, ticel, hsno, tsno, fuwindi, fvwindi,
     &                  uocn, vocn, ficeu, ficex, hicex, stor)

c Main driver of the dynamical sea-ice model of Flato and Hibler (1990),
c Annals of Glaciology,14,72-77. Two cavitating-fluid iteration methods
c are provided, selected by logical flag flato in comdice. In one
c method (flato=.true.), the free-drift equations are solved just once,
c then the velocities themselves are iteratively adjusted to remove or 
c reduce convergence, incrementing a "pressure" field up to a maximum
c strength. In the other method (flato=.false.), the "free-drift +
c pressure" equations are solved iteratively, again incrementing the
c pressure field where there is convergence up to a maximum strength.

c A staggered Arakawa-B grid is used, with thickness and concentration
c on the "h-grid", and velocities on the "u-grid" at the corners of the
c h-grid. The u-grid pt (i,j) is at the SW corner of h-grid pt (i,j).
c Land-ocean grid boxes from lsx coincide with h-grid boxes, so that
c h-grid edges are at the poles. (So u-grid centers are at the poles;
c u and v are set to zero at the poles.) A no-slip boundary condition
c is applied at land boundaries.

c All velocities are zeroed over land by a local "umask", so 
c field values on land are never changed by this package. This is
c necessary since lsx uses the same thickness and temperature arrays 
c for soil and seaice, and for vectorization this package computes
c results both over land and ocean. (Also see comments below loop 500
c in advect. u-grid arrays (including those
c supplied and returned to the ocean model) are wrapped in longitude
c but not in latitude, so are dimensioned (0:plon+1, plat+1), with
c lat # 1 at the south pole and lat # plat+1 at the north pole.

c Ice velocities at the poles are fixed at zero, so some u-grid fields
c including wind stresses and ocean currents at j=1 and plat+1 are
c never used. However for uniformity all u-grid arrays have the same 
c dimensions.

c The include-file comdice contains variables unique to this model.
c The files compar and comgrd contain basic commons shared with the
c lsx surface package of the agcm. MKS units are used throughout.

c All arguments are supplied and on the h-grid unless otherwise noted.

c *** The "flato=.false." iterative method is relatively untested ***
c *** Ocean currents uocn,vocn are currently modified by polefilt ***

c fice    = seaice areal fraction, ie, concentration (modified)
c hicelay = thickness of each seaice layer, meters (modified)
c ticelay = temperature of each seaice layer, deg K (modified)
c hsnolay = thickness of each snow layer, meters (modified)
c tsnolay = temperature of each snow layer, deg K (modified)
c fuwindi = eastward wind stress on seaice, Nl/m**2 (modified)
c fvwindi = northward wind stress on seaice, Nl/m**2 (modified)
c uocn    = eastward  surface ocean current, m/s (u-grid)
c vocn    = northward surface ocean current, m/s (u-grid)

      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'COMDICE.f' !

      dimension
     &  fice(plon,plat),             stor(plon,plat),
     &  hice(plon,plat),             ticel(plon,plat,2),
     &  hsno(plon,plat),             tsno(plon,plat),
     &  fuwindi(plon,plat),          fvwindi(plon,plat),
     &  uocn(0:plon+1,plat+1),       vocn(0:plon+1,plat+1),
     &  ficex(plon,plat),hicex(plon,plat)

c Local h-grid arrays

      dimension 
     &  ficeold(plon,plat),          hiceold(plon,plat),
     &  hsnoold(plon,plat),          div(plon,plat),
     &  pres(plon,plat),             streng(plon,plat)

c Local u-grid arrays

      dimension
     &  hiceu(0:plon+1,plat+1),      ficeu(0:plon+1,plat+1),
     &  fuwindis(0:plon+1,plat+1),   fvwindis(0:plon+1,plat+1)

      save fuwindis, fvwindis

      common/icocstr/octau(0:plon+1,plat+1),octav(0:plon+1,plat+1)

      logical ifdone

c Statement function for ice strength (Nl m-1)
      fstreng(h,f) = 2.5e4 * h * h * exp (-20.*(1.-f))

c Initial setup, done once only at the (re)start of each run

c Initialize some local arrays globally, since ice model will operate on
c reduced latitudinal regions of non-zero ice, but sometimes refer to
c adjacent points. (Also hice, ficeold will be needed globally for 
c history field computations.)

      do 50 j=1,plat
      do 50 i=1,plon
        adgain(i,j)=0.0
        ficeold(i,j) = fice(i,j)
        icolaps(i,j) = -1
        workp(i,j)=0.0
   50 continue

      do 52 j=1,plat+1
      do 52 i=0,plon+1
        u(i,j) = 0.
        v(i,j) = 0.
        hiceu(i,j) = 0.
        ficeu(i,j) = 0.
        octau(i,j) = 0.
        octav(i,j) = 0.
   52 continue

c Set latha,etc, to global extent for calls to timefilt,polefilt below

      latha = 1
      lathb = plat
      latua = 1
      latub = plat+1

c Time-smooth and polar-filter the seaice wind stresses

      call timefilt (fuwindi, fvwindi, fuwindis, fvwindis)
      call polefilt (fuwindis, fvwindis, 83, 1)

c Overall loop for each hemisphere (1 for SH, 2 for NH)

      do 1000 ihem = 1,2

c Find the reduced latitudinal band for this hemisphere, bounded by
c latha,lathb (h-grid) and latua,latub (u-grid). Skip this hemisphere
c if no ice, flagged by latha=0. (latha,etc, are in common dicegrid.)

      call icebound (fice, ihem)
      if (latha.eq.0) goto 1000

c Interpolate hice and fice to u-grid (hiceu, ficeu)
chal use extended ice values over land points
      call ltou (hicex, hiceu)
      call ltou (ficex, ficeu)

c Dynamical section for u and v.
c First, initialize pressure to zero and set strength.

      do 200 j=latha,lathb
      do 200 i=1,plon
        pres(i,j) = 0.
        streng(i,j) = fstreng ( hice(i,j), fice(i,j) )
  200 continue

c If "flato" method, solve for free-drift velocities. If "not flato"
c method, initialize u,v to zero for the first iteration.

      if (flato) then
        call icefree (hiceu,fuwindis,fvwindis,uocn,vocn,pres)
      else
        do 210 j=latua,latub
        do 210 i=0,plon+1
          u(i,j) = 0.
          v(i,j) = 0.
  210   continue
      endif

c Begin cavitating-fluid iteration

      do 250 iloop=1,nloop

c          Save pre-iteration u,v for convergence test in icediag.
c          (So cavit[2] and icefree must not change worku[2]).

        do 252 j=latua,latub
        do 252 i=0,plon+1
          worku(i,j) = u(i,j)
          worku2(i,j) = v(i,j)
  252   continue
 
c          If "flato" method, apply cavitating-fluid increment to u,v
c          in cavit. If "not flato" method, apply cavitating-fluid
c          increment to pressure pres in cavit2, then solve for the
c          "free-drift + pressure" u,v in icefree.

        if (flato) then
          call cavit (div, pres, streng)
        else
          call cavit2 (div, pres, streng)
          call icefree (hiceu,fuwindis,fvwindis,uocn,vocn,pres)
        endif

      if(iloop.eq.2)then
        do j=1,plat
        do i=1,plon
          work2nd(i,j)=workp(i,j)
        enddo
        enddo
      endif
      call polefilt(u,v,80,4)

c          Test if converged - if so, compute diagnostics and skip out

        call icediag (fice,ficeu, pres,streng, iloop, ifdone, ihem)

chbg Prevent two step de-coupling of iterating solution
        if(flato)then
          do j=latua,latub
          do i=0,plon+1
            u(i,j) = 0.7*u(i,j)+0.3*worku(i,j)
            v(i,j) = 0.7*v(i,j)+0.3*worku2(i,j)
          enddo
          enddo
        endif
        if (ifdone) goto 260

  250 continue
  260 continue

c set up divergence arrays.
      do j=1,plat
      do i=1,plon
        workl(i,j)=workp(i,j)
        workdf(i,j)=workl(i,j)-work2nd(i,j)
      enddo
      enddo

c Compute the net stress on the ocean (f[u,v]ocn, on u-grid) for use
c by the dynamic ocean model

      do 270  j=latua,latub
      do 270  i=1,plon
      octau(i,j)=csubw*(coswturn(j)*(u(i,j)-uocn(i,j))-sinwturn(j)*
     &(v(i,j)-vocn(i,j)))
      octav(i,j)=csubw*(coswturn(j)*(v(i,j)-vocn(i,j))+sinwturn(j)*
     &(u(i,j)-uocn(i,j)))
 270  continue 

c Advection section.
c Fields are advected in a particular order so already-advected fields
c can be used to recover thickness from mass and temperature from heat.
c Each layer of the multi-layer thermodynamic ice and snow models is
c advected totally independently of the other layers (!). Any illogical
c results (eg, snow with no ice) should be taken care of by logical
c adjustments in the thermodynamic model(s), but I don't think such
c cases can arise here except possibly in pathological situations(?).

c        advect ice concentration

      do 400 j=latha,lathb
      do 400 i=1,plon
        ficeold(i,j) = fice(i,j)
        work(i,j) = fice(i,j)
        work2(i,j) = 1.
  400 continue
      call advect (fice, work, work2, 0,ihem,pres,streng)

c          advect ice layer mass, recovering thickness

        do 422 j=latha,lathb
        do 422 i=1,plon
          hiceold(i,j) = hice(i,j)
          work(i,j) = ficeold(i,j) * hice(i,j)
          work2(i,j) = fice(i,j)
  422   continue
        call advect (hice(1,1), work, work2, 1,ihem,pres,streng)

c adjust for concentrations greater than 1 
       do 460 j=latha,lathb
       do 460 i=1,plon
       if(fice(i,j).ge.1.00)then
       hice(i,j)=hice(i,j)*fice(i,j)
        fice(i,j)=1.00
       endif
  460  continue

c sum advective gain in ice thickness
c multiply by the change in ice concentration to avoid bias to ice edge points
        do 432 j=latha,lathb
        do 432 i=1,plon
        if(imask(i,j).eq.1.and.fice(i,j).gt.1.e-3)then
        adgain(i,j)=(hice(i,j)-hiceold(i,j))
        endif
  432   continue

c          advect ice layer heat, recovering temperature

        do 425 k=1,2
        do 424 j=latha,lathb
        do 424 i=1,plon
          work(i,j)=0.5*ficeold(i,j)*hiceold(i,j)*ticel(i,j,k)
          work2(i,j) = 0.5*fice(i,j) * hice(i,j)
  424   continue
        call advect (ticel(1,1,k), work, work2, 2,ihem,pres,streng)

  425   continue

c        advect ice brine resevoir

      do 440 j=latha,lathb
      do 440 i=1,plon
        work(i,j) = ficeold(i,j)*stor(i,j)
        work2(i,j) = fice(i,j)
  440 continue
      call advect (stor, work, work2, 3,ihem,pres,streng)

c          advect snow layer mass, recovering thickness

        do 522 j=latha,lathb
        do 522 i=1,plon
          hsnoold(i,j) = hsno(i,j)
          work(i,j) = ficeold(i,j) *  hsno(i,j)
          work2(i,j) = fice(i,j) 
  522   continue
        call advect (hsno(1,1), work, work2, 11,ihem,pres,streng)

c          advect snow layer heat, recovering temperature

        do 524 j=latha,lathb
        do 524 i=1,plon
          work(i,j)  = ficeold(i,j) *  hsnoold(i,j)
     &                 * tsno(i,j)
          work2(i,j) = fice(i,j)  * hsno(i,j)
  524   continue
        call advect (tsno(1,1), work, work2, 12,ihem,pres,streng)

c Kill off tiny amounts of ice or snow. ficemin and fsnomin (supplied
c by lsx) should be larger than what can typically be advected into a
c previously ice/snow-free grid box in one time step, otherwise the ice
c margin can never advance.

      do 600 j = latha,lathb
      do 600 i = 1,plon

         if(hsno(i,j).lt.0.0001)then
         hsno(i,j)=0.0
         tsno(i,j)=0.0
         endif
    
  600 continue

c Loop back for other hemisphere

      call wrapu (octau)
      call wrapu (octav)

 1000 continue

      return
      end
