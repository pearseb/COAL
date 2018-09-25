c $Log: icediag.f,v $
c Revision 1.4  1996/10/24 01:02:54  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1996/03/21  03:18:49  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.2  1993/12/17  15:32:49  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  92/06/16  11:54:51  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/ocwind in COMDICE.f
c                  u - eastward  velocity on u-grid 
c                  v - northward velocity on u-grid 
c
c              from common/dicework in COMDICE.f
c                  worku,worku2 - velocity work arrays
c
c              from common/dicegrid in COMDICE.f
c                  latha, lathb - h-grid)latitudinal indices bounding
c                  latua, latub - u-grid)the extent of ice
c                  srfarea - surface area of grid boxes
c
c              from arguments
c                  fice  - ice concentration on h-grid  
c                  ficeu - ice concentration on u-grid 
c                  ihem  - 1 for SH, 2 for NH
c                  iloop - iteration number 
c                  pres  - internal ice "pressure" on h-grid 
c                  streng- internal ice strength on h-grid 
c
c     Output:  from common/dicediag in COMDICE.f
c                  icolaps - flag for system collapse
c                  loop    - loop counter
c
c              from arguments
c                  ifdone - flag, set to .true. if last iteration
c
c     In/Out:  from common/dicediag in COMDICE.f
c                  rmsinc - root mean square u,v increment
c                  rmsvel - root mean square velocities
c
c 
      subroutine icediag (fice,ficeu,pres,streng,iloop,ifdone,ihem)

c (i)   Calculates rmsinc, the rms u,v increment for current iteration 
c (ii)  Tests if last iteration
c (iii) If last iteration, calculates diagnostics (in common comdiag)

c u      = eastward  velocity on u-grid (supplied)
c v      = northward velocity on u-grid (supplied)
c fice   = ice concentration on h-grid (supplied)
c ficeu  = ice concentration on u-grid (supplied)
c pres   = internal ice "pressure" on h-grid (supplied)
c streng = internal ice strength on h-grid (supplied)
c iloop  = iteration number (supplied)
c ifdone = flag, set to .true. if last iteration (returned)
c ihem   = 1 for SH, 2 for NH

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension 
     &   fice(plon,plat),              ficeu(0:plon+1,plat+1),
     &   pres(plon,plat),              streng(plon,plat)

      logical ifdone


      ifdone = .false.


c Calculate rms ice-area-weighted velocity increment for this iteration

      rmsinc(ihem) = 0.
      tareau = 0.
      do 100 j=latua,latub
        do 102 i=1,plon
          rmsinc(ihem) = rmsinc(ihem)
     &           + ((u(i,j)-worku(i,j))**2 + (v(i,j)-worku2(i,j))**2)
     &             * srfareau(j)*ficeu(i,j)
          tareau = tareau + srfareau(j)*ficeu(i,j)
  102   continue
  100 continue
      if (tareau.gt.0.) rmsinc(ihem) = sqrt(rmsinc(ihem)/tareau)


c     If last loop, set ifdone flag and calc diagnostics (in comdiag)

      if ( (rmsinc(ihem).le.convinc .and. iloop.ge.nloopmin)
     &     .or. iloop.eq.nloop ) then

        ifdone = .true.

        loop(ihem) = iloop

        rmsvel(ihem) = 0.
        do 200 j=latua,latub
        do 200 i=1,plon
          rmsvel(ihem) = rmsvel(ihem)
     &           + (u(i,j)**2 + v(i,j)**2) * srfareau(j)*ficeu(i,j)
  200   continue
        if (tareau.gt.0.) rmsvel(ihem) = sqrt(rmsvel(ihem)/tareau)
c       write(6,*)'rms  ',rmsvel(ihem)

        do 210 j=latha,lathb
        do 210 i=1,plon
          icolaps(i,j) = -1                                    !masscomp
          if (fice(i,j).gt.0.) then                            !masscomp
            icolaps(i,j) = 0                                   !masscomp
            if (pres(i,j).gt.0.) icolaps(i,j) = 1              !masscomp
            if (pres(i,j).ge.streng(i,j)) icolaps(i,j) = 2     !masscomp
          endif                                                !masscomp
cray      icolaps(i,j) = cvmgt (0, -1,
cray *                          fice(i,j).gt.0.) 
cray      icolaps(i,j) = cvmgt (1, icolaps(i,j),
cray *                          fice(i,j).gt.0. .and. pres(i,j).gt.0.)
cray      icolaps(i,j) = cvmgt (2, icolaps(i,j),
cray *                   fice(i,j).gt.0..and.pres(i,j).ge.streng(i,j))
  210   continue


c       Warn only if failed to get to 2*convinc...below that but 
c       above convinc is not too bad

c       if (rmsinc(ihem) .gt. 2.*convinc)
c    *  write(ioterm,220) ihem, iloop, nloop, rmsinc(ihem), convinc
c 220   format(' *** Warning. icediag 220. Bad cavit ice convergence.'
c    *        /' ihem,iloop,nloop,rmsinc,convinc=',3i4,2e15.5)

      endif

      return
      end
