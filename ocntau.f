# 1 "ocntau.F"
c Removing the declaration of the COMMON block /GRADNS/ to a header file.
c SJP 2009/04/25
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Minor enhancements to the coupling of the new, high-resolution version of the
c ocean model to the atmosphere.
c SJP 2007/12/21
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Declaration of /FCORR/ moved to FCORR.f.
c SJP 2004/09/10
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c Modified to include the gridpoint at the tip of the Antarctic Peninsula that,
c at R21 resolution, is treated as land by the AGCM, but as ocean by the OGCM.
c The other 10 such gridpoints were already correctly treated by this
c subroutine, but this additional gridpoint was inexplicably omitted.
c SJP 2003/06/10
c
c $Log: ocntau.f,v $
c Revision 1.16  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.15  2001/02/28 04:36:40  rot032
c Further tidy ups from HBG
c
c Revision 1.14  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.13  2000/06/20 02:08:34  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.12  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.11  1998/12/10  00:55:51  ldr
c HBG changes to V5-1-21
c
c Revision 1.10  1997/12/19  02:03:17  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.9  1996/10/24  01:03:30  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.8  1996/03/21  03:19:15  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.7  1994/08/08  16:22:51  ldr
c Fix up RCS header at top of file.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  leads - if T and semice=T allow leads in sea-ice model
c                  lcouple - flag to run model in coupled mode, default F
c                  savefcor - if true, save stats for flux-corrections for
c                             coupled ocean model
c
c              from common/ficeug in this subroutine
c                  ficeu - ice concentration on u-grid
c
c              from common/icocstr in this subroutine
c                  octau, octav - ice ocean stress
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from argument
c                  ijdyn - ice dynamics time step factor
c
c     In/Out:  from common/atocstr in this subroutine
c             otaux,otauy - wind stresses, zonal and meridional (dynes/cm**2)
c
c              from common/a2o in A2O.f
c                  otx,oty - coupled AGCM and ICE stress onto ocean
c
c              from common/fcorr in this subroutine
c                  sfluxes - surface flxues for flux corrections
c
      subroutine ocntau(ijdyn)

c Filler routine for ocean stresses from atmos model
c Set all land stresses to zero.
c Then average points with stresses to get expansion
c of stresses over the adjacent land points.
c This is done so that the change to the ocean U/V-grid
c (which is different to the atmos model U/V-grid)
c gets sufficient surrounding values with realistic magnitude stresses.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'OPARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ijdyn

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
c.. Ocean forcing terms
      include 'A2O.f'
      include 'FCORR.f'
      include 'FEWFLAGS.f'
      include 'LSMI.f'

      real ficeu,fice
      common/ficeug/ficeu(0:plon+1,plat+1),fice(plon,plat)

c.. atmos-ocean stresses
      real otaux,otauy
      common/atocstr/otaux(ln2,lat),otauy(ln2,lat)

c.. ice-ocean stresses
      real octau,octav
      common/icocstr/octau(0:plon+1,plat+1),octav(0:plon+1,plat+1)

C Local work arrays and variables
      real staux(0:lon+1,0:lat2+1),stauy(0:lon+1,0:lat2+1)
      real rtaux(0:lon+1,0:lat2+1),rtauy(0:lon+1,0:lat2+1)
      integer mskl(0:lon+1,0:lat2+1),msklsum(lon,lat2)

      integer i, ia, iap1, isum, j, ja, jap1, lg, ma, mg, ns

      real taux_agcm
      real tauy_agcm
      real xfice
      real ustar3

# 143


C Local data, functions etc

# 179


c...  When using the high-resolution version of the ocean model, there is no
c...  mismatch between the positions of the coastlines on the AGCM and OGCM
c...  grids. As a result, no grid-matching operations need to be performed.



      logical first
      data first/.true./
      save first,msklsum

C Start code : ----------------------------------------------------------

      if(first)then

c.... Create a mask of 1 for non-land, else 0
      do 21 ns=1,2
      do 21 lg=1,lat
      j=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 31 mg=1,lon
      ma=mg+(ns-1)*lon
      i=mg
      mskl(i,j)=1
   31 if(imsl(ma,lg).eq.4)mskl(i,j)=0
      mskl(0,j)=mskl(lon,j)
   21 mskl(lon+1,j)=mskl(1,j)
      do 41 i=0,lon+1
      mskl(i,0)=mskl(i,1)
   41 mskl(i,lat2+1)=mskl(i,lat2)

# 225


      do 47 j=1,lat2
      do 47 i=1,lon
      msklsum(i,j)=0
      if(mskl(i,j).eq.0)then
      msklsum(i,j)=
     &  mskl(i+1,j  )+mskl(i+1,j+1)+mskl(i  ,j+1)+mskl(i-1,j+1)
     & +mskl(i-1,j  )+mskl(i-1,j-1)+mskl(i  ,j-1)+mskl(i+1,j-1)
      endif
   47 continue

        first=.false.
      end if

c.... Set all land point stresses to zero
      do 10 lg=1,lat
      do 10 mg=1,ln2
      if(imsl(mg,lg).eq.4)then
        otaux(mg,lg)=0.0
        otauy(mg,lg)=0.0
      endif
   10 continue

c.... Make global (i,j) arrays from (mg,lg,ns) arrays with
c.... wrap around and extend N & S.
      DO 20 ns=1,2
      do 20 lg=1,lat
      j=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 30 mg=1,lon
      ma=mg+(ns-1)*lon
      i=mg
      staux(i,j)=otaux(ma,lg)
   30 stauy(i,j)=otauy(ma,lg)
      staux(0,j)=staux(lon,j)
      stauy(0,j)=stauy(lon,j)
      staux(lon+1,j)=staux(1,j)
   20 stauy(lon+1,j)=stauy(1,j)
      do 40 i=0,lon+1
      staux(i,0)=staux(i,1)
      stauy(i,0)=stauy(i,1)
      staux(i,lat2+1)=staux(i,lat2)
   40 stauy(i,lat2+1)=stauy(i,lat2)

# 292


c.... Keep a copy of staux,stauy

      do 45 j=1,lat2
      do 45 i=1,lon
      rtaux(i,j)=staux(i,j)
   45 rtauy(i,j)=stauy(i,j)

c.... Run through the mask checking all points that have a
c.... zero in the mask (i.e. a land point).
c.... If zero, try to create a stress at that point by averaging
c.... the 8 surrounding points, using only values that are
c.... non-zero (non-land values)
      do 50 j=1,lat2
      do 50 i=1,lon
      isum=msklsum(i,j)
      if(isum.gt.0)then
      rtaux(i,j)=(     staux(i+1,j  )+staux(i+1,j+1)
     & +staux(i  ,j+1)+staux(i-1,j+1)+staux(i-1,j  )
     & +staux(i-1,j-1)+staux(i  ,j-1)+staux(i+1,j-1) )/isum
      rtauy(i,j)=(     stauy(i+1,j  )+stauy(i+1,j+1)
     & +stauy(i  ,j+1)+stauy(i-1,j+1)+stauy(i-1,j  )
     & +stauy(i-1,j-1)+stauy(i  ,j-1)+stauy(i+1,j-1) )/isum
      endif
   50 continue
     
c.... wrap around and extend

      do 60 j=1,lat2
      rtaux(0,j)=rtaux(lon,j)
      rtauy(0,j)=rtauy(lon,j)
      rtaux(lon+1,j)=rtaux(1,j)
   60 rtauy(lon+1,j)=rtauy(1,j)
      do 70 i=0,lon+1
      rtaux(i,0)=rtaux(i,1)
      rtauy(i,0)=rtauy(i,1)
      rtaux(i,lat2+1)=rtaux(i,lat2)
   70 rtauy(i,lat2+1)=rtauy(i,lat2)

c.... Now transfer the data to the ocean U/V grid (a NE shift)

      do 80 j=1,lat2
      do 80 i=1,lon
      staux(i,j)=0.25*( rtaux(i,j)+rtaux(i+1,j)
     &  +rtaux(i,j+1)+rtaux(i+1,j+1)  )
      stauy(i,j)=0.25*( rtauy(i,j)+rtauy(i+1,j)
     &  +rtauy(i,j+1)+rtauy(i+1,j+1)  )
   80 continue

c Before these ocean stresses can be used, the ice-ocean stresses
c from the dynamical ice model (if in use) must be added.
c They must be weighted by the leads area. There is a fractional
c ice cover given by ficeu on the ice-model U/V grid. Note that the
c ice-model U/V grid has a SW shift (unlike the ocean model NE shift)
c in the U/V index.

      If(leads)Then

      do 140 j=1,lat2
      do 140 i=1,lon
c..  ficeu is the fractional coverage on the U/V grid obtained by
c..  enlarging the ice cover over adjacent land areas, before
c..  interpolating. See icedrive.f which call icefhx.f to spread
c..  fice into ficex (extended over land points). Then dynice.f
c..  uses ltou/ltoh to interpolate fice onto the U/V grid (ficeu).
c..  This means that ficeu will exist for both ice points and surrounding
c..  pseudo land areas. However, the ocean model will only use stresses
c..  for the ocean points and will zero out stresses over land points.
      xfice=ficeu(i+1,j+1)
      if(xfice.gt.0.0)then
        staux(i,j)=(1.0-xfice)*staux(i,j)+xfice*octau(i+1,j+1)*10.0
        stauy(i,j)=(1.0-xfice)*stauy(i,j)+xfice*octav(i+1,j+1)*10.0
      endif
  140 continue

      End If

c.... Save these for forcing the ocean model.
c.... The ocean model assumes that forcing terms have been computed
c.... and summed for each atmos step.
c.... There is a factor ijdyn since the ice dynamics is only carried
c.... out every ijdyn atmos steps, and the atmos-ocean stresses, although
c.... computed every atmos step, can only be used to create the total
c.... stress when the ice dynamics stresses are available. ijdyn is usually
c.... 2 when the ice model is in use, else it will be 1.

c.... For calculating flux corrections :
      if(savefcor)then
      do 110 j=1,lat2
      do 110 i=1,lon
        sfluxes(i,j,1)=sfluxes(i,j,1)+staux(i,j)*ijdyn
        sfluxes(i,j,2)=sfluxes(i,j,2)+stauy(i,j)*ijdyn
c-- compute (u*)**3 in units (m/sec)**3
c-- note that stresses (staux and stauy) were set up in radstres.f
c--  with a factor of 10 (to get from Kgm/m/sec**2 to dyne/cm**2)
          taux_agcm=staux(i,j)/10.0      ! Kgm/m/sec**2
          tauy_agcm=stauy(i,j)/10.0      ! Kgm/m/sec**2
c-- ustar3 units (m/sec)**3
c-- 1025 is sea water density in Kgm/m**3
          ustar3=( sqrt(taux_agcm**2+tauy_agcm**2) /1025.0 )**1.5
        sfluxes(i,j,6)=sfluxes(i,j,6)+    ustar3*ijdyn
  110 continue
      endif

c.... For Ocean model forcing :
      if(lcouple)then

c.. MOM1

# 411


c...  Interpolate the wind stresses onto the new ocean model grid

      do j = 1, jmt-2
        ja = j  / 2
        jap1 = (j + 1) / 2

        do i = 1, imt-2
          ia = i / 2
          iap1 = (i + 1) / 2
          if (i .eq. 1) ia = (imt - 1)/2

          if (mod(i, 2) .eq. 0) then
            if (mod(j, 2) .eq. 0) then

c...  i even, j even
              otx(i, j) = ijdyn * staux(ia, ja)
              oty(i, j) = ijdyn * stauy(ia, ja)

            else

c...  i even, j odd
              if (j .eq. 1) then
                otx(i, j) = ijdyn * staux(ia, jap1)
                oty(i, j) = ijdyn * stauy(ia, jap1)
              else
                otx(i, j) = ijdyn *
     &            (cos(alatu(ja+1)) * staux(ia, ja) +
     &             cos(alatu(ja+2)) * staux(ia, jap1)) /
     &            (cos(alatu(ja+1)) + cos(alatu(ja+2)))
                oty(i, j) = ijdyn *
     &            (cos(alatu(ja+1)) * stauy(ia, ja) +
     &             cos(alatu(ja+2)) * stauy(ia, jap1)) /
     &            (cos(alatu(ja+1)) + cos(alatu(ja+2)))
              end if

            endif
          else
            if (mod(j, 2) .eq. 0) then

c...  i odd, j even
              otx(i, j) = 0.5 * ijdyn * (staux(ia, ja) +
     &                                   staux(iap1, ja))
              oty(i, j) = 0.5 * ijdyn * (stauy(ia, ja) +
     &                                   stauy(iap1, ja))

            else

c...  i odd, j odd
              if (j .eq. 1) then
                otx(i, j) = 0.5 * ijdyn * (staux(ia, jap1) +
     &                                     staux(iap1, jap1))
                oty(i, j) = 0.5 * ijdyn * (stauy(ia, jap1) +
     &                                     stauy(iap1, jap1))
              else
                otx(i, j) = ijdyn *
     &            (cos(alatu(ja+1)) * staux(ia, ja) +
     &             cos(alatu(ja+1)) * staux(iap1, ja) +
     &             cos(alatu(ja+2)) * staux(ia, jap1) +
     &             cos(alatu(ja+2)) * staux(iap1, jap1)) /
     &            (2.0 * (cos(alatu(ja+1)) + cos(alatu(ja+2))))
                oty(i, j) = ijdyn *
     &            (cos(alatu(ja+1)) * stauy(ia, ja) +
     &             cos(alatu(ja+1)) * stauy(iap1, ja) +
     &             cos(alatu(ja+2)) * stauy(ia, jap1) +
     &             cos(alatu(ja+2)) * stauy(iap1, jap1)) /
     &            (2.0 * (cos(alatu(ja+1)) + cos(alatu(ja+2))))
              end if

            end if
          end if
        end do
      end do



      endif

      return
      end
c---------------------------------------------------------------------
c
c Similar to above, but to suit MOM2 model
c
      subroutine ocntauMOM2(ijdyn)

c Filler routine for ocean stresses from atmos model
c
c  This is for the MOM2 code
c
c  The collection of stresses for data/flux corrections is prepared
c  on an equivalent AGCM UV grid.
c
c  The AGCM stress data to be passed directly to the MOM2 model
c  (when coupled) is prepared on the AGCM T grid. The MOM2 code
c  (see gosbc.F) will transfer this to the OGCM UV grid
c
c Set all land stresses to zero.
c Then average points with stresses to get expansion
c of stresses over the adjacent land points.
c This is done so that the change to the ocean U/V-grid
c (which is different to the atmos model U/V-grid)
c gets sufficient surrounding values with realistic magnitude stresses.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ijdyn

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
      include 'FCORR.f'
      include 'FEWFLAGS.f'
      include 'LSMI.f'
c      include 'GRADNS.f'

      real ficeu,fice
      common/ficeug/ficeu(0:plon+1,plat+1),fice(plon,plat)

c     real x90

c.. atmos-ocean stresses
      real otaux,otauy
      common/atocstr/otaux(ln2,lat),otauy(ln2,lat)

c.. ice-ocean stresses
      real octau,octav
      common/icocstr/octau(0:plon+1,plat+1),octav(0:plon+1,plat+1)

c.. MOM2 data
      real atmvar
      common/CSIROAGCM/atmvar(lon,lat2,6) ! ima=lon, jma=lat2  to match MOM2

C Local work arrays and variables
      real staux(0:lon+1,0:lat2+1),stauy(0:lon+1,0:lat2+1)
      real rtaux(0:lon+1,0:lat2+1),rtauy(0:lon+1,0:lat2+1)
      integer mskl(0:lon+1,0:lat2+1),msklsum(lon,lat2)
      real ttaux(lon,lat2),ttauy(lon,lat2)

      integer i
      integer isum
      integer j
      integer lg
      integer ma
      integer mg
      integer ns

      real avtau
      real avtav
      real taux_agcm
      real tauy_agcm
      real xfice
      real ustar3

C Local data, functions etc
      logical first
      data first/.true./
      save first,msklsum

C Start code : ----------------------------------------------------------

      if(first)then

c.... Create a mask of 1 for non-land, else 0
      do 21 ns=1,2
      do 21 lg=1,lat
      j=lg
      if(ns.eq.1)j=lat2+1-lg
      do 31 mg=1,lon
      ma=mg+(ns-1)*lon
      i=mg
      mskl(i,j)=1
   31 if(imsl(ma,lg).eq.4)mskl(i,j)=0
      mskl(0,j)=mskl(lon,j)
   21 mskl(lon+1,j)=mskl(1,j)
      do 41 i=0,lon+1
      mskl(i,0)=mskl(i,1)
   41 mskl(i,lat2+1)=mskl(i,lat2)

      do 47 j=1,lat2
      do 47 i=1,lon
      msklsum(i,j)=0
      if(mskl(i,j).eq.0)then
      msklsum(i,j)=
     &  mskl(i+1,j  )+mskl(i+1,j+1)+mskl(i  ,j+1)+mskl(i-1,j+1)
     & +mskl(i-1,j  )+mskl(i-1,j-1)+mskl(i  ,j-1)+mskl(i+1,j-1)
      endif
   47 continue

        first=.false.
      end if

c.... Set all land point stresses to zero
      do 10 lg=1,lat
      do 10 mg=1,ln2
      if(imsl(mg,lg).eq.4)then
        otaux(mg,lg)=0.0
        otauy(mg,lg)=0.0
      endif
   10 continue

c----------------------------------------------------------
c.... The next part involving staux,stauy is to compute
c....  stresses for ocean model on the ocean model U,V grid
c.... These may be used for forcing the ocean model
c....  or for computing flux corrections.
c.... (They will not be used for driving MOM2 because the
c....  AGCM data for that is assumed to be on the AGCM T
c....  grid, not the (NE shift) UV grid. See next section
c....  for that input to the MOM2 model.)
c----------------------------------------------------------

c.... Make global (i,j) arrays from (mg,lg,ns) arrays with
c.... wrap around and extend N & S.
      DO 20 ns=1,2
      do 20 lg=1,lat
      j=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 30 mg=1,lon
      ma=mg+(ns-1)*lon
      i=mg
      staux(i,j)=otaux(ma,lg)
   30 stauy(i,j)=otauy(ma,lg)
      staux(0,j)=staux(lon,j)
      stauy(0,j)=stauy(lon,j)
      staux(lon+1,j)=staux(1,j)
   20 stauy(lon+1,j)=stauy(1,j)
      do 40 i=0,lon+1
      staux(i,0)=staux(i,1)
      stauy(i,0)=stauy(i,1)
      staux(i,lat2+1)=staux(i,lat2)
   40 stauy(i,lat2+1)=stauy(i,lat2)

c.... Keep a copy of staux,stauy

      do 45 j=1,lat2
      do 45 i=1,lon
      rtaux(i,j)=staux(i,j)
   45 rtauy(i,j)=stauy(i,j)

c.... Run through the mask checking all points that have a
c.... zero in the mask (i.e. a land point).
c.... If zero, try to create a stress at that point by averaging
c.... the 8 surrounding points, using only values that are
c.... non-zero (non-land values)
      do 50 j=1,lat2
      do 50 i=1,lon
      isum=msklsum(i,j)
      if(isum.gt.0)then
      rtaux(i,j)=(     staux(i+1,j  )+staux(i+1,j+1)
     & +staux(i  ,j+1)+staux(i-1,j+1)+staux(i-1,j  )
     & +staux(i-1,j-1)+staux(i  ,j-1)+staux(i+1,j-1) )/isum
      rtauy(i,j)=(     stauy(i+1,j  )+stauy(i+1,j+1)
     & +stauy(i  ,j+1)+stauy(i-1,j+1)+stauy(i-1,j  )
     & +stauy(i-1,j-1)+stauy(i  ,j-1)+stauy(i+1,j-1) )/isum
      endif
   50 continue
     
c.... wrap around and extend

      do 60 j=1,lat2
      rtaux(0,j)=rtaux(lon,j)
      rtauy(0,j)=rtauy(lon,j)
      rtaux(lon+1,j)=rtaux(1,j)
   60 rtauy(lon+1,j)=rtauy(1,j)
      do 70 i=0,lon+1
      rtaux(i,0)=rtaux(i,1)
      rtauy(i,0)=rtauy(i,1)
      rtaux(i,lat2+1)=rtaux(i,lat2)
   70 rtauy(i,lat2+1)=rtauy(i,lat2)

c.... Now transfer the data to the ocean U/V grid (a NE shift)

      do 80 j=1,lat2
      do 80 i=1,lon
      staux(i,j)=0.25*( rtaux(i,j)+rtaux(i+1,j)
     &  +rtaux(i,j+1)+rtaux(i+1,j+1)  )
      stauy(i,j)=0.25*( rtauy(i,j)+rtauy(i+1,j)
     &  +rtauy(i,j+1)+rtauy(i+1,j+1)  )
   80 continue

c Before these ocean stresses can be used, the ice-ocean stresses
c from the dynamical ice model (if in use) must be added.
c They must be weighted by the leads area. There is a fractional
c ice cover given by ficeu on the ice-model U/V grid. Note that the
c ice-model U/V grid has a SW shift (unlike the ocean model NE shift)
c in the U/V index.

      If(leads)Then

      do 140 j=1,lat2
      do 140 i=1,lon
c..  ficeu is the fractional coverage on the U/V grid obtained by
c..  enlarging the ice cover over adjacent land areas, before
c..  interpolating. See icedrive.f which call icefhx.f to spread
c..  fice into ficex (extended over land points). Then dynice.f
c..  uses ltou/ltoh to interpolate fice onto the U/V grid (ficeu).
c..  This means that ficeu will exist for both ice points and surrounding
c..  pseudo land areas. However, the ocean model will only use stresses
c..  for the ocean points and will zero out stresses over land points.
      xfice=ficeu(i+1,j+1)
      if(xfice.gt.0.0)then
        staux(i,j)=(1.0-xfice)*staux(i,j)+xfice*octau(i+1,j+1)*10.0
        stauy(i,j)=(1.0-xfice)*stauy(i,j)+xfice*octav(i+1,j+1)*10.0
      endif
  140 continue

      End If

c.... Save these for forcing the ocean model.
c.... The ocean model assumes that forcing terms have been computed
c.... and summed for each atmos step.
c.... There is a factor ijdyn since the ice dynamics is only carried
c.... out every ijdyn atmos steps, and the atmos-ocean stresses, although
c.... computed every atmos step, can only be used to create the total
c.... stress when the ice dynamics stresses are available. ijdyn is usually
c.... 2 when the ice model is in use, else it will be 1.

c.... For data or calculating flux corrections or driving MOM2 in stand
c....  alone mode.(Note that the data "X" is on a NE shifted grid
c....  relative to the AGCM A points)
c
c        A...A
c        .....
c        ..X..
c        .....
c        A...A
c
c     X longitudes are thus ((mg-0.5)*360/lon,mg=1,lon)
c     X latitudes are (0.5*(glats(j)+glats(j+1)),j=1,lat2)

c     print *,lat2
c     x90=90.0
c     write(6,1616)(0.5*(glats(j)+glats(j+1))*180/pi,j=1,lat2-1),x90
c     print *,lon
c     write(6,1616)((mg-0.5)*360/lon,mg=1,lon)
c1616 format(8(1x,f8.4))

      if(savefcor)then
      do 110 j=1,lat2
      do 110 i=1,lon
        sfluxes(i,j,1)=sfluxes(i,j,1)+staux(i,j)*ijdyn
        sfluxes(i,j,2)=sfluxes(i,j,2)+stauy(i,j)*ijdyn
c-- compute (u*)**3 in units (m/sec)**3
c-- note that stresses (staux and stauy) were set up in radstres.f
c--  with a factor of 10 (to get from Kgm/m/sec**2 to dyne/cm**2)
          taux_agcm=staux(i,j)/10.0        ! Kgm/m/sec**2
          tauy_agcm=stauy(i,j)/10.0        ! Kgm/m/sec**2
c-- ustar3 units (m/sec)**3
c-- 1025 is sea water density in Kgm/m**3
          ustar3=( sqrt(taux_agcm**2+tauy_agcm**2) /1025.0 )**1.5
        sfluxes(i,j,6)=sfluxes(i,j,6)+    ustar3*ijdyn
  110 continue
      endif

c----------------------------------------------------------
c.... The next part involving staux,stauy is to compute
c....  stresses on the AGCM T grion the AGCM T grid.
c.... These will be used for driving MOM2 because the
c.... MOM2 code (see gosbc.F) assunes that these stresses
c.... are on the AGCM T grid, not the (NE shift) UV grid.
c----------------------------------------------------------

c.... Make global (i,j) arrays from (mg,lg,ns) arrays
      do ns=1,2
      do lg=1,lat
        j=(ns-1)*lg+(lat2+1-lg)*(2-ns)
        do mg=1,lon
          ma=mg+(ns-1)*lon
          i=mg
          ttaux(i,j)=otaux(ma,lg)
          ttauy(i,j)=otauy(ma,lg)
        enddo
      enddo
      enddo

c Before these ocean stresses can be used, the ice-ocean stresses
c from the dynamical ice model (if in use) must be added.
c They must be weighted by the leads area. There is a fractional
c ice cover given by fice on the ice-model T grid. Note that the
c ice-model U/V grid has a SW shift (unlike the ocean model NE shift)
c in the U/V index.

      If(leads)Then

      do j=1,lat2
      do i=1,lon
c..  fice is the fractional coverage on the AGCM T grid
c..  The ice stresses are on the ice-model U/V points
c
c      ITau(i,j+1) ------------ ITau(i+1,j+1)
c          |                        |
c          |                        |
c          |                        |
c          |         Ice(i,j)       |
c          |                        |
c          |                        |
c          |                        |
c      ITau(i,j) -------------- ITau(i+1,j)
c
c
      xfice=fice(i,j)
      if(xfice.gt.0.0)then
        avtau=0.25*(octau(i,j)+octau(i+1,j)+octau(i,j+1)+octau(i+1,j+1))
        avtav=0.25*(octav(i,j)+octav(i+1,j)+octav(i,j+1)+octav(i+1,j+1))
        ttaux(i,j)=(1.0-xfice)*ttaux(i,j)+xfice*avtau*10.0
        ttauy(i,j)=(1.0-xfice)*ttauy(i,j)+xfice*avtav*10.0
      endif
      enddo
      enddo

      End If

c.... Save these for forcing the ocean model.
c.... The ocean model assumes that forcing terms have been computed
c.... and summed for each atmos step.
c.... There is a factor ijdyn since the ice dynamics is only carried
c.... out every ijdyn atmos steps, and the atmos-ocean stresses, although
c.... computed every atmos step, can only be used to create the total
c.... stress when the ice dynamics stresses are available. ijdyn is usually
c.... 2 when the ice model is in use, else it will be 1.

c.... For Ocean model forcing :
      if(lcouple)then

c.. MOM2 data (averaging done by MOM2 code)
      do j=1,lat2
      do i=1,lon
         atmvar(i,j,1)=ttaux(i,j)*ijdyn ! Stress in dyne/cm**2
         atmvar(i,j,2)=ttauy(i,j)*ijdyn ! Stress in dyne/cm**2
c-- compute (u*)**3 in units (m/sec)**3
c-- note that stresses (ttaux and ttauy) were set up in radstres.f
c--  with a factor of 10 (to get from Kgm/m/sec**2 to dyne/cm**2)
          taux_agcm=ttaux(i,j)/10.0      ! Kgm/m/sec**2
          tauy_agcm=ttauy(i,j)/10.0      ! Kgm/m/sec**2
c-- ustar3 units (m/sec)**3
c-- 1025 is sea water density in Kgm/m**3
          ustar3=( sqrt(taux_agcm**2+tauy_agcm**2) /1025.0 )**1.5
         atmvar(i,j,6)=    ustar3*ijdyn
      enddo
      enddo

      endif

      return
      end
