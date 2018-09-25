c Integrating Jess Trevena's freshwater hosing source code into the model.
c SJP 2009/08/06
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Fix a bug whereby DOWNSCALE was being called even when the model was running
c in stand-alone atmosphere mode.
c SJP 2008/02/03
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Replaced the line "include 'OCEANPARS.f'" with "include 'OPARAMS.f'",
c enabling the header file OCEANPARS.f to be removed from the model source
c code.
c SJP 2007/05/31
c
c Modified to divide DSIS, DSISB and DSFW by DT when saving for diagnostic
c purposes.
c SJP 2004/09/24
c
c Modified to calculate the salinity tendency from actual ocean salinities,
c applying a uniform global correction in order to ensure conservation of
c freshwater. For diagnostic purposes, a number of additional fields are now
c saved when SAVEFCOR=T.
c SJP 2004/09/13
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c Modified for the changes to the freshwater fluxes, whereby the ice water flux
c is split into its two components.
c SJP 2003/06/19
c
c Modifications to treatment of runoff at R21 resolution, as follows:
c (1) runoff is estimated for the 18 gridpoints that are treated as ocean by
c     the AGCM, but as land by the OGCM.
c (2) the runoff data is then passed to subroutine CONSERVE.
c SJP 2003/06/10
c
c Modified so that, at R21 resolution, the fields on the OGCM grid are passed
c to subroutine CONSERVE. This routine applies a correction to the values, in
c order to ensure that the global integral is conserved. Also, in loop 15, we
c ensure that kk loops over the values 1 to 6 - previously, the upper limit was
c set to 5, which appears to be an error.
c SJP 2003/04/28
c
c $Log: ocforce.f,v $
c Revision 1.17  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.16  2001/02/28 04:36:39  rot032
c Further tidy ups from HBG
c
c Revision 1.15  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.14  2001/02/12 05:39:53  rot032
c New river routing scheme and tidy-ups from HBG
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
c Revision 1.10  1997/12/23  00:23:38  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.9  1997/12/19  01:25:40  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.8.1.1  1997/12/19  02:03:17  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.8  1996/10/24  01:03:29  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.7  1994/08/08  16:22:48  ldr
c Fix up RCS header at top of file.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  lcouple - flag to run model in coupled mode, default F
c                  savefcor - if true, save stats for flux-corrections for
c                             coupled ocean model
c
c              from common/gausl in GAUSL.f
c                  w - normal guass weights
c
c              from common/rungrid in this subroutine
c                  landlg, landmg, minlgx  ) runoff relocation
c                  minmgx, spread, topland ) data
c
c              from common/timex in TIMEX.f
c                  dt - time step in seconds
c                  mstep - timestep in minutes
c                  ratlm - fraction of month elapsed, runs 0-1
c
c              from common/tm in this subroutine
c                  smo - flux correction from ocean to climatology
c
c              from arguments
c                  ind - if 0 = first pass
c
c     In/Out:  from common/aforce in this subroutine
c                  athf - heat flux to ocean
c                  athfa - heat flux to ocean with zero sub-ice heat input
c                  dsfw - P-E : fresh water to ocean (m/step)
c                  dsis - ice change : water flux to ocean (m/step)
c                  dsisb - sublimation from ice : water flux to ocean (m/step)
c                  drunof - runoff
c
c              from common/a2o in A2O.f
c                  osalf - ocean surface salt flux
c                  osurh - ocean surface heat flux
c     		     osrad - ocean surface incident radiation
c		     oswnd - ocean surface wind at 10
c                  osice - ice concentration 
c
c              from common/fcorr in this subroutine
c                  sfluxes - surface flxues for flux corrections
c
      subroutine ocforce(ind)
c
c  To compute ocean forcing terms - Heat and salinity.
c

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'OPARAMS.f'

C Argument list
      integer ind

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'BCOGCM.f'
      include 'A2O.f'
      include 'FCORR.f'
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'HOSING.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'
      include 'TTFC.f'

      real drunof,dsis,dsisb,dsfw,athf,athfa,atsi
      common/aforce/drunof(ln2,lat),dsis(ln2,lat),dsisb(ln2,lat)
     & ,dsfw(ln2,lat),athf(ln2,lat),athfa(ln2,lat),atsi(ln2,lat)

      include 'gas.h'

      real topland,spread
      integer landlg,landmg,minlgx,minmgx
      common/rungrid/topland(nlanp),spread(nlanp,8)
     & ,landlg(nlanp),landmg(nlanp),minlgx(ndips),minmgx(ndips)

      real atmvar
      common/CSIROAGCM/atmvar(lon,lat2,6) ! ima=lon, jma=lat2  to match MOM2

      real fwfpme,fwficea,fwficeb,fwfriv
      common/fwfdat/fwfpme(ln2,lat),fwficea(ln2,lat),fwficeb(ln2,lat),
     &              fwfriv(ln2,lat)

C Local work arrays and variables
      real atsf(ln2,lat)
      real rimsl(lon,0:lat2)
      real work1(ln2,lat)
      real salo(ln2, lat)
      real atsf1(ln2, lat)
      real atsf2(ln2, lat)
      real atsf3(ln2, lat)
      real atsfold(ln2, lat)
      real oflux(imt-2, jmt-2)

      integer i, it, j, kk, ktr, lg, lga, lx, ma, mg, mga, ns

      real dz1, runof, salice, saloold, topt

      parameter (salice = 0.01)
      parameter (saloold = 0.035)

#ifdef OCEAN_LOW

      real athfao(ln2, lat)
      real athfo(ln2, lat)
      real dsfwo(ln2, lat)
      real dsisbo(ln2, lat)
      real dsiso(ln2, lat)
      real drunofo(ln2, lat)

      integer mx

#endif

C Local data, functions etc
      integer mgsp(8),lgsp(8)
      data mgsp/1,1,0,-1,-1,-1,0,1/
      data lgsp/0,1,1,1,0,-1,-1,-1/

C Start code : ----------------------------------------------------------

#ifdef OCEAN_LOW

c.... R21 grid matching data :-

c.... 6 land points to be treated as ocean points, and 
c.... which get values by EW averaging.
c.... 1=Svarlbard, 2=Iceland, 3&4=Indon, 5=end of Chile, 6=Antarctic
c.... These are the (mg,lg,ns) coordinates
          integer islmg(6),isllg(6),islns(6)
          data islmg/ 3,62,23,20,52,54/
          data isllg/ 4, 8,26,27,12, 8/
          data islns/ 1, 1, 1, 2, 2, 2/

C20  *********. ***********       .         .     **  .         .  **
C21  ***********  ********        .         .     **  .         .  **
C22  ***********..***.***..........................*...............**
C23  **********   **  ***         .         .       **.         .  **
C24  *********.    *   **         .         .        **         .  **
C25  **********    *   L*         .         .         * **      .  **
C26  *********.        L.  S      .         .         .****     .  **
C27  ..*******.........L................................****.........
C28    ****** .        *LL*       .         .         .******   .
C-------------------------------EQUATOR------------------------------
C28    ****** .        **LSS  SS  .         .         .*******  .
C27    L***** .         S     *S  .         .         .******** .
C26     ***** .         .     **  .         .         .******** .
C25  ...*****..................L........................*******......
C24     ****  .         .   ** *  .         .         . ******* .
C23     ****L*.         .  *****  .         .         . ******  .
C22     ****L*.         . ******  .         .         . ******  .
C21     ****L*.         . ******* .         .         . *****   .
C20  ...***...............*******.......................*****........
C19     ***   .         . ******* .         .         . ****    .

c.... There are 5 land points to be treated as ocean points, but
c.... cannot get values by EW averaging. They will take values from 
c.... adjacent points. They are all in the Indonesian region.
c.... These are the (mg,lg,ns) coordinates :
          integer indmg(5),indlg(5),indns(5)
          data indmg/22,23,26,27,27/
          data indlg/28,28,28,28,27/
          data indns/ 2, 2, 2, 2, 2/
c.... Indonesian point 1 gets data from point to the South,
c.... point 2 from point to the East,
c.... point 3 from point to the West,
c.... point 4 from point to the East,
c.... point 5 from point to the East :
          integer indmgx(5),indlgx(5)
          data indmgx/ 0,+1,-1,+1,+1/
          data indlgx/-1, 0, 0, 0, 0/
c.... End of R21 grid matching data ....

#else

c...  When using the high-resolution version of the ocean model, there is no
c...  mismatch between the positions of the coastlines on the AGCM and OGCM
c...  grids. As a result, no grid-matching operations need to be performed.

#endif

      real wg(lat2)
      logical first
      data first/.true./
      save wg,first

      IF(ind.eq.1)THEN

c---- Arrays for calculating salinity tendencies for ocean model
c---- ind=1 : set arrays to zero
      do 10 lg=1,lat
      do 10 mg=1,ln2
        drunof(mg,lg)=0.0
        dsis  (mg,lg)=0.0
        dsisb (mg,lg)=0.0
        dsfw  (mg,lg)=0.0
        athf  (mg,lg)=0.0
        athfa (mg,lg)=0.0
   10 continue

      return

      END IF

      IF(ind.eq.2)THEN

#ifdef OCEAN_LOW

        if(mw.eq.22)then ! R21 only

c.... Set initial values on OGCM grid equal to values on AGCM grid
          do lg = 1, lat
            do mg = 1, ln2
              dsiso(mg, lg) = dsis(mg, lg)
              dsisbo(mg, lg) = dsisb(mg, lg)
              dsfwo(mg, lg) = dsfw(mg, lg)
              athfo(mg, lg) = athf(mg, lg)
              athfao(mg, lg) = athfa(mg, lg)
            end do
          end do

c.... Create values for the ocean model points which are
c.... land points in the agcm.
CSJP          do 15 kk=1,5
          do 15 kk=1,6
          ns=islns(kk)
          ma=islmg(kk)+(ns-1)*lon
          lg=isllg(kk)
          dsiso (ma,lg)=0.5*(dsis (ma+1,lg)+dsis (ma-1,lg))
          dsisbo(ma,lg)=0.5*(dsisb(ma+1,lg)+dsisb(ma-1,lg))
          dsfwo (ma,lg)=0.5*(dsfw (ma+1,lg)+dsfw (ma-1,lg))
          athfo (ma,lg)=0.5*(athf (ma+1,lg)+athf (ma-1,lg))
   15     athfao(ma,lg)=0.5*(athfa(ma+1,lg)+athfa(ma-1,lg))
          do 16 kk=1,5
          ns=indns(kk)
          ma=indmg(kk)+(ns-1)*lon
          mx=ma+indmgx(kk)
          lg=indlg(kk)
          lx=lg+indlgx(kk)
          dsiso (ma,lg)=dsis (mx,lx)
          dsisbo(ma,lg)=dsisb(mx,lx)
          dsfwo (ma,lg)=dsfw (mx,lx)
          athfo (ma,lg)=athf (mx,lx)
   16     athfao(ma,lg)=athfa(mx,lx)

c.... Apply conservation to these fields
          call conserve(dsis, dsiso, 1)
          call conserve(dsisb, dsisbo, 1)
          call conserve(dsfw, dsfwo, 1)
          call conserve(athf, athfo, 1)
          call conserve(athfa, athfao, 1)

        endif ! R21 only

C---- Compute Salinity Forcing

c.... At R21 resolution, there are 18 gridpoints that are treated as ocean by
c.... the AGCM, but as land by the OGCM. Estimate runoff for each of these
c.... points, using E-W averaging, N-S averaging or by taking values from
c.... adjacent points, as appropriate.
      if (mw .eq. 22) then

c.... Set initial values on OGCM grid equal to values on AGCM grid
        do lg = 1, lat
          do mg = 1, ln2
            drunofo(mg, lg) = drunof(mg, lg)
          end do
        end do

c.... Estimate runoff for additional land points
        drunofo(91, 16) = 0.5 * (drunof(91, 15) + drunof(91, 17))
        drunofo(72, 21) = 0.5 * (drunof(71, 21) + drunof(73, 21))
        drunofo(72, 22) = 0.5 * (drunof(71, 22) + drunof(73, 22))
        drunofo(72, 23) = 0.5 * (drunof(71, 23) + drunof(73, 23))
        drunofo(91, 25) = 0.5 * (drunof(91, 24) + drunof(91, 26))
        drunofo(67, 27) = drunof(68, 27)
        drunofo(85, 28) = 0.5 * (drunof(84, 28) + drunof(86, 28))
        drunofo(20, 28) = (2.0 * drunof(19, 28) + drunof(22, 28)) / 3.0
        drunofo(21, 28) = (drunof(19, 28) + 2.0 * drunof(22, 28)) / 3.0
        drunofo(19, 27) = 0.75 * drunof(19, 28) + 0.25 * drunof(19, 24)
        drunofo(19, 26) = 0.5 * (drunof(19, 28) + drunof(19, 24))
        drunofo(19, 25) = 0.25 * drunof(19, 28) + 0.75 * drunof(19, 24)
        drunofo(23, 18) = 0.5 * (drunof(22, 18) + drunof(24, 18))
        drunofo(64, 17) = 0.5 * (drunof(64, 18) + drunof(64, 16))
        drunofo(52, 9) = 0.5 * (drunof(52, 10) + drunof(52, 8))
        drunofo(53, 9) = 0.5 * (drunof(53, 10) + drunof(53, 8))
        drunofo(35, 8) = 0.5 * (drunof(34, 8) + drunof(36, 8))
        drunofo(52, 5) = 0.5 * (drunof(52, 6) + drunof(52, 4))

c.... Apply conservation to this field
        call conserve(drunof, drunofo, 2)

      end if

#endif

c---- Relocate land based runoff to ocean points
c---- Add to Freshwater flux at each ocean point
c----
c---- Transfer runoff data to a global array
      do 30 ns=1,2
      do 30 lg=1,lat
      lx=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 30 mg=1,lon
      ma=mg+(ns-1)*lon
   30 rimsl(mg,lx)=drunof(ma,lg)

c---- create global G-weights
      if(first)then
        do 28 lg=1,lat
        wg(lg)=w(lg)
   28   wg(lat2+1-lg)=w(lg)
        first=.false.
      end if

      IF(newriver)THEN

c---- Send runoff down rivers/reservoirs with time delay to the oceans.
        call newriv(rimsl,wg)

      ELSE

c---- Relocate the runoff by downslope method working from top down.
      ktr=0
      do 32 kk=1,nlanp
      topt=topland(kk)
      lg=landlg(kk)
      mg=landmg(kk)
      runof=rimsl(mg,lg)*wg(lg)
      if(topt.lt.0.0)then
c....   for troughs - relocate runoff to lower nearby point 
        ktr=ktr+1
        mga=minmgx(ktr)
        lga=minlgx(ktr)
        rimsl(mga,lga)=rimsl(mga,lga)+runof/wg(lga)
      else
c....  otherwise spread runoff downslope to any lower surrounding points
        do 34 it=1,8
        if(spread(kk,it).gt.0.0)then
         mga=mg+mgsp(it)
         lga=lg+lgsp(it)
         if(mga.eq.0)mga=lon
         if(mga.eq.lon+1)mga=1
         if(lga.eq.0)then
          lga=1
          mga=mga+lon/2
          if(mga.gt.lon)mga=mga-lon
         end if
        rimsl(mga,lga)=rimsl(mga,lga)+runof*spread(kk,it)/wg(lga)
        endif
   34   continue
      endif
      rimsl(mg,lg)=0
   32 continue

      ENDIF

c---- Save the Fresh water flix P-E stats
      if(fwf_sflg)then
        do lg=1,lat
        do mg=1,ln2
          fwfpme(mg,lg)=fwfpme(mg,lg)+dsfw(mg,lg) ! m of water per step
          work1(mg,lg)=dsfw(mg,lg) ! Save dsfw (before runoff added)
        enddo
        enddo
      endif

c...  Apply freshwater hosing if required
      if (hosing_flag) call hose(rimsl)

c---- Add the relocated runoff to the freshwater flux
      do 36 ns=1,2
      do 36 lg=1,lat
      lx=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 36 mg=1,lon
      ma=mg+(ns-1)*lon
   36 dsfw(ma,lg)=dsfw(ma,lg)+rimsl(mg,lx)

c---- Now compute the salinity tendency due to atmos forcing
c
c.... We now use actual ocean SSSs to calculate the salinity tendency, applying
c.... a uniform global correction in order to ensure conservation of freshwater
c
c.... The salinity of ice (salice) is 1% fixed.
c.... The ocean salinity (salo) is set at 0.035 to ensure conservation
c.... of ocean salinity given a zero annual mean fresh water input
c
c.... The salinity tendency is given by (Mk3) DS/Dt = atsf =
c.... = {Dis*(Sice-So) + DIsb*Sice - Dfw*So}/Dt/dz1
c....
c.... All depths in meters of ice/water
c.... Ocean salinities are in units # 0.035

c.... Calculate the actual ocean SSSs for this timestep
      do ns = 1, 2
        do lg = 1, lat
          lx = (ns-1) * lg + (lat2+1-lg) * (2-ns)
          do mg = 1, lon
            ma = mg + (ns-1) * lon
            if (lcouple) then
              salo(ma, lg) = 0.5 * (asal(mg, lx, 1) +
     &                              asal(mg, lx, 2))
            else
              salo(ma, lg) = (1.0-ratlm) * asal(mg, lx, 1) +
     &                       ratlm * asal(mg, lx, 2)
            end if
          end do
        end do
      end do
 
c     dzlev1 : depth of 1st ocean level (meters) : see TTFC.f
c      (If not coupled then dzlev1 must be prespecified here
c       for relevant ocean to be coupled to) :
      if(.not.lcouple)then
        if(mw.eq.64)then
          dzlev1=10.0        ! Mk3 T63/2 31 level ocean model
c  NOTE: The MOM2 code (see atmos.F) has now been coded so that
c   dzlev1=10 must always be used HERE for Mk3.
c   (regardless of the thickness of the 1st ocean model level)
        elseif(mw.eq.22)then
          dzlev1=25.0        ! Mk2 R21 12 level ocean model
        else
c         dzlev1=??.??
          print *,'Please specify dzlev1 for resolution mw=',mw,
     &            ' in ocforce.f'
          stop
        endif
      endif
      dz1=dzlev1

c.... Calculate the components of the salinity tendency, and sum them in order
c.... to obtain initial values for the total salinity tendency. We also
c.... calculate the total salinity tendency the old way (using a uniform
c.... surface salinity of 35 psu) for comparison.
      do lg = 1, lat
        do mg = 1, ln2
          atsf1(mg, lg) = dsis(mg, lg) * (salice - salo(mg, lg)) / 
     &                    (dt*dz1)
          atsf2(mg, lg) = dsisb(mg, lg) * salice / (dt*dz1)
          atsf3(mg, lg) = 0.0 - dsfw(mg, lg) * salo(mg, lg) / (dt*dz1)
          atsf(mg, lg) = atsf1(mg, lg) + atsf2(mg, lg) + atsf3(mg, lg)
          atsfold(mg, lg) = (dsis(mg, lg) * (salice - saloold) +
     &                       dsisb(mg, lg) * salice -
     &                       dsfw(mg, lg) * saloold) / (dt*dz1)
        end do
      end do

c.... Apply a uniform global correction, in order to ensure conservation of
c.... freshwater
#ifdef OCEAN_LOW
      call conserve_fw(dsis, dsisb, dsfw, atsf)
#else
      call conserve_fw2(dsis, dsisb, dsfw, atsf)
#endif

      if(fwf_sflg)then

c.... These statistics are redundant, so we calculate them the old way, using
c.... a uniform value for SALO. This avoids any possibility of divison by zero
c.... taking place.
c
c---- Save the Fresh water flux stats from Ice and from River outflow 
c----  (see above for equivalent components used in the salinity forcing)
c---- units are m of water per step

        do lg = 1, lat
          do mg = 1, ln2
            fwficea(mg, lg) = fwficea(mg, lg) + dsis(mg, lg) *
     &                        (saloold - salice) / saloold
            fwficeb(mg, lg) = fwficeb(mg, lg) - dsisb(mg,lg) * salice /
     &                                                         saloold
            fwfriv(mg, lg) = fwfriv(mg, lg) + (dsfw(mg, lg) -
     &                                         work1(mg, lg))
          end do
        end do

      end if

c---- Gather Ocean Heat and Salinity Forcing, and Solar

      if(savefcor)then
      do 20 ns=1,2
      do 20 lg=1,lat
       lx=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 20 mg=1,lon
        ma=mg+(ns-1)*lon
        sfluxes(mg,lx,3)=sfluxes(mg,lx,3)+athf(ma,lg)
        sfluxes(mg,lx,4)=sfluxes(mg,lx,4)+atsf(ma,lg)
        sfluxes(mg,lx,5)=sfluxes(mg,lx,5)+atsi(ma,lg)
        sfluxes(mg,lx,7)=sfluxes(mg,lx,7)+athfa(ma,lg)
        sfluxes(mg,lx,8)=sfluxes(mg,lx,8)+dsis(ma,lg)/dt
        sfluxes(mg,lx,9)=sfluxes(mg,lx,9)+dsisb(ma,lg)/dt
        sfluxes(mg,lx,10)=sfluxes(mg,lx,10)+dsfw(ma,lg)/dt
        sfluxes(mg,lx,11)=sfluxes(mg,lx,11)+atsf1(ma,lg)
        sfluxes(mg,lx,12)=sfluxes(mg,lx,12)+atsf2(ma,lg)
        sfluxes(mg,lx,13)=sfluxes(mg,lx,13)+atsf3(ma,lg)
        sfluxes(mg,lx,14)=sfluxes(mg,lx,14)+atsfold(ma,lg)
   20 continue
      endif

      if(lcouple)then

      do 50 ns=1,2
      do 50 lg=1,lat
       lx=(ns-1)*lg+(lat2+1-lg)*(2-ns)
      do 50 mg=1,lon
        ma=mg+(ns-1)*lon
c.. MOM2 data (averaging done by MOM2 code)
        atmvar(mg,lx,3)=athf(ma,lg)
        atmvar(mg,lx,4)=atsf(ma,lg)
        atmvar(mg,lx,5)=atsi(ma,lg)
   50 continue

#ifdef OCEAN_LOW

      do j = 1, lat2
        if (j .gt. lat) then
          lg = lat2 + 1 - j
        else
          lg = j
        end if
        do i = 1, lon
          if (j .gt. lat) then
            ma = i
          else
            ma = i + lon
          end if
          osurh(i, j) = osurh(i, j) + athf(ma, lg)
          osalf(i, j) = osalf(i, j) + atsf(ma, lg)
          osrad(i, j) = osrad(i, j) + atsi(ma, lg)
          oswnd(i, j) = oswnd(i, j) + wind10(ma, lg)
	   osice(i, j) = osice(i, j) + sice(ma, lg)
        end do
      end do

#else

c...  Use bilinear interpolation to downscale the surface heat flux and surface
c...  salinity tendency to the resolution of the new ocean model grid

      call downscale(athf, oflux)
      do j = 1, jmt-2
        do i = 1, imt-2
          osurh(i, j) = osurh(i, j) + oflux(i, j)
        end do
      end do

      call downscale(atsf, oflux)
      do j = 1, jmt-2
        do i = 1, imt-2
          osalf(i, j) = osalf(i, j) + oflux(i, j)
        end do
      end do

      call downscale(atsi, oflux)
      do j = 1, jmt-2
        do i = 1, imt-2
          osrad(i, j) = osrad(i, j) + oflux(i, j)
        end do
      end do

      call downscale(wind10, oflux)
      do j = 1, jmt-2
        do i = 1, imt-2
          oswnd(i, j) = oswnd(i, j) + oflux(i, j)
        end do
      end do

      call downscale(sice, oflux)
      do j = 1, jmt-2
        do i = 1, imt-2
          osice(i, j) = osice(i, j) + oflux(i, j)
        end do
      end do

#endif

      end if   ! lcouple

      END IF   ! ind

      return
      end
C---------------------------------------------------------------------
      subroutine newriv(runoff,wg)

C This is the new river routing scheme ; Currently for T63 only

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real runoff(lon,0:lat2) ! == rimsl() == drunof() im m water/step
      real wg(lat2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'LSMI.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'

      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

      real gradwe(lon,lat2),gradsn(lon,0:lat2)
C At T63, there are 0:309 points treated as troughs
      integer ntrpnt(0:309)
      integer ntrmg(0:309),ntrlg(0:309)
      integer ntre(0:309),ntrn(0:309),ntrw(0:309),ntrs(0:309)
      logical basin(0:309)
      real distf(0:309)
      common/rivflow/gradwe,gradsn
     & ,ntrpnt,ntrmg,ntrlg,ntre,ntrn,ntrw,ntrs,basin,distf

C Local work arrays and variables
      real resvr2(0:lon+1,0:lat2)
      real factr(0:309)
      character*1 ch1(lon,lat2)

      integer k
      integer kk
      integer lg
      integer lgH
      integer lgO
      integer lgn
      integer lgns
      integer lgs
      integer loops
      integer ma
      integer maxnum
      integer mg
      integer mge
      integer mgH
      integer mgO
      integer mgp
      integer mgw
      integer ns

      real flowsn
      real flowwe
      real resx
      real srfarea
      real sumla
      real sumwg

      logical mapriv

C Local data, functions etc
      integer Lakes(4)
      data Lakes/194,195,185,166/ ! Great Lakes trough indicators
      real flowrate
      data flowrate/0.001/
      character*1 chset(20)
      character*20 ch20
      equivalence (ch20,chset)
      data ch20/'+123456789ABCDEFGHIJ'/

C Start code : ----------------------------------------------------------

c Routine to move the runoff downstream with time delay

c Extend runoff array "below" South Pole using opposite points
c  and the reservoir array
         do mg=1,lon
           mgp=mg+lon/2
           if(mgp.gt.lon)mgp=mgp-lon
           runoff(mg,0)=runoff(mgp,1)
           resvr1(mg,0)=resvr1(mgp,1)
         enddo

c Add runoff to the reservoir (1)
        do lg=0,lat2
        do mg=1,lon
          resvr1(mg,lg)=resvr1(mg,lg)+runoff(mg,lg) ! m/step
c  Note that for the Lakes, the reservoir value may be < 0
          resx=resvr1(mg,lg) ! Prevent very small river values
          if((resx.gt.0.0).and.(resx.lt.1.0e-15))resvr1(mg,lg)=0.0
        enddo
        enddo
C Cyclic continuity
        do lg=0,lat2
          resvr1(0,lg)=resvr1(lon,lg)
          resvr1(lon+1,lg)=resvr1(1,lg)
        enddo

C Optional river mapping (set mapriv=.false. if not wanted)
c       mapriv=.true.
        mapriv=.false.
        if(mapriv.and.(nsteps.eq.1))then

          do lgns=1,lat2
           ns=2-(lgns-1)/lat
           lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
           do mg=1,lon
            ma=mg+(ns-1)*lon
            if(imsl(ma,lg).eq.4)then
             if(resvr1(mg,lgns).gt.0.0)then
               kk=min(1+nint(10.0*resvr1(mg,lgns)),20)
               ch1(mg,lgns)=chset(kk)
             endif
             if(resvr1(mg,lgns).eq.0.0)ch1(mg,lgns)='.'
             if(resvr1(mg,lgns).lt.0.0)ch1(mg,lgns)='-'
            else
             ch1(mg,lgns)=' '
            endif
           enddo
          enddo

c---- print out has max 128 characters per line
          maxnum=128
          loops=(lon-1)/maxnum+1
          maxnum=lon/loops

          do kk=1,loops
           mgw=kk*maxnum
           mge=mgw-maxnum+1
           if(kk.eq.loops)mgw=lon
           print *,'River Reservoirs: +,0(.),-'
           print *,'+=trace, -=Lake below sill, 1:J = 10*(metres),',
     &      ' J is >=2m'
           do lg=lat2,1,-1
            write(6,1617)lg,(ch1(mg,lg),mg=mge,mgw)
 1617 format(i2,1x,128a1) 
           enddo
          enddo

        endif
c--> End Mapping

C Copy contents to temporary reservoir (2)
        do lg=0,lat2
        do mg=0,lon+1
          resvr2(mg,lg)=resvr1(mg,lg)
          resvr1(mg,lg)=max(resvr1(mg,lg),0.0) ! No flow if < 0
        enddo
        enddo

c Compute the gradient flow :
C Reservoir (2) will be updated by using the gradient flow
C  between the same grid box values held in reservoir(1)

c W->E gradients between triangular regions
c Each T63 box has 15x15 1/8 degree pixels

c      +123456789012345+123456789012345+
c      5..............x|x..............|
c      4.............xx|xx.............|
c      3............xxx|xxx............|
c      2...........xxxx|xxxx...........|
c      1..........xxxxx|xxxxx..........|
c      0.........xxxxxx|xxxxxx.........|
c      9........xxxxxxx|xxxxxxx........|
c      8.....W.xxxxxxxx|xxxxxxxx.E.....|
c      7........xxxxxxx|xxxxxxx........|
c      6.........xxxxxx|xxxxxx.........|
c      5..........xxxxx|xxxxx..........|
c      4...........xxxx|xxxx...........|
c      3............xxx|xxx............|
c      2.............xx|xx.............|
c      1..............x|x..............|
c      +++++++++++++++++++++++++++++++++
c           (mg,lg)    |   (mg+1,lg)
c                      |
c                      |
c             gradsn(mg,lg)

c      gridT63=1.875*pi*erad/180.0 ! grid length in m (# 208Km)
c       sumw = sum of topographic heights in W triangle etc
c         gradwe(mg,lg)=(sumw-sume)/(gridT63*coslat(lg))
c           then converted to Miller et al flow velociity, then scaled
c       For W->E flow between (mg) and (mg+1)

c S->N gradients between triangular regions

c      +123456789012345+
c      5...............|
c      4...............|
c      3...............|
c      2...............|
c      1...............|
c      0.......N.......|
c      9...............|
c      8.......x.......| (mg,lg+1)
c      7......xxx......|
c      6.....xxxxx.....|
c      5....xxxxxxx....|
c      4...xxxxxxxxx...|
c      3..xxxxxxxxxxx..|
c      2.xxxxxxxxxxxxx.|
c      1xxxxxxxxxxxxxxx|
c      +123456789012345+ <-- gradsn(mg,lg)
c      5xxxxxxxxxxxxxxx|
c      4.xxxxxxxxxxxxx.|
c      3..xxxxxxxxxxx..|
c      2...xxxxxxxxx...|
c      1....xxxxxxx....|
c      0.....xxxxx.....|
c      9......xxx......|
c      8.......x.......| (mg,lg)
c      7...............|
c      6.......S.......|
c      5...............|
c      4...............|
c      3...............|
c      2...............|
c      1...............|
c      +++++++++++++++++

c       sumn = sum of topographic heights in N triangle etc
c         gradsn(mg,lg)=(sums-sumn)/gridT63,
c           then converted to Miller et al flow velociity, then scaled
c       For S->N flow between (lg) and (lg+1)

        do lg=1,lat2-1 ! lat2 is N Pole

c-- Next loops which were (1,lon) now broken into (1,lon-1) + final bit 
c       so as to vectorise -- >>
         do mg=1,lon-1
          flowwe=gradwe(mg,lg)
          if(flowwe.gt.0.0)then
            resvr2(mg+1,lg)=resvr2(mg+1,lg)+flowwe*resvr1(mg,lg)
            resvr2(mg,lg)  =resvr2(mg,lg)  -flowwe*resvr1(mg,lg)
          endif
         enddo
         do mg=1,lon-1
          flowwe=gradwe(mg,lg)
          if(flowwe.lt.0.0)then
            resvr2(mg+1,lg)=resvr2(mg+1,lg)+flowwe*resvr1(mg+1,lg)
            resvr2(mg,lg)  =resvr2(mg,lg)  -flowwe*resvr1(mg+1,lg)
          endif
         enddo
          flowwe=gradwe(lon,lg)
          if(flowwe.gt.0.0)then
            resvr2(1,lg)  =resvr2(1,lg)  +flowwe*resvr1(lon,lg)
            resvr2(lon,lg)=resvr2(lon,lg)-flowwe*resvr1(lon,lg)
          endif
          if(flowwe.lt.0.0)then
            resvr2(1,lg)  =resvr2(1,lg)  +flowwe*resvr1(1,lg)
            resvr2(lon,lg)=resvr2(lon,lg)-flowwe*resvr1(1,lg)
          endif
c-- >>

         do mg=1,lon
          flowsn=gradsn(mg,lg)
          if(flowsn.gt.0.0)then
            resvr2(mg,lg+1)=resvr2(mg,lg+1)+flowsn*resvr1(mg,lg)
     &                                          *wg(lg)/wg(lg+1)
            resvr2(mg,lg)  =resvr2(mg,lg)  -flowsn*resvr1(mg,lg)
          endif
          if(flowsn.lt.0.0)then
            resvr2(mg,lg+1)=resvr2(mg,lg+1)+flowsn*resvr1(mg,lg+1)
            resvr2(mg,lg)  =resvr2(mg,lg)  -flowsn*resvr1(mg,lg+1)
     &                                            *wg(lg+1)/wg(lg)
          endif
         enddo

        enddo
c
c Take care of possible flow across South Pole (lg=0 and lg=1)
c
        lg=0
        do mg=1,lon
          flowsn=gradsn(mg,lg)
          if(flowsn.gt.0.0)then
            resvr2(mg,lg+1)=resvr2(mg,lg+1)+flowsn*resvr1(mg,lg)
          endif
          if(flowsn.lt.0.0)then
            resvr2(mg,lg+1)=resvr2(mg,lg+1)+flowsn*resvr1(mg,lg+1)
          endif
        enddo
c Update lg=0 row (using "opposite" points)
        do mg=1,lon
          mgp=mg+lon/2
          if(mgp.gt.lon)mgp=mgp-lon
          resvr2(mg,lg)=resvr2(mgp,lg+1)
        enddo

C Now replace reservoir(1) values by those held 
C  (temporarily) in reservoir(2)
        do lg=0,lat2
        do mg=1,lon
          resvr1(mg,lg)=resvr2(mg,lg)
        enddo
        enddo
C Add cyclic continuity
        do lg=0,lat2
          resvr1(0,lg)=resvr1(lon,lg)
          resvr1(lon+1,lg)=resvr1(1,lg)
          resvr2(0,lg)=resvr2(lon,lg)
          resvr2(lon+1,lg)=resvr2(1,lg)
        enddo
 

C Transfer from troughs to another precalculated point
C (Update reservoir(1) values by means of (same) reservoir(2) values)

c Compute amount removed from each trough - place in factr()
        do kk=0,309 ! Number of troughs
         mg=ntrmg(kk)
         lg=ntrlg(kk)
c Note: distf() = distance factor is pre-computed in landrun
c The "flowrate" can be adjusted as required.
c (If G.Lakes are below sill depth (i.e. resvr2<0) then prevent any flow)
         factr(kk)=(flowrate/distf(kk))*max(0.0,resvr2(mg,lg))
        enddo

c Now adjust the troughs (mg,lg) and the downstream points
        do kk=0,309 ! Number of troughs
         mg=ntrmg(kk)
         lg=ntrlg(kk)
         lgn=lg+ntrn(kk)
         lgs=lg-ntrs(kk)
         mge=mg+ntre(kk)
         mgw=mg-ntrw(kk)

         resvr1(mg,lg)=resvr1(mg,lg)-factr(kk)

         if(ntrpnt(kk).eq.1)resvr1(mge,lg)=resvr1(mge,lg)+
     &                      factr(kk)                        ! E flow

         if(ntrpnt(kk).eq.2)resvr1(mge,lgn)=resvr1(mge,lgn)+
     &                      factr(kk)*wg(lg)/wg(lgn)         ! NE flow

         if(ntrpnt(kk).eq.3)resvr1(mg,lgn)=resvr1(mg,lgn)+
     &                      factr(kk)*wg(lg)/wg(lgn)         ! N flow

         if(ntrpnt(kk).eq.4)resvr1(mgw,lgn)=resvr1(mgw,lgn)+
     &                      factr(kk)*wg(lg)/wg(lgn)         ! NW flow

         if(ntrpnt(kk).eq.5)resvr1(mgw,lg)=resvr1(mgw,lg)+
     &                      factr(kk)                        ! W flow

         if(ntrpnt(kk).eq.6)resvr1(mgw,lgs)=resvr1(mgw,lgs)+
     &                      factr(kk)*wg(lg)/wg(lgs)         ! SW flow

         if(ntrpnt(kk).eq.7)resvr1(mg,lgs)=resvr1(mg,lgs)+
     &                      factr(kk)*wg(lg)/wg(lgs)         ! S flow

         if(ntrpnt(kk).eq.8)resvr1(mge,lgs)=resvr1(mge,lgs)+
     &                      factr(kk)*wg(lg)/wg(lgs)         ! SE flow

        enddo

c Adjust the Great Lakes
c
c   Great Lakes (L)
c      144 145 146 147 148 149 150 151 152 153
c     .........................................
c (lg).   .  Superior .   .   .   .   .   .   .(lgns)
c  23 .   . L . L . L .   .   .   .   .   .   .  74
c     .   .194.   .195.   .   .   .   .   .   .
c     .........................................
c     .   .   .   M   .   Huron   .   .   .   .
c  24 .   .   .   i L . L . L . L .   .   .   .  73
c     .   .   .   c   .   .   .185.   .   .   .
c     ............h............................
c     .   .   .   i   .   .   .   .   Ontario .
c  25 .   .   .   g L .   . L .   .   . L .   .  72
c     .   .   .   a   .   .   .   .   .176.   .
c     ............n............................
c     .   .   .   .   .   . Erie  .   .   .   .
c  26 .   .   .   .   .   . L . L .   .   .   .  71
c     .   .   .   .   .   .166.   .   .   .   .
c     .........................................
c
c Trough 194=Lake Superior(W)
c Trough 195=Lake Superior(E)
c Trough 185=Lake Huron
c Trough 166=Lake Erie
c   Niagara Falls
c Trough 176=Lake Ontario
c (See landrun.f for slower transfer rates out of troughs 166 and 176)
c
c First remove any water used to top up the G.Lakes "puddle" depth
c (done in surfa). This will be removed initially from L.Huron trough.
c  (and then spread to other Great Lakes troughs - see below)
c  This may make the G.Lakes attain a negative value.
c  This is equivalent to the G.Lakes being below sill depth
c  and outflow is prevented (see above). This also implies that the 
c  G.Lakes have unlimited water available for evaporation.
         kk=Lakes(3) ! L.Huron = 185
         mgH=ntrmg(kk)
         lgH=ntrlg(kk)
         sumla=0.0
         do lg=23,26
         lgns=lat2+1-lg
         do mg=145,150
           if(lakesg(mg,lg))then
             sumla=sumla+resrem(mg,lgns)*wg(lgns)/wg(lgH)
           endif
         enddo
         enddo
         resvr1(mgH,lgH)=resvr1(mgH,lgH)-sumla

c Lake Ontario also needs water removed (not included in the above)
         mgO=ntrmg(176)
         lgO=ntrlg(176)
         resvr1(mgO,lgO)=resvr1(mgO,lgO)-resrem(mgO,lgO)

c The reservoirs (troughs) 194,195,185,166 are now averaged to
c  give the same amount in each
         sumla=0.0
         sumwg=0.0
         do k=1,4
          kk=Lakes(k)
          mg=ntrmg(kk)
          lg=ntrlg(kk)
          sumla=sumla+resvr1(mg,lg)*wg(lg)
          sumwg=sumwg+wg(lg)
         enddo
         sumla=sumla/sumwg

         do k=1,4
          kk=Lakes(k)
          mg=ntrmg(kk)
          lg=ntrlg(kk)
          resvr1(mg,lg)=sumla
         enddo

c Adjust for Lake Winnipeg evap
c
c      139 140 141 142
c     .................
c (lg).   .   .   .   .(lgns)
c  19 .   .   .tr .   .  78
c     .   .   .221.   .
c     .................
c     .   Winnipeg.   .
c  20 .   . L . L .   .  77
c     .   .   .   .   .
c     .................
c       (puddle water out of trough 221)
         sumla=(resrem(140,77)+resrem(141,77))*wg(77)
         resvr1(141,78)=resvr1(141,78)-sumla/wg(78)
         
c Adjust for "Baker Lakes" evap
c
c      139 140 141 142 143
c     .....................
c (lg).   .Baker Lakes.   .(lgns)
c  14 .   . L . L .tr .   .  83
c     .   .   .   .273.   .
c     .....................
c       (puddle water out of trough 273)
         sumla=resrem(140,83)+resrem(141,83)
         resvr1(142,83)=resvr1(142,83)-sumla

c Adjust for southern Hudson Bay lakes
c
c      147 148 149 150 151 152 153
c     ........|___.....___|........
c (lg).   .   .   |   |   .   .   .(lgns)
c  20 .   .   . L |   | L .   .   .  77
c     .   .   .   |   |   .   .   .
c     ............|___|............
c     .   .   .   .   .   .   .   .
c  21 .   .tr .   . L .   .   .   .  76
c     .   .209------> .   .   .   .
c     .............................
c       (trough 209 flows to middle reservoir/lake point)
c Take puddle water out of each lake. Also average the
c  3 Hudson lakes to give the same amount in each.
         resvr1(149,77)=resvr1(149,77)-resrem(149,77)
         resvr1(151,77)=resvr1(151,77)-resrem(151,77)
         resvr1(150,76)=resvr1(150,76)-resrem(150,76)
         sumla=(resvr1(149,77)+resvr1(151,77))*wg(77)
     &                   +resvr1(150,76)*wg(76)
         sumwg=2*wg(77)+wg(76)
         sumla=sumla/sumwg

         resvr1(149,77)=sumla
         resvr1(151,77)=sumla
         resvr1(150,76)=sumla
         
c End of lakes adjustments.

c Return the non-land resvr1() in runoff.
c (units of m water per step to match drunof == rimsl)
c This is the runoff that has finally reached the ocean from the 
c  rivers/reservoirs.
c Then set the ocean part of resvr1() to zero
      do lgns=1,lat2
       ns=2-(lgns-1)/lat
       lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
       do mg=1,lon
        runoff(mg,lgns)=resvr1(mg,lgns)
        ma=mg+(ns-1)*lon
        if(imsl(ma,lg).eq.4)then
          runoff(mg,lgns)=0.0
        else
          resvr1(mg,lgns)=0.0
        endif
       enddo
      enddo


      if(rsv_sflg)then
        do lg=1,lat2
        do mg=1,lon
c Collect stats on ocean outflow 
c (m of water per step converted to 1000 m**3/sec)
          srfarea=erad*erad*wg(lg)*2.0*pi/lon
          srunoff(mg,lg)=srunoff(mg,lg)+runoff(mg,lg)*srfarea/dt/1000
C Collect stats on the average reservoir amount (m of water)
          sresvr1(mg,lg)=sresvr1(mg,lg)+resvr1(mg,lg)
        enddo
        enddo
      endif

      return
      end
