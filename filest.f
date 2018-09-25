c Modified for five-character experiment names.
c SJP 2009/04/22
c
c nsib comparsion is modified to lsm_type = "nsib "
c AJA 2009/01/22
c
c SJP 2004/09/29
c
c (1) Declaration of /FCORR/ moved to FCORR.f.
c (2) Modified so that 14, rather than 7, fields are saved when SAVEFCOR=T.
c SJP 2004/09/10
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c Modified for the changes to the freshwater fluxes, whereby the ice water flux
c is split into its two components. The masking of the fluxes over land is
c also removed - this was erroneous, as the freshwater fluxes are calculated
c on the OGCM grid.
c SJP 2003/06/18
c
c Modified so that DTM can be saved without having to set SAVEFCOR=T.
c SJP 2003/05/29
c
c $Log: filest.f,v $
c Revision 1.54  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.53  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.52.1.1  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.52  2001/03/07 04:28:58  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.51  2001/02/22 05:34:38  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.50  2001/02/12 05:39:47  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.49.1.1  2001/02/22 05:56:37  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.49  2000/12/08 03:58:52  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.48  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.47  2000/06/20 02:27:09  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.46  1999/06/30 05:29:37  rot032
c New comments from HBG.
c
c Revision 1.45  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.44  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.43  1998/12/10  00:55:33  ldr
c HBG changes to V5-1-21
c
c Revision 1.42  1997/12/19  02:03:11  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.41  1997/07/24  05:42:49  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.40  1997/06/11  02:21:29  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.39  1996/06/13  02:06:29  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.38  1996/04/17  01:32:52  ldr
c Changes from TIE to allow initialization of monthly mean netcdf files
c from within model.
c
c Revision 1.37  1996/03/21  03:18:40  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.36  1995/11/23  06:03:28  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.35  1995/08/08  02:02:14  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.34  1995/06/30  02:44:40  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.33  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.32  1994/08/08  17:21:15  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.31  94/08/08  13:13:01  ldr
c Add flag c_sflg for cloud stats.
c 
c Revision 1.30  94/08/04  17:06:57  ldr
c Merge of LDR and HBG changes.
c 
c Revision 1.29  94/08/04  10:37:07  ldr
c Comment out diagnostic of convective cloud.
c 
c Revision 1.28.1.1  94/08/04  16:55:15  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.28.1.2  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.28  94/03/30  12:34:20  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.27  93/12/17  15:32:29  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.26  93/11/03  13:06:52  ldr
c Replace silly flag zfl_sflg with savezflux.
c 
c Revision 1.25  93/11/03  11:53:36  ldr
c Merge of HBG with LDR changes since V4-4-24l.
c 
c Revision 1.24  93/10/19  10:31:40  ldr
c Add extra diagnostics for qcloud scheme.
c 
      subroutine filest(str)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      character str*13

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'BGRND.f'
      include 'ECPARM.f'
      include 'ECTRAC.f'
      include 'FCORR.f'
      include 'FEWFLAGS.f'
      include 'LSMI.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'
      include 'VERTV.f'        !Coupled model

      real hlat_a,clat_a
      common/latheat/hlat_a(ln2,lat,nl),clat_a(ln2,lat,nl)

      real tgmax,tgmin,tsmax,tsmin,htsmax,htsmin
     & ,extsmax,extsmin,tggmax,tggmin,tgfmax,tgfmin
     & ,htggmax,htggmin,htgfmax,htgfmin,hpmc
      common/maxmin/tgmax(ln2,lat),tgmin(ln2,lat),
     &              tsmax(ln2,lat),tsmin(ln2,lat),
     &              htsmax(ln2,lat),htsmin(ln2,lat),
     &              extsmax(ln2,lat),extsmin(ln2,lat),
     &              tggmax(ln2,lat),tggmin(ln2,lat),
     &              tgfmax(ln2,lat),tgfmin(ln2,lat),
     &              htggmax(ln2,lat),htggmin(ln2,lat),
     &              htgfmax(ln2,lat),htgfmin(ln2,lat)
     &             ,hpmc(ln2,lat)

      real fwfpme,fwficea,fwficeb,fwfriv
      common/fwfdat/fwfpme(ln2,lat),fwficea(ln2,lat),fwficeb(ln2,lat),
     &              fwfriv(ln2,lat)

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

C Local work arrays and variables
      real work1(ln2,lat),work2(ln2,lat),work3(ln2,lat),work4(ln2,lat)

      integer k
      integer ke
      integer kx
      integer lg
      integer lgns
      integer ma
      integer mg
      integer ns
      integer nv

      real avg
      real consd
      real consr
      real consx
      real conx

      character char2*2
      character varname*3

C Local data, functions etc

C Start code : ----------------------------------------------------------

C
C FUNCTION : TO ARCHIVE MONTHLY PHYSICAL DATA OF CSIRO MODEL
C         AND DYNAMICAL DATA ON P-SURFACES
C

C****
C**** TO OUTPUT PHYSICAL + DYNAMICAL STATS TO AT END OF MONTH
C****
      consd=1.0/(nsteps*mstep/1440.0) ! 1/(no. of days in run)
      consx=1.0/nsteps
      consr=float(nrad)/nsteps
C**** con?? ARE SCALINGS FOR THE PHYSICAL DATA GATHERED
C**** (PHYSICAL DATA FROM 14 TO 50 IN RADSTM ARRAYS,
C****  WHILST 1-13 ARE NOT STATS

c---- The data collected is for the full global field unless the
c---- resolution is T63 when some fields will have a 9:1 compression
c---- See ncinit.f where compression is set at 1 or 9.

c WRITE SURFACE FLUXES TO BE USED TO CALCULATE FLUX CORRECTIONS FOR
c COUPLED MODEL.

      if(savefcor)then
        write(70) nsteps
        do nv=1,14

cxxxx Note : dimension of fort.70 changed for Mk3.1 model [ Now has (u*)**3 ]
c....        Increased again to include ATHFA. SJP 2003/06/20
c....        Increased to add a further 7 fields. SJP 2004/09/10
c
c     sfluxes(:, :,  1) = taux (N/m**2)
c     sfluxes(:, :,  2) = tauy
c     sfluxes(:, :,  3) = Total ocean heat flux (w/m**2)
c     sfluxes(:, :,  4) = Total ocean salinity tendency
c     sfluxes(:, :,  5) = Solar into ocean (w/m**2)
c     sfluxes(:, :,  6) = (u*)**3=(sqrt(taux**2+tauy**2)/1025.0)**(3/2)
c     sfluxes(:, :,  7) = Total ocean heat flux with zero sub-ice heat input
c     sfluxes(:, :,  8) = DSIS (water flux arising from ice growth/melt)
c     sfluxes(:, :,  9) = DSISB (water flux arising from ice sublimation)
c     sfluxes(:, :, 10) = DSFW (P-E plus runoff)
c     sfluxes(:, :, 11) = ATSF1 (the salinity tendency arising from DSIS)
c     sfluxes(:, :, 12) = ATSF2 (the salinity tendency arising from DSISB)
c     sfluxes(:, :, 13) = ATSF3 (the salinity tendency arising from DSFW)
c     sfluxes(:, :, 14) = ATSFOLD (ATSF calculated using S_o = 35psu)

          write(70)((consx*sfluxes(mg,lg,nv),mg=1,lon),lg=1,lat2)

        enddo
CSJP        if (dtm_sflg) call collst(52,consx,str)         ! (Tmlo-Tintp)
      endif

c.... It's safe to save DTM, even if SAVEFCOR=F.
      if (dtm_sflg) call collst(52,consx,str)         ! (Tmlo-Tintp)

c---- Monthly statistical files that may be generated (listed alphabetically)
c      sacd[yymm.???] = Average drag                     N/m2
c      sals[yymm.???] = Surface albedo                   0 - 1
c      sc[01:nl][yymm.???] = Cloud per model level       0 - 1
c      sclc[yymm.???] = Convective cloud                 0 - 1
c      scld[yymm.???] = Total cloud                      0 - 1
c      scli[yymm.???] = Liquid cloud fraction            0 - 1
c      sclh[yymm.???] = High cloud                       0 - 1
c      scll[yymm.???] = Low cloud                        0 - 1
c      sclm[yymm.???] = Middle cloud                     0 - 1
c      sdtm[yymm.???] = MLO temp deviation               K
c      sevp[yymm.???] = Evaporation                      mm/day
c      sfw1[yymm.???] = Fresh water flux : P-E           mm/day
c      sfw2[yymm.???] = Fresh water flux : Ice water A   mm/day
c      sfw3[yymm.???] = Fresh water flux : Ice water B   mm/day
c      sfw4[yymm.???] = Fresh water flux : River outflow mm/day
c      sg[01:nl][yymm.???] = GPH per level               m       (P-lev)
c      shfl[yymm.???] = Sensible heat flux               W/m2
c      sicd[yymm.???] = Ice concentration * depth        m
c      sich[yymm.???] = Ice advection                    m
c      sico[yymm.???] = Ice concentration                0 - 1
c      sicu[yymm.???] = Ice zonal velocity               m/s
c      sicv[yymm.???] = Ice merid velocity               m/s
c      sinr[yymm.???] = Interception of rain by foliage  mm/day
c      sire[yymm.???] = Ice redistribution               m
c      sisf[yymm.???] = Ice to ocean salt flux           mm/day
c      sitf[yymm.???] = Ice to ocean heat flux           W/m2
c      siwp[yymm.???] = Frozen (atmos) water path        Kgm/m2
c      sl[01:nl][yymm.???] = Latent heating per model level
c      slwp[yymm.???] = Liquid water path                Kgm/m2
c      sper[yymm.???] = Soil percolation
c      spev[yymm.???] = Potential evaporation            mm/day
c      spmc[yymm.???] = Moisture in puddles              mm
c      spsl[yymm.???] = Sea level pressure               hPa
c      spwc[yymm.???] = Precipitable water content       mm
c      sq[01:nl][yymm.???] = Moisture per level          Kgm/Kgm (P-lev)
c      sr[01:nl][yymm.???] = RH per level                0 - 100 (P-lev)
c      sref[yymm.???] = Effective radius for liquid clouds 0 - 20
c      srev[yymm.???] = Rain evaporation                 mm/day
c      srgc[yymm.???] = Net clear sky LW at ground       W/m2
c      srgd[yymm.???] = Downward LW at ground            W/m2
c      srgn[yymm.???] = Net LW at ground                 W/m2
c      srnc[yymm.???] = Convective rain                  mm/day
c      srnd[yymm.???] = Total rain                       mm/day
c      srsv[yymm.???] = River/reservoirs                 m
c      srtc[yymm.???] = Clear sky LW at top              W/m2
c      srtu[yymm.???] = LW outwards at top               W/m2
c      srun[yymm.???] = Runoff                           mm/day
c      srvo[yymm.???] = River outflow (as rainfall)      m/month
c      ssev[yymm.???] = Scaling evaporation              mm/day
c      ssgc[yymm.???] = Net clear sky SW at ground       W/m2
c      ssgd[yymm.???] = Downward SW at ground            W/m2
c      ssgn[yymm.???] = Net SW at ground                 W/m2
c      ssid[yymm.???] = Sea ice depth                    m
c      ssnd[yymm.???] = Snow depth                       m
c      ssno[yymm.???] = Snowfall                         mm/day
c      ssoc[yymm.???] = Clear sky SW out at top          W/m2
c      ssot[yymm.???] = SW outwards at top               W/m2
c      sssb[yymm.???] = Snow sublimation                 mm/day
c      st[01:nl][yymm.???] = Temperature per level       K      (P-lev)
c      stax[yymm.???] = Surface stress Eastward          N/m2
c      stay[yymm.???] = Surface stress Northward         N/m2
c      stb2[yymm.???] = Soil temperature level 2         K
c      stb3[yymm.???] = Soil temperature level 3 (deep)  K
c      stgf[yymm.???] = Vegetated ground temperature     K
c      stgg[yymm.???] = Bare ground temperature          K
c      sthd[yymm.???] = Daily max high screen temp       K
c      sthf[yymm.???] = Daily max vegetated ground temp  K
c      sthg[yymm.???] = Daily max bare ground temp       K
c      sthm[yymm.???] = Extreme high screen temp/month   K
c      stld[yymm.???] = Daily min low sceen temp         K
c      stlf[yymm.???] = Daily min vegetated ground temp  K
c      stlg[yymm.???] = Daily min bare ground temp       K
c      stlm[yymm.???] = Extreme low  screen temp/month   K
c      stsc[yymm.???] = Screen temperature               K
c      stsu[yymm.???] = Surface temperature              K
c      su[01:nl][yymm.???] = Zonal wind per level        m/s   (P-lev)
c      sv[01:nl][yymm.???] = Meridional wind per level   m/s   (P-lev)
c      svmo[yymm.???] = Surface wind speed               m/s
c      swdf[yymm.???] = Ice divergence removed by rheology +-100
c      swfb[yymm.???] = Lower level soil moisture        Fraction
c      swfg[yymm.???] = Upper level soil moisture        Fraction
c      swls[yymm.???] = Ice residual divergence          0 - 100
c     IF (sltrace) THEN
c       tt[01:nl][yymm.???] = Tracer t average conc      ppmv (sig-lev)
c       tf[01:nl][yymm.???] = Tracer f average conc      ppmv (sig-lev)
c       tv[01:nl][yymm.???] = Tracer v average conc      ppmv (sig-lev)
c     ENDIF
c     IF (coupled_aero) THEN
c       sdms[yymm.???] = dimethyl sulfide vertical integral mgS/m2
c       sso2[yymm.???] = sulfur dioxide vertical integral   mgS/m2
c       sso4[yymm.???] = sulfate vertical integral          mgS/m2
c       ss2d[yymm.???] = SO2 dry deposition                 kg/m2/s
c       ss4d[yymm.???] = SO4 dry deposition                 kg/m2/s
c       ss2w[yymm.???] = SO2 wet deposition                 kg/m2/s
c       ss4w[yymm.???] = SO4 wet deposition                 kg/m2/s
c       ssem[yymm.???] = Sulfur emissions                   kg/m2/s
c       ss2o[yymm.???] = Oxidation of SO2 by OH             kg/m2/s
c       ss2h[yymm.???] = Oxidation of SO2 by H202           kg/m2/s
c       ss23[yymm.???] = Oxidation of SO2 by O3             kg/m2/s
c       sdm1[yymm.???] = dimethyl sulfide level 1           ug/m3
c       ss21[yymm.???] = sulfur dioxide level 1             ug/m3
c       ss41[yymm.???] = sulfate level 1                    ug/m3
c     ENDIF
c----

c Sea-ice model stuff

      if(semice)then
        if ( ico_sflg ) then
           call collstx(10,consx,str) !ico
           call collstx(11,consx,str) !icd
        endif
        if ( itf_sflg ) call collstx(12,consx,str) !itf
        if ( isf_sflg ) call collstx(13,consx,str) !isf
      endif

C**
C**   OUTPUT PHYSICAL STATISTICS WITH EACH STATISTIC COLLECTED
C**    INTO SINGLE ARRAY FORM, WITH DATA ORDERED AS FOLLOWS :
C**     HSND,HSID,HTST,HTSCRN
C**    ,HWFG,HWFB,HRN,HEVAP,HFLUX
C**    ,HCLD,HCLL,HCLM,HCLH
C**    ,HRG,HRT,HSIN,HSG
C**     NEW RADIATION STATS: HRTCLR, HSOUTCLR at end (LDR 9/90)
C**
C      4 GENERAL SURFACE CONDITION VALUES : HSND=SNOW DEPTH;
C        HSID=SEA ICE DEPTH; HTST ; HTSCR=TSCREEN
      if ( snd_sflg ) call collst(29,consx,str)
      if ( sid_sflg ) call collst(30,consx,str)
      if ( tsu_sflg ) call collst(25,consx,str)
      if ( tsc_sflg ) call collst(26,consx,str)
C      5 SURFACE HYDROLOGY RELATED VALUES : WFCG;WFCB;RAINFALL;
C       HEVAP;HFLUX
      if ( wfg_sflg ) then
         call collst(27,consx,str)
         varname='pmc'             ! Puddle depth
         call collph(varname,consx,str,hpmc)
      endif
      if ( wfb_sflg ) call collst(28,consx,str)
      if ( rnd_sflg ) call collst(15,consd,str)
      if ( evp_sflg ) call collst(14,consd,str)
      if ( hfl_sflg ) call collst(16,consx,str)
C     TOTAL CLOUD + 3 CLOUD TYPES.
      if ( cld_sflg ) call collst(17,consr,str)
      if ( cll_sflg ) call collst(18,consr,str)
      if ( clm_sflg ) call collst(19,consr,str)
      if ( clh_sflg ) call collst(20,consr,str)

C Diagnostics for LDR cloud scheme.

      if ( qcloud )then
        if ( clc_sflg )call collst(51,consx,str) !Conv cloud
        if ( lwp_sflg )then
          call collst(53,consx,str) !Liquid water path 
          call collst(71,consx,str) !Ice water path 
        endif
        if ( pwc_sflg )call collst(54,consx,str) !Precipitable water
        if ( rev_sflg )call collst(55,consd,str) !Evap of rain
        if ( ssb_sflg )call collst(56,consd,str) !Subl of snow
        if ( sno_sflg )call collst(57,consd,str) !Snowfall
        if ( ref_sflg )then
          call collst(58,consr,str) !Reff for liq. clouds
          call collst(59,consr,str) !Liq. cloud fraction to scale Reff
        endif
      endif

      if (coupled_aero) then
        call collst(60,consx,str) !DMS vertically integrated
        call collst(61,consx,str) !SO2 vertically integrated
        call collst(62,consx,str) !SO4 vertically integrated
        call collst(63,consx,str) !SO2 dry deposition, 10^(-12) kg/m2/s
        call collst(64,consx,str) !SO4 dry deposition, ditto
        call collst(65,consx,str) !SO2 wet deposition, ditto
        call collst(66,consx,str) !SO4 wet deposition, ditto
        call collst(67,consx,str) !S   emissions,      ditto
        call collst(68,consx,str) !SO2 to SO4 ox. by OH, ditto
        call collst(69,consx,str) !SO2 to SO4 ox. by H202, ditto
        call collst(70,consx,str) !SO2 to SO4 ox. by O3, ditto
        call collst(72,consx,str) !DMS level 1 (ug/m3)
        call collst(73,consx,str) !SO2 level 1 (ug/m3)
        call collst(74,consx,str) !SO4 level 1 (ug/m3)
        call collst(75,consx,str) !DMS to SO2 ox. by OH,  10^(-12) kg/m2/s
        call collst(76,consx,str) !DMS to SO2 ox. by NO3, ditto
        call collst(77,consx,str) !Total convective wet deposition of S
      endif

C     LWR VALUES(RG,RT);SOLAR INPUT;SG.
      if ( rgn_sflg ) call collst(21,consx,str)
      if ( rtu_sflg ) call collst(22,consr,str)
      if ( sot_sflg ) call collst(23,consr,str)
      if ( sgn_sflg ) call collst(24,consx,str)

c NEW PHYSICAL STATS (clear sky LW and SW outward fluxes, 
c soil temperatures tb2,tb3)

      if ( rtc_sflg ) call collst(31,consr,str) !hrtclr
      if ( soc_sflg ) call collst(32,consr,str) !hsoutclr
      if ( tb2_sflg ) call collst(33,consx,str) !htb2
      if ( tb3_sflg ) call collst(34,consx,str) !htb3

c MONTHLY MEAN SURFACE WIND SPEED AND SURFACE STRESSES
      
      if ( vmo_sflg ) call collst(35,consx,str) !vmod
      if ( tax_sflg ) call collst(36,consx,str) !taux
      if ( tay_sflg ) call collst(37,consx,str) !tauy

c CONVECTIVE RAIN

      if ( rnc_sflg ) call collst(38,consd,str)

c Some NSIB fields

!!bp      if( lsm_type .eq. "nsib " ) then
       if ( run_sflg ) call collst(39,consd,str) ! runoff
       if ( int_sflg ) call collst(40,consd,str) ! interception
       if ( als_sflg ) call collst(41,consx,str) ! albedo
       if ( tgg_sflg ) call collst(42,consx,str) ! bare ground temp
       if ( tgf_sflg ) call collst(43,consx,str) ! vegetated ground temp
       if ( pev_sflg ) call collst(44,consd,str) ! potential evaporation
       if ( sev_sflg ) call collst(45,consd,str) ! scaling evaporation
       if ( per_sflg ) call collst(50,consd,str) ! deep soil moisture percolation
!!bp      endif

c SURFACE CLOUD FORCING
      if ( rgc_sflg ) call collst(46,consr,str) !hrgclr
      if ( sgc_sflg ) call collst(47,consr,str) !hsgclr

c SURFACE DOWNWARD LW AND SW
      if ( rgd_sflg ) call collst(48,consr,str) !hrgdn
      if ( sgd_sflg ) call collst(49,consr,str) !hsgdn

c OCEAN FRESH WATER FLUX COMPONENTS (mms/day)
      if ( fwf_sflg ) then
       do lg=1,lat
       do mg=1,ln2
CSJP         if(imsl(mg,lg).eq.4)then  !   mask out land points
CSJP           work1(mg,lg)=-7777777.0
CSJP           work2(mg,lg)=-7777777.0
CSJP           work3(mg,lg)=-7777777.0
CSJP         else                      ! m/month to mms/day
           work1(mg,lg)=fwfpme(mg,lg)*consd*1000.0
           work2(mg,lg)=fwficea(mg,lg)*consd*1000.0
           work3(mg,lg)=fwficeb(mg,lg)*consd*1000.0
           work4(mg,lg)=fwfriv(mg,lg)*consd*1000.0
CSJP         endif
       enddo
       enddo
       varname='fw1'             ! P-E
       call collph(varname,1.0,str,work1)
       varname='fw2'             ! Ice water A
       call collph(varname,1.0,str,work2)
       varname='fw3'             ! Ice water B
       call collph(varname,1.0,str,work3)
       varname='fw4'             ! River outflow
       call collph(varname,1.0,str,work4)
      endif

C
C     OUTPUT MONTHLY SUM OF MERIDIONAL TRANSPORT STATS TO TAPE51.
C       (THESE ARE ONLY COLLECTED 4 TIMES A DAY)
C     Done using nifty equivalenced array PVBTBX

      if ( savezflux ) then
        open(51,file='zflx'//str,form='formatted')
        ke=lat*2*(7*nl+nl-1)
        avg=1.0/kdynm
        do 30 k=1,ke
 30     pvbtbx(k)=pvbtbx(k)*avg
        write(51,9040)pvbtb,pvtb,pvbub,pvub,pvbqb,pvqb,psv,pssd
 9040   format(1p6e12.5)
        close(51)
      endif
C****
C**** PROCESS THE DYNAMICAL DATA (GATHERED EVERY 1/4 DAY)
C****
      write(6,100)kdynm
  100 format(1h ,'averaging p-surface dynamical data -kdynm=',i4)
      avg=1.0/kdynm
      do 31 k=nbelow+1,nl
      do 31 lg=1,lat
      do 31 mg=1,ln2
      tmonth(mg,lg,k)=tmonth(mg,lg,k)*avg
      umonth(mg,lg,k)=umonth(mg,lg,k)*avg
      vmonth(mg,lg,k)=vmonth(mg,lg,k)*avg
      rmonth(mg,lg,k)=rmonth(mg,lg,k)*avg
      gmonth(mg,lg,k)=gmonth(mg,lg,k)*avg
   31 qmonth(mg,lg,k)=qmonth(mg,lg,k)*avg
      do 33 lg=1,lat
      do 33 mg=1,ln2
   33 pmonth(mg,lg)=pmonth(mg,lg)*avg
C**   SPECIAL CHECKS FOR LEVEL 1,2,3,4 DATA : MAY BE
C**   NON-EXISTENT (OR NEAR NON EXISTENT) VALUES DUE TO
C**   BELOW GROUND LEVEL VALUES. IBETA() HOLDS THE NUMBER
C**   OF ADDITIONS OF ABOVE SURFACE VALUES. IF IBETA()=KDYNM
C**   THEN VALUES EXISTED ALWAYS. SELECT KDYNM/4 (ARBITRARY)
C**   AS CUT OFF FOR INCLUSION OF DATA . REPLACE UNWANTED
C**   VALUES WITH -7777777 TO IDENTIFY.
      do 35 k=1,nbelow
      do 35 lg=1,lat
      do 35 mg=1,ln2
      if(ibeta(mg,lg,k).lt.(kdynm/4)) then
        tmonth(mg,lg,k)=-7777777.0
        umonth(mg,lg,k)=-7777777.0
        vmonth(mg,lg,k)=-7777777.0
        rmonth(mg,lg,k)=-7777777.0
        gmonth(mg,lg,k)=-7777777.0
        qmonth(mg,lg,k)=-7777777.0
        ibeta(mg,lg,k)=1
      end if
   35 continue
      do 36 k=1,nbelow
      do 36 lg=1,lat
      do 36 mg=1,ln2
      tmonth(mg,lg,k)=tmonth(mg,lg,k)/ibeta(mg,lg,k)
      umonth(mg,lg,k)=umonth(mg,lg,k)/ibeta(mg,lg,k)
      vmonth(mg,lg,k)=vmonth(mg,lg,k)/ibeta(mg,lg,k)
      rmonth(mg,lg,k)=rmonth(mg,lg,k)/ibeta(mg,lg,k)
      gmonth(mg,lg,k)=gmonth(mg,lg,k)/ibeta(mg,lg,k)
   36 qmonth(mg,lg,k)=qmonth(mg,lg,k)/ibeta(mg,lg,k)

c Write dynamical stats (T,U,V,Rh,Q,GPH,Pmsl)
c Allow two different formats: traditional ASCII and new NetCDF.

      if ( t_sflg ) then
        do 82 k=1,nl
        write(char2,'(i2.2)')k
        varname='t'//char2
 82     call colldy(varname,1.0,str,tmonth(1,1,k))
      endif
      
      if ( u_sflg ) then
        do 84 k=1,nl
        write(char2,'(i2.2)')k
        varname='u'//char2
 84     call colldy(varname,1.0,str,umonth(1,1,k))
      endif
      
      if ( v_sflg ) then
        do 86 k=1,nl
        write(char2,'(i2.2)')k
        varname='v'//char2
 86     call colldy(varname,1.0,str,vmonth(1,1,k))
      endif
      
      if ( rh_sflg ) then
        do 87 k=1,nl
        write(char2,'(i2.2)')k
        varname='r'//char2
 87     call colldy(varname,1.0,str,rmonth(1,1,k))
      endif

      if ( q_sflg ) then
        do 88 k=1,nl
        write(char2,'(i2.2)')k
        varname='q'//char2
 88     call colldy(varname,1.0,str,qmonth(1,1,k))
      endif
      
      if ( g_sflg ) then
        do 89 k=1,nl
        write(char2,'(i2.2)')k
        varname='g'//char2
 89     call colldy(varname,1.0,str,gmonth(1,1,k))
      endif
      
      varname='psl'
c.... must use k=1 for PS
      if ( psl_sflg ) call colldy(varname,1.0,str,pmonth(1,1))
      
c     Daily max and min tscreen
      
      varname='thd'             ! T High Daily
      if ( thd_sflg ) call collph(varname,consd,str,htsmax)
      varname='tld'             ! T Low Daily
      if ( tld_sflg ) call collph(varname,consd,str,htsmin)
      
c     Daily max and min tgg (bare ground temperature)
      
      varname='thg'             ! TGG High Daily
      if ( thg_sflg ) call collph(varname,consd,str,htggmax)
      varname='tlg'             ! TGG Low Daily
      if ( tlg_sflg ) call collph(varname,consd,str,htggmin)
      
c     Daily max and min tgf (vegetated ground temperature)
      
      varname='thf'             ! TGF High Daily
      if ( thf_sflg ) call collph(varname,consd,str,htgfmax)
      varname='tlf'             ! TGF Low Daily
      if ( tlf_sflg ) call collph(varname,consd,str,htgfmin)
      
c     Monthly extreme max and min tscreen
      
      varname='thm'             ! T High Month
      if ( thm_sflg ) call collph(varname,1.0,str,extsmax)
      varname='tlm'             ! T Low Month
      if ( tlm_sflg ) call collph(varname,1.0,str,extsmin)

c     New river routing
      if(newriver.and.rsv_sflg)then

        do lg=1,lat
        do ns=1,2
         lgns=(ns-1)*lg+(lat2+1-lg)*(2-ns)
         ma=(ns-1)*lon
         do mg=1,lon
c work1 = River outflow (1000 m**3/s) ; work2 = Avgerage reservoir (m)
          work1(mg+ma,lg)=srunoff(mg,lgns)*consx
          work2(mg+ma,lg)=sresvr1(mg,lgns)*consx
         enddo
        enddo
        enddo

        do lg=1,lat
        do mg=1,ln2
         if(imsl(mg,lg).eq.4)then
           work1(mg,lg)=-7777777.0 ! mask out land points
         else
           work2(mg,lg)=-7777777.0 ! mask out sea points
         endif
        enddo
        enddo

        varname='rvo'             ! River outflow (1000 m**3/sec)
        call collph(varname,1.0,str,work1)

        varname='rsv'             ! Average rivers/reservoirs
        call collph(varname,1.0,str,work2)

      endif
      
c     Latent heating
      
      if ( l_sflg ) then
        conx=0.5*consd          !For latent heating
        do 90 k=1,nl
        write(char2,'(i2.2)')k
        varname='l'//char2
 90     call colldy(varname,conx,str,hlat_a(1,1,k))
      endif

c Cloud amounts

      if ( c_sflg ) then
c Global cloud is stored for all levels for consistency with qcloud scheme.
c Total cloud, low, middle, and high cloud collected earlier
c in this routine.

        do 95 k=1,nl
        write(char2,'(i2.2)')k
        varname='c'//char2
 95     call colldy(varname,consr,str,clat_a(1,1,k))
      endif

      if ( coupled_aero ) then

        do k=1,nl
          write(char2,'(i2.2)')k
          varname='d'//char2  !DMS
          call colldy(varname,consx,str,hxtg(1,1,k,1))
        enddo

        do k=1,nl
          write(char2,'(i2.2)')k
          varname='o'//char2  !SO2
          call colldy(varname,consx,str,hxtg(1,1,k,2))
        enddo

        do k=1,nl
          write(char2,'(i2.2)')k
          varname='s'//char2  !SO4
          call colldy(varname,consx,str,hxtg(1,1,k,3))
        enddo

      endif

      kx=nsteps*mstep/1440
      write(6,120)kx,nsteps
  120 format(1h ,i2,' day averaged physical stats and dynamical'
     &,' fields (on p surfaces) written to file51, nsteps=',i6)
      return
      end
