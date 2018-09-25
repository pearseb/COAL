c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Updated for the new version of the subroutine hstring.
c SJP 2008/11/20
c
c Modified for five-digit year numbers.
c SJP 2004/09/29
c
c Preset the pressure levels on which to output dagnostics, in accordance with
c the pressure levels on which data must be provided for inclusion in PMIP2.
c SJP 2004/02/02
c
c (1) Modified for the changes to the freshwater fluxes, whereby the ice water
c flux is split into its two components. Valid range of FW1 increased from
c (-40, 100) to (-100, 100).
c (2) Valid range of RUN increased from (0, 25) to (0, 100).
c SJP 2003/06/19
c
c Valid range of DTM increased from (-10, 10) to (-30, 30).
c SJP 2003/06/07
c
c Modified so that DTM can be saved without having to set SAVEFCOR=T.
c SJP 2003/05/29
c
c Specify the precision of all integer arrays and real variables that are
c used as arguments to netCDF routines. This ensures portability of the code
c across 32- and 64-bit platforms.
c SJP 2001/12/13
c
c In mattrib and tmattrib, ind replaced with ind2 in the argument list. ind
c is then set equal to ind2 at the start of the subroutine. Otherwise, a
c segmentation fault occurs when the value of ind is changed.
c SJP 2001/11/23
c
c Updated from v2.4 to v3 of netCDF.
c SJP 2001/11/21
c
c $Log: ncinit.f,v $
c Revision 1.31  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.30  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.29.1.1  2001/10/12 02:13:44  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.29  2001/06/04 02:26:52  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.28  2001/03/07 04:28:57  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.27  2001/02/22 06:46:47  rot032
c Merge LDR and HBG changes.
c
c Revision 1.26  2001/02/12 05:39:44  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.25.1.1  2001/02/22 05:56:37  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.25  2000/12/08 03:58:52  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.24  2000/11/14 07:00:17  rot032
c Add details for sulfur-cycle fields.
c
c Revision 1.23  2000/11/14 03:55:01  rot032
c Merge HBG changes with ncshort change.
c
c Revision 1.22  2000/08/28 06:11:08  rot032
c Change nclong to ncshort for SGI, as per HBG.
c
c Revision 1.21.1.1  2000/11/14 03:11:35  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.21  2000/08/17 03:02:48  rot032
c Add initialization for IWP and ECHAM aerosol fields.
c
c Revision 1.20  2000/06/21 05:35:22  rot032
c Fix for SGI.
c
c Revision 1.19  2000/06/20 02:18:38  rot032
c Merge of HBG changes with LDR's LWP fixup.
c
c Revision 1.18  1999/07/20 02:37:29  rot032
c Correct units for LWP.
c
c Revision 1.17.1.1  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.17  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.16  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.15  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.14  1998/12/10  00:55:31  ldr
c HBG changes to V5-1-21
c
c Revision 1.13  1998/05/26  05:10:49  ldr
c Final 10k run changes (mainly ACH salinity stuff) merged into V5-1-9.
c
c Revision 1.12  1998/02/02  06:33:18  ldr
c Increase the limits on ire and ich by factor of 10 (suggested by SPO).
c
c Revision 1.6.1.4  1998/05/26  04:48:55  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.11  1998/01/30  04:43:20  ldr
c Merge SPO's Cd diagnostic fix into latest version.
c
c Revision 1.10  1998/01/08  05:46:17  ldr
c Correct the isf fix.
c
c Revision 1.9  1998/01/08  02:51:11  ldr
c MSLP fix.
c
c Revision 1.8  1998/01/08  01:44:40  ldr
c Correct units and range for ice-ocean salt flux (isf).
c
c Revision 1.7.1.1  1998/01/30  04:36:14  ldr
c Fix from SPO for Cd diagnostic.
c
c Revision 1.7  1997/12/19  02:03:10  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.6  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.5  1997/07/24  22:58:55  mrd
c Fix zero length string problem.
c
c Revision 1.4  1997/07/22  05:04:40  mrd
c Get rid of rubbish characters at end of string from hstring.
c
c Revision 1.3  1997/06/11  02:21:28  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.2  1996/08/02  00:29:40  mrd
c Changed check of file pre-exisitence to use correct inquire specifier.
c
c Revision 1.1  1996/04/17  01:31:11  ldr
c Initial revision
c
      subroutine ncinit
c
c     Create netcdf files for monthly means (TIE 4/96).
c     

      implicit none
C Global parameters
      include 'PARAMS.f'
      integer numcdf
      parameter (numcdf=500) ! Must be >= total number of NetCDF files

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'STFLAGS.f'

      real tpa,tpb,tpc,wpa,wpb,wpc,qpa,qpb,qpc,plev,fred
      common/dystb/tpa(nl),tpb(nl),tpc(nl),wpa(nl),wpb(nl),wpc(nl)
     &,qpa(nl),qpb(nl),qpc(nl),plev(nl),fred(nl)

      integer ncdfcom(numcdf),nncdf
      character*3 ncdfchr(numcdf)
      common/ncdfiles/nncdf,ncdfcom,ncdfchr

C Local work arrays and variables
      character  lname*50, vname*3, pname*6
      character*2 lev(nl)

      integer icomp
      integer itcom
      integer k

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Set icomp=compression of data factor - It is only used for T63,
c  and then for some variables only (see mattrib calls below).
c At T63, icomp=9 . If mattrib call uses icomp, then T63 data arrays 
c  are taken from (192,96) down to (64,32)  - this is a  9:1  compression.
c If no compression required, replace the icomp by 1 in mattrib call.
c This routine is the only place where icomp has to be set (Mk3 model)
      icomp=1
      if(lw.eq.64)icomp=9
c
c     Check if data collection flag is set to true
c
c     Single level Variables
c
      nncdf=0 ! Counter for number of NetCDF/ASCII files generated

C Surface statistics flags
      if (evp_sflg) then
         lname='evaporation' 
         call mattrib('evp', lname, 'mm/day', 1)
      endif
      if (pev_sflg) then
         lname='potential evap'
         call mattrib('pev', lname, 'mm/day', 1)
      endif
      if (sev_sflg) then
         lname='scaling evap'
         call mattrib('sev', lname, 'mm/day', 1)
      endif
      if (rnd_sflg) then
         lname='precipitation'
         call mattrib('rnd', lname, 'mm/day', 1)
      endif
      if (rnc_sflg) then
         lname='convective rainfall'
         call mattrib('rnc', lname, 'mm/day', 1)
      endif
      if (hfl_sflg) then
         lname='sensible heat flux'
         call mattrib('hfl', lname, 'W/m2', 1)
      endif
      if (wfg_sflg) then
         lname='soil moisture upper'
         call mattrib('wfg', lname, ' ', 1)
         lname='moisture puddles'
         call mattrib('pmc', lname, 'mm', 1)
      endif
      if (wfb_sflg) then
         lname='soil moisture lower'
         call mattrib('wfb', lname, ' ', 1)
      endif
      if (run_sflg) then
         lname='runoff'
         call mattrib('run', lname, 'mm/day', 1)
      endif
      if (per_sflg) then
         lname='Soil percolation'
         call mattrib('per', lname, 'mm/day', 1)
      endif
      if (int_sflg) then
         lname='interception'
         call mattrib('inr', lname, 'mm/day', 1)
      endif
      if (psl_sflg) then
         lname='mean sea-level pressure'
         call mattrib('psl', lname, 'hPa', 1)
      endif
      if (vmo_sflg) then
         lname='surface wind speed'
         call mattrib('vmo', lname, 'm/s', icomp)
      endif
      if (tax_sflg) then
         lname='surface stress east'
         call mattrib('tax', lname, 'N/m2', 1)
      endif
      if (tay_sflg) then
         lname='surface stress nth'
         call mattrib('tay', lname, 'N/m2', 1)
      endif
      if (tsu_sflg) then
         lname='surface temperature'
         call mattrib('tsu', lname, 'K', 1)
      endif
      if (tsc_sflg) then
         lname='screen temperature'
         call mattrib('tsc', lname, 'K', 1)
      endif
      if (tb2_sflg) then
         lname='soil temp level 2'
         call mattrib('tb2', lname, 'K', 1)
      endif
      if (tb3_sflg) then
         lname='soil temp lowest level'
         call mattrib('tb3', lname, 'K', 1)
      endif
      if (tgg_sflg) then
         lname='bare ground temp'
         call mattrib('tgg', lname, 'K', 1)
      endif
      if (tgf_sflg) then
         lname='veg ground temp'
         call mattrib('tgf', lname, 'K', 1)
      endif
      if (thd_sflg) then
         lname='mean daily max temp'
         call mattrib('thd', lname, 'K', 1)
      endif
      if (tld_sflg) then
         lname='mean daily min temp'
         call mattrib('tld', lname, 'K', 1)
      endif
      if (thg_sflg) then
         lname='bare ground Tmax'
         call mattrib('thg', lname, 'K', 1)
      endif
      if (tlg_sflg) then
         lname='bare ground Tmin'
         call mattrib('tlg', lname, 'K', 1)
      endif
      if (thf_sflg) then
         lname='Veg ground Tmax'
         call mattrib('thf', lname, 'K', 1)
      endif
      if (tlf_sflg) then
         lname='Veg ground Tmin'
         call mattrib('tlf', lname, 'K', 1)
      endif
      if (thm_sflg) then
         lname='extreme max temp'
         call mattrib('thm', lname, 'K', 1)
      endif
      if (tlm_sflg) then
         lname='extreme min temp'
         call mattrib('tlm', lname, 'K', 1)
      endif
CSJP      if (savefcor.and.dtm_sflg) then
      if (dtm_sflg) then
         lname='Tmlo error'
         call mattrib('dtm', lname, 'K', 1)
      endif

      if(newriver)then
       if(rsv_sflg)then
         lname='reservoirs'! river/resrvr amount stored in m (not mm)
         call mattrib('rsv', lname, 'm', 1)
         lname='river outflow' ! river outflow stored as 1000 m**3/sec
         call mattrib('rvo', lname, '1000 m3/s', 1)
       endif
      endif

C Cloud statistics
      if (cld_sflg) then
         lname='total cloud'
         call mattrib('cld', lname, ' ', icomp)
      endif
      if (cll_sflg) then
         lname='low cloud'
         call mattrib('cll', lname, ' ', icomp)
      endif
      if (clm_sflg) then
         lname='middle cloud'
         call mattrib('clm', lname, ' ', icomp)
      endif
      if (clh_sflg) then
         lname='high cloud'
         call mattrib('clh', lname, ' ', icomp)
      endif

C Radiation statistics
      if (rgn_sflg) then
         lname='Net LW at ground'
         call mattrib('rgn', lname, 'W/m2', 1)
      endif
      if (rgd_sflg) then
         lname='Downward LW ground'
         call mattrib('rgd', lname, 'W/m2', icomp)
      endif
      if (rgc_sflg) then
         lname='Net LW ground clear'
         call mattrib('rgc', lname, 'W/m2', icomp)
      endif
      if (sgn_sflg) then
         lname='Net SW at ground'
         call mattrib('sgn', lname, 'W/m2', 1)
      endif
      if (sgd_sflg) then
         lname='Downward SW ground'
         call mattrib('sgd', lname, 'W/m2', icomp)
      endif
      if (sgc_sflg) then
         lname='Net SW ground clear'
         call mattrib('sgc', lname, 'W/m2', icomp)
      endif
      if (rtu_sflg) then
         lname='LW out at top'
         call mattrib('rtu', lname, 'W/m2', icomp)
      endif
      if (rtc_sflg) then
         lname='LW out clear sky'
         call mattrib('rtc', lname, 'W/m2', icomp)
      endif
      if (sot_sflg) then
         lname='SW out at top'
         call mattrib('sot', lname, 'W/m2', icomp)
      endif
      if (soc_sflg) then
         lname='SW out clear sky'
         call mattrib('soc', lname, 'W/m2', icomp)
      endif
      if (als_sflg) then
         lname='surface albedo'
         call mattrib('als', lname, ' ', 1)
      endif

C Snow & ice statistics
      if (snd_sflg) then
         lname='snow depth'
         call mattrib('snd', lname, 'cm', 1)
      endif
      if (sid_sflg) then
         lname='sea-ice depth'
         call mattrib('sid', lname, 'm', 1)
      endif
      if (ico_sflg) then
         lname='ice concentration'
         call mattrib('ico', lname, ' ', 1)
         lname='ice depth * conc'
         call mattrib('icd', lname, 'm', 1)
      endif
      if (itf_sflg) then
         lname='ice-ocean heat flux'
         call mattrib('itf', lname, 'W/m2', 1)
      endif
      if (isf_sflg) then
         lname='ice-ocean salt flux'
         call mattrib('isf', lname, 'm', 1)
      endif
      if (icu_sflg) then
         lname='ice zonal velocity'
         call mattrib('icu', lname, 'm/s', 1)
      endif
      if (icv_sflg) then
         lname='ice merid velocity'
         call mattrib('icv', lname, 'm/s', 1)
      endif
      if (div_sflg) then
         lname='ice residual divergence'
         call mattrib('wls', lname, 's-1', 1)
         lname='divergence removed by rheology'
         call mattrib('wdf', lname, 's-1', 1)
      endif
      if (gro_sflg) then
         lname='Monthly ice growth'
         call mattrib('gro', lname, 'm', 1)
      endif
      if (ire_sflg) then
         lname='Ice redistribution'
         call mattrib('ire', lname, 'm', 1)
      endif
      if (ich_sflg) then
         lname='Ice advection'
         call mattrib('ich', lname, 'm', 1)
      endif

C Fresh water (fw) flux into ocean statistics
      if (fwf_sflg) then
         lname='Fresh water flux : Precip-Evap'
         call mattrib('fw1', lname, 'mm/day', 1)
         lname='Fresh water flux : Ice water A'
         call mattrib('fw2', lname, 'mm/day', 1)
         lname='Fresh water flux : Ice water B'
         call mattrib('fw3', lname, 'mm/day', 1)
         lname='Fresh water flux : River outflow'
         call mattrib('fw4', lname, 'mm/day', 1)
      endif

C Qcloud statistics
      IF (qcloud) THEN
      if (sno_sflg) then
         lname='Snowfall'
         call mattrib('sno', lname, 'mm/day', icomp)
      endif
      if (rev_sflg) then
         lname='Rain evaporation'
         call mattrib('rev', lname, 'mm/day', icomp)
      endif
      if (ssb_sflg) then
         lname='Snow sublimation'
         call mattrib('ssb', lname, 'mm/day', icomp)
      endif
      if (clc_sflg) then
         lname='convective cloud'
         call mattrib('clc', lname, ' ', icomp)
      endif
      if (lwp_sflg) then
         lname='liquid water path'
         call mattrib('lwp', lname, 'kg/m2', icomp)
         lname='ice water path'
         call mattrib('iwp', lname, 'kg/m2', icomp)
      endif
      if (pwc_sflg) then
         lname='precipitable water'
         call mattrib('pwc', lname, 'mm', icomp)
      endif
      if (ref_sflg) then
         lname='Effective radius for liquid clouds'
         call mattrib('ref', lname, 'um', icomp)
         lname='Liquid cloud fraction'
         call mattrib('cli', lname, ' ', icomp)
      endif
      ENDIF

c****************************************************************
c**** Next - Model data for levels in the atmosphere
c****************************************************************
      do k=1,nl
         write(lev(k),'(i2.2)')k
      enddo

c****************************************************************
c**** Data collected on pre-set pressure levels given by plev(nl)
c****************************************************************

c**** The pressure levels (pre-set) to which the model dynamics variables
c**** of temperature, winds, water vapour (U,V,T,Q) and the derived
c**** relative humidity (RH) are interpolated are set up here.
c**** Geopotential height (G) is also set up

c**** The data levels MAY be given by 1000mb*(Sigma-level-value) of
c**** the model being used (usual strategy)
CSJP      do k=1,nl
CSJP        plev(k)=1000.0*sig(k)
CSJP      enddo
c****    OR
c**** they may be preset to any sensible values BUT NOTE THAT
c**** when specifying arbitrary data levels (must have Nl of them)
c**** then the required pressure levels given by plev()
c**** must be in order and decreasing in pressure (mbs). 
CSJP  Set the levels to the standard 17 WMO pressure levels, as required for
CSJP  participation in PMIP2. Add an additional level at 5 hPa, as the model
CSJP  has 18 vertical levels.
      plev(1)  = 1000.0
      plev(2)  =  925.0
      plev(3)  =  850.0
      plev(4)  =  700.0
      plev(5)  =  600.0
      plev(6)  =  500.0
      plev(7)  =  400.0
      plev(8)  =  300.0
      plev(9)  =  250.0
      plev(10) =  200.0
      plev(11) =  150.0
      plev(12) =  100.0
      plev(13) =   70.0
      plev(14) =   50.0
      plev(15) =   30.0
      plev(16) =   20.0
      plev(17) =   10.0
      plev(18) =    5.0

c     Zonal Wind Variable (p-level)
c
      if (u_sflg) then
         do k=1,nl
            vname='u'//lev(k)
            write(pname,'(f6.1)')plev(k)
            lname='Zonal wind at pressure='//pname//' mbs'
            call mattrib(vname, lname, 'm/s', icomp)
         enddo   
      endif

c     Meridional Wind Variable (p-level)
c
      if (v_sflg) then
         do k=1,nl
            vname='v'//lev(k)
            write(pname,'(f6.1)')plev(k)
            lname='Meridional wind at pressure='//pname//' mbs'
            call mattrib(vname, lname, 'm/s', icomp)
         enddo   
      endif

c     Temperature Variable (p-level)
c
      if (t_sflg) then
         do k=1,nl
            vname='t'//lev(k)
            write(pname,'(f6.1)')plev(k)
            lname='Temperature at pressure='//pname//' mbs'
            call mattrib(vname, lname, 'K', icomp)
         enddo   
      endif

c     Specific Humidity Variable (p-level)
c
      if (q_sflg) then
         do k=1,nl
            vname='q'//lev(k)
            write(pname,'(f6.1)')plev(k)
            lname='Specific humidity at pressure='//pname//' mbs'
            call mattrib(vname, lname, 'kg/kg', icomp)
         enddo   
      endif  

c     Relative Humidity Variable (p-level)
c
      if (rh_sflg) then
         do k=1,nl
            vname='r'//lev(k)
            write(pname,'(f6.1)')plev(k)
            lname='Relative humidity at pressure='//pname//' mbs'
            call mattrib(vname, lname, ' ', icomp)
         enddo
      endif

c     Geopotential Height Variable (p-level)
c
      if (g_sflg) then
         do k=1,nl
            vname='g'//lev(k)
            write(pname,'(f6.1)')plev(k)
            lname='Geopotential height at pressure='//pname//' mbs'
            call mattrib(vname, lname, 'm', icomp)
         enddo
      endif

c****************************************************************
c**** Data collected on model 'sigma' levels or 'hybrid' levels
c****************************************************************

c     Cloud Variable (Sigma-level)
c
      if (c_sflg) then
         do k=1,nl
            vname='c'//lev(k)
            lname='Cloud at sigma level '//lev(k)
            call mattrib(vname, lname, ' ', icomp)
         enddo
      endif

c     Latent Heat Variable (Sigma-level)
c
      if (l_sflg) then
         do k=1,nl
            vname='l'//lev(k)
            lname='Latent heat at sigma level '//lev(k)
            call mattrib(vname, lname, 'K/day', icomp)
         enddo
      endif

c     Tracer (if any) Variables (Sigma-level)
c
      if (sltrace) then
c
c Data output within filewr note. T63 data has compression.
c
c Month average tracer concentrations:-tt??, tf??, tv??
c
         if(mw.eq.64)then
           itcom=9
         else
           itcom=1
         endif
         do k=1,nl
            vname='t'//lev(k)
            lname='Tracer t average conc at sigma level '//lev(k)
            call tmattrib(vname, lname, 'ppmv', itcom)
         enddo
         do k=1,nl
            vname='f'//lev(k)
            lname='Tracer f average conc at sigma level '//lev(k)
            call tmattrib(vname, lname, 'ppmv', itcom)
         enddo
         do k=1,nl
            vname='v'//lev(k)
            lname='Tracer v average conc at sigma level '//lev(k)
            call tmattrib(vname, lname, 'ppmv', itcom)
         enddo
c
c End of month tracer concentrations:-tT??, tF??, tV??
c
c        do k=1,nl
c           vname='T'//lev(k)
c           lname='Tracer t final conc at sigma level '//lev(k)
c           call tmattrib(vname, lname, 'ppmv', itcom)
c        enddo
c        do k=1,nl
c           vname='F'//lev(k)
c           lname='Tracer f final conc at sigma level '//lev(k)
c           call tmattrib(vname, lname, 'ppmv', itcom)
c        enddo
c        do k=1,nl
c           vname='V'//lev(k)
c           lname='Tracer v final conc at sigma level '//lev(k)
c           call tmattrib(vname, lname, 'ppmv', itcom)
c        enddo

      endif

c****************************************************************
c**** aerosol data
c****************************************************************

      IF (coupled_aero) THEN
         lname='dimethyl sulfide vertical integral'
         call mattrib('dms', lname, 'mgS/m2', 1)
         lname='sulfur dioxide vertical integral'
         call mattrib('so2', lname, 'mgS/m2', 1)
         lname='sulfate vertical integral'
         call mattrib('so4', lname, 'mgS/m2', 1)
         lname='SO2 dry deposition'
         call mattrib('s2d', lname, 'kg/m2/s', 1)   !Really 10^-12 kg/m2/s for all these
         lname='SO4 dry deposition'
         call mattrib('s4d', lname, 'kg/m2/s', 1)
         lname='SO2 wet deposition'
         call mattrib('s2w', lname, 'kg/m2/s', 1)
         lname='SO4 wet deposition'
         call mattrib('s4w', lname, 'kg/m2/s', 1)
         lname='Total S convective wet deposition'
         call mattrib('scw', lname, '10^-12 kg/m2/s', 1)
         lname='Sulfur emissions'
         call mattrib('sem', lname, 'kg/m2/s', 1)
         lname='Oxidation of SO2 by OH'
         call mattrib('s2o', lname, 'kg/m2/s', 1)
         lname='Oxidation of SO2 by H202'
         call mattrib('s2h', lname, 'kg/m2/s', 1)
         lname='Oxidation of SO2 by O3'
         call mattrib('s23', lname, 'kg/m2/s', 1)
         lname='Oxidation of DMS by OH'
         call mattrib('dmo', lname, '10^(-12) kg/m2/s', 1)
         lname='Oxidation of DMS by NO3'
         call mattrib('dmn', lname, '10^(-12) kg/m2/s', 1)
         lname='dimethyl sulfide level 1'
         call mattrib('dm1', lname, 'ug/m3', 1)
c         call mattrib('dm1', lname, 'pptv', 1)
         lname='sulfur dioxide level 1'
         call mattrib('s21', lname, 'ug/m3', 1)
c         call mattrib('s21', lname, 'pptv', 1)
         lname='sulfate level 1'
         call mattrib('s41', lname, 'ug/m3', 1)

         do k=1,nl
           vname='d'//lev(k)
           lname='DMS at sigma level '//lev(k)
           call mattrib(vname, lname, 'pptv', icomp)
         enddo

         do k=1,nl
           vname='o'//lev(k)
           lname='Sulfur dioxide at sigma level '//lev(k)
           call mattrib(vname, lname, 'pptv', icomp)
         enddo

         do k=1,nl
           vname='s'//lev(k)
           lname='Sulfate at sigma level '//lev(k)
           call mattrib(vname, lname, 'pptv', icomp)
         enddo

      ENDIF

c****************************************************************
c**** Some extra data (not "sflg" controlled)
c****************************************************************

      if (cdmap) then
c        lname='Drag - non leads'
c        call mattrib('cdm', lname, 'N/m2', 1)
c        lname='Leads drag'
c        call mattrib('pcd', lname, 'N/m2', 1)
         lname='Average drag'
         call mattrib('acd', lname, 'N/m2', 1)
      endif

c see prtcd.f
c        lname='Cd * 1000'
c        call mattrib('cdm', lname, ' ', 1)

c****************************************************************

      if(nncdf.gt.numcdf)then
         print*, 'Number of NetCDF/ASCII data files to be generated'
         print*, 'is ',nncdf
         print*, 'Increase parameter numcdf to at least ',nncdf
         print*, 'in subr ncinit + mattrib and 3 times in collst'
         stop
      end if

      return
      end
c
      subroutine mattrib(vname, lname, units, ind2)
c
c     Sets up an array to show compression per data file.
c     Checks if monthly mean netcdf file already exists and creates it 
c     if necessary
c
      implicit none
C Global parameters
      include 'netcdf.inc'
      include 'PARAMS.f'
      integer ndims
      parameter(ndims=4)
      real radtodeg
      parameter(radtodeg=57.29577951)
      integer numcdf
      parameter (numcdf=500) ! Must be >= total number of NetCDF files

C Argument list
      character vname*3
      character lname*(*)
      character units*(*)
      integer ind2

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
      include 'GAUSL.f'
      include 'TIMEX.f'
      include 'VERSION.f' 
      include 'FILES.f'

      integer ncdfcom(numcdf),nncdf
      character*3 ncdfchr(numcdf)
      common/ncdfiles/nncdf,ncdfcom,ncdfchr

C Local work arrays and variables
      real*4 ylat(lat2)
      real*4 xlon(lon)
      real*4 glat9(32),glon9(64) ! Note : for T63 compression only
      integer   dims(ndims)
      integer   vid
      integer   lngstr
      integer   ierr, mhistid
      integer   londim, latdim, mthdim, yrdim
      integer   latid, lonid, mthid, yrid
      logical used
      character filename*14
      character hist*82

      integer ind
      integer i
      integer j
      integer lat9
      integer lon9

C Local data, functions etc
      integer mth(12)
      data mth / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /

C Start code : ----------------------------------------------------------

c Set ind equal to value of argument ind2
      ind = ind2

c---- ind controls whether the full field is dumped (ind.eq.1), or
c----     if a compressed (area averaged) form is created (ind.eq.9).
c---- (ind=9 for T63 only). If T63, then data from T63 in arrays (192,96)
c---- may be compressed (averaged) to arrays (64,32). The use of "64" and
c---- "32" in glat9 and glon9 definitions above 
c---- is NOT a hangover from the R21 model!!
c----

      if(ind.ne.9)ind=1
c     Check if T63 is in use
      if (ind.eq.9) then 
        if (mw.ne.64) then
         print *,'ncinit/mattrib'
         print *,'wrong resolution for 9:1 compression arrays'
         stop
        end if
        lat9=lat2/3
        lon9=lon/3
      end if
c----------------------------------

c     Save index number of NetCDF file, and its character*3 name
      nncdf=nncdf+1
      if(nncdf.le.numcdf)then
        ncdfchr(nncdf)=vname ! NetCDF file name and
        ncdfcom(nncdf)=ind   ! compression factor (1 or 9)
      endif

c     Create NetCDF filename and check if file exists

      filename='s'//vname//'_'//str(9:13)//'.nc' !Harvey's suggestion      
      inquire(file=filename, exist=used)

c     Do nothing if the file already exists
      if (used) return

c     Set up the netcdf file
      ierr =  nf_create(filename, nf_clobber, mhistid)
      if (ierr .ne. nf_noerr) then
         print*, 'Error in creating file', ierr
         stop
      end if

c     Create dimensions: lon, lat, mth, yr
      IF(ind.eq.1)THEN
        ierr = nf_def_dim(mhistid, 'longitude', lon, londim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
        ierr = nf_def_dim(mhistid, 'latitude', lat2, latdim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ELSE
        ierr = nf_def_dim(mhistid, 'longitude', lon9, londim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
        ierr = nf_def_dim(mhistid, 'latitude', lat9, latdim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ENDIF
      ierr = nf_def_dim(mhistid, 'month', 12, mthdim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ierr = nf_def_dim(mhistid, 'year', nf_unlimited, yrdim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      dims(1)=londim
      dims(2)=latdim
      dims(3)=mthdim
      dims(4)=yrdim

c     Model Version
      ierr = nf_put_att_text(mhistid, nf_global, 'version',
     &                       len(trim(version)), trim(version))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     Experiment description
      ierr = nf_put_att_text(mhistid, nf_global, 'runtype', 5,
     &                       runtype)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     History
      call hstring(hist)
      ierr = nf_put_att_text(mhistid, nf_global, 'history',
     &                       len(trim(hist)), trim(hist))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     Define attributes for the dimensions
      ierr = nf_def_var(mhistid, 'latitude', nf_float, 1, dims(2),
     &                  latid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ierr = nf_put_att_text(mhistid, latid, 'units', 13,
     &                       'degrees_north')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      ierr = nf_def_var(mhistid, 'longitude', nf_float, 1, dims(1),
     &                  lonid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ierr = nf_put_att_text(mhistid, lonid, 'units', 12,
     &                       'degrees_east')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      ierr = nf_def_var(mhistid, 'month', nf_int, 1, dims(3), mthid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ierr = nf_put_att_text(mhistid, mthid, 'units', 6, 'months')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      ierr = nf_def_var(mhistid, 'year', nf_int, 1, dims(4), yrid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      ierr = nf_put_att_text(mhistid, yrid, 'units', 5, 'years')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     Define variable
      ierr = nf_def_var(mhistid, vname, nf_float, ndims, dims, vid) 
      if (ierr .ne. nf_noerr) then
         print*, ' Error in variable declaration ', vname
         stop
      end if

c     Define attributes for variable
      ierr = nf_put_att_text(mhistid, vid, 'long_name', lngstr(lname),
     &                       lname)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      if ( lngstr(units).ne.0 ) then
         ierr = nf_put_att_text(mhistid, vid, 'units', lngstr(units),
     &                          units)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"
      end if        

c     Leave define mode
      ierr = nf_enddef(mhistid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     Write the current year
      ierr = nf_put_var1_int(mhistid, yrid, 1, iyear)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Write the months
      ierr = nf_put_var_int(mhistid, mthid, mth)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c---- Files with or without data compression

c     Model latitudes
      do j=1,lat
         ylat(j) = real((-radtodeg * acos(sia(j))), 4)
         ylat(lat2+1-j) = -ylat(j)
      end do
c     Model longitudes
      do i=1,lon
         xlon(i) = real((float(i-1) * 360./float(lon)), 4)
      end do

      if((vname.eq.'icu').or.(vname.eq.'icv')) then
         do i=1,plon
            xlon(i) = real(blonu(i), 4)
         enddo
         do j=1,plat
            ylat(j) = real(blatu(j), 4)
         enddo
      endif

      IF(ind.eq.1)THEN

c---- No compression
c     Write the latitudes
      ierr = nf_put_var_real(mhistid, latid, ylat)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     Write the longitudes 
      ierr = nf_put_var_real(mhistid, lonid, xlon)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      ELSE

c---- 9:1 compression for T63
c     Write the latitudes
      do j=1,lat9
        glat9(j)=ylat(3*j-1)
      enddo
      ierr = nf_put_var_real(mhistid, latid, glat9)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

c     Write the longitudes 
      do i=1,lon9
        glon9(i)=xlon(3*i-2)
      enddo
      ierr = nf_put_var_real(mhistid, lonid, glon9)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      ENDIF

c     Close the netCDF file
      ierr = nf_close(mhistid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in mattrib"

      return
      end

      subroutine tmattrib(vname, lname, units, ind2)
c
c     Sets up an array to show compression per data file.
c     Checks if monthly mean netcdf file already exists and creates it 
c     if necessary
c
c     For tracer data files starting with "t" (not "s")
c
c     tt01 - tt??; tf01 - tf??; tv01 - tv?? (nl levels of month average data)
c     tT01 - tT??; tF01 - tF??; tV01 - tV?? (nl levels of end of month data)
c
      implicit none
C Global parameters
      include 'netcdf.inc'
      include 'PARAMS.f'
      integer ndims
      parameter(ndims=4)
      real radtodeg
      parameter(radtodeg=57.29577951)
      integer numcdf
      parameter (numcdf=500) ! Must be >= total number of NetCDF files

C Argument list
      character vname*3
      character lname*(*)
      character units*(*)
      integer ind2

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
      include 'GAUSL.f'
      include 'TIMEX.f'
      include 'VERSION.f' 
      include 'FILES.f'

      integer ncdfcom(numcdf),nncdf
      character*3 ncdfchr(numcdf)
      common/ncdfiles/nncdf,ncdfcom,ncdfchr

C Local work arrays and variables
      real*4 ylat(lat2)
      real*4 xlon(lon)
      real*4 glat9(32),glon9(64) ! Note : for T63 compression only
      integer   dims(ndims)
      integer   vid
      integer   lngstr
      integer   ierr, mhistid
      integer   londim, latdim, mthdim, yrdim
      integer   latid, lonid, mthid, yrid
      logical used
      character filename*14
      character hist*82

      integer ind
      integer i
      integer j
      integer lat9
      integer lon9

C Local data, functions etc
      integer mth(12)
      data mth / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /

C Start code : ----------------------------------------------------------

c Set ind equal to value of argument ind2
      ind = ind2

c---- ind controls whether the full field is dumped (ind.eq.1), or
c----     if a compressed (area averaged) form is created (ind.eq.9).
c---- (ind=9 for T63 only). If T63, then data from T63 in arrays (192,96)
c---- may be compressed (averaged) to arrays (64,32). The use of "64" and
c---- "32" in glat9 and glon9 definitions above 
c---- is NOT a hangover from the R21 model!!
c----

      if(ind.ne.9)ind=1
c     Check if T63 is in use
      if (ind.eq.9) then 
        if (mw.ne.64) then
         print *,'ncinit/mattrib'
         print *,'wrong resolution for 9:1 compression arrays'
         stop
        end if
        lat9=lat2/3
        lon9=lon/3
      end if
c----------------------------------

c     Save index number of NetCDF file, and its character*3 name
      nncdf=nncdf+1
      if(nncdf.le.numcdf)then
        ncdfchr(nncdf)=vname ! NetCDF file name and
        ncdfcom(nncdf)=ind   ! compression factor (1 or 9)
      endif

c     Create NetCDF filename and check if file exists

      filename='t'//vname//'_'//str(9:13)//'.nc' !Harvey's suggestion      
      inquire(file=filename, exist=used)

c     Do nothing if the file already exists
      if (used) return

c     Set up the netcdf file
      ierr =  nf_create(filename, nf_clobber, mhistid)
      if (ierr .ne. nf_noerr) then
         print*, 'Error in creating file', ierr
         stop
      end if

c     Create dimensions: lon, lat, mth, yr
      IF(ind.eq.1)THEN
        ierr = nf_def_dim(mhistid, 'longitude', lon, londim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
        ierr = nf_def_dim(mhistid, 'latitude', lat2, latdim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ELSE
        ierr = nf_def_dim(mhistid, 'longitude', lon9, londim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
        ierr = nf_def_dim(mhistid, 'latitude', lat9, latdim)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ENDIF
      ierr = nf_def_dim(mhistid, 'month', 12, mthdim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ierr = nf_def_dim(mhistid, 'year', nf_unlimited, yrdim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      dims(1)=londim
      dims(2)=latdim
      dims(3)=mthdim
      dims(4)=yrdim

c     Model Version
      ierr = nf_put_att_text(mhistid, nf_global, 'version',
     &                       len(trim(version)), trim(version))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Experiment description
      ierr = nf_put_att_text(mhistid, nf_global, 'runtype', 5,
     &                       runtype)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     History
      call hstring(hist)
      ierr = nf_put_att_text(mhistid, nf_global, 'history',
     &                       len(trim(hist)), trim(hist))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Define attributes for the dimensions
      ierr = nf_def_var(mhistid, 'latitude', nf_float, 1, dims(2),
     &                  latid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ierr = nf_put_att_text(mhistid, latid, 'units', 13,
     &                       'degrees_north')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      ierr = nf_def_var(mhistid, 'longitude', nf_float, 1, dims(1),
     &                  lonid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ierr = nf_put_att_text(mhistid, lonid, 'units', 12,
     &                       'degrees_east')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      ierr = nf_def_var(mhistid, 'month', nf_int, 1, dims(3), mthid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ierr = nf_put_att_text(mhistid, mthid, 'units', 6, 'months')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      ierr = nf_def_var(mhistid, 'year', nf_int, 1, dims(4), yrid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      ierr = nf_put_att_text(mhistid, yrid, 'units', 5, 'years')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Define variable
      ierr = nf_def_var(mhistid, vname, nf_float, ndims, dims, vid) 
      if (ierr .ne. nf_noerr) then
         print*, ' Error in variable declaration ', vname
         stop
      end if

c     Define attributes for variable
      ierr = nf_put_att_text(mhistid, vid, 'long_name', lngstr(lname),
     &                       lname)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      if ( lngstr(units).ne.0 ) then
         ierr = nf_put_att_text(mhistid, vid, 'units', lngstr(units),
     &                          units)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"
      end if        

c     Leave define mode
      ierr = nf_enddef(mhistid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Write the current year
      ierr = nf_put_var1_int(mhistid, yrid, 1, iyear)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Write the months
      ierr = nf_put_var_int(mhistid, mthid, mth)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c---- Files with or without data compression

c     Model latitudes
      do j=1,lat
         ylat(j) = real((-radtodeg * acos(sia(j))), 4)
         ylat(lat2+1-j) = -ylat(j)
      end do
c     Model longitudes
      do i=1,lon
         xlon(i) = real((float(i-1) * 360./float(lon)), 4)
      end do

      if((vname.eq.'icu').or.(vname.eq.'icv')) then
         do i=1,plon
            xlon(i) = real(blonu(i), 4)
         enddo
         do j=1,plat
            ylat(j) = real(blatu(j), 4)
         enddo
      endif

      IF(ind.eq.1)THEN

c---- No compression
c     Write the latitudes
      ierr = nf_put_var_real(mhistid, latid, ylat)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Write the longitudes 
      ierr = nf_put_var_real(mhistid, lonid, xlon)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      ELSE

c---- 9:1 compression for T63
c     Write the latitudes
      do j=1,lat9
        glat9(j)=ylat(3*j-1)
      enddo
      ierr = nf_put_var_real(mhistid, latid, glat9)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

c     Write the longitudes 
      do i=1,lon9
        glon9(i)=xlon(3*i-2)
      enddo
      ierr = nf_put_var_real(mhistid, lonid, glon9)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      ENDIF

c     Close the netCDF file
      ierr = nf_close(mhistid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in tmattrib"

      return
      end
