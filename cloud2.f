c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Moved REFAC1 and REFAC2 to /CLOUDPAR/, so that they can be read from
c namelist input. The default values are set in READNML1 to 0.85 and 0.95 
c respectively.
c SJP 2004/09/24
c
c Values of refac1 and refac2 reduced from 0.85/0.95 to 0.665/0.885, as these
c values have been diagnosed, in conjunction with an increase in the value of
c rcrits in newcloud.f, to give surface energy balance for pre-industrial
c (~1850) conditions. Note that, if it is desirable to achieve surface energy
c balance for present-day (~1992) conditions instead, then refac1/refac2 should
c be set to 0.7/0.9.
c SJP 2003/05/24
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: cloud2.f,v $
c Revision 1.40  2001/11/07 06:56:28  rot032
c Some further minor tuning from HBG
c
c Revision 1.39  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.38  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.36  2001/02/28 05:13:28  rot032
c Declare diffk.
c
c Revision 1.37.1.1  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.37  2001/06/04 02:27:02  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.35  2001/02/22 05:34:44  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.34  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.33  1999/06/30 05:29:37  rot032
c Final tuing of refac from HBG.
c
c Revision 1.32  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.31  1998/12/10  00:55:53  ldr
c HBG changes to V5-1-21
c
c Revision 1.30  1998/05/26  05:42:26  ldr
c Qcloud changes from LDR: new icefall Kessler option, Platt's (1996)
c optical properties and emissivity "fix" for low clouds.
c
c Revision 1.29  1997/12/17  23:22:56  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.28  1997/10/03  05:45:49  ldr
c Changes for sulfates from LDR
c
c Revision 1.27  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.26  1997/06/11  02:21:31  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.25  1997/02/21  00:26:10  ldr
c Go back to rcrits/l=0.85/0.75 for 18L version, tidy cloud2.f and make
c clddia.f, cloud.f general for NL.
c
c Revision 1.24  1996/11/22  03:42:44  ldr
c Comments and tidy-ups from LDR.
c
c******************************************************************************
c
c This is the interface between the Fels-Schwarzkopf radiation scheme and
c LDR's prognostic cloud water scheme. It is called by radfs.
c
c INPUT/OUTPUT:
c Input:
c
c See also include files PARAMS.f, PHYSPARAMS.f, CPARAMS.f, RDPARM.f, HCON.f,
c which contain parameters and physical constants.
c
c from common/fewflags in FEWFLAGS.f
c      debug - namelist flag to control single column debugging
c      lgdebug - latitude index for single column debugging
c      insdebug - hemisphere index for single column debugging
c      mgdebug  - longitude index for single column debugging
c
c from common/hybrpr in HYBRPR.f
c      dprf - pressure thickness at each sigma level
c      prf  - pressure at full levels
c
c from common/levdata in this subroutine
c      nlow - top level for low cloud
c      nmid - top level for mid cloud
c
c from common/lsmi in LSMI.f
c      imsl - land-sea mask ( = 4 for land points )
c
c from common/radisw in RADISW.f
c      coszro - zenith angle at grid pt. (used in SW scheme)
c
c from arguments
c      cldoff - flag set to T for clear sky radiation calculations
c      lg     - latitude index
c      ttg    - temperature
c      qlg    - cloud liquid water mixing ratio (kg/kg)
c      qfg    - cloud ice mixing ratio (kg/kg)
c      cfrac  - total cloud fraction (stratiform + convective)
c      clcon  - convective cloud fraction
c      qccon  - cloud water mixing ratio of convective clouds (kg/kg)
c
c Output:
c
c in common/radisw in RADISW.f (passed back to radiation scheme)
c      camt   - cloud amounts (locations specified by ktop/kbtm indices) 
c      cirab  - absorptivity of clouds in the near IR band (used in SW scheme)
c      cirrf  - reflectivity of clouds in the near IR band (used in SW scheme)
c      cuvrf - reflectivity of clouds in the visible band (used in SW scheme)
c      emcld - cloud emissivity (used in LW scheme)
c      kbtm  - index of (data level) pressure of cloud bottom (used in LW)
c      kbtmsw- index of (flux level) pressure of cloud bottom (used in SW)
c      ktop  - index of (data level) pressure of cloud top (used in LW)
c      ktopsw- index of (flux level) pressure of cloud top (used in SW)
c      nclds - no. clouds at each grid point
c
c in arguments
c      clat - cloud amount diagnostic (upside-down version of cfrac array)
c      clh - high level cloud diagnostic
c      cll - low level cloud diagnostic
c      clm - mid level cloud diagnostic
c
c******************************************************************************
 
      subroutine cloud2(cldoff,lg,ttg,qlg,qfg,cfrac,clcon,qccon,
     &                  cdso4,                                 !Inputs
     &                  clat,cll,clm,clh,Refflm,cldliq)        !Outputs

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /RADISW/ )

C Global parameters
      include 'PARAMS.f'     !Input model grid dimensions
      include 'PHYSPARAMS.f' !Input physical constants
      include 'CPARAMS.f'    !Input cloud scheme parameters
      include 'RDPARM.f'     !Input radiation scheme parameters
      include 'HCON.f'       !Input radiation physical constants

C Argument list
      logical cldoff
      integer lg
      real ttg(ln2,l)
      real qlg(ln2,l)
      real qfg(ln2,l)
      real cfrac(ln2,l)
      real clcon(ln2,l)
      real qccon(ln2,l)
      real cdso4(ln2,nl)
      real clat(ln2,l)
      real cll(ln2)
      real clm(ln2)
      real clh(ln2)
      real Refflm(ln2)
      real cldliq(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'CLOUDPAR.f'
      include 'HYBRPR.f'
      include 'RADISW.f'     !Output various things (see above)

C Global data blocks
      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.
      include 'LSMI.f'       !Input imsl (land sea mask)

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

C Local work arrays and variables
      real taul(ln2,l), taui(ln2,l) !Visible optical depths
      real Reffl(ln2,l), Reffi(ln2,l) !Effective radii (um)
      real Emw(ln2,l), Abw(ln2,l), Rew1(ln2,l), Rew2(ln2,l) !Water clouds
      real Emi(ln2,l), Abi(ln2,l), Rei1(ln2,l), Rei2(ln2,l) !Ice clouds
      real qlptot(ln2),taultot(ln2)
      real fice(ln2,l),tau_sfac(ln2,l)
      real rk(ln2),Cdrop(ln2,nl)

      integer i
      integer k
      integer mg
      integer nc
      integer ns

      real ab
      real cfl
      real cldht
      real deltai
      real diffk
      real dz
      real em
      real fcon
      real qlpath
      real refac
CSJP      real refac1
CSJP      real refac2
      real re1
      real re2
      real rhoa
      real sigmai
      real tciwc
      real tclwc
      real temp_correction
      real tmid
      real trani
      real tranw
      real wice
      real wliq

C Local data, functions etc

C Start code : ----------------------------------------------------------

c***INITIALIZE THE CLOUD AND CLOUD INDEX FIELDS
c     Except for the ground layer (nc=1) the assumption is that
c     no cloud exists. also, for this purpose, the cloud index is set
c     at one (p=0)
c     Don't set cirrf, cuvrf, cirab at the surface because these are set
c     the albedo by radfs.
      do 130 i=1,ln2
         camt(i,1)=zero
         emcld(i,1)=one
         ktop(i,1)=1
         kbtm(i,1)=1
         ktopsw(i,1)=1
         kbtmsw(i,1)=1
 130  continue
      do 140 k=2,lp1
         do 140 i=1,ln2
         camt(i,k)=zero
         emcld(i,k)=one
         cirrf(i,k)=0.
         cuvrf(i,k)=0.
         cirab(i,k)=0.
         ktop(i,k)=1
         kbtm(i,k)=1
         ktopsw(i,k)=1
         kbtmsw(i,k)=1
 140  continue
c***NOW SET CLOUD AND CLOUD INDEX FIELDS DEPENDING ON THE NO. OF CLOUDS
      nc=1
      do 150 mg=1,ln2
c---FIRST, THE ground layer (nc=1)
         emcld(mg,nc)=one
         camt(mg,nc)=one
         ktop(mg,nc)=lp1
         kbtm(mg,nc)=lp1
         ktopsw(mg,nc)=lp1
         kbtmsw(mg,nc)=lp1
 150  continue

      If (cldoff) Then

         do 10 mg=1,ln2
            cll(mg)=0.
            clm(mg)=0.
            clh(mg)=0.
            nclds(mg)=0
 10      continue

      Else

        if(ukconv.and..not.coupled_aero)then
c--- No lowest level cloud for radiation with ukconv
          do mg=1,ln2
            cfrac(mg,1)=0.0
          enddo
        endif

        do k=1,nl
          do mg=1,ln2
            clat(mg,nlp-k)=cfrac(mg,k) !To pass back for plotting
          enddo
        enddo
        
        do mg=1,ln2
          nclds(mg)=0
          cll(mg)=0.
          clm(mg)=0.
          clh(mg)=0.
          cldliq(mg)=0.
          qlptot(mg)=0.
          taultot(mg)=0.
        enddo
        
        
c Diagnose low, middle and high clouds; nlow,nmid are set up in initax.f
        
        do k=1,nlow
          do mg=1,ln2
            cll(mg)=cll(mg)+cfrac(mg,k)-cll(mg)*cfrac(mg,k)
          enddo
        enddo
        do k=nlow+1,nmid
          do mg=1,ln2
            clm(mg)=clm(mg)+cfrac(mg,k)-clm(mg)*cfrac(mg,k)
          enddo
        enddo
        do k=nmid+1,nl-1
          do mg=1,ln2
            clh(mg)=clh(mg)+cfrac(mg,k)-clh(mg)*cfrac(mg,k)
          enddo
        enddo

c Set up rk and Cdrop
        do mg=1,ln2
          if(imsl(mg,lg).lt.4)then !sea
            rk(mg)=0.8
          else            !land
            rk(mg)=0.67
          endif
        enddo
        
        if(naerosol_i(1).gt.0)then
          do k=1,nl-1
            do mg=1,ln2
              Cdrop(mg,k)=cdso4(mg,k)
            enddo
          enddo
        else
          do mg=1,ln2
            if(imsl(mg,lg).lt.4)then !sea
              do k=1,nl-1
                Cdrop(mg,k)=Cdrops
              enddo
            else            !land
              do k=1,nl-1
                Cdrop(mg,k)=Cdropl
              enddo
            endif
          enddo
        endif

c Define the emissivity (Em), and the SW properties (Re, Ab) for liquid (w)
c and ice (i) clouds respectively.
        
        do k=1,nl-1
          do mg=1,ln2
            taul(mg,k)=0.
            taui(mg,k)=0.
            Reffl(mg,k)=0.
            Reffi(mg,k)=0.
            Emw(mg,k)=0.
            Emi(mg,k)=0.
            if(cfrac(mg,k).gt.0)then
C***              fcon=min(1.,qccon(mg,k)/(qlg(mg,k)+qfg(mg,k)))
C***              tau_sfac_c=0.5
C***              fmin=0.6 !tau_sfac_s at C=0.2
C***              fmax=0.9 !tau_sfac_s at C=0.8
C***              Cstrat=(cfrac(mg,k)-clcon(mg,k))/(1-clcon(mg,k))
C***              tau_sfac_s=fmin+(fmax-fmin)*(Cstrat-0.2)/0.6
C***              tau_sfac_s=max(min(fmax,tau_sfac_s),fmin)
C***                    !Scale fac for tau
C***              tau_sfac(mg,k)=fcon*tau_sfac_c + (1-fcon)*tau_sfac_s
              tau_sfac(mg,k)=1.
              fice(mg,k) = qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))
            endif         !cfrac
          enddo
        enddo
              
c Liquid water clouds
        do k=1,nl-1
          do mg=1,ln2
            if((cfrac(mg,k).gt.0).and.(qlg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              cfl=cfrac(mg,k)*(1-fice(mg,k))
              Wliq=rhoa*qlg(mg,k)/cfl     !kg/m^3
              cldliq(mg)=cldliq(mg)+cfl-cldliq(mg)*cfl !Liquid cloud cover
                
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                
              Reffl(mg,k)=
     &             (3*2*Wliq/(4*rhow*pi*rk(mg)*Cdrop(mg,k)))**(1./3)
c              Reffl(mg,k)=
c    &             (3*Wliq/(4*rhow*pi*rk(mg)*Cdrop(mg,k)))**(1./3)
              qlpath=Wliq*dz
              taul(mg,k)=tau_sfac(mg,k)*1.5*qlpath/(rhow*Reffl(mg,k))
              qlptot(mg)=qlptot(mg)+qlpath*cfl
              taultot(mg)=taultot(mg)+taul(mg,k)*cfl

c Water cloud emissivity according to Martin Platt

C***          deltvl=taul(mg,k)*1.26 !Mult by 2^(1/3) so using mid cloud Reff
C***          deltal=min(0.4*deltvl, 45.) !IR optical depth for liq.
C***          if(deltvl.gt.0.4)then
C***            diffk=1.6
C***          else
C***            diffk=1.8
C***          endif
C***          Emw(mg,k) = 1.0 - exp(-diffk*deltal) !em of strat water clouds

c Or, water-cloud emissivity following the Sunshine scheme

              cldht=dz/1000.       !in km
              tclwc=Wliq*1000.     !in g/m**3
              tranw = exp( -1.66 * cldht * 50.885 * tclwc ** 0.769917)
              Emw(mg,k) = 1.0 - tranw    ! em is (1 - transmittance)
            endif         !cfrac
          enddo
        enddo
              
c Ice clouds : Choose scheme according to resolution
       IF(lw.eq.22)THEN

        do k=1,nl-1
          do mg=1,ln2
            if((cfrac(mg,k).gt.0).and.(qfg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*fice(mg,k)) !kg/m**3
              Reffi(mg,k)=3.73e-4*Wice**0.216 !Lohmann et al. (1999)
              taui(mg,k)=1.5*Wice*dz/(rhoice*Reffi(mg,k))
              deltai=min(0.5*taui(mg,k), 45.) !IR optical depth for ice.
              taui(mg,k)=tau_sfac(mg,k)*taui(mg,k)

c Ice-cloud emissivity following Platt

              if(taui(mg,k).gt.0.4)then
                diffk=1.6
              else
                diffk=1.8
              endif
              Emi(mg,k) = 1.0 - exp(-diffk*deltai)
            endif         !cfrac
          enddo
        enddo

       ELSE

        do k=1,nl-1
          do mg=1,ln2
            if((cfrac(mg,k).gt.0).and.(qfg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*fice(mg,k)) !kg/m**3
              sigmai = aice*Wice**bice !visible ext. coeff. for ice
              taui(mg,k)=sigmai*dz !visible opt. depth for ice
              Reffi(mg,k)=1.5*Wice*dz/(rhoice*taui(mg,k))
              taui(mg,k)=tau_sfac(mg,k)*taui(mg,k)

c Ice-cloud emissivity following the Sunshine scheme

              tmid=ttg(mg,k)-273.15  !in celsius
              cldht=dz/1000.       !in km
              tciwc=Wice*1000.     !in g/m**3

              temp_correction=1.047E+00+tmid
     &        *(-9.13e-05+tmid
     &        *(2.026e-04-1.056e-05*tmid))
              temp_correction=max(1.0, temp_correction)
              trani = exp(-1.66*temp_correction * tciwc*cldht
     &              / (0.630689d-01+0.265874*tciwc))
c-- Limit ice cloud emessivities
              trani=min(0.70,trani)
              Emi(mg,k) = 1.0 - trani    ! em is (1 - transmittance)
            endif         !cfrac
          enddo
        enddo

       ENDIF

c Calculate the effective radius of liquid water clouds seen from above

        do mg=1,ln2
          if(taultot(mg).gt.0.)then
            Refflm(mg)=1.5*qlptot(mg)/(rhow*taultot(mg))
          else
            Refflm(mg)=0.
          endif
        enddo
          
c Calculate the SW cloud radiative properties for liquid water and ice clouds
c respectively, following Tony Slingo's (1989) Delta-Eddington scheme.
        
        call slingo (Reffl, taul, coszro, !inputs
     &       Rew1, Rew2, Abw )  !outputs
        
        call slingi (Reffi, taui, coszro, !inputs
     &       Rei1, Rei2, Abi )  !outputs

CSJP        if(lw.eq.22)then

C Values diagnosed to give surface energy balance for pre-industrial conditions
CSJP          refac1=0.655
CSJP          refac2=0.885

C Values diagnosed to give surface energy balance for present-day conditions
CSJP          refac1=0.7
CSJP          refac2=0.9

C Original values
CSJP          refac1=0.85
CSJP          refac2=0.95

CSJP        else
CSJP          refac1=0.90
CSJP          refac2=1.00
CSJP        endif

        do k=1,nl-1
          do mg=1,ln2
            if(cfrac(mg,k).gt.0.)then
              fcon=min(1.,qccon(mg,k)/(qlg(mg,k)+qfg(mg,k)))
c Original refac :
c             refac=0.7*fcon+0.9*(1-fcon)
c Mk3 with no direct aerosol effect :
c             refac=0.7*fcon+0.85*(1-fcon)
c Mk3 with direct aerosol effect :
              refac=refac1*fcon+refac2*(1-fcon)

              Rei1(mg,k)=min(refac*Rei1(mg,k),1.)
              Rei2(mg,k)=min(refac*Rei2(mg,k),1.-2*Abi(mg,k))
              Rew1(mg,k)=min(refac*Rew1(mg,k),1.)
              Rew2(mg,k)=min(refac*Rew2(mg,k),1.-2*Abw(mg,k))
            endif
          enddo
        enddo

c Weight cloud properties by liquid/ice fraction
        
        do k=1,nl-1
          do mg=1,ln2
            if(cfrac(mg,k).gt.0.)then
              Re1 = fice(mg,k)*Rei1(mg,k) + (1-fice(mg,k))*Rew1(mg,k)
              Re2 = fice(mg,k)*Rei2(mg,k) + (1-fice(mg,k))*Rew2(mg,k)
              Em = fice(mg,k)*Emi(mg,k) + (1-fice(mg,k))*Emw(mg,k)
              if(prf(mg,k).gt.800.) Em = 1.
              Ab = fice(mg,k)*Abi(mg,k) + (1-fice(mg,k))*Abw(mg,k)
              
              nclds(mg)=nclds(mg)+1
              nc=nclds(mg)+1
              camt(mg,nc)=cfrac(mg,k)
              ktop(mg,nc)=nlp-k
              kbtm(mg,nc)=nlp-k
              emcld(mg,nc)=Em
              ktopsw(mg,nc)=nlp-k
              kbtmsw(mg,nc)=nlp-k+1
              cuvrf(mg,nc)=Re1
              cirrf(mg,nc)=Re2
              cirab(mg,nc)=2*Ab
            endif
          enddo
        enddo

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            write(25,'(a,i1)')'After cloud2'
            write(25,9)'cdso4 ',(cdso4(mg,k),k=1,nl)
            write(25,91)'Rew ',((Rew1(mg,k)+Rew2(mg,k))/2,k=1,nl)
            write(25,91)'Rei ',((Rei1(mg,k)+Rei2(mg,k))/2,k=1,nl)
            write(25,91)'Abw ',(Abw(mg,k),k=1,nl)
            write(25,91)'Abi ',(Abi(mg,k),k=1,nl)
            write(25,91)'Emw ',(Emw(mg,k),k=1,nl)
            write(25,91)'Emi ',(Emi(mg,k),k=1,nl)
            write(25,9)'Reffl ',(Reffl(mg,k),k=1,nl)
            write(25,9)'Reffi ',(Reffi(mg,k),k=1,nl)
            write(25,*)
          endif
        endif
 9      format(a,30g10.3)
 91     format(a,30f10.3)
      End If   !cldoff
      

      return
      end

c******************************************************************************

c Calculate SW radiative properties for water clouds using delta-Eddington
c scheme as given by Slingo (1989) JAS 46, 1419-1427.
c Coefficients di and fi modified to use Reff in SI units.

      subroutine slingo(reff, tau, mu0,       !inputs
     &                  refl1, refl2, abso )  !outputs

      implicit none
C Global parameters
      include 'PARAMS.f'
      integer nbands
      parameter (nbands=4)

C Argument list
      real reff(ln2,nl)
      real tau(ln2,nl)
      real mu0(ln2)
      real refl1(ln2,nl)
      real refl2(ln2,nl)
      real abso(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      double precision exparg,denom,epsilon,omega,f,omwf

      integer i
      integer k
      integer mg

      real absband
      real alpha1
      real alpha2
      real alpha3
      real alpha4
      real beta
      real beta0
      real e
      real gam1
      real gam2
      real gi
      real rdif
      real rm
      real rdir
      real tdb
      real tdif
      real tdir
      real ttot
      real u1
      real u2

C Local data, functions etc
      real ci(nbands)
      data ci / -5.62e-8, -6.94e-6, 4.64e-4,  2.01e-1  /

      real di(nbands)
c      data di / 1.63e-7,  2.35e-5,  1.24e-3,  7.56e-3  /
      data di / 1.63e-1,  2.35e+1,  1.24e+3,  7.56e+3  / !Si units

      real ei(nbands)
      data ei / 0.829,    0.794,    0.754,    0.826    /

      real fi(nbands)
c      data fi / 2.482e-3, 4.226e-3, 6.560e-3, 4.353e-3 /
      data fi / 2.482e+3, 4.226e+3, 6.560e+3, 4.353e+3 / !SI units

      real wi(nbands)
      data wi / 0.459760, 0.326158, 0.180608, 0.033474 /
 
      real wi24
      data wi24 / 0.540240 / !Sum of wi, i=2 to 4

C Start code : ----------------------------------------------------------

      do k=1,nl
        do mg=1,ln2
          refl1(mg,k)=0.
          refl2(mg,k)=0.
          abso(mg,k)=0.
        enddo
      enddo

      do i=1,nbands
        do k=1,nl-1
          do mg=1,ln2
            if(tau(mg,k).gt.0..and.mu0(mg).gt.0.)then
               reff(mg,k)=min(20.e-6,max(4.e-6,reff(mg,k)))
               omega=1-(ci(i)+di(i)*reff(mg,k))     
               gi=ei(i)+fi(i)*reff(mg,k)
               beta0=(3./7.)*(1-gi)
               beta=0.5-0.75*mu0(mg)*gi/(1+gi)
               f=gi**2
               U1=7./4.
               U2=(7./4)*(1.-(1-omega)/(7*omega*beta0))
               alpha1=U1*(1.-omega*(1-beta0))
               alpha2=U2*omega*beta0
               alpha3=(1-f)*omega*beta
               alpha4=(1-f)*omega*(1-beta)
               epsilon=sqrt(alpha1**2-alpha2**2)
               rM=alpha2/(alpha1+epsilon)
               E=exp(-epsilon*tau(mg,k))
               omwf=1-omega*f
               denom=omwf**2-epsilon**2*mu0(mg)**2
               gam1=(omwf*alpha3-mu0(mg)*(alpha1*alpha3+alpha2*alpha4))
     &               /denom 
               gam2=(-omwf*alpha4-mu0(mg)*(alpha1*alpha4+alpha2*alpha3))
     &               /denom
               exparg=dmin1(70.0d0,omwf*tau(mg,k)/mu0(mg))
               Tdb=exp(-exparg)
               Rdif=rM*(1-E**2)/(1-(E*rM)**2)
               Tdif=E*(1-rM**2)/(1-(E*rM)**2)
               Rdir=-gam2*Rdif-gam1*Tdb*Tdif+gam1
               Tdir=-gam2*Tdif-gam1*Tdb*Rdif+gam2*Tdb
               Ttot=Tdb+Tdir
               Absband=1-Rdir-Ttot
               Absband=max(0., Absband) !Needed for 32 bit
c               if(absband.gt.1..or.absband.lt.0)then
c                 print*,'Warning slingi: band, abs =',i,absband
c               endif
               abso(mg,k)=abso(mg,k)+Absband*wi(i)
               if(i.eq.1)then
                 refl1(mg,k)=Rdir
               else
                 refl2(mg,k)=refl2(mg,k)+Rdir*wi(i)/wi24
               endif
             endif
           enddo
         enddo
       enddo

       return
       end

c******************************************************************************

c Slingo type scheme for ice cloud SW properties (similar to water clouds)
c Single scattering albedo follows Francis etal (1994) QJRMS 120, 809--848.
c Lambdas modified to use Reff in SI units.

      subroutine slingi(reff, tau, mu0,       !inputs
     &                  refl1, refl2, abso )  !outputs

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      integer nbands
      parameter (nbands=4)

C Argument list
      real reff(ln2,nl)
      real tau(ln2,nl)
      real mu0(ln2)
      real refl1(ln2,nl)
      real refl2(ln2,nl)
      real abso(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      double precision vm !For 32 bit
      double precision exparg,denom,epsilon,omega,f,omwf

      integer i
      integer k
      integer mg

      real absband
      real alpha1
      real alpha2
      real alpha3
      real alpha4
      real beta
      real beta0
      real e
      real gam1
      real gam2
      real gi
      real rdif
      real rm
      real rdir
      real tdb
      real tdif
      real tdir
      real ttot
      real u1
      real u2

C Local data, functions etc
      real lambda(nbands)
      data lambda / 0.470e-6, 0.940e-6, 1.785e-6, 3.190e-6 /

      real nprime(nbands)
      data nprime / 5.31e-9,  8.74e-7,  3.13e-4,  1.07e-1  /

      real wi(nbands)
      data wi     / 0.459760, 0.326158, 0.180608, 0.033474 /

      real wi24
      data wi24   / 0.540240 / !Sum of wi, i=2 to 4

C Start code : ----------------------------------------------------------

      do k=1,nl
        do mg=1,ln2
          refl1(mg,k)=0.
          refl2(mg,k)=0.
          abso(mg,k)=0.
        enddo
      enddo

      do i=1,nbands
        do k=1,nl-1
          do mg=1,ln2
            if(tau(mg,k).gt.0..and.mu0(mg).gt.0.)then
               vm=2*pi*reff(mg,k)*nprime(i)/lambda(i)
               omega=0.5+(vm+3)/(6*(vm+1)**3)
               gi=0.8  !Constant asymmetry parameter
               beta0=(3./7.)*(1-gi)
               beta=0.5-0.75*mu0(mg)*gi/(1+gi)
               f=gi**2
               U1=7./4.
               U2=(7./4)*(1.-(1-omega)/(7*omega*beta0))
               alpha1=U1*(1.-omega*(1-beta0))
               alpha2=U2*omega*beta0
               alpha3=(1-f)*omega*beta
               alpha4=(1-f)*omega*(1-beta)
               epsilon=sqrt(alpha1**2-alpha2**2)
               rM=alpha2/(alpha1+epsilon)
               E=exp(-epsilon*tau(mg,k))
               omwf=1-omega*f
               denom=omwf**2-epsilon**2*mu0(mg)**2
               gam1=(omwf*alpha3-mu0(mg)*(alpha1*alpha3+alpha2*alpha4))
     &               /denom
               gam2=(-omwf*alpha4-mu0(mg)*(alpha1*alpha4+alpha2*alpha3))
     &               /denom
               exparg=dmin1(70.0d0,omwf*tau(mg,k)/mu0(mg))
               Tdb=exp(-exparg)
               Rdif=rM*(1-E**2)/(1-(E*rM)**2)
               Tdif=E*(1-rM**2)/(1-(E*rM)**2)
               Rdir=-gam2*Rdif-gam1*Tdb*Tdif+gam1
               Tdir=-gam2*Tdif-gam1*Tdb*Rdif+gam2*Tdb
               Ttot=Tdb+Tdir
               Absband=1-Rdir-Ttot
               Absband=max(0., Absband) !Needed for 32 bit
c               if(absband.gt.1..or.absband.lt.0)then
c                 print*,'Warning slingi: band, abs =',i,absband
c               endif
               abso(mg,k)=abso(mg,k)+Absband*wi(i)
               if(i.eq.1)then
                 refl1(mg,k)=Rdir
               else
                 refl2(mg,k)=refl2(mg,k)+Rdir*wi(i)/wi24
               endif
             endif
           enddo
         enddo
       enddo

       return
       end
