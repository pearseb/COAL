#if defined dmr_cfc && !defined dmr_oline 
* **NO CODE TO COMPILE FOR CFC RUNS

#else


c  These routines are required to run the BGC module 
c  To compile, first apply preprocessing
c
c 25-11-2015 PJB
c   general aesthetic clean of the code
c
c 
c
c
c
c


c=====================================================================c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c.....................................................................c

      subroutine bio_i 

c subroutine "bio_i" does the following:
c
c 1. reads in the BGC model arrays/variables defined in bio.h
c             the ocean model arrays/variables contained in obgc.h
c 2. reads namelist variables defined in "input" file
c 3. defines a number of initial values 
c       - these may be overwriten by the values provided in
c         the "input" file
c 4. Sets air-sea flux algorithms for different tracers with "igas"
c       - igas == 1     Oxygen
c       - igas == 2     Carbon Dixoide
c       - igas == 3     F12
c       - igas == 4     F11
c       - igas == 5     Carbon 13
c       - igas == 6     Carbon 14 
c       - igas == 7     Nitrous Oxide
c
c
c
c
c......................................................................c

*   INITIALISE THE VARIABLES OF THE BIOGEOCHEMCIAL MODEL

      include "obgc.h"
      include "extra.h"
      include "bio.h"

      parameter (ntmaxp1 = nbio2d)          ! ntmaxp1 is...
      character *12 trc_name(ntmaxp1)        
      
      ! reads variables stated in the input file on model run initiation
      namelist /bgc_names/ s_npp,
     + ratio_pop, pwl_pop,
     + rain_ratio, ratio_pic, pwl_pic, terr_input,
     + pmodel, temp, time, xmixl, pcon, pwr_sat,
     + tr_off, 
     + igas,
     + n_pho, n_alk, n_dic, n_c14, n_si, n_oxy, n_c13, n_fe,
     + n_no3, n_n15, n_tou, n_age, n_tomz, n_n2o,
     + trc_name,
     + idump, jdump, kdump, xdump, pdump,
     + pwl_opal, ratio_opal,
     + geo_seq,
     + tglobal, ivredc, ivredo, ivredpic,
     + ipic, den, fix, atmdep, 
     + sedfluxes, Boh_a, subgrid, Sden_acl,
     + vary_stoich, biomol, solarcarb,
     + Den_o2lim, Den_nlim, Den_type,
     + NO3diagnostics, Fe_kd, p_kd, s_npp_fix,
     + MichaelisMenton, p_k, n_k, fe_k,
     + OptimalUptake, VoA,
     + ReminRelax, ReminTemp_Marsay, ReminTemp_Q10, Q10, T_rem, 
     + ReminPico, ReminDiagnostics,
     + pN2O, N2Otemp, N2Odepth,
     + omz_min, omz_mean, omz_max,
     + N15diagnostics, latj, loni,
     + Eass, Erem, Efix, Eatm, Ewc, Esed,
     + vary_frac13C, diff_out
     
      ! Set default initial values for BGC model prior to reading input file
      data pwl_pop /-.836/
      data pwl_pic /-2./
      data pwl_opal /12e3/  ! 12 km
      data pwr_sat /0./
      data temp,time,xmixl,trad,pcon /24.,76.,50.,0.,0./
      data geo_seq /ntmax*0/
      data ivredc /0/
      data ivredo /0/, ipic /0/, ivredpic /0/
      
      ! Define the lowest possible value for nutrients in the model 
      data c0 /0./
      
      ! Nfixation is off before reading namelist
      data fixswitch /0./
      
      ! Tracer names
      data trc_name /ntmaxp1*' '/
      print*, trc_name,ntmaxp1
      data ratio_opal /ntmax*0./, 
     +     ratio_pop /ntmax*0./, ratio_pic /ntmax*0./
      data tglobal /ntmax*1./
      
      ! Is silicic acid simulated by the BGC model?
      if (n_si .ne.0 ) then
          ratio_opal(n_si)=10.
      endif

      ! Read the input file
      read  (5, bgc_names)
      write (6, bgc_names)
      
      do n=1,nt
         if (trc_name(n) .ne. ' ') then
             trname(n)=trc_name(n)
         endif
      enddo

      trname(1) = "Temperature"
      trname(2) = "Salinity"
      if (n_oxy.ne.0) then
          trname(n_oxy)="Oxygen"
          igas(n_oxy-2)=1
      endif
      if (n_pho.ne.0) then
          trname(n_pho)="Phosphate"
          igas(n_pho-2)=0
      endif
      if (n_dic.ne.0) then
          trname(n_dic)="DIC"
          igas(n_dic-2)=2
      endif
      if (n_alk.ne.0) then
          trname(n_alk)="Alk"
          igas(n_alk-2)=0
      endif
      if (n_c13.ne.0) then
          trname(n_c13)="C13"
          igas(n_c13-2)=2
      endif
      if (n_c14.ne.0) then 
          trname(n_c14)="C14"
          igas(n_c14-2) =2 
      endif
      if (n_si.ne.0) then
          trname(n_si)="Silicate"
          igas(n_si-2)=0
      endif
      if (n_fe.ne.0) then
          trname(n_fe)="Iron"
          igas(n_fe-2)=0
      endif
      if (n_no3.ne.0) then
          trname(n_no3)="Nitrate"
          igas(n_no3-2)=0
      endif
      if (n_n15.ne.0) then
          trname(n_n15)="Nitrogen15"
          igas(n_n15-2)=0
      endif
      if (n_tou.ne.0) then
          trname(n_tou)="True OU"
          igas(n_tou-2)=1
      endif
      if (n_age.ne.0) then
          trname(n_age)="Ideal Age"
          igas(n_age-2)=0
      endif
      if (n_tomz.ne.0) then
          trname(n_tomz)="Age of OMZ"
          igas(n_tomz-2)=0
      endif
         
      if (n_n2o.ne.0) then
          trname(n_n2o)="NitrousOxide"
          igas(n_n2o-2)=7
      endif

      print*, "OBGC: Tracer names = ",trname
      print*, "OBGC: gas indexes = ",igas

        ntr_req = max(n_oxy,n_dic,n_c13,n_c14,n_alk,n_pho,n_si,
     +                n_fe,n_no3,n_n15,n_tou,n_age,n_tomz,n_n2o)
      if (nt.ne. ntr_req) then 
          print *,"OBGC: nt does not equal number of requested tracers"
          print*,"nt = ",nt,"  ntr_req = ",ntr_req
          stop
      endif


#ifdef rjm_co2_dumping
      write(6,*) 'CO2 dumping site',xt(idump),yt(jdump),zt(kdump),
     +           ' Amount (GT/yr)',xdump
      xdump = xdump*1e9*1e15/3.73248e8    ! convert to pmol/s
#endif

#ifdef rjm_geo_seq
      gCtppm = 0.48
      write(6,*) 'Prescribe geological sequestration',geo_seq
      do l=1,nt
         geo_seq(l) =geo_seq(l)/(365.*24.*3600)  ! GtC/yr to GtC/sec
         ttflux(l)=0.
      enddo
#endif

      write(6,*) 'Phosphate addition in the southern subtropics',pdump


*   INITIALISE SEDIMENTS AND AIR-SEA FLUXES TO ZERO
      do j=1,jmt
         do i=1,imt
            sediments(i,j)=c0
            do n=1,nt-2
               fluxgas(i,j,n)=c0
            enddo
        enddo
      enddo

*   DISPLAY FRACTIONATION FACTORS FOR NITROGEN ISOTOPES
      if (n_n15.ne.0) then 
         print*, " "
         print*, " "
         print*, "NITROGEN ISOTOPE FRACTIONATION FACTORS "
         print*, " "
         print*, " ASSIMILATION = ", Eass, "per mil"
         print*, " REMINERALISATION = ", Erem, "per mil"
         print*, " N2 FIXATION = ", Efix, "per mil"
         print*, " AEOLIAN DEPOSITION = ", Eatm, "per mil"
         print*, " PELAGIC DENITRIFICATION = ", Ewc, "per mil"
         print*, " SEDIMENT DENITRIFICATION = ", Esed, "per mil"
         print*, " "
         print*, " "
      endif

c SBC for the BGC module come from atmosphere
*   READ SEASONAL FORCINGS AT SURFACE OCEAN FOR BGC MODEL
*       - windspeed squared
*       - shortwave radiation
*       - surface pressure
*       - sea ice area fraction
*       - aeolian dust flux
*       - aeolian Nr flux (reactive nitrogen)
*   READ subgrid scale bathymetry file
* need windspeed, shortwave radiation and surface pressure
      call nread_ocmip
      print*,' '
      print*,' '


      if (.not.subgrid) then
         print*, " "
         print*, " Erasing subgrid scale bathymetry from calculations "
         print*, " "
         do j=1,jmt
         do i=1,imt
         do k=1,km
           sgb_frac(i,j,k) = 0.0
         enddo
         enddo
         enddo
      endif

*   INITIALISE PARTICULATE ORGANIC PHOSPHATE REMINERALISATION FUNCTIONS - for pop and pic

      ! constants
      zpop=1e-4     ! 1/10000cm or 1/100m
      zdia=4e-4     ! 1/2500 cm or 1/25 m
      zpic=1e+2     ! m to cm
      i_pop=4 !2   !not used
      i_pic=12 !6
      i_opal=12     ! minimum depth to start remineralization
      
      do j=1,jmt
         do i=1,imt
      
            fmin_pop(i,j,1)=1.   ! fraction of organic matter conserved at surface
            fmin_pop(i,j,km)=0.  ! fraction of organic matter conserved at ocean bottom
            fmin_dia(i,j,1)=1.   ! fraction of organic matter conserved at surface
            fmin_dia(i,j,km)=0.  ! fraction of organic matter conserved at ocean bottom
            fmin_pic(i,j,1)=1.   ! fraction of inorganic carbon conserved at surface
            fmin_pic(i,j,km)=0.  ! fraction of inorganic carbon conserved at ocean bottom
            fmin_opal(i,j,1)=1.
            fmin_opal(i,j,km)=0.

            sedtot(i,j) = 0.0
            
            do k=1,km
               
               dizt(k) = 1./dzt(k) 
               dz_t=dzt(k)*.5
               sedtot(i,j) = sedtot(i,j) + sgb_frac(i,j,k)
               
               fmin_pop(i,j,k) = min(1.0, ((zt(k)+dz_t)*zpop)**pwl_pop)
     +                           *(1.0-sedtot(i,j))
               fmin_dia(i,j,k) = min(1.0, ((zt(k)+dz_t)*zdia)**pwl_pop)
     +                           *(1.0-sedtot(i,j))
               fmin_pic(i,j,k) = min(1.0, exp( -(zt(k)+dz_t) /
     +                           (zpic*pwl_pic)) ) *(1.0-sedtot(i,j))
               fmin_opal(i,j,k)= min(1.0, exp( -(zt(k)+dz_t) /
     +                           (zpic*pwl_opal)) ) *(1.0-sedtot(i,j))
            enddo
            fmin_pic(i,j,1) = 1.0
            if (.not.sedfluxes) then
               fmin_pop(i,j,km) = 0.0
               fmin_dia(i,j,km) = 0.0
               fmin_pic(i,j,km) = 0.0
               fmin_opal(i,j,km) = 0.0
            endif

         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c




c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine obgc_atmco2 (atmco2)


c get the co2 from rrvco2 (RADISW.f) for coupled run
c requires this interface to get access to rrvco2 without messing up 
c ocean code
c......................................................................c

      include "PARAMS.f"
      include "RDPARM.f"
      include "RADISW.f"    ! holds the atmospheric co2 variable

      atmco2 = rrvco2 *1e6   !ppm units
 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c======================================================================c










c======================================================================c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c......................................................................c

      subroutine bio_geo_chem (j, n, EPmax)

* The ocean bgc module
c
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c    primary author:  r.j.matear  
c                     richard.matear@ml.csiro.au
c
c    secondary author: p.j.buchanan
c                      pearse.buchanan@utas.edu.au
c
c            -  Nitrogen cycle (no3, d15n)
c            -  Sedimentary chemistry
c      
c 	The matrix source will hold the source and sink terms derived 
c   from the the biogeochemical model
c
c   The bio_geo_chem subroutine does the following:
c
c       0. Calculate carbon chemistry terms
c       1. Biological uptake 
c       2. Biological remineralization
c       3. Check Oxygen and sediments 
c       4. Air-sea exchange
c
c   Setup for the following tracer order
c   O2, PO4, CO2
c
c......................................................................c

      use msingledouble
      use mvars
      include "obgc.h"
      include "extra.h"
      include "bio.h"
      include "OCEAN_NML.f"
      include "SCALAR.f"
      include "bulkf_variables.f"
      
      !	 include "co2.h"
      include "FEWFLAGS.f"  ! logical flags of which lcouple is used here

      REAL(kind=rx1), DIMENSION(1) :: CCco2flux, CCco2ex, CCdpCO2
      REAL(kind=rx1), DIMENSION(1) :: CCph,CCpco2,CCfco2,CCco2,CChco3,
     +   CCco3, CCOmegaA, CCOmegaC, CCBetaD, CCrhoSW, CCp, CCtempis
      REAL(kind=rx1), DIMENSION(1) :: CCtem, CCsal, CCalk, CCdic, CCsil,
     +   CCpo4, CCpatm, CCdep, CClat, CClon, kw660, xco2
      character(10) :: optCON = 'mol/m3'
      character(10) :: optT = 'Tpot'
      character(10) :: optP = 'm'

      parameter (istrt=2, iend=imt-1)
      save pbtime,itimes,sum1,sum2,sum3,sum4

      data pbtime,itimes /0.,-1/
      c0=0.

* Index to use for biology
      ibdt=taum1

      if (n.eq.3 .and. j.eq.2) then
         itimes = itimes+1
      endif
      if (n.eq.3) then
         call tx_growth(j)
      endif

      ! rjm Needs to be fixed to do c14

c C14 Tracer
      if (n.eq. n_c14) then
         c14_atm=1e2
         c14_lambda = -3.8560619e-12
         c14_gamma = 3.215e-7/dzt(1)*1e2
         
         do i=istrt,iend
            source(i,1) = (c14_atm - tb(i,1,n) )*c14_gamma
            fluxgas(i,j,n-2) = source(i,1)
         enddo
         do k=2,km
            do i=istrt,iend
               source(i,k) = c14_lambda*tb(i,k,n)
            enddo
         enddo
         return
      endif

c --------------------------------------------

************************************************************************
* ****  STEP 1 -->  MODIFY ATMOSPHERIC PROPERTIES IF APPROPRIATE  **** *  
* ****                     & COMPUTE CARBON CHEMISTRY             **** *                    
************************************************************************

      if (n.eq.3) then ! --> bgc is on

*   1.1  --  Modify pCO2 and C13 in atmosphere if appropriate
      
         ! compute year of integration 
         ttt = itt/float(knitd*365)
         
         ! Do the following calcs at the very start of the timestep
         if (j.eq.2) then

            ! Get atmospheric partial pressure of CO2
            pco2a=co2_history(ttt)
            if (lcouple) then
               call obgc_atmco2(pco2a)
            endif
            
            ! Get the fraction of c13 to total co2 in the atmosphere
            fc13a = 0.011164382 - 2.0225714e-7*max(0.,(pco2a-280.) ) 
#ifdef dmr_atm_box
            fc13a=rc_atm(pco2a)
#endif

            print *,'Timestep number =',itt,'   Year =',ttt,
     +              '   Timesteps per day =',knitd
            print *,'Atmospheric CO2 =',pco2a,'d13C frac =',fc13a

         endif ! --> j.eq.2 (means first latitude)

         tott = ttt   ! start the total time from 2004


*   1.2  --  Compute carbon chemistry terms if dic tracer is used

         if (n_dic.gt.0) then  
            do i = istrt,iend 
               if (kmt(i,j).gt.0) then

                  ! Collect tracers required for the calculation
                  CCtem(1)=max(-2.0, min(40.0,tb(i,1,1)))  ! TEMPERATURE
                  CCsal(1)=max(0.0,(tb(i,1,2)+.035)*1e3) ! SALINITY
                  
                  if (n_alk .le. 0) then        ! ALKALINITY
                     CCalk(1) = 2377  ! 2431 old value
                  else
                     CCalk(1) = max(0.0,tb(i,1,n_alk))
                  endif
                  
                  CCdic(1) = max(0.0,tb(i,1,n_dic))         ! TOTAL CARBON
                  
                  if (n_pho.gt.0) then         ! PHOSPHATE 
                     CCpo4(1) = max(tb(i,1,n_pho),0.)
                  else
                     CCpo4(1) = 0.0
                  endif
                  
                  CCsil(1) = max(silicate(i,j),0.)    ! SILICATE
                  
                  ! Collect other necessary terms
                  CCpatm(1) = sbcbio(i,j,6)
                  CClat(1) = radian*phi(j)
                  CClon(1) = (i-1)*2.8125
                  CCdep(1) = 12.5

                  ! Convert input tracers to mol/m3
                  CCalk(1) = CCalk(1)*1e-3
                  CCdic(1) = CCdic(1)*1e-3
                  CCpo4(1) = CCpo4(1)*1e-3
                  CCsil(1) = CCsil(1)*1e-3
  
                  
!     OUTPUT variables:
!     =================
!     ph   = pH on total scale
!     pco2 = CO2 partial pressure (uatm)
!     fco2 = CO2 fugacity (uatm)
!     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
!     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
!     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
!     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
!     OmegaC = Omega for calcite, i.e., the   calcite saturation state
!     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
!     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
!     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
!     tempis  = in-situ temperature [degrees C]

                  call vars(CCph, CCpco2, CCfco2, CCco2, CChco3, CCco3, 
     +                      CCOmegaA, CCOmegaC, CCBetaD, CCrhoSW, CCp, 
     +                      CCtempis, 
     +                      CCtem, CCsal, CCalk, CCdic, CCsil, CCpo4,
     +                      CCpatm, CCdep, CClat, N=1,             
     +                      optCON='mol/m3',optT='Tpot   ',optP='m ')

                  ! store particular parameters for later use
                  pco2o(i,j,n_dic-2) = CCpco2(1)
                  co2sw(i,j) = CCco2(1)
                  omegaar(i,j) = CCOmegaA(1)
                  omegaca(i,j) = CCOmegaC(1)
                  fco3(i,j) = CCco3(1)/(CCco2(1)+CChco3(1)+CCco3(1))
                  rhosw(i,j) = CCrhoSW(1)


               endif  ! masks
            enddo  ! i - loop
         endif  ! end n_dic if statement

* ****  END STEP 1  **** *


         
**************************************************
* ****  STEP 2 -->  COMPUTE NEW PRODUCTION  **** *
**************************************************
         
*       Units are mmol/m3/s  cm
*       = 1e-2 mmol/m2/s

         
         do i=istrt,iend
         if (kmt(i,j).gt.0) then
       
*   2.1.  --  Calculate the C:H:O:N:P stoichiometries of organic matter
            
            if (vary_stoich) then
               call stoich(i,j, max(0.0,tb(i,1,n_pho)), 
     +                          max(0.0,tb(i,1,n_no3)) )  ! see Galbraith & Martiny 2015
            else
               carb2P(i,j) = 106
               NtoP(i,j) = 16
               
               ! Implicit assumption that all organic matter is carbohydrate
               HtoP(i,j) = 2.0*carb2P(i,j) + 3.0*NtoP(i,j) + 3.0
               OtoP(i,j) = carb2P(i,j) + 4.0

               ! Calculate O2 and NO3 demand for oxic and suboxic remineralisation
               !    Both assume "complete" remineralisation, whereby
               !    organics are converted instantly to end products
               !    PO4, CO2 and NO3
            endif
            o2_rem(i,j) = -(carb2P(i,j) + 0.25*HtoP(i,j) - 
     +                      0.50*OtoP(i,j) - 0.75*NtoP(i,j) + 1.25)
     +                      - 2.0*NtoP(i,j)
               
            no3_rem(i,j) = -(0.8 * (carb2P(i,j) + 0.25*HtoP(i,j) - 
     +                      0.50*OtoP(i,j) - 0.75*NtoP(i,j) + 1.25)
     +                      + 0.6*NtoP(i,j) )
            
            ! Calculate fractionation of c13
            if (vary_frac13C .and. n_dic.gt.0) then ! uses growth rate and [CO2]aq (Laws et al., 1995, Geo. et Cos. Acta)
               c13toc(i,j) = carb2P(i,j) * min(0.985, max(0.97,
     +                       (1.0 - ((0.371 - (Vmax(i,j)*86400.0 / 
     +                       (co2sw(i,j)*1e3)))/0.015)*1e-3) ))
            else
               c13toc(i,j) = carb2P(i,j)*103.8/106.0
            endif
            
            ! read in phosphate concentrations
            if (n_pho.gt.0) then
               phs = tb(i,1,n_pho) !phosphate at the surface  
            else
               phs = c0
            endif
            if (phs .le. 0.02) phs=c0
         
            ! read in nitrate concentrations 
            if (n_no3.gt.0) then
               no3 = max(0.0, tb(i,1,n_no3))
            else
               no3 = c0
            endif
            if (no3 .le. 0.1) no3=c0

            
            
            if (fix) then
*********************************************
***   2.2  --  BEGIN NITROGEN FIXATION    ***

*   2.21  --  Determine oxygen and nitrate consumption of Nfixer OM            
               CP_fix = 331.0 !Karl & Letelier 2008 MEPS
               NP_fix = 50.0  ! Mills & Arrigo (2010) Nature Geoscience
               HP_fix = 2.0*CP_fix + 3.0*NP_fix + 3.0
               OP_fix = CP_fix + 4.0
               
               o2_rem_f = -(CP_fix + 0.25*HP_fix - 
     +                         0.50*OP_fix - 0.75*NP_fix + 1.25)
     +                         - 2.0*NP_fix     ! -431 O2:P
               no3_rem_f = -(0.8 * (CP_fix + 0.25*HP_fix - 
     +                         0.50*OP_fix - 0.75*NP_fix + 1.25)
     +                         + 0.6*NP_fix )   ! -294.8 NO3:P
                    

*   2.23  --  Calculate temperature limitation (0-1)
               temp = max(-2.0, min(40.0,tb(i,1,1)))
               !Fmax is in doublings per day
               Fmax = max(0.01,  ( (-0.0042 * temp**2.0) +
     +                   (0.2253 * temp) - 2.7819 ) ) !from Kreist2015
                
*   2.24  --  Calculate nutrient limitation
               if (n_fe.gt.0) then
                  dFe = tb(i,1,n_fe)
                  R_fix(i,j) = max(0.01, ! n fixation happens everywhere
     +                         min(exp(-no3),max(0.,tanh(2*dFe-Fe_kd))))
               else
                  R_fix(i,j) = max(0.01,
     +                         min(exp(-no3),(phs/(phs+p_kd)) ))
               endif

                  
*   2.25  --  Calculate PONfix production

               pop_fix(i,1,j) = s_npp_fix*Fmax*dzt(1)*R_fix(i,j)
     +                          *(1.-sbcbio(i,j,3))*(1./86400.)
               
         
*   2.26  --  Determine percentage of available phosphate for diazotrophs            

               if (pop_fix(i,1,j).gt.0.0) then
               ! 1. check if enough PO4 is available (convert to mmol PO4 m-3 per timestep) 
                  ratio = phs/(pop_fix(i,1,j)*dizt(1)*c2dtts)
                  pop_fix(i,1,j) = min(1.0,ratio) * pop_fix(i,1,j)
                  phs = phs-pop_fix(i,1,j)
               endif

            endif
            

***   2.2  --  END NITROGEN FIXATION    ***
*******************************************
          


*   2.3  --  Calculate Export Production (Particulate Organic Phosphate)            

     
            if (MichaelisMenton) then
               if (n_fe.gt.0 .and. n_no3.gt.0) then 
                  dfe = max(tb(i,1,n_fe),0.)
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min( phs/(phs+p_k), dfe/(dfe+fe_k),
     +                              no3/(no3+n_k), 
     +                              fpgr_temp(i,j) )
        
               elseif (n_fe.gt.0 .and. n_no3.eq.0) then
                  dfe = max(tb(i,1,n_fe),0.)
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min( phs/(phs+p_k), dfe/(dfe+fe_k),
     +                              fpgr_temp(i,j) )
               
               elseif (n_fe.eq.0 .and. n_no3.gt.0) then 
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min( phs/(phs+p_k),
     +                              no3/(no3+n_k), 
     +                              fpgr_temp(i,j) )
               
               else
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min( phs /(phs+p_k), 
     +                              fpgr_temp(i,j) )
               
               endif
            endif


            
            if (OptimalUptake) then

               phs = max(1e-10, phs) !phosphate at the surface  
               no3 = max(1e-10, no3) !nitrate at the surface  
               temp = tb(i,1,1)
               
               if (n_no3.gt.0) then 
                  fa = max( (1.+(no3/VoA)**0.5)**(-1.), 
     +                      (1.+((phs*NtoP(i,j))/VoA)**0.5)**(-1.)  )
               else
                  fa = (1.+((phs*NtoP(i,j))/VoA)**0.5)**(-1.)
               endif
             
               if (NtoP(i,j).gt.0.0) then 
                  Pou = phs/( (phs/(1.-fa)) + (VoA/(NtoP(i,j)*fa)) ) 
                  Nou = no3/( (no3/(1.-fa)) + (VoA/fa) )
               else
                  Pou = 0.0 
                  Nou = 0.0 
               endif
              
               if (n_fe.gt.0 .and. n_no3.gt.0) then 
                  dfe = max(tb(i,1,n_fe),0.)
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                       min(Pou,Nou,dfe/(dfe+fe_k),fpgr_temp(i,j))
        
               elseif (n_fe.gt.0 .and. n_no3.eq.0) then
                  dfe = max(tb(i,1,n_fe),0.)
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min(Pou, dfe/(dfe+fe_k), fpgr_temp(i,j))
               
               elseif (n_fe.eq.0 .and. n_no3.gt.0) then 
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min( Pou, Nou, fpgr_temp(i,j) )
               
               else
                  pop(i,1,j) = s_npp*vmax(i,j)*dzt(1)* 
     +                         min( Pou, fpgr_temp(i,j) )
               
               endif
            endif 

            
*   2.4  --  Ensure that pop can only remove what is available
            
            if (pop(i,1,j).gt.0.0) then
               ! 1. check if enough PO4 is available (convert to mmol PO4 m-3 per timestep) 
               pop(i,1,j) = min(1.0,
     +                     (phs/(pop(i,1,j)*dizt(1)*c2dtts)) )
     +                      * pop(i,1,j)
               ! 2. check if enough NO3 is available
               if (n_no3.gt.0) then
                  pop(i,1,j) = min(1.0, max(0., ( 
     +                         (tb(i,1,n_no3)-0.1) / (pop(i,1,j)*
     +                         dizt(1)*c2dtts*NtoP(i,j)) ) ))*pop(i,1,j)
               endif
            endif

            
*   2.5  --  Determine a measure of C:P based on variable pCO2 concentrations
            
            if (ivredc.eq.0) then
               vredf = 1.
            endif
            if (ivredc.eq.1) then 
               vredf= (1. + (pco2o(i,j,n_dic-2) - 280.)*2./700./6.6)
            endif
            if (ivredc.eq.2) then
               vredf= (1. + (pco2a - 280.)*2./700./6.6)
            endif


*   2.6  --  COMPUTE POC and PON FROM POP 
         !Note poc is not used in source equation only in pic scaling below
            poc(i,1,j) = pop(i,1,j)*carb2P(i,j)*vredf  

*   2.7  --  COMPUTE CaCO3 SHELL PRODUCTION AS A CONSTANT VALUE OF POP
            if (ivredpic .eq.1) then
               pic(i,1,j) = rain_ratio*poc(i,1,j)
            endif
            if (ivredpic .eq.0) then 
               pic(i,1,j) = rain_ratio*pop(i,1,j)*carb2P(i,j)
            endif


*   2.8  --  CALCULATE CaCO3 PRODUCTION SCALING DEPENDENT ON omegaca
            if (ipic.ge.1) then
               if (omegaca(i,j).gt.1) then
                  var_satst(i,j) = (omegaca(i,j)-1)**pwr_sat ! from Ridgwell 2007 Biogeosciences
               else
                  var_satst(i,j) = 0  !no pic export
               endif  !end test for omega

               ! scaled based on poc which includes the vred term !
               if (ivredpic.eq.1) then
                  pic(i,1,j) = var_satst(i,j)*rain_ratio*poc(i,1,j)
               endif
  
               ! scaled based on pop with constant c/p 
               if (ivredpic .eq. 0) then
                  pic(i,1,j) = var_satst(i,j)*rain_ratio*pop(i,1,j)*
     +                         carb2P(i,j)
               endif
            endif ! end omegaca dependency

*   2.9  --  COMPUTE THE OPAL PRODUCTION
            if (n_si.ne.0) then
               opal(i,1,j) = pop(i,1,j)*tb(i,1,n_si)/(tb(i,1,n_si)+4)
            endif

         endif
         enddo  ! end looping over the ith dimension (longitude)

* ****  END STEP 2  **** *
         

**********************************************************************
* ****  STEP 3 -->  COMPUTE REMINERALISATION OF ORGANIC MATTER **** **
*                       - specifically: POP, PIC, OPAL               *
**********************************************************************
         
         if (ReminTemp_Marsay) then
*   3.1  --  Apply temperature dependence to exponent (Marsay 2015 PNAS)
            
         do i=istrt,iend
         if (kmt(i,j).ne.0) then

            ! find depth interval over which temp will be averaged
            dep = 0.0
            do k=1,9
               dep = dep + dzt(k)
            enddo
            ! find average temperature over mesopelagic zone
            temp = 0.0
            do k=1,9
               temp = temp + tb(i,k,1)*(dzt(k)/dep)
            enddo

            ! apply temp-dependence to power law exponent
            Mcurve_b(i,j) = max(0.6, 0.062*temp + 0.303)*(-1.)
            
            ! recalculate the Martin Curve for pop and pop_fix
            sedtot(i,j) = 0.0
            do k=1,km
               sedtot(i,j) = sedtot(i,j) + sgb_frac(i,j,k)
               fmin_pop(i,j,k) = min(1.0, ((zt(k)+dzt(k)*0.5)*zpop)
     +                           **Mcurve_b(i,j)) * (1.0-sedtot(i,j))
               fmin_dia(i,j,k) = min(1.0, ((zt(k)+dzt(k)*0.5)*zdia)
     +                           **Mcurve_b(i,j)) * (1.0-sedtot(i,j))
            enddo
            if (.not.sedfluxes) then
               fmin_pop(i,j,km) = 0.0
               fmin_dia(i,j,km) = 0.0
            endif
         
         endif
         enddo

         endif


         if (ReminPico) then
*   3.1  --  Apply picoplankton dependence to exponent
            
            if (EPmax.gt.0.0) then
               do i=istrt,iend
               if (kmt(i,j).ne.0) then
                  
                  ! find export production (mg C m-2 hr-1)
                  EP(i,j) = poc(i,1,j) * (1e-2*12.0*(c2dtts/2.0)) /
     +                      ((c2dtts/2.0)/3600.0)
               
                  ! find fraction of picoplankton based on export production
                  Fpico = 0.51 - 0.26*(EP(i,j)/(EP(i,j) + EPmax/2.0)) 
              
                  ! find transfer efficiency
                  Teff = 0.47 - 0.81*Fpico 
                 
                  if (ReminTemp_Q10) then
                  ! 3.3  --  Apply temperature dependence
            
                     ! find depth interval over which temp will be averaged
                     dep = 0.0
                     do k=1,10
                        dep = dep + dzt(k)
                     enddo
                     ! find average temperature over mesopelagic zone
                     temp = 0.0
                     do k=1,10
                        temp = temp + tb(i,k,1)*(dzt(k)/dep)
                     enddo

                     ! apply temp-dependence based on Q10 (Matsumoto 2007)
                     Teff = max(1e-6, Teff-0.02*Q10**((temp-T_rem)/10.))

                  endif
            
                  ! find b exponent
                  Mcurve_b(i,j) = log10(Teff)/log10(10.0)
            
                  ! recalculate the Martin Curve
                  sedtot(i,j) = 0.0
                  do k=1,km
                     sedtot(i,j) = sedtot(i,j) + sgb_frac(i,j,k)
                     fmin_pop(i,j,k) =min(1.0, ((zt(k)+dzt(k)*0.5)*zpop)
     +                                **Mcurve_b(i,j))*(1.0-sedtot(i,j))
                     fmin_dia(i,j,k) =min(1.0, ((zt(k)+dzt(k)*0.5)*zdia)
     +                                **Mcurve_b(i,j))*(1.0-sedtot(i,j))
                  enddo
                  if (.not.sedfluxes) then
                     fmin_pop(i,j,km) = 0.0
                     fmin_dia(i,j,km) = 0.0
                  endif

               endif
               enddo
            endif

         endif
            
         
*   3.3  --  COMPUTE REMINERALISATION OF POP, PIC, OPAL
         do i=istrt,iend
         if (kmt(i,j).gt.0) then
            do k=2,kmt(i,j)
               pop(i,k,j) = -pop(i,1,j)*
     +                      (fmin_pop(i,j,k-1)-fmin_pop(i,j,k) )
               pic(i,k,j) = -pic(i,1,j) *
     +                      (fmin_pic(i,j,k-1) - fmin_pic(i,j,k) )
               if (fix) then
               pop_fix(i,k,j) = -pop_fix(i,1,j)*
     +                          (fmin_dia(i,j,k-1)-fmin_dia(i,j,k) )
               endif
               if (n_si.ne.0) then
                  opal(i,k,j) = -opal(i,1,j)*
     +                          (fmin_opal(i,j,k-1) - fmin_opal(i,j,k) )
               endif
          
            enddo

*   3.4  --  TREAT ORGANIC MATTER AT OCEAN BOTTOM (SEDIMENT OR WATER)
       
            if (sedfluxes) then
            
               ! Remineralise remaining organic matter in sediments
               op_sed(i,kmt(i,j),j) = -pop(i,1,j)*
     +                                (fmin_pop(i,j,kmt(i,j))-0.0)
               ic_sed(i,kmt(i,j),j) = -pic(i,1,j)*
     +                                (fmin_pic(i,j,kmt(i,j))-0.0)
               if (fix) then
                  op_sed_f(i,kmt(i,j),j) = -pop_fix(i,1,j)*
     +                                     (fmin_dia(i,j,kmt(i,j))-0.0)
               endif
               if (n_si.ne.0) then
                  si_sed(i,kmt(i,j),j) =-opal(i,1,j)*
     +                                  (fmin_opal(i,j,kmt(i,j))-0.0)
               endif

            else ! not sedfluxes
        
               ! Remineralise all remaining organic matter in water column
               pop(i,kmt(i,j),j) = -pop(i,1,j)*
     +                            (fmin_pop(i,j,kmt(i,j)-1)-0.0)
               op_sed(i,kmt(i,j),j) = 0.0
               pic(i,kmt(i,j),j) = -pic(i,1,j) *
     +                             (fmin_pic(i,j,kmt(i,j)-1)-0.0)
               ic_sed(i,kmt(i,j),j) = 0.0
               if (fix) then
                  pop_fix(i,kmt(i,j),j) = -pop_fix(i,1,j)*
     +                                    (fmin_dia(i,j,kmt(i,j)-1)-0.0)
                  op_sed_f(i,kmt(i,j),j) = 0.0
               endif
               if (n_si.ne.0) then
                  opal(i,kmt(i,j),j) = -opal(i,1,j)*
     +                                (fmin_opal(i,j,kmt(i,j)-1)-0.0)
                  si_sed(i,kmt(i,j),j) = 0.0
               endif

            endif ! sedfluxes


*   3.5  --  CHECK CONSERVATION OF MASS THROUGH WATER COLUMN
            
         tot_b = pop(i,1,j)
         tot_a = 0.0
         do k=2,kmt(i,j)-1
            tot_a = tot_a + pop(i,k,j)
         enddo
         tot_a = tot_a + pop(i,kmt(i,j),j) + op_sed(i,kmt(i,j),j)
            
         dif = tot_b+tot_a
         tot = tot_b-tot_a
         if (ABS(dif).gt.epsilon(tot)) then
            print*, " "
            print*, "   REMINERALISATION (POP) SCHEME NOT     "
            print*, "       CONSERVING MATTER                  "
            print*, " "
            print*," POP at surface = ",tot_b
            print*," POP sub-surface = ",tot_a
            print*," Fraction POP to sediment = ",sedtot(i,j)
            print*, " "
            print*, 'i = ',i,' j = ',j,'k = ',kmt(i,j)
            print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
            print*, " "
            print*, " pop(i,1,j) =", pop(i,1,j)
            do k=2,km
              print*, "k = ", k," pop(i,k,j) =", pop(i,k,j)
            enddo
            print*, " "
            print*, "              STOPPING MODEL                   "
            print*, " "
            print*, " "
            stop
         endif
         
         tot_b = pic(i,1,j)
         tot_a = 0.0
         do k=2,kmt(i,j)-1
            tot_a = tot_a + pic(i,k,j)
         enddo
         tot_a = tot_a + pic(i,kmt(i,j),j) + ic_sed(i,kmt(i,j),j)
            
         dif = tot_b+tot_a
         tot = tot_b-tot_a
         if (ABS(dif).gt.epsilon(tot)) then
            print*, " "
            print*, "   REMINERALISATION (PIC) SCHEME NOT     "
            print*, "       CONSERVING MATTER                  "
            print*, " "
            print*," PIC at surface = ",tot_b
            print*," PIC sub-surface = ",tot_a
            print*," Fraction PIC to sediment = ", sedtot(i,j)
            print*, " "
            print*, 'i = ',i,' j = ',j,'k = ',kmt(i,j)
            print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
            print*, " "
            print*, " pic(i,1,j) =", pic(i,1,j)
            do k=2,km
              print*, "k = ", k," pic(i,k,j) =", pic(i,k,j)
            enddo
            print*, " "
            print*, "              STOPPING MODEL                   "
            print*, " "
            print*, " "
            stop
         endif

         endif ! kmt(i,j).ne.0
         enddo ! i loop
      
      endif  ! end calculations on first tracer


* ****  END STEP 3  **** *


      
********************************************************
* ****  STEP 4 -->  COMPUTE AIR-SEA GAS EXCHANGE  **** *
********************************************************

*   4.1  --  Oxygen tracers
      if (n.eq.n_oxy .or. n.eq.n_tou) then    !setup for oxygen but the code can do other gases
         ng = igas(n-2)

         do i = istrt,iend
            if (kmt(i,j).ne.0) then 
         
               ! Units are 1e-2 mmol/m2/s
               temp = max(-2.0, min(40.0,tb(i,1,1)))
               saln = max(0.0, min(45.0,(tb(i,1,2)+.035)*1e3))

               ! Gas transfer velocity [cm/s] 
               ! Updated according to OMIP6 biogeochemical protocols
               x_flux = 0.251 * (schmidt_no(ng,temp)/660)**(-0.5)
     +                  *sbcbio(i,j,1)*(1./3600) ! conversion
               
               !atmosphere and surface ocean O2 [mmol/cm3]
               oxs = (tb(i,1,n)+tr_off(n))
               oxa = (o2sat(temp,saln)*sbcbio(i,j,6))
            
               ! exchange in mmol/cm2/s
               fluxgas(i,j,n-2) = x_flux*(oxa-oxs)

               ! for OMIP6 output calculate partial pressure of O2 in seawater
               !pO2(i,j) = oxs/(oxa/(0.20946*(1.-pH2O(temp,saln))))

            endif
         enddo

      else

*   4.2  --  Phosphate
         if (n.eq.n_pho) then
            do i = istrt,iend
               fluxgas(i,j,n-2) = terr_input
               ! phosphate fertilization
               rtd = 180./3.14159
               spy = 365.*3600*24
               yy = phit(j)*rtd
               if (yy.ge.-50 .and. yy.le.-20 .and. fm(i,1).ne.0) then
                  fluxgas(i,j,n-2) = fluxgas(i,j,n-2) + pdump / spy 
               endif
            enddo
         endif

*   4.3  --  Nitrate         
         if (n.eq.n_no3) then
            do i = istrt,iend
               if (atmdep) then
                  fluxgas(i,j,n-2) = sbcbio(i,j,5)
               else
                  fluxgas(i,j,n-2) = 0.0
               endif
            enddo
         endif

*   4.4  --  Nitrogen 15          
         if (n.eq.n_n15) then
            do i = istrt,iend
               if (atmdep) then
                  fluxgas(i,j,n-2) = (1.-(1./(1. + Eatm*1e-3 + 1.)))
     +                                *sbcbio(i,j,5)
               else
                  fluxgas(i,j,n-2) = 0.0
               endif
            enddo
         endif

         
*   4.5  --  Alkalinity         
         if (n.eq.n_alk) then
            do i = istrt,iend
               if (atmdep) then
                fluxgas(i,j,n-2) = (-terr_input*NtoP(i,j)-sbcbio(i,j,5))
               else
                fluxgas(i,j,n-2) = -terr_input*NtoP(i,j)
               endif
            enddo
         endif

*   4.6  --  Iron         
         if (n.eq.n_fe) then
            do i = istrt,iend
               fluxgas(i,j,n-2) = sbcbio(i,j,4) 
            enddo
         endif

*   4.7  --  Age Tracers         
         if (n.eq.n_age) then
            do i = istrt,iend
               fluxgas(i,j,n-2) = 0.0 
            enddo
         endif
         if (n.eq.n_tomz) then
            do i = istrt,iend
               fluxgas(i,j,n-2) = 0.0 
            enddo
         endif

*   4.8  --  Nitrous Oxide         
         if (n.eq.n_n2o) then
            ng = igas(n-2)
            do i = istrt,iend
               if (kmt(i,j).gt.0) then   ! Only calculate for ocean points
                  temp = max(-2.0, min(40.0, tb(i,1,1)))
                  saln = max(0., min(45.0,(tb(i,1,2)+.035)*1e3))
                  
                  !1. Calculate Gas transfer velocity (Kw) [cm/s]
                  x_N2O(i,j) = 0.251*(schmidt_no(ng,temp)/660.0)**(-0.5)
     +                         *sbcbio(i,j,1)*(1./3600) !cm/hr --> cm/s

                  !2. Get surface water concentration (umol/cm3)
                  N2Oo = (tb(i,1,n))
                 
                  !3. Calculate N2O concentration in moist air (umol/cm3)
                  N2Oa = (pN2O*1e-9 * phizero(temp,saln,ng)*1e3*
     +                   sbcbio(i,j,6))
                  
                  ! exchange in umol/m2/s
                  fluxgas(i,j,n-2) = x_N2O(i,j)*(N2Oa - N2Oo)
                  
               endif
            enddo
         endif
         
*   4.9  --  Dissolved Inorganic Carbon         
         ! For CO2 exchange the units are umol/l *cm/s  ==>  1e-2 mmol/m2/s
         if (n.eq.n_dic .or. n.eq.n_c13 .or. n.eq.n_c14) then 
            ng = igas(n-2)
            do i = istrt,iend
         
               if (kmt(i,j).gt.0) then   ! Only calculate for ocean points
                  temp = max(-2.0, min(40.0,tb(i,1,1)))
                  saln = max(0., min(45.0,(tb(i,1,2)+.035)*1e3))
                  
                  ! 1. - calculate gas transfer velocity [cm/s]
                     x_co2(i,j) = 0.251 * (660.0/scco2(temp))**0.5 *
     +                            sbcbio(i,j,1)*(1.0/3600.0)

                  ! 2. - calculate the atmopsheric mole fraction CO2
                     tk = 273.15 + temp
                     pH2O = exp(24.4543 - 67.4509*(100.0/tk) -
     +                      4.8489*log(tk/100) - 0.000544d0*saln)
                     pco2atm = (sbcbio(i,j,6)-pH2O)*pco2a 

                  ! 3. - calculate fCO2atm [uatm] from mole fraction CO2
                     fco2atm = fco2(pco2atm, temp, sbcbio(i,j,6),12.5)

                  ! 4. surface K0 [(mol/kg) / atm] at T, S of surface water
                     invtk = 1.0/tk
                     tmp = 9345.17*invtk-60.2409+23.3585*log(tk/100.0)
                     K0 = exp(tmp + saln*(0.023517 - 
     +                        0.00023656*tk + 0.0047036e-4*tk*tk) )

                  ! 5. compute "Atmospheric" [CO2*], air-sea CO2 flux, & Delta pCO2
                     
                     !Equil. [CO2*] for atm CO2 at Patm & sfc-water T,S 
                     ! from [mol/m3] to [mmol/cm3]
                     co2starair = K0*fco2atm*1.0e-6*rhosw(i,j)*1e3
                     
                     !Oceanic [CO2*] from [mol/m3] to [mmol/cm3] from vars.f90
                     co2star = co2sw(i,j)*1e3
                    
                     dpco2 = pco2o(i,j,n-1) - pco2atm
                     !Delta pCO2 (oceanic - atmospheric pCO2) [uatm]


                  if (n.eq.n_dic ) then
                  

                     ! 6. compute co2flux across air-sea interface in mmol/cm2/s
                     fluxgas(i,j,n-2) = x_co2(i,j)* 
     +                                 (co2starair - co2star)


                  endif ! end dic
         
*   4.11  --  C13 Exchange (must include CO2 exchange)
                  if (n.eq.n_c13 ) then
                  
                     ! Temperature Dependent fractionation
                     f_dis = (0.0049*temp-1.31)/1000.0 + 1.0
                     f_spe = (0.0144*temp*fco3(i,j) 
     +                        - 0.107*temp + 10.53)/1000.0 + 1.0

                     ratio = tb(i,1,n)/tb(i,1,n_dic)

                     fluxgas(i,j,n-2) = x_co2(i,j)*0.99912*f_dis *
     +                                  (fc13a*co2starair - 
     +                                   co2star*ratio/f_spe)

                  endif  ! endif for c13  

**   4.12  --  C14 Exchange (must include CO2 exchange)
*                  if (n.eq.n_c14 ) then
*                  
*                     ratio = tb(i,1,n)/tb(i,1,n_dic)
*                     c14o = co2star*ratio  ! mmol/m3
*                     del14c_atm = 0.0
*                     c14a = co2starair * (1.0 + del14c_atm*1e-3) ! mmol/m3
*
*                     fluxgas(i,j,n-2) = x_co2(i,j)*(c14a-c14o)
*
*                  endif  ! endif for c14  


               endif  ! endif fm>0

            enddo  ! end i loop
         endif  ! endif for carbon
      endif  ! endif airsea flux test
 
* ****  END STEP 4  **** *


******************************************************
* ****  STEP 5 -->  COMPUTE THE VIRTUAL FLUXES  **** *
******************************************************
      ! assumes fluxgas has been updated prior to the following
      ! for n_pho and n_alk the terrerial input set fluxgas

      do i = istrt,iend
         if (.not.lcouple) then
            if (kmt(i,j).ne.0) then 
            salf = flux(i,j,2) ! salf = change in salinity in top layer per second
            
               fluxgas(i,j,n-2) = fluxgas(i,j,n-2) + salf 
     +                            * tglobal(n)/0.035 * dzt(1)
         
            endif
         else
            if (kmt(i,j).ne.0) then           
            salf = flux(i,j,2)! salf = change in salinity * rhow * dz(1) *0.01
            WMTOP = 1.  ! v1.2 fixes units  RHOW*DZ(1)*0.01
            
               fluxgas(i,j,n-2) = fluxgas(i,j,n-2) + salf 
     +                            * tglobal(n)/0.035*dzt(1)/wmtop
         
            endif
         endif
      enddo  
     
* ****  END STEP 5  **** *

      
*****************************************************************
* ****  STEP 6 --> COMPUTE OXIC & SUBOXIC REMINERALISATION **** *
* ****                      IN THE WATER COLUMN            **** *
*****************************************************************

      if (n.eq.3) then ! --> if biogeochemistry is on
      if (ReminRelax) then 
      
* -------------------------------------------------------------------- *
*                                                                      *
*   The following section of code ensures that where oxygen is         *
*   limiting, and denitrification cannot account for the reminer-      *
*   alisation of organic matter (OM), remineralisation is also         *
*   limited.                                                           *
*                                                                      *
*   Oxygen should therefore not removed past zero. OM is passed to     *
*   the next depth for consumption. If this cannot be achieved to      *
*   completion, then the remaining OM is passed on again until         *
*   either:                                                            *
*       - all OM can be remineralised                                  *
*       - the last depth is encounted (kmt(i,j)), where all OM is      *
*         remineralised                                                *
*                                                                      *
* -------------------------------------------------------------------- *

      ! beforehand, select only organics not hitting sediment and
      ! collect the amount of organics hitting sediment
      do i = istrt,iend
      if (kmt(i,j).ne.0) then
         do k = 1,kmt(i,j)
    
            if (k.lt.kmt(i,j)) then

            op_tot(i,k,j) = pop(i,k,j) * (1.0-sgb_frac(i,j,k)) 
            op_sed(i,k,j) = pop(i,k,j) * sgb_frac(i,j,k) 

            if (fix) then
               op_tot_f(i,k,j) = pop_fix(i,k,j) * (1.0-sgb_frac(i,j,k))
               op_sed_f(i,k,j) = pop_fix(i,k,j) * sgb_frac(i,j,k)
            else
               op_tot_f(i,k,j) = 0.0
               op_sed_f(i,k,j) = 0.0
            endif

            else
            ! Calculate additional organics to sediment in final depth
            op_tot(i,k,j) = pop(i,k,j) * (1.0-sgb_frac(i,j,k)) 
            op_sed(i,k,j) = op_sed(i,k,j) + (pop(i,k,j)*sgb_frac(i,j,k))

            if (fix) then
               op_tot_f(i,k,j) = pop_fix(i,k,j) * (1.0-sgb_frac(i,j,k))
               op_sed_f(i,k,j) = op_sed_f(i,k,j) + 
     +                           (pop_fix(i,k,j) * sgb_frac(i,j,k))
            else
               op_tot_f(i,k,j) = 0.0
               op_sed_f(i,k,j) = 0.0
            endif
            
            endif
              
         enddo
      endif
      enddo


      do i = istrt,iend
      if (kmt(i,j).ne.0) then
         do k = 1,kmt(i,j)
      
            oxy = max(0.0, tb(i,k,n_oxy))
      
            if (n_no3.gt.0 .and. den) then

               no3 = max(0.0, tb(i,k,n_no3))

*   6.1  --  Using oxygen, calculate the strength of denitrification

               if (oxy.lt.Den_o2lim) then
                     
                  if (Den_type.eq."Linear ") then
                     R_den(i,k,j) = min(1., max(0.,
     +                              (Den_o2lim - oxy)/Den_o2lim ))
                  endif
                        
                  if (Den_type.eq."QuadPos") then
                     R_den(i,k,j) = min(1.0, max(0.0,
     +                             ( ((Den_o2lim - oxy)/
     +                                 Den_o2lim)**2. ) ) )
                  endif

                  if (Den_type.eq."QuadNeg") then
                     R_den(i,k,j) = min(1.0, max(0.0,
     +                             ( -(oxy/Den_o2lim)**2 + 1.0)))
                  endif
                    
                  if (Den_type.eq."Sigmoid") then
                     R_den(i,k,j) = min(1.0, max(0.0,
     +                             (1./(1.0 - exp(-Den_o2lim/2.0) +
     +                              exp(oxy-(Den_o2lim/2.)) )) ))
                  endif
                     
               else
                  R_den(i,k,j) = 0.0
               endif

          
*   6.2  --  To ensure that nitrate consumption slows as no3 becomes
*            less available, and therefore ensure that no3 never really
*            reaches zero in the ocean interior, I here introduce a
*            reduction in denitrification as NO3 becomes limiting
               
               if (no3.le.Den_nlim .and. R_den(i,k,j).gt.0.0) then
                  denrelax = .5 + .5*tanh(0.25*no3-Den_nlim*0.25-2.5)
               else
                  denrelax = 1.0
               endif
            
               if (denrelax.lt.R_den(i,k,j)) then
                  R_den(i,k,j) = denrelax
               endif
            
          
            else ! if den is FALSE
               R_den(i,k,j) = 0.0
            endif

            
*   6.3  --  Calculate amount of pop to be remineralised by oxic and
*            suboxic (denitrification) remineralisation
 
            op_rem(i,k,j) = op_tot(i,k,j)*(1.0-R_den(i,k,j))
            op_den(i,k,j) = op_tot(i,k,j)*R_den(i,k,j)

            if (fix) then
               op_rem_f(i,k,j) = op_tot_f(i,k,j)*(1.0-R_den(i,k,j))
               op_den_f(i,k,j) = op_tot_f(i,k,j)*R_den(i,k,j)
            else
               op_rem_f(i,k,j) = 0.0
               op_den_f(i,k,j) = 0.0
            endif


*   6.4  --  Now, determine whether the organic matter allocated to
*            oxygen can in fact be remineralised by the oxygen that is 
*            available

            if (-op_rem(i,k,j).gt.0.0) then
               R_Olim(i,k,j) = min(1.0, max(0.0, ( (oxy-0.5) / (
     +                         op_rem(i,k,j)*dizt(k)*c2dtts*
     +                         o2_rem(i,j) ) )))
            else
               R_Olim(i,k,j) = 1.0    ! full oxic remineralisation possible of OM
            endif
 
            
*   6.5  --  Now, do the same for nitrate and denitrification

            if (-op_den(i,k,j).gt.0.0) then
               R_Nlim(i,k,j) = min(1.0, max(0.0, ( no3 / (
     +                         op_den(i,k,j)*dizt(k)*c2dtts*
     +                         no3_rem(i,j) ) )))
            else
               R_Nlim(i,k,j) = 1.0    ! full suboxic remineralisation possible of OM
            endif
            
            
*   6.6  --  If R_Olim and R_Nlim are both 1.0, then no OM is required
*            to be passed to the next depth, and the total OM can be
*            simply calculated as...

            if (R_Olim(i,k,j).eq.1.0 .and. R_Nlim(i,k,j).eq.1.0) then
               op_tot(i,k,j) = op_rem(i,k,j) + op_den(i,k,j)

               ! If all OM from general phytoplankton can be remineralised, 
               ! then the OM from nitrogen fixers can now be calculated as...
            
               if (fix) then

                  oxy = max(0.0, tb(i,k,n_oxy))
                  no3 = max(0.0, tb(i,k,n_no3))
                 
                  oxy=oxy-op_rem(i,k,j)*dizt(k)*c2dtts*o2_rem(i,j)
                  no3=no3-op_den(i,k,j)*dizt(k)*c2dtts*no3_rem(i,j)
                  
                  if (-op_rem_f(i,k,j).gt.0.0) then
                     R_Olim_f(i,k,j) = min(1.0, max(0.0, ( 
     +                                 (oxy-0.5) / (op_rem_f(i,k,j)*
     +                                 dizt(k)*c2dtts*o2_rem_f) )))
                  else
                     R_Olim_f(i,k,j) = 1.0    ! full oxic remineralisation possible of OM
                  endif
                  if (-op_den_f(i,k,j).gt.0.0) then
                     R_Nlim_f(i,k,j) = min(1.0, ( no3 / ( 
     +                                 op_den_f(i,k,j)*dizt(k)*
     +                                 c2dtts*no3_rem_f  )))
                  else
                     R_Nlim_f(i,k,j) = 1.0    ! full suboxic remineralisation possible of OM
                  endif
            
                  if (R_Olim_f(i,k,j).eq.1.0 .and. 
     +                R_Nlim_f(i,k,j).eq.1.0) then
                      op_tot_f(i,k,j) = op_rem_f(i,k,j)+op_den_f(i,k,j)
                  else
                     if (k.lt.kmt(i,j)) then
                        !1
                        op_rem_f(i,k,j) =op_rem_f(i,k,j)*R_Olim_f(i,k,j)
                        op_den_f(i,k,j) =op_den_f(i,k,j)*R_Nlim_f(i,k,j)
                        !2
                        op_unrem_f(i,k,j) = op_tot_f(i,k,j) -
     +                                 (op_rem_f(i,k,j)+op_den_f(i,k,j))
                        !3
                        op_tot_f(i,k,j) =op_rem_f(i,k,j)+op_den_f(i,k,j)
                        op_sed_f(i,k+1,j) = op_sed_f(i,k+1,j)+
     +                                         op_unrem_f(i,k,j)*
     +                                         sgb_frac(i,j,k+1)
                        op_tot_f(i,k+1,j) = op_tot_f(i,k+1,j)+
     +                                      op_unrem_f(i,k,j)*
     +                                      (1.0-sgb_frac(i,j,k+1))
                     else
                        if (sedfluxes) then
                           !1
                           op_rem_f(i,k,j) = op_rem_f(i,k,j)*
     +                                       R_Olim_f(i,k,j)
                           op_den_f(i,k,j) = op_den_f(i,k,j)*
     +                                       R_Nlim_f(i,k,j)
                           !2
                           op_unrem_f(i,k,j) = op_tot_f(i,k,j) -
     +                                         (op_rem_f(i,k,j) + 
     +                                          op_den_f(i,k,j))
                           !3
                           op_tot_f(i,k,j) = op_rem_f(i,k,j) +
     +                                       op_den_f(i,k,j)
                           op_sed_f(i,k,j) = op_sed_f(i,k,j) +
     +                                       op_unrem_f(i,k,j)
                        else
                           op_tot_f(i,k,j) = op_rem_f(i,k,j) +
     +                                       op_den_f(i,k,j)
                        endif ! --> sedfluxes

                     endif ! --> k.lt.kmt(i,j)
                  endif ! --> limitation factors for remin of f-org
               
               else
                  
                  op_rem_f(i,k,j) = 0.0
                  op_den_f(i,k,j) = 0.0
                  op_unrem_f(i,k,j) = 0.0
                  op_tot_f(i,k,j) = 0.0
                  op_sed_f(i,k,j) = 0.0
                
               endif ! --> nitrogen fixation on


*   6.7  --  However, if this is not the case, then we must do a few
*            things:
*               1. alter the OM of the current box so that neither oxygen
*                  nor nitrate is removed past zero
*               2. calculate the total OM that is not remineralised by
*                  oxic and suboxic processes
*               3. if this is not the deepest box at the given lat and
*                  lon, shift this quantity of unremineralised OM into the
*                  next box. Otherwise, remineralise everything            
            
            else
               
               if (k.lt.kmt(i,j)) then
                 
                  !1
                  op_rem(i,k,j) = op_rem(i,k,j) * R_Olim(i,k,j)
                  op_den(i,k,j) = op_den(i,k,j) * R_Nlim(i,k,j)

                  !2
                  op_unrem(i,k,j) = op_tot(i,k,j) -
     +                              (op_rem(i,k,j)+op_den(i,k,j))
                  
                  !3
                  op_tot(i,k,j) = op_rem(i,k,j) + op_den(i,k,j)
                  op_sed(i,k+1,j) = op_sed(i,k+1,j)+
     +                              op_unrem(i,k,j)*sgb_frac(i,j,k+1)
                  op_tot(i,k+1,j) = op_tot(i,k+1,j)+op_unrem(i,k,j)
     +                              *(1.-sgb_frac(i,j,k+1))
               else
                  if (sedfluxes) then
                     !1
                     op_rem(i,k,j) = op_rem(i,k,j) * R_Olim(i,k,j)
                     op_den(i,k,j) = op_den(i,k,j) * R_Nlim(i,k,j)

                     !2
                     op_unrem(i,k,j) = op_tot(i,k,j) -
     +                                 (op_rem(i,k,j)+op_den(i,k,j))
                     
                     !3
                     op_tot(i,k,j) = op_rem(i,k,j) + op_den(i,k,j)
                     op_sed(i,k,j) = op_sed(i,k,j) + op_unrem(i,k,j)
                  else 
                     op_tot(i,k,j) = op_rem(i,k,j) + op_den(i,k,j)
                  endif ! --> sedfluxes
               endif ! --> k.lt.kmt(i,j)

            
               if (fix) then
               
                  if (k.lt.kmt(i,j)) then
                     R_Olim_f(i,k,j) = 0.0
                     R_Nlim_f(i,k,j) = 0.0
               
                     !1
                     op_rem_f(i,k,j) = op_rem_f(i,k,j) * R_Olim_f(i,k,j)
                     op_den_f(i,k,j) = op_den_f(i,k,j) * R_Nlim_f(i,k,j)

                     !2
                     op_unrem_f(i,k,j) = op_tot_f(i,k,j) -
     +                                 (op_rem_f(i,k,j)+op_den_f(i,k,j))
               
                     !3
                     op_tot_f(i,k,j) = op_rem_f(i,k,j) + op_den_f(i,k,j)
                     op_sed_f(i,k+1,j) = op_sed_f(i,k+1,j)+
     +                                   op_unrem_f(i,k,j)*
     +                                   sgb_frac(i,j,k+1)
                     op_tot_f(i,k+1,j) = op_tot_f(i,k+1,j)+
     +                                   op_unrem_f(i,k,j)*
     +                                   (1.0-sgb_frac(i,j,k+1))
                  else
                    if (sedfluxes) then
                       !1
                       op_rem_f(i,k,j) = op_rem_f(i,k,j)*R_Olim_f(i,k,j)
                       op_den_f(i,k,j) = op_den_f(i,k,j)*R_Nlim_f(i,k,j)
                       !2
                       op_unrem_f(i,k,j) = op_tot_f(i,k,j) -
     +                                 (op_rem_f(i,k,j)+op_den_f(i,k,j))
                       !3
                       op_tot_f(i,k,j) = op_rem_f(i,k,j)+op_den_f(i,k,j)
                       op_sed_f(i,k,j)=op_sed_f(i,k,j)+op_unrem_f(i,k,j)

                    else
                       op_tot_f(i,k,j) = op_rem_f(i,k,j)+op_den_f(i,k,j)
                    endif ! --> sedfluxes
                  endif ! --> k.lt.kmt(i,j)
               
               else
                  
                  op_rem_f(i,k,j) = 0.0
                  op_den_f(i,k,j) = 0.0
                  op_unrem_f(i,k,j) = 0.0
                  op_tot_f(i,k,j) = 0.0
                  op_sed_f(i,k,j) = 0.0
                
               endif ! --> nitrogen fixation turned on
                
            endif ! oxygen and nitrogen limitation factors

         enddo ! --> loop over k
      endif ! --> kmt(i,j).ne.0
      enddo ! --> loop over i

      else  ! --> ReminRelax

         do i = istrt,iend
         if (kmt(i,j).ne.0) then
            do k = 1,kmt(i,j)
        
               if (k.lt.kmt(i,j)) then
               
                  op_tot(i,k,j) = pop(i,k,j) * (1.0-sgb_frac(i,j,k)) 
                  op_sed(i,k,j) = pop(i,k,j) - op_tot(i,k,j) 

                  if (fix) then
                     op_tot_f(i,k,j) = pop_fix(i,k,j)*
     +                                 (1.0-sgb_frac(i,j,k))
                     op_sed_f(i,k,j) = pop_fix(i,k,j) - op_tot_f(i,k,j)
                  else
                     op_tot_f(i,k,j) = 0.0
                     op_sed_f(i,k,j) = 0.0
                  endif

               else

                  op_tot(i,k,j) = pop(i,k,j) * (1.0-sgb_frac(i,j,k)) 
                  op_sed(i,k,j) = op_sed(i,k,j) + 
     +                            (pop(i,k,j) - op_tot(i,k,j))

                  if (fix) then
                     op_tot_f(i,k,j) = pop_fix(i,k,j) * 
     +                                 (1.0-sgb_frac(i,j,k))
                     op_sed_f(i,k,j) = op_sed_f(i,k,j) + 
     +                                 (pop_fix(i,k,j)-op_tot_f(i,k,j))
                  else
                     op_tot_f(i,k,j) = 0.0
                     op_sed_f(i,k,j) = 0.0
                  endif

               endif ! --> k.lt.kmt(i,j)


               ! Get strength of denitrification

               oxy = max(0.0, tb(i,k,n_oxy))
      
               if (n_no3.gt.0 .and. den) then

                  no3 = max(0.0, tb(i,k,n_no3))


                  if (oxy.lt.Den_o2lim) then
                        
                     if (Den_type.eq."Linear ") then
                        R_den(i,k,j) = min(1., max(0.,
     +                                 (Den_o2lim - oxy)/Den_o2lim ))
                     endif
                           
                     if (Den_type.eq."QuadPos") then
                        R_den(i,k,j) = min(1.0, max(0.0,
     +                                ( ((Den_o2lim - oxy)/
     +                                    Den_o2lim)**2. ) ) )
                     endif

                     if (Den_type.eq."QuadNeg") then
                        R_den(i,k,j) = min(1.0, max(0.0,
     +                                ( -(oxy/Den_o2lim)**2 + 1.0)))
                     endif
                       
                     if (Den_type.eq."Sigmoid") then
                        R_den(i,k,j) = min(1.0, max(0.0,
     +                                (1./(1.0 - exp(-Den_o2lim/2.0) +
     +                                 exp(oxy-(Den_o2lim/2.)) )) ))
                     endif
                        
                  else
                     R_den(i,k,j) = 0.0
                  endif
               
                  if (no3.le.Den_nlim .and. R_den(i,k,j).gt.0.0) then
                     denrelax = .5 + .5*tanh(0.25*no3-Den_nlim*0.25-2.5)
                  else
                     denrelax = 1.0
                  endif
            
                  if (denrelax.lt.R_den(i,k,j)) then
                     R_den(i,k,j) = denrelax
                  endif
          
                  
               else ! if den is FALSE
                  R_den(i,k,j) = 0.0
               endif
               
               op_rem(i,k,j) = op_tot(i,k,j)*(1.0-R_den(i,k,j))
               op_den(i,k,j) = op_tot(i,k,j)*R_den(i,k,j)

               if (fix) then
                  op_rem_f(i,k,j) = op_tot_f(i,k,j)*(1.0-R_den(i,k,j))
                  op_den_f(i,k,j) = op_tot_f(i,k,j)*R_den(i,k,j)
               else
                  op_rem_f(i,k,j) = 0.0
                  op_den_f(i,k,j) = 0.0
               endif

            enddo ! --> k loop
         endif ! --> kmt(i,j).gt.0
         enddo ! --> i loop
      
      endif  ! --> ReminRelax
      endif  ! --> n.eq.n_oxy (i.e. that bgc is on)
     

*   6.8  --  Check conservation of matter
            
      if (n.eq.3) then
         do i = istrt,iend
         if (kmt(i,j).ne.0) then
            tot_b = 0.0
            do k = 2,kmt(i,j)

               tot_b = tot_b + op_tot(i,k,j) + op_sed(i,k,j) +
     +                 op_tot_f(i,k,j) + op_sed_f(i,k,j)

            enddo

            tot_a = op_tot(i,1,j) + op_tot_f(i,1,j)

            dif = ABS(tot_b)-tot_a
            tot = ABS(tot_b)+tot_a

            if (ABS(dif).gt.epsilon(tot)) then
               print*, " "
               print*, "   REMINERALISATION SCHEME NOT         "
               print*, "       CONSERVING MATTER                  "
               print*, " "
               print*," Org P in total m-2 s-1 = ",op_tot(i,1,j)
               print*,"F-Org P in total m-2 s-1 = ",op_tot_f(i,1,j)
               print*," sum of reactions below surface = ",tot_b
               print*, " "
               print*, 'i = ',i,' j = ',j
               print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
               print*, " "
               print*, " "
               print*, "        STOPPING MODEL                   "
               print*, " "
               print*, " "
               stop
            endif


            tot = 0.0
            dif = 0.0
            do k = 1,kmt(i,j)

            ! Check that all organic matter is accounted for at each k
      
            tot_b = op_tot(i,k,j) + op_tot_f(i,k,j)
            tot_a = op_rem(i,k,j) + op_rem_f(i,k,j) + 
     +              op_den(i,k,j) + op_den_f(i,k,j)
               
            dif = tot_b-tot_a
            tot = tot_b+tot_a
            if (ABS(dif).gt.epsilon(tot)) then
               print*, " "
               print*, "   PELAGIC REMINERALISATION SCHEME NOT         "
               print*, "            CONSERVING MATTER                  "
               print*, " "
               print*," Org P in total m-2 s-1 = ",tot_b
               print*," Org P sum of parts m-2 s-1 = ",tot_a
               print*, " "
               print*, 'i = ',i,' j = ',j,'k = ',k
               print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
               print*, " "
               print*, " no3 = ", no3, " oxy = ", oxy
               print*, " optot =", op_tot(i,k,j),
     +                 " optot_f = ",op_tot_f(i,k,j)
               print*, " op_rem =", op_rem(i,k,j),
     +                 " op_den = ", op_den(i,k,j)
               print*, " op_rem_f =", op_rem_f(i,k,j),
     +                 " op_den_f = ", op_den_f(i,k,j)
               print*, " "
               print*, "              STOPPING MODEL                   "
               print*, " "
               print*, " "
               stop
            endif

            enddo
         endif
         enddo
      endif

      if (N15diagnostics) then
      do i = istrt,iend
         if (kmt(i,j).ne.0) then
            do k = 1,km
               if (n.eq.n_n15 .and. j.eq.latj .and. i.eq.loni) then
               print*, " "
               print*, "Pelagic Org", "oxic rem", "suboxic rem"
               print*, op_tot(i,k,j),op_rem(i,k,j),op_den(i,k,j)
               print*, " "
               print*, "Pelagic F-Org", "oxic rem", "suboxic rem"
               print*, op_tot_f(i,k,j),op_rem_f(i,k,j),op_den_f(i,k,j)
               print*, " "
               print*, "Sediment Org", "Sediment F-Org"
               print*, op_sed(i,k,j),op_sed_f(i,k,j)
               print*, " "
               endif
            enddo
         endif
      enddo
      endif

* ****  END STEP 6  **** *

      
***********************************************************
* ****  STEP 7 -->  DO SEDIMENT PROCESSES  **** *
***********************************************************

      if (n.eq.3) then  ! --> make sure calcs only happen once
      if (sedfluxes) then

*   7.1  --  Calculate and save flux of organics to sediment

      ! Here, we calculate what happens to the organics in sediments 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c      !!!   NOTE :: please read if subgrid==.TRUE.   !!!              c
c                                                                      c
c                                                                      c
c   The implementation of subgrid scale bathymetry will increase       c
c   the amount of organic matter that is treated via sedimentary       c
c   processes.                                                         c
c                                                                      c
c   In Step 3, the remineralisation curve was altered using the        c
c   subgrid scale bathymetry fraction of sediments for each water      c
c   column.                                                            c
c                                                                      c
c   In step 6 above, the remineralisation of organic matter            c
c   through the water column was recalculated to account for the       c
c   conservation of oxygen and nitrate concentrations. If O2 or        c
c   NO3 are not available, then oxic and suboxic remineralisation      c
c   don't take place. The un-remineralised organic matter is           c
c   passed to the next depth.                                          c
c                                                                      c
c   A fraction of the unremineralised organic matter that is passed    c
c   to the next depth is allocated to the sediments, while the rest    c
c   is allocated to the water column.                                  c
c                                                                      c
c   The introduction of additional sediment at shallower depths        c
c   therefore requires remineralisation to occur in some places        c
c   where it otherwise would not (i.e. where O2 and NO3 are not        c
c   available.                                                         c
c                                                                      c
c   To account for water columns where O2 and NO3 are limiting,        c
c   but where organic matter is hitting the sediments and must         c
c   be remineralised, we assume this process occurs via sulfate        c
c   reduction:                                                         c
c                                                                      c
c    59 H2SO4 + Org --> 106 CO2 + 16 NH3 + H3PO4 + 59 H2S + 62 H2O     c
c                                                                      c
c   However, because SO4 is not YET explicitly simulated, the          c
c   following calculations are made:                                   c
c      1. We find the rain rate of organics hitting the sediment       c
c      2. We calculate the amount of O2 and NO3 required to            c
c         remineralise the organics via the Bohlen parameterisation    c
c      3. We find the amount of O2 and NO3 that are available,         c
c         and determine if these concentrations are limiting           c
c      4. Finally, if both O2 and NO3 are deemed to be limiting,       c
c         we remineralise all remaining organic matter without         c
c         removing O2 or NO3. This step assumes that the organics      c
c         are remineralised by sulfate reduction.                      c
c      5. We save the amount of organics remineralised via             c
c         sulfate reduction to an output array, so that it can         c
c         be compared to estimates in the literature that suggest      c
c         roughly 18% of total sediment remineralisation.              c
c                                                                      c
c                                                                      c
c    The following code executes these steps.                          c
c                                                                      c
c                                                                      c
c    Author: Pearse J Buchanan        July 2017                        c
c    email: pearse.buchanan@utas.edu.au  /                             c
c           pearse.buchanan@hotmail.com                                c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      do i = istrt,iend
         if (kmt(i,j).gt.0) then
         do k = 1,kmt(i,j)


*   7.1  --  Find the rain rate of organics hitting sediments

            ! We have already saved the sedimentary rain of organics for
            ! all depth levels in the steps above. This information is  
            ! held in the op_sed, op_sed_f, ic_sed and si_sed arrays.

*   7.2  --  Calculate sedimentary remineralisation via either oxic or
*             suboxic means (i.e. sedimentary denitrification)
*
*                   Bohlen2012 - parameterised function based on rain
*                                rate of organic carbon, and the ambient
*                                oxygen and nitrate concentrations
*          
            if (n_no3.gt.0) then ! requires no3 cycle to be on
            if (op_sed(i,k,j).lt.0.0 .or. op_sed_f(i,k,j).lt.0.0) then ! Org must be hitting sediment 
            
               oxy = max(0.0, tb(i,k,n_oxy)) 
               no3 = max(0.0, tb(i,k,n_no3))

               o_n = oxy-no3

               ! positive units of mmol N m-3 per leap timestep using
               ! Bohlen parameterisation, where units of N loss is found
               ! based on amount of org C hitting sediments
               sdenTg = (0.04+Boh_a*0.98**o_n) * 
     +                 (-op_sed(i,k,j)*carb2P(i,j)*c2dtts*dizt(k)) ! 0.0 if no org
               sdenTf = (0.04+Boh_a*0.98**o_n) * 
     +                  (-op_sed_f(i,k,j)*CP_fix*c2dtts*dizt(k)) ! 0.0 if no org

               ! Accelerate denitrification in upper 250 metres
               if (sgb_frac(i,j,k).gt.0.0 .and. zt(k).lt.25000) then
                  sdenTg = min((1.0-sgb_frac(i,j,k))*Sden_acl+1.0,20.0)*
     +                     sdenTg             
                  sdenTf = min((1.0-sgb_frac(i,j,k))*Sden_acl+1.0,20.0)*
     +                     sdenTf
               endif
                  
               ! Is enough NO3 available to remineralise?
               ! positive units of mmol N m-3 per leap timestep
               if (sdenTg.gt.0.0) then
                  sdenG(i,k,j) = min(1.0, (no3*0.66)/sdenTg)*sdenTg ! remineralise general org first
                  sdenG(i,k,j) = (0.5+0.5*tanh(10.*no3-5.))*sdenG(i,k,j)! apply relaxation at low no3 (< 1) values
               else
                  sdenG(i,k,j) = 0.0
               endif
               if (sdenTf.gt.0.0) then
                  sdenF(i,k,j) = min(1.0, (no3*0.66-sdenG(i,k,j)) /
     +                           sdenTf)*sdenTf
                  sdenF(i,k,j) = (0.5+0.5*tanh(10.*no3-5.))*sdenF(i,k,j)! apply relaxation at low no3 values
               else
                  sdenF(i,k,j) = 0.0
               endif
              
               ! convert back to 1e-2 mmol PO4 m-2 s-1 (negative units) 
               ! and make sure that they are not > op_sed or op_sed_f
               sdenG(i,k,j) = sdenG(i,k,j)/ 
     +                        (no3_rem(i,j)*c2dtts*dizt(k))
               sdenG(i,k,j) = max(op_sed(i,k,j), sdenG(i,k,j))
               sdenF(i,k,j) = sdenF(i,k,j)/
     +                        (no3_rem_f*c2dtts*dizt(k))
               sdenF(i,k,j) = max(op_sed_f(i,k,j), sdenF(i,k,j))
               
               ! Assign remaining organic matter to oxygen
               ! negative units of 1e-2 mmol PO4 m-2 s-1 
               if (sdenG(i,k,j).gt.op_sed(i,k,j)) then
                  sremG(i,k,j) = (1.-sdenG(i,k,j)/op_sed(i,k,j))
     +                           *op_sed(i,k,j)
                  sremG(i,k,j) = min(1.0, (oxy*0.66)/(sremG(i,k,j)*
     +                         o2_rem(i,j)*c2dtts*dizt(k)))*sremG(i,k,j)
                  sremG(i,k,j) = (0.5+0.5*tanh(0.2*oxy-5.))*sremG(i,k,j)! apply relaxation at low oxy values
               else
                  sremG(i,k,j) = 0.0
               endif
               if (sdenF(i,k,j).gt.op_sed_f(i,k,j)) then
                  sremF(i,k,j) = (1.-sdenF(i,k,j)/op_sed_f(i,k,j))
     +                           *op_sed_f(i,k,j)
                  sremF(i,k,j) = min(1.0,
     +                           (oxy*0.66)/((sremF(i,k,j)+sremG(i,k,j))
     +                           *o2_rem_f*c2dtts*dizt(k)))*sremF(i,k,j)
                  sremF(i,k,j) = (0.5+0.5*tanh(0.2*oxy-5.))*sremF(i,k,j)! apply relaxation at low oxy values
               else
                  sremF(i,k,j) = 0.0
               endif
               
               ! Assign remaining organic matter to sulfate
               ! negative units of 1e-2 mmol PO4 m-2 s-1 
               if (op_sed(i,k,j).lt.0.0) then
                  ssulG(i,k,j) = (1.-(sdenG(i,k,j)+sremG(i,k,j))
     +                           /op_sed(i,k,j))*op_sed(i,k,j)
               else
                  ssulG(i,k,j) = 0.0
               endif
               if (op_sed_f(i,k,j).lt.0.0) then
                  ssulF(i,k,j) = (1.-(sdenF(i,k,j)+sremF(i,k,j))
     +                           /op_sed_f(i,k,j))*op_sed_f(i,k,j)
               else
                  ssulF(i,k,j) = 0.0
               endif
            
            ! Check that all organic matter is accounted for at each k
            tot_b = op_sed(i,k,j) + op_sed_f(i,k,j)
            tot_a = sdenG(i,k,j) + sdenF(i,k,j) +
     +              sremG(i,k,j) + sremF(i,k,j) +
     +              ssulG(i,k,j) + ssulF(i,k,j)
               
            dif = tot_b-tot_a
            tot = tot_b+tot_a
            if (ABS(dif).gt.epsilon(tot)) then
               print*, " "
               print*, "   SEDIMENTARY REMINERALISATION SCHEME NOT     "
               print*, "            CONSERVING MATTER                  "
               print*, " "
               print*," Org P arriving m-2 s-1 = ",tot_b
               print*," Org P removed m-2 s-1 = ",tot_a
               print*, " "
               print*, 'i = ',i,' j = ',j,'k = ',k
               print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
               print*, " "
               print*, " no3 = ", no3, " oxy = ", oxy
               print*, " sdenTg =", sdenTg," sdenTf = ", sdenTf
               print*, " sdenG =", sdenG(i,k,j)," sdenF = ",sdenF(i,k,j)
               print*, " sremG =", sremG(i,k,j)," sremF = ",sremF(i,k,j)
               print*, " ssulG =", ssulG(i,k,j)," ssulF = ",ssulF(i,k,j)
               print*, " "
               print*, "              STOPPING MODEL                   "
               print*, " "
               print*, " "
               stop
            endif
      
            else 
         
               sdenG(i,k,j) = 0.0
               sdenF(i,k,j) = 0.0
               sremG(i,k,j) = 0.0
               sremF(i,k,j) = 0.0
               ssulG(i,k,j) = 0.0
               ssulF(i,k,j) = 0.0
               
            endif ! --> if stuff is arriving to seafloor 
            
            ! Calculate sulfate reduction in sediments for output (units
            ! of mmol C s-1)
            sulS(i,k,j) = ( -ssulG(i,k,j)*carb2P(i,j)
     +                      -ssulF(i,k,J)*CP_fix )*dizt(k)
            
            if (N15diagnostics .and. j.eq.latj .and. i.eq.loni) then
            print*, " "
            print*, " "
            print*, " CHECK SEDIMENT FLUX FIELDS "
            print*,"depth=",k,"lat=",j*1.59,"lon=",i*2.8125
            print*,op_rem(i,k,j), op_den(i,k,j)
            print*,sremG(i,k,j), sdenG(i,k,j), ssulG(i,k,j)
            print*,op_rem_f(i,k,j), op_den_f(i,k,j)
            print*,sremF(i,k,j), sdenF(i,k,j), ssulF(i,k,j)
            endif
           

            else ! --> N cycle is not on
            
            if (op_sed(i,k,j).lt.0.0) then ! Org must be hitting sediment 
            
               oxy = max(0.0, tb(i,k,n_oxy))
               
               sdenG(i,k,j) = 0.0 
               sdenF(i,k,j) = 0.0 
               sremF(i,k,j) = 0.0 
               ssulF(i,k,j) = 0.0 

               ! Assign organic matter to oxygen
               ! negative units of 1e-2 mmol PO4 m-2 s-1 
               if (op_sed(i,k,j).lt.0.0) then
                  sremG(i,k,j) = (1.-sdenG(i,k,j)/op_sed(i,k,j))
     +                           *op_sed(i,k,j)
                  sremG(i,k,j) = min(1.0, (oxy*0.66)/(sremG(i,k,j)*
     +                         o2_rem(i,j)*c2dtts*dizt(k)))*sremG(i,k,j)
                  sremG(i,k,j) = (0.5+0.5*tanh(0.2*oxy-5.))*sremG(i,k,j)! apply relaxation at low oxy values
               else
                  sremG(i,k,j) = 0.0
               endif
               
               ! Assign remaining organic matter to sulfate
               ! negative units of 1e-2 mmol PO4 m-2 s-1 
               if (op_sed(i,k,j).lt.0.0) then
                  ssulG(i,k,j) = (1.-(sdenG(i,k,j)+sremG(i,k,j))
     +                           /op_sed(i,k,j))*op_sed(i,k,j)
               else
                  ssulG(i,k,j) = 0.0
               endif
               
            ! Check that all organic matter is accounted for at each k
            tot_b = op_sed(i,k,j)
            tot_a = sdenG(i,k,j) + sremG(i,k,j) + ssulG(i,k,j) 
               
            dif = tot_b-tot_a
            tot = tot_b+tot_a
            if (ABS(dif).gt.epsilon(tot)) then
               print*, " "
               print*, "   SEDIMENTARY REMINERALISATION SCHEME NOT     "
               print*, "            CONSERVING MATTER                  "
               print*, " "
               print*," Org P arriving m-2 s-1 = ",tot_b
               print*," Org P removed m-2 s-1 = ",tot_a
               print*, " "
               print*, 'i = ',i,' j = ',j,'k = ',k
               print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
               print*, " "
               print*, " oxy = ", oxy
               print*, " sdenTg =", sdenTg
               print*, " sdenG =", sdenG(i,k,j)
               print*, " sremG =", sremG(i,k,j)
               print*, " ssulG =", ssulG(i,k,j)
               print*, " "
               print*, "              STOPPING MODEL                   "
               print*, " "
               print*, " "
               stop
            endif
      
            else ! --> if stuff isn't arriving to sediment 
         
               sdenG(i,k,j) = 0.0
               sdenF(i,k,j) = 0.0
               sremG(i,k,j) = 0.0
               sremF(i,k,j) = 0.0
               ssulG(i,k,j) = 0.0
               ssulF(i,k,j) = 0.0
               
            endif ! --> if stuff is arriving to seafloor 
               
            ! Calculate sulfate reduction in sediments for output
            sulS(i,k,j) = ( -ssulG(i,k,j)*carb2P(i,j)
     +                      -ssulF(i,k,J)*CP_fix )*dizt(k)

            endif ! --> N cycle to be on (n_no3.gt.0)
          
            enddo ! --> k loop
            endif ! --> kmt(i,j).gt.0
         enddo ! --> i loop
            
      else  ! --> if sedfluxes is FALSE 
         
         do i = istrt,iend
            if (kmt(i,j).ne.0) then
            do k = 1,km
               sdenG(i,k,j) = 0.0
               sdenF(i,k,j) = 0.0
               sremG(i,k,j) = 0.0
               sremF(i,k,j) = 0.0
               ssulG(i,k,j) = 0.0
               ssulF(i,k,j) = 0.0
      
            enddo
            endif
         enddo
      
      endif  ! --> sedfluxes 

      endif  ! --> n.eq.3

* ****  END STEP 7  **** *


      
***********************************************************
* ****  STEP 8 -->  Nitrogen Isotopic Fractionation  **** *
***********************************************************

      if (n.eq.n_n15) then
        do i = istrt,iend
        if (kmt(i,j).ne.0) then
        
        ! Initialise arrays to zero before iteration over k
        tot = 0.0
        dif = 0.0

        do k = 1,km
           pon_cons(k) = 0.0
           ponf_cons(k) = 0.0
           on_rem(k) = 0.0
           on_rem_f(k) = 0.0
           no3_den(k) = 0.0
           no3_den_f(k) = 0.0
           pon15(i,k,j) = 0.0
           pon14(i,k,j) = 0.0
           ponf15(i,k,j) = 0.0
           ponf14(i,k,j) = 0.0
        enddo
        
        do k = 1,km
           
           no3 = tb(i,k,n_no3)
           n15 = tb(i,k,n_n15)
        
           if (no3.lt.0.0 .and. n15.lt.0.0) then
              no3 = ABS(tb(i,k,n_no3))
              n15 = ABS(tb(i,k,n_n15))
           endif
        
           n14 = no3-n15


*   BEGIN FRACTIONATION OF NITROGEN *

c
c   Routine to fractionate 15N during processes affecting NO3
c        
c       INPUTS:
c           i           =   longitudinal index for grid
c           j           =   latitudinal index for grid
c           k           =   depth index for grid
c           no3         =   nitrate concentration at grid box
c           n15         =   15nitrate concentration at grid box 
c           n14         =   14nitrate concentration at grid box 
c           N15diagnostics  =   Boolean for output to screen for 15N 
c           loni / lonj      =   Write output for 15N at this coord 
c       OUTPUTS:
c           pon15(i,k,j)    =   amount of 15NO3 in Org (mmol/m3)
c           pon14(i,k,j)    =   amount of 14NO3 in Org (mmol/m3)
c           pon_cons(k)     =   amount of Org consumed as depth inc
c           ponf15(i,k,j)   =   amount of 15NO3 in F-Org (mmol/m3)
c           ponf14(i,k,j)   =   amount of 14NO3 in F-Org (mmol/m3)
c           ponf_cons(k)    =   amount of F-Org consumed as depth inc
c           d_n15(i,k,j)    =   Change in 15NO3 due to assim/remin Org
c           d_n15d(i,k,j)   =   Change in 15NO3 due to denitrif Org
c           d_n15f(i,k,j)   =   Change in 15NO3 due to assim/remin F-Org
c           d_n15df(i,k,j)  =   Change in 15NO3 due to denitrif F-Org
c           d15org(i,j)     =   The del15N value or organics at sediment
c                    
c
c        This routine fractionates the nitrogen isotopes held within
c        nitrate during assimilation by phytoplankton, denitrification
c        in the water column, and denitrification in the sediments.
c        Nitrogen fixation and atmospheric deposition of reactive
c        nitrogen introduce nitrate at set isotopic signatures, and
c        therefore do not fractionate.
c
c        The fractionation occurs stages 4 and 5 of each section of the
c        following code, where the epsilon fractionation factor (Exxx) is
c        applied against the isotopic signature of the nitrate in the
c        water column (IUC), using the accumulated product equation (see
c        Altabet & Francois, 2001, Deep Sea Research).
c
c        AUTHOR:
c           Pearse James Buchanan          
c           pearse.buchanan@utas.edu.au
c           pearse.buchanan@hotmail.com
c           pearse.buchanan@princeton.edu
c           
c
c......................................................................c


*   8.1  --  Calculate fractionation due to assimilation
        if (k.eq.1) then
               
           if (op_tot(i,k,j).gt.0.0) then ! If organic matter is created 

              ! 1. get NO3 uptake
              on_rem(k) = op_tot(i,k,j)*NtoP(i,j)*c2dtts*dizt(k)
              ! 2. get utilisation
              NUC = min(0.999, max(0.001, on_rem(k)/no3 ))
              ! 3. get n15:n14 ratio
              IUC = min(2.0, max(0.5, n15/n14 ))
              ! 4. calculate fractionation factor
              f_ass = IUC+Eass*(1-NUC)/NUC*log(1.0-NUC)*1e-3*IUC
              ! 5. get change in 15NO3 
              d_n15(i,k,j) = (1. - (1./(f_ass+1.)))*on_rem(k)
                     
              ! 6. store the assimilation 15N in organic matter
              pon15(i,k,j) = d_n15(i,k,j)
              pon14(i,k,j) = on_rem(k) - pon15(i,k,j)
              pon_cons(k) = 0.0
              
           else
              
              on_rem(k) = 0.0
              d_n15(i,k,j) = 0.0
              pon14(i,k,j) = 0.0
              pon15(i,k,j) = 0.0
              pon_cons(k) = 0.0

           endif

           ! 5. calculate storage of 15N due to nitrogen fixation
           if (op_tot_f(i,k,j).gt.0.0) then

              ! N2 fixation introduces NO3 with 15N:14N of 0.999
              on_rem_f(k) = op_tot_f(i,k,j)*NP_fix*c2dtts*dizt(k)
              ponf15(i,k,j) = (1.-(1./(1.+(1.+Efix*1e-3))))*on_rem_f(k)
              ponf14(i,k,j) = on_rem_f(k) - ponf15(i,k,j)
              ponf_cons(k) = 0.0
              
              ! set d_n15f to zero to account for no 15N uptake from water
              d_n15f(i,k,j) = 0.0
              
           else  ! if no N-fixation - no storage of 15N in organics

              on_rem_f(k) = 0.0
              ponf15(i,k,j) = 0.0
              ponf14(i,k,j) = 0.0
              ponf_cons(k) = 0.0
              d_n15f(i,k,j) = 0.0
               
           endif

        endif ! k.eq.1     
              

*   8.2  --  Calculate below surface processes:
*            OXIC REMINERALISATION = release of 15N from organic matter              
*            PELAGIC DENITRIFICATION = consumption of 15N from water column
*            SEDIMENTARY PROCESSES = release/consumption of 15N from sediment


        if (k.gt.1) then   
          
           ! track changes in isotopes as depth increases
           pon15(i,k,j) = pon15(i,k-1,j)*(dizt(k)/dizt(k-1))
           pon14(i,k,j) = pon14(i,k-1,j)*(dizt(k)/dizt(k-1))
              
           ! OXIC AND SUBOXIC REMINERALISATION
           ! get the release of NO3 due to remineralisation
           on_rem(k) = (op_tot(i,k,j) + op_sed(i,k,j))
     +                 *dizt(k)*c2dtts*NtoP(i,j)
           
           ! track change in consumption as depth increases
           pon_cons(k) = pon_cons(k-1) + 
     +                   (-on_rem(k)*(dizt(1)/dizt(k)))
        
           if (on_rem(k).lt.0.0) then ! if stuff was remineralised
        
              if (k.lt.kmt(i,j)) then
                 ! calculate the release of 15N (value negative)
                 d_n15(i,k,j) = (1. - (1./(pon15(i,k,j) / 
     +                           pon14(i,k,j)+1.)))*on_rem(k)
                 
                 ! calculate change in pon15 and pon14 due to frac
                 pon15(i,k,j)=pon15(i,k,j)-(-d_n15(i,k,j))
                 pon14(i,k,j)=pon14(i,k,j)-(-on_rem(k)-(-d_n15(i,k,j)))
              else
                 ! release all remaining 15N (value negative) 
                 d_n15(i,k,j) = -pon15(i,k,j)
                 pon15(i,k,j) = 0.0
                 pon14(i,k,j) = 0.0
              endif
           
           else
              
              on_rem(k) = 0.0
              d_n15(i,k,j) = 0.0
           
           endif

               
           ! DENITRIFICATION
           if (op_den(i,k,j).lt.0.0 .or. sdenG(i,k,j).lt.0.0) then
              
              ! get the consumption of NO3 due to suboxic remineralisation
              
              ! 1. get NO3 uptake
              tot = op_den(i,k,j) + sdenG(i,k,j) 
              no3_den(k) = tot*dizt(k)*c2dtts*no3_rem(i,j)
              ratio = op_den(i,k,j)/tot
              ! 2. get utilisation
              NUCpel = min(0.999, max(0.001,
     +                 no3_den(k)*ratio/no3 ))
              NUCsed = min(0.999, max(0.001,
     +                 no3_den(k)*(1.0-ratio)/no3 ))
              ! 3. get n15:n14 ratio
              IUC = min(2.0, max(0.5, n15/n14 ))
              ! 4. calculate fractionation factors
              f_pel =IUC+Ewc*(1-NUCpel)/NUCpel*log(1.0-NUCpel)*IUC*1e-3
              f_sed =IUC+Esed*(1-NUCsed)/NUCsed*log(1.0-NUCsed)*IUC*1e-3
              ! 5. get change in 15NO3 with fractionation
              d_n15d(i,k,j) = (1.-(1./(f_pel+1.)))*no3_den(k)*ratio +
     +                      (1.-(1./(f_sed+1.)))*no3_den(k)*(1.0-ratio)

           else
             
              tot = 0.0
              no3_den(k) = 0.0
              d_n15d(i,k,j) = 0.0
              ratio = 0.0
            
           endif ! --> op_den(i,k,j).lt.0.0

           
           ! 2. Do the same for N-fixation
           
           ! get the release of NO3 due to remineralisation
           on_rem_f(k) = (op_tot_f(i,k,j) + op_sed_f(i,k,j)) *
     +                    dizt(k)*c2dtts*NP_fix
              
           ! track changes in organics and isotopes as depth increases
           ponf_cons(k) = ponf_cons(k-1) + 
     +                    (-on_rem_f(k)*(dizt(1)/dizt(k)))
           ponf15(i,k,j) = ponf15(i,k-1,j)*(dizt(k)/dizt(k-1))
           ponf14(i,k,j) = ponf14(i,k-1,j)*(dizt(k)/dizt(k-1))
         
           if (on_rem_f(k).lt.0.0) then
          
              if (k.lt.kmt(i,j)) then
                 ! calculate the release of 15N (value negative)
                 d_n15f(i,k,j) = (1.-(1./(ponf15(i,k,j)/
     +                           ponf14(i,k,j)+1.)))*on_rem_f(k)

                 ! calculate change in pon15 and pon14 due to frac
                 ponf15(i,k,j) = ponf15(i,k,j) - (-d_n15f(i,k,j))
                 ponf14(i,k,j) = ponf14(i,k,j) - 
     +                           (-on_rem_f(k)-(-d_n15f(i,k,j)))
              else
                 ! release all remaining 15N (value negative) 
                 d_n15f(i,k,j) = -ponf15(i,k,j)
                 ponf15(i,k,j) = 0.0
                 ponf14(i,k,j) = 0.0
              endif

           else
             
              on_rem_f(k) = 0.0
              d_n15f(i,k,j) = 0.0

           endif ! --> on_rem_f.lt.0.0


           ! DENITRIFICATION
           if (op_den_f(i,k,j).lt.0.0 .or. sdenF(i,k,j).lt.0.0) then
              
              ! get the consumption of NO3 due to suboxic remineralisation
              
              ! 1. get NO3 uptake
              tot = op_den_f(i,k,j) + sdenF(i,k,j) 
              no3_den_f(k) = tot*dizt(k)*c2dtts*no3_rem_f ! no3_den_f now positive
              ratio = op_den_f(i,k,j)/tot
              ! 2. get utilisation
              NUCpel = min(0.999, max(0.001, (no3_den(k)+no3_den_f(k))
     +                 *ratio/no3 ))
              NUCsed = min(0.999, max(0.001, (no3_den(k)+no3_den_f(k))
     +                 *(1.0-ratio)/no3 ))
              ! 3. get n15:n14 ratio
              IUC = min(2.0, max(0.5, n15/n14 ))
              ! 4. calculate fractionation factors
              f_pel =IUC+Ewc*(1-NUCpel)/NUCpel*log(1.0-NUCpel)*1e-3*IUC
              f_sed =IUC+Esed*(1-NUCsed)/NUCsed*log(1.0-NUCsed)*1e-3*IUC
              ! 5. get change in 15NO3 with fractionation
              d_n15df(i,k,j) = (1.-(1./(f_pel+1.)))*no3_den_f(k)*ratio +
     +                     (1.-(1./(f_sed+1.)))*no3_den_f(k)*(1.0-ratio)
              
           else
            
              tot = 0.0
              no3_den_f(k) = 0.0
              d_n15df(i,k,j) = 0.0
              ratio = 0.0

           endif ! --> op_den(i,k,j).lt.0.0



**   8.3  --  Undertake important checks of conservation
*
*           ! Check that on_rem = release of 15N + 14N 
*           dif=0.0
*           tot=0.0
*           dif = on_rem(k) + (
*     +           ( pon15(i,k-1,j)*(dizt(k)/dizt(k-1)) - pon15(i,k,j) ) +
*     +           ( pon14(i,k-1,j)*(dizt(k)/dizt(k-1)) - pon14(i,k,j) ) )
*           tot = on_rem(k) - (
*     +           ( pon15(i,k-1,j)*(dizt(k)/dizt(k-1)) - pon15(i,k,j) ) +
*     +           ( pon14(i,k-1,j)*(dizt(k)/dizt(k-1)) - pon14(i,k,j) ) )
*           if (ABS(dif).gt.epsilon(tot)) then
*           print*, " "
*           print*, " RELEASE OF NO3 NOT TRACKING THE COMBINED    "
*           print*, "        ISOTOPES OF ORG MATTER "
*           print*, " "
*           print*, " Tracer number =", n
*           print*," PON released (mmol) = ",on_rem(k)
*           print*," pon15(i,1,j) =",pon15(i,1,j)
*           print*," pon15(i,k-1,j) =",pon15(i,k-1,j)
*           print*," pon15(i,k,j) =",pon15(i,k,j)
*           print*," pon14(i,1,j) =",pon14(i,1,j)
*           print*," pon14(i,k-1,j) =",pon14(i,k-1,j)
*           print*," pon14(i,k,j) =",pon14(i,k,j)
*           print*," dizt(k-1) =", dizt(k-1),
*     +            " dizt(k) =", dizt(k)
*           print*, " "
*           print*, " total = ", tot
*           print*, " difference = ", dif
*           print*, " "
*           print*, 'i = ',i,' j = ',j,'k = ',k
*           print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
*           print*, 'kmt = ', kmt(i,j)
*           print*, " "
*           print*, " "
*           print*, "          STOPPING MODEL            "
*           print*, " "
*           print*, " "
*           stop
*           endif
*           
*           dif = on_rem_f(k)+((ponf15(i,k-1,j)*dizt(k)/dizt(k-1)
*     +           - ponf15(i,k,j)) +
*     +           (ponf14(i,k-1,j)*dizt(k)/dizt(k-1)
*     +           - ponf14(i,k,j)))
*           tot = on_rem_f(k)-((ponf15(i,k-1,j)*dizt(k)/dizt(k-1)
*     +           - ponf15(i,k,j)) +
*     +           (ponf14(i,k-1,j)*dizt(k)/dizt(k-1)
*     +           - ponf14(i,k,j)))
*           if (ABS(dif).gt.epsilon(tot)) then
*           print*, " "
*           print*, " RELEASE OF NO3 NOT TRACKING THE COMBINED    "
*           print*, "        ISOTOPES OF F-ORG MATTER "
*           print*, " "
*           print*, " Tracer number =", n
*           print*," PONf released (mmol) = ",on_rem_f(k)
*           print*," ponf15(i,1,j) =",ponf15(i,1,j)
*           print*," ponf15(i,k-1,j) =",ponf15(i,k-1,j)
*           print*," ponf15(i,k,j) =",ponf15(i,k,j)
*           print*," ponf14(i,1,j) =",ponf14(i,1,j)
*           print*," ponf14(i,k-1,j) =",ponf14(i,k-1,j)
*           print*," ponf14(i,k,j) =",ponf14(i,k,j)
*           print*," dizt(i,k-1,j) =", dizt(k-1),
*     +            " dizt(i,k,j) =", dizt(k)
*           print*, " "
*           print*, " total = ", tot
*           print*, " difference = ", dif
*           print*, " "
*           print*, 'i = ',i,' j = ',j,'k = ',k
*           print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
*           print*, 'kmt = ', kmt(i,j)
*           print*, " "
*           print*, " "
*           print*, "          STOPPING MODEL            "
*           print*, " "
*           print*, " "
*           stop
*           endif
           
           
*   8.4  --  Save delta 15N of organic matter arriving at the depth 

           if ((pon14(i,k-1,j)+ponf14(i,k-1,j)).gt.0.0) then
              d15Norg(i,k,j) = ( (pon15(i,k-1,j)+ponf15(i,k-1,j)) /
     +                         (pon14(i,k-1,j)+ponf14(i,k-1,j))-1.0)*1e3
           endif
           
        endif  ! --> k.gt.1  



*   8.5  --  Print output of interest to screen if N15diagnostics==TRUE

        if (no3.gt.0.0) then
        if (N15diagnostics .and. j.eq.latj .and. i.eq.loni) then
           print*, " "
           print*, " "
           print*, " RESULTS FOR SOURCE EQUATION "
           print*,"depth=",k,"lat=",j*1.59,"lon=",i*2.8125
           print*, "NO3 = ",no3,"15N = ",n15
           print*, "Water delta15N = ",(( (n15 / 
     +                                 (no3-n15))
     +                              / 1.0) -1.)*1000.
           print*, "NO3 assim = ",on_rem(k), on_rem_f(k)
           print*, "NO3 denit = ",no3_den(k), no3_den_f(k)
           print*, "15N assim = ",d_n15(i,k,j),
     +              d_n15f(i,k,j)
           print*, "15N denit = ",d_n15d(i,k,j),
     +              d_n15df(i,k,j)
           print*, "New delta15N = ", (( 
     +((n15-d_n15(i,k,j)-d_n15f(i,k,j)-d_n15d(i,k,j)-d_n15df(i,k,j))
     +      / ( (no3-on_rem(k)-on_rem_f(k)-no3_den(k)-no3_den_f(k)) -
     +(n15-d_n15(i,k,j)-d_n15f(i,k,j)-d_n15d(i,k,j)-d_n15df(i,k,j)) 
     +                  )) / 1.0) -1.)*1000.
        endif
        endif


        enddo ! do k = 1,kmt(i,j)   
        endif ! if kmt(i,j).gt.0
        enddo ! do i = istrt,iend
      endif !if n.eq.n_n15


* ****  END STEP 8  **** *


****************************************************************
* ****  STEP 9 -->  NITROUS OXIDE PRODUCTION/CONSUMPTION  **** *
****************************************************************

      if (n.eq.n_n2o) then 
      if (n_tou.gt.0 .and. n_oxy.gt.0 .and. 
     +    n_age.gt.0 .and. n_tomz.gt.0) then
         do i = istrt,iend

         if (kmt(i,j).ne.0) then  ! --> select for only ocean points
         
*   9.1  --  Shelf N2O Production (Freing et al. 2012 Phil Trans R Soc B)
         if (kmt(i,j).ge.1 .and. kmt(i,j).le.6) then
            do k = 3,km

            age = max(1e-3,tb(i,k,n_age))
            temp = max(-2.0, min(40.0, tb(i,k,1)))
            tour = (tb(i,k,n_tou)-tb(i,k,n_oxy))/age
            ! units umol/m3 yr-1
            d_n2o(i,k,j) = tour*0.1793*exp(-zt(k)/350.0) + 0.5374

            enddo
         endif
         
*   9.2  --  Open Ocean N2O Production (Freing et al. 2012 Phil Trans R Soc B)
         if (kmt(i,j).gt.6) then
            do k = 3,km

            age = max(1e-3,tb(i,k,n_age))
            
            if (N2Otemp) then

               temp = max(-2.0, min(40.0, tb(i,k,1)))
               tour =(tb(i,k,n_tou)-tb(i,k,n_oxy))/age
               ! units umol/m3 yr-1
               d_n2o(i,k,j) = tour*0.0665*exp(-temp/20000.0) - 0.0032
               

            elseif (N2Odepth) then

               tour =(tb(i,k,n_tou)-tb(i,k,n_oxy))/age
               ! units umol/m3 yr-1
               d_n2o(i,k,j) = tour*0.0658*exp(-zt(k)/20000.0) - 0.0065

            endif
         
            enddo
         endif
       
*   9.3 --  Production and consumption due to denitrification
         if (den) then
            do k = 3,kmt(i,j)

               if (op_den(i,k,j).ne.0.0 .or. op_den_f(i,k,j).ne.0.0)then

                  age = max(1e-3, tb(i,k,n_age))
                  tomz = max(1e-3, tb(i,k,n_tomz))

                  if (tb(i,k,n_oxy).gt.1.0) then
                  
                     n2o_Pden = max(0.0, 1.0/tomz * (
     +                          (omz_max - (d_n2o(i,k,j)+tb(i,k,n_n2o)))
     +                          - 0.5*(omz_max - omz_mean) ))
                     d_n2o(i,k,j) = d_n2o(i,k,j) + n2o_Pden
                  
                  else
                     
                     n2o_Cden = max(0.0, 1.0/tomz * (
     +                          ((d_n2o(i,k,j)+tb(i,k,n_n2o)) - omz_min)
     +                          - 0.5*(omz_mean - omz_min) ))
                     d_n2o(i,k,j) = d_n2o(i,k,j) - n2o_Cden
                  
                  endif

               endif
               
            enddo
         endif

         endif  ! --> kmt(i,j).gt.0

         ! No N2O production or consumption above 50 m
         d_n2o(i,1,j) = 0.0
         d_n2o(i,2,j) = 0.0

         enddo

      else

        print*, " "
        print*, " PLEASE TURN ON OXYGEN, TOU, AGE and OMZ-AGE TRACERS "
        print*, "   They are required for the calculation    "
        print*, "             of Nitrous Oxide.              "
        print*, " "
        stop
        
      endif  
      endif  ! --> n.eq.n_n2o
      

* ****  END STEP 9  **** *
      
*****************************************************
* ****  STEP 10 -->  COMPUTE SOURCES AND SINKS  **** *
*****************************************************

      do i = istrt,iend
         
         isurf = 1.
         fixer = 0.
      
      if (kmt(i,j).ne.0) then
         do k = 1,kmt(i,j)

*   10.1  --  Modify the P:C:O ratio since pop is used in source calculation
            vredf = 1.
            if (n.eq.n_dic .or. n.eq.n_c13 .or. n.eq.n_oxy) then 
               if (ivredc.eq.0) then 
                  vredf = 1.
               endif
               if (ivredc.eq.1) then
                  vredf = ( 1. + ( pco2o(i,j,n_dic-2) - 280. ) 
     +                      * 2. / 700. / 6.6 )
               endif
               if (ivredc.eq.2) then
                  vredf = (1. + (pco2a - 280.) * 2. / 700. / 6.6)
               endif
               if (ivredo.eq.0 .and. n.eq.n_oxy) then
                  vredf = 1.
               endif
            endif

*   10.2  --  Net source/sink for each tracer

               if (n.eq.n_oxy) then
                  source(i,k) = ( -op_rem(i,k,j)*o2_rem(i,j) 
     +                            -op_rem_f(i,k,j)*o2_rem_f
     +                            -sremG(i,k,j)*o2_rem(i,j)
     +                            -sremF(i,k,j)*o2_rem_f
     +                            -ic_sed(i,k,j)*ratio_pic(n)
     +                            -pic(i,k,j) * ratio_pic(n) +
     +                            fluxgas(i,j,n-2) * isurf )*dizt(k)
               endif
               if (n.eq.n_dic) then
                  source(i,k) = (( -op_tot(i,k,j)*carb2P(i,j) 
     +                             -op_tot_f(i,k,j)*CP_fix
     +                             -op_sed(i,k,j)*carb2P(i,j)
     +                             -op_sed_f(i,k,j)*CP_fix )*vredf  
     +                             -ic_sed(i,k,j)*ratio_pic(n)
     +                             -pic(i,k,j) * ratio_pic(n) +
     +                              fluxgas(i,j,n-2) * isurf ) * dizt(k)
               endif
               if (n.eq.n_c13) then
                  if (k.eq.1) then
                     ratio = ( tb(i,1,n) + tr_off(n) ) /
     +                        ( tb(i,1,n_dic)+tr_off(n_dic) )
                  endif
                  source(i,k) = (((-op_tot(i,k,j)*c13toc(i,j)
     +                             -op_tot_f(i,k,j)*CP_fix*0.988
     +                             -op_sed(i,k,j)*c13toc(i,j)
     +                             -op_sed_f(i,k,j)*CP_fix*0.988
     +                             )*vredf 
     +                             -ic_sed(i,k,j)*ratio_pic(n)*0.998
     +                             -pic(i,k,j)*ratio_pic(n)*0.998)*ratio
     +                           + fluxgas(i,j,n-2)*isurf ) * dizt(k)
               endif
*               if (n.eq.n_c14) then  ! mol/m3
*                 source(i,k) = -(log(2.0)/(5700*365*86400.))*tb(i,k,n)
*     +                          + fluxgas(i,j,n-2)*isurf*dizt(k)
*               endif
               if (n.eq.n_pho) then
                  source(i,k) = ( -op_tot(i,k,j) - op_tot_f(i,k,j)
     +                            -op_sed(i,k,j) - op_sed_f(i,k,j)
     +                            -pic(i,k,j)*ratio_pic(n) 
     +                            -ic_sed(i,k,j)*ratio_pic(n)
     +                            +fluxgas(i,j,n-2)*isurf ) * dizt(k)
               endif
               if (n.eq.n_alk) then
*               !See Wolf-Gladrow 2007 Marine Chem for explanation
                  source(i,k) = (( -op_tot(i,k,j)*NtoP(i,j) 
     +                             -op_sed(i,k,j)*NtoP(i,j)
     +                             -op_tot_f(i,k,j)*NP_fix*fixer
     +                             -op_sed_f(i,k,j)*NP_fix*fixer
     +                             -op_den(i,k,j)*no3_rem(i,j)
     +                             -op_den_f(i,k,j)*no3_rem_f
     +                             -sdenG(i,k,j)*no3_rem(i,j)
     +                             -sdenF(i,k,j)*no3_rem_f)*(-1)
     +                             -pic(i,k,j)*ratio_pic(n)
     +                             -ic_sed(i,k,j)*ratio_pic(n)
     +                             +fluxgas(i,j,n-2)*isurf)*dizt(k)
               endif
               if (n.eq.n_fe) then   
                  source(i,k) = ( (-op_tot(i,k,j)*ratio_pop(n) 
     +                             -op_tot_f(i,k,j)*0.64
     +                             -op_sed(i,k,j)*ratio_pop(n)
     +                             -op_sed_f(i,k,j)*0.64)*vredf !0.64 is Fe:P of Kustka et al 2003 LO
     +                             -pic(i,k,j)*ratio_pic(n)
     +                             -ic_sed(i,k,j)*ratio_pic(n)
     +                            + fluxgas(i,j,n-2)*isurf)*dizt(k)
     +                          - max(0.0, (tb(i,k,n)) - 0.60) 
     +                          / (86400. * 365)
               endif
               if (n.eq.n_no3) then
                  ! Capture pelagic denitrification rate 
                  denP(i,k,j) = ( -op_den(i,k,j)*no3_rem(i,j) 
     +                            -op_den_f(i,k,j)*no3_rem_f )
     +                            *dizt(k)

                  ! Capture sedimentary denitrification rate 
                  denS(i,k,j) = ( -sdenG(i,k,j)*no3_rem(i,j)
     +                            -sdenF(i,k,j)*no3_rem_f )
     +                            *dizt(k)
                  ! Capture nitrogen fixation rate
                  Nfix(i,k,j) = (-op_tot_f(i,k,j)*NP_fix
     +                           -op_sed_f(i,k,j)*NP_fix)*fixer*dizt(k)
                  Nexp(i,k,j) = (-op_tot(i,k,j)*NtoP(i,j)
     +                           -op_sed(i,k,j)*NtoP(i,j))*dizt(k)

                  source(i,k) = ( -op_tot(i,k,j)*NtoP(i,j) ! water column
     +                            -op_sed(i,k,j)*NtoP(i,j)  ! sediments
     +                            -pic(i,k,j)*ratio_pic(n)
     +                            -ic_sed(i,k,j)*ratio_pic(n)
     +                            + (fluxgas(i,j,n-2)*isurf) )*dizt(k)
     +                            + denP(i,k,j)+Nfix(i,k,j)+denS(i,k,j)

               endif
               if (n.eq.n_n15) then
                   source(i,k) = (( -d_n15(i,k,j) - d_n15f(i,k,j)
     +                              -d_n15d(i,k,j) - d_n15df(i,k,j) )
     +                           / (c2dtts)) +
     +                           (fluxgas(i,j,n-2)*isurf)*dizt(k)
               endif
               if (n.eq.n_tou) then
                  source(i,k) = (fluxgas(i,j,n-2)*isurf) * dizt(k)
               endif
               if (n.eq.n_n2o) then
                  source(i,k) = (d_n2o(i,k,j)*(1.0/86400.0) +
     +                           fluxgas(i,j,n-2)*isurf) * dizt(k)
               endif
               
            isurf = c0
            fixer = 1.

         enddo
      endif
      enddo


#ifdef rjm_geo_seq
      if (igas(n-2).eq.2) then 
         if (n.eq.n_dic .and. j.eq.2) then
            ttflux(1) = ttflux(1) + 1
         endif
         do i = 2,imt-1
            da = dxt(i)*dyt(j)*cst(j)
            ttflux(n) = ttflux(n) + fluxgas(i,j,n-2) * da
         enddo

         ! write fluxes
         if (j.eq.jmt-1) then
            fc = 1e-5 * 1e-4 * 12e-15  ! convert flux to GtC/s 
            fcc = gCtppm * tott * spy
            stflux(n) = ttflux(n) * fc / max(ttflux(1), 1.)
            if (ifix(tott).eq.tott ) then
               write(91,*) tott, n, dpco2a, geo_seq(n)*fcc,
     +                     stflux(n)*fcc, stflux(n_dic)*fcc,
     +                     ( geo_seq(n)-stflux(n_dic)+stflux(n) )*fcc,
     +                     ttflux(1)
            endif
         endif
      endif
#endif

* ****  END STEP 10  **** *
      
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c




c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine tx_growth (j)

c	Author: r.j.matear
c 	
c    Compute adjustment to export production based on variable 
c    temperature and mixed layer depth
c
c......................................................................c

      include "obgc.h"
      include "bio.h"
#ifdef dmr_tconstant
      include "oline.h"
#endif

* Time terms
      parameter (istrt=2, iend=imt-1)
      c0=0.
      p5=.5
      par=pmodel(9)
      alpha=pmodel(10)
      dps = 1/86400.   ! inverse of seconds per day

* Compute mixed layer depth
      do i = istrt,iend 

         k = 1
         mld = k
         
         do k = 1,km 
            if (mld.eq.k .and. k.le.kmt(i,j)-1) then
               if (abs (tb(i,1,1) - tb(i,k,1)) .lt. .01) then
                  mld = k
               endif
            endif
         enddo
        
         xmld(i,j) = (zt(mld) + dzt(mld)*p5)   !cm
         xmld(i,j) = dzt(1)

* factor for the mld variability
         fmld = dzt(1)/xmld(i,j)  ! reduced growth with mld

#ifdef dmr_tconstant
         tsurf = 10
#endif

         fpgr_temp(i,j) = c0
         if (kmt(i,j) .ne. 0 ) then 
            ! Compute Temperature and mixed layer depth adjustments
            tsurf = max(-2.0, min(40.0, tb(i,1,1)))
            fmax =  0.6*(1.066)**tsurf
            f1 = sbcbio(i,j,2)*alpha*par

            ! Evans and Parslow PI equation
            !	 	fi = fmax*f1/sqrt(fmax**2 + f1**2 )
            ! Brian Griffiths PI equation
            fi = fmax * (1 - exp( -f1/fmax) )

            Vmax(i,j) = fmax*dps  ! per second
            ! old   fpgr_temp(i,j) = fi*fmld*s_npp*dps  !per second
            fpgr_temp(i,j) = fi*fmld*dps / Vmax(i,j) !unitless


         endif

      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c



c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine stoich (i,j,phs,no3)

c
c   Routine to determine the stoichiometry of Carbon, Hydrogen, Oxygen
c   and Nitrogen to Phosphate based on the amount of Phosphate and
c   Nitrogen in any given surface grid box.
c       INPUT:
c           i   =   longitudinal index for grid
c           j   =   latitudinal index for grid
c           phs =   phosphate concentration at grid(i,j)
c           no3 =   nitrate concentration at grid (i,j)
c           EP  =   export production at grid (i,j) mg C m-2 hr-1
c       OUTPUT:
c           carb2P(i,j) =   Carbon:Phosphorus ratio 
c           NtoP(i,j)   =   Nitrogen:Phosphorus ratio
c           HtoP(i,j)   =   Hydrogen:Phosphorus ratio
c           OtoP(i,j)   =   Oxygen:Phosphorus ratio
c           
c   Remin_o2(i,j)   =   moles of O2 removed during aerobic reminen
c   Remin_no3(i,j)  =   moles of NO3 removed during anaerobic remin
c
c        AUTHOR:
c           Pearse James Buchanan          
c           pearse.buchanan@utas.edu.au
c
c......................................................................c

        include "obgc.h"
        include "bio.h"
        
        INTEGER :: mld
        REAL, DIMENSION(i,j) :: NtoC  
        REAL :: phospho, nucleic, ATP, protein, carbs,
     +          lipids, PUFA, SFA, MUFA, lipidC, lipidH, lipidO,
     +          nucleic_g, protein_g, phospho_g, carb_g, ATP_g, 
     +          lipid_g,
     +          nucleicC, nucleicH, nucleicO, nucleicN, nucleicP,
     +          proteinC, proteinH, proteinO, proteinN, 
     +          phosphoC, phosphoH, phosphoO, phosphoN, phosphoP,
     +          ATPC, ATPH, ATPO, ATPN, ATPP,
     +          carbC, carbH, carbO,
     +          molarC, molarH, molarO, molarN, molarP,
     +          decay, solar
        

* 1. determine the C:P and N:P ratios based on the nutrient limitation
*    scheme defined by Galbraith & Martiny (2015) PNAS, which first
*    caluclates C:P and N:C, after which we can calculate N:P.        
        if (kmt(i,j).gt.0) then
        carb2P(i,j) = 1.0/((6.9*phs + 6.0)/1000.0) * fm(i,1)
        
        if (n_no3.gt.0) then
           NtoC(i,j) = 0.125 + 0.03*no3/(0.32+no3) * fm(i,1)
           NtoP(i,j) = carb2P(i,j)*NtoC(i,j) * fm (i,1)
        else
           NtoP(i,j) = 16.0
        endif



        if (biomol) then
* 2. calculate % by weight of phosphoglycerides
           phospho = 5. + 10. * (phs/(phs+1.0))

* 3. prescribe nucleic acid and ATP %           
           nucleic = 9.
           ATP = 0.1

* 4. calculate % of proteins
           molarP = ( nucleic*(31./341.) + phospho*(31./714.72) +
     +                ATP*(93./507.) )/31.
           molarN = NtoP(i,j)*molarP
           protein = ( molarN*14. - nucleic*(52.5/341.) - 
     +                 phospho*(6.02/714.72) - ATP*(70./507.) )
     +                / (560./3320.)

* 5. calculate carbs and lipid %
           if (solarcarb) then
              ! find mixed layer depth
              k = 1
              mld = k
       
              do k = 1,km 
                 if (mld.eq.k .and. k.le.kmt(i,j)-1) then
                    if (abs (tb(i,1,1) - tb(i,k,1)) .lt. .01) then
                       mld = k
                    endif
                 endif
              enddo
              xmld(i,j) = (zt(mld) + dzt(mld)*0.5)/100.0   !cm --> m
       
              ! find solar saturation
              decay = 0.01148   !(Morel & Maritorena, 2001)
              solar = (obc(i,j,m,2)/10.) * exp(-decay*xmld(i,j))
        
              ! find carbohydrate %
              carbs = (100. - protein - phospho - nucleic - ATP)*0.5
     +                *(0.5 + 1.4*(solar/10.0))
           else
              carbs = (100. - protein - phospho - nucleic - ATP)*0.5*1.4
           endif
           
           lipids = 100. - protein - phospho - nucleic - ATP - carbs
        
* 6. calculate lipid elemental composition
           PUFA = (41.2 - 0.433*tb(i,1,1))/100.0
           SFA = (27.3 + 0.279*tb(i,1,1))/100.0
           MUFA = 1. - PUFA - SFA

           lipidC = PUFA*21. + MUFA*18. + SFA*15.
           lipidH = PUFA*31. + MUFA*30. + SFA*34.
           lipidO = 2.

* 7. calculate the grams of each element in the organic matter
           nucleic_g = 12*9.625 + 14.0 + 16*8.0 + 14*3.75 + 31.0
           protein_g = 12*147.0 + 228.0 + 16*46.0 + 14*40.0 + 32.0
           phospho_g = 12*37.9 + 72.5 + 16*9.4 + 14*0.43 + 31.0
           carb_g = 12*6.0 + 10.0 + 16*2.0 
           ATP_g = 12*10.0 + 16.0 + 16*13.0 + 14*5.0 + 3*31.0
           lipid_g = 12.0*lipidC + lipidH + lipidO*16.0
           
           nucleicC = (12.0*9.625)/nucleic_g
           nucleicH = 14.0/nucleic_g
           nucleicO = (16.0*8.0)/nucleic_g
           nucleicN = (14.0*3.75)/nucleic_g 
           nucleicP = 31.0/nucleic_g

           proteinC = (12.0*147.0)/protein_g
           proteinH = 228.0/protein_g
           proteinO = (16.0*46.0)/protein_g
           proteinN = (14.0*40.0)/protein_g 

           phosphoC = (12.0*37.9)/phospho_g
           phosphoH = 72.5/phospho_g
           phosphoO = (16.0*9.4)/phospho_g
           phosphoN = (14.0*0.43)/phospho_g 
           phosphoP = 31.0/phospho_g

           carbC = (12.0*6.0)/carb_g
           carbH = 10.0/carb_g
           carbO = (16.0*5.0)/carb_g

           ATPC = (12.0*10.0)/ATP_g
           ATPH = 16.0/ATP_g
           ATPO = (16.0*13.0)/ATP_g
           ATPN = (14.0*5.0)/ATP_g 
           ATPP = (31.0*3.0)/ATP_g

           lipidC = (lipidC*12.0)/lipid_g
           lipidH = lipidH/lipid_g
           lipidO = (lipidO*16.0)/lipid_g

* 8. calculate the molar mass of each element
           molarC = (nucleicC*nucleic + proteinC*protein + 
     +               phosphoC*phospho + carbC*carbs + ATPC*ATP + 
     +               lipidC*lipids)/12.
           molarH = (nucleicH*nucleic + proteinH*protein + 
     +               phosphoH*phospho + carbH*carbs + ATPH*ATP + 
     +               lipidH*lipids)
           molarO = (nucleicO*nucleic + proteinO*protein + 
     +               phosphoO*phospho + carbO*carbs + ATPO*ATP + 
     +               lipidO*lipids)/16.
           molarN = (nucleicN*nucleic + proteinN*protein + 
     +               phosphoN*phospho + ATPN*ATP)/14.
           molarP = (nucleicP*nucleic +  phosphoP*phospho +ATPP*ATP)/31.
           
* 9. calculate the new C:H:O:N:P ratios
           carb2P(i,j) = molarC/molarP
           HtoP(i,j) = molarH/molarP
           OtoP(i,j) = molarO/molarP
           NtoP(i,j) = molarN/molarP
           
        else

           HtoP(i,j) = 2*carb2P(i,j) + 3*NtoP(i,j) + 3
           OtoP(i,j) = carb2P(i,j) + 4

        endif
        endif


        
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c


c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine atm_source (joff, js, je, is, ie, n)


#ifdef rjm_atm_source
c
c
c  Routine to determine the surface sources for the co2 perturbation
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c    author:  r.j.matear  e-mail richard.matear@ml.csiro.au
c
c  Setup for the following tracer order
c  dCO2
c......................................................................c

      include "param.h"
      include "mw.h"

      include "bio.h"
      include "bio_sbc.h"
      
      ! grid info
      include "coord.h"
      include "grdvar.h"
      
      ! c2dtts
      include "scalar.h"

      ! kmt
      include "levind.h"

      ! Time terms
      include "tmngr.h"

      parameter (istrt=2, iend=imt-1)
      save pbtime,itimes
      data pbtime,itimes /0.,0/

*   Index to use for biology
      ibdt=taum1
#ifdef ios_upstrm
      ibdt=tau
#endif
    
      if (n.eq.3 .and. joff.eq.0) then
         itimes=itimes+1
      endif

*   Modify pCO2 and C13 in atmosphere
#ifdef rjm_anthropogenic
      xspin=0.
      if (pbtime.eq.0) then
         pbtime=prelyr
      endif
      
      ttt=prelyr-pbtime
      
      if (ttt .le. xspin ) then
         pco2a=280.
      else if (ttt .le.100+xspin) then
         pco2a = 280 + 0.28048302*(ttt-xspin)
      else
         pco2a= 28.48302*exp( (ttt-100.+xspin)/41.936498)+280
      endif
      
      ! old value	fc13a = 0.01 - 1.818571e-7*(pco2a-280.)
      fc13a = 0.011164382 - 2.0225714e-7*(pco2a-280.)
      pco2a=pco2a-280 
      print *,'Pco2 evolution',ttt,prelyr,pbtime,pco2a,fc13a

#else

      pco2a=280 -280  ! set to current value
      pco2a=pco2_atm(ibdt)
      fc13a =0.011164382  ! old value .01
      print *,'Pco2 evolution',ttt,prelyr,pco2a
#endif

*   Calculate air-sea flux for other tracers
      if (igas(n-2).gt.0) then 
         do j = js,je
            jrow = j +joff
            do i = istrt,iend

               ! Only calculate for ocean points
               if (tmask(i,1,j) .ne. 0) then
                  temp=t(i,1,j,1,ibdt)
                  saln=max (0., (t(i,1,j,2,ibdt)+.035)*1e3 )

               ! For CO2 exchange
               !   umol/l *cm/s
               !   or 1e-2 mmol/m2/s
                  if (igas(n-2) .eq. 2) then

                     tdic = t(i,1,j,n,ibdt)+tr_off(n)
                     tdic_2 = tdic*tdic
                     tdic_3 = tdic_2*tdic
                     tdic_4 = tdic_2*tdic_2
                     tdic_5 = tdic*tdic_4
                 
                     pco2o(i,jrow,n_dic-2) = (1.5568 - 1.3993e-2*temp) 
     +                                       *tdic 
     +                               +(7.4706-0.20207*temp)*1e-3*tdic_2 
     +                               -(1.2748-0.12015*temp)*1e-5*tdic_3 
     +                               +(2.4491-0.12639*temp)*1e-7*tdic_4 
     +                               -(1.5468-0.15326*temp)*1e-10*tdic_5

                     !Compute CO2 solubility once and store in pop
                     if (itimes .eq. 1) then
                        tk = temp+273.15
                        tk2 = tk/100.0
                        pop(i,1,jrow) = exp( -60.2409 + 9345.17/tk + 
     +                                       23.3585*Alog(tk2) +
     +                                       saln*(2.3517E-2 - 
     +                                       2.3656E-2*tk2 + 4.7036E-3 *
     +                                       tk2*tk2) )
                     endif

                     x_co2(i,jrow) = 7.964 * sbcbio(i,jrow,1) /
     +                               schmidt_no(2,temp)**.5 * 
     +                               pop(i,1,jrow) * 2.7777778e-4
                     fluxgas(i,jrow,n-2) = x_co2(i,jrow) *
     +                                     ( pco2a - 
     +                                     pco2o(i,jrow,n_dic-2) )

#ifdef rjm_co2_dumping
                     if (i.eq.idump .and. jrow.eq.jdump) then
                        print*,'pco2',i,jrow,tdic,pco2o(i,jrow,n_dic-2),
     +                      fluxgas(i,jrow,n-2)
                     endif
#endif

                     ! C13 exchange (must include CO2 exchange)
                  elseif (igas(n-2) .eq. 5) then
                     !   cas=0.98940
                     !   csa=0.9984

                     ! Temperature Dependent fractionation
                     temk = 1./(temp+273.15)
                     cas = (-0.373*temk + 1.00019) *.9995
                     csa = (-9.866*temk + 1.02412) *.9995

                     fluxgas(i,jrow,n-2) = cas * x_co2(i,jrow) * pco2a
     +                                     * fc13a - csa * x_co2(i,jrow)
     +                                     * pco2o(i,jrow,n_dic-2)
     +                                     * ( t(i,1,j,n,ibdt) 
     +                                     + tr_off(n) )
     +                                     / ( t(i,1,j,n_dic,ibdt) 
     +                                     + tr_off(n_dic) )

                  endif    ! end tracer type 
               endif    ! end mask condition
            enddo    ! end j loop
         enddo    ! end i loop
      endif    ! end atm-source condition


* Compute sources and sinks
      do j = js,je
         jrow = j+joff

         do i = istrt,iend

#ifdef rjm_co2_dumping
            if (jrow.eq.jdump) then
               source(idump,kdump,j) = xdump * cstdxtr(i,j) * 
     +                                 dytr(jrow) * dztr(kdump)
            endif
#endif

            isurf=1.
            if (igas(n-2) .ne. 5) then
               source(i,1,j) = (fluxgas(i,jrow,n-2)*isurf )*dizt(1)
            else 
               ratio = ( t(i,1,j,n,ibdt) + tr_off(n) ) /
     +                 ( t(i,1,j,n_dic,ibdt) + tr_off(n_dic) )
               source(i,1,j) = ( fluxgas(i,jrow,n-2) * isurf ) * dizt(1)
            endif
         enddo
      enddo


#endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c





c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine sediment(i,k,j,n)
      
c                                                                      c
c                                                                      c
c                                                                      c
c......................................................................c
      
      include "obgc.h"
      include "bio.h"
      include "extra.h"
      
      dimension twodt(km)

      twodt(k) = c2dtts

c Check that oxygen is not less than zero
      if (n.eq.3 .and. ta(i,k,n).lt.0) then 

         rr = ratio_pop(n)
         
         if (ratio_pop(n).eq.0) then
            rr = 1.
         endif 
     
         rr=1./rr

         if (k.gt.1 .and. k.lt.kmt(i,j) ) then
            sred = min( ta(i,k,n), 0. ) / twodt(k) !*dzt(k)/twodt(k)
            ta(i,k,n) = max( ta(i,k,n), 0. )
            source(i,k+1) = source(i,k+1) + sred * dzt(k) * dizt(k+1)
            sred = sred * rr * dzt(k)
            pop(i,k,j) = pop(i,k,j) + sred
            pop(i,k+1,j) = pop(i,k+1,j) - sred
         elseif (k.eq.kmt(i,j) ) then
            sediments(i,j) = min( ta(i,k,n), 0. ) * dzt(k) 
     +                            * rr / twodt(k)
            ta(i,k,n) = max( ta(i,k,n), 0. )
            pop(i,k,j) = pop(i,k,j) + sediments(i,j)
            sum = sum + sediments(i,j)
         endif
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c




c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      function rc_atm(pco2)

c                                                                      c
c                                                                      c
c computes the c13/c12 ratio in the atmosphere                         c
c......................................................................c

      include "obgc.h"
      include "bio.h"
      
      ri=0.011164382
      ci=280./.48
      rc_atm=(ci*ri - atm_box_tot_c13)/(ci-atm_box_tot)
      !	print*,'Atm Box',ttt
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c




c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      function co2_history(ttt)

c
c
c......................................................................c

      common /co2h/ ccc(501),xxx(501)

#ifdef dmr_atm_box
      include "obgc.h"
      include "bio.h"
#endif

      data icall /0/
      save icall,ico2 

* Time to 1881 used in the co2 history
      iyr = 0  ! no need to delay start to 1880
      tstart = 1880-iyr

#ifdef dmr_atm_box
      if (icall.eq.0) then
         open(12,file='atm_box.dat',status='old')
         read(12,*) xt,pc,tot,ann,rc,tot_c13
         close(12)
         atm_box_tot=tot
         atm_box_ann=0
         atm_box_ann_c13=0
         atm_box_tot_c13=tot_c13
         atm_box_rc=rc
         icall=1
      endif
      co2_history=280 - atm_box_tot*.48
      !	print*,'Atm Box',ttt,co2_history,atm_box_tot,atm_box_ann
#else
      if (icall.eq.0) then
         atm_box_tot=0
         open(11,file='co2_history.dat',status='old')
         do i = 1,501
            read(11,*,end=91) xxx(i),ccc(i)
         enddo

 91      print*,'rjm number of co2 points is ', i-1
         ico2 = i-1
         icall=1
         close(11)
      endif

      tstart = xxx(1)-iyr   ! use the first year in the file
      if (ttt.lt.iyr) then
         co2_history=280
      else
         i1 = ttt+1-iyr
         if (i1.le.ico2-1) then 
            i2 = i1+1
            co2_history = ( ccc(i2) - ccc(i1) ) * 
     +                    (ttt - xxx(i1) + tstart) + ccc(i1)
         else
            co2_history = ccc(ico2)
         endif
      endif
#endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c



c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine file_name2(i1,a7)

c......................................................................c

      character *10 a7
      character *1 a1
      character *2 a2
      character *3 a3
      a7 = 'svmo.nc000'

#ifdef cexp4
      i2 = i1-1
      i2a = i2/50
      i3 = i2a*50 +1
#else
      i2 = i1-56
      i2a = i2/15
      i3 = i2a*15 +56
#endif

      if (i3.le.9) then
         write(a1,'(i1)') i3
         a7(10:10)=a1
      elseif (i3.le.99) then
         write(a2,'(i2)') i3
         a7(9:10)=a2
      else
         write(a3,'(i3)') i3
         a7(8:10)=a3
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c======================================================================c




c======================================================================c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......................................................................c

      subroutine file_name(i1,a7)

c......................................................................c

      character *18 a7
      character *1 a1
      character *2 a2
      character *3 a3
      a7 = 'archive/com0000.nc'
      if (i1.le.9) then
         write(a1,'(i1)') i1
         a7(15:15) = a1
      elseif (i1.le.99) then
         write(a2,'(i2)') i1
         a7(14:15) = a2
      else
         write(a3,'(i3)') i1
         a7(13:15) = a3
      endif
      return
      end

#endif
