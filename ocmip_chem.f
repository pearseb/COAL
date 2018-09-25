* routine to read in the forcing information 
* gas exchange coefficient, fice, and surface pressure
      subroutine nread_ocmip
* Read forcing data from a netcdf file
      include "netcdf.inc"

      include "obgc.h"
      include "bio.h"
      real*4 tvar1(imt,jmt,12),tvar2(imt,jmt,12),
     +       tvar3(imtm2,jmtm2,12),tvar4(imtm2,jmtm2,km),
     +       tvar5(imt,jmt,km)
      integer s1(3), c1(3), c2(3)
      data s1 /1,1,1/, c1 /imtm2,jmtm2,12/, c2 /imtm2,jmtm2,km/

* zero arrays
      do m=1,12
         do i=1,imt
            do j=1,jmt
               tvar1(i,j,m) = 0
               tvar2(i,j,m) = 0
            enddo
         enddo
      enddo
      do k=1,km
         do i=1,imt
            do j=1,jmt
               tvar5(i,j,k) = 0
            enddo
         enddo
      enddo

* Get file 
      
      id_f=ncopn("SGBfrac.nc",ncnowrit,icode)
      write(6,*) "Reading SGBfrac.nc file"
      id_v=ncvid(id_f, "SGBfrac", icode)
      call ncvgt(id_f, id_v, s1, c2,tvar4,icode)

      tvar5(2:imtm1,2:jmtm1,1:km) = tvar4(:,:,:)

      write(6,*) "Saving sub-grid scale bathymetry for sediment"
      do k=1,km
         do j=1,jmt
            do i=1,imt
               sgb_frac(i,j,k) = tvar5(i,j,k) 
            enddo
         enddo
      enddo

      print*, " "
      print*, " "


      id_f=ncopn("gasx_ocmip2.nc",ncnowrit,icode)
      write(6,*) "Reading gasx_ocmip2.nc file"
      
      id_v=ncvid(id_f, "WND2", icode)
      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3(:,:,:)
      
      id_v=ncvid(id_f, "FICE", icode)
      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
      tvar2(2:imtm1,2:jmtm1,1:12) = tvar3(:,:,:)
      
      write(6,*) "Preparing surface wind sqaured"
      do m=1,12
         do j=1,jmt
            do i=1,imt
               obc(i,j,m,1) = max(0.0,tvar1(i,j,m)*(1-tvar2(i,j,m)))
            enddo
         enddo
      enddo
      
      
      id_v=ncvid(id_f, "NSWRS", icode)
      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3(:,:,:)
      
      write(6,*) "Preparing Net Short Wave Radiation"
      do m=1,12
         do j=1,jmt
            do i=1,imt
               obc(i,j,m,2) = max(0.0, tvar1(i,j,m)*(1-tvar2(i,j,m))) !scale by seaice [] 
            enddo
         enddo
      enddo
     
       ! PRESSURE NOT USED ! 
*      id_v=ncvid(id_f, "P", icode)
*      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
*      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3
      
      write(6,*) "Storing Sea ice field"
      do m=1,12
         do j=1,jmt
            do i=1,imt
               obc(i,j,m,3) = max(0.0, min(1.0,tvar2(i,j,m)))
            enddo
         enddo
      enddo
      
      
      id_v=ncvid(id_f, "DUST", icode)
      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3(:,:,:)
      
      write(6,*) "Storing Aeolian Dust field"
      do m=1,12
         do j=1,jmt
            do i=1,imt
               ! input umol/m2/s= pmol/kg/s m but need to convert to pmol/kg/s cm
               obc(i,j,m,4) = max(0.0,tvar1(i,j,m)*100*(1-tvar2(i,j,m)))! convert to cm 
            enddo
         enddo
      enddo
      
      
      id_v=ncvid(id_f, "NOX", icode)
      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3(:,:,:)
      
      write(6,*) "Storing NOx deposition field"
      do m=1,12
         do j=1,jmt
            do i=1,imt
               ! input mmol/m2/s= umol/kg/s m but need to convert to umol/kg/s cm
               obc(i,j,m,5) = max(0.0,tvar1(i,j,m)*100*(1-tvar2(i,j,m)))! convert to cm 
            enddo
         enddo
      enddo
      
      
      id_v=ncvid(id_f, "P", icode)
      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3(:,:,:)
      
      write(6,*) "Storing Atmospheric Pressure field"
      do m=1,12
         do j=1,jmt
            do i=1,imt
               obc(i,j,m,6) = max(0.0,tvar1(i,j,m)) ! units in atm 
            enddo
         enddo
      enddo
      
*      id_v=ncvid(id_f, "DayLength", icode)
*      call ncvgt(id_f, id_v, s1, c1,tvar3,icode) 
*      tvar1(2:imtm1,2:jmtm1,1:12) = tvar3
*      
*      write(6,*) "Storing DayLength field"
*      do m=1,12
*         do j=1,jmt
*            do i=1,imt
*               obc(i,j,m,5) = tvar1(i,j,m) 
*            enddo
*         enddo
*      enddo
       
      return
      end

c_
c_ ---------------------------------------------------------------------
c_ 
	subroutine co2flux (t,s,kw660,ppo,dic,alk,po4,si,dic14,
	1	xco2,xc14, kwco2,  co2ex,c14ex,pco2surf,dpco2)
c
c**********************************************************************
c
c  Computes the time rate of change of DIC in the surface
c  layer due to air-sea gas exchange in mol/m^3/s.
c
c  Inputs:
c    t        model surface temperature (deg C)
c    s        model surface salinity (permil)
c    kw660    gas transfer velocity at a Schmidt number of 660, accounting
c               for sea ice fraction (m/s)
c    ppo      surface pressure divided by 1 atm
c    dic      surface DIC concentration (mol/m^3)
c    alk      surface alkalinity (eq/m^3)
c    po4      surface phosphate concentration (mol/m^3)
c    si       surface silicate concentration (mol/m^3)
c    xco2     atmospheric CO2 mixing ratio (ppm)
c  Output:
c    kwco2    gas exchange coefficient
c    co2ex    time rate of change of DIC in the surface layer due
c               to air-sea exchange (mol/m^3/s)
c
c  Two functions are called:
c    scco2    Schmidt number of CO2
cc   co2sato  CO2 saturation concentration at 1 atm (mol/m^3)
c  One subroutine is called:
c    co2calc  Computes actual aqueous CO2 concentration in surface
c             layer (mol/m^3), the saturation concentration at a
c             given atmospheric pressure, and there difference. 
c
c  Numbers in brackets refer to equation numbers in simulation design
c  document.
c
c  Ray Najjar 1/29/99
c  Modifications: James Orr 4/4/99 
c                 (in particular, see changes near "cc" prefix)
c
c**********************************************************************
c
      real kwco2,kw660
c
c  Compute the transfer velocity for CO2 in m/s [4]
c
cc    Somethings seems to be missing here, and units change is unneccessary
cc    because OCMIP gas exchange fields (i.e., xKw) is already in m/s!
cc    kwco2 = (scco2(t)/660)**-0.5*0.01/3600.0
      kwco2 =          Kw660   * (660/ScCO2(t))**0.5 
c
c  ---------------------------------------------------------------------------
c  >>> RGN: Do the carbonate equilibrium calculations.  Pick values of
c  >>>      ph range.  You may want to change these per suggestions of Sabine
c  >>>      and Key.  Note also that co2calc subroutine expects 
c  >>>      concentrations in mol/kg and eq/kg, so average surface density 
c  >>>      is used to convert.
C
c  >>> JCO: Yes, but one can pass arguments in mol/m^3 directly,
c  >>>      IF input and output arguments are consistent (i.e., always in
c  >>>      per volume basis, not per mass). See my comments in modified 
c  >>>      co2calc.f.  Related errors are small.

      phlo = 3.0
      phhi = 12.0

cc    dicin = dic/rhoavg
cc    alkin = alk/rhoavg
cc    po4in = po4/rhoavg
cc    siin = si/rhoavg
      dicin = dic
      alkin = alk
      po4in = po4
      dic14in = dic14
      siin = si
c  ---------------------------------------------------------------------------

cc    call co2calc(t,s,dicin,alkin,po4in,siin,phlo,phi,ph,co2star)
      CALL co2calc(t,s,dicin,alkin,po4in,siin,phlo,phhi,ph
	1	,xco2,ppo
     &            ,co2star,dco2star,pCO2surf,dpco2)

c  Compute the saturation concentration for CO2 [3]
c
cc    co2sat = co2sato(t,s,xco2)*ppo
c
c  Compute time rate of change of CO2 due to gas exchange [1]
c
c dropped the 1/dz1
cc    co2ex = kwco2*(co2sat-co2star)/dz1
      co2ex = kwco2*(dco2star)

c  Compute time rate of change of c14 due to gas exchange
	ratm = (1+xc14*1e-3)
	rocn = dic14in/(dicin+1e-10)
      c14ex = kwco2*( (dco2star+co2star)*ratm  - co2star*rocn)

cc    Further handling of output arguments (1) pCO2surf (pCO2ocn) and 
cc    (2) dpCO2 (pCO2ocn - pCO2atm) are left to the discretion of the 
cc    the modeler.  Both arguments are needed for standard OCMIP-2 output. 
c
      return
      END







c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /cs/home/csdmr/csrjm/mark2/code/RCS/ocmip_chem.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 2000/01/04 03:43:01 $   ;  $State: Exp $
c_ $Author: csrjm $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: ocmip_chem.f,v $
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c_ Revision 1.8  1999/07/16 11:40:33  orr
c_ Modifications by Keith Lindsay to fix inconsistency with common block
c_ "species" (not the same in "ta_iter_1.f").
c_ Also comment lines changed/added by J. Orr.
c_
c_ Revision 1.7  1999/04/26  13:04:54  orr
c_ Modified USAGE comment, to include new arguments
c_
c_ Revision 1.6  1999/04/14 12:55:52  orr
c_ Changed units for input arguments for tracers:
c_ formerly in mol/metric ton (T); now in mol/m^3.
c_ Used 1024.5 kg/m^3 as a constant conversion factor.
c_ Modelers can now pass tracers on a per volume basis, as carried in models.
c_
c_ Revision 1.5  1999/04/06 16:57:37  orr
c_ Changed calc of dpCO2 to account for diff atm pressure
c_
c_ Revision 1.4  1999/04/06 13:17:58  orr
c_ Added 2 output arguments: pCO2surf and dpCO2
c_ (see section 4 of Biotic HOWTO)
c_
c_ Revision 1.3  1999/04/05 15:59:11  orr
c_ Cleaned up comments regarding units
c_
c_ Revision 1.2  1999/04/04 02:35:00  orr
c_ Changed units for input and output tracers:
c_ previously in umol/kg; now in mol/T
c_ Can also pass input and output in mol/m^3 with little error.
c_
c_ Revision 1.1  1999/04/03  21:59:56  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      subroutine co2calc(t,s,dic_in,ta_in,pt_in,sit_in
     &                  ,phlo,phhi,ph,xco2_in,atmpres
     &                  ,co2star,dco2star,pCO2surf,dpco2)
C
C-------------------------------------------------------------------------
C SUBROUTINE CO2CALC
C
C PURPOSE
C	Calculate delta co2* from total alkalinity and total CO2 at
C temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
C
C USAGE
C       call co2calc(t,s,dic_in,ta_in,pt_in,sit_in
C    &                  ,phlo,phhi,ph,xco2_in,atmpres
C    &                  ,co2star,dco2star,pCO2surf,dpco2)
C
C INPUT
C	dic_in = total inorganic carbon (mol/m^3) 
C                where 1 T = 1 metric ton = 1000 kg
C	ta_in  = total alkalinity (eq/m^3) 
C	pt_in  = inorganic phosphate (mol/m^3) 
C	sit_in = inorganic silicate (mol/m^3) 
C	t      = temperature (degrees C)
C	s      = salinity (PSU)
C	phlo   = lower limit of pH range
C	phhi   = upper limit of pH range
C	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
C	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
C
C       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
C             used to initialize variables dic, ta, pt, sit, and xco2.
C             * Variables dic, ta, pt, and sit are in the common block 
C               "species".
C             * Variable xco2 is a local variable.
C             * Variables with "_in" suffix have different units 
C               than those without.

C OUTPUT
C	co2star  = CO2*water (mol/m^3)
C	dco2star = delta CO2 (mol/m^3)
c       pco2surf = oceanic pCO2 (ppmv)
c       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
C
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
c     - Models carry tracers in mol/m^3 (on a per volume basis)
c     - Conversely, this routine, which was written by observationalists 
c       (C. Sabine and R. Key), passes input arguments in umol/kg  
c       (i.e., on a per mass basis)
c     - I have changed things slightly so that input arguments are in mol/m^3,
c     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
c       should be given in mol/m^3; output arguments "co2star" and "dco2star"  
c       are likewise in mol/m^3.

C FILES and PROGRAMS NEEDED
C	drtsafe
C	ta_iter_1
C
C--------------------------------------------------------------------------
C
	real invtk,is,is2
        real k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,kcal,karag
	common /const/k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,k0,
	1	kcal,karag,ff,htotal
	common /species/bt,st,ft,sit,pt,dic,ta,calcium
	external ta_iter_1
C

c       ---------------------------------------------------------------------
C       Change units from the input of mol/m^3 -> mol/kg:
c       (1 mol/m^3)  x (1 m^3/1024.5 kg)
c       where the ocean's mean surface density is 1024.5 kg/m^3
c       Note: mol/kg are actually what the body of this routine uses 
c       for calculations.  
c       ---------------------------------------------------------------------
	permil = 1.0 / 1024.5
c       To convert input in mol/m^3 -> mol/kg 
	pt=pt_in*permil
	sit=sit_in*permil
	ta=ta_in*permil
	dic=dic_in*permil

c       ---------------------------------------------------------------------
C       Change units from uatm to atm. That is, atm is what the body of 
c       this routine uses for calculations.
c       ---------------------------------------------------------------------
	permeg=1.e-6
c       To convert input in uatm -> atm
	xco2=xco2_in*permeg
c       ---------------------------------------------------------------------
C
C*************************************************************************
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code. 
C Once calculated, the constants are
C stored and passed in the common block "const". The original version of this
C code was based on the code by Dickson in Version 2 of "Handbook of Methods
C for the Analysis of the Various Parameters of the Carbon Dioxide System
C in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
C
C Derive simple terms used more than once
C
	tk = 273.15 + t
	tk100 = tk/100.0
	tk1002=tk100*tk100
	invtk=1.0/tk
	dlogtk=log(tk)
	is=19.924*s/(1000.-1.005*s)
	is2=is*is
	sqrtis=sqrt(is)
	s2=s*s
	sqrts=sqrt(s)
	s15=s**1.5
	scl=s/1.80655
C
C f = k0(1-pH2O)*correction term for non-ideality
C
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
C
	ff = exp(-162.8301 + 218.2968/tk100  +
     &		90.9241*log(tk100) - 1.47696*tk1002 +
     &		s * (.025695 - .025225*tk100 + 
     &		0.0049867*tk1002))
C
C K0 from Weiss 1974
C
	k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +
     &		s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

C
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
C
	k1=10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
     &		0.0118 * s + 0.000116*s2))
C
	k2=10**(-1*(1394.7*invtk + 4.777 - 
     &		0.0184*s + 0.000118*s2))
C
C kb = [H][BO2]/[HBO2]
C
C Millero p.669 (1995) using data from Dickson (1990)
C
	kb=exp((-8966.90 - 2890.53*sqrts - 77.942*s +
     &		1.728*s15 - 0.0996*s2)*invtk +
     &		(148.0248 + 137.1942*sqrts + 1.62142*s) +
     &		(-24.4344 - 25.085*sqrts - 0.2474*s) *
     &		dlogtk + 0.053105*sqrts*tk)
C
C k1p = [H][H2PO4]/[H3PO4]
C
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
C
	k1p = exp(-4576.752*invtk + 115.525 - 18.453 * dlogtk +
     &		(-106.736*invtk + 0.69171) * sqrts +
     &		(-0.65643*invtk - 0.01844) * s)
C
C k2p = [H][HPO4]/[H2PO4]
C
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
C
	k2p = exp(-8814.715*invtk + 172.0883 - 27.927 * dlogtk +
     &		(-160.340*invtk + 1.3566) * sqrts +
     &		(0.37335*invtk - 0.05778) * s)
C
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
C
C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
C
	k3p = exp(-3070.75*invtk - 18.141 + 
     &		(17.27039*invtk + 2.81197) *
     &		sqrts + (-44.99486*invtk - 0.09984) * s)
C
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C
C Millero p.671 (1995) using data from Yao and Millero (1995)
C
	ksi = exp(-8904.2*invtk + 117.385 - 19.334 * dlogtk +
     &		(-458.79*invtk + 3.5913) * sqrtis +
     &		(188.74*invtk - 1.5998) * is +
     &		(-12.1652*invtk + 0.07871) * is2 +
     &		log(1.0-0.001005*s))
C
C------------------------------------------------------------------------
C kw = [H][OH]
C
C Millero p.670 (1995) using composite data
C
	kw = exp(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk +
     &		(118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     &		sqrts - 0.01615 * s)
C
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4]
C
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C
	ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     &		(-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
     &		(35474*invtk - 771.54 + 114.723*dlogtk) * is -
     &		2698*invtk*is**1.5 + 1776*invtk*is2 +
     &		log(1.0 - 0.001005*s))
C
C------------------------------------------------------------------------
C kf = [H][F]/[HF]
C
C Dickson and Riley (1979) -- change pH scale to total
C
	kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     &		log(1.0 - 0.001005*s) + 
     &		log(1.0 + (0.1400/96.062)*(scl)/ks))
C
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, calcium and fluoride
C
C Uppstrom (1974)
	bt = 0.000232 * scl/10.811
C Morris & Riley (1966)
	st = 0.14 * scl/96.062
C Riley (1965)
	ft = 0.000067 * scl/18.9984
C*************************************************************************
C
C Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
C The solution converges to err of xacc. The solution must be within
C the range x1 to x2.
C
C If DIC and TA are known then either a root finding or iterative method
C must be used to calculate htotal. In this case we use the Newton-Raphson
C "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
C error trapping removed).
C
C As currently set, this procedure iterates about 12 times. The x1 and x2
C values set below will accomodate ANY oceanographic values. If an initial
C guess of the pH is known, then the number of iterations can be reduced to
C about 5 by narrowing the gap between x1 and x2. It is recommended that
C the first few time steps be run with x1 and x2 set as below. After that,
C set x1 and x2 to the previous value of the pH +/- ~0.5. The current
C setting of xacc will result in co2star accurate to 3 significant figures
C (xx.y). Making xacc bigger will result in faster convergence also, but this
C is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
C
C Parentheses added around negative exponents (Keith Lindsay)
C
	x1 = 10.0**(-phhi)
	x2 = 10.0**(-phlo)
	xacc = 1.e-10
	htotal = drtsafe(ta_iter_1,x1,x2,xacc)
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
	htotal2=htotal*htotal
	co2star=dic*htotal2/(htotal2 + k1*htotal + k1*k2)
	co2starair=xco2*ff*atmpres
	dco2star=co2starair-co2star
	ph=-log10(htotal)

c
c       ---------------------------------------------------------------
cc      Add two output arguments for storing pCO2surf
cc      Should we be using K0 or ff for the solubility here?
c       ---------------------------------------------------------------
        pCO2surf = co2star / ff
        dpCO2    = pCO2surf - xco2*atmpres
C
C  Convert units of output arguments
c      Note: co2star and dco2star are calculated in mol/kg within this routine 
c      Thus Convert now from mol/kg -> mol/m^3
       co2star  = co2star / permil
       dco2star = dco2star / permil
	co2starair= co2starair/permil

c      Note: pCO2surf and dpCO2 are calculated in atm above. 
c      Thus convert now to uatm
       pCO2surf = pCO2surf / permeg
       dpCO2    = dpCO2 / permeg
!	print *,'xxx',co2star,co2starair,dco2star

C
	return
	end
c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /cs/home/csdmr/csrjm/mark2/code/RCS/ocmip_chem.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 2000/01/04 03:43:01 $   ;  $State: Exp $
c_ $Author: csrjm $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: ocmip_chem.f,v $
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c_ Revision 1.1  1999/03/22 12:57:48  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      function scco2(t)
c
c*********************************************************************
c
c  Computes the Schmidt number of CO2 in seawater using the 
c  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
c  7373-7382).  Input is temperature in deg C.
c
c*********************************************************************
c
      scco2 = 2073.1 - 125.62*t + 3.6276*t**2 - 0.043219*t**3
      if (scco2 .lt. 0 ) print*,'schmidt no',scco2,t
c
      return
      end
c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /cs/home/csdmr/csrjm/mark2/code/RCS/ocmip_chem.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 2000/01/04 03:43:01 $   ;  $State: Exp $
c_ $Author: csrjm $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: ocmip_chem.f,v $
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c_ Revision 1.2  1999/09/01 17:55:41  orr
c_ Fixed sign error in dfn/dx following remarks of C. Voelker (10/Aug/1999)
c_
c_ Revision 1.1  1999/04/03 22:00:42  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      subroutine ta_iter_1(x,fn,df)
      real k12,k12p,k123p
      real k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      common /const/k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,k0,
     +        kcal,karag,ff,htotal
      common /species/bt,st,ft,sit,pt,dic,ta,calcium
C
C This routine expresses TA as a function of DIC, htotal and constants.
C It also calculates the derivative of this function with respect to 
C htotal. It is used in the iterative solution for htotal. In the call
C "x" is the input value for htotal, "fn" is the calculated value for TA
C and "df" is the value for dTA/dhtotal
C
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0 + st/ks
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0*x2 + 2.0*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0*x + k1
C
C	fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
C
      fn = k1*x*dic/b +
     &     2.0*dic*k12/b +
     &     bt/(1.0 + x/kb) +
     &     kw/x +
     &     pt*k12p*x/a +
     &     2.0*pt*k123p/a +
     &     sit/(1.0 + x/ksi) -
     &     x/c -
     &     st/(1.0 + ks/x/c) -
     &     ft/(1.0 + kf/x) -
     &     pt*x3/a -
     &     ta
C
C	df = dfn/dx
C
      df = ((k1*dic*b) - k1*x*dic*db)/b2 -
     +     2.0*dic*k12*db/b2 -
     +     bt/kb/(1.0+x/kb)**2. - kw/x2 +
     +     (pt*k12p*(a - x*da))/a2 - 2.0*pt*k123p*da/a2 -
     +     sit/ksi/(1.0+x/ksi)**2. - 1.0/c +
     +     st*(1.0 + ks/x/c)**(-2.0)*(ks/c/x2) +
     +     ft*(1.0 + kf/x)**(-2.)*kf/x2 -
     +     pt*x2*(3.0*a-x*da)/a2

      return
      end
c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /cs/home/csdmr/csrjm/mark2/code/RCS/ocmip_chem.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 2000/01/04 03:43:01 $   ;  $State: Exp $
c_ $Author: csrjm $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: ocmip_chem.f,v $
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c Revision 1.1  2000/01/04  03:43:01  csrjm
c Initial revision
c
c_ Revision 1.1  1999/04/03 22:00:42  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      REAL FUNCTION DRTSAFE(FUNCD,X1,X2,XACC)
C
C	File taken from Numerical Recipes. Modified  R.M.Key 4/94
C
      MAXIT=100
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      IF(FL .LT. 0.0) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(DRTSAFE,F,DF)
      DO 100, J=1,MAXIT
        IF(((DRTSAFE-XH)*DF-F)*((DRTSAFE-XL)*DF-F) .GE. 0.0 .OR.
     &	      ABS(2.0*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          DRTSAFE=XL+DX
          IF(XL .EQ. DRTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE
          DRTSAFE=DRTSAFE-DX
          IF(TEMP .EQ. DRTSAFE)RETURN
	END IF
        IF(ABS(DX) .LT. XACC)RETURN
        CALL FUNCD(DRTSAFE,F,DF)
        IF(F .LT. 0.0) THEN
          XL=DRTSAFE
          FL=F
        ELSE
          XH=DRTSAFE
          FH=F
        END IF
  100  CONTINUE
      RETURN
      END
c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/orr/WWW/Biotic/boundcond/RCS/o2sato.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 1999/03/22 12:57:48 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: o2sato.f,v $
c_ Revision 1.1  1999/03/22 12:57:48  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      function o2sato(T,S)
c
C ********************************************************************
C                                     
C Computes the oxygen saturation concentration at 1 atm total pressure
c in mol/m^3 given the temperature (t, in deg C) and the salinity (s,
c in permil). 
C
C FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
C THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
c
C *** NOTE: THE "A3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
C *** IT SHOULDN'T BE THERE.                                ***
C
C o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
c 0 permil <= S <= 42 permil
C C
C CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
c o2sato = 0.282015 mol/m^3
C
C ********************************************************************
c
      DATA A0/ 2.00907   /,A1/ 3.22014   /, A2/ 4.05010 /,
     $     A3/ 4.94457   /,A4/-2.56847E-1/, A5/ 3.88767 /
      DATA B0/-6.24523E-3/,B1/-7.37614E-3/,
     $     B2/-1.03410E-2/,B3/-8.17083E-3/
      DATA C0/-4.88682E-7/
C
      TT  = 298.15-T
      TK  = 273.15+T
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5
     $     + S*(B0 + B1*TS + B2*TS2 + B3*TS3)
     $     + C0*(S*S)
      o2sato = EXP(CO)
c
c  Convert from ml/l to mol/m^3
c
      o2sato = o2sato/22391.6*1000.0
C
      RETURN
      END
c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/orr/WWW/Biotic/boundcond/RCS/o2flux.f,v $ 
c_ $Revision: 1.2 $
c_ $Date: 1999/04/04 02:43:12 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: o2flux.f,v $
c_ Revision 1.2  1999/04/04 02:43:12  orr
c_ Fixed calc of )2 piston velocity.
c_
c_ Revision 1.1  1999/03/22 12:57:48  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      subroutine o2flux(t,s,kw660,ppo,o2,dz,o2ex)
c
c**********************************************************************
c
c  Computes the time rate of change of oxygen in the surface
c  layer due to air-sea gas exchange in mol/m^/s.
c
c  Inputs:
c    t       model surface temperature (deg C)
c    s       model surface salinity (permil)
c    kw660   gas transfer velocity at a Schmidt number of 660, accounting
c              for sea ice fraction (m/s)
c    ppo     surface pressure divided by 1 atm.
c    o2      surface ocean O2 concentration (mol/m^3)
c    dz      thickness of surface grid box (m)
c  Output:
c    o2ex    time rate of change of oxygen in the surface layer due
c              to air-sea exchange (mol/m^3/s)
c
c  Two functions are called:
c    sco2    Schmidt number of oxygen
c    o2sato  oxygen saturation concentration at 1 atm (mol/m^3)
c
c  Numbers in brackets refer to equation numbers in simulation design
c  document.
c
c  Ray Najjar 1/29/99
c
c**********************************************************************
c
      real kwo2,kw660
c
c  Compute the transfer velocity for O2 in m/s [4]
c
cc    kwo2 = (sco2(t)/660)**-0.5*0.01/3600.0
      kwo2 = Kw660 * (660/sco2(t))**0.5
c
c  Compute the saturation concentrations for O2 [3]
c
      o2sat = o2sato(t,s)*ppo
c
c  Compute time rate of change of O2 due to gas exchange [1]
c
      o2ex = kwo2*(o2sat-o2)/dz
c
      return
      end
c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/orr/WWW/Biotic/boundcond/RCS/sco2.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 1999/03/22 12:57:48 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: sco2.f,v $
c_ Revision 1.1  1999/03/22 12:57:48  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      function sco2(t)
c
c*********************************************************************
c
c  Computes the Schmidt number of oxygen in seawater using the
c  formulation proposed by Keeling et al. (1998, Global Biogeochem.
c  Cycles, 12, 141-163).  Input is temperature in deg C.
c
c*********************************************************************
c
      sco2 = 1638.0 - 81.83*t + 1.483*t**2 - 0.008004*t**3
c
      return
      end
