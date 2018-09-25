c  Routines called by biogeochemical module
c =====================================================================
c  Function to calculate oxygen saturation in mmol/m3
      
      
      function o2sat(t,s)

! THE FOLLOWING DOCUMENTATION IS TAKEN DIRECTLY FROM MOCSY PACKAGE *pjb*

!    Purpose:
!    Compute O2 saturation concentration of surface seawater
!    (mol/m3) at 1 atm (Garcia & Gordon, 1992)
!
!    ********************************************************************
!    Computes the oxygen saturation concentration at 1 atm total pressure
!    in mol/m^3 given sea surface temperature T (deg C) and salinity S (permil) 
!
!    From: Garcia & Gordon (1992) Oxygen solubility in
!    seawater: better fitting equations,
!          Limnol. Oceanogr., 37(6), 1307-1312.
!          This routine uses:
!          - equation (8) on page 1310
!          - coefficients from Table 1, column 2 (Benson & Krause, [cm3/dm3], i.e, same as [ml/L])
!
!    *** NOTE: The "A3*ts^2" term in equation (8) in the paper is a TYPO.
!    *** It shouldn't be there. It is not used in this routine.
!
!    'o2sat' is fit between T(freezing) <= T <= 40(deg C) and  0 <= S <= 42 permil
!
!    CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, o2sat = 282.015  mmol/m^3
!    ********************************************************************

      REAL :: a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, c0,
     +        tt, tk, ts, ts2, ts3, ts4, ts5, o2sat
      DATA a0/ 2.00907 /, a1/ 3.22014     /, a2/ 4.05010 /,  
     +     a3/ 4.94457 /, a4/ -2.56847e-1 /, a5/ 3.88767 /
      DATA b0/ -6.24523e-3 /, b1/ -7.37614e-3 /, b2/ -1.03410e-2 /, 
     +     b3/ -8.17083e-3 /
      DATA c0/ -4.88682e-7 /
      
      tt = 298.15 - t
      tk = 273.15 + t
      ts = alog(tt/tk)

      ts2 = ts**2
      ts3 = ts**3
      ts4 = ts**4
      ts5 = ts**5

      o2sat = exp(a0 + a1*ts + a2*ts2 + a3*ts3 + a4*ts4 + a5*ts5 +
     +            s*(b0 + b1*ts + b2*ts2 + b3*ts3) + c0*(s*s)) 
     +        / 22391.6*1e6 ! convert from ml/L to mmol/m3

      return
      end

c =====================================================================

c =====================================================================
c	Function to calculation the fCO2atm value from T,Patm,depth,pCO2
c
      function fco2(pco2, temp, Patm, dep)
      
      IMPLICIT NONE

      REAL :: pco2, temp, Patm, dep
      REAL :: tk, prb, Ptot, Rgas_atm, B, Del, xCO2approx, xc2,
     +        fugcoeff, fco2

        tk = 273.15 + temp
        prb = dep/10.0
        Ptot = Patm + prb/1.01325      !Total pressure (atmospheric + hydrostatic) [atm]
        Rgas_atm = 82.05736            !R in (cm3 * atm) / (mol * K) from CODATA (2006)
      
        ! To compute fugcoeff, we need 3 other terms (B, Del, xc2)
        ! as well as 3 others above (tk, Ptot, Rgas_atm)

        B = -1636.75 + 12.0408*tk - 0.0327957*(tk*tk) +
     +        0.0000316528*(tk*tk*tk)
        Del = 57.7 - 0.118*tk

        ! "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
        ! x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
        ! Let's assume that xCO2 = pCO2. Resulting fugcoeff is identical to 8th digit after the decimal.


        xCO2approx = pco2 * 1.e-6
        xc2 = (1.0 - xCO2approx)**2
        fugcoeff = exp( Ptot*(B + 2.0*xc2*Del)/(Rgas_atm*tk) )
        fco2 = pco2 * fugcoeff

      return
      end
C=======================================================================


c =====================================================================
c	Function to calculation the pCO2 value from T,S,Alk, Tco2
c
	function pco2(t,s,ta,tc,d0,iflag)
	save ah
	tsi=0.
	tp=0.
	s=max(0.0001,s)
	call peng(t+273.15,s )
*	call unesco51(t+273.15,s)
        call tcta2pc(s,t,ta,tc,tp,tsi,
     1               ah,pc,co,c1,c2,d0,iflag)
*	write(6,*) 't tc ta ph pco2 H2CO2 HCO3 CO3'
*	write(6,'(8f10.3)') t,tc,ta,-alog10(ah),pc,co,c1,c2
	pco2=pc
	return
	end
C=======================================================================
C
C   PURPOSE: 
C      To estimate the concentration of the carbonate species:
C         -  aH  (Activity of Hydrogen ion)
C         -  pCO2
C         -  [CO2] + [H2CO3]
C         -  [HCO3-]
C         -  [CO3=]
C      in seawater from measurements of:
C         -  TOTAL ALKALINITY
C         -  TOTAL CO2
C         -  TOTAL PHOSPHATE
C         -  TOTAL SILICATE
C
C
C   PROGRAMMER: Y.-H. CHAN  (original code)
C               J.Page      (minor modifications)
C   
C   COMPUTER:   IOS VAX 6410
C   LANGUAGE:   VAX FORTRAN  (Enhanced ANSI-77)
C
C
C   REFERENCES:
C     1. Peng, T.-H., T. Takahashi, W.S. Broecker & J. Olafsson (1987)          <= primary reference
C        SEASONAL VARIABILITY OF CARBON DIOXIDE, NUTRIENTS AND OXYGEN
C        IN THE NORTHERN  ATLANTIC SURFACE WATER: OBSERVATIONS AND A MODEL.
C        TELLUS 39B,  439--458.
C     2. Catherine Goyet and Alain Poisson                                      <= not applied as yet
C        New determination of carbonic acid dissociation constants in 
C        seawater as a function of temperature and salinity.
C        Deep-Sea Research, Vol 36, No. 11, pp1635-1654, 1989
C
C-----------------------------------------------------------------------
       subroutine tcta2pc( sal,swt,talk,tco2,tp,tsi,       !input
     &   ah,pco2,cca,hco3,co3,d0,iflag)       ! output
      
      implicit none
      integer n,iflag
      real  pco2         !Partial pressure of CO2 in seawater (uatm)
      real  cca          !Conc of CO2+H2CO3  (uM/kg)
      real  hco3         !Conc of Bicarbonate (uM/kg)
      real  co3          !Conc of Carbonate   (uM/kg)
      real  sal          !Salinity
      real  swt          !Sea Water Temperature (C)
      real  tp           !Total Phosphate in seawater  (uM/kg)
      real  tsi          !Total Silicate in seawater   (uM/kg)
      real  talk         !Total Alkalinity  (uEq/kg)
      real  tco2         !Total CO2 in seawater  (uM/kg)
      real  ah           !Activity of Hydrogen ion
      real  ph           !pH of seawater
      real  d0           !see kc0	
      real  kc0          !1st apparent dissociation constant of Carbonic Acid
      REAL  KC1          !2nd ...
      REAL  KC2          !3rd ...
      REAL  KB           !First dissociation constant of Boric Acid
      REAL  KP1          !1st apparent dissociation constant of Phosphoric Acid
      REAL  KP2
      REAL  KP3
      REAL  KSI          !First dissociation constant of Silicic Acid
      REAL  KW           !Dissociation constant of water
      REAL  FH
      REAL  S            !Salinity
      REAL  T            !Sea water temp (C)
      REAL  TA           !Total Alkalinity
      REAL  TB           !Total Borate
      REAL  TC           !Total CO2
      REAL  TK           !Sea water temp (K)
      REAL  TK2
      REAL  AA
      REAL  AC           !Carbonate Alkalinity
      REAL  AHN
      REAL  AHI
      REAL  AB           !Borate Alkalinity
      REAL  AP           !Phosphate Alkalinity
      REAL  ASI          !Silicate Alkalinity
      REAL  AW           !Water Alkalinity
      REAL  DEL
	COMMON /KVALUES/ KC0,KC1,KC2,KB,KP1,KP2,KP3,KSI,KW,FH
      S   = (SAL)
      T   = (SWT)
      TA  = (TALK)
      TC  = (TCO2)
      TK  = T + 273.15
      TK2 = TK /100.0
      IF (S.LE.0.0 .or. TK.LE.0.0) THEN
*          PRINT 1006,SWT,SAL,TALK,TCO2
         STOP
      ENDIF
*	CALL PENG(TK,S) 
      TB = 4.106e2 * (S/35.)    !Culkin (1965)
*          PRINT *,SWT,SAL,TALK,TCO2
C      ****************************************
C      *  ITERATION TO FIND BEST VALUE OF AH  *
C      ****************************************
       DEL = 1.0
       AHN = 1.0e-8
	 iflag=0
	if (ah .ne. 0) ahn=ah
       N = 0
       DO WHILE (DEL.GE.1.0e-5 .AND. N.LT.20)
          N = N + 1
          AHI = AHN  
          AB  = TB  * (KB / (AHI + KB))
          ASI = 0 ! (TSI) * (KSI / (AHI + KSI))
          AP  = 0 ! (1.0/(1.0 + KP2/AHI + KP2*KP3/(AHI*AHI))
*     &          + 2.0/(1.0 + AHI/KP2 + KP3/AHI)
*     &          + 3.0/(1.0 + AHI/KP3 + AHI*AHI/(KP2*KP3))) * (TP)
          AW  = FH * (KW * 1.0e6 / AHI) - (AHI * 1.0e6 / FH)
          AC  = TA - AB - ASI - AP - AW
          AA  = TC - AC
          AHN = 0.5D0 * (KC1/AC)
     &          *(AA + SQRT(AA*AA + 4.0*AC*(KC2/KC1)*(2.0*TC-AC)))
          DEL = ABS(1.0 - AHI/AHN) 
       END DO
       CCA  = (2.0*TC - AC) / (2.0 + KC1/AHN)
       PCO2 = CCA / KC0 
       HCO3 = CCA * (KC1 / AHN)
       CO3  = HCO3 * (KC2 / AHN)
       AH   = (AHN)
c If the maximum number of iterations are reach assume
c pco2 is zero and reset ahn
       IF (N.GE.20 .AND. DEL.GE.1.0e-5) then 
		PRINT 1001, N, DEL
		print*,'PH = ',ah,tco2,talk,swt,sal
		iflag=-1
		pco2=0
		ah=1e-8
		ahn=1e-8
	 endif
*	  WRITE(6, 1003) AH,PH,PCO2,CCA,HCO3,CO3,AC,AB,ASI,AP,AW
	d0=kc0
 1001  FORMAT ('0CONVERGENCE FAILED IN SUBROUTINE TCTA2PC ITERATION'/ 
     & ' NO. OF ITERATION =',I3, '   DEL =',D15.8)
 1002  FORMAT (' ', 10D13.6)
 1003  FORMAT ('0', 6E15.8/' ',5e15.8/)
 1004  FORMAT ('0', 6E15.8,2e15.8/' ',8e15.8)
 1005  FORMAT('0E R R O R  in subroutine TCTA2PC()  =>  ',
     &    'value of AH outside legal range => cannot compute PH'/
     & '   SWT =',F12.3,' sea water temp   (C)'/
     & '   SAL =',F12.3,' salinity'/
     & '   TALK=',F12.3,' total Alkalinity  (uE/kg)'/
     & '   TCO2=',F12.3,' total CO2         (uM/kg)'/
     & '   TP  =',F12.3,' phosphate         (uM/kg)'/
     & '   TSI =',F12.3,' silicate          (uM/kg)'/
     & '   AH  =',F12.3,' hydrogen ion activity')
 1006  FORMAT('0E R R O R  in subroutine TCTA2PC()  =>  '/
     & '   SWT =',F12.3,' sea water temp   (C)'/
     & '   SAL =',F12.3,' salinity'/
     & '   TALK=',F12.3,' total Alkalinity  (uE/kg)'/
     & '   TCO2=',F12.3,' total CO2         (uM/kg)')
       RETURN
       END
	SUBROUTINE PENG(TK,S)
	IMPLICIT REAL  (A-H,K,O-Z)
	COMMON /KVALUES/ KC0,KC1,KC2,KB,KP1,KP2,KP3,KSI,KW,FH
      TK2 = TK /100.0D0
	T=TK-273.15D0
      KC0 = EXP( -60.2409 + 9345.17/TK + 23.3585*aLOG(TK2)          !Weiss (1974)
     &   + S*(2.3517e-2 - 2.3656e-2*TK2 + 4.7036e-3*TK2*TK2))
      KC1 = 10.0**(13.7201 - 3.1334e-2*TK - 3235.76/TK               !Mehrbach et al (1973)
     &   - 1.30D-5*S*TK + 0.1032*SQRT(S) )
      KC2 = 10.0**( -5371.9645 - 1.671221*TK + 128375.28/TK          !Mehrback et al (1973)
     &   + 2194.3055*aLOG10(TK) - 0.22913*S - 18.3802*aLOG10(S)
     &   + 8.0944e-4*S*TK + 5617.11*aLOG10(S)/TK - 2.136*S/TK)   
      KB  = 10.0**( -9.26 + 8.86e-3*S + 1.0e-2*T)                    !Lyman (1956)
C###  KP1 = 2.0e-2              !not used (negligible contrib.)      
      KP2 = EXP(-9.039 - 1450.0/TK)                                 !Kester & Pytkowicz (1967)
      KP3 = EXP( 4.466 - 7276.0/TK)                                 !Kester & Pytkowicz (1967)
      KSI = 4.0e-10                                                  !Ingri (1959)
      KW = EXP(148.9802- 13847.26/TK - 23.6521*aLOG(TK) - 0.019813*S      !Culberson & Pytkowicz (1973)
     &   + SQRT(S)*( -79.2447 + 3298.72/TK + 12.0408*aLOG(TK)))
      FH = 1.29 - 2.04e-3*TK + (4.61e-4 - 1.48e-6*TK)*S*S             !Culberson & Pytkowicz (1973) / Takahashi (1982a)
*      print *, KC0,KC1,KC2,KB,KSI,KP2,KP3,KW,FH
	RETURN
	END
	SUBROUTINE UNESCO51(TK,S)
C*  THE NEW KC1,K2,KB COMES FROM UNESCO REPORT 51
	IMPLICIT REAL   (A-H,K,O-Z)
	COMMON /KVALUES/ KC0,KC1,KC2,KB,KP1,KP2,KP3,KSI,KW,FH
        TK2 = TK /100.0e0
	T=TK-273.15e0
      KC0 = EXP( -60.2409 + 9345.17/TK + 23.3585*aLOG(TK2)          !Weiss (1974)
     &   + S*(2.3517e-2 - 2.3656e-2*TK2 + 4.7036e-3*TK2*TK2))
	KB01=148.0248 - 8966.90/TK -24.4344* aLOG(TK)
	KB= EXP ( KB01 +  (0.5998-75.25/TK)*S**.5 - 0.01767* S)
	
	KC101=6320.81/TK - 126.3405 + 19.568 * aLOG(TK)
	KC201=5143.69/TK - 90.1833 + 14.613 *aLOG(TK)
	
	kjunk= (KC101 +  (19.894-840.39/TK -3.0189*aLOG(TK))*S**.5
	1        + 0.0068*S )
	KC1= 10**(-kjunk)
	kjunk=(KC201 +   (17.176-690.59/TK- 2.6719*aLOG(TK))*S**.5
	1        + 0.0217*S )
	KC2= 10**(-kjunk)
C###  KP1 = 2.0e-2              !not used (negligible contrib.)      
      KP2 = EXP(-9.039 - 1450.0/TK)                                 !Kester & Pytkowicz (1967)
      KP3 = EXP( 4.466 - 7276.0/TK)                                 !Kester & Pytkowicz (1967)
      KSI = 4.0e-10                                                  !Ingri (1959)
      KW = EXP(148.9802- 13847.26/TK - 23.6521*aLOG(TK) - 0.019813*S      !Culberson & Pytkowicz (1973)
     &   + SQRT(S)*( -79.2447 + 3298.72/TK + 12.0408*aLOG(TK)))
      FH = 1.29 - 2.04e-3*TK + (4.61e-4 - 1.48e-6*TK)*S*S             !Culberson & Pytkowicz (1973) / Takahashi (1982a)
	
*      print *, KC0,KC1,KC2,KB,KSI,KP2,KP3,KW,FH
	RETURN
	END
* ****************************************************
*	Calculate the schmidt number for various gases 
*  in seawater (Wanninkhof 1992, JGR 97c).
*
      function  schmidt_no(igas,t)
* igas 	- the gas to calculate, 
* 1) O2, 2) CO2, 3) SF6, 4) DMS
* 5) CFC-12, 6) CFC-12, 7) N2O

* t		- temperature degrees C

      parameter (ngas=7)
      dimension a(5,ngas)      
      common /c_sch/ a
      data a /1920.4, 135.6,  5.2122, 0.10939, 0.00093777, ! O2
     1        2116.8, 136.25, 4.7353, 0.092307, 0.0007555, ! Co2
     1        3177.5, 200.57, 6.8865, 0.13335, 0.0010877,  ! SF6
     1        2855.7, 177.63, 6.0438, 0.11645, 0.00094743, ! DMS
     1        3828.1, 249.86, 8.7603, 0.1716, 0.001408,    ! CCl2F2 (CFC-12)
     1        3579.2, 222.63, 7.5749, 0.14595, 0.0011874,  ! CCl3F (CFCF-11)
     1        2356.2, 166.38, 6.3952, 0.13422, 0.0011506/  ! N2O

      schmidt_no = a(1,igas) - a(2,igas)*t + a(3,igas)*t*t 
     +              - a(4,igas)*t*t*t + a(5,igas)*t*t*t*t
      return
      end


* ****************************************************
*	Calculate the solubility of for various gases in 
*  seawater (Wanninkhof 1992, JGR 97c).
*	See original references for CO2 and CFC
	function solu_gas (igas,t,s)
* igas - the gas to calculate, 1) oxygen, 2) CO2, 3) f-12, 4) f-12
* t	 - temperature degrees C
* s	 - salinity
	parameter(ngas=4)
	dimension a(7,ngas)
	common /c_solu/ a
c2345678901234567890123456789312345678941234567895123456789612345678971234567898
	data a /-58.3877,85.8079,23.8439, -0.034892, 0.015568, -0.0019387,0.,
	1	 -58.0931,90.5069,22.2940,  0.027766, -0.025888,  0.0050578,0.,
	2 -218.0971, 298.9702, 113.8049, -1.39165, -0.143566, 0.091015, 
	2 -0.0153924,
	3 -229.9261, 319.6552, 119.4471, -1.39165, -0.142382, 0.091459, 
	3 -0.0157274/
	tk=t+273.15
	td=tk*1e-2
	tdi=1e2/tk
	if (igas .eq.1) then
	 solu_gas=0
	else if (igas.eq.2) then
	 b=a(1,igas)+ a(2,igas)*tdi + a(3,igas)*alog(td) + 
	1	s*(a(4,igas) + a(5,igas)*td + a(6,igas)*td*td )
	 solu_gas= exp(b) 
	else if (igas .gt.2) then
* includes the effect of r.h. and is converted to concentration by
* multiplying by atmospheric dry air mole fraction
	 b=a(1,igas)+ a(2,igas)*tdi + a(3,igas)*alog(td) + 
	1	a(4,igas)*td*td +
	1	s*(a(5,igas) + a(6,igas)*td + a(7,igas)*td*td )
	solu_gas= exp(b) 
	endif
	return
	end
*	Chemical enhancement of CO2 exchange
	function chem_co2_x(t)
	chem_co2_x= 2.5*(0.5246+1.6256e-2*t+4.9946e-4*t*t)
	return
	end


c =====================================================================
c  Function to calculate surface humidity (pH2O)
      
      function pH2O(t,s)

        tk = t+273.15
        pH2O = 24.4543 - 67.4509*(100./tk) - 4.8489*alog(tk/100.0) - 
     +         0.000544*s
        
      return
      end


c =====================================================================
c  Function to calculate gas saturation in umol/kg
      function phizero(t,s,igas)

c Constants come from Orr et al., 2017 Geoscientific Model Development 

      select case (igas)
      case (1)
c Oxygen
      a1 = -177.7888
      a2 =255.5907
      a3 =146.4813
      a4 =-22.2040
      b1 =-0.037362
      b2 =0.016504
      b3 =-0.0020564

      case (2)
c CO2
      a1 = -160.7333
      a2 = 215.4152
      a3 = 89.8920
      a4 = -1.47759
      b1 = 0.029941
      b2 = -0.027455
      b3 = 0.0053407
      
      case (3)
c SF6
      a1 = -80.0343
      a2 = 117.232
      a3 = 29.5817
      a4 = 0.0
      b1 = 0.0335183
      b2 = -0.0373942
      b3 = 0.00774862

      case (5)
c CFC-12
      a1 = -218.0971
      a2 = 298.9702
      a3 = 113.8049
      a4 = -1.39165
      b1 = -0.143566
      b2 = 0.091015
      b3 = -0.0153924

      case (6)
c CFC-11
      a1 = -229.9261
      a2 = 319.6552
      a3 = 119.4471
      a4 = -1.39165
      b1 = -0.142382
      b2 = 0.091459
      b3 = -0.0157274

      case (7)
c N2O
      a1 = -165.8806 
      a2 = 222.8743 
      a3 = 92.0792  
      a4 = -1.48425 
      b1 = -0.056235 
      b2 = 0.031619 
      b3 = -0.0048472

      case default
       print*, ' no gas defined for igas=',igas
 
      end select

      t100=(t+273.15)/100.
      
      ! Following gives solubility function (mol/L/atm)
      phizero = exp(a1 + a2/t100 + a3*alog(t100) + a4*t100*t100 + 
     +              s*(b1 + b2*t100 + b3*t100*t100))
      
      phizero = phizero*1e6 !mol/L --> mmol/m3

      return
      end


      real function kquad(c,t)
      real  c(3),r,t
      r=10**(-3.0)
      kquad=(c(1) + c(2)*t + c(3)* t**2.0 * r)*r
      return
      end

      subroutine pmolvol(t,s,p)
	include "co2.h"
      real dvfresh(3,13),dvsalt(3,13), dkfresh(3,13),dksalt(
     *3,13), kp(13)
      real vquad,kquad
      equivalence (kp,k1)

      data dvfresh/ -30.54,0.1849,-2.3366, -29.81,0.1150,-1.816, -25.60,
     *0.2324,-3.6246, -38.75,0.1991,-2.813, -23.10,0.0990,-1.560, -12.81
     *,-0.0028,-1.254, -19.48,0.1991,-2.813, -27.59,0.1150,-1.816, -38.3
     *8,0.1600,-2.526, 0.0,0.0,0.0, 32.3,0.0,0.0, -65.28,0.397,-5.155, -
     *62.50, 0.397,-5.155/
      data dvsalt/ -26.69,0.1976,-3.111, -17.32,0.1758,-2.647, -20.02,0.
     *1119,-1.409, -28.56,0.1211,-0.321, -18.03,0.0466,0.316, -9.78,-0.0
     *09,-0.942, -14.51,0.1211,-0.321, -23.12,0.1758,-2.647, -26.57,0.20
     *20,-3.042, 0.0,0.0,0.0, 32.3,0.0,0.0, -45.464,0.3529,-4.985, -42.6
     *80,0.3529,-4.985/
      data dkfresh/ -6.22,0.1368,-1.233, -5.74,0.093,-1.896, -7.33,0.136
     *8,-1.233, -10.13,0.2396,-3.961, -5.39,0.093,-1.896, -5.57,0.131,-1
     *.160, -6.51,0.2396,-3.961, -5.24,0.0930,-1.896, -10.10,0.1792,-3.6
     *53, 0.0,0.0,0.0, 0.0,0.0,0.0, 18.47,0.1956,-2.212, 18.47,0.1956,-2
     *.212/
      data dksalt/ -3.90,0.1700,0.0, -4.87,0.0900,0.0, -5.13,0.0794,0.0,
     * -3.0,0.0427,0.0, -4.53,0.0900,0.0, -3.91,0.054,0.0, -2.67,0.0427,
     *0.0, -5.15,0.0900,0.0, -4.08,0.0714,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0,
     * -13.70,0.1245,0.0, -13.70,0.1245,0.0/
      pbar=p/10.0
      rt=(273.15+t)*83.145
      saltrat=sqrt(s)/sqrt(35.0)
      i=1
23212 if(.not.(i.le.13))goto 23214
      dvf = vquad(dvfresh(1,i),t)
      dkf = kquad(dkfresh(1,i),t)
      dv = dvf + (vquad(dvsalt(1,i),t) - dvf) * saltrat
      dk = dkf + (kquad(dksalt(1,i),t) - dkf) * saltrat
      kp(i) = kp(i)*exp(-dv*pbar/rt + 0.5*dk*pbar**2.0/rt)
23213 i=i+1
      goto 23212
23214 continue
      return
      end

	function vquad(c,t)
	real c(3),r,t
	r=10**(-3.0)
	vquad=c(1) + c(2)*t + c(3)* t**2.0 * r
	return
	end



C -------------------------------------------------------------------- C

      FUNCTION rho(salt, temp, pbar)

C     Compute in-situ density from salinity (psu), in situ temperature
C     (C), & pressure (bar)

          REAL :: s, t, p, rhow, rho0 
          REAL :: a, b, c 
          REAL :: Ksbmw, Ksbm0, Ksbm
          REAL :: rho

          s = salt
          t = temp
          p = pbar

          ! Density of pure water
          rhow = 999.842594 + 6.793952e-2*t          
     +           - 9.095290e-3*t*t + 1.001685e-4*t**3  
     +           - 1.120083e-6*t**4 + 6.536332e-9*t**5

          a = 8.24493e-1 - 4.0899e-3*t + 7.6438e-5*t*t 
     +        - 8.2467e-7*t**3 + 5.3875e-9*t**4
          b= -5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t*t
          c = 4.8314e-4

          ! Density of seawater at 1 atm, P=0
          rho0 = rhow + a*s + b*s*SQRT(s) + c*s**2.0


          ! Secant bulk modulus of pure water
          ! The secant bulk modulus is the average change in pressure
          ! divided by the total change in volume per unit of initial
          ! volume.
          Ksbmw = 19652.21 + 148.4206*t - 2.327105*t*t 
     +            + 1.360477e-2*t**3 - 5.155288e-5*t**4

          ! secant bulk modulus of seawater at 1 atm
          Ksbm0 = Ksbmw + s*( 54.6746 - 0.603459*t + 
     +            1.09987e-2*t**2  - 6.1670e-5*t**3) +
     +            s*SQRT(s)*( 7.944e-2 + 1.6483e-2*t - 5.3009e-4*t**2)

          ! secant bulk modulus of seawater at s,T,p
          Ksbm = Ksbm0 + p*(3.239908 + 1.43713e-3*t +
     +           1.16092e-4*t**2 - 5.77905e-7*t**3) + 
     +           p*s*(2.2838e-3 - 1.0981e-5*t - 1.6078e-6*t**2) +
     +           p*s*SQRT(s)*1.91075e-4 + 
     +           p*p*(8.50935e-5 - 6.12293e-6*t + 5.2787e-8*t**2) +
     +           p*p*s*(-9.9348e-7 + 2.0816e-8*t + 9.1697e-10*t**2)

          ! Density of seawater at S,T,P
          rho = rho0/(1.0 - p/Ksbm)

      RETURN
      END FUNCTION rho

C -------------------------------------------------------------------- C



C -------------------------------------------------------------------- C
C    THIS IS THE OLD CARBONATE CHEMISTRY CALCUALTION                   C

      subroutine co2calc1(t,s,dicin,talkin,po4in,siin,
     +                    phlo,phhi,ph,xco2in,atmpres,co2star,dco2star,
     +                    omegaar,omegaca,co3,fco2,p,ff_out)



C-------------------------------------------------------------------------
C SUBROUTINE CO2CALC
C
C PURPOSE
C	Calculate delta co2* from total alkalinity and total CO2 at
C temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
C
C USAGE
C	call co2calc(t,s,phlo,phhi,ph,xco2,atmpres,co2star,dco2star)
C
C INPUT
C	dic = total inorganic carbon (umol/kg) [common]
C	ta  = total alkalinity (umol/kg) [common]
C	pt  = inorganic phosphate (umol/kg) [common]
C	sit = inorganic silicate (umol/kg) [common]
C	t   = temperature (degrees C)
C	s   = salinity (PSU)
C	phlo= lower limit of pH range
C	phhi= upper limit of pH range
C	xco2=atmospheric mole fraction CO2 in dry air (ppmv) 
C	atmpres = atmospheric pressure in atmospheres (1 atm==1013.25mbar)
C OUTPUT
C	dco2star = delta CO2 ()
C
C FILES and PROGRAMS NEEDED
C	drtsafe
C	ta_iter_1
C
C--------------------------------------------------------------------------
C
      include "co2.h"
      external ta_iter_1
C
C Fix units from the input of umol/kg to mol/kg
C
c       To convert input in m mol/m^3 -> mol/kg 
      permil = 1.0 / 1024.5*1e-3
!	permil=1.e-6
      pt=po4in*permil
      sit=siin*permil
      ta=talkin*permil
      dic=dicin*permil

c       To convert input in uatm -> atm
      xco2=xco2in*1e-6
C
C*************************************************************************
C	Calculate all constants needed to convert between various measured
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
      is=19.924*s/(1000.-1.005*s) ! Ionic strength (see Millero 1995)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5
      scl=s/1.80655
C
C f = k0(1-pH2O)*correction term for non-ideality
C
C Solubility of atmospheric CO2       
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
C
      ff = exp( -162.8301 + 218.2968/tk100 + 90.9241*log(tk100) - 
     +          1.47696*tk1002 +
     +          s*(.025695 - .025225*tk100 + 0.0049867*tk1002) )
C
C K0 from Weiss 1974, where K0 is moles of CO2 per kg per atm 
C       That is, the solubility of CO2 in SEAWATER     
C
      k0 = exp(-60.2409 + 93.4517/tk100 + 23.3585*log(tk100) +
     +         s * (0.023517 - 0.023656*tk100 + 0.0047036*tk1002))

C
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
C
      k1 = 10**(-1 * (-62.008 + 3670.7*invtk + 9.7944*dlogtk -
     +                0.0118*s + 0.000116*s2) )
C
      k2 = 10**(-1 * (4.777 + 1394.7*invtk - 
     +                0.0184*s + 0.000118*s2))
C
C kb = [H][B(OH)4]/[B(OH)3]
C
C Millero p.669 (1995) using data from Dickson (1990)
C
      kb = exp( (-8966.90 - 2890.53*sqrts - 77.942*s +
     &           1.728*s15 - 0.0996*s2)*invtk +
     &          (148.0248 + 137.1942*sqrts + 1.62142*s) +
     &          (-24.4344 - 25.085*sqrts - 0.2474*s)*dlogtk
     &          + 0.053105*sqrts*tk)
C
C k1p = [H][H2PO4]/[H3PO4]
C
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
C
      k1p = exp(-4576.752*invtk + 115.525 - 18.453*dlogtk +
     &          (-106.736*invtk + 0.69171)*sqrts +
     &          (-0.65643*invtk - 0.01844)*s)
C
C k2p = [H][HPO4]/[H2PO4]
C
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
C
      k2p = exp(-8814.715*invtk + 172.0883 - 27.927*dlogtk +
     &          (-160.340*invtk + 1.3566)*sqrts +
     &          (0.37335*invtk - 0.05778)*s)
C
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
C
C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
C
      k3p = exp(-3070.75*invtk - 18.141 + 
     &          (17.27039*invtk + 2.81197)*sqrts +
     &          (-44.99486*invtk - 0.09984)*s)
C
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C
C Millero p.671 (1995) using data from Yao and Millero (1995)
C
      ksi = exp(-8904.2*invtk + 117.385 - 19.334 * dlogtk +
     &          (-458.79*invtk + 3.5913) * sqrtis +
     &          (188.74*invtk - 1.5998) * is +
     &          (-12.1652*invtk + 0.07871) * is2 +
     &          log(1.0-0.001005*s))
C
C------------------------------------------------------------------------
C kw = [H][OH]
C
C Millero p.670 (1995) using composite data
C
      kw = exp(-13847.26*invtk + 148.9652 - 23.6521*dlogtk +
     &         (118.67*invtk - 5.977 + 1.0495*dlogtk)*sqrts -
     &         0.01615*s)
C
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4]
C
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C
      ks = exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     &       (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
     &       (35474*invtk - 771.54 + 114.723*dlogtk) * is -
     &       2698*invtk*is**1.5 + 1776*invtk*is2 +
     &       log(1.0 - 0.001005*s))
C
C------------------------------------------------------------------------
C kf = [H][F]/[HF]
C
C Dickson and Riley (1979) -- change pH scale to total
C
      kf = exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     &         log(1.0 - 0.001005*s) + 
     &         log(1.0 + (0.1400/96.062)*(scl)/ks))
C


C kcal and karagonite
c ksp0cal = gamma*[Ca]gamma*[CO3]/[CaCO3 calcite]
c  ksp0arag = gamma*[Ca]gamma*[CO3]/[CaCO3 aragonite]
c  Mucci (1983)


      lnksp0cal = -395.8293 + 6537.773/tk + 71.595*log(tk) - 
     +            0.17959*tk
      lnksp0arag = -395.9180 + 6685.079/tk + 71.595*log(tk) -
     &             0.17959*tk

c  kcal = [Ca][CO3]/[CaCO3 calcite]
      kcal = exp(lnksp0cal +
     +           (-1.78938 + 410.64/tk + .0065453 * tk) * sqrt(s) -
     +            0.17755 * s + .0094979 * s**1.5)

c  karag = [Ca][CO3]/[CaCO3 aragonite]
      karag = exp(lnksp0arag +
     &            (-.157481 + 202.938/tk + .0039780 * tk) * sqrt(s) -
     &             0.23067 * s + .0136808 * s**1.5)

* correct for the pressure effects
      call pmolvol(t,s,p)

C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate,calcium, and fluoride
C
C Uppstrom (1974)
      bt = 0.000232 * scl/10.811
C Morris & Riley (1966)
      st = 0.14 * scl/96.062
c Millero (1982)
      calcium = .010280d0*s/35. 
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
      x1 = 10.0**(-phhi)
      x2 = 10.0**(-phlo)
      xacc = 1.e-10
      htotal = drtsafe(ta_iter_1,x1,x2,xacc)
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
      htotal2 = htotal*htotal
      co2star = (dic*htotal2) / (htotal2 + k1*htotal + k1*k2)
      ! use ff instead of k0 to include effect of water vapour in air
      co2starair = xco2*ff*atmpres
      dco2star = co2starair - co2star
      ph = -log10(htotal)
      co3 = dic*(k1*k2/(htotal**2 + k1*htotal + k1*k2))
      omegaca = calcium*co3/kcal
      omegaar = calcium*co3/karag
C
! use ff to be consistent with the calculation of co2starair
! I think ff is used to convert mol/kg to mol flux
! while k0 converts mol/m3 to mol flux
! I will use k0 since the fluxes are computed with pco2 to give mol/m3
      fco2 = co2star*1e6/k0 ! using k0 fixes constant loss of DIC
      !fco2 = co2star*1e6/ff
      ff_out = k0      !k0  !ff
      
C  Convert units of output arguments
c      Note: co2star and dco2star are calculated in mol/kg within this routine 
c      Thus Convert now from mol/kg -> m mol/m^3
      co2star  = co2star / permil
      dco2star = dco2star / permil
      co2starair=co2starair/permil
      co3 = co3 /permil

!	print *,'yyy', co2star,co2starair,dco2star
      return
      end

