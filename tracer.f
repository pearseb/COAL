# 1 "tracer.F"
c Included capability for bulk formula of heat fluxes (shortwave,
c longwave, sensible and latent) for running ocean-only model
c Pearse J Buchanan 2018/07/17
c
c Calculation for total negative oxygen concentration due to
c remineralisation in absence of oxygen added
c   - oxylow(i,k)
c Pearse J Buchanan 2015/12/09
c
c (1) Removing the ITM=1 option from the ocean model source code.
c (2) Updating the default configuration of the high-resolution version of the
c     ocean model.
c SJP 2009/06/23
c
c Further modifications to the ocean model output routines to improve memory
c usage. The features that were causing stack overflows have now been resolved.
c SJP 2009/05/11
c
c Modifying the ocean model output routines for improved memory usage.
c SJP 2009/05/06
c
c Enabling user control over the ocean model statistics that are saved to file.
c SJP 2009/04/19
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Implementing the equation of state of McDougall et al (2003).
c SJP 2009/04/08
c
c Added code which checks for negative salinities, and which aborts the
c simulation if these are detected.
c SJP 2009/04/06
c
c Modified so that flux adjustments are only applied when FLUXADJ=T.
c SJP 2008/03/08
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c (1) Modify the mixing across unresolved straits, when using the new ocean
c     model grid, to reflect the reduction in the minimum depth of the ocean
c     from 160m to 50m.
c (2) Add an imposed mixing between the Caspian and Black Seas, when using the
c     new ocean model grid, at the rate of 0.01 Sv.
c SJP 2007/11/25
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Rationalise the diagnostic information written to standard output.
c SJP 2007/11/21
c
c Modified to remove the imposed mixing across Hudson Strait and the Strait of
c Gibraltar when using the new ocean model grid, as these straits are now
c resolved.
c SJP 2007/10/06
c
c Manually inlined the subroutine TRACER1, which comprised the ENTRY points
c TRACER4, TRACER5 and TRACER7, into this routine.
c SJP 2007/06/19
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/18
c
c Modified to implement mixing across unresolved straits for the new ocean
c model grid.
c SJP 2007/06/13
c
c Modified to use OCEAN_LOW preprocessor macro, rather than OCEAN_DBL.
c SJP 2007/06/04
c
c (1) Modified to enable the new ocean model grid.
c (2) Removed TRACRI, which is redundant, plus the call to DENCAL, the
c     references to the header files ICEOC.f, RHOJMIX.f and DEPL.f, and the
c     declarations of variables DIFF, SUM, KL, KN and JZX.
c SJP 2007/06/02
c
c Added the line "include 'PARAMS.f'", as this line is no longer included via
c the header file OPARAMS.f.
c SJP 2007/05/31
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Modified to enable relaxation of the coupled model SSTs and SSSs towards
c prescribed values.
c SJP 2006/01/05
c
c (1) Removed the section of code which resets ocean temperatures less than
c     -1.85 degC in the coupled model, but not in the stand-alone OGCM. As well
c     as creating a very significant energy conservation error, the OGCM should
c     have identical physics, whether it is running in stand-alone mode or as
c     part of the coupled model. Any differences in the physics are strongly
c     undesirable, as they will represent a source of drift in the coupled
c     model.
c (2) Added support for a variable relaxation timescale for the stand-alone
c     OGCM. The relaxation timescale, TRELAX, is obtained via /TRELAX/, and
c     GAMMA is calculated accordingly. If TRELAX=20.0, the same value for GAMMA
c     is used as previously, ensuring bit-for-bit reproduciblity with the
c     previous version of the model.
c SJP 2005/03/19
c
c Completely updated the code for mixing across unresolved straits. The
c previous code left both the Baltic Sea and Black Sea disconnected from the
c world ocean, and contained serious errors with regard to the mixing between
c the world ocean and both Hudson Bay and the Mediterranean Sea.
c
c The modified code fixes these errors, and mixes across five unresolved
c straits as follows:
c
c (1) Baltic Sea/North Sea                 0.02 Sv
c (2) Black Sea/Mediterranean Sea          0.01 Sv
c (3) Hudson Bay/Atlantic Ocean            0.1 Sv    [was 1.0 Sv]
c (4) Mediterranean Sea/Atlantic Ocean     0.75 Sv   [was 1.0 Sv]
c (5) Persian Gulf/Indian Ocean            0.2 Sv    [was 0.05 Sv]
c
c Only the Caspian Sea remains disconnected from the world ocean. However, as
c this is a landlocked sea in reality, this is appropriate.
c
c Note also that the modified code assumes a vertical resolution of 21 levels
c and that, if the vertical resolution were to be changed, then at the very
c least the values for k_baltic, k_black etc would have to be updated.
c SJP 2004/01/15
c
c Calls to TRACER1, TRACER2, TRACER3, TRACER8 and TRACER9 removed, as these are
c redundant.
c SJP 2004/01/06
c
c Array bounds violation in loop 900 fixed.
c SJP 2004/01/04
c
c TRACER1 split off into a separate file, tracer1.f.
c SJP 2003/12/30
c
c Adjacent loops fused, where possible.
c SJP 2003/12/19
c
c (1) Modified so that the Gent-McWilliams eddy parameterisation scheme is
c always used, enabling some optimisation. The parameter IGM is thus rendered
c irrelevant.
c (2) Some tidying up.
c SJP 2003/12/18
c
c Modified for changes to /orestart/.
c SJP 2003/09/04
c
c Modified to make use of /orestart/.
c SJP 2003/09/03
c
c Some calls to MATRIX commented out, as they were requesting that values
c out of array bounds be plotted. Some loops restructured to avoid array
c bounds violations.
c SJP 2003/05/02
c
c Parameter definitions moved to include file OPARAMS.f.
c SJP 2003/04/29
c
c $Log: tracer.f,v $
c Revision 1.21  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.20  1998/05/26  05:10:49  ldr
c Final 10k run changes (mainly ACH salinity stuff) merged into V5-1-9.
c
c Revision 1.19  1997/12/23  00:23:35  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.18.1.1  1998/05/26  04:48:56  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.18  1997/12/19  01:25:36  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Include mixing Persian Gulf - Indian Ocean at 0.05 Sv rate
c Include array for storing depth of convective mixed layer
c Change to 21 levels in ocean model, and insertion of
c eddy-induced transport (major changes delineated)
c
c Revision 1.17.1.1  1997/12/19  02:03:13  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.17  1997/03/06  23:33:59  ldr
c Corrections from HBG to the recent tidy-ups to ocean routines.
c
c Revision 1.16  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.15  1994/08/08  17:22:56  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.14  94/08/04  16:56:47  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c
c Revision 1.14.1.1  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.13  94/03/30  12:35:37  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.12  93/12/20  16:22:09  ldr
c Minor changes to V4-4-53l from HBG for coupled model
c
c Revision 1.11  93/12/17  15:34:16  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.10  93/12/17  11:51:53  ldr
c Changes to V4-4-45l from HBG for coupled model
c
c Revision 1.9  93/11/29  11:38:44  ldr
c Changes to V4-4-32l from HBG for coupled model
c
c Revision 1.8  93/11/03  11:44:39  ldr
c Changes to V4-4-24l from HBG for coupled model
c
c Revision 1.7  93/10/05  13:07:43  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.6  93/08/19  15:11:06  ldr
c Minor cosmetic changes.
c
      SUBROUTINE TRACER(J, EPmax)
C
C=======================================================================
C                                                                    ===
C  TRACER COMPUTES, FOR ONE ROW, THE NT TRACERS, WHERE:              ===
C           J=THE ROW NUMBER                                         ===
C                                                                    ===
C=======================================================================
C
      implicit none
C
C---------------------------------------------------------------------
C  DEFINE GLOBAL DATA
C---------------------------------------------------------------------
C
      include 'OPARAMS.f'
      include 'OCEAN_NML.f'
      include 'ORESTART.f'
      include 'FULLWD.f'
      include 'SCALAR.f'
      include 'ONEDIM.f'
      include 'WORKSP.f'
      include 'COEFFS.f'
      include 'ACHEXT.f'
      include 'SFC4.f'
      include 'ETRANS1.f'
      include 'ETRANS2.f'
      include 'ETRANS3.f'
      include 'LEVD.f'
      include 'CONVECTD.f'
      include 'FLXBLOK.f'
      include 'FEWFLAGS.f'
      include 'FCOR.f'
      include 'A2O.f'    !agcm
      include 'O2A.f'    !agcm
      include 'CRELAX.f'
      include 'OHIST.f'
      include 'bulkf_variables.f' !pjb

      integer i, j, k, l, m, kz, jj, isave, ieave, is, ie, iredo,
     &        im, n, mm, idx, ism1, iea, ieb, ii, kp1
      real fx, fxa, fxb, fxc, fxd, fxe, gdrho, a0, qsfc, boxvol
      REAL AA(IMT,KM),BB(IMT,KM),CC(IMT,KM),
     &     EE(IMT,KMP1),FF(IMT,KMP1,NT),EFDR(IMT),xn(imt,km)
      EQUIVALENCE (AA,E),(BB,FK1),(CC,FK3),(EE,TDIF)
      real atv(imt,km)
      real sden(imt), tden(imt), rden(imt)
      real rxue(imt,km), ryue(imt,km), rzue(imt,km), rxun(imt,km),
     &     ryun(imt,km), rzun(imt,km)

c...  These values are currently OK for both the old and new ocean model grids.
c...  They might need to be changed, however, if there was any change to the
c...  positions of the coastlines, or to the vertical resolution.
c...
c...  Ideally, the mixing across unresolved straits should be fully
c...  generalised, with all the associated parameters being read from an input
c...  file.
c...
c...  SJP 2007/06/

c...  Parameters and variables associated with mixing across unresolved
c...  straits. SJP 2004/01/15
c...
c...  Define the mean rates of exchange (in Sv)
# 301

      real f_baltic, f_black, f_caspian, f_persian, f_indo
      parameter (f_baltic  = 0.02)
      parameter (f_black   = 0.01)
      parameter (f_caspian = 0.01)
      parameter (f_persian = 0.2)
      parameter (f_indo    = 0.1)


c...  Define the sill depths (in m)
# 322

      real z_baltic, z_black, z_caspian, z_persian, z_indo
      parameter (z_baltic  = 160.0)
      parameter (z_black   = 160.0)
      parameter (z_caspian = 160.0)
      parameter (z_persian = 160.0)
      parameter (z_indo    = 80.0)


c...  Define the maximum levels at which mixing occurs
# 343

      integer k_baltic, k_black, k_caspian, k_persian, k_indo
      parameter (k_baltic  = 5)
      parameter (k_black   = 5)
      parameter (k_caspian = 5)
      parameter (k_persian = 5)
      parameter (k_indo    = 3)


c...  Declare variable to hold the exchange coefficient (in s^-1)
      real alpha

c rjm ...
      include "bio.h"
      include "extra.h"
c rjm

C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
C=======================================================================
C  BEGIN INTRODUCTORY SECTION, PREPARING VARIOUS  ======================
C  ARRAYS FOR THE COMPUTATION OF THE TRACERS      ======================
C=======================================================================
C
C-------------------------------------------------------------
C  COMPUTE FIVE COMPONENTS OF ISOPYCNAL DIFF. TENSOR FK
C  (SCALE DEN. GRADIENTS BY FXE TO KEEP EXPONENTS IN RANGE LATER)
C-------------------------------------------------------------
C
      FXA=1.E-25
      FXB=2.
      FXC=.5
      FXD=1.
      FXE=1.E10
      IF(J.EQ.2) THEN

        DO 202 K=1,KM
        DO 202 I=1,IMT
          RY(I,K)=0.
 202    CONTINUE

        DO 204 I=1,IMT
          RZP(I,1)=0.
 204    CONTINUE

        if (m2003_eos) then

          do k = 1, km, 2
            kp1 = k + 1
            if (kp1 .gt. km) kp1 = km
            if (omask(j, k)) then
              do i = 1, imt
                sden(i) = 1000.0*tb(i, k, 2)+35.0
                tden(i) = tb(i, k, 1)
              end do
              call m2003(sden, tden, zdb(kp1), rden)
              do i = 1, imt
                tempa(i, k) = rden(i)
              end do
              call m2003(sden, tden, zdb(k), rden)
              do i = 1, imt
                tempb(i, k) = rden(i)
              end do
            else
              do i = 1, imt
                tempa(i, k) = 1.0
                tempb(i, k) = 1.0
              end do
            end if
          end do

          do k = 2, km, 2
            kp1 = k + 1
            if (kp1 .gt. km) kp1 = km
            if (omask(j, k)) then
              do i = 1, imt
                sden(i) = 1000.0*tb(i, k, 2)+35.0
                tden(i) = tb(i, k, 1)
              end do
              call m2003(sden, tden, zdb(k), rden)
              do i = 1, imt
                tempa(i, k) = rden(i)
              end do
              call m2003(sden, tden, zdb(kp1), rden)
              do i = 1, imt
                tempb(i, k) = rden(i)
              end do
            else
              do i = 1, imt
                tempa(i, k) = 1.0
                tempb(i, k) = 1.0
              end do
            end if
          end do

        else

          CALL STATEC(TB,TB(1,1,2),TEMPA,TDIF,TDIF(1,1,2),1)
          CALL STATEC(TB,TB(1,1,2),TEMPB,TDIF,TDIF(1,1,2),2)

        end if

        DO K=1,KM
          DO I=1,IMTM1
            RXP(I,K)=FM(I,K)*FM(I+1,K)*(TEMPB(I+1,K)-TEMPB(I,K))*FXE
          END DO
          RXP(IMT,K)=RXP(2,K)
        END DO

        DO 208 K=2,KM,2
        DO 208 I=1,IMT
          RZP(I,K)=TEMPA(I,K-1)-TEMPA(I,K)
          RZP(I,K+1)=TEMPB(I,K)-TEMPB(I,K+1)
 208    CONTINUE

        DO 210 K=2,KM
        DO 210 I=1,IMT
          RZP(I,K)=FM(I,K-1)*FM(I,K)*RZP(I,K)*FXE
 210    CONTINUE

        DO 212 I=1,IMT
          RZP(I,KMP1)=0.
 212    CONTINUE

        DO 214 M=1,NT
        DO 214 K=1,KM
        DO 214 I=1,IMT
          ESAV(I,K,M)=0.
 214    CONTINUE

        DO 216 L=1,3
        DO 216 K=1,KMP1
        DO 216 I=1,IMT
          E(I,K,L)=0.
 216    CONTINUE

        FK3(IMT,KM,3)=0.
      ENDIF

      do k = 1, km
        do i = 1, imt
          RX(I,K)=RXP(I,K)
          RYM(I,K)=RY(I,K)
          RZ(I,K)=RZP(I,K)
        end do
      end do

      do i = 1, imt
        RZ(I,KMP1)=RZP(I,KMP1)
      end do

      if (m2003_eos) then

        do k = 1, km, 2
          kp1 = k + 1
          if (kp1 .gt. km) kp1 = km
          if (omask(j+1, k)) then
            do i = 1, imt
              sden(i) = 1000.0*tbp(i, k, 2)+35.0
              tden(i) = tbp(i, k, 1)
            end do
            call m2003(sden, tden, zdb(kp1), rden)
            do i = 1, imt
              tempa(i, k) = rden(i)
            end do
            call m2003(sden, tden, zdb(k), rden)
            do i = 1, imt
              tempb(i, k) = rden(i)
            end do
          else
            do i = 1, imt
              tempa(i, k) = 1.0
              tempb(i, k) = 1.0
            end do
          end if
        end do

        do k = 2, km, 2
          kp1 = k + 1
          if (kp1 .gt. km) kp1 = km
          if (omask(j+1, k)) then
            do i = 1, imt
              sden(i) = 1000.0*tbp(i, k, 2)+35.0
              tden(i) = tbp(i, k, 1)
            end do
            call m2003(sden, tden, zdb(k), rden)
            do i = 1, imt
              tempa(i, k) = rden(i)
            end do
            call m2003(sden, tden, zdb(kp1), rden)
            do i = 1, imt
              tempb(i, k) = rden(i)
            end do
          else
            do i = 1, imt
              tempa(i, k) = 1.0
              tempb(i, k) = 1.0
            end do
          end if
        end do

      else

        CALL STATEC(TBP,TBP(1,1,2),TEMPA,TDIF,TDIF(1,1,2),1)
        CALL STATEC(TBP,TBP(1,1,2),TEMPB,TDIF,TDIF(1,1,2),2)

      end if

      DO 222 K=2,KM,2
      DO 222 I=1,IMT
        RZP(I,K)=TEMPA(I,K-1)-TEMPA(I,K)
        RZP(I,K+1)=TEMPB(I,K)-TEMPB(I,K+1)
 222  CONTINUE

      DO 224 K=2,KM
      DO 224 I=1,IMT
        RZP(I,K)=FMP(I,K-1)*FMP(I,K)*RZP(I,K)*FXE
 224  CONTINUE

      DO 226 I=1,IMT
        RZP(I,KMP1)=0.
 226  CONTINUE

      if (m2003_eos) then

        do k = 1, km, 2
          if (omask(j, k)) then
            do i = 1, imt
              sden(i) = 1000.0*tb(i, k, 2)+35.0
              tden(i) = tb(i, k, 1)
            end do
            call m2003(sden, tden, zdb(k), rden)
            do i = 1, imt
              tempa(i, k) = rden(i)
            end do
          else
            do i = 1, imt
              tempa(i, k) = 1.0
            end do
          end if
        end do

        do k = 2, km, 2
          kp1 = k + 1
          if (kp1 .gt. km) kp1 = km
          if (omask(j, k)) then
            do i = 1, imt
              sden(i) = 1000.0*tb(i, k, 2)+35.0
              tden(i) = tb(i, k, 1)
            end do
            call m2003(sden, tden, zdb(kp1), rden)
            do i = 1, imt
              tempa(i, k) = rden(i)
            end do
          else
            do i = 1, imt
              tempa(i, k) = 1.0
            end do
          end if
        end do

      else

        CALL STATEC(TB,TB(1,1,2),TEMPA,TDIF,TDIF(1,1,2),2)

      end if

      DO K=1,KM
        DO I=1,IMTM1
          RXP(I,K)=FMP(I,K)*FMP(I+1,K)*(TEMPB(I+1,K)-TEMPB(I,K))*FXE
          RY (I,K)=FMP(I,K)*FM(I,K)*(TEMPB(I,K)-TEMPA(I,K))*FXE
        END DO
        RXP(IMT,K)=RXP(2,K)
        RY(IMT,K)=RY(2,K)
        DO I=1,IMTM1
          E(I,K,1)=RX(I,K)*FXB*DXU2RQ(I,K)*CSTR(J)
          E(I,K,2)=(RY(I+1,K)+RY(I,K)+RYM(I+1,K)+RYM(I,K))*DYT4R(J)
          E(I,K,3)=(RZ(I,K)+RZ(I+1,K)+RZ(I,K+1)+RZ(I+1,K+1))*FXC*
     &     DZ2RQ(I,K)
          TEMPA(I,K)=-SQRT(E(I,K,1)**2+E(I,K,2)**2)*SLMXRF
     &       *SQRT(DTXQ(I,K))
        END DO
        E(IMT,K,1)=E(2,K,1)
        E(IMT,K,2)=E(2,K,2)
        E(IMT,K,3)=E(2,K,3)
        TEMPA(IMT,K)=TEMPA(2,K)
      END DO

c             store density gradients for use in computing uedd

        do k=2,km
        do i=2,imtm1
          rxue(i,k) = (e(i,k,1) + e(i,k-1,1)) * 0.5
        enddo
        enddo

      do k = 1, km
        IF(E(1,K,3).GT.TEMPA(1,K)) E(1,K,3)=TEMPA(1,K)
        TEMPB(1,K)=-E(1,K,1)/(E(1,K,3)**2+FXA)
        FK1(1,K,3)=E(1,K,3)*TEMPB(1,K)
        do i = 2, imt
          IF(E(I,K,3).GT.TEMPA(I,K)) E(I,K,3)=TEMPA(I,K)
          TEMPB(I,K)=-E(I,K,1)/(E(I,K,3)**2+FXA)
          FK1(I,K,3)=E(I,K,3)*TEMPB(I,K)
          E(I,K,1)=(RXP(I,K)+RXP(I-1,K)+RX(I,K)+RX(I-1,K))*DXT4RQ(I,K)*
     &     CSR(J)
          E(I,K,2)=RY(I,K)*DYUR(J)
          E(I,K,3)=(RZP(I,K)+RZ(I,K)+RZP(I,K+1)+RZ(I,K+1))*FXC*
     &      DZ2RQ(I,K)
          TEMPA(I,K)=-SQRT(E(I,K,1)**2+E(I,K,2)**2)*SLMXRF
     &       *SQRT(DTXQ(I,K))
        END DO
        E(1,K,1) = E(IMTM1,K,1)
        E(1,K,2) = E(IMTM1,K,2)
        E(1,K,3) = E(IMTM1,K,3)
        TEMPA(1,K) = TEMPA(IMTM1,K)
      end do

c             store density gradients for use in computing uedd

        do k=2,km
        do i=2,imtm1
          rxun(i,k) = (e(i,k,1) + e(i,k-1,1)) * 0.5
          ryun(i,k) = (e(i,k,2) + e(i,k-1,2)) * 0.5
        enddo
        enddo

      do i = 1, imt
        IF(E(I,1,3).GT.TEMPA(I,1)) E(I,1,3)=TEMPA(I,1)
        TEMPB(I,1)=-E(I,1,2)/(E(I,1,3)**2+FXA)
        FK2(I,1,3)=E(I,1,3)*TEMPB(I,1)
      end do

      DO K=2,KM
        IF(E(1,K,3).GT.TEMPA(1,K)) E(1,K,3)=TEMPA(1,K)
        TEMPB(1,K)=-E(1,K,2)/(E(1,K,3)**2+FXA)
        FK2(1,K,3)=E(1,K,3)*TEMPB(1,K)
        DO I=2,IMT
          IF(E(I,K,3).GT.TEMPA(I,K)) E(I,K,3)=TEMPA(I,K)
          TEMPB(I,K)=-E(I,K,2)/(E(I,K,3)**2+FXA)
          FK2(I,K,3)=E(I,K,3)*TEMPB(I,K)
          E(I,K,1)=(RX(I,K-1)+RX(I-1,K-1)+RX(I,K)+RX(I-1,K))*
     &     DXT4RQ(I,K)*CSTR(J)
          E(I,K,2)=(RY(I,K-1)+RYM(I,K-1)+RY(I,K)+RYM(I,K))*DYT4R(J)
          E(I,K,3)=RZ(I,K)*DZZ2RQ(I,K)*FXB
          TEMPB(I,K)=E(I,K,1)**2+E(I,K,2)**2
          TEMPA(I,K)=-SQRT(TEMPB(I,K))*SLMXRF
     &       *SQRT(DTXQ(I,K))
          slope(i,k) = sqrt(tempb(i,k)/(e(i,k,3)**2+fxa)) * fm(i,k)
        END DO
        E(1,K,1)=E(IMTM1,K,1)
        E(1,K,2)=E(IMTM1,K,2)
        E(1,K,3)=E(IMTM1,K,3)
        TEMPB(1,K)=TEMPB(IMTM1,K)
        TEMPA(1,K)=TEMPA(IMTM1,K)
        slope(1,k) = slope(imtm1,k)
      END DO

      DO I=1,IMT
        E(I,1,1)=0.0
        E(I,1,2)=0.0
        E(I,1,3)=0.0
        TEMPB(I,1)=0.0
        TEMPA(I,1)=0.0
        slope(i,1) = 0.0
      END DO

      do k = 1, km
        if (k .gt. 1) then
          do i = 1, imtm1
            ryue(i,k) = (e(i,k,2) + e(i+1,k,2)) * 0.5
            rzue(i,k) = (e(i,k,3) + e(i+1,k,3)) * 0.5
            rzun(i,k) = (rz(i,k) + rzp(i,k)) * dzz2rq(i,k)
          end do
        end if
        do i = 1, imt
          IF(E(I,K,3).GT.TEMPA(I,K)) E(I,K,3)=TEMPA(I,K)
          TEMPA(I,K)=FXD/(E(I,K,3)**2+FXA)
          FK3(I,K,3)=TEMPA(I,K)*TEMPB(I,K)
          TEMPB(I,K)=-E(I,K,3)*TEMPA(I,K)
          FK3(I,K,1)=E(I,K,1)*TEMPB(I,K)
          FK3(I,K,2)=E(I,K,2)*TEMPB(I,K)
          FK1(I,K,3)=FK1(I,K,3)*AHIQ(I,K)
          FK2(I,K,3)=FK2(I,K,3)*AHIQ(I,K)
        end do
      end do

      DO 249 L=1,3
      DO 249 K=1,KM
      DO 249 I=1,IMT
        FK3(I,K,L)=FK3(I,K,L)*AHIQ(I,K)
 249  CONTINUE

c --------------------------------------------------------------------------
c         begin calculation of isopycnal-mixing transport velocity
c                according to gent and mcwilliams (1990, j.p.o.)

c                   uedd = d(ahi*rx/rz)/dz,  vedd = d(ahi*ry/rz)/dz,
c                      wedd =d(ahi*rx/rz)/dx + d(ahi*ry/rz)/dy

c  (i.e., the horizontal component is proportional to the vertical gradient
c     of the neutral surface slope, the vertical component is such as to
c     achieve continuity.  note that the condition of no normal velocity
c     at the surface and bottom/side boundaries requires that
c     ahe = 0 at all the above boundaries.
c impose slope limitation to keep transport velocities within reasonable limits
c -----------------------------------------------------------------------------

c       first,  impose slope constraint by adjusting vertical density gradient
c
c             second, compute diffusivity (ahs(k)) times slope of surface
c             note that ahs(k) is already defined at top of the respective
c             grid box.
c             (first invert rz, noting that rz here is always zero or negative)

      do k = 2, km
        do i = 2, imtm1
          tempa(i,k)=rxue(i,k)**2+ryue(i,k)**2
          tempa(i,k)=-sqrt(tempa(i,k))*slopemx(k)
          tempb(i,k)=rxun(i,k)**2+ryun(i,k)**2
          tempb(i,k)=-sqrt(tempb(i,k))*slopemx(k)
          if(rzue(i,k).gt.tempa(i,k)) rzue(i,k)=tempa(i,k)
          if(rzun(i,k).gt.tempb(i,k)) rzun(i,k)=tempb(i,k)
          tempa(i,k) = -1./(fxa-rzue(i,k))
          tempb(i,k)=  -1./(fxa-rzun(i,k))
          arxrz(i,k)=aheq(i,k)*rxue(i,k)*tempa(i,k)*fm(i,k)*fm(i+1,k)
          aryrz(i,k)=aheq(i,k)*ryun(i,k)*tempb(i,k)*fm(i,k)*fmp(i,k)
        end do
        arxrz(imt,k)=arxrz(2,k)
        aryrz(imt,k)=aryrz(2,k)
        arxrz(1,k)=arxrz(imtm1,k)
        aryrz(1,k)=aryrz(imtm1,k)
      end do

c   define bottom and surface values of ahi*rx/rz, to satisfy wedd=0 condition

      do i=1,imt                                                        
        arxrz(i,1)    = 0.                                              
        aryrz(i,1)    = 0.                                              
        arxrz(i,kmp1) = 0.                                              
        aryrz(i,kmp1) = 0.                                              
      enddo                                                             

c           third, take vertical gradient to give i.t. velocity

      do k=1,km                                                       
      do i=1,imt                                                     
        uedd(i,k) = fxb*dz2rq(i,k)*(arxrz(i,k)-arxrz(i,k+1))        
        vedd(i,k) = fxb*dz2rq(i,k)*(aryrz(i,k)-aryrz(i,k+1))       
      enddo                                                       
      enddo                                                      

c---------------------------------------------------------------------
c                   find wedd exactly as per the physical w
c---------------------------------------------------------------------

c  find advective coefficient 'fuwedd' for west  face of t box
c                           & 'fvnedd' for north face of t box
c  Following code modified 25.02.98 to ensure definition of all variables
c  Note: assumes j=1 row to be land (typically Antarctica)
c---------------------------------------------------------------------

      if(j.le.2)then

        do k=1,km
        do i=2,imt
          fvnedd(i,k) = 0.
        enddo
        enddo

      endif

      fxa=cstr(j)*dytr(j)                                            
      fxb=fxa*cs(j)                                                 

      do k = 1, km
        do i = 2, imt
          fuwedd(i,k)= uedd(i-1,k) * cstr(j) * 2.
          fvstedd(i,k)=fvnedd(i,k)*cst(j-1)*dyt(j-1)*fxa
          fvnedd(i,k)= vedd(i,k)*(dxuq(i,k)+dxuq(i-1,k))*fxb*dxt4rq(i,k)
        end do
        fuwedd(imt,k)=fuwedd(2,k)
        fuwedd(1,k)  =fuwedd(imtm1,k)
        fvnedd(imt,k)=fvnedd(2,k)
        fvnedd(1,k)  =fvnedd(imtm1,k)
        fvstedd(imt,k)=fvstedd(2,k)
        fvstedd(1,k)  =fvstedd(imtm1,k)
      end do

c---------------------------------------------------------------------
c  compute vertical velocity in t columns
c---------------------------------------------------------------------
                                                                      
c  1st, set vertical velocity at the surface to zero (rigid-lid)
                                                                    
      fx=0.0                                                       

      do i=1,imt                                              
        wedd(i,1)=fx                                             
      enddo                                                   

c
c  2nd, compute change of w between levels
c

      do k=1,km                                         
      do i=2,imtm1                                     
        wedd(i,k+1)=c2dzq(i,k)*((fuwedd(i+1,k)-fuwedd(i,k))*dxt4rq(i,k) 
     *                       +fvnedd(i  ,k)-fvstedd(i,k))               
      enddo                                                   
      enddo                                                  
                                                                       
c  3rd, integrate downward from the surface
                                                                     
      do k=1,km                                                 
      do i=2,imtm1                                             
        wedd(i,k+1)=wedd(i,k)+wedd(i,k+1)                         
      enddo                                                   
      enddo                                                  

c      eliminate roundoff error in bottom wedd (wedd(i,kmp1) < 1e-14 in test)

      do i=1,imt                                                        
        kz=kmt(i,j)                                                     
        wedd(i,kz+1)=0.                                                 
      enddo                                                  

      do k=2,kmp1                                                   
        wedd(imt,k)=wedd(2,k)                                          
        wedd(1,k)  =wedd(imtm1,k)                                     
      enddo                                                   

      if(nergy.eq.1)then                                           

       if (j .eq. jplot) then
        print 8055, j,itt,0.1                                       
 8055   format(/12x,'u edd    velocity for j =',i3,10x,'itt=',i7,  
     &                                              '  scl=',f6.2)
        call matrix(uedd,imt,1,imt,0,km,.1)                      
        print 8056, j,itt,0.1                                   
 8056   format(/12x,'v edd    velocity for j =',i3,10x,'itt=',i7,   
     &                                              '  scl=',f6.2) 
        call matrix(vedd,imt,1,imt,0,km,.1)                       
        print 8057, j,itt,1e-4                                   
 8057   format(/12x,'w edd    velocity for j =',i3,10x,'itt=',i7,       
     &                                              '  scl=',1pe10.2)   
        call matrix(wedd,imt,1,imt,0,kmp1,1.e-4)                        
       endif                                            

       do k=1,km                                 
       do i=2,imtm1
        if(sqrt(uedd(i,k)**2 + vedd(i,k)**2).gt.16.)print *,         
     *   'warn edd  u=',uedd(i,k),'  v=',vedd(i,k),'  at i,j,k:',i,j,k  
       enddo                                                  
       enddo                                                  

       do i=2,imtm1
        if(abs(wedd(i,kmp1)).gt.1e-15)print *,                      
     *   'warn -- wedd(kmp1)=',wedd(i,kmp1),'  at i,j,k:',i,j,k     
       enddo                                                  

      endif                                                       

c              end computation of eddy transport velocity
c -----------------------------------------------------------------------------

C  SET ATV TO BE FUNCTION OF STATIC STABILITY (GARGETT, 1984)

         GDRHO = 960.
         A0 = 1.E-3

         DO 962 K = 2, KM
         DO 962 I = 1, IMT
            IF(RZ(I,K).GE.0.) THEN
              ATV(I,K) = 1.E6
            ELSE
              XN(I,K)  = SQRT (-GDRHO*RZ(I,K)/(DZZ(K)*FXE))
              ATV(I,K) = A0/XN(I,K)
              IF(ATV(I,K).LT.0.3) ATV(I,K) = 0.3
            ENDIF
 962     CONTINUE

C  SET ATV TO BE LARGE IN SURFACE LAYERS, TO MIMIC WIND-FORCED MIXING

         DO 963 I = 1, IMT
              IF(ATV(I,2).LT.20.) ATV(I,2) = 20.
              IF(ATV(I,3).LT.1.5) ATV(I,3) = 1.5
 963     CONTINUE

         do k = 2, km
         do i = 1, imt
              FK3(I,K,3) = FK3(I,K,3) + ATV(I,K)
              Xdif(i,j,k) = FK3(i,k,1)
              Ydif(i,j,k) = FK3(i,k,2)
              Wdif(i,j,k) = FK3(i,k,3)
              Gar84(i,j,k) = ATV(i,k)
         enddo
         enddo

C
C---------------------------------------------------------------------
C  FIND ADVECTIVE COEFFICIENT 'FUW' FOR WEST  FACE OF T BOX
C                           & 'FVN' FOR NORTH FACE OF T BOX
C---------------------------------------------------------------------
C
      FXA=CSTR(J)*DYTR(J)
      FXB=FXA*CS(J)

      DO K=1,KM
        DO I=2,IMT
          FUW(I,K)=(U(I-1,K)*DYU (J  )+UM(I-1,K)*DYU (J-1  ))*FXA
          FVN(I,K)=(V(I  ,K)*DXUQ(I,K)+V (I-1,K)*DXUQ(I-1,K))*FXB
     &             *DXT4RQ(I,K)
        END DO
        FUW(1,K)=FUW(IMTM1,K)
        FVN(1,K)=FVN(IMTM1,K)
      END DO

C
C---------------------------------------------------------------------
C  COMPUTE VERTICAL VELOCITY IN T COLUMNS
C---------------------------------------------------------------------
C
C  1ST, SET VERTICAL VELOCITY AT THE SURFACE TO ZERO (RIGID-LID)
C

      FX=0.0

      DO 700 I=1,IMT
        W(I,1)=FX
 700  CONTINUE

C
C  2ND, COMPUTE CHANGE OF W BETWEEN LEVELS
C

      DO K=1,KM
        DO I=1,IMTM1
          W(I,K+1)=C2DZQ(I,K)*((FUW(I+1,K)-FUW (I,K))*DXT4RQ(I,K)
     &                         +FVN(I  ,K)-FVST(I,K))
        END DO
        W(IMT,K+1)=W(2,K+1)
      END DO

C
C  3RD, INTEGRATE DOWNWARD FROM THE SURFACE
C

      DO 712 K=1,KM
      DO 712 I=1,IMT
        W(I,K+1)=W(I,K)+W(I,K+1)
 712  CONTINUE

C
C---------------------------------------------------------------------
C  SET BOUNDARY CONDITIONS FOR VERTICAL DIFFUSION OF TRACERS
C---------------------------------------------------------------------
C
      FX=0.0

      DO 730 I=1,IMT
        KZ=KMT(I,J)
        W(I,KZ+1)=FX
 730  CONTINUE

C
C=======================================================================
C  END INTRODUCTORY SECTION  ===========================================
C=======================================================================
C
C=======================================================================
C  BEGIN COMPUTATION OF THE TRACERS.               =====================
C  THE NEW VALUES "TA", WILL FIRST BE LOADED WITH  =====================
C  THE TIME RATE OF CHANGE, AND THEN UPDATED.      =====================
C=======================================================================
C

      DO 855 M=1,NT
C
C---------------------------------------------------------------------
C  COMPUTE TOTAL ADVECTION OF TRACERS
C---------------------------------------------------------------------

c     Note lower case additions in the following 40 lines for eddy transport

C
C  1ST, COMPUTE FLUX THROUGH WEST FACE OF T BOX
C
C  2ND, COMPUTE ZONAL FLUX DIVERGENCE
C
C  3RD, ADD IN MERIDIONAL FLUX DIVERGENCE
C
C  4TH, COMPUTE FLUX THROUGH TOP OF T BOX
C

      DO K=1,KM
        DO I=2,IMT
          TEMPA(I,K)=FUW(I,K)*(T(I,K,M)+T(I-1,K,M))
     *                   + fuwedd(i,k)*(t(i,k,m)+t(i-1,k,m))
        ENDDO    
        TEMPA(1,K)=TEMPA(IMTM1,K)
        DO I=1,IMTM1
          TA(I,K,M)=(TEMPA(I,K)-TEMPA(I+1,K))*DXT4RQ(I,K)
          
          TA(I,K,M)=TA(I,K,M)-FVN (I,K)*(TP(I,K,M)+T (I,K,M))
     *                   - fvnedd(i,k)*(tp(i,k,m)+t (i,k,m))
     &                       +FVST(I,K)*(T (I,K,M)+TM(I,K,M))
     *                   + fvstedd(i,k)*(t(i,k,m)+tm(i,k,m))

        END DO
        TA(IMT,K,M)=TA(2,K,M)
        if (k .gt. 1) then
          do i = 1, imt
            TEMPB(I,K)=W(I,K)*(T(I,K-1,M)+T(I,K,M))
     *                     + wedd(i,k)*(t(i,k-1,m)+t(i,k,m))
          end do
        end if
      ENDDO    

      DO I=1,IMT
        TEMPB(I,1)=0.0
        TEMPB(I,KMP1)=0.0
      ENDDO

C
C  5TH, ADD IN VERTICAL FLUX DIVERGENCE
C

      DO 824 K=1,KM
      DO 824 I=1,IMT
        TA(I,K,M)=TA(I,K,M)+(TEMPB(I,K+1)-TEMPB(I,K))*DZ2RQ(I,K)
 824  CONTINUE

C ---------- GAMMA = 1.1574E-6  corresponds to 10 day relaxation time --
C ---------- GAMMA = 5.787E-7   corresponds to 20 day relaxation time --
C ---------- GAMMA = 2.315E-7   corresponds to 50 day relaxation time --

      if (abs(trelax-20.0) .lt. 0.001) then
        GAMMA = 5.787e-7
      else
        gamma = 1.0 / (86400.0 * trelax)
      end if

      IF(M.EQ.1)THEN  ! Temperature
c.... RHOCW = RHOW * Cw = [1025 Kg/M**3] * [Cw J/Kg/K] = 4.1E6 J/M**3/K
c.... CTOP is RHOCW(J/M**3/K) * DZ(cm) * 0.01 = J/M**2/K
        CTOP    = RHOCW*DZ(1)*0.01

        If(.not.lcouple)Then

         DO 955 I=2,IMTM1
           if (.not. bulk_force) then
              ! Restore SSTs and SSSs to prescribed climatology
              TA(I,1,1) = TA(I,1,1) + GAMMA*(SST(I,J)-TB(I,1,1))*FM(I,1)
              FLUX(I,J,1) = CTOP*GAMMA*(SST(I,J)-TB(I,1,1))*FM(I,1)
           else
              ! 1. Solve for radiative and turbulent heat fluxes
              call bulkf_formula(
     I          uwnd(i,j),vwnd(i,j),TB(i,1,1),spechum(i,j),airtem(i,j),
     I          swdown(i,j),fclt(i,j),fice(i,j),lwdown(i,j),
     I          rain(i,j),snow(i,j),roff(i,j),
     O          swflx(i,j),lwup(i,j),lwflx(i,j),senflx(i,j),latflx(i,j),
     O          hflx(i,j),fflx(i,j),WSX(I,J),WSY(I,J),
     I          i,j, ncar_bulk_method, mit_bulk_method)
!              hflx(i,j) = hflx(i,j) + rain(i,j)*2.5e6 ! --> latent heat of vaporisation
!              hflx(i,j) = hflx(i,j) - snow(i,j)*3.337e5 ! --> latent heat of fusion

              ! 3. Account for heating in sea ice affected surface points,
              !    such that temperatures warm when temperatures dip below
              !    freezing.
*             call bulkf_seaice(
*     I            hflx(i,j),ctop,dttsf,sbcbio(i,j,3),T(i,1,1),TB(i,1,1),
*     O            frzflx(i,j),fwiflx(i,j),
*     I            pnu,i,j)
*             !    Update the relevant flux fields
*             latflx(i,j) = latflx(i,j) + frzflx(i,j)
*             hflx(i,j) = hflx(i,j) + frzflx(i,j)
*             fflx(i,j) = fflx(i,j) + fwiflx(i,j)

              ! 4. Apply the heat fluxes to the temperature variable
              !    (freshwater later)
              TA(I,1,1) = TA(I,1,1) + hflx(i,j)/ctop * fm(i,1)
              FLUX(I,J,1) = hflx(i,j)*fm(i,1)
 
              ! 5. Apply restoring to SSTs where atmospheric
              ! temperatures are less than 1.0 degC, basically where sea
              ! ice exists. This approach is warranted given that we do
              ! not account for latent heating during sea ice formation
              ! and latent cooling during sea ice melting.
              if (airtem(i,j).lt.274.15) then
                 TA(I,1,1) = TA(I,1,1) + 
     +                       GAMMA*(SST(I,J)-TB(I,1,1))*FM(I,1)
                 FLUX(I,J,1) = FLUX(I,J,1) + 
     +                       GAMMA*(SST(I,J)-TB(I,1,1))*FM(I,1)
              endif
            endif

  955   CONTINUE

          Else ! --> if in coupled mode

C ***********   FORCE OCEAN WITH AGCM FLUXES (IN MKS UNITS)*************

        do i = 2, imtm1
          if (fluxadj) then
            qsfc = osurh(i-1, j-1) - wt1 * hfcor(i, j, lmon1)
     &                             - wt2 * hfcor(i, j, lmon2)
          else
            qsfc = osurh(i-1, j-1)
          end if
          ta(i, 1, 1) = ta(i, 1, 1) + (qsfc / ctop) * fm(i, 1)
          flux(i, j, 1) = qsfc * fm(i, 1)
        end do

        if (crelax_flag) then
          crelax_gamma = 1.0 / (86400.0 * crelax_tau)
          do i = 2, imtm1
            ta(i, 1, 1) = ta(i, 1, 1) + crelax_gamma * (sst(i, j) -
     &                                  tb(i, 1, 1)) * fm(i, 1)
            crelax_qnf(i, j) = ctop * crelax_gamma * (sst(i, j) -
     &                         tb(i, 1, 1)) * fm(i, 1)
            flux(i, j, 1) = flux(i, j, 1) + crelax_qnf(i, j)
          end do
        end if

        End If

c.....  Save the surface heat flux if required
        if (save_stfht .and. (mix .eq. 0)) then
          do i = 2, imtm1
            ohist_stfht(i, j) = ohist_stfht(i, j) + flux(i, j, 1)
          end do
        end if

c rjm ...
      ELSE if (m.eq.2) then ! Salinity

        If(.not.lcouple)Then

        DO 956 I = 2,IMTM1
         if (.not.bulk_force) then
            TA(I,1,2) = TA(I,1,2) + GAMMA*(SAL(I,J)-TB(I,1,2))*FM(I,1)
            FLUX(I,J,2)= GAMMA*(SAL(i,j)-TB(I,1,2))*FM(I,1)
         endif
         if (bulk_force) then
            ! convert fw flux to salinity tendency (psu/s) and then
            ! back to model units (psu/1000 - 0.035)
            ! Apply relaxation term x2 on top of freshwater flux,
            ! accelerating the relaxation when seaice is present
            TA(I,1,2) = TA(I,1,2) + 
     +                  (fflx(i,j)*(TB(I,1,2)+0.035)/(DZ(1)*0.01) + 
     +                  GAMMA*(SAL(I,J)-TB(I,1,2)))*FM(I,1)
            FLUX(I,J,2)= (fflx(i,j)*(TB(I,1,2)+0.035)/(DZ(1)*0.01) + 
     +                   GAMMA*(SAL(i,j)-TB(I,1,2)))*FM(I,1)
         endif
         !TA(I,1,2) = TA(I,1,2) +  GAMMA*(SAL(I,J)-TB(I,1,2))*FM(I,1)
         !FLUX(I,J,2)= GAMMA*(SAL(i,j)-TB(I,1,2))*FM(I,1)
  956   CONTINUE

        Else

c........perform salinity calculations.........
        do i = 2, imtm1
          if (fluxadj) then
            qsfc = osalf(i-1, j-1) - wt1 * sfcor(i, j, lmon1)
     &                             - wt2 * sfcor(i, j, lmon2)
          else
            qsfc = osalf(i-1, j-1)
          end if
          ta(i, 1, 2) = ta(i, 1, 2) + qsfc * fm(i, 1)
          flux(i, j, 2) = qsfc * fm(i, 1)
        end do

        if (crelax_flag) then
          crelax_gamma = 1.0 / (86400.0 * crelax_tau)
          do i = 2, imtm1
            ta(i, 1, 2) = ta(i, 1, 2) + crelax_gamma * (sal(i, j) -
     &                                  tb(i, 1, 2)) * fm(i, 1)
            crelax_salf(i, j) = crelax_gamma * (sal(i, j) -
     &                          tb(i, 1, 2)) * fm(i, 1)
            flux(i, j, 2) = flux(i, j, 2) + crelax_salf(i, j)
          end do
        end if

        End If

c.....  Save the surface salinity tendency if required
        if (save_stfsal .and. (mix .eq. 0)) then
          do i = 2, imtm1
            ohist_stfsal(i, j) = ohist_stfsal(i, j) + flux(i, j, 2)
          end do
        end if

      ENDIF

C
C  ADD IN ISOPYCNAL DIFFUSION OF TRACERS (EXCEPT ZZ TERM)
C-------------------------------------------------------------
C
      FXA=.5
      FXB=2.

      do k = 1, km
        do i = 1, imtm1
          ! First, apply diffusivities to size of grid in X and Y directions
          TEMPA(I,K)=AHIQ(I,K)+AHHQ(I,K) !isopycnal + horizontal diffusivity
          E(I,K,1)=
     &      TEMPA(I,K)*(TB(I+1,K,M)-TB(I,K,M))*FXB*DXU2RQ(I,K)*CSTR(J)
          E(I,K,2)=TEMPA(I,K)*(TBP(I,K,M)-TB(I,K,M))*DYUR(J)
        end do
        TEMPA(IMT,K)=AHIQ(IMT,K)+AHHQ(IMT,K)
        E(IMT,K,1)=E(2,K,1)
        E(IMT,K,2)=E(2,K,2)
        ! Calculate diffusive transport along isopycnals
        if ((k .gt. 1) .and. (k .lt. km)) then
          do i = 1, imt
            TEMPA(I,K)=FM(I,K-1)*(TB(I,K-1,M)-TB(I,K  ,M))
            TEMPB(I,K)=FM(I,K+1)*(TB(I,K  ,M)-TB(I,K+1,M))
          end do
          DO I=1,IMTM1
            E(I,K,1)=E(I,K,1)
     &        +FK1(I,K,3)*( TEMPA(I,K)+TEMPA(I+1,K)
     &                     +TEMPB(I,K)+TEMPB(I+1,K))*DZ2RQ(I,K)*FXA
            E(I,K,2)=E(I,K,2)
     &        +FK2(I,K,3)*( FM (I,K-1)*(TB (I,K-1,M)-TB (I,K  ,M))
     &                     +FM (I,K+1)*(TB (I,K  ,M)-TB (I,K+1,M))
     &                     +FMP(I,K-1)*(TBP(I,K-1,M)-TBP(I,K  ,M))
     &                     +FMP(I,K+1)*(TBP(I,K  ,M)-TBP(I,K+1,M)))
     &                     *DZ2RQ(I,K)*FXA
          END DO
          E(IMT,K,1)=E(2,K,1)
          E(IMT,K,2)=E(2,K,2)
        end if
      end do

      DO 260 I=1,IMTM1
        E(I, 1,1)=E(I, 1,1)
     &    +FK1(I,1,3)*( FM(I  ,2)*(TB (I  ,1,M)-TB (I  ,2,M))
     &                 +FM(I+1,2)*(TB (I+1,1,M)-TB (I+1,2,M)))
     &                 *DZ2RQ(I,1)*FXA
        E(I,KM,1)=E(I,KM,1)
     &    +FK1(I,KM,3)*( FM(I  ,KMM1)*(TB(I  ,KMM1,M)-TB(I  ,KM,M))
     &                  +FM(I+1,KMM1)*(TB(I+1,KMM1,M)-TB(I+1,KM,M)))
     &                  *DZ2RQ(I,KM)*FXA
        E(I, 1,2)=E(I, 1,2)
     &    +FK2(I,1,3)*( FM (I,2)*(TB (I,1,M)-TB (I,2,M))
     &                 +FMP(I,2)*(TBP(I,1,M)-TBP(I,2,M)))
     &                 *DZ2RQ(I,1)*FXA
        E(I,KM,2)=E(I,KM,2)
     &    +FK2(I,KM,3)*( FM (I,KMM1)*(TB (I,KMM1,M)-TB (I,KM,M))
     &                  +FMP(I,KMM1)*(TBP(I,KMM1,M)-TBP(I,KM,M)))
     &                  *DZ2RQ(I,KM)*FXA
        E(I,1,3)=0.
        E(I,KMP1,3)=0.
 260  CONTINUE

      E(IMT, 1,1)=E(2, 1,1)
      E(IMT,KM,1)=E(2,KM,1)
      E(IMT, 1,2)=E(2, 1,2)
      E(IMT,KM,2)=E(2,KM,2)
      E(IMT,1,3)=0.0
      E(IMT,KMP1,3)=0.0

      DO K=1,KM
        DO I=2,IMTM1
          TEMPA(I,K)=FM(I-1,K)*(TB(I  ,K,M)-TB(I-1,K,M))
          TEMPB(I,K)=FM(I+1,K)*(TB(I+1,K,M)-TB(I  ,K,M))
        END DO
        TEMPA(1,K)=TEMPA(IMTM1,K)
        TEMPB(1,K)=TEMPB(IMTM1,K)
        TEMPA(IMT,K)=TEMPA(2,K)
        TEMPB(IMT,K)=TEMPB(2,K)
      END DO

      DO 264 K=2,KM
      DO 264 I=1,IMT
        E(I,K,3)=
     &    FK3(I,K,1)*( TEMPA(I,K)+TEMPA(I,K-1)
     &                +TEMPB(I,K)+TEMPB(I,K-1))*DXT4RQ(I,K)*CSTR(J)
 264  CONTINUE

      DO 266 K=1,KM
      DO 266 I=1,IMT
        TEMPA(I,K)=FMM(I,K)*(TB (I,K,M)-TBM(I,K,M))
        TEMPB(I,K)=FMP(I,K)*(TBP(I,K,M)-TB (I,K,M))
 266  CONTINUE

      do k = 1, km
        if (k .gt. 1) then
          do i = 1, imt
            E(I,K,3)=E(I,K,3)
     &       +FK3(I,K,2)*( TEMPA(I,K)+TEMPA(I,K-1)
     &                    +TEMPB(I,K)+TEMPB(I,K-1))*DYT4R(J)
            E(I,K,3)=FM(I,K-1)*FM(I,K)*E(I,K,3)
          end do
        end if
        DO I=1,IMTM1
          E(I,K,1)=FM (I+1,K)*FM(I,K)*E(I,K,1)
          E(I,K,2)=FMP(I  ,K)*FM(I,K)*E(I,K,2)*CS(J)
        END DO
        E(IMT,K,1)=E(2,K,1)
        E(IMT,K,2)=E(2,K,2)
      end do

      DO K=1,KM
        DO I=2,IMT
           TA(I,K,M) = TA(I,K,M) + 
     &                 ((E(I,K,1)-E(I-1,K,1))*4.*DXT4RQ(I,K) +
     &                 (E(I,K,2)-ESAV(I,K,M))*DYTR(J))*CSTR(J) +
     &                 (E(I,K,3)-E(I,K+1,3))*DZ2RQ(I,K)*FXB
           ESAV(I,K,M)=E(I,K,2)
        END DO
        TA(1,K,M)=TA(IMTM1,K,M)
        ESAV(1,K,M)=ESAV(IMTM1,K,M)
      END DO

C
C---------------------------------------------------------------------
C  COMPUTE NEW TRACERS, RESETTING LAND POINTS TO ZERO
C---------------------------------------------------------------------
C
c rjm ...
      do k = 1,km
         do i = 1,imt
            source(i,k) = 0.
         enddo
      enddo

c Add call to bgc component
      if (m.gt. 2) then
         call bio_geo_chem (j, m, EPmax)
      endif
# 1380


      ! Begin calculation of tracers based on master equation
      do k = 1,km
         do i = 1,imt

* -------------------------------------------------------------------- *
*                                                                      *
*       The following "leap-frog" calculation takes the total physical *
*           and biogeochemical changes to each tracer for the current  *
*           timestep. The arrays and variables within the equation     *
*           can be defined as follows:                                 *
*               - TB(i,k,m) is the total concentration of the "m"th    *
*                 tracer at the previous timestep                      *
*               - C2DTTS = DTTSF*2.0                                   *
*                   where DTTSF is equal to the number of seconds per  *
*                   timestep for the ocean model                       *
*               - TA(i,k,m) is the total change in the "m"th tracer's  *
*                 concentration due to physical processes              *
*               - source(i,k,m) is the total change in the "m"th       *
*                 tracer's concentration due to biogeochemical         *
*                 processes calcualted by csiro_obgc.F                 *
*               - FM(i,k) is a land-ocean mask, where ocean grids      *
*                 are equal to 1.0, while land is equal to 0.0         *
*               - DTXQ(i,k) is an acceleration factor for each level   *
*                 of the model, and is set in the OCEAN_NML.F file     *
*                 to 1.0                                               *
*
         
            if (M.eq.n_age) then

               if (K.eq.1) then
                  TA(I,K,M) = 0.0
               else
                  TA(I,K,M) = (TB(I,K,M) + 
     +                        TA(I,K,M)*C2DTTS*DTXQ(I,K))*FM(I,K) +
     +                        (c2dtts/(86400.0*365.0))*fm(i,k)
               endif

            elseif (M.eq.n_tomz) then

               if (tb(i,k,n_oxy).gt.10.0) then
                  TA(I,K,M) = 0.0
               else
                  TA(I,K,M) = (TB(I,K,M) + 
     +                        TA(I,K,M)*C2DTTS*DTXQ(I,K))*FM(I,K) +
     +                        (c2dtts/(86400.0*365.0))*fm(i,k)
               endif
            
            else

               TA(I,K,M) = (TB(I,K,M) + 
     +                     TA(I,K,M)*C2DTTS*DTXQ(I,K))*FM(I,K) +
     +                     source(i,k)*c2dtts*fm(i,k)
               
            endif
*                                                                      *
* -------------------------------------------------------------------- *

* -------------------------------------------------------------------- *
               
*   Undertake various corrections to tracers caused by numerical       *
*   instabilities in the physical transport of tracers and prepare     *
*   important summary statistics for the Nitrogen Cycle for model      *
*   output and balancing the Nitrate budget                            *

               
            if (m.eq.n_oxy) then
               oxyneg(i,k,j) = min(ta(i,k,m),0.0)
               ta(i,k,m) = max(ta(i,k,m),0.)
               if (sedfluxes) sulS(i,k,j) = sulS(i,k,j)*86400.0
            endif
          
            if (m.eq.n_no3) then
               
               if (den) then   ! put the units in day-1
                  denP(i,k,j) = -denP(i,k,j)*86400.0
                  if (sedfluxes) denS(i,k,j) = -denS(i,k,j)*86400.0
               endif

               Nexp(i,k,j) = Nexp(i,k,j)*86400.0
              
               if (fix) then
                  Nfix(i,k,j) = Nfix(i,k,j)*86400.0
               endif
               
               ! Capture change in NO3 between timesteps
               no3change(i,k,j) = (ta(i,k,m)-tb(i,k,m))/2.

            endif
            
*            if (m.eq.n_n15 .and. k.le.kmt(i,j) ) then
*
*               if (abs(ta(i,k,n_no3)).lt.0.1) then
*
*                  ta(i,k,m) = ta(i,k,n_no3)*0.5
*
*               endif
*
*            endif
            
            if (m.eq.n_age) then
               ta(i,k,m) = max(0.0, ta(i,k,m))
            endif
            if (m.eq.n_tomz) then
               ta(i,k,m) = max(0.0, ta(i,k,m))
            endif
               
* -------------------------------------------------------------------- *
               
            call sediment(i,k,j,m)
c rjm
         enddo
      enddo

c rjm ...
      ! fix fe at bottom
*      do i = 1,imt
*         if (m.eq.n_fe) then
*            ta(i,kmt(i,j),m) = 0.1
*         endif
*      enddo
* rjm

C
 855  CONTINUE

c...  ***************************************************************
c...  ***  MIX ACROSS UNRESOLVED STRAITS AFTER CONVECTIVE MIXING  ***
c...  ***************************************************************
c...
c...  There are five unresolved straits across which we mix:
c...
c...  (1) Baltic Sea/North Sea
c...  (2) Black Sea/Mediterranean Sea
c...  (3) Hudson Bay/Atlantic Ocean          [OLD OCEAN MODEL GRID ONLY]
c...  (4) Mediterranean Sea/Atlantic Ocean   [OLD OCEAN MODEL GRID ONLY]
c...  (5) Persian Gulf/Indian Ocean
c...
c...  On the new ocean model grid, mixing is also performed between the Caspian
c...  Sea and the Black Sea. This is not physically realistic, as the Caspian
c...  Sea is landlocked in reality, but is intended to reduce drift within
c...  the coupled model. For reasons of backwards compatibility, this mixing is
c...  not currently performed on the old ocean model grid.
c...
c...  The sill depths, and the rates of exchange across each of these straits,
c...  are defined as parameters above.
c...
c...  If F represents the rate of exchange, dz the sill depth and dx(1), dx(2),
c...  dy(1), dy(2) the zonal and meridional dimensions of the gridboxes on
c...  either side of each strait, then the exchange coefficients are given by
c...
c...  alpha(1) = F / (dx(1) * dy(1) * dz)
c...  alpha(2) = F / (dx(2) * dy(2) * dz)
c...
c...  and the tracers are updated as follows:
c...
c...  t'(1) = t(1) + alpha(1) * c2dtts * (t(2) - t(1))
c...  t'(2) = t(2) + alpha(2) * c2dtts * (t(1) - t(2))
c...
c...  where t(1), t(2) are the previous tracer values and t'(1), t'(2) are the
c...  new values.
c...
c...  Note that, as a result of the leap-frog timestepping, we use c2dtts
c...  rather than dtts.

c...  (1) Mix between the Baltic Sea and the North Sea

# 1574


      if (j .eq. 93) then
        alpha = 1.0e10 * f_baltic / (cst(93) * dxt(9) * dyt(93) *
     &                               z_baltic)
        do m = 1, nt
          do k = 1, k_baltic
            ta(9, k, m) = ta(9, k, m) + alpha * c2dtts *
     &                                  (tbp(10, k, m) - tb(9, k, m))
          end do
        end do
      end if

      if (j .eq. 94) then
        alpha = 1.0e10 * f_baltic / (cst(94) * dxt(10) * dyt(94) *
     &                               z_baltic)
        do m = 1, nt
          do k = 1, k_baltic
            ta(10, k, m) = ta(10, k, m) + alpha * c2dtts *
     &                                    (tbm(9, k, m) - tb(10, k, m))
          end do
        end do
      end if



c...  (2) Mix between the Black Sea and the Mediterranean Sea

# 1619


c...  On the new ocean model grid, mix between points (9, 84:85) and
c...  (12, 84:85)

      if (j .eq. 84) then
        alpha = 1.0e10 * f_black / ((cst(84) * dxt(9) * dyt(84) +
     &                               cst(85) * dxt(9) * dyt(85)) *
     &                              z_black)
        do m = 1, nt
          do k = 1, k_black
            ta(9, k, m) = ta(9, k, m) + alpha * c2dtts *
     &                    (0.5 * (tb(12, k, m) + tbp(12, k, m)) - 
     &                     0.5 * (tb(9, k, m) + tbp(9, k, m)))
            ta(12, k, m) = ta(12, k, m) + alpha * c2dtts *
     &                     (0.5 * (tb(9, k, m) + tbp(9, k, m)) - 
     &                      0.5 * (tb(12, k, m) + tbp(12, k, m)))
          end do
        end do
      end if

      if (j .eq. 85) then
        alpha = 1.0e10 * f_black / ((cst(84) * dxt(9) * dyt(84) +
     &                               cst(85) * dxt(9) * dyt(85)) *
     &                              z_black)
        do m = 1, nt
          do k = 1, k_black
            ta(9, k, m) = ta(9, k, m) + alpha * c2dtts *
     &                    (0.5 * (tbm(12, k, m) + tb(12, k, m)) - 
     &                     0.5 * (tbm(9, k, m) + tb(9, k, m)))
            ta(12, k, m) = ta(12, k, m) + alpha * c2dtts *
     &                     (0.5 * (tbm(9, k, m) + tb(9, k, m)) -  
     &                      0.5 * (tbm(12, k, m) + tb(12, k, m)))
          end do
        end do
      end if



c...  (3) Mix between Hudson Bay and the Atlantic Ocean

# 1675


c...  This is not necessary on the new ocean model grid, as Hudson Strait is
c...  explicitly resolved



c...  (4) Mix between the Mediterranean Sea and the Atlantic Ocean

# 1698


c...  This is not necessary on the new ocean model grid, as the Strait of
c...  Gibraltar is explicitly resolved



c... (5) Mix between the Persian Gulf and the Indian Ocean - following the
c...     original approach by Tony Hirst, we mix between _both_ Persian Gulf
c...     gridpoints and the Indian Ocean on the old ocean model grid, in order
c...     to avoid "severe grid Pechlet violation giving negative salinity"

# 1776


c...  On the new ocean model grid, mix between points (23, 74) and (24, 73)
c...  only

      if (j .eq. 73) then
        alpha = 1.0e10 * f_persian / (cst(73) * dxt(24) * dyt(73) *
     &                                z_persian)
        do m = 1, nt
          do k = 1, k_persian
            ta(24, k, m) = ta(24, k, m) + alpha * c2dtts *
     &                                    (tbp(23, k, m) - tb(24, k, m))
          end do
        end do
      end if

      if (j .eq. 74) then
        alpha = 1.0e10 * f_persian / (cst(74) * dxt(23) * dyt(74) *
     &                                z_persian)
        do m = 1, nt
          do k = 1, k_persian
            ta(23, k, m) = ta(23, k, m) + alpha * c2dtts *
     &                                    (tbm(24, k, m) - tb(23, k, m))
          end do
        end do
      end if



c...  (6) Mix between the Caspian Sea and the Black Sea

# 1813


c...  On the new ocean model grid, mix between points (17, 84:85) and
c...  (20, 84:85)

      if (j .eq. 84) then
        alpha = 1.0e10 * f_caspian / ((cst(84) * dxt(17) * dyt(84) +
     &                                 cst(85) * dxt(17) * dyt(85)) *
     &                                z_caspian)
        do m = 1, nt
          do k = 1, k_caspian
            ta(17, k, m) = ta(17, k, m) + alpha * c2dtts *
     &                     (0.5 * (tb(20, k, m) + tbp(20, k, m)) -
     &                      0.5 * (tb(17, k, m) + tbp(17, k, m)))
            ta(20, k, m) = ta(20, k, m) + alpha * c2dtts *
     &                     (0.5 * (tb(17, k, m) + tbp(17, k, m)) -
     &                      0.5 * (tb(20, k, m) + tbp(20, k, m)))
          end do
        end do
      end if

      if (j .eq. 85) then
        alpha = 1.0e10 * f_caspian / ((cst(84) * dxt(17) * dyt(84) +
     &                                 cst(85) * dxt(17) * dyt(85)) *
     &                                z_caspian)
        do m = 1, nt
          do k = 1, k_caspian
            ta(17, k, m) = ta(17, k, m) + alpha * c2dtts *
     &                     (0.5 * (tbm(20, k, m) + tb(20, k, m)) -
     &                      0.5 * (tbm(17, k, m) + tb(17, k, m)))
            ta(20, k, m) = ta(20, k, m) + alpha * c2dtts *
     &                     (0.5 * (tbm(17, k, m) + tb(17, k, m)) -
     &                      0.5 * (tbm(20, k, m) + tb(20, k, m)))
          end do
        end do
      end if



c...  (7) Mix within the waters northeast of Sumatra, Indonesia

# 1894


c...  On the new ocean model grid, mix between points (40, 58:59) and
c...  (41, 58:59)

      if (j .eq. 58) then
        alpha = 1.0e10 * f_indo / ((cst(58) * dxt(40) * dyt(58) +
     &                              cst(59) * dxt(40) * dyt(59)) *
     &                              z_indo)
        do m = 1, nt
          do k = 1, k_indo
            ta(40, k, m) = ta(40, k, m) + alpha * c2dtts *
     &                     (0.5 * (tb(41, k, m) + tbp(41, k, m)) -
     &                      0.5 * (tb(40, k, m) + tbp(40, k, m)))
            ta(41, k, m) = ta(41, k, m) + alpha * c2dtts *
     &                     (0.5 * (tb(40, k, m) + tbp(40, k, m)) -
     &                      0.5 * (tb(41, k, m) + tbp(41, k, m)))
          end do
        end do
      end if

      if (j .eq. 59) then
        alpha = 1.0e10 * f_indo / ((cst(58) * dxt(40) * dyt(58) +
     &                              cst(59) * dxt(40) * dyt(59)) *
     &                              z_indo)
        do m = 1, nt
          do k = 1, k_indo
            ta(40, k, m) = ta(40, k, m) + alpha * c2dtts *
     &                     (0.5 * (tbm(41, k, m) + tb(41, k, m)) -
     &                      0.5 * (tbm(40, k, m) + tb(40, k, m)))
            ta(41, k, m) = ta(41, k, m) + alpha * c2dtts *
     &                     (0.5 * (tbm(40, k, m) + tb(40, k, m)) -
     &                      0.5 * (tbm(41, k, m) + tb(41, k, m)))
          end do
        end do
      end if



C  COMPUTE ZZ COMPONENT OF DIFFUSION IMPLICITLY
C-------------------------------------------------------------
C
      FXA=4.0
      FXB=1.0
      FXC=0.0

      DO 861 I=1,IMT
        AA(I,   1)=FXC
        CC(I,  KM)=FXC
        EE(I,KM+1)=FXC
 861  CONTINUE

      DO 862 M=1,NT
      DO 862 I=1,IMT
        FF(I,KM+1,M)=FXC
 862  CONTINUE

      DO 863 K=1,KM
      DO 863 I=1,IMT
        TEMPA(I,K)=C2DTTS*DZ2RQ(I,K)*FXA*FM(I,K)
     &    *DTXQ(I,K)
        TEMPB(I,K)=FK3(I,K,3)*DZZ2RQ(I,K)
 863  CONTINUE

      do i = 1, imt
        AA(I,KM)=TEMPA(I,KM)*TEMPB(I,KM)
        CC(I,1)=TEMPA(I,1)*TEMPB(I,2)*FM(I,2)
      end do

      do k = 2, kmm1
        do i = 1, imt
          AA(I,K)=TEMPA(I,K)*TEMPB(I,K)
          CC(I,K)=TEMPA(I,K)*TEMPB(I,K+1)*FM(I,K+1)
          BB(I,K)=AA(I,K)+CC(I,K)+FXB
        end do
      end do

      do i = 1, imt
        BB(I,1)=AA(I,1)+CC(I,1)+FXB
        BB(I,KM)=AA(I,KM)+CC(I,KM)+FXB
      end do

      DO 875 K=KM,1,-1
        DO 872 I=2,IMTM1
          EFDR(I)=FXB/(BB(I,K)-CC(I,K)*EE(I,K+1))
 872    CONTINUE
        DO 874 I=2,IMTM1
          EE(I,K)=AA(I,K)*EFDR(I)
 874    CONTINUE
      DO 875 M=1,NT
      DO 875 I=2,IMTM1
        FF(I,K,M)=(TA(I,K,M)+CC(I,K)*FF(I,K+1,M))*EFDR(I)
 875  CONTINUE

      DO 886 M=1,NT
        DO 880 I=2,IMTM1
          BB(I,1)=FF(I,1,M)
 880    CONTINUE
        DO 884 K=2,KM
        DO 884 I=2,IMTM1
          BB(I,K)=EE(I,K)*BB(I,K-1)+FF(I,K,M)
 884    CONTINUE
        DO 885 K=1,KM
        DO 885 I=1,IMT
          TA(I,K,M)=BB(I,K)*FM(I,K)
           if (N15diagnostics) then
           if (m.eq.n_n15 .and. j.eq.latj .and. i.eq.loni 
     +         .and. k.le.kmt(i,j)) then 
               print*, " "
               print*, " After Vertical Diffusion "
               print*, "NO3 (tb) = ", tb(i,k,n_no3)
               print*, "NO3 (t) = ", t(i,k,n_no3)
               print*, "NO3 (ta) = ", ta(i,k,n_no3)
               print*, "15N (tb) = ", tb(i,k,m)
               print*, "15N (t) = ", t(i,k,m)
               print*, "15N (ta) = ", ta(i,k,m)
               print*,"delta 15N = ",(((ta(i,k,m)/
     +          (ta(i,k,n_no3)-ta(i,k,m)))/1.0)-1.)*1000
            endif 
            endif
 885    CONTINUE
 886  CONTINUE

C
C---------------------------------------------------------------------
C  DO ANALYSIS OF TRACER FORCING ON ENERGY TIMESTEP
C---------------------------------------------------------------------
C
      IF(NERGY.EQ.0.OR.MXP.EQ.1) GO TO 920

      DO 910 I=2,IMTM1
        KZ=KMT(I,J)
        IF (KZ.EQ.0) GOTO 910
        DO 900 M=1,NT
        DO 900 K=1,KZ
          BOXVOL = CST(J)*DXT(I)*DYT(J)*DZ(K)
C
C  1ST, COMPUTE TRACER CHANGE DUE TO ADVECTION
C
          TTDTOT(2,M)=TTDTOT(2,M)+BOXVOL*
     &              ((-FUW (I+1,K)*(T (I+1,K,M)+T (I  ,K,M))
     &                +FUW (I  ,K)*(T (I  ,K,M)+T (I-1,K,M)))*DXT4R(I)
     &                -FVN (I  ,K)*(TP(I  ,K,M)+T (I  ,K,M))
     &                +FVST(I  ,K)*(T (I  ,K,M)+TM(I  ,K,M)))
          if (k .eq. 1) then
            TTDTOT(3,M)=TTDTOT(3,M)+BOXVOL*
     &                  (W(I,2)*(T(I,1,M)+T(I,2,M))
     &                  -W(I,1)*T(I,1,M))*DZ2R(1)
          else if (k .eq. km) then
            TTDTOT(3,M)=TTDTOT(3,M)+BOXVOL*
     &                  (W(I,KM+1)*T(I,KM,M)
     &                  -W(I,KM)*(T(I,KM-1,M)+T(I,KM,M)))*DZ2R(KM)
          else
            TTDTOT(3,M)=TTDTOT(3,M)+BOXVOL*
     &                  (W(I,K+1)*(T(I,K  ,M)+T(I,K+1,M))
     &                  -W(I,K  )*(T(I,K-1,M)+T(I,K  ,M)))*DZ2R(K)
          end if
 900    CONTINUE
C
C  4TH, COMPUTE TOTAL ENERGY EXCHANGE BETWEEN POTENTIAL AND KINETIC
C
        IF(KZ.LT.2) GO TO 910
        FX=CST(J)*DXT(I)*DYT(J)*GRAV*0.5
        DO 905 K=2,KZ
          BUOY=BUOY-FX*DZZ(K)*W(I,K)*(RHOS(I,K-1)+RHOS(I,K))
 905    CONTINUE
 910  CONTINUE

 920  CONTINUE
C
C=======================================================================
C  END COMPUTATION OF THE TRACERS  =====================================
C=======================================================================
C
C---------------------------------------------------------------------
C  INTEGRATE TOTAL CHANGES IN T,S AND SQUARED T,S ON ENERGY TIMESTEP
C---------------------------------------------------------------------
C
      IF(NERGY.EQ.1 .AND. MXP.NE.1) THEN

        DO 970 M=1,NT
        DO 970 K=1,KM
         FX=CST(J)*DYT(J)*DZ(K)/C2DTTS
        DO 970 I=2,IMTM1
          TTDTOT(1,M)=TTDTOT(1,M)+(TA(I,K,M)   -TB(I,K,M)   )*FX*DXT(I)
          TVAR(M)    =TVAR(M)    +(TA(I,K,M)**2-TB(I,K,M)**2)*FX*DXT(I)
 970    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  FOURIER FILTER TRACERS AT HIGH LATITUDES
C---------------------------------------------------------------------
C
      IF((J.GT.JFT1.AND.J.LT.JFT2).OR.J.LT.JFRST)GO TO 1190
      JJ=J-JFRST+1 ! jj varies from 107-112 as j varies from 108-113
      IF (J.GE.JFT2) JJ=JJ-JSKPT+1  ! turns jj into 1-5
C
C  IF PREVIOUS STRIPS WERE OF SAME LENGTH, DONT RECOMPUTE FOURIER COEFFS
C
      ISAVE=0
      IEAVE=0

      DO 1140 L=1,LSEGF
        DO 1135 K=1,KM
          IF(ISTF(JJ,L,K).EQ.0) GO TO 1135
          IS=ISTF(JJ,L,K)
          IE=IETF(JJ,L,K)
          IREDO=0
          IF(IS.NE.ISAVE .OR. IE.NE.IEAVE) THEN
            IREDO=-1
            ISAVE=IS
            IEAVE=IE
            IM=IE-IS+1
            IF(IM.NE.IMTM2.OR.KMT(1,J).LT.K) THEN
              M=1
              N=NINT(IM*CST(J)*CSTR(JFT0))
            ELSE
              M=3
              N=NINT(IM*CST(J)*CSTR(JFT0)*0.5)
            ENDIF
          ENDIF
          DO 1130 MM=1,NT
            IDX=IREDO+MM
            ISM1=IS-1
            IEA=IE
            IF(IE.GE.IMT) IEA=IMTM1
            DO 1100 I=IS,IEA
              TDIF(I-ISM1,K,1)=TA(I,K,MM)
 1100       CONTINUE
            IF(IE.GE.IMT) THEN
              IEB=IE-IMTM2
              II=IMTM1-IS
              DO 1102 I=2,IEB
                TDIF(I+II,K,1)=TA(I,K,MM)
 1102         CONTINUE
            ENDIF
C
            CALL FILTER(TDIF(1,K,1),IM,M,N,IDX)
C
            DO 1120 I=IS,IEA
              TA(I,K,MM)=TDIF(I-ISM1,K,1)
 1120       CONTINUE
            IF(IE.GE.IMT) THEN
              DO 1122 I=2,IEB
                TA(I,K,MM)=TDIF(I+II,K,1)
 1122         CONTINUE
            ENDIF
 1130     CONTINUE
 1135   CONTINUE
 1140 CONTINUE

 1190 CONTINUE

c...  Check for negative salinities, and abort if these have arisen. Negative
c...  salinities indicate that something has gone seriously wrong with the
c...  simulation, and it should not be allowed to continue. SJP 2009/04/06

      do k = 1, km
        do i = 2, imtm1
          if (ta(i, k, 2) .lt. -0.035) then
            write (*, *)
            write (*, *) "ABORTING: Negative salinity in ocean"
            write (*, *)
            write (*, *) "i = ", i
            write (*, *) "j = ", j
            write (*, *) "k = ", k
            write (*, *)
            write (*, *) "Salinity = ", 1000.0*ta(i, k, 2)+35.0, " psu"
            write (*, *)
            stop
          end if
        end do
      end do

      if (lcouple) then

      DO 9001 I=2,IMT-1
C*************FILL OSAL ARRAY FOR USE BY THE AGCM******************
      OSAL(I,J,2)=TA(I,1,2)+0.035
      OSAL(I,J,1)=T(I,1,2)+0.035
C*************FILL OSST ARRAY FOR USE BY THE AGCM******************
C*************CONVERT TO KELVIN AND MULTIPLY BY -1 ****************
      OSST(I,J,2)=-1.*(TA(I,1,1)+273.15)
 9001 OSST(I,J,1)=-1.*(T(I,1,1)+273.15)

      end if

C
C---------------------------------------------------------------------
C  ACCUMULATE INTEGRATED ABSOLUTE CHANGES IN T EVERY NTSI TIMESTEPS
C---------------------------------------------------------------------
C
      IF(MOD(ITT,NTSI).EQ.0) THEN
        FX=0.5*CST(J)*DYT(J)/C2DTTS

        do m = 1, nt
          do k = 1, km
            TDIF(1,K,M)=ABS(TA(1,K,M)-TB(1,K,M))*C2DZQ(1,K)*FX*
     &                                DXTQ(1,K)
            TDIF(1,K,M)=TDIF(1,K,M)/DTXQ(1,K)
            TDIF(IMT,K,M)=ABS(TA(IMT,K,M)-TB(IMT,K,M))*C2DZQ(IMT,K)*FX*
     &                                DXTQ(IMT,K)
            TDIF(IMT,K,M)=TDIF(IMT,K,M)/DTXQ(IMT,K)
            do i = 2, imtm1
              TDIF(I,K,M)=ABS(TA(I,K,M)-TB(I,K,M))*C2DZQ(I,K)*FX*
     &                                  DXTQ(I,K)
              TDIF(I,K,M)=TDIF(I,K,M)/DTXQ(I,K)
              DTABS(M)=DTABS(M)+TDIF(I,K,M)
            end do
          end do
        end do

      ENDIF
C
C---------------------------------------------------------------------
C  TRANSFER QUANTITIES COMPUTED TO THE NORTH OF THE PRESENT ROW
C  TO BE DEFINED TO THE SOUTH IN THE COMPUTATION OF THE NEXT ROW
C---------------------------------------------------------------------
C
      FX=CST(J)*DYT(J)*CSTR(J+1)*DYTR(J+1)

      DO 990 K=1,KM
      DO 990 I=1,IMT
        FVST(I,K)=FVN(I,K)*FX
        RHOS(I,K)=RHON(I,K)
 990  CONTINUE

C
C---------------------------------------------------------------------
C  SET CYCLIC BOUNDARY CONDITIONS ON NEWLY COMPUTED TRACERS
C---------------------------------------------------------------------
C

      DO 992 M=1,NT
      DO 992 K=1,KM
        TA(1  ,K,M)=TA(IMTM1,K,M)
        TA(IMT,K,M)=TA(2    ,K,M)
 992  CONTINUE

C
C---------------------------------------------------------------------
C  SET NEW VELOCITIES AT NORTHERN WALL TO ZERO SINCE NO PASS THROUGH
C  CLINIC IS MADE FOR THIS ROW
C---------------------------------------------------------------------
C
      IF(J.EQ.JMTM1) THEN
        FX=0.0

        DO 680 K=1,KM
        DO 680 I=1,IMT
          UA(I,K)=FX
          VA(I,K)=FX
 680    CONTINUE

      ENDIF

      RETURN
      END
