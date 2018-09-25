# 1 "OPARAMS.F"
c Incorporating the code which calculates the meridional overturning
c streamfunctions into the core model source code.
c SJP 2009/04/17
c
c Removed several obsolete parameters.
c SJP 2008/03/08
c
c Reduce the value of NISLE from 11 to 10 for the new ocean model grid.
c SJP 2008/02/03
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Increse the value of NISLE from 3 to 11 for the new ocean model grid.
c SJP 2007/10/18
c
c Modify the default parameter values for the new ocean model grid, reducing
c both JFT1 and JFU1 from 5 to 2.
c SJP 2007/07/06
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/18
c
c Modified to use OCEAN_LOW preprocessor macro, rather than OCEAN_DBL.
c SJP 2007/06/04
c
c (1) Removed the following parameters, which are not referenced anywhere
c     within the model: ILONG, ILAT, NAT, NAQ, NAX, NPX, NIX, NATP and NRIV.
c (2) Removed the line "include 'PARAMS.f'", as the contents of this file are
c     no longer referenced.
c (3) Modified to enable the new ocean model grid.
c SJP 2007/05/31
c
c Change the following parameters to the values used by Dave Bi:
c
c  JFT0   51 -> 54
c  JFT1    7 ->  3
c  JFT2   52 -> 55
c  JFU0   50 -> 54
c  JFU1    7 ->  3
c  JFU2   52 -> 55
c
c These changes remove all Fourier filtering in the Southern Hemisphere, and
c restrict the filtering in the Northern Hemisphere to latitudes greater than
c or equal to 81.2 degN (tracers) and 82.8 degN (velocities).
c SJP 2004/03/30
c
c Deleted NLAND, as this is not used.
c SJP 2003/12/29
c
c Parameters used by the OGCM.
c SJP 2003/04/29

      integer imt, jmt, km, nt, lseg, nisle, lsegf, jfrst, jft0, jft1,
     &        jft2, jfu0, jfu1, jfu2, imu, imtp1, imtm1, imtm2, imum1,
     &        imum2, jmtp1, jmtm1, jmtm2, jscan, kmp1, kmp2, kmm1,
     &        jskpt, jskpu, njtbft, njtbfu, nieven, ntmin2, nbasin
# 66

      PARAMETER (IMT=130,JMT=114,KM=21,NT=11,LSEG=5,NISLE=10
     &,JFRST=2,JFT0=107,JFT1=2,JFT2=108,JFU0=107,JFU1=2,JFU2=108
     &,NBASIN=4

     &,LSEGF=LSEG
     &,IMU=IMT
     &,IMTP1=IMT+1,IMTM1=IMT-1,IMTM2=IMT-2,IMUM1=IMU-1,IMUM2=IMU-2
     &,JMTP1=JMT+1,JMTM1=JMT-1,JMTM2=JMT-2,JSCAN=JMTM2
     &,KMP1=KM+1,KMP2=KM+2,KMM1=KM-1
     &,JSKPT=JFT2-JFT1,JSKPU=JFU2-JFU1
     &,NJTBFT=(JFT1-JFRST+1)+(JMTM1-JFT2+1)
     &,NJTBFU=(JFU1-JFRST+1)+(JMTM1-JFU2+1)
     &,NIEVEN=2*((NISLE+1)/2)
     &,NTMIN2=NT+1/NT)
