c Enabling user control over the ocean model statistics that are saved to file.
c SJP 2009/04/19
c
c Incorporating the code which calculates the meridional overturning
c streamfunctions into the core model source code.
c SJP 2009/04/17
c
c Implementing the equation of state of McDougall et al (2003).
c SJP 2009/04/08
c
c Implementing a Robert time filter.
c SJP 2008/12/17
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/18
c
c Defines the COMMON block /WORKSP/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      logical omask, vmask, veddmask
      integer basin
      real ta, ua, va, tbp, ubp, vbp, tp, up, vp, tb, ub, vb, wsy,
     &     t, u, v, wsx, tbm, ubm, vbm, tm, um, vm, uclin, vclin,
     &     usav, vsav, rhon, rhos, fuw, fvn, fvsu, fvst, fmm, fm,
     &     fmp, gm, uover, udif, uunder, vover, vdif, vunder, w,
     &     tempa, tempb, tdif, zuneng, zuseng, zvneng, zvseng, fk1,
     &     fk2, fk3, rxp, rx, ry, rym, rzp, rz, e, esav, dxtq, dxuq,
     &     dxt4rq, dxu2rq, dz2rq, dzzq, dzz2rq, c2dzq, eemq, ffmq,
     &     ahiq, ahifac, dtxq, cq, toq, soq, ciq, toiq, soiq,
     &     tf, uf, vf, ssfub, ssfvb, ssfubp, ssfvbp
      COMMON /WORKSP/
     &TA (IMT,KM,NT),UA (IMT,KM),VA (IMT,KM),
     &TBP(IMT,KM,NT),UBP(IMT,KM),VBP(IMT,KM),
     &TP (IMT,KM,NT),UP (IMT,KM),VP (IMT,KM),
     &TB (IMT,KM,NT),UB (IMT,KM),VB (IMT,KM),WSY (IMT,JMT),
     &T  (IMT,KM,NT),U  (IMT,KM),V  (IMT,KM),WSX (IMT,JMT),
     &TBM(IMT,KM,NT),UBM(IMT,KM),VBM(IMT,KM),
     &TM (IMT,KM,NT),UM (IMT,KM),VM (IMT,KM),
     & tf(imt, km, nt), uf(imt, km), vf(imt, km),
     & ssfub(imt), ssfvb(imt), ssfubp(imt), ssfvbp(imt)
      COMMON /WORKSP/
     & UCLIN(IMT,KM),VCLIN(IMT,KM),USAV (IMT,KM),VSAV (IMT,KM),
     & RHON (IMT,KM),RHOS (IMT,KM),FUW  (IMT,KM),FVN  (IMT,KM),
     & FVSU (IMT,KM),FVST (IMT,KM),
     & FMM  (IMT,KM),FM   (IMT,KM),FMP (IMT,KM),
     &               GM   (IMT,KM),
     &    UOVER(IMT),UDIF (IMT,KM),UUNDER(IMT),
     &    VOVER(IMT),VDIF (IMT,KM),VUNDER(IMT),
     & W(IMT,KMP1),TEMPA(IMT,KMP1),TEMPB(IMT,KMP1),
     & TDIF(IMT,KMP2,NTMIN2),
     & ZUNENG(IMT,8),ZUSENG(IMT,8),ZVNENG(IMT,8),ZVSENG(IMT,8)
     & ,FK1(IMT,KM,3:3),FK2(IMT,KM,3:3),FK3(IMT,KM,3)
     & ,RXP(IMT,KM),RX(IMT,KM),RY(IMT,KM),RYM(IMT,KM)
     & ,RZP(IMT,KMP1),RZ(IMT,KMP1),E(IMT,KMP1,3),ESAV(IMT,KM,NT)
      COMMON /WORKSP/
     & DXTQ  (IMT,KM),DXUQ  (IMT,KM),DXT4RQ(IMT,KM),DXU2RQ(IMT,KM),
     & DZ2RQ (IMT,KM),DZZQ  (IMT,KM),DZZ2RQ(IMT,KM),C2DZQ (IMT,KM),
     & EEMQ  (IMT,KM),FFMQ  (IMT,KM),AHIQ  (IMT,KM), AHIFAC (JMT,KM),
     & DTXQ   (IMT,KM),
     & CQ (IMT,KM,9  ),TOQ (IMT,KM  ),SOQ (IMT,KM  ),
     & CIQ(IMT,KM,9,2),TOIQ(IMT,KM,2),SOIQ(IMT,KM,2),
     & basin(2:imtm1, 2:jmtm1), omask(jmt, km),
     & vmask(2:jmtm1, km, nbasin), veddmask(2:jmtm1, km, nbasin)
