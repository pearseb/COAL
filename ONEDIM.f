c Implementing the equation of state of McDougall et al (2003).
c SJP 2009/04/08
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/18
c
c Defines the COMMON block /ONEDIM/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real dxt, dxtr, dxt2r, dxu, dxur, dxu2r, dxu4r, dxt4r, sfu, sfub,
     &     sfv, sfvb, zun, zus, zvn, zvs, dyt, dytr, dyt2r, dyu, dyur,
     &     dyu2r, dyu4r, dyt4r, cs, csr, cst, cstr, phi, phit, sine,
     &     tng, c2dz, dz, dz2r, eem, ffm, zdz, ahi, dzz, dzz2r, zdzz,
     &     tinit, zdb
      COMMON /ONEDIM/
     & DXT  (IMT),DXTR (IMT),DXT2R(IMT),DXU  (IMT),DXUR (IMT),DXU2R(IMT)
     &,DXU4R(IMT),DXT4R(IMT),SFU  (IMT),SFUB (IMT),SFV  (IMT),SFVB (IMT)
     &,ZUN  (IMT),ZUS  (IMT),ZVN  (IMT),ZVS  (IMT)
     &,DYT  (JMT),DYTR (JMT),DYT2R(JMT),DYU  (JMT),DYUR (JMT),DYU2R(JMT)
     &,DYU4R(JMT),DYT4R(JMT),CS   (JMT),CSR  (JMT),CST  (JMT),CSTR (JMT)
     &,PHI  (JMT),PHIT (JMT),SINE (JMT),TNG  (JMT)
     &,C2DZ ( KM),DZ   ( KM),DZ2R ( KM),EEM  ( KM)
     &,FFM  ( KM),ZDZ  ( KM),AHI  ( KM)
     &,DZZ (KMP1),DZZ2R(KMP1),ZDZZ(KMP1),TINIT(KM,NT), zdb(km)
