c Purpose
c -------
c Defines the COMMON block /OHIST/, which contains the ocean model history
c arrays.
c
c Contains
c --------
c OHIST_SMFZON		Zonal wind stress
c OHIST_SMFMER		Meridional wind stress
c OHIST_STFHT		Surface heat flux
c OHIST_STFSAL		Surface surface salinity
c OHIST_TEMP		Potential temperature
c OHIST_SAL		Salinity
c OHIST_RHO		Density
c OHIST_U		Zonal velocity
c OHIST_V		Meridional velocity
c OHIST_W		Vertical velocity
c OHIST_UEDD		Zonal eddy-induced velocity
c OHIST_VEDD		Meridional eddy-induced velocity
c OHIST_WEDD		Vertical eddy-induced velocity
c OHIST_RES		Barotropic streamfunction
c OHIST_CDEPTHM		Maximum depth of convection
c
c History
c -------
c 2009 May 5	Steven Phipps	Original version

      real ohist_smfzon, ohist_smfmer, ohist_stfht, ohist_stfsal,
     &     ohist_temp, ohist_sal, ohist_rho, ohist_u, ohist_v, ohist_w,
     &     ohist_uedd, ohist_vedd, ohist_wedd, ohist_res, ohist_cdepthm
      common /ohist/ ohist_smfzon(2:imtm1, 2:jmtm1),
     &               ohist_smfmer(2:imtm1, 2:jmtm1),
     &               ohist_stfht(2:imtm1, 2:jmtm1),
     &               ohist_stfsal(2:imtm1, 2:jmtm1),
     &               ohist_temp(2:imtm1, 2:jmtm1, km),
     &               ohist_sal(2:imtm1, 2:jmtm1, km),
     &               ohist_rho(2:imtm1, 2:jmtm1, km),
     &               ohist_u(2:imtm1, 2:jmtm1, km),
     &               ohist_v(2:imtm1, 2:jmtm1, km),
     &               ohist_w(2:imtm1, 2:jmtm1, km),
     &               ohist_uedd(2:imtm1, 2:jmtm1, km),
     &               ohist_vedd(2:imtm1, 2:jmtm1, km),
     &               ohist_wedd(2:imtm1, 2:jmtm1, km),
     &               ohist_res(2:imtm1, 2:jmtm1),
     &               ohist_cdepthm(2:imtm1, 2:jmtm1)
