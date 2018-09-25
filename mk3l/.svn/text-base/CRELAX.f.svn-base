c Defines common block CRELAX, which contains the variables relating to the
c relaxation of the coupled model SSTs and SSSs towards prescribed values.
c
c SJP 2006/01/05

c Common block CRELAX contains:
c
c CRELAX_FLAG	If .TRUE., relaxation is employed within the coupled model
c
c CRELAX_TAU	The relaxation timescale (days)
c
c CRELAX_GAMMA	The relaxation constant (1/s)
c
c CRELAX_QNF	The additional heat flux (W/m^2)
c
c CRELAX_SALF	The additional salinity tendency (1/s)
c
c CRELAX_STFHT	The monthly-mean additional heat flux (W/m^2)
c
c CRELAX_STFSAL	The monthly-mean additional salinity tendency (1/s)

      logical crelax_flag
      real crelax_tau, crelax_gamma, crelax_qnf, crelax_salf,
     &     crelax_stfht, crelax_stfsal
      common /crelax/ crelax_flag, crelax_tau, crelax_gamma,
     &                crelax_qnf(imt, jmt), crelax_salf(imt, jmt),
     &                crelax_stfht(imt, jmt), crelax_stfsal(imt, jmt)
