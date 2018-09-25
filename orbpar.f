c This subroutine derives the three orbital parameters required by
c the model. It was obtained from the GISS atmosphere-ocean model
c website at http://aom.giss.nasa.gov/solar.html.
c
c It has been modified to
c (1) remove some unused data
c (2) improve the use of arrays, so that the code does not reference
c     outside array bounds
c (3) return the longitude of the perihelion, rather than the
c     longitude of the sun at perihelion
c (4) return angles in degrees, rather than radians.
c
c SJP 2001/12/23

      SUBROUTINE ORBPAR (YEAR,ECCEN,OBLIQ,OMEGVP)
C****
C**** ORBPAR calculates the three orbital parameters as a function of
C**** YEAR.  The source of these calculations is: Andre L. Berger,
C**** 1978, "Long-Term Variations of Daily Insolation and Quaternary
C**** Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
C**** Berger, May 1978, "A Simple Algorithm to Compute Long Term
C**** Variations of Daily Insolation", published by Institut
C**** D'Astronomie de Geophysique, Universite Catholique de Louvain,
C**** Louvain-la Neuve, No. 18.
C****
C**** Tables and equations refer to the first reference (JAS).  The
C**** corresponding table or equation in the second reference is
C**** enclosed in parentheses.
C****
C**** Input:  YEAR   = years A.D. are positive, B.C. are negative
C**** Output: ECCEN  = eccentricity of orbital ellipse
C****         OBLIQ  = latitude of Tropic of Cancer in radians
C****         OMEGVP = longitude of perihelion =
C****                = spatial angle from vernal equinox to perihelion
C****                  in radians with sun as angle vertex
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TWOPI=6.283185307179586477d0, PI180=TWOPI/360.d0)
      REAL*8 TABL10(3,47),TABL40(3,19),TABL50(3,10)
C**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
      DATA TABL10/-2462.22D0,  31.609970D0,  251.9025D0,
     2             -857.32D0,  32.620499D0,  280.8325D0,
     3             -629.32D0,  24.172195D0,  128.3057D0,
     4             -414.28D0,  31.983780D0,  292.7251D0,
     5             -311.76D0,  44.828339D0,   15.3747D0,
     6              308.94D0,  30.973251D0,  263.7952D0,
     7             -162.55D0,  43.668243D0,  308.4258D0,
     8             -116.11D0,  32.246689D0,  240.0099D0,
     9              101.12D0,  30.599442D0,  222.9725D0,
     O              -67.69D0,  42.681320D0,  268.7810D0,
     1               24.91D0,  43.836456D0,  316.7998D0,
     2               22.58D0,  47.439438D0,  319.6023D0,
     3              -21.16D0,  63.219955D0,  143.8050D0,
     4              -15.65D0,  64.230484D0,  172.7351D0,
     5               15.39D0,   1.010530D0,   28.9300D0,
     6               14.67D0,   7.437771D0,  123.5968D0,
     7              -11.73D0,  55.782181D0,   20.2082D0,
     8               10.27D0,   0.373813D0,   40.8226D0,
     9                6.49D0,  13.218362D0,  123.4722D0,
     O                5.85D0,  62.583237D0,  155.6977D0,
     1               -5.49D0,  63.593765D0,  184.6277D0,
     2               -5.43D0,  76.438309D0,  267.2771D0,
     3                5.16D0,  45.815262D0,   55.0196D0,
     4                5.08D0,   8.448301D0,  152.5268D0,
     5               -4.07D0,  56.792709D0,   49.1382D0,
     6                3.72D0,  49.747849D0,  204.6609D0,
     7                3.40D0,  12.058272D0,   56.5233D0,
     8               -2.83D0,  75.278214D0,  200.3284D0,
     9               -2.66D0,  65.241013D0,  201.6651D0,
     O               -2.57D0,  64.604294D0,  213.5577D0,
     1               -2.47D0,   1.647247D0,   17.0374D0,
     2                2.46D0,   7.811584D0,  164.4194D0,
     3                2.25D0,  12.207832D0,   94.5422D0,
     4               -2.08D0,  63.856659D0,  131.9124D0,
     5               -1.97D0,  56.155991D0,   61.0309D0,
     6               -1.88D0,  77.448837D0,  296.2073D0,
     7               -1.85D0,   6.801054D0,  135.4894D0,
     8                1.82D0,  62.209412D0,  114.8750D0,
     9                1.76D0,  20.656128D0,  247.0691D0,
     O               -1.54D0,  48.344406D0,  256.6113D0,
     1                1.47D0,  55.145462D0,   32.1008D0,
     2               -1.46D0,  69.000534D0,  143.6804D0,
     3                1.42D0,  11.071350D0,   16.8784D0,
     4               -1.18D0,  74.291306D0,  160.6835D0,
     5                1.18D0,  11.047742D0,   27.5932D0,
     6               -1.13D0,   0.636717D0,  348.1074D0,
     7                1.09D0,  12.844549D0,   82.6496D0/
C**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
      DATA TABL40/ .01860798D0,   4.207205D0,   28.620089D0,
     2             .01627522D0,   7.346091D0,  193.788772D0,
     3            -.01300660D0,  17.857263D0,  308.307024D0,
     4             .00988829D0,  17.220546D0,  320.199637D0,
     5            -.00336700D0,  16.846733D0,  279.376984D0,
     6             .00333077D0,   5.199079D0,   87.195000D0,
     7            -.00235400D0,  18.231076D0,  349.129677D0,
     8             .00140015D0,  26.216758D0,  128.443387D0,
     9             .00100700D0,   6.359169D0,  154.143880D0,
     O             .00085700D0,  16.210016D0,  291.269597D0,
     1             .00064990D0,   3.065181D0,  114.860583D0,
     2             .00059900D0,  16.583829D0,  332.092251D0,
     3             .00037800D0,  18.493980D0,  296.414411D0,
     4            -.00033700D0,   6.190953D0,  145.769910D0,
     5             .00027600D0,  18.867793D0,  337.237063D0,
     6             .00018200D0,  17.425567D0,  152.092288D0,
     7            -.00017400D0,   6.186001D0,  126.839891D0,
     8            -.00012400D0,  18.417441D0,  210.667199D0,
     9             .00001250D0,   0.667863D0,   72.108838D0/
C**** Table 5 (3).  General precession in longitude: psi
      DATA TABL50/ 7391.02D0,  31.609970D0,  251.9025D0,
     2             2555.15D0,  32.620499D0,  280.8325D0,
     3             2022.76D0,  24.172195D0,  128.3057D0,
     4            -1973.65D0,    .636717D0,  348.1074D0,
     5             1240.23D0,  31.983780D0,  292.7251D0,
     6              953.87D0,   3.138886D0,  165.1686D0,
     7             -931.75D0,  30.973251D0,  263.7952D0,
     8              872.38D0,  44.828339D0,   15.3747D0,
     9              606.35D0,    .991874D0,   58.5749D0,
     O             -496.03D0,    .373813D0,   40.8226D0/
C****
      YM1950 = YEAR-1950.
C****
C**** Obliquity from Table 1 (2):
C****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
C****   OBLIQD = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
C****
      SUMC = 0.
      DO 110 I=1,47
      ARG    = PI180*(YM1950*TABL10(2,I)/3600.+TABL10(3,I))
  110 SUMC   = SUMC + TABL10(1,I)*COS(ARG)
      OBLIQ  = 23.320556D0 + SUMC/3600.
C****
C**** Eccentricity from Table 4 (1):
C****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
C****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
C****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
C****
      ESINPI = 0.
      ECOSPI = 0.
      DO 210 I=1,19
      ARG    = PI180*(YM1950*TABL40(2,I)/3600.+TABL40(3,I))
      ESINPI = ESINPI + TABL40(1,I)*SIN(ARG)
  210 ECOSPI = ECOSPI + TABL40(1,I)*COS(ARG)
      ECCEN  = SQRT(ESINPI*ESINPI+ECOSPI*ECOSPI)
C****
C**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
C****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
C****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
C****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
C****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
C****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
C****
      PIE = ATAN2(ESINPI,ECOSPI)
      FSINFD = 0.
      DO 310 I=1,10
      ARG    = PI180*(YM1950*TABL50(2,I)/3600.+TABL50(3,I))
  310 FSINFD = FSINFD + TABL50(1,I)*SIN(ARG)
      PSI    = PI180*(3.392506D0+(YM1950*50.439273D0+FSINFD)/3600.)
      OMEGVP = MOD(PIE+PSI,TWOPI)
      IF(OMEGVP.lt.0.)  OMEGVP = OMEGVP + TWOPI
      OMEGVP = OMEGVP / PI180
C****
      RETURN
      END
