c Purpose
c -------
c Calculates the density of seawater, using the equation of state of:
c
c McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel (2003), Accurate
c and Computationally Efficient Algorithms for Potential Temperature and
c Density of Seawater, Journal of Atmospheric and Oceanic Technology, 20,
c 730-741.
c
c Note that, for incorporation into the CSIRO Mk3L OGCM, the values of the
c coefficients have been modified so that the density is returned in g cm-3,
c rather than kg m-3.
c
c Inputs
c ------
c S		Salinity (psu)
c T		Potential temperature (degC)
c P		Pressure (db)
c
c Outputs
c -------
c RHO		Density (g cm-3)
c
c History
c -------
c 2009 May 8	Steven Phipps	Original version

      subroutine m2003(s, t, p, rho)

      implicit none

C Global parameters
      include 'OPARAMS.f' 

C Argument list
      real s(imt), t(imt), p, rho(imt)

C Local work arrays and variables
      integer i
      real p2, p3, s15(imt), s2(imt), t2(imt), t3(imt), t4(imt)

C Local data, functions etc
      real a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11,
     &     b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12
      parameter (a0 = 9.99843699e-1)
      parameter (a1 = 7.35212840e-3)
      parameter (a2 = -5.45928211e-5)
      parameter (a3 = 3.98476704e-7)
      parameter (a4 = 2.96938239e-3)
      parameter (a5 = -7.23268813e-6)
      parameter (a6 = 2.12382341e-6)
      parameter (a7 = 1.04004591e-5)
      parameter (a8 = 1.03970529e-10)
      parameter (a9 = 5.18761880e-9)
      parameter (a10 = -3.24041825e-11)
      parameter (a11 = -1.23869360e-14)
      parameter (b0 = 1.0)
      parameter (b1 = 7.28606739e-3)
      parameter (b2 = -4.60835542e-5)
      parameter (b3 = 3.68390573e-7)
      parameter (b4 = 1.80809186e-10)
      parameter (b5 = 2.14691708e-3)
      parameter (b6 = -9.27062484e-6)
      parameter (b7 = -1.78343643e-10)
      parameter (b8 = 4.76534122e-6)
      parameter (b9 = 1.63410736e-9)
      parameter (b10 = 5.30848875e-6)
      parameter (b11 = -3.03175128e-16)
      parameter (b12 = -1.27934137e-17)

C Start code : ------------------------------------------------------------

c...  Calculate density
      p2 = p*p
      p3 = p*p2
      do i = 1, imt
        s15(i) = s(i)*sqrt(s(i))
        s2(i) = s(i)*s(i)
        t2(i) = t(i)*t(i)
        t3(i) = t(i)*t2(i)
        t4(i) = t(i)*t3(i)
        rho(i) = (a0 + a1*t(i) + a2*t2(i) + a3*t3(i) + a4*s(i) +
     &            a5*s(i)*t(i) + a6*s2(i) + a7*p + a8*p*t2(i) +
     &            a9*p*s(i) + a10*p2 + a11*p2*t2(i)) /
     &           (b0 + b1*t(i) + b2*t2(i) + b3*t3(i) + b4*t4(i) +
     &            b5*s(i) + b6*s(i)*t(i) + b7*s(i)*t3(i) + b8*s15(i) +
     &            b9*s15(i)*t2(i) + b10*p + b11*p2*t3(i) + b12*p3*t(i))
      end do

      return
      end subroutine
