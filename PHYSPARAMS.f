c $Log: PHYSPARAMS.f,v $
c Revision 1.4  1996/11/28 01:03:35  mrd
c Add type declaration for pi.
c
c Revision 1.3  1996/10/24  01:03:32  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.2  1996/04/18  06:37:00  ldr
c Added space to RCS header to avoid annoying warning on Seca.
c
c Revision 1.1  1996/03/21  03:21:09  ldr
c Initial revision
c
c Physical constants
c
c
      real stefbo, erad, eradsq, cp, rdry, epsil, rvap, hlf, hl
     &, hls, hlcp, grav, hlfcp, hlscp, sq2, cappa, tomg, pi

      parameter (stefbo=5.67e-8) !Stefan-Boltzmann constant
      parameter (erad=6.37122e6, eradsq=erad*erad) !Radius of earth
      parameter (cp=1004.64) ! Specific heat of dry air at const P
      parameter (rdry=287.04) ! Specific gas const for dry air
      parameter (epsil=0.622) ! Ratio molec wt of H2O vapour to dry air
      parameter (rvap=461.) ! Gas constant for water vapour
      parameter (hlf=3.35e5) ! Latent heat of fusion (at 0 deg C)
      parameter (hl=2.5e6) !Latent heat of vaporization (at 0 deg. C)
      parameter (hls=hl+hlf) !  "      "   " sublimation
      parameter (hlcp=hl/cp, hlfcp=hlf/cp, hlscp=hlcp+hlfcp)
      parameter (grav=9.80616) ! Acceleration of gravity

      parameter (sq2=1.414213562373092) ! Square root of 2
      parameter (cappa=rdry/cp) 
      parameter (tomg=2*7.2921233e-5) ! 2*omega
      parameter (pi=3.14159265) !Good ol' pi
