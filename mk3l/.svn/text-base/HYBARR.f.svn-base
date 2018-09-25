c $Log: HYBARR.f,v $
c Revision 1.2  1997/06/13 06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.1  1995/08/18  06:10:24  ldr
c Initial revision
c
c Hybrid prressure level array indicators : p = an.P* + bn.P00
      real anh    ! Half level an
      real bnh    ! Half level bn
      real anf    ! Full level an
      real bnf    ! Full level bn
      real dadnf  ! Full level d(an)/d(eta)
      real dbdnf  ! Full level d(bn)/d(eta)
      common/hybarr/anh(nlp),bnh(nlp),anf(nl),bnf(nl)
     &,dadnf(nl),dbdnf(nl)
