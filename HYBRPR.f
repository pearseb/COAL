c $Log: HYBRPR.f,v $
c Revision 1.3  2000/11/14 03:11:37  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.2  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.1  1995/08/18  06:10:24  ldr
c Initial revision
c
      real prf,prh,dprf,pdpsk,muf,pdp00k
      common/hybrpr/prf(ln2,nl),prh(ln2,nlp),dprf(ln2,nl)
     &,pdpsk(ln2,nl),muf(ln2,nl),pdp00k(ln2,nl)
