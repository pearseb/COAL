c $Log: ECFIELD.f,v $
c Revision 1.1  2000/06/20 02:22:22  rot032
c Initial revision
c

      common / ecfield / field(ln2e,numfl2,late)
      real zoxidant(ln2e,numfl2,late)
      equivalence(field, zoxidant)
