c $Log: VERTV.f,v $
c Revision 1.6  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.5  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1994/08/08  16:22:32  ldr
c Fix up RCS header at top of file.
c
      integer nvertv
      parameter (nvertv=lat*2*(7*nl+nlm))

      real pvbtb(lat,2,nl)
      real pvtb (lat,2,nl)
      real pvbub(lat,2,nl)
      real pvub (lat,2,nl)
      real pvbqb(lat,2,nl)
      real pvqb (lat,2,nl)
      real psv  (lat,2,nl)
      real pssd (lat,2,nlm)

      common /vertv/ pvbtb, pvtb, pvbub, pvub, pvbqb, pvqb, psv, pssd

      real pvbtbx(nvertv)
      equivalence (pvbtbx(1),pvbtb(1,1,1))
