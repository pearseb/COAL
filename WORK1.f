c $Log: WORK1.f,v $
c Revision 1.5  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.4  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.3  1993/12/17  15:31:38  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  92/12/09  14:42:43  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  91/03/14  09:49:55  ldr
c Initial revision
c 
      real un(ln2,nl)
      real vn(ln2,nl)
      real pn(ln2)
      real ten(ln2,nl)
      real fun(ln2,nl)
      real fvn(ln2,nl)
      real rmn(ln2,nl)
      real rmnt(ln2,nl)

      common/work1/un,vn,pn,ten,fun,fvn,rmn,rmnt
