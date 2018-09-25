c $Log: UVPGD.f,v $
c Revision 1.4  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.3  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/12/17  15:31:36  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  93/03/15  15:49:30  ldr
c Initial revision
c 
      real ugd(lon,lat2,nl)  ! indexed for horiz SLT
      real vgd(lon,lat2,nl)
      real pgd(lon,lat2)
      real dpsdt(lon,lat2,nl)
      real sdot(lon,nl,lat2) ! indexed for vert  SLT

      common /uvpgd/ ugd, vgd, pgd, dpsdt, sdot

      real ureal(lon,lat2,nl),
     &     vreal(lon,lat2,nl),
     &     pstar(lon,lat2)
      equivalence (ugd,ureal),(vgd,vreal),(pgd,pstar)
