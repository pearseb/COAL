c Added bulk formula (heat and freshwater fluxes) functionality
c Pearse J Buchanan 16/07/2018

c $Log: SSTSAL.f,v $
c Revision 1.2  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1993/10/06  11:25:34  ldr
c Initial revision
c
c  
      real sstm(imt,jmt,2)
      real salm(imt,jmt,2)
      real spechumm(imt,jmt,2)
      real airtemm(imt,jmt,2)
      real slpm(imt,jmt,2)
      real swdownm(imt,jmt,2)
      real lwdownm(imt,jmt,2)
      real fcltm(imt,jmt,2)
      real rainm(imt,jmt,2)
      real snowm(imt,jmt,2)
      real roffm(imt,jmt,2)
      real ficem(imt,jmt,2)
      real uwndm(imt,jmt,2)
      real vwndm(imt,jmt,2)
      real dustm(imt,jmt,2)
      real noxm(imt,jmt,2)

      common /sstsal/ sstm, salm, spechumm, airtemm, slpm, swdownm, 
     +                lwdownm, fcltm, rainm, snowm, roffm, ficem, uwndm,
     +                vwndm, dustm, noxm
