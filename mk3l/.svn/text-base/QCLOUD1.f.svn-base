c $Log: QCLOUD1.f,v $
c Revision 1.11  2001/02/22 05:34:44  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.10  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.9  1998/12/10  00:55:54  ldr
c HBG changes to V5-1-21
c
c Revision 1.8  1997/02/21  00:33:29  ldr
c Make nc=1 the default.
c
c Revision 1.7  1995/08/08  02:11:54  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.6  1995/06/30  02:44:41  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.5  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1994/08/08  13:10:33  ldr
c Implement memory-saving trick: NX for tracers and NC for qcloud.
c
c Revision 1.3  94/06/20  09:50:36  ldr
c Add global array for cloud fraction cfb.
c 
c Revision 1.2  94/05/13  14:55:44  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.1  93/09/24  11:46:27  ldr
c Initial revision
c 
      integer nc, lonc, latc
      parameter (nc=1)          !Memory saving trick - set NC=1 to use qcloud.
      parameter (lonc=nc*(ln2-1)+1, latc=nc*(lat-1)+1) ! = ln2, lat if nc=1
      real   qlb(lonc,nl,latc)!Cloud liquid water at t
      real  qlbm(lonc,nl,latc)!Cloud liquid water at t-1
      real   qfb(lonc,nl,latc)!Cloud ice at t
      real  qfbm(lonc,nl,latc)!Cloud ice at t-1
      real   cfb(lonc,nl,latc)!Cloud fraction (not a prognostic variable)
      real cfrad(lonc,nl,latc)!Time averaged cloud cover passed to radiation
      real qlrad(lonc,nl,latc)!Time averaged cld liq watr passed to radiation
      real qfrad(lonc,nl,latc)!Time averaged cloud ice passed to radiation
      real qcrad(lonc,nl,latc)!Time averaged conv cloud watr passed to radatn

      common /qcloud1/ qlb, qlbm, qfb, qfbm, cfb,
     &                 cfrad, qlrad, qfrad, qcrad
