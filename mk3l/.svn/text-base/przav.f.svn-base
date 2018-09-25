c Removing the declaration of the COMMON block /GRADNS/ to a header file.
c SJP 2009/04/25
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/11
c
c $Log: przav.f,v $
c Revision 1.45  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.44  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.43  2001/02/12 05:39:42  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.42  2000/06/20 02:08:31  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.41  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.40  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.39  1998/12/10  00:55:28  ldr
c HBG changes to V5-1-21
c
c Revision 1.38  1997/05/19  07:31:53  ldr
c Implement zmean cfrac diagnostic for dcloud scheme.
c
c Revision 1.37  1996/10/24  01:03:08  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.36  1996/01/09  06:18:16  ldr
c This is the version used for the long run H06 and write-up.
c
c Revision 1.35  1995/11/30  04:40:16  ldr
c Final updates to qcloud scheme to bring it to long run h63 used for
c BMRC workshop and write-up.
c
c Revision 1.34  1995/08/18  06:14:46  ldr
c LDR's changes to V4-7-12l to bring cloud scheme to run g57.
c
c Revision 1.33  1995/08/14  05:29:54  ldr
c HBG's new improved printing routines.
c
c Revision 1.32  1995/08/08  02:02:17  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.31  1994/09/12  16:32:39  ldr
c Fix up cloud forcing diagnostic if clforflag=F.
c
c Revision 1.30  94/08/09  12:35:31  ldr
c Rationalize zonal mean diagnostics: Use conzp flag to control HBG's
c cloud map and introduce savezonu, savezonv for dump of zmean u,v fields.
c Also remove savezcls which is fairly irrelevant.
c 
c Revision 1.29  94/08/08  17:22:07  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.28  94/08/04  10:33:54  ldr
c Split zonal mean cloud water diagnostic into separate liquid and frozen
c components.
c 
c Revision 1.27  94/06/30  09:53:18  ldr
c Add more zonal mean diagnostic plot options.
c 
c Revision 1.26  94/06/13  11:01:20  ldr
c Fix up 18 level temperature header alignment.
c 
c Revision 1.25  94/05/26  14:08:25  ldr
c Add new diagnostic for zonal mean cloud cover over sea.
c 
c Revision 1.24  94/03/30  10:20:59  ldr
c Tidy up 18 level diagnostics and add zqlpath (liq water path) diagnostic.
c 
c Revision 1.23  93/12/17  15:33:29  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.22  93/11/23  17:03:05  ldr
c Fix for 18 level version: Average levels to squash 18 levs into 9, except
c for sdot, for which alternate levels only are printed.
c 
c Revision 1.21  93/11/05  10:18:58  ldr
c Correct format of cloud amout printouts to use nint rather than int.
c 
      subroutine przav(ndadd,itt)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ndadd
      integer itt

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'BGRND.f'
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'GAUSL.f'
      include 'LSMI.f'
      include 'RADAV.f'
      include 'PRINTT.f'
      include 'TIMEX.f'
      include 'VERTV.f'        !Coupled model
      include 'GRADNS.f'

      real u,v,t,ww,ps
      common/dynav/u(lat,2,nl),v(lat,2,nl),t(lat,2,nl)
     &,ww(lat,2,nlm),ps(lat,2)

      real gev,grn,ghf,grh,ght,gcl,gcll,gclm,gclh,grg
     & ,grt,gsg,gsin,gsou,galb
     & ,gta,gts,gtsl,gtss,gsbl,ggw,grns,grnc,grtclr,gsouclr
     & ,gclc,gasbal,gatbal,gsfbal
     & ,gochfx,gflxic,gflxml,gflxse,gflxmj,grebal
     & ,gsbl_l,gsbl_i,gsbl_s,gcolwv
     & ,guz,gvz,gtz,gww
      common/maty/gev,grn
     & ,ghf,grh(nl),ght(nl),gcl,gcll,gclm,gclh,grg
     & ,grt,gsg,gsin,gsou,galb
     & ,gta,gts,gtsl,gtss,gsbl,ggw,grns,grnc,grtclr,gsouclr
     & ,gclc,gasbal,gatbal,gsfbal
     & ,gochfx,gflxic,gflxml,gflxse,gflxmj,grebal
     & ,gsbl_l,gsbl_i,gsbl_s,gcolwv
     & ,guz(nl),gvz(nl),gtz(nl),gww(nl)
      integer ngevx
      parameter (ngevx=2*nl+37+4*nl)
      real gevx(ngevx)
      equivalence (gevx(1),gev)

      real zisum
      integer jice
      common/iceprt/zisum(lat,2),jice(lat,2)

C Local work arrays and variables
      real vt(lat2,nl),ut(lat,2,nlp)
      real umx(10),umxt(10)
      real specstat(35+nl)            
      real tslzx(lat,2),tsszx(lat,2),sblzx(lat,2),gwzx(lat,2)
      real qout(lat2,nl)

      integer lgu(10),lgut(10),kpnt(10)
      integer jetex(lat2,nl)
      integer knl(lat,2)

      character*5 uchnl(nl),vchnl(nl),wchnl(nl),tchnl(nl)
     & ,rhch(nl),htch(nl)

      character*50 header

      logical testmax,testmin,tester

      integer i
      integer ii
      integer inl
      integer j
      integer k
      integer kbm
      integer kbt
      integer ki
      integer kis
      integer kk
      integer kl
      integer kpass
      integer ks
      integer ksk
      integer ksmin
      integer ksum
      integer ktm
      integer ktt
      integer kx
      integer lath
      integer latmin
      integer lg
      integer lgns
      integer lnor
      integer lnx
      integer lsou
      integer lsx
      integer lx
      integer ma
      integer mg
      integer ns

      real avg
      real consdr
      real consrad
      real conser
      real daysec
      real gicsum
      real gicsumn
      real gicsums
      real gswclf
      real pmin
      real rainmx
      real swclf
      real umax
      real vmax
      real vmin
      real wicsum
      real wicsumn
      real wicsums
      real wland
      real wsea
      real xx

C Local data, functions etc

C Start code : ----------------------------------------------------------

c****
c**** TO PRINT ZONAL AVERAGES (DYNAMICS + PHYSICS)
c****

c**** DETERMINE THE NUMBER OF NON-LAND POINTS PER LATITUDE
      do 20 ns=1,2
         do 20 lg=1,lat
         kk=0
         do 10 mg=1,lon
            ma=mg+(ns-1)*lon
 10         if (imsl(ma,lg).ne.4) kk=kk+1
 20      knl(lg,ns)=kk
c**** DETERMINE THE SEA/LAND WEIGHTING FACTORS
      wsea=0.0
      do 30 lg=1,lat
 30      wsea=wsea+w(lg)*(knl(lg,1)+knl(lg,2))
      wsea=wsea*0.5
      wland=lon-wsea
c**** TO COMPUTE THE AVERAGES
      consdr=mstep/1440.0/ndadd   ! (1/NSTEPS)
c     Scaling for quantities accumulated once per radiation step
      consrad=consdr*nrad
      conser=1.0/ndadd
      avg=1.0/kdynm
c      avg=(float(nrunday)/float(ndadd))/float(kdynm)

c**** liquid water cloud stats

      if (qcloud) then
         do 60 ns=1,2
            do 60 lg=1,lat
            ks=max(1,knl(lg,ns))
            qlpz(lg,ns)=qlpz(lg,ns)*consdr
            qlpz(lg,ns)=qlpz(lg,ns)/ks !Cloud liq water path over sea
 60      continue

c                convert to g/kg for plotting
CSJP         do 80 j=1,lat2*nl
CSJP            qlz(j,1,1)=1000*qlz(j,1,1)*consdr !cloud liquid water
CSJP            qfz(j,1,1)=1000*qfz(j,1,1)*consdr !cloud ice
CSJP            cfracz(j,1,1)=cfracz(j,1,1)*consdr
CSJP            qsublz(j,1,1)=1.e9*qsublz(j,1,1)*consdr !sublimation of snow
CSJP                                                    !(10^-6 kg/kg/s) 
CSJP            qevapz(j,1,1)=1.e9*qevapz(j,1,1)*consdr !evaporation of rain
CSJP            qautoz(j,1,1)=1.e9*qautoz(j,1,1)*consdr !autoconversion of ql
CSJP            qcollz(j,1,1)=1.e9*qcollz(j,1,1)*consdr !collection of ql by r
CSJP 80         qaccrz(j,1,1)=1.e9*qaccrz(j,1,1)*consdr !accretion  of ql by s
        do 80 k = 1, nl
          do 80 j = 1, 2
            do 80 i = 1, lat
              qlz(i,j,k)=1000*qlz(i,j,k)*consdr !cloud liquid water
              qfz(i,j,k)=1000*qfz(i,j,k)*consdr !cloud ice
              cfracz(i,j,k)=cfracz(i,j,k)*consdr
              qsublz(i,j,k)=1.e9*qsublz(i,j,k)*consdr !sublimation of snow
                                                      !(10^-6 kg/kg/s)
              qevapz(i,j,k)=1.e9*qevapz(i,j,k)*consdr !evaporation of rain
              qautoz(i,j,k)=1.e9*qautoz(i,j,k)*consdr !autoconversion of ql
              qcollz(i,j,k)=1.e9*qcollz(i,j,k)*consdr !collection of ql by r
 80           qaccrz(i,j,k)=1.e9*qaccrz(i,j,k)*consdr !accretion  of ql by s

      end if

      if(.not.qcloud)then
CSJP        do j=1,lat2*nl
CSJP          cfracz(j,1,1)=cfracz(j,1,1)*consdr
CSJP        enddo
        do k = 1, nl
          do j = 1, 2
            do i = 1, lat
              cfracz(i,j,k)=cfracz(i,j,k)*consdr
            end do
          end do
        end do
      endif

c---------------------------------------------------------------------
c--------------------- Check number of levels ------------------------
c---------------------------------------------------------------------
c----  skip levels printed if nl>24
      ksk=1+(nl-1)/24

      do 85 j=1,ngevx
 85      gevx(j)=0.0

c---------------------------------------------------------------------
c--------------------- Dynamics Maps ---------------------------------
c---------------------------------------------------------------------

      do 90 k=1,nl
         do 90 ns=1,2
         do 90 lg=1,lat
         u(lg,ns,k)=u(lg,ns,k)*avg
         v(lg,ns,k)=v(lg,ns,k)*avg
         t(lg,ns,k)=t(lg,ns,k)*avg
         if (k .eq. nl) then
           ps(lg, ns) = ps(lg, ns) * avg
         else
           ww(lg, ns, k) = ww(lg, ns, k) * avg
         end if
 90   continue
CSJP c.... ww + ps has vertical extent nlm+1 = nl
CSJP 90      ww(lg,ns,k)=ww(lg,ns,k)*avg

      if (dynzp) then

         do 100 k=1,nl
            inl=nint(sig(k)*1000.0)
            write (uchnl(k),"(' U',i3.3)") inl
            write (vchnl(k),"(' V',i3.3)") inl
            write (tchnl(k),"(' T',i3.3)") inl
            inl=nint(sigh(k+1)*1000.0)
 100        write (wchnl(k),"(' W',i3.3)") inl

      do 102 ns=1,2
         do 102 lg=1,lat
         do 102 k=1,nl
         guz(k)=guz(k)+u(lg,ns,k)*w(lg)
         gvz(k)=gvz(k)+v(lg,ns,k)*w(lg)
         gtz(k)=gtz(k)+t(lg,ns,k)*w(lg)
         if (k .eq. nl) then
           gww(k) = gww(k) + ps(lg, ns) * w(lg)
         else
           gww(k) = gww(k) + ww(lg, ns, k) * w(lg)
         end if
 102  continue
CSJP 102     gww(k)=gww(k)+ww(lg,ns,k)*w(lg)
      do 104 k=1,nl
         guz(k)=guz(k)*0.5
         gvz(k)=gvz(k)*0.5
         gtz(k)=gtz(k)*0.5
 104     gww(k)=gww(k)*0.5

 105  format(3x,24a5)
 107  format(1x,i2,24f5.1)
c---- U map  -----
         print *,'  UZ-MAP  in m/sec'
         write (6,105) (uchnl(k),k=1,nl,ksk)
         do 110 ns=1,2
            do 110 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
 110        write (6,107) lx,(u(lx,ns,k),k=1,nl,ksk)
         write (6,105) (uchnl(k),k=1,nl,ksk)
         write (6,107) lx,(guz(k),k=1,nl,ksk)
c---- V map  -----
         print *,'  VZ-MAP  in m/sec'
         write (6,105) (vchnl(k),k=1,nl,ksk)
         do 120 ns=1,2
            do 120 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
 120        write (6,107) lx,(v(lx,ns,k),k=1,nl,ksk)
         write (6,105) (vchnl(k),k=1,nl,ksk)
         write (6,107) lx,(gvz(k),k=1,nl,ksk)
c---- T map  -----
         print *,'  TZ-MAP  in degrees C'
         write (6,105) (tchnl(k),k=1,nl,ksk)
         do 130 ns=1,2
            do 130 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
 130        write (6,107) lx,(t(lx,ns,k)-273.0,k=1,nl,ksk)
         write (6,105) (tchnl(k),k=1,nl,ksk)
         write (6,107) lx,(gtz(k)-273.0,k=1,nl,ksk)
c---- Sigma-dot map  -----
         print *,' SDZ-MAP  * 1.0e+05'
         write (6,105) (wchnl(k),k=1,nlm,ksk)
         do 140 ns=1,2
            do 140 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
 140        write (6,107) lx,(ww(lx,ns,k)*1.0e+05,k=1,nlm,ksk)
         write (6,105) (wchnl(k),k=1,nlm,ksk)
         write (6,107) lx,(gww(k)*1.0e+05,k=1,nlm,ksk)

      end if

c---- General energy balances, Evap, Rain, Surface heat flux etc.
c---- (Room for more stats here)

      do 150 ns=1,2
      do 150 lg=1,lat
         asbalz(lg,ns)=asbalz(lg,ns)*consdr
         atbalz(lg,ns)=atbalz(lg,ns)*consdr
         sfbalz(lg,ns)=sfbalz(lg,ns)*consdr
         sblz  (lg,ns)=sblz  (lg,ns)*consdr
         flxicz(lg,ns)=flxicz(lg,ns)*consdr
         flxmlz(lg,ns)=flxmlz(lg,ns)*consdr
         flxsez(lg,ns)=flxsez(lg,ns)*consdr
         flxmjz(lg,ns)=flxmjz(lg,ns)*consdr
         rebalz(lg,ns)=rebalz(lg,ns)*consdr
         sblz_l(lg,ns)=sblz_l(lg,ns)*consdr
         sblz_i(lg,ns)=sblz_i(lg,ns)*consdr
         sblz_s(lg,ns)=sblz_s(lg,ns)*consdr
         colwvz(lg,ns)=colwvz(lg,ns)*consdr
         evz(lg,ns)=evz(lg,ns)*conser
         rnz(lg,ns)=rnz(lg,ns)*conser
 150     hfz(lg,ns)=hfz(lg,ns)*consdr

      rainmx=-1.e30
      do 180 ns=1,2
         do 180 lg=1,lat
         rainmx=max(rainmx,rnz(lg,ns))
         gasbal=gasbal+asbalz(lg,ns)*w(lg)
         gatbal=gatbal+atbalz(lg,ns)*w(lg)
         gsfbal=gsfbal+sfbalz(lg,ns)*w(lg)
         gochfx=gochfx+sblz  (lg,ns)*w(lg)
         gflxic=gflxic+flxicz(lg,ns)*w(lg)
         gflxml=gflxml+flxmlz(lg,ns)*w(lg)
         gflxse=gflxse+flxsez(lg,ns)*w(lg)
         gflxmj=gflxmj+flxmjz(lg,ns)*w(lg)
         grebal=grebal+rebalz(lg,ns)*w(lg)
         gsbl_l=gsbl_l+sblz_l(lg,ns)*w(lg)
         gsbl_i=gsbl_i+sblz_i(lg,ns)*w(lg)
         gsbl_s=gsbl_s+sblz_s(lg,ns)*w(lg)
         gcolwv=gcolwv+colwvz(lg,ns)*w(lg)
         gev=gev+evz(lg,ns)*w(lg)
         grn=grn+rnz(lg,ns)*w(lg)
 180     ghf=ghf+hfz(lg,ns)*w(lg)
      gasbal=gasbal*0.5
      gatbal=gatbal*0.5
      gsfbal=gsfbal*0.5
      gochfx=gochfx*0.5/lon
      gflxic=gflxic*0.5
      gflxml=gflxml*0.5
      gflxse=gflxse*0.5
      gflxmj=gflxmj*0.5
      grebal=grebal*0.5
      gsbl_l=gsbl_l*0.5
      gsbl_i=gsbl_i*0.5
      gsbl_s=gsbl_s*0.5
      gcolwv=gcolwv*0.5
      gev=gev*0.5
      grn=grn*0.5
      ghf=ghf*0.5

      if (zavgp) then

c---- Columns of stats here
         print *,' General diagnostics'
         print *,'       psl     asbal  atbal  sfbal lwpath evap rain',
     &' hflux colwv  ohflx  flxic  flxml  flxse  flxmj  rebal',
     &'  sbl_l  sbl_i  sbl_s'
         do 200 ns=1,2
            do 200 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
            write (6,210) lx,ps(lx,ns),asbalz(lx,ns),atbalz(lx,ns),
     &       sfbalz(lx,ns),qlpz(lx,ns),evz(lx,ns),rnz(lx,ns),hfz(lx,ns)
     &      ,colwvz(lx,ns),sblz(lx,ns)/lon
     &      ,flxicz(lx,ns),flxmlz(lx,ns),flxsez(lx,ns),flxmjz(lx,ns)
     &      ,rebalz(lx,ns),sblz_l(lx,ns),sblz_i(lx,ns),sblz_s(lx,ns)
 200     continue
 210     format (1x,i2,2x,f9.3,3(1x,f6.1),1x,f6.3,2f5.1,2f6.1,9f7.1)
         print *,' global mean : asbal  atbal  sfbal    evap    rain ',
     &' hflux colwv  ohflx  flxic  flxml  flxse  flxmj  rebal',
     &'  sbl_l  sbl_i  sbl_s'
         write (6,220) gasbal,gatbal,gsfbal,gev,grn,ghf
     &                ,gcolwv,gochfx,gflxic,gflxml,gflxse,gflxmj,grebal
     &                ,gsbl_l,gsbl_i,gsbl_s
 220     format (14x,3(1x,f6.2),1x,2f8.5,f6.1,f6.2,9f7.2)
         write (6,230)
 230     format (2x,'asbal=(Sin-Sout-Rt); atbal=(Sin-Sout-Sg)',
     &    '+ (Rg-Rt)+Fg+Rn; sfbal=Sg-Rg-Fg-Eg')
         print *,'ohflx= Global mean heat flux down into ocean (model)'
         print *,' at the sea-atmosphere interface and below the ice'
     &    ,' points'
         print *,'flxse= Global mean heat flux down into sea+mlo points'
         print *,'flxic= Global mean heat flux up into ice points'
         print *,'  ==> ohflx= flxse - flxic'
         print *,'flxml= Global mean heat flux up into mlo points'
         print *,'flxmj= Global mean heat flux for MLO:SEA jump'
         print *,'rebal= Heat flux due to Rain-Evap difference'
         print *,'  ==> asbal= atbal + sfbal - rebal'
         print *,'sbl_l= Global mean heat flux down into land'
         print *,'sbl_i= Global mean heat flux down into ice'
         print *,'sbl_s= Global mean heat flux down into sea+mlo'
         print *,'  ==> sfbal= sbl_s + sbl_i + sbl_l'
         print *,' Near ocean heat balance when sbl_s = flxic '

      end if

c---------------------------------------------------------------------
c-------------- PRINT-PLOT HEATING FOR TOP OF ATMOS (R), -------------
c-----------------   ATMOSPHERE (A), AND OCEAN (O) -------------------
      if (plotheat) call aplota (1)

c---------------------------------------------------------------------
c---------------- PHYSICAL STATS -------------------------------------
c---------------------------------------------------------------------
c****
c**** RH,RAD HEATINGS
c****
CSJP      do 240 j=1,lat2*nl
CSJP 240     rhz(j,1,1)=rhz(j,1,1)*consdr
      do 240 k = 1, nl
        do 240 j = 1, 2
          do 240 i = 1, lat
 240        rhz(i,j,k)=rhz(i,j,k)*consdr
      daysec=86400.0
CSJP      do 250 j=1,lat2*nl
CSJP 250     htz(j,1,1)=htz(j,1,1)*(consdr*daysec)
      do 250 k = 1, nl
        do 250 j = 1, 2
          do 250 i = 1, lat
 250        htz(i,j,k)=htz(i,j,k)*(consdr*daysec)
      do 260 ns=1,2
         do 260 lg=1,lat
         do 260 k=1,nl
         grh(k)=grh(k)+rhz(lg,ns,k)*w(lg)
 260     ght(k)=ght(k)+htz(lg,ns,k)*w(lg)
      do 270 k=1,nl
         grh(k)=grh(k)*0.5
 270     ght(k)=ght(k)*0.5

c---------------------------------------------------------------------
c--------------------- Physics 1 Maps --------------------------------
c---------------------------------------------------------------------

      if (phz1p) then

         do 280 k=1,nl
            inl=nint(sig(k)*1000.0)
            write (rhch(k),"(' R',i3.3)") inl
 280        write (htch(k),"(' H',i3.3)") inl
c---- RH map  -----
         print *,' RHZ-MAP  in %'
         write (6,105) (rhch(k),k=1,nl,ksk)
         do 290 ns=1,2
            do 290 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
 290        write (6,107) lx,(rhz(lx,ns,k),k=1,nl,ksk)
         write (6,105) (rhch(k),k=1,nl,ksk)
         write (6,107) lx,(grh(k),k=1,nl,ksk)
c---- HT map  -----
         print *,' HTZ-MAP  in Deg/Day'
         write (6,105) (htch(k),k=1,nl,ksk)
         do 300 ns=1,2
            do 300 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
 300        write (6,'(1x,i2,24f5.2)') lx,(htz(lx,ns,k),k=1,nl,ksk)
         write (6,105) (htch(k),k=1,nl,ksk)
         write (6,'(1x,i2,24f5.2)') lx,(ght(k),k=1,nl,ksk)

      end if

c---------------------------------------------------------------------
c----- CLOUD,RG,SG,TA,T*,T*L,T*S,OCBAL,GWET,RAINS,RAINC --------------
c---------------------------------------------------------------------
CSJP      do 310 j=1,lat2
CSJP         clz(j,1)=clz(j,1)*consrad
CSJP         cllz(j,1)=cllz(j,1)*consrad
CSJP         clmz(j,1)=clmz(j,1)*consrad
CSJP         clhz(j,1)=clhz(j,1)*consrad
CSJP         clcz(j,1)=clcz(j,1)*consdr
CSJP         rgz(j,1)=rgz(j,1)*consdr
CSJP         rtz(j,1)=rtz(j,1)*consrad
CSJP         rnsz(j,1)=rnsz(j,1)*conser
CSJP         rncz(j,1)=rncz(j,1)*conser
CSJP         sgz(j,1)=sgz(j,1)*consdr
CSJP         sinz(j,1)=sinz(j,1)*consrad
CSJP         souz(j,1)=souz(j,1)*consrad
CSJP         albz(j,1)=albz(j,1)*consrad
CSJP         taz(j,1)=taz(j,1)*consdr
CSJP         tsz(j,1)=tsz(j,1)*consdr
CSJP         tslz(j,1)=tslz(j,1)*consdr
CSJP         tssz(j,1)=tssz(j,1)*consdr
CSJPc        sblz(j,1)=sblz(j,1)*consdr
CSJP         gwz(j,1)=gwz(j,1)*consdr
CSJP         rtclrz(j,1)=rtclrz(j,1)*consrad
CSJP         souclrz(j,1)=souclrz(j,1)*consrad
CSJP 310  continue
      do 310 j = 1, 2
        do 310 i = 1, lat
          clz(i,j)=clz(i,j)*consrad
          cllz(i,j)=cllz(i,j)*consrad
          clmz(i,j)=clmz(i,j)*consrad
          clhz(i,j)=clhz(i,j)*consrad
          clcz(i,j)=clcz(i,j)*consdr
          rgz(i,j)=rgz(i,j)*consdr
          rtz(i,j)=rtz(i,j)*consrad
          rnsz(i,j)=rnsz(i,j)*conser
          rncz(i,j)=rncz(i,j)*conser
          sgz(i,j)=sgz(i,j)*consdr
          sinz(i,j)=sinz(i,j)*consrad
          souz(i,j)=souz(i,j)*consrad
          albz(i,j)=albz(i,j)*consrad
          taz(i,j)=taz(i,j)*consdr
          tsz(i,j)=tsz(i,j)*consdr
          tslz(i,j)=tslz(i,j)*consdr
          tssz(i,j)=tssz(i,j)*consdr
c         sblz(i,j)=sblz(i,j)*consdr
          gwz(i,j)=gwz(i,j)*consdr
          rtclrz(i,j)=rtclrz(i,j)*consrad
          souclrz(i,j)=souclrz(i,j)*consrad
 310  continue

      do 320 ns=1,2
         do 320 lg=1,lat
         gcl=gcl+clz(lg,ns)*w(lg)
         gcll=gcll+cllz(lg,ns)*w(lg)
         gclm=gclm+clmz(lg,ns)*w(lg)
         gclh=gclh+clhz(lg,ns)*w(lg)
         gclc=gclc+clcz(lg,ns)*w(lg)
         grg=grg+rgz(lg,ns)*w(lg)
         grt=grt+rtz(lg,ns)*w(lg)
         grtclr=grtclr+rtclrz(lg,ns)*w(lg)
         gsg=gsg+sgz(lg,ns)*w(lg)
         gsin=gsin+sinz(lg,ns)*w(lg)
         gsou=gsou+souz(lg,ns)*w(lg)
         gsouclr=gsouclr+souclrz(lg,ns)*w(lg)
         galb=galb+albz(lg,ns)*w(lg)
         grns=grns+rnsz(lg,ns)*w(lg)
         grnc=grnc+rncz(lg,ns)*w(lg)
         gta=gta+taz(lg,ns)*w(lg)
         gts=gts+tsz(lg,ns)*w(lg)
c**    GSBL,GTSS FOR SEA ONLY
         gsbl=gsbl+sblz(lg,ns)*w(lg)
         gtss=gtss+tssz(lg,ns)*w(lg)
c**   GGW,GTSL FOR LAND ONLY
         ggw=ggw+gwz(lg,ns)*w(lg)
         gtsl=gtsl+tslz(lg,ns)*w(lg)
 320  continue

      do 330 j=2*nl+4,2*nl+24
 330     gevx(j)=gevx(j)*0.5
      gtsl=gtsl/wland
      gtss=gtss/wsea
      gsbl=gsbl/wsea
      ggw=ggw/wland

      do 360 ns=1,2
         do 360 lg=1,lat
         lx=(ns-1)*(lat+1)+lg*(3-2*ns)
         kl=max(1,lon-knl(lx,ns))
         gwzx(lx,ns)=gwz(lx,ns)/kl
         tslzx(lx,ns)=tslz(lx,ns)/kl
         ks=max(1,knl(lx,ns))
         sblzx(lx,ns)=sblz(lx,ns)/ks
         tsszx(lx,ns)=tssz(lx,ns)/ks
 360  continue

c---------------------------------------------------------------------
c--------------------- Physics 2 Maps --------------------------------
c---------------------------------------------------------------------

      if (phz2p) then

         write (6,370)
 370     format ('  cloud  cl  cm  ch  cc   rg     rt     sg   '
     &     ,'solin  solrf  lwclf  swclf   alb'
     &     ,'   tscr    t*     t*l    t*s   ocbal  gwt  rns   rnc')
         swclf=0.0
         do 380 ns=1,2
            do 380 lg=1,lat
            lx=(ns-1)*(lat+1)+lg*(3-2*ns)
            if (clforflag) swclf=souclrz(lx,ns)-souz(lx,ns)
            write (6,390)lx,nint(100.*clz(lx,ns)),nint(100.*cllz(lx,ns))
     &       ,nint(100.*clmz(lx,ns)),nint(100.*clhz(lx,ns))
     &       ,nint(100.*clcz(lx,ns)),rgz(lx,ns),rtz(lx,ns),sgz(lx,ns)
     &       ,sinz(lx,ns),souz(lx,ns),(rtclrz(lx,ns)-rtz(lx,ns)),swclf,
     &       albz(lx,ns),taz(lx,ns),tsz(lx,ns),tslzx(lx,ns),tsszx(lx,ns)
     &       ,sblzx(lx,ns),gwzx(lx,ns),rnsz(lx,ns),rncz(lx,ns)
 380     continue
 390     format (1x,i2,5(1x,i3),7(f6.1,1x),f5.2,5(1x,f6.1),3(1x,f5.2))
         write (6,370)
         gswclf=0.
         if (clforflag) gswclf=(gsouclr-gsou)
         write (6,390) lx,nint(100.*gcl),nint(100.*gcll),nint(100.*gclm)
     &    ,nint(100.*gclh),nint(100.*gclc),grg,grt,gsg,gsin,gsou,
     &    (grtclr-grt),gswclf,galb,gta,gts,gtsl,gtss,gsbl,ggw,grns,grnc

      end if

      if(gwicm)then
        gicsum=0.0
        wicsum=0.0
        gicsumn=0.0
        wicsumn=0.0
        gicsums=0.0
        wicsums=0.0
        do lg=1,lat
        if(jice(lg,1).gt.0)then
          gicsum=gicsum+zisum(lg,1)*w(lg)
          gicsumn=gicsumn+zisum(lg,1)*w(lg)
          wicsum=wicsum+jice(lg,1)*w(lg)
          wicsumn=wicsumn+jice(lg,1)*w(lg)
          zisum(lg,1)=zisum(lg,1)/jice(lg,1)
        endif
        if(jice(lg,2).gt.0)then
          gicsum=gicsum+zisum(lg,2)*w(lg)
          gicsums=gicsums+zisum(lg,2)*w(lg)
          wicsum=wicsum+jice(lg,2)*w(lg)
          wicsums=wicsums+jice(lg,2)*w(lg)
          zisum(lg,2)=zisum(lg,2)/jice(lg,2)
        endif
        enddo
        gicsum=gicsum/wicsum
        gicsumn=gicsumn/wicsumn
        gicsums=gicsums/wicsums
        write(6,1800)gicsum
 1800 format('ICE GLOBAL=',f7.2)
        write(6,1801)gicsumn
 1801 format('ICE NORTH =',f7.2)
c       do lg=1,lat
c         if(zisum(lg,1).gt.0.0)write(6,1802)lg,zisum(lg,1)
c       enddo
c 1802 format(i2,2x,f7.2)
        write(6,1803)gicsums
 1803 format('ICE SOUTH =',f7.2)
c       do lg=lat,1,-1
c         if(zisum(lg,2).gt.0.0)write(6,1802)lg,zisum(lg,2)
c       enddo
      endif

c---------------------------------------------------------------------
c-------------- PRINT-PLOT CLOUD  PROFILES ---------------------------
      if (plotclds) call aplota (2)
c---------------------------------------------------------------------
c---- PRINT-PLOT Top of Atmos radiation (Nett sol in, Nett OLR out) --
      if (plotnetr) call aplota (3)
c---------------------------------------------------------------------
c----------- PRINT-PLOT EVAP/RAIN PROFILES ---------------------------
      if (plotevrn) call aplota (4)
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c---------------- gbrfpr.f OUTPUT SCANNING ---------------------------
c---------------------------------------------------------------------

      if (savegbrf) then

         specstat(13)=gasbal
         specstat(14)=gatbal
         specstat(15)=gsfbal
         specstat(16)=gev
         specstat(17)=grn
         specstat(18)=ghf
         specstat(33)=rainmx
         do 400 k=1,nl
 400        specstat(33+k)=grh(k)
         specstat(19)=gcl*100.
         specstat(20)=gcll*100.
         specstat(21)=gclm*100.
         specstat(22)=gclh*100.
         specstat(23)=grt
         specstat(24)=gsou
         specstat(25)=grtclr-grt
         specstat(26)=gsouclr-gsou
         specstat(27)=gts
         specstat(28)=gtsl
         specstat(29)=gtss
         specstat(30)=gsbl
         specstat(31)=ggw
         specstat(32)=grnc

c    find min sub Antarctic pressure

         pmin=1.e30
         latmin=0
c---- lsx,lnx are G-lat range for required feature - selected via R21 mo
         lsx=int((7.1/28)*lat)
         lnx=int((12.1/28)*lat)
c     print *,' searching fpr pmin between SH lats ',lsx,lnx
         do 410 lg=lsx,lnx
            if (pmin.lt.ps(lg,2)) go to 410
            pmin=ps(lg,2)
            latmin=-lg
 410     continue
         specstat(nl+34)=pmin
         specstat(nl+35)=float(latmin)

c    select peak jets and easterly band limits

c    for nhem find single peak
c---- check for NH  Zonal Wind Jets - max of 10
c---- at tropopause level
         lnor=lat/2
c---- determine tropopause levels (sigma values 0.35 to 0.15)
c---- kbt(bottom),ktt (top)
         do 420 k=1,nl
            if (sig(k).lt.0.35) go to 430
 420     continue
 430     kbt=k
         do 440 k=nl,1,-1
            if (sig(k).gt.0.15) go to 450
 440     continue
 450     ktt=k
         ns=1
         do 460 k=1,10
            umxt(k)=-1.0e+30
 460        umx(k)=-1.0e+30
         ki=0
c     do 140 k=2,nl-1
c        do 140 lg=2,lat-1
c---- lsx,lnx are G-lat range for required feature - selected via
c---- R21 model
         lsx=int((12.1/28)*lat)
         lnx=int((22.1/28)*lat)
         do 470 k=kbt,ktt
            do 470 lg=lsx,lnx
            xx=u(lg,ns,k)
            if ((xx.gt.u(lg+1,ns,k)).and.(xx.gt.u(lg-1,ns,k)).and.
     &       (xx.gt.u(lg,ns,k+1)).and.(xx.gt.u(lg,ns,k-1))) then
               ki=ki+1
               umx(ki)=xx
               umxt(ki)=xx
               lgu(ki)=(3-2*ns)*lg
               lgut(ki)=lgu(ki)
               if (ki.eq.10) go to 480
            end if
 470     continue
c---- sort the U max values (if more than 1 (quasi) jet)
 480     if (ki.gt.1) then
            do 500 kis=1,ki
               kx=0
               umax=-1.0e+30
               do 490 k=1,10
                  if (umxt(k).gt.umax) then
                     umax=umxt(k)
                     kx=k
                  end if
 490           continue
               kpnt(kis)=kx
               umxt(kx)=-1.0e+30
 500        continue
            do 510 kis=1,ki
               kx=kpnt(kis)
               umxt(kis)=umx(kx)
               lgut(kis)=lgu(kx)
 510        continue
         end if
c     print *,' There are ',ki,' NH jets'
c---- choose max NH jet only
         specstat(1)=umxt(1)
         specstat(2)=float(lgut(1))

c    for shem find one or two peaks

c---- check for SH  Zonal Wind Jets - max of 10
c---- between sigma=0.35 and top of atmos
         do 520 k=1,nl
            do 520 ns=1,2
            do 520 lg=1,lat
 520        ut(lg,ns,k)=u(lg,ns,k)
c---- checking up to level nl : need boundary value to get strat jet in
         do 530 ns=1,2
            do 530 lg=1,lat
 530        ut(lg,ns,nlp)=u(lg,ns,nl)*0.9
         ns=2
         do 540 k=1,10
            umxt(k)=-1.0e+30
 540        umx(k)=-1.0e+30
         ki=0
c---- lsx,lnx are G-lat range for required feature - selected via
c---- R21 model
         lsx=int((7.1/28)*lat)
         lnx=int((26.1/28)*lat)
         do 550 k=kbt,nl
            do 550 lg=lsx,lnx
            xx=u(lg,ns,k)
            if ((xx.gt.ut(lg+1,ns,k)).and.(xx.gt.ut(lg-1,ns,k)).and.
     &       (xx.gt.ut(lg,ns,k+1)).and.(xx.gt.ut(lg,ns,k-1))) then
               ki=ki+1
               umx(ki)=xx
               umxt(ki)=xx
               lgu(ki)=(3-2*ns)*lg
               lgut(ki)=lgu(ki)
               if (ki.eq.10) go to 560
            end if
 550     continue
c---- sort the U max values (if more than 1 (quasi) jet)
 560     if (ki.gt.1) then
            do 580 kis=1,ki
               kx=0
               umax=-1.0e+30
               do 570 k=1,10
                  if (umxt(k).gt.umax) then
                     umax=umxt(k)
                     kx=k
                  end if
 570           continue
               kpnt(kis)=kx
               umxt(kx)=-1.0e+30
 580        continue
            do 590 kis=1,ki
               kx=kpnt(kis)
               umxt(kis)=umx(kx)
               lgut(kis)=lgu(kx)
 590        continue
         end if
c     print *,' There are ',ki,' SH jets'
         specstat(7)=0.0
         specstat(8)=0.0
         ki=min(2,ki)
         do 600 kis=1,ki
            specstat(3+kis*2)=umxt(kis)
 600        specstat(4+kis*2)=float(lgut(kis))

c     find peak v

c---- check for Meridional Wind maximums and minimums
c---- at tropopause and between 45N and 45S  (max of 10)
         lath=lat/2
         lnor=lath
         lsou=lat+lath+1
c     print *,' V(max,min) between lats ',lnor,lsou,' at levs ',kbt,ktt
         do 610 k=1,nl
            do 610 lg=1,lat
            vt(lg,k)=v(lg,1,k)
 610        vt(lat2p-lg,k)=v(lg,2,k)
c---- Do two passes : 1 for maximums, 2 for minimums
         do 690 kpass=1,2
            do 620 k=1,10
               umxt(k)=0.0
 620           umx(k)=0.0
            ki=0
            do 630 k=kbt,ktt
               do 630 lg=lnor,lsou
               xx=vt(lg,k)
               testmax=.false.
               if ((xx.gt.vt(lg+1,k)).and.(xx.gt.vt(lg-1,k)).and.(xx.gt.
     &          vt(lg,k+1)).and.(xx.gt.vt(lg,k-1))) testmax=.true.
               testmin=.false.
               if ((xx.lt.vt(lg+1,k)).and.(xx.lt.vt(lg-1,k)).and.(xx.lt.
     &          vt(lg,k+1)).and.(xx.lt.vt(lg,k-1))) testmin=.true.
               tester=testmax
               if (kpass.eq.2) tester=testmin
               if (tester) then
                  ki=ki+1
                  umx(ki)=xx
                  umxt(ki)=xx
                  lgu(ki)=lg
                  if (lg.gt.lat) lgu(ki)=lg-lat2p
                  lgut(ki)=lgu(ki)
                  if (ki.eq.10) go to 640
               end if
 630        continue
c---- sort the V max/min values (if more than 1 (quasi) jet)
 640        if (ki.gt.1) then
               do 670 kis=1,ki
                  kx=0
                  if (kpass.eq.1) then
                     vmax=0.0
                     do 650 k=1,10
                        if (umxt(k).gt.vmax) then
                           vmax=umxt(k)
                           kx=k
                        end if
 650                 continue
                  end if
                  if (kpass.eq.2) then
                     vmin=0.0
                     do 660 k=1,10
                        if (umxt(k).lt.vmin) then
                           vmin=umxt(k)
                           kx=k
                        end if
 660                 continue
                  end if
                  kpnt(kis)=kx
                  umxt(kx)=0.0
 670           continue
               do 680 kis=1,ki
                  kx=kpnt(kis)
                  umxt(kis)=umx(kx)
                  lgut(kis)=lgu(kx)
 680           continue
            end if
c---- store the max V jet
            if (kpass.eq.1) then
               specstat(9)=umxt(1)
               specstat(10)=float(lgut(1))
            else
               specstat(11)=umxt(1)
               specstat(12)=float(lgut(1))
            end if
 690     continue

c---- search for easterly jet minimum extent between 45N and 45S
c---- in the mid levels of the atmosphere (near equator)

         lath=lat/2
         lnor=lath
         lsou=lat+lath+1
c---- determine mid levels of atmos (sigma values 0.7 to 0.15)
         do 700 k=1,nl
            if (sig(k).lt.0.7) go to 710
 700     continue
 710     kbm=k
         do 720 k=nl,1,-1
            if (sig(k).gt.0.15) go to 730
 720     continue
 730     ktm=k
c     print *,' V extent between lats ',lnor,lsou,' at levs ',kbm,ktm
         do 750 k=kbm,ktm
            do 740 lg=1,lat2
 740           jetex(lg,k)=0
            do 750 lg=1,lat
            if (u(lg,1,k).lt.0.0) jetex(lg,k)=lg
 750        if (u(lg,2,k).lt.0.0) jetex(lat2p-lg,k)=-lg
c---- check if any level has no easterly wind in tropics
         ksmin=lat
         kx=0
         lnx=0
         lsx=0
         do 770 k=kbm,ktm
            ksum=0
            do 760 lg=lnor,lsou
 760           if (jetex(lg,k).ne.0) ksum=ksum+1
c---- if ksum=0, then there is no zone of Easterly wind
c----  at this level - jump out of checking
            if (ksum.eq.0) go to 790
            if (ksum.lt.ksmin) then
               kx=k
               ksmin=ksum
            end if
 770     continue
c---- level kx has minimum extent of equatorial Easterly jet
         do 780 lg=lnor,lsou
            if ((jetex(lg,kx).eq.0).and.(jetex(lg+1,kx).ne.0)) lnx=
     &       jetex(lg+1,kx)
            if ((jetex(lg,kx).ne.0).and.(jetex(lg+1,kx).eq.0)) lsx=
     &       jetex(lg,kx)
 780     continue
 790     specstat(3)=float(lnx)
         specstat(4)=float(lsx)

         write (69,800) (specstat(ii),ii=1,12)
 800     format (1x,f7.2,f5.0,2f5.0,4(f7.2,f5.0))
         write (69,810) (specstat(ii),ii=13,22),specstat(33)
 810     format (1x,3f7.2,2f10.5,f7.2,4f5.1,f7.2)
         write (69,820) (specstat(ii),ii=23,32)
 820     format (1x,7f7.2,3f8.3)
         write (69,830) (specstat(ii),ii=34,33+nl)
 830     format (9(1x,f5.1))
         write (69,840) (specstat(ii),ii=34+nl,35+nl)
 840     format (4x,f7.2,f5.0)

      end if

c---------------------------------------------------------------------
c------------ Special diagnostics for coupled run --------------------
c---------------------------------------------------------------------

      if (lcouple) then

c        call ocwrit9 (pvbtb,nl,'pvbtb   ',itt)
c        call ocwrit9 (pvtb, nl,'pvtb    ',itt)
c        call ocwrit9 (pvbub,nl,'pvbub   ',itt)
c        call ocwrit9 (pvub, nl,'pvub    ',itt)
c        call ocwrit9 (pvbqb,nl,'pvbqb   ',itt)
c        call ocwrit9 (pvqb, nl,'pvqb    ',itt)
c        call ocwrit9 (psv,  nl,'psv     ',itt)
c        call ocwrit9 (pssd,nlm,'pssd    ',itt)
c        call ocwrit9 (u,    nl,'u       ',itt)
c        call ocwrit9 (v,    nl,'v       ',itt)
c        call ocwrit9 (t,    nl,'t       ',itt)
c        call ocwrit9 (ww,  nlm,'sd      ',itt)

c        call ocwrite (ps,      'psl     ',itt)
c        call ocwrite (clz,     'cloud   ',itt)
c        call ocwrite (cllz,    'cll     ',itt)
c        call ocwrite (clmz,    'clm     ',itt)
c        call ocwrite (clhz,    'clh     ',itt)
c        call ocwrite (evz,     'evap    ',itt)
c        call ocwrite (rnz,     'rain    ',itt)
c        call ocwrite (hfz,     'hflux   ',itt)

c        call ocwrit9 (rhz,  nl,'rh      ',itt)
c        call ocwrit9 (htz,  nl,'heating ',itt)

c        call ocwrite (rgz,     'rg      ',itt)
c        call ocwrite (rtz,     'rt      ',itt)
c        call ocwrite (rnsz,    'rains   ',itt)
c        call ocwrite (rncz,    'rainc   ',itt)
c        call ocwrite (sgz,     'sg      ',itt)
c        call ocwrite (sinz,    'solin   ',itt)
c        call ocwrite (souz,    'solrf   ',itt)
c        call ocwrite (albz,    'alb     ',itt)
c        call ocwrite (taz,     'tscrn   ',itt)
c        call ocwrite (tsz,     't*      ',itt)
c        call ocwrite (tslzx,   't*l     ',itt)
c        call ocwrite (tsszx,   't*s     ',itt)
c        call ocwrite (sblzx,   'ocbal   ',itt)
c        call ocwrite (gwzx,    'gwet    ',itt)
c        call ocwrite (asbalz  ,'energyin',itt)
c        call ocwrite (atbalz  ,'energyat',itt)
c        call ocwrite (sfbalz  ,'energysf',itt)
c        call ocwrite (rtclrz  ,'rtclr   ',itt)
c        call ocwrite (souclrz ,'souclr  ',itt)

      end if

c---------------------------------------------------------------------
c--------------------- Cif Dump --------------------------------------
c---------------------------------------------------------------------
c Dump various zonal mean diagnostics in convenient form for plotting
c with con_cif. These are equivalent to HBG's zonal mean prints.
      if (savezont) then
         header='Zonal mean temperature                            '
         call dumpzcif (t,header,'zont'//str,glats,sig)
      end if
      if (savezonr) then
         header='Zonal mean relative humidity                      '
         call dumpzcif (rhz,header,'zonr'//str,glats,sig)
      end if
      if (savezonu) then
         header='Zonal mean zonal wind                             '
         call dumpzcif (u,header,'zonu'//str,glats,sig)
      end if
      if (savezonv) then
         header='Zonal mean meridional wind                        '
         call dumpzcif (v,header,'zonv'//str,glats,sig)
      end if
      if (qcloud.and.savezqcl) then
         header='Zonal mean cloud liquid water mixing ratio (g/kg) '
         call dumpzcif (qlz,header,'zqcl'//str,glats,sig)
         header='Zonal mean cloud ice mixing ratio (g/kg)          '
         call dumpzcif (qfz,header,'zqcf'//str,glats,sig)
         header='Zonal mean evaporation of rain (10**-9 kg/kg/s)   '
         call dumpzcif (qevapz,header,'zqev'//str,glats,sig)
         header='Zonal mean sublimation of snow (10**-9 kg/kg/s)   '
         call dumpzcif (qsublz,header,'zqsb'//str,glats,sig)
         header='Zonal mean autoconversion of ql (10**-9 kg/kg/s)  '
         call dumpzcif (qautoz,header,'zaut'//str,glats,sig)
         header='Zonal mean collection of ql by rain 10**-9 kg/kg/s'
         call dumpzcif (qcollz,header,'zcol'//str,glats,sig)
         header='Zonal mean accretion of ql by snow 10**-9 kg/kg/s '
         call dumpzcif (qaccrz,header,'zacr'//str,glats,sig)
      end if
      if (savezcfr) then
         header='Zonal mean cloud fraction                         '
         call dumpzcif (cfracz,header,'zcfr'//str,glats,sig)
      end if
      if (savezonw) then
         open(51,file='zonw'//str)
         write (51,870) nl
         write (51,880) (sigh(k)*1000,k=1,nl)
         write (51,870) lat2
         write (51,880) (glats(lgns)*180.0/pi,lgns=1,lat2)
         write (51,'(a)') 'Zonal mean d(sigma)/dt'
         do 850 ns=1,2
            do 850 lg=1,lat
            lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
 850        qout(lgns,1)=0. !Pad k=1 with zeroes
         do 860 k=2,nl
            do 860 ns=1,2
            do 860 lg=1,lat
            lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
 860        qout(lgns,k)=ww(lg,ns,k-1)
         write (51,880) qout
         close (51)
      end if
 870  format (i4)
 880  format (1p6e12.5)
c**** RESET ZONAL STATS ARRAYS TO ZERO IN SUB "ZEROST"
      return
      end
C---------------------------------------------------------------------
      subroutine dumpzcif(zfield,header,filename,glats,sig)

c Dump zonal mean cif files.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real zfield(lat,2,nl)
      character*50 header
      character*17 filename
      real glats(lat2)
      real sig(nl)

C Local work arrays and variables
      real zout(lat2,nl)

      integer k
      integer lg
      integer lgns
      integer ns

C Start code : ----------------------------------------------------------

      open(51,file=filename)
      write(51,933)nl
      write(51,934)(sig(k)*1000,k=1,nl)
      write(51,933)lat2
      write(51,934)(glats(lgns)*180.0/pi,lgns=1,lat2)
      write(51,'(a50)')header

      do k=1,nl
        do ns=1,2
          do lg=1,lat
            lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
            zout(lgns,k)=zfield(lg,ns,k)
          enddo
        enddo
      enddo
      write(51,934)zout
      close(unit=51)

 933  format(i4)
 934  format (1p6e12.5)
      
      return
      end
