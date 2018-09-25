c (1) Declaration of /FCORR/ moved to FCORR.f.
c (2) Modified to reflect the fact that SFLUXES now contains 14 fields, rather
c     than 7.
c SJP 2004/09/10
c
c Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c heat input, to be saved.
c SJP 2003/06/20
c
c Modified for the changes to the freshwater fluxes, whereby the ice water flux
c is split into its two components.
c SJP 2003/06/19
c
c $Log: zerost.f,v $
c Revision 1.27  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.26  2001/03/07 04:28:57  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.25  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.24  2001/02/12 05:39:43  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.23  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.22  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.21  1998/12/10  00:55:29  ldr
c HBG changes to V5-1-21
c
c Revision 1.20  1998/01/30  04:36:13  ldr
c Fix from SPO for Cd diagnostic.
c
c Revision 1.19  1997/12/23  00:33:38  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.18  1997/12/19  02:03:10  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.17  1994/08/08  17:23:23  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.16  94/08/04  16:56:50  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.15  94/03/30  12:35:46  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.14  93/12/17  15:34:28  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.13  93/11/03  11:44:41  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.12  93/10/05  13:07:51  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.11  93/09/21  12:43:10  ldr
c Moved initialization of max and min temperatures from filerd to zerost
c where they belong.
c 
c Revision 1.10  93/09/17  12:18:26  mrd
c Added diagnostic for monthly mean atmospheric latent heating
c 
c Revision 1.9  93/08/06  16:14:19  ldr
c Introduce new flag savefcor to control saving of flux correction data.
c 
      subroutine zerost(kzst)

      implicit none
C Global parameters
      include 'PARAMS.f'
      integer ndynav,nprtar,nprtar1
      parameter (ndynav=lat*2*(3*nl+nlm+1))
      parameter (nprtar=lon*lat*2*3)
      parameter (nprtar1=lon*lat*2*4)

C Argument list
      integer kzst

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'BGRND.f'
      include 'ECPARM.f'
      include 'ECTRAC.f'
      include 'FCORR.f'
      include 'FEWFLAGS.f'
      include 'MASIV2.f'
      include 'MASIV4.f'
      include 'RADAV.f'
      include 'STFLAGS.f'
      include 'SURF1.f'
      include 'VERTV.f'        !Coupled model

      real cdavg,pcdavg
      common/cdmapx/cdavg(ln2,lat),pcdavg(ln2,lat)

      real ochf,occur,g50dt
      integer khf
      common/curent/ochf(ln2,lat),occur(ln2,lat),khf(ln2,lat)
     &,g50dt(ln2,lat)

      real u,v,t,ww,ps
      common/dynav/u(lat,2,nl),v(lat,2,nl),t(lat,2,nl)
     &,ww(lat,2,nlm),ps(lat,2)
        real ux(ndynav)
        equivalence (ux(1),u(1,1,1))

      real fwfpme,fwficea,fwficeb,fwfriv
      common/fwfdat/fwfpme(ln2,lat),fwficea(ln2,lat),fwficeb(ln2,lat),
     &              fwfriv(ln2,lat)

      real hlat_a,clat_a
      common/latheat/hlat_a(ln2,lat,nl),clat_a(ln2,lat,nl)

      real tgmax,tgmin,tsmax,tsmin,htsmax,htsmin
     & ,extsmax,extsmin,tggmax,tggmin,tgfmax,tgfmin
     & ,htggmax,htggmin,htgfmax,htgfmin,hpmc
      common/maxmin/tgmax(ln2,lat),tgmin(ln2,lat),
     &              tsmax(ln2,lat),tsmin(ln2,lat),
     &              htsmax(ln2,lat),htsmin(ln2,lat),
     &              extsmax(ln2,lat),extsmin(ln2,lat),
     &              tggmax(ln2,lat),tggmin(ln2,lat),
     &              tgfmax(ln2,lat),tgfmin(ln2,lat),
     &              htggmax(ln2,lat),htggmin(ln2,lat),
     &              htgfmax(ln2,lat),htgfmin(ln2,lat)
     &             ,hpmc(ln2,lat)

      real statsice
      common/masiv5/statsice(ln2,14,lat)

      real devap,drain,dvgg
      common/prtar/devap(ln2,lat),drain(ln2,lat),dvgg(ln2,lat)
        real devapx(nprtar)
        equivalence (devapx(1),devap(1,1))

      real sshf,slhf,sswr,slwr
      common/prtar1/sshf(ln2,lat),slhf(ln2,lat),sswr(ln2,lat)
     & ,slwr(ln2,lat)
        real sshfx(nprtar1)
        equivalence (sshfx(1),sshf(1,1))

      real brain
      common/rblock/brain(ln2,lat)

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

C Local work arrays and variables

      integer k
      integer lg
      integer mg
      integer nv

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** TO ZERO PHYSICAL STATISTIC ARRAYS FOR DATA FILE 21

      if(kzst.eq.0)then
        do 10 lg=1, lat
        do 10 k=14,max_radstm
        do 10 mg=1,ln2
          radstm(mg,k,lg)=0.0
   10   continue

        if(semice)then
        do 145 lg=1,lat
        do 145 k=10,14
        do 145 mg=1,ln2
          statsice(mg,k,lg)=0.0
  145   continue
        endif

c ZERO VARIOUS ARRAYS

        kdynm=0
        do 15 lg=1,lat
        do 15 mg=1,ln2
          brain(mg,lg)=0.0
          cdavg(mg,lg)=0.0
          pcdavg(mg,lg)=0.0
          pmonth(mg,lg)=0.0
          ochf(mg,lg)=0.0
          occur(mg,lg)=0.0
          khf(mg,lg)=0
          tsmax(mg,lg)=0.               
          htsmax(mg,lg)=0.
          htsmin(mg,lg)=0.
          extsmax(mg,lg)=0.
          tsmin(mg,lg)=999.             
          extsmin(mg,lg)=999.
c These ones are from NSIB scheme
          tgfmax(mg,lg)=tgf(mg,lg)
          tgfmin(mg,lg)=tgf(mg,lg)
          htggmax(mg,lg)=0.
          htgfmax(mg,lg)=0.
          htggmin(mg,lg)=0.
          htgfmin(mg,lg)=0.
          hpmc(mg,lg)=0.
   15   continue 
        do 151 k=1,nbelow
        do 151 lg=1,lat
        do 151 mg=1,ln2
          ibeta(mg,lg,k)=0
  151   continue 
        do 17 k=1,nl
        do 17 lg=1,lat
        do 17 mg=1,ln2
          tmonth(mg,lg,k)=0.0
          umonth(mg,lg,k)=0.0
          vmonth(mg,lg,k)=0.0
          qmonth(mg,lg,k)=0.0
          rmonth(mg,lg,k)=0.0
          gmonth(mg,lg,k)=0.0
          hlat_a(mg,lg,k)=0.0
          clat_a(mg,lg,k)=0.0
   17   continue 

C**** TO ZERO MERIDIONAL TRANSPORT ARRAYS FOR DATA FILE 21
        do 23 k=1,nvertv
 23     pvbtbx(k)=0.0

C**** ZERO ARRAY HOLDING ATMOSPHERIC FLUXES USED TO COMPUTE FLUX CORRECTIONS
C**** FOR COUPLED MODEL
        
        if(savefcor)then
          do nv=1,14
            do lg=1,lat2
              do mg=1,lon
                sfluxes(mg,lg,nv)=0.
              enddo
            enddo
          enddo
        endif

C**** River routing stats
      if(newriver)then
        do lg=1,lat2
        do mg=1,lon
          srunoff(mg,lg)=0.0
          sresvr1(mg,lg)=0.0
        enddo
        enddo
      endif

c**** Fresh water flux stats
      if(fwf_sflg)then
        do lg=1,lat
        do mg=1,ln2
          fwfpme(mg,lg)=0.0
          fwficea(mg,lg)=0.0
          fwficeb(mg,lg)=0.0
          fwfriv(mg,lg)=0.0
        enddo
        enddo
      endif

c SET MAX AND MIN TEMPS TO VALUES OF SURFACE TEMP READ FROM RESTART FILE

        do 96 lg=1,lat
        do 96 mg=1,ln2
          tgmax(mg,lg)=savegrid(mg,3,lg)
          tgmin(mg,lg)=savegrid(mg,3,lg)
c These ones are from NSIB scheme
          tggmax(mg,lg)=tggsl(mg,1,lg)
          tggmin(mg,lg)=tggsl(mg,1,lg)
   96   continue

c For aerosol scheme...

        if(coupled_aero)then
          hxtg(:,:,:,:)=0.
        endif

      endif ! end of if(kst.eq.0)

C**** THE FOLLOWING DONE IN ALL CALLS TO ZEROST
C**** TO ZERO ZONAL STATS ARRAYS
      do 24 k=1,ndynav
 24   ux(k)=0.0
      do 25 k=1,nradav
 25   clzx(k)=0.0
      do 26 k=1,nprtar
 26   devapx(k)=0.0
      do 28 k=1,nprtar1
 28   sshfx(k)=0.0

      return
      end
