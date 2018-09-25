c Modified for five-character experiment names.
c SJP 2009/04/22
c
c $Log: prtt.f,v $
c Revision 1.10  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.9  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.8  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.7  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.6  93/12/17  15:33:28  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.5  93/06/18  13:07:50  ldr
c Modify Hal's printed map of PMSL so that it can be plotted with con_cif.
c 
c Revision 1.4  92/11/20  12:14:06  ldr
c Changes for extra printouts/ corrections.
c 
c Revision 1.3  92/04/23  12:12:29  ldr
c Diagnostic routines generalized by HBG for R42.
c 
c Revision 1.2  91/03/13  12:59:38  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:51  ldr
c Initial release V3-0
c 
      subroutine prtt(ndadd)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ndadd

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FILES.f'
      include 'PRINTT.f'
      include 'TIMEX.f'

      real devap,drain,dvgg
      common/prtar/devap(ln2,lat),drain(ln2,lat),dvgg(ln2,lat)

      real sshf,slhf,sswr,slwr
      common/prtar1/sshf(ln2,lat),slhf(ln2,lat),sswr(ln2,lat),
     &slwr(ln2,lat)

C Local work arrays and variables
      real tx(lon,lat2)
      character*17 fname

      integer lg
      integer lgns
      integer loops
      integer lx
      integer ma
      integer maxnum
      integer mg
      integer mga
      integer mgb
      integer ns

      real conx
      real cony

C Local data, functions etc

C Start code : ----------------------------------------------------------

C---- Global maps - Optional - Rain,Evap, Pmsl (or Upper level divergence)
C---- Global maps - Optional - SSHF,SLHF,SSWR,SLWR
      conx=mstep/1440.0/ndadd
      cony=1.0/ndadd

      if (.not.gmap2)return
c---- maximum number of integers/characters per line is 32
      maxnum=32
      loops=(lon-1)/maxnum+1

      if (.not.rainm) go to 30
c---- rainmap
      do 20 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         write (6,10) ndadd,mga,mgb
 10      format (1x,i2,' day average rain : mg=',i3,' to ',i3)
 20      call prtt1 (drain,cony,mga,mgb)

 30   if (.not.evapm) go to 60
c---- evap map
      do 50 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         write (6,40) ndadd,mga,mgb
 40      format (1x,i2,' day average evap(mms*10) : mg=',i3,' to ',i3)
 50      call prtt1 (devap,10.0*cony,mga,mgb)

c---- surface pressure map or 350mb divergence map
c Replace old printed map with plotable file in normal "cif" format
c Old code is commented out below.
 60   if (.not.pmslm) go to 120
      do 70 lgns=1,lat2
         ns=2-(lgns-1)/lat
         lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
           do 70 mg=1,lon
             ma=mg+(ns-1)*lon
             tx(mg,lgns)=dvgg(ma,lg)
 70        continue
           fname='pmsl'//str
           print*,'Writing instantaneous pmsl to file ',fname
           open(51,file=fname)
             write(51,9040)tx
           close(51)
 9040 format(1p6e12.5)
           
c---- If surface pressure, remove mean values
c**   PMSL MAPS (MBS WITH 1000 OR 900 MBS REMOVED)
C***         zzz=dvgg(mg,lg)
C***         if (zzz.ge.999.99) zzz=zzz-1000.0
C***         if (zzz.ge.900.0) zzz=zzz-900.0
C*** 70      dvgg(mg,lg)=zzz
C***      do 90 lx=1,loops
C***         mgb=lx*maxnum
C***         mga=mgb-maxnum+1
C***         mgb=min(lon,mgb)
C***         write (6,80)mga,mgb
C*** 80      format (1x,'instantaneous p*(-1000 or -900)  : mg=',i3,' to ',
C***     1    i3)
C*** 90      call prtt1 (dvgg,1.0,mga,mgb)
c**   DIVERGENCE MAP (UNITS 1.0E-06/SEC) AT 350MB
c     do 110 lx=1,loops
c        mgb=lx*maxnum
c        mga=mgb-maxnum+1
c        mgb=min(lon,mgb)
c        write (6,100) ndadd,mga,mgb
c100     format (1x,'instantaneous p*(-1000 or -900)  : mg=',i3,' to ',
c    1    i3)
c110     call prtt1 (dvgg,conx*1.0e+06,mga,mgb)

 120  if (.not.surfm) return
c---- print surface sensible heat flux
c---- some with coded output : e.g. A49=149, C65=365, F 1=501
      do 140 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         write (6,130) ndadd,mga,mgb
 130     format (1x,i2,' day average surface sensible heat flux : mg='
     & ,i3,' to ',i3)
 140     call prtt3 (sshf,conx,mga,mgb)
c---- print surface latent heat flux
      do 160 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         write (6,150) ndadd,mga,mgb
 150     format (1x,i2,' day average latent heat flux : mg=',i3,' to '
     & ,i3,
     & ' : coded output : e.g. A49=149, C65=365, F 1=501')
 160     call prtt2 (slhf,conx,mga,mgb)
c---- print surface solar heat flux
      do 180 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         write (6,170) ndadd,mga,mgb
 170     format (1x,i2,' day average surface solar flux : mg=',i3,
     & ' to ',i3,
     & ' : coded output : e.g. A49=149, C65=365, F 1=501')
 180     call prtt2 (sswr,conx,mga,mgb)
c---- print surface long wave flux
      do 200 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         write (6,190) ndadd,mga,mgb
 190     format (1x,i2,' day average surface LW flux:  mg=',i3,' to '
     & ,i3,
     & ' : coded output : e.g. A49=149, C65=365, F 1=501')
 200     call prtt2 (slwr,conx,mga,mgb)
      return
      end
C---------------------------------------------------------------------
      subroutine prtt1 (arrx,con,mga,mgb)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real arrx(ln2,lat)
      real con
      integer mga
      integer mgb

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'LSMI.f'

C Local work arrays and variables
      integer iar(32)
      character*1 car(32)

      integer mg
      integer lg
      integer mx

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- maxnumber for this routine is 32
c---- Mapping routine for PRTT - Global maps with "L" to denote land
c---- (A1,I2)
      write (6,10) (mg,mg=mga,mgb)
 10   format (1x,' mg=',32i3)
      do 30 lg=1,lat
         mx=0
         do 20 mg=mga,mgb
            mx=mx+1
            car(mx)=' '
            if (imsl(mg,lg).eq.4) car(mx)='L'
 20         iar(mx)=nint(arrx(mg,lg)*con)
 30      write (6,40) lg,(car(mg),iar(mg),mg=1,mx)
 40   format (1x,i2,2x,32(a1,i2))
      do 60 lg=lat,1,-1
         mx=0
         do 50 mg=mga,mgb
            mx=mx+1
            car(mx)=' '
            if (imsl(mg+lon,lg).eq.4) car(mx)='L'
 50         iar(mx)=nint(arrx(mg+lon,lg)*con)
 60      write (6,40) lg,(car(mg),iar(mg),mg=1,mx)
      write (6,10) (mg,mg=mga,mgb)

      return
      end
C---------------------------------------------------------------------
      subroutine prtt2 (arrx,con,mga,mgb)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real arrx(ln2,lat)
      real con
      integer mga
      integer mgb

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer iar(32)
      character*1 car(32)

      integer iarx
      integer ihun
      integer lg
      integer mg
      integer mx

C Local data, functions etc
      character*1 cint(10)
      data cint/' ','A','B','C','D','E','F','G','H','I'/

C Start code : ----------------------------------------------------------

c---- maxnumber for this routine is 32
c---- Mapping routine for PRTT - Global maps  (I3)
      write (6,10) (mg,mg=mga,mgb)
 10   format (1x,' mg=',32i3)
      do 30 lg=1,lat
         mx=0
         do 20 mg=mga,mgb
            mx=mx+1
            iarx=nint(arrx(mg,lg)*con)
            ihun=min(9,iarx/100)
            car(mx)=cint(ihun+1)
 20         iar(mx)=iarx-ihun*100
 30      write (6,40) lg,(car(mg),iar(mg),mg=1,mx)
 40   format (1x,i2,2x,32(a1,i2))
      do 60 lg=lat,1,-1
         mx=0
         do 50 mg=mga,mgb
            mx=mx+1
            iarx=nint(arrx(mg+lon,lg)*con)
            ihun=min(9,iarx/100)
            car(mx)=cint(ihun+1)
 50         iar(mx)=iarx-ihun*100
 60      write (6,40) lg,(car(mg),iar(mg),mg=1,mx)
         write (6,10) (mg,mg=mga,mgb)

      return
      end
C---------------------------------------------------------------------
      subroutine prtt3 (arrx,con,mga,mgb)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real arrx(ln2,lat)
      real con
      integer mga
      integer mgb

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer iar(32)

      integer lg
      integer mg
      integer mx

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- maxnumber for this routine is 32
c---- Mapping routine for PRTT - Global maps  (I3)
      write (6,10) (mg,mg=mga,mgb)
 10   format (1x,' mg=',32i3)
      do 30 lg=1,lat
         mx=0
         do 20 mg=mga,mgb
            mx=mx+1
 20         iar(mx)=nint(arrx(mg,lg)*con)
 30      write (6,40) lg,(iar(mg),mg=1,mx)
 40   format (1x,i2,2x,32i3)
      do 60 lg=lat,1,-1
         mx=0
         do 50 mg=mga,mgb
            mx=mx+1
 50         iar(mx)=nint(arrx(mg+lon,lg)*con)
 60      write (6,40) lg,(iar(mg),mg=1,mx)
         write (6,10) (mg,mg=mga,mgb)

      return
      end
