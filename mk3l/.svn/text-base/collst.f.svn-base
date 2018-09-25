c Fixing a bug within the subroutine COLLDYU, whereby the array ARR, rather
c than the array TX, was being passed to NC_PUT.
c SJP 2009/05/11
c
c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Modified for five-digit year numbers.
c SJP 2004/09/29
c
c Move definitions of LON9 and LAT9 out of IF blocks in order to
c avoid compiler warnings.
c SJP 2001/12/13
c
c $Log: collst.f,v $
c Revision 1.36  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.35  2001/02/28 04:36:39  rot032
c Further tidy ups from HBG
c
c Revision 1.34  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.33  2000/12/08 03:58:53  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.32  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.31  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.30  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.29  1999/06/16 06:21:54  rot032
c HBG changes to V5-3
c
c Revision 1.28  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.27  1997/06/11  02:21:31  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.26  1996/10/24  01:02:33  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.25  1996/06/13  02:05:51  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.24  1996/03/21  03:18:32  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.23  1995/11/23  06:03:27  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.22  1995/08/08  02:02:13  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.21  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.20  95/03/14  05:02:17  ldr
c Merge of INS changes for new ncput into LDR's latest version.
c 
c Revision 1.19  1995/03/14  04:58:52  ldr
c Changes for new version of ncput which handles reduced T63 arrays.
c
c Revision 1.18  94/08/08  17:20:55  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.17.1.1  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.17  94/03/30  12:34:05  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.16  93/12/17  15:31:54  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.15  93/11/29  14:40:23  ldr
c Merge of HBG and SPO changes.
c 
c Revision 1.14  93/10/19  11:14:42  ldr
c Add extra diagnostics for qcloud scheme.
c 
c Revision 1.13.1.1  93/11/29  12:33:10  ldr
c SPO's changes to rationalize printing of icuv statistics.
c 
c Revision 1.13  93/10/05  13:05:30  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.12  93/09/13  14:35:13  mrd
c Addded deep soil moisture percolation diagnostic
c 
c Revision 1.11  93/07/28  12:45:21  ldr
c Change varname for interception to inr so as not to confuse ncgen.
c 
c Revision 1.10  93/07/22  10:17:25  mrd
c Add diagnostics for clear sky net surface radiation and downward radiation
c at surface.
c
c     INPUT/OUTPUT
c     Input:   from common/gausl in GAUSL.f
c                  w - normal guass weights
c
c              from common/masiv2 in MASIV2.f
c                  radstm - global statistics array
c
c              from arguments
c                  con - factoring constant
c                  str - character string YYYYYMM.XXXXX where XXXXX is runtype
c                  kvar- variable number
c
c     In/Out:  from common/timex in TIMEX.f
c                  iyear - year counter for model run
c                  month - month counter
c
c 
      subroutine collst (kvar,con,str)

      implicit none
C Global parameters
      include 'PARAMS.f'
      integer numcdf
      parameter (numcdf=500) ! Must be >= total number of NetCDF files
      integer lon9
      parameter (lon9=lon/3)
      integer lat9
      parameter (lat9=lat2/3)

C Argument list
      integer kvar
      real con
      character*13 str

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'GAUSL.f'
      include 'MASIV2.f'   !Input stats
      include 'TIMEX.f'    !Input iyear, month

      integer ncdfcom(numcdf),nncdf
      character*3 ncdfchr(numcdf)
      common/ncdfiles/nncdf,ncdfcom,ncdfchr

C Local work arrays and variables
      real tx(lon,lat2)
      integer ind ! compression factor (1 or 9)
      character*13 filename
c---- These are the 9:1 compression arrays
      real txmid(64,lat2)
      real txcomp(64,32)    ! this is the minimal field for plotting.
      real wg(lat2)
c----

      integer kcc
      integer lg
      integer lid
      integer mg
      integer mga
      integer mid

      real sumw

C Local data, functions etc
      character*3 varname(14:max_radstm)
      data varname /'evp','rnd','hfl','cld','cll','clm','clh','rgn',
     &     'rtu','sot','sgn','tsu','tsc','wfg','wfb','snd','sid',
     &     'rtc','soc','tb2','tb3','vmo','tax','tay','rnc','run',
     &     'inr','als','tgg','tgf','pev','sev','rgc','sgc','rgd',
     &     'sgd','per','clc','dtm','lwp','pwc','rev','ssb','sno',
     &     'ref','cli','dms','so2','so4','s2d','s4d','s2w','s4w',
     &     'sem','s2o','s2h','s23','iwp','dm1','s21','s41','dmo',
     &     'dmn','scw'/

C Start code : ----------------------------------------------------------

c**** ROUTINE TO COLLECT PHYSICAL STATS FROM THE N1,S1,N2,S2...
c**** METHOD OF STORAGE IN /MASIV2/RADSTM INTO ARRAYS SUITABLE FOR
c**** PLOTTING,ETC : EG. (LON,LAT*2) WHERE 1,1 IS (SOUTH POLE,GM).
c----
c---- Note : Data is collected from arrays arranged (ln2,kvar,lat)
c----

c---- ind controls whether the full field is dumped (ind.eq.1), or
c----     if a compressed (area averaged) form is created (ind.eq.9).
c---- (ind=9 for T63 only). If T63, then
c---- this routine takes data from T63 in arrays (192,96) and may
c---- compress (average) to arrays (64,32). The use of "64" etc is
c---- NOT a hangover from the R21 model!!
c----
c        lat k+1  :  x - x - x  :              : X :
c        lat k    :  x - x - x  :  reduced to  : X :
c        lat k-1  :  x - x - x  :              : X :
c 
c----      by 3 point averaging (simple) in EW direction followed
c----      by a 3 point weighted average in NS direction onto lat k
c----      i.e. 9:1 reduction in data storage.

c Determine which statistic, and the associated compression (if any)
      ind=1
      do kcc=1,nncdf
        if(varname(kvar).eq.ncdfchr(kcc))ind=ncdfcom(kcc)
      enddo
      
      if(ind.eq.9)then
c---- check that T63 is in use
      if (mw.ne.64) then
         print *,' wrong resolution for collst 9:1 compression'
         stop
      end if
      endif
 
c Work out filename:
c Character string STR*13 contains YYYYYMM.XXXXX where XXXXX is runtype
      filename='s'//varname(kvar)//'_'//str(9:13)//'.nc' !Harvey's suggestion

c Transfer data to plotting array : (1,1) is SP,GM
c and Write data
       do 20 lg=1,lat
          do 20 mg=1,lon
          tx(mg,lg)=radstm(mg+lon,kvar,lg)*con
 20       tx(mg,lat2p-lg)=radstm(mg,kvar,lg)*con
      if (ind.eq.1) then
c---- write full field
        call nc_put(lon, lat2, filename, varname(kvar), iyear, month,
     $              tx)
      else
c---- form 9:1 compressed field and then write
         do 30 lg=1,lat2
            do 30 mg=1,lon-2,3
            mid=(mg/3)+1
            mga=mg-1
            if (mga.eq.0) mga=lon
 30         txmid(mid,lg)=(tx(mga,lg)+tx(mg,lg)+tx(mg+1,lg))/3.0
c---- create global G-weights
            do 38 lg=1,lat
            wg(lg)=w(lg)
 38         wg(lat2+1-lg)=w(lg)
         do 40 lg=2,lat2-1,3
            lid=(lg/3)+1
            sumw=wg(lg-1)+wg(lg)+wg(lg+1)
            do 40 mg=1,64
            txcomp(mg,lid)=(wg(lg-1)*txmid(mg,lg-1)+
     &       wg(lg)*txmid(mg,lg)+wg(lg+1)*txmid(mg,lg+1))/sumw
 40      continue
         call nc_put(lon9, lat9, filename, varname(kvar), iyear, month,
     $               txcomp)
      end if

      return
      end
C---------------------------------------------------------------------
      subroutine colldy (varname,con,str,arr)

      implicit none
C Global parameters
      include 'PARAMS.f'
      integer numcdf
      parameter (numcdf=500) ! Must be >= total number of NetCDF files
      integer lon9
      parameter (lon9=lon/3)
      integer lat9
      parameter (lat9=lat2/3)

C Argument list
      character*3 varname
      real con
      character*13 str
      real arr(ln2,lat)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'GAUSL.f'
      include 'TIMEX.f'    !Input iyear, month

      integer ncdfcom(numcdf),nncdf
      character*3 ncdfchr(numcdf)
      common/ncdfiles/nncdf,ncdfcom,ncdfchr

C Local work arrays and variables
      real tx(lon,lat2)
      integer ind ! compression factor (1 or 9)
      character*13 filename
c---- These are the 9:1 compression arrays
      real txmid(64,lat2)
      real txcomp(64,32)    ! this is the minimal field for plotting.
      real wg(lat2)
      real flag(lon,lat2)
c----

      integer kcc
      integer lg
      integer lid
      integer mg
      integer mga
      integer mid

      real sumf
      real sumwf
      real wf
      real wfm1
      real wfp1

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** ROUTINE TO OUTPUT DYNAMICAL STATS FROM THE N1,S1,N2,S2...
c**** METHOD OF STORAGE INTO ARRAYS SUITABLE FOR
c**** PLOTTING,ETC : EG. (LON,LAT*2) WHERE 1,1 IS (SOUTH POLE,GM).
c----
c---- Note : Data is collected from arrays arranged (ln2,lat)
c----

c Determine which statistic, and the associated compression (if any)
      ind=1
      do kcc=1,nncdf
        if(varname.eq.ncdfchr(kcc))ind=ncdfcom(kcc)
      enddo

      if(ind.eq.9)then
c---- check that T63 is in use
      if (mw.ne.64) then
         print *,' wrong resolution for collst 9:1 compression'
         stop
      end if
      endif

c Work out filename
      filename='s'//varname//'_'//str(9:13)//'.nc'

c Transfer data to plotting array : (1,1) is SP,GM
c and Write data
      do 20 lg=1,lat
         do 20 mg=1,lon
         tx(mg,lg)=arr(mg+lon,lg)*con
         if(arr(mg+lon,lg).eq.-7777777.0)tx(mg,lg)=-7777777.0
         tx(mg,lat2p-lg)=arr(mg,lg)*con
         if(arr(mg,lg).eq.-7777777.0)tx(mg,lat2p-lg)=-7777777.0
 20   continue
      if (ind.eq.1) then
c---- write full field
        call nc_put(lon, lat2, filename, varname, iyear, month, tx)
      else
c---- write 9:1 compressed field
c---- Note - The dynamical fields may contain dummy data inside
c----        mountains flagged by -7777777
         do 30 lg=1,lat2
            do 30 mg=1,lon
            flag(mg,lg)=1.0
            if (tx(mg,lg).eq.-7777777.0) flag(mg,lg)=0.0
 30      continue
         do 40 lg=1,lat2
            do 40 mg=1,lon-2,3
            mid=(mg/3)+1
            mga=mg-1
            if (mga.eq.0) mga=lon
            sumf=flag(mga,lg)+flag(mg,lg)+flag(mg+1,lg)
            txmid(mid,lg)=-7777777.0
            if (sumf.gt.0.5) txmid(mid,lg)=
     &       (tx(mga,lg)*flag(mga,lg)+tx(mg,lg)*flag(mg,lg)+
     &        tx(mg+1,lg)*flag(mg+1,lg))/sumf
 40      continue
         do 50 lg=1,lat2
            do 50 mg=1,64
            flag(mg,lg)=1.0
 50         if (txmid(mg,lg).eq.-7777777.0) flag(mg,lg)=0.0
c---- create global G-weights
            do 58 lg=1,lat
            wg(lg)=w(lg)
 58         wg(lat2+1-lg)=w(lg)
         do 60 lg=2,lat2-1,3
            lid=(lg/3)+1
            do 60 mg=1,64
            wfm1=wg(lg-1)*flag(mg,lg-1)
            wf  =wg(lg)  *flag(mg,lg)
            wfp1=wg(lg+1)*flag(mg,lg+1)
            sumwf=wfm1+wf+wfp1
            txcomp(mg,lid)=-7777777.0
            if (sumwf.gt.0.0) txcomp(mg,lid)=(txmid(mg,lg-1)*wfm1+
     &       txmid(mg,lg)*wf+txmid(mg,lg+1)*wfp1)/sumwf
 60      continue
        call nc_put(lon9, lat9, filename, varname, iyear, month,
     $              txcomp)
      end if

      return
      end
C---------------------------------------------------------------------
      subroutine collph (varname,con,str,arr)

C Global parameters
      include 'PARAMS.f'
      integer numcdf
      parameter (numcdf=500) ! Must be >= total number of NetCDF files
      integer lon9
      parameter (lon9=lon/3)
      integer lat9
      parameter (lat9=lat2/3)

C Argument list
      character*3 varname
      real con
      character*13 str
      real arr(ln2,lat)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'GAUSL.f'
      include 'TIMEX.f'    !Input iyear, month

      integer ncdfcom(numcdf),nncdf
      character*3 ncdfchr(numcdf)
      common/ncdfiles/nncdf,ncdfcom,ncdfchr

C Local work arrays and variables
      real tx(lon,lat2)
      character*13 filename
c---- These are the 9:1 compression arrays
      real txmid(64,lat2)
      real txcomp(64,32)    ! this is the minimal field for plotting.
      real wg(lat2)
c---- 

      integer kcc
      integer lg
      integer lid
      integer mg
      integer mga
      integer mid

      real sumw

C Local data, functions etc

C Start code : ----------------------------------------------------------
c**** ROUTINE TO OUTPUT PHYSICAL STATS FROM THE N1,S1,N2,S2...
c**** METHOD OF STORAGE IN /MAXMIN/ INTO ARRAYS SUITABLE FOR
c**** PLOTTING,ETC : EG. (LON,LAT*2) WHERE 1,1 IS (SOUTH POLE,GM).
c----
c---- Note : Data is collected from arrays arranged (ln2,lat)
c----

c Determine which statistic, and the associated compression (if any)
      ind=1
      do kcc=1,nncdf
        if(varname.eq.ncdfchr(kcc))ind=ncdfcom(kcc)
      enddo

      if(ind.eq.9)then
c---- check that T63 is in use
      if (mw.ne.64) then
         print *,' wrong resolution for collst 9:1 compression'
         stop
      end if
      endif

c Work out filename
      filename='s'//varname//'_'//str(9:13)//'.nc'

c Transfer data to plotting array : (1,1) is SP,GM
c and Write data
      do 20 lg=1,lat
         do 20 mg=1,lon
         tx(mg,lg)=arr(mg+lon,lg)*con
 20      tx(mg,lat2p-lg)=arr(mg,lg)*con
      if (ind.eq.1) then
c---- write full field
        call nc_put(lon, lat2, filename, varname, iyear, month, tx)
      else
c---- write 9:1 compressed field
         do 30 lg=1,lat2
            do 30 mg=1,lon-2,3
            mid=(mg/3)+1
            mga=mg-1
            if (mga.eq.0) mga=lon
            txmid(mid,lg)=(tx(mga,lg)+tx(mg,lg)+tx(mg+1,lg))/3.0
 30      continue
c---- create global G-weights
            do 38 lg=1,lat
            wg(lg)=w(lg)
 38         wg(lat2+1-lg)=w(lg)
         do 40 lg=2,lat2-1,3
            lid=(lg/3)+1
            sumw=wg(lg-1)+wg(lg)+wg(lg+1)
            do 40 mg=1,64
            txcomp(mg,lid)=(wg(lg-1)*txmid(mg,lg-1)+
     &       wg(lg)*txmid(mg,lg)+wg(lg+1)*txmid(mg,lg+1))/sumw
 40      continue
        call nc_put(lon9, lat9, filename, varname, iyear, month,
     $              txcomp)
      end if

      return
      end
C---------------------------------------------------------------------
      subroutine colldyu (varname,con,str,arr)

c---- Collect statistics on (U,V) grid for ice model  ----------

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      character*3 varname
      real con
      character*13 str
      real arr(lon,lat2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
      include 'TIMEX.f'    !Input iyear, month

C Local work arrays and variables
      real tx(lon,lat2)
      character*13 filename
      integer lg
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** ROUTINE TO OUTPUT ICE MODEL STATS FROM THE N1,S1,N2,S2...
c**** METHOD OF STORAGE INTO ARRAYS SUITABLE FOR
c**** PLOTTING,ETC : EG. (LON,LAT*2) WHERE 1,1 IS (SOUTH POLE,GM).
c----
c---- Note : Data is collected from arrays arranged (lon,lat2)
c----

c Work out filename
      filename='s'//varname//'_'//str(9:13)//'.nc'

c and Write data
c Transfer data to plotting array : (1,1) is SP,GM
c---- write full field (No compression for ice model data)
         do 20 lg=1,lat2
            do 20 mg=1,lon
 20         tx(mg,lg)=arr(mg,lg)*con

      call nc_put(lon, lat2, filename, varname, iyear, month, tx)

      return
      end
