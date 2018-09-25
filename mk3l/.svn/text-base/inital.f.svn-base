c Removing the declaration of the COMMON block /GRADNS/ to a header file.
c SJP 2009/04/25
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Transferred COMMON blocks to separate header files, as follows:
c /MDAY/  ->  MDAY.f
c SJP 2007/05/29
c
c Modified to enable five-digit year numbers.
c SJP 2004/09/22
c
c $Log: inital.f,v $
c Revision 1.61  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.60  2001/02/12 05:39:45  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.59  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.58  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.57  1998/12/10  00:55:32  ldr
c HBG changes to V5-1-21
c
c Revision 1.56  1998/05/27  02:07:36  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.55  1997/12/23  00:23:34  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.54  1997/12/19  01:25:33  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.53.3.1  1997/12/19  02:03:11  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.53  1996/10/24  01:27:40  ldr
c Merge of TIE changes with LDR/MRD changes since V5-0-12.
c
c Revision 1.52  1996/08/12  01:51:22  mrd
c Generalise for triangular truncations other than T63
c
c Revision 1.51.1.1  1996/10/24  01:02:57  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.51  1996/06/13  02:06:56  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.50  1996/04/04  02:16:16  mrd
c Move code to get a cleaner split between spectral initialisation in gauleg
c and other initialisation in inital. gauleg can now be called before inital.
c
c Revision 1.49  1996/03/21  03:18:52  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.48  1996/02/19  23:23:07  ldr
c Do vertical SLT properly for 24 level version.
c
c Revision 1.47  1995/11/23  05:52:11  ldr
c Fix for T63 header and new nsstop=-2 option.
c
c Revision 1.46  1995/08/29  01:30:09  ldr
c HBG's corrections to his V4-7-13h, to make results with hybrid=F agree
c with previous version.
c
c Revision 1.44  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.43  94/08/08  17:21:36  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.42.1.1  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.42  93/12/17  15:32:57  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.41  93/10/05  13:06:31  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.40  93/08/31  11:14:01  ldr
c Changes from HBG to fix problems with RS files in multi-month runs.
c 
c Revision 1.39  93/08/19  15:08:23  ldr
c Minor cosmetic changes.
c 
c Revision 1.38  93/08/19  11:50:07  ldr
c Replace INCD with NDSTOP in namelist and add required housekeeping code.
c 
c Revision 1.37  93/08/17  11:53:54  ldr
c Changed MSTEP=0 option for R42 and added print of ratlm.
c 
c Revision 1.36  93/08/10  16:13:48  ldr
c Merge of HBG changes to V4-3-28l with other changes since then.
c 
c Revision 1.35  93/08/03  11:35:49  ldr
c  Replace intrinsic functions with generic forms where possible.
c 
c Revision 1.34  93/07/29  12:46:10  ldr
c Corrected calculation of ratlm from previous revision.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  ifwds  - default F, only T when doing forward step from
c                           new restart file
c                  ltrace - if T, calls Langrangian tracer routine, tracera
c                  ndstop - if non-zero, stop model run after ndstop days
c                  nsstop - if > zero, stop model run after nsstop steps
c                  rainflag  - if T write daily rainfall to file
c                  saveqflux - if T,save fields required for running model
c                              in its qflux form to the file qflxyymm.xxx
c
c              from common/files in FILES.f
c                  irfilename - name of input restart file
c
c              from arguments
c                  mms - counter of model run months, not calendar months
c
c     In/Out:  from common/cnsta in CNSTA.f
c                  sig  - sigma values    sigh - sigma values for half levels
c                  algf - log(sig(k))     algh - log(sigh(k))
c                  dsk  - sigma thicknesses (1 to nl)
c
c              from common/si in SI.f
c                  am  - geopotential height matrix
c                  amh -      "         "      "   , half levels
c
c              from common/timex in TIMEX.f
c                  iyear - year counter for model run
c                  kday  - first day of run   kdayp - last day of run
c                  ldays - month count of days   month - month counter
c                  mstep - timestep in minutes 
c                  ndays - day counter for model run
c                  kdays - year count of days, runs 1 to 365
c                  nstep - no. of time steps for this run
c                  ratlm - fraction of month elapsed, runs 0-1
c                  nrad - no timesteps between calls to the radiation scheme,
c                         set to zero for no radiation, default 4
c
c              from common/worka in WORKA.f
c                  phisb - global mean geopotential height
c                  phnb  - approx global mean model level geopotential height
c                  tmean - preset isothermal temperature
c
c     Output:  from common/cnsta in CNSTA.f
c                  dshdk - d(sigh)/dk at bottom of layer k
c                  sigk - sigma**cappa      sighk - sigh(k)**cappa
c
c              from common/cnste in CNSTE.f
c                  toph - top half of diffusion indicator matrix
c
c              from common/gradns in GRADNS.f
c                  sigi - 1.0-sig(k)
c
c              from common/timex in TIMEX.f
c                  dt   - time step in seconds
c                  mins - current model time in mins
c
c              from common/worka in WORKA.f
c                  emean  - preset global mean energy const/model level
c                  tmeanh - preset isothermal temperature at half levels
c
c
      subroutine inital(mms)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer mms

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'CNSTE.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'RKLF.f'
      include 'SI.f'
      include 'TIMEX.f'
      include 'WORKA.f'
      include 'MDAY.f'
      include 'GRADNS.f'

C Local work arrays and variables
      real atm(nlm),atp(nlm),btm(nlm),btp(nlm)
      character*20 filename

      character*50 exptyp
      save exptyp

      integer i
      integer ifiler
      integer ierr
      integer ios
      integer j
      integer k
      integer km
      integer kp
      integer ksum
      integer ll
      integer lls
      integer minl
      integer mm
      integer mstepx
      integer nday
      integer ndayp

      real alf
      real den
      real rk
      real sigt

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** ROUTINE TO SET UP INITIAL VALUES

c Set up sigma half levels

      if(nl.gt.18)then !Use general cubic Smagorinsky formula
        sigt=0.004458  !Select top sigma level
        den=.75*( 1./nl-.5 - .5*(2./nl-1.)**3 )
        alf=(sigt-.25 -.25*(2./nl-1.)**3 )/den
      else   !Use old simpler formula
        alf=1.
      endif

      do k=1,nl+1
        rk=(.5*nl-(k-1))/nl
        sigh(k)=.5 + 1.5*alf*rk - 2.*(3.*alf-2.)*rk**3
c Calculate dsigh/dk at bottom of layer k (for SLT scheme)
        dshdk(k)=-1.5*(alf/nl - (3.*alf-2)*(nl-2*k+2.)**2/nl**3)
      enddo
      sigh(1)=1.0
      sigh(nl+1)=0.0

      do k=1,nl
         kp=k+1
         sig(k)=0.5*(sigh(k)+sigh(kp))
         sigk(k)=sig(k)**cappa
         sighk(k)=sigh(k)**cappa
         algf(k)=log(sig(k))
         algh(k)=log(sigh(k))
         dsk(k)=sigh(k)-sigh(kp)
      end do
      sighk(nlp)=0.0
      do k=1,nl
         sigi(k)=1.0-sig(k)
      end do
      do 9900 k=1,nl-1
      btm(k)=1.0/(algf(k)-algf(k+1))
      btp(k)=-btm(k)
      atm(k)=-algf(k+1)*btm(k)
 9900 atp(k)=algf(k)*btm(k)
      do 331 k=1,nl
  331 write(6,1612)k,sig(k),sigh(k)
 1612 format(1x,'k=',i2,' sig=',f7.4,' sigh=',f7.4)

c**** Using (sigma) OR (hybrid) vertical coordinate 
      call vertc


C**   CREATE DIFFUSION MULTIPLIER ARRAY FOR UPPER SPECTRUM
C**   (HAS TO BE SET FOR EACH TRUNCATION)
C**   AT T63, THE UPPER PART IS DIFFUSED.
C**   BUT IF jsfdiff=.true., THEN USE JSF'S SPECTRAL DIFFUSION.

      do 1513 mm=1,mw
      do 1513 ll=1,lw
      if(jsfdiff)then
        toph(ll,mm)=1.0
      else
        toph(ll,mm)=0.0
      endif
 1513 continue

      if(jsfdiff)then

       filename = 'rklfred.t63.nec'
       open (99,file=filename,status='old',form='unformatted'
     &  ,iostat=ios)
       if (ios .ne. 0) 
     &  write (6,*) 'open rklfred failed, ios = ',ios
       read (99) rklfred
       close (99)

      else

       minl=lw
c      if ( trunc.eq.'T' ) minl=lw-mw/16
       if ( trunc.eq.'T' ) minl=lw-mw/8
       do 152 mm=1,mw
       lls=max(1,minl-mm+1)
       if(lls.gt.lw)go to 152
       do 153 ll=lls,lw
  153  toph(ll,mm)=1.0
  152  continue

      endif

c****  temp/geopotl coupling for phi=am*t+phist
c***  based upon t=a+b.log(sigma) : energy conserving -
c****  changes made to dynamical terms (in DYNMNL).
      do 532 j=1,nl
      do 531 i=1,nl
      am(i,j)=0.0
      amh(i,j)=0.0
  531 um(i,j)=0.0
  532 um(j,j)=1.0
c**** define temporary full level geopotential matrix values
      am(1,1)=rdry*(atm(1)*(algh(1)-algf(1))+btm(1)*(algh(1)**2
     & -algf(1)**2)*0.5)
      am(1,2)=rdry*(atp(1)*(algh(1)-algf(1))+btp(1)*(algh(1)**2
     & -algf(1)**2)*0.5)
      am(2,1)=rdry*(atm(1)*(algh(1)-algf(2))+btm(1)*(algh(1)**2
     & -algf(2)**2)*0.5)
      am(2,2)=rdry*(atp(1)*(algh(1)-algf(2))+btp(1)*(algh(1)**2
     & -algf(2)**2)*0.5)
      do 9925 i=3,nl
      do 9925 j=1,i
      if(j.eq.i-1)go to 9915
      if(j.eq.i)go to 9920
      am(i,j) = am(i-1,j)
      go to 9925
 9915 am(i,j) = am(i-1,j)
     &+rdry*(atm(i-1)*(algf(i-1)-algf(i))+btm(i-1)*(algf(i-1)**2
     & -algf(i)**2)*0.5)
      go to 9925
 9920 am(i,j) = am(i-1,j)
     &+rdry*(atp(i-1)*(algf(i-1)-algf(i))+btp(i-1)*(algf(i-1)**2
     & -algf(i)**2)*0.5)
 9925 continue
c**** create half level geopotential matrix values
      amh(1,1)=rdry*(atm(1)*(algh(1)-algh(2))+btm(1)*(algh(1)**2
     & -algh(2)**2)*0.5)
      amh(1,2)=rdry*(atp(1)*(algh(1)-algh(2))+btp(1)*(algh(1)**2
     & -algh(2)**2)*0.5)
      do 9950 i=2,nl-1
      do 9945 j=1,nl
 9945 amh(i,j)=am(i,j)
c     j=i
      amh(i,i) = amh(i,i)
     &+rdry*(atm(i)*(algf(i)-algh(i+1))+btm(i)*(algf(i)**2
     & -algh(i+1)**2)*0.5)
c     j=i+1
      amh(i,i+1) = amh(i,i+1)
     &+rdry*(atp(i)*(algf(i)-algh(i+1))+btp(i)*(algf(i)**2
     & -algh(i+1)**2)*0.5)
 9950 continue
      do 9955 j=1,nl
 9955 amh(nl,j)=2.0*am(nl,j)-amh(nl-1,j)
c**** create final full level geopotential matrix values
      do 9975 j=1,nl
 9975 am(1,j)=0.5*amh(1,j)
      do 9980 i=2,nl
      do 9980 j=1,nl
 9980 am(i,j)=0.5*(amh(i-1,j)+amh(i,j))

c     write(6,4475)
c4475 format(1x,' t=a+b.log(sigma) energy cons am(,)')
c     do 4476 j=1,nl
c4476 write(6,4477)(am(i,j),i=1,nl)
c4477 format(1x,9f7.1)

C****
C**** SET ISOTHERMAL MEAN TEMPS
C****
      do 123 k=1,nl
  123 tmean(k)=290.0
      write(6,124)tmean(1)
  124 format(1h ,'isothermal mean temp=',f6.1)
      tmeanh(1)=0.0
      do 636 k=2,nl
      km=k-1
  636 tmeanh(k)=(dsk(km)*tmean(k)+dsk(k)*tmean(km))/(dsk(km)+dsk(k))
      tmeanh(nlp)=0.0
c  Value for the new, rougher topography.
      phisb=2281.0
C**** set mean k.e. term at each level
C****   note - these are not final - temporary values only
      do 637 k=1,nl
  637 emean(k)=(nl/9.0)*310.0/eradsq
C**** set mean geopotential term at each level
      do 122 k=1,nl
      phnb(k)=phisb
      do 122 j=1,nl
  122 phnb(k)=phnb(k)+am(k,j)*tmean(j)

c----  OPEN THE MAIN ATMOSPHERIC RESTART FILE TO GET DAYNUMBER
c----  Only do this at first month of a multiple month run
      If (mms.eq.1) Then

      ifiler=8  
      open(unit=8,file=irfilename,form='unformatted',status='old'
     &     ,iostat=ierr)
      call errcheck(ierr,'restart file     ','inital    ')
c----  READ IN HEADER FROM RESTART FILE : EXPTYP IS CHARACTER*50
c----  AND GET INITIAL DAY FROM RESTART FILE
      read(ifiler)exptyp
      read(ifiler)ndays,mstepx
      close(unit=ifiler)

      if(ifwds.or.nsstop.eq.-2)then
        ndays=0
        write(6,*)' Setting NDAYS to 0 in inital (New model)'
      endif

      End IF

      kday=ndays
C**   CHECK FOR END OF MONTH AFTER INCD DAYS
      if(ndstop.eq.0)then
        incd=31
      else
        incd=ndstop
      endif
      nday=kday-(kday/365)*365
      month=0
   21 month=month+1
      nday=nday-mdays(month)
      if(nday.ge.0)go to 21
      nday=nday+mdays(month)
      ndayp=nday+incd
      if(ndayp.gt.mdays(month))ndayp=mdays(month)
      incd=ndayp-nday
      kdayp=kday+incd
      incd=kdayp-kday
      if(incd.gt.0)go to 72
      write(6,630)kday
  630 format(1h ,'model stopped at day ',i5)
      stop
   72 write(6,410)kday,kdayp,incd
  410 format (1h ,'start from day',i9,' up to day',i9,' incd=',i2)
      if(mstep.eq.0)then  ! set mstep automatically
        if(lw.eq.22)then
          mstep=30
        elseif((lw.eq.43).or.(lw.eq.64))then
          mstep=15
        endif
        nrad=120/mstep          !Allow 2 hours between radiation calls
        print*,'MSTEP has been automatically set in inital to ',mstep
        print*,'NRAD has been correspondingly set to ',nrad        
      endif
      dt=mstep*60.0
      nstep=incd*(24*60/mstep)
      write(6,411)nstep
  411 format (1h ,'number of atmos timesteps to be carried out=',i5)
      mins = int(kday, 8) * int(1440, 8)
C**   KDAYS : YEAR COUNT OF DAYS , RUNS 1 TO 365 , UPDATED AT START OF
C**   EACH DAY.
C**   LDAYS : MONTH COUNT OF DAYS , RUNS 0 TO MDAYS(MONTH)-1 , I.E.
C**   IT COUNTS THE NUMBER OD DAYS INTEGRATED SO FAR IN A GIVEN MONTH.
C**   REQUIRED FOR T*,A* INTERPOLATION.
C**   MONTH : GIVES CURRENT MONTH NUMBER
      ndays=kday
      kdays=ndays+1
      iyear=kdays/365
      kdays=kdays-iyear*365
C**   DETERMINE MONTH
      ksum=0
      month=0
   20 month=month+1
      ksum=ksum+mdays(month)
      if(kdays.gt.ksum)go to 20
      ldays=kdays-ksum+mdays(month)-1
      write(6,412)iyear,month,ldays+1
  412 format(1h ,'model restarted at year ',i5.5,
     & ' month ',i2,' day ',i2)
CJAN  MONTH=1
CJAN  LDAYS=0
      ratlm=float(ldays)/mdays(month)
      write(6,6676)kdays,ratlm
 6676 format(1h ,' inital : kdays=',i5,' ratlm=',f7.4)

c----
c---- OPEN DATA STORAGE FILES
c----
      call openfl(exptyp)

      return
      end
C---------------------------------------------------------------------
      subroutine vertc

c***********************************************************************
c**** To compute HYBRID vertical coordinate A(n), B(n) and dA/dn, dB/dn
c**** values for CSIRO model. A(n),B(n) - half
c**** levels and dA/dn,dB/dn,B(n) - full levels. 
c****
c**** The pressure at any level (n) is defined by
c****
c****  p(n,Pst) = Pooeta*A(n) + Pst*B(n); Pst=surface pressure
c****                                     Pooeta=1013.2
c**** 
c**** This also applies for SIGMA vertical coordinates if A(n)=0,
c**** B(n)=sigma, dA/dn=0 and dB/dn=1
c**********************************************************************

      implicit none
C Global parameters
      include 'PARAMS.f'
c**** neq = number of conditions/matrix equations = 10
      integer neq,neq2
      parameter (neq=10,neq2=neq/2)

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f' 
      include 'FEWFLAGS.f'
      include 'HYBARR.f'

C Local work arrays and variables
      real c(neq,1)
      real uk2(neq,neq)

      integer i
      integer j
      integer k

      real alpha
      real apb
      real ex
      real pooeta
      real ptm2
      real p05
      real p500
      real suma
      real sumb

      logical oldhybrid

C Local data, functions etc

C Start code : ----------------------------------------------------------

      IF (hybrid) THEN

c****=============== FOR HYBRID COORDINATES ===============**** 

c**** Create the elements of matrics (equations) and solve
      Pooeta=1013.2
c**** Top two half levels to be constant pressure
      ptm2=Pooeta*sigh(nlp-2)
c**** Lowest two levels thickness at 1013.2mb = alpha*(lowest two levels
c**** thickness at 500mb)
      alpha=(500.0-ptm2)/(Pooeta-ptm2)
      alpha=alpha*0.5    
      p500=alpha*Pooeta-500.0
      p05=(1.0-alpha)*Pooeta

      do 10 j=1,neq
      c(j,1)=0.0
      do 10 i=1,neq
   10 uk2(i,j)=0.0

c**** condition 1
      do 12 j=1,neq2
   12 uk2(1,j)=1.0
c**** condition 2
      do 14 j=1,neq2
   14 uk2(2,j+neq2)=1.0
      c(2,1)=1.0
c**** condition 3
      do 16 j=1,neq2
   16 uk2(3,j+neq2)=sigh(nlp-2)**j
c**** condition 4
      do 18 j=1,neq2
   18 uk2(4,j+neq2)=sigh(nlp-1)**j
c**** condition 5
      do 22 j=1,neq2
   22 uk2(5,j)=sigh(nlp-2)**j
      c(5,1)=sigh(nlp-2)
c**** condition 6
      do 24 j=1,neq2
   24 uk2(6,j) =sigh(nlp-1)**j
      c(6,1)=sigh(nlp-1)
c**** condition 7
      do 26 j=1,neq2
      uk2(7,j)=p05*(sigh(1)**j-sigh(2)**j)
   26 uk2(7,j+neq2)=-p500*(sigh(1)**j-sigh(2)**j)
c**** condition 8
      do 28 j=1,neq2
      uk2(8,j)=p05*(sigh(2)**j-sigh(3)**j)
   28 uk2(8,j+neq2)=-p500*(sigh(2)**j-sigh(3)**j)
c**** conditions 9 and 10 (ok for nl>=7)
      if(nl.lt.7)then
        print *,'hybrid 9 & 10 not ok for nl<7'
        stop
      endif
      i=(4*nl)/9
      do 32 j=1,neq2
      uk2(9,j)=sigh(i)**j
   32 uk2(9,j+neq2)=sigh(i)**j
      c(9,1)=sigh(i)
      i=(6*nl)/9
      do 34 j=1,neq2
      uk2(10,j)=sigh(i)**j
   34 uk2(10,j+neq2)=sigh(i)**j
      c(10,1)=sigh(i)

c**** FIND THE MATRIX SOLUTION: uk2(i,j) * xx(i,1) = c(i,1) ****
       call gaussj(neq,neq,1,1,uk2,c)
c**** MATRIX c(i,1) HAS BEEN REPLACED BY THE SOLUTION xx(i,1) 

c      oldhybrid=.true. ! This can be made false - see below
      oldhybrid=.false. ! For Mk3.1 etc

      If (oldhybrid) Then 
c----
c---- Default (original) hybrid level calculation (Mk2 & Mk3)
c----

c**** A=anh(k) and B=bnh(k) at each half level 
      do 582 k=1,nlp
      suma=0.0
      sumb=0.0
      do 584 j=neq2,1,-1
      ex=0.0
      if(sigh(k).ne.0.0)ex=sigh(k)**j
      suma=suma+c(j,1)*ex
  584 sumb=sumb+c(j+neq2,1)*ex
      anh(k)=suma
      bnh(k)=sumb
      apb=suma+sumb  
  582 continue
c**** Set A (half level 1) = B (half level 8,9) = 0.0 (exactly) as 
c**** defined in conditions 1,3 and 4 (see notes)
      anh(1)=0.0
      bnh(nlp-1)=0.0
      bnh(nlp-2)=0.0

c**** For consisentcy, and energy conservation, the values
c**** at FULL levels must be computed from HALF level values
      do 782 k=1,nl
      anf(k)=0.5*(anh(k)+anh(k+1))
      bnf(k)=0.5*(bnh(k)+bnh(k+1))
      dadnf(k)=(anh(k)-anh(k+1))/(sigh(k)-sigh(k+1))
      dbdnf(k)=(bnh(k)-bnh(k+1))/(sigh(k)-sigh(k+1))
  782 continue

      Else
c----
c---- More precise calculations for hybrid levels
c----  taking into account machine round off errors (Mk3++)
c----

c**** A=anh(k) and B=bnh(k) at each half level 
      do 5821 k=1,nlp
      suma=0.0
      sumb=0.0
CDIR$ NEXTSCALAR
*PDIR NOVECTOR
      do 5841 j=neq2,1,-1
      ex=0.0
      if(sigh(k).ne.0.0)ex=sigh(k)**j
      suma=suma+c(j,1)*ex
 5841 sumb=sumb+c(j+neq2,1)*ex
      anh(k)=suma
      bnh(k)=sumb
 5821 continue
*PDIR VECTOR

C**** Certain conditions should be met.
C**** These are :
c****  (1) A (half level 1) = B (half level nl & nl-1) = 0.0 (exactly)
c****      as defined in conditions 1,3 and 4 (see notes)
c****  (2) A (half level) + B (half level) = sigh (half level)
c****  (3) A (full level) + B (full level) = sig (full level)
c****  (4) dAdNf (full level) + dBdNf (full level) = 1.0
c****     with dAdNf(nl)=dAdNf(nl-1)=1.0 & dBdNf(nl)=dBdNf(nl-1)=0.0
C****
c**** However, due to the calculations involved, and round off,
c****  these criteria do not follow identically. 
c**** So adjustments will be undertaken to guarantee the above.
c**** These adjustments are minor (last decimal places), and
c****  will be applied first to both numbers equally, and if
c****  that is not sufficient, then the smallest number will
c****  be adjusted to guarantee the above conditions. See below.

c****  (1) A (half level 1) = B (half level nl & nl-1) = 0.0 (exactly)
      anh(1)=0.0
      bnh(nlp-1)=0.0
      bnh(nlp-2)=0.0
c****  (2) A (half level) + B (half level) = sigh (half level)
c****   Adjust (if needed) so that anh+bnh = sigh
      anh(nlp-1)=sigh(nlp-1)
      anh(nlp-2)=sigh(nlp-2)
      do 7801 k=2,nlp-3
        apb=(anh(k)+bnh(k))-sigh(k)
        anh(k)=anh(k)-apb/2
        bnh(k)=bnh(k)-apb/2
        apb=anh(k)+bnh(k)
        if(apb.ne.sigh(k))then
          if(anh(k).gt.bnh(k))then
            bnh(k)=sigh(k)-anh(k)
          else
            anh(k)=sigh(k)-bnh(k)
          endif
        endif
 7801 continue

c****  (3) A (full level) + B (full level) = sig (full level)
c**** For consisentcy, and energy conservation, the values
c**** at FULL levels must be computed from HALF level values
      do 7821 k=1,nl
        anf(k)=0.5*(anh(k)+anh(k+1))
        bnf(k)=0.5*(bnh(k)+bnh(k+1))
        apb=(anf(k)+bnf(k))-sig(k)
        anf(k)=anf(k)-apb/2
        bnf(k)=bnf(k)-apb/2
        apb=anf(k)+bnf(k)
        if(apb.ne.sig(k))then
          if(anf(k).gt.bnf(k))then
            bnf(k)=sig(k)-anf(k)
          else
            anf(k)=sig(k)-bnf(k)
          endif
        endif
 7821 continue

c****  (4) dAdNf (full level) + dBdNf (full level) = 1.0
c****     with dAdNf(nl)=dAdNf(nl-1)=1.0 & dBdNf(nl)=dBdNf(nl-1)=0.0
      do 7841 k=1,nl-2
        dadnf(k)=(anh(k)-anh(k+1))/(sigh(k)-sigh(k+1))
        dbdnf(k)=(bnh(k)-bnh(k+1))/(sigh(k)-sigh(k+1))
        apb=(dadnf(k)+dbdnf(k))-1.0
        dadnf(k)=dadnf(k)-apb/2
        dbdnf(k)=dbdnf(k)-apb/2
        apb=dadnf(k)+dbdnf(k)
        if(apb.ne.1.0)then
          if(dadnf(k).gt.dbdnf(k))then
            dbdnf(k)=1.0-dadnf(k)
          else
            dadnf(k)=1.0-dbdnf(k)
          endif
        endif
 7841 continue
      dadnf(nl)=1.0
      dbdnf(nl)=0.0
      dadnf(nl-1)=1.0
      dbdnf(nl-1)=0.0

      EndIf ! (oldhybrid)

c      write(88,*)"HYBRID COORDINATES"
c      do k=1,nl
c      write(88,*)(dadnf(k)+dbdnf(k))
c      enddo

c---- For use in model, multiply anf,anh,dadnf by factor Pooeta
      do k=1,nl
        anf(k)=anf(k)*Pooeta
        dadnf(k)=dadnf(k)*Pooeta
      enddo
      do k=1,nlp
        anh(k)=anh(k)*Pooeta
      enddo

      ELSE          ! end of Hybrid coordinate section

c****=============== FOR SIGMA COORDINATES ===============**** 
c**** A=anh(k), B=bnh(k), at each half level 
      do k=1,nlp
       anh(k)=0.0
       bnh(k)=sigh(k)
      enddo  
c**** A=anf(k), B=bnf(k),dadnf(k),dbdnf(k) at each full level 
      do k=1,nl
       anf(k)=0.0
       bnf(k)=sig(k)
       dadnf(k)=0.0
       dbdnf(k)=1.0
      enddo  

c      write(88,*)"SIGMA COORDINATES"
c      do k=1,nl
c      write(88,*)(dadnf(k)+dbdnf(k))
c      enddo

      ENDIF          ! end of Sigma coordinate section
 
      return
      end
C---------------------------------------------------------------------
      subroutine gaussj(n,np,m,mp,a,b)

      implicit none
C Global parameters
      integer nmax
      parameter (nmax=50)

C Argument list
      integer n
      integer np
      integer m
      integer mp
      real a(np,np)
      real b(np,mp)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer indxc(nmax),indxr(nmax),ipiv(nmax)

      integer i
      integer icol
      integer irow
      integer j
      integer k
      integer l
      integer ll

      real big
      real dum
      real pivinv

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do 11 j=1,n
        ipiv(j)=0
   11 continue

      do 22 i=1,n

      big=0.
      do 13 j=1,n
      if(ipiv(j).ne.1)then
      do 12 k=1,n
      if (ipiv(k).eq.0) then
        if (abs(a(j,k)).ge.big)then
          big=abs(a(j,k))
          irow=j
          icol=k
        endif
      else if (ipiv(k).gt.1) then
        print *,'singular matrix in gaussj'
        stop
      endif
   12 continue
      endif
   13 continue

      ipiv(icol)=ipiv(icol)+1

      if (irow.ne.icol) then

      do 14 l=1,n
        dum=a(irow,l)
        a(irow,l)=a(icol,l)
        a(icol,l)=dum
   14 continue
      do 15 l=1,m
        dum=b(irow,l)
        b(irow,l)=b(icol,l)
        b(icol,l)=dum
   15 continue

      endif

      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.) then
        print *,'singular matrix in gaussj'
        stop
      endif
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.

      do 16 l=1,n
        a(icol,l)=a(icol,l)*pivinv
   16 continue
      do 17 l=1,m
        b(icol,l)=b(icol,l)*pivinv
   17 continue

      do 21 ll=1,n

      if(ll.ne.icol)then

      dum=a(ll,icol)
      a(ll,icol)=0.
      do 18 l=1,n
        a(ll,l)=a(ll,l)-a(icol,l)*dum
   18 continue
      do 19 l=1,m
        b(ll,l)=b(ll,l)-b(icol,l)*dum
   19 continue

      endif

   21 continue

   22 continue

      do 24 l=n,1,-1

      if(indxr(l).ne.indxc(l))then

      do 23 k=1,n
        dum=a(k,indxr(l))
        a(k,indxr(l))=a(k,indxc(l))
        a(k,indxc(l))=dum
   23 continue

      endif

   24 continue

      return
      end
