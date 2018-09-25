c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c nsib comparison is modified to lsm_type = "nsib "
c AJA 2009/01/22
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/08
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Changes to add machine type 'ALPH', which makes use of the super-fast FFTW
c FFT library.
c SJP 2001/11/22
c
c $Log: phys.f,v $
c Revision 1.58  2001/02/28 04:36:37  rot032
c Further tidy ups from HBG
c
c Revision 1.57  2001/02/22 05:34:38  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.56  2001/02/12 05:39:47  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.55  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.54  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.53  1998/12/10  00:55:35  ldr
c HBG changes to V5-1-21
c
c Revision 1.52  1997/12/23  00:23:34  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.51  1997/12/17  23:22:45  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.50.1.1  1997/12/19  02:03:11  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.50  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.49  1997/03/06  04:35:28  mrd
c Force appropriate R or T truncation in call to ftospec, no matter how
c LLMAX is defined.
c
c Revision 1.48  1996/11/19  06:05:53  ldr
c Parallel Legendre transforms for kaos.
c
c Revision 1.47  1996/10/24  01:59:28  ldr
c Fix up merge of MRD/TIE changes; plmr not needed any more.
c
c Revision 1.46  1996/10/24  01:27:40  ldr
c Merge of TIE changes with LDR/MRD changes since V5-0-12.
c
c Revision 1.45  1996/08/12  01:51:22  mrd
c Generalise for triangular truncations other than T63
c
c Revision 1.44.1.1  1996/10/24  01:03:06  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.44  1996/06/13  02:07:15  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.43  1996/03/21  03:18:57  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.42  1995/10/17  05:08:27  ldr
c Made common block /curent/ the right length.
c
c Revision 1.41  1995/10/04  06:38:16  mrd
c Changes for more efficient writing of history files.
c
c Revision 1.40  1995/09/07  06:38:04  ldr
c Corrected padding on /legnd/.
c
c Revision 1.39  1995/08/31  04:30:46  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.38  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c

*VOCL TOTAL,REPEAT(999999)

      subroutine phys(tdt,exptyp)

      implicit none

!$OMP THREADPRIVATE ( /ENERX/ )
!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /PHYSICAL/ )
!$OMP THREADPRIVATE ( /STATS/ )
!$OMP THREADPRIVATE ( /VEGDAT/ )
!$OMP THREADPRIVATE ( /WORK1/ )
!$OMP THREADPRIVATE ( /WORKNS/ )

c Multiprocessor version of PHYS

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real tdt
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real uzon
      common/enerx/uzon(nl,2)

      real plmx,plm,cplm,pad
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)

      include 'PHYSICAL.f'

      real radst
      common/stats/radst(ln2,14:max_radstm)

      real veg
      common/vegdat/veg(ln2,4)

      include 'WORK1.f'

      real tfnr,tfni
      common/workns/tfnr(mw,nl,2),tfni(mw,nl,2)

C Global data blocks
      include 'CHMAP.f'
      include 'FLDRI.f'
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'HIST.f'
      include 'MASIV2.f'
      include 'MASIV4.f'
      include 'PRINTT.f'
      include 'RMGRID.f'
      include 'TIMEX.f'
      include 'UVPGD.f'

      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel

      integer lwh
      parameter (lwh=(lw+1)/2)
      real plmg,cplmg,plmgo,plmge,cplmgo,cplmge,rampm
      common/glegnd/plmg(lw1,mw,lat),cplmg(lw,mw,lat)
     &           ,plmgo(lat,lwh,mw),plmge(lat,lwh,mw)
     &           ,cplmgo(lat,lwh,mw),cplmge(lat,lwh,mw)
     &           ,rampm(mw,nl)

      real ronmx,sonmx
      integer iphysm,iphys
      common/masiv1/ronmx(ln2,nl,lat,2),sonmx(ln2,nl,lat,2)
     &             ,iphysm,iphys

      real utr,ur
      common/ubplm/utr(nl,lw),ur(lw,nl)

      real vegm
      common/vegdatm/vegm(ln2,4,lat)

C Local work arrays and variables
c-- Dimension following at mw+1 to prevent Bank Conflicts at T63
      real tfor(mw+1,nl,lat),tfoi(mw+1,nl,lat),tfer(mw+1,nl,lat),
     &     tfei(mw+1,nl,lat)
      real ufor(nl),ufer(nl)
      real utrl(nl,lw,lat)

      logical radstep,modmin

      integer ihist
      integer k
      integer lg
      integer lgn
      integer ll
      integer ll_lim
      integer mg
      integer nv
      integer i, j

      real x3

C Local data, functions etc
      include 'LLMAX.f'

C Start code : ----------------------------------------------------------

C**** ROUTINE TO CALCULATE ALL PHYSICAL ADJUSTMENTS,
C****   SURFACE FLUXES,RADIATIVE HEATING,ETC.

C***** FILL MAPPING ARRAYS WITH BLANKS BEFORE PHYSICAL LOOP
C***** IF LAST TIMESTEP BEFORE END OF DAY
*PDIR SERIAL
      if(nrad.eq.0)then
        radstep=.false.
      else
        radstep = mod(mins/int(mstep, 8),int(nrad, 8)).eq.0_8
      endif

      modmin=mod(mins+int(mstep, 8),1440_8).eq.0_8
      if(modmin.and.mlomap)call chset(i10d)
      if(modmin.and.gwicm)then
        call chset(chwb)
        call chset(chwg)
        call chset(snmap)
        call chset(snwice)
        call chset(snicel)
      endif
      if(modmin.and.cvrnm)then
        call chset(cnmlp)
        call chset(rainp)
      endif
      if((mod(mins+int(nrad*mstep, 8),1440_8).eq.0_8).and.cldm)then
        call chset(chlowc)
        call chset(chmidc)
        call chset(chhic)
        call chset(chtotc)
      endif

c Call initfs to calculate solar parameters for radiation calculation
c Argument 1 means this is called once each timestep
      if(nrad.ne.0)call initfs(1) 

c Initialize logical arrays to keep track of points that have just
c frozen or just melted. These arrays are used after the physics loop to
c adjust points at the equatorward latitude between SEA and MLO.
      if((nsteps.eq.0).and.(.not.qflux))call just_fm(1)
*PDIR ENDSERIAL

!$OMP  PARALLEL
!$OMP& PRIVATE (i, j, k, lg, lgn, ll, ll_lim, mg, mm, nv, ufer, ufor,
!$OMP&          x3)      
!$OMP& SHARED  (exptyp, impvor, iphys, lsm_type, pgd, plmg,
!$OMP&          plmge, plmgo, radstep, radstm, rgt, rmg, ronmx,
!$OMP&          savegrid, savehist, sonmx, tdt, tei, ter, tfei, tfer,
!$OMP&          tfoi, tfor, ugd, utr, utrl, vegm, vgd, w)

C****
C**** 1111 : MAIN GAUSS LOOP FOR PHYSICAL ADJUSTMENTS
C****

!$OMP DO SCHEDULE(DYNAMIC)

      do 1111 lg=1,lat
C***  THE NORTH & SOUTH GRID POINT VALUES ARE COMPUTED TOGETHER
c----
c     extract the appropriate legndre polys and store in plm()
c     for (lw,mw) and in plmx() for (lw1,mw).
      do 33 mm=1,mw 
      do 33 ll=1,llmax(mm)+1
   33 plmx(ll,mm) = plmg(ll,mm,lg) 
      do 330 mm=1,mw
      do 330 ll=1,llmax(mm)
  330 plm(ll,mm)=plmx(ll,mm)
      if(trunc.eq.'T')then
         do mm=2,mw
            plm(mm,lw+2-mm)=0.0
         end do
      end if

C**** TO INPUT ONE ROW(NORTH AND SOUTH) OF ACCUMULATED RADIATION DATA
C**** PLUS SURFACE DATA

!      print *, "before loop 39"

      do 39 nv=14,max_radstm
      do 39 mg=1,ln2
   39 radst(mg,nv)=radstm(mg,nv,lg)

!      print *, "after loop 39"

      do 392 nv=1,ngrid
        do 391 mg=1,ln2
          griddata(mg,nv)=savegrid(mg,nv,lg)
 391    continue
 392  continue

      if( lsm_type .eq. "nsib ") then
       do 400 k=1,4
       do 400 mg=1,ln2
       veg(mg,k)=vegm(mg,k,lg)
 400   continue
      endif

C**** TRANSFORM PHYSICS QUANTITIES ONTO GRID VIA SUB PTOG()

C     CALL TIMER('PTOG    ',3)
      if(mw.eq.64)then
        call ptogcr63    !T63 version (optimized for Cray)
      else
        call ptogcray    !R21/R42 version (optimized for Cray)
      endif
C     CALL TIMER('PTOG    ',4)

      do 58 k=1,nl
      do 58 mg=1,ln2
      rmnt(mg,k)=rgt(mg,k,lg)
   58 rmn(mg,k)=rmg(mg,k,lg)

C**** RADIN() : RADIATION , CONVECTIVE CONTRIBUTIONS
C     CALL TIMER('RADIN   ',3)
      call radin(lg,tdt,exptyp)
C     CALL TIMER('RADIN   ',4)

      do 59 k=1,nl
      do 59 mg=1,ln2
   59 rmg(mg,k,lg)=rmn(mg,k)

C**** SAVE UN,VN,PN FOR DYNAMICS STEP

      lgn=lat2p-lg
      do 40 k=1,nl
      do 40 mg=1,lon
      ugd(mg,lgn,k)=un(mg,k)
      ugd(mg,lg,k)=un(mg+lon,k)
      vgd(mg,lgn,k)=vn(mg,k)
   40 vgd(mg,lg,k)=vn(mg+lon,k)
      do 42 mg=1,lon
      pgd(mg,lgn)=pn(mg)
   42 pgd(mg,lg)=pn(mg+lon)

C**** TO OUTPUT ONE ROW OF ACCUMULATED RADIATION DATA TO STORAGE ARRAY
C**** PLUS SURFACE DATA

      do 66 nv=14,max_radstm
      do 66 mg=1,ln2
   66   radstm(mg,nv,lg)=radst(mg,nv)
      do 662 nv=1,ngrid
      do 662 mg=1,ln2
  662   savegrid(mg,nv,lg)=griddata(mg,nv)


C**   WRITE CURRENT GRID POINT FRICTION TO STORAGE

      do 67 k=1,nl
      do 67 mg=1,ln2
      ronmx(mg,k,lg,iphys)=ron(mg,k)
   67 sonmx(mg,k,lg,iphys)=son(mg,k)

c     If history is to be saved on this step, put
c      some data into global arrays :
      if(savehist(1).or.savehist(2))call hist_wlat(lg,radstep)

C**** SPECTRAL RESYNTHESIS :
C**    FFT TEMPS & MOISTURE FLUXES
C     CALL TIMER('PH MFFTM',3)
c   Physics FFT transforming 1 variables (nl) and for NH,SH at
c    same time : (nl)*2
c   Next argument (1) => Physics variables
c
c   The input data is stacked (ten etc) in a common block (work1)
c   The return data (tfnr,tfni etc) is also stacked in a common block
c    (workns). However, to facilitate vectorizing below, a change of
c    position of the last two arguments in (tfnr,tfni) relative to (ten)
c    etc is required. The change from [ten(ln2,nl) to tfnr(mw,nl,2) etc]
c    is achieved through precomputed pointer "igtos_p" :
c    see gauleg.f and mfftm.f

CSJP  Former machine dependence at this point
      call mfftma(ten, tfnr, tfni, nl*2, 1)

C     CALL TIMER('PH MFFTM',4)

C**** CREATE THE PHYSICALLY ADJUSTED TEMP FIELD AND MOISTURE FLUXES
C**   SPECTRALLY BY SEQUENTIALLY ADDING ODD AND EVEN COMPONENTS
C     CALL TIMER('PH SPADD',3)

      x3=w(lg)
      do 69 k=1,nl
      do 69 mm=1,mw
      tfor(mm,k,lg)=(tfnr(mm,k,1)+tfnr(mm,k,2))*x3
      tfoi(mm,k,lg)=(tfni(mm,k,1)+tfni(mm,k,2))*x3
      tfer(mm,k,lg)=(tfnr(mm,k,1)-tfnr(mm,k,2))*x3
      tfei(mm,k,lg)=(tfni(mm,k,1)-tfni(mm,k,2))*x3
 69   continue

      if(impvor)then

      do 691 k=1,nl
      ufor(k)=(uzon(k,1)+uzon(k,2))*x3/lon
      ufer(k)=(uzon(k,1)-uzon(k,2))*x3/lon
  691 continue
      do 731 ll=1,lw,2
      do 731 k=1,nl
      utrl(k,ll,lg)=ufor(k)*plm(ll,1)
  731 continue
      do 741 ll=2,lw,2
      do 741 k=1,nl
      utrl(k,ll,lg)=ufer(k)*plm(ll,1)
  741 continue

      end if

 1111 continue

!$OMP END DO

C**** END OF MAIN GAUSS LOOP


c Zero ter, tei arrays

!$OMP DO SCHEDULE(DYNAMIC)

      do 82 k=1,nl
CSJP        do ll=1,lw*mw
CSJP          ter(ll,1,k)=0.
CSJP          tei(ll,1,k)=0.
CSJP        enddo
        do 82 j = 1, mw
          do 82 i = 1, lw
            ter(i, j, k) = 0.0
            tei(i, j, k) = 0.0
 82   continue

!$OMP END DO

      if(impvor)then

!$OMP DO SCHEDULE(DYNAMIC)

      do ll=1,lw
        do lg=1,lat
          do k=1,nl
            utr(k,ll)=utr(k,ll)+utrl(k,ll,lg)
          enddo
        enddo
      enddo

!$OMP END DO

      endif

c Do the Legendre transform to re-create the spectral temperature field
c Parallelize over mm

c In a triangular model the truncation is enforced here, no matter how
c LLMAX is defined. This is necessary to maintain the truncation because
c plmg has an extra meridional wave

!$OMP DO SCHEDULE(DYNAMIC)

      do mm=1,mw
         ll_lim = lw1-mm
         if ( trunc .eq. 'R' ) ll_lim = lw

         call ftospec(ll_lim,mm,
     &            plmgo(1,1,mm),plmge(1,1,mm),
     &            tfor,tfoi,tfer,tfei,ter,tei)
      enddo

!$OMP END DO

!$OMP END PARALLEL

c----

*PDIR SERIAL

c Adjust IMSL between MLO and SEA if poleward point has just melted or frozen
      if(.not.qflux)call just_fm(2)

c Save history files
      do ihist=1,2
       if(savehist(ihist))then
        if(mod(mins+int(mstep, 8),int(hist_interval(ihist), 8)).eq.0_8)
     &    call hist_save(ihist, mstep, nrad)
        if(radstep .and. (mod(mins+int(nrad*mstep, 8),
     &                        int(hist_interval(ihist), 8)).eq.0_8))
     &    call hist_cld(ihist)
       endif
      enddo

*PDIR ENDSERIAL

      return
      end
C---------------------------------------------------------------------
      subroutine chset(charr)
      include 'PARAMS.f'
      character*1 charr(ln2,lat)
c.... Fill character array with blanks
      do 520 lg=1,lat
      do 520 mg=1,ln2
  520   charr(mg,lg)=' '
      return
      end
