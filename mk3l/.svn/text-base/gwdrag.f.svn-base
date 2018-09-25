c Inserted code to detect floating point exceptions before they occur. If
c the square root of a negative number is about to be taken, the run is aborted
c and the values of key variables are written to standard output.
c Note that this is only done for the calculation of the Brunt-Vaisala
c frequency at the surface and for levels 1 to NL-1, as these are the only
c points at which I have ever encountered FPEs.
c SJP 2004/01/05
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: gwdrag.f,v $
c Revision 1.17  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.16  1998/12/10  00:55:41  ldr
c HBG changes to V5-1-21
c
c Revision 1.15  1997/12/17  23:22:52  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.14  1996/10/24  01:02:49  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.13  1996/06/13  02:06:36  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.12  1996/03/21  03:18:44  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.11  1996/02/19  04:09:49  ldr
c Generalize for 24 levels.
c
c Revision 1.10  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.9  1994/08/08  13:16:22  ldr
c Push the debug prints down into individual routines.
c
c Revision 1.8  94/08/04  16:55:29  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.7  93/12/17  15:32:40  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.6  93/07/06  16:30:55  ldr
c      delt(nl) now defined for nl levels - was only for 9 via data statement
c 
c Revision 1.5  92/12/09  14:43:28  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.4  92/06/16  12:02:37  ldr
c Rewritten to pass physical variables as arguments rather than in common,
c in order to implement sea-ice model.
c 
c Revision 1.3  91/05/02  17:06:23  igw
c limit topographic variance
c 
c Revision 1.2  91/03/13  12:58:27  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:28  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT:
c     Input:   from common/cnsta in CNSTA.f
c                  dsk - sigma thicknesses (1 to nl)
c
c              from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  insdebug - flag to control hemisphere debugging
c                  lgdebug  - flag to control latitude debugging
c                  mgdebug  - flag to control longitude debugging
c
c               from common/hybrpr in HYBRPR.f
c                   dprf  - pressure thickness at each sigma level
c                   prf   - pressure at full levels
c                   muf   - pressure derivative with respect to coordinates, at
c                           full levels
c                   pdpsk - (pressure/surface pressure)**cappa 
c
c               from arguments
c                   he  - variance of sub-grid scale topographys
c                   ins - hemisphere index       lg - latitude index
c                   pg  - surface pressure(mbs)  tg - surface temperature
c                   ttg - physical temperature   tdt- 2 X timestep (dt)
c                   u   - zonal wind             v  - meridional wind 
c
c     In/Out:   from arguments
c                   bvnf - Brunt-Vaisala frequency, full levels
c                   ron  - measure of east-west stress at surface
c                   son  - measure of north-south stress at surface
c                   dkegw - change in KE  
c 
      subroutine gwdrag(lg,tdt,he,tg,ttg,pg,u,v,
     &                  ron,son,bvnf,dkegw)

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real tdt
      real he(ln2)
      real tg(ln2)
      real ttg(ln2,nl)
      real pg(ln2)
      real u(ln2,nl)
      real v(ln2,nl)
      real ron(ln2,nl)
      real son(ln2,nl)
      real bvnf(ln2,nl)
      real dkegw(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'

C Local work arrays and variables
      real wmag(ln2)
      real uu(ln2,nl),thf(ln2,nl),dthdz(ln2,nl),fni(ln2,nl)
      real bvng(ln2),temp(ln2),fnii(ln2),unst(ln2),apuw(ln2),apvw(ln2)
     & ,alam(ln2)
      real dzi(ln2,nl),delt(nl)
      real hel(ln2)
      real temp_sjp

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do 400 k=1,nl-1
  400 delt(k)=1.0
      delt(nl)=0.0
C****
C**** GRAVITY WAVE DRAG - VECTORISABLE - FOR ALL POINTS.
C****
C**** FULL LEVEL VELOCITIES HAVE BEEN UPDATED BY HVERTMX
C**** FORM NEW WMAG AT SURFACE

      do 681 mg=1,ln2
  681 wmag(mg)=max(sqrt(u(mg,1)*u(mg,1)+v(mg,1)*v(mg,1)),1.0)

c**** limit launching height : Palmer et al use limit on variance of
c**** (400 M)**2. Launching height is 2*Std Dev. This limit on 
c**** launching height is then 2*400=800 M. This may be a bit severe.
c**** According to Palmer this prevents 2-grid noise at steep edge of
c**** Himalayas etc.
      do 682 mg=1,ln2
  682 hel(mg)=min(he(mg),800.0)

c---- This dz is being computed elsewhere!
      dzx=grav/rdry
      do 683 k=1,nl
      do 683 mg=1,ln2
  683 dzi(mg,k)=(prf(mg,k)/dprf(mg,k))*dzx/ttg(mg,k)

      do 700 k=1,nl
      do 700 mg=1,ln2
      uumg=(u(mg,k)*u(mg,1)+v(mg,k)*v(mg,1))/wmag(mg)
  700 uu(mg,k)=max(0.0,uumg)
C**** SET UU() TO ZERO ABOVE IF UU() ZERO BELOW
C**** UU>0 AT K=1, UU>=0 AT K=2 - ONLY SET FOR K=3 TO Nl
      do 702 k=3,nl
      do 702 mg=1,ln2
  702 if(uu(mg,k-1).eq.0.0)uu(mg,k)=0.0
C**** CALCULATE BRUNT-VAISALA FREQUENCY
C**** SURFACE BVNG() AND FULL LEVELS BVNF(,)
C**      PUT THETA IN THF()
      do 7031 k=1,nl
      do 7031 mg=1,ln2
 7031 thf(mg,k)=ttg(mg,k)/pdpsk(mg,k)
      do 703 mg=1,ln2
C**** IF UNSTABLE (SURFACE TO BL) THEN NO GWD (SET UNST()=0)
      unst(mg)=1.0
  703 fnii(mg)=tg(mg)+0.1
      do 704 mg=1,ln2
      if(tg(mg).gt.thf(mg,1))unst(mg)=0.0
C**      CALC D(TH)/DZ
  704 if(thf(mg,1).lt.fnii(mg))thf(mg,1)=fnii(mg)
      do 25 k=2,nl
      do 705 mg=1,ln2
  705 fnii(mg)=thf(mg,k-1)+0.1
      do 706 mg=1,ln2
  706 if(thf(mg,k).lt.fnii(mg))thf(mg,k)=fnii(mg)
   25 continue
      do 707 mg=1,ln2
  707 dthdz(mg,1)=(thf(mg,1)-tg(mg))*dzi(mg,1)
      do 30 k=2,nl
      do 7071 mg=1,ln2
 7071 dthdz(mg,k)=(thf(mg,k)-thf(mg,k-1))*dzi(mg,k)
   30 continue
C**    CALC BVNG()
      do mg = 1, ln2
        temp_sjp = grav * dthdz(mg, 1) / tg(mg) 
        if (temp_sjp .lt. 0.0) then
          write (*, *)
          write (*, *) "ABORTING: Fatal error in GWDRAG"
          write (*, *)
          write (*, *) "mg = ", mg
          write (*, *) "lg = ", lg
          write (*, *)
          write (*, *) "thf(mg, 1)  = ", thf(mg, 1)
          write (*, *) "tg(mg)      = ", tg(mg)
          write (*, *) "prf(mg, 1)  = ", prf(mg, 1)
          write (*, *) "dprf(mg, 1) = ", dprf(mg, 1)
          write (*, *) "ttg(mg, 1)  = ", ttg(mg, 1)
          write (*, *)
          stop
        end if
        bvng(mg) = sqrt(temp_sjp)
      end do
CSJP      do 708 mg=1,ln2
CSJP  708 bvng(mg)=sqrt(grav*dthdz(mg,1)/tg(mg))
C**** LIMIT VALUE OF BVNG
      do 709 mg=1,ln2
c-709 if(bvng(mg).gt.0.01)bvng(mg)=0.01
  709 if(bvng(mg).gt.0.1)bvng(mg)=0.1
C**    CALC BVNF(,)
      do k = 1, nl-1
        do mg = 1, ln2
          temp_sjp = 0.5 * grav * (dthdz(mg, k) + dthdz(mg, k+1))
     &                            / thf(mg, k)
          if (temp_sjp .lt. 0.0) then
            write (*, *)
            write (*, *) "ABORTING: Fatal error in GWDRAG"
            write (*, *)
            write (*, *) "mg = ", mg
            write (*, *) "lg = ", lg
            write (*, *) "k  = ", k
            write (*, *)
            write (*, *) "thf(mg, k-1)  = ", thf(mg, k-1)
            write (*, *) "thf(mg, k)    = ", thf(mg, k)
            write (*, *) "thf(mg, k+1)  = ", thf(mg, k+1)
            write (*, *) "prf(mg, k)    = ", prf(mg, k)
            write (*, *) "prf(mg, k+1)  = ", prf(mg, k+1)
            write (*, *) "dprf(mg, k)   = ", dprf(mg, k)
            write (*, *) "dprf(mg, k+1) = ", dprf(mg, k+1)
            write (*, *) "ttg(mg, k)    = ", ttg(mg, k)
            write (*, *) "ttg(mg, k+1)  = ", ttg(mg, k+1)
            write (*, *)
            stop
          end if
          bvnf(mg, k) = sqrt(temp_sjp)
        end do
      end do
CSJP      do 40 k=1,nl-1
CSJP      do 710 mg=1,ln2
CSJP      bvnf(mg,k)=sqrt(  (0.5*grav)*(dthdz(mg,k)+dthdz(mg,k+1))/
CSJP     &  thf(mg,k) )
CSJP  710 continue
CSJP   40 continue
      fc2=0.5
      do 711 mg=1,ln2
      bvnf(mg,nl)=sqrt(  grav*dthdz(mg,nl)/thf(mg,nl) )
C**** FROUDE NUMBER CALCS
C****   CALCULATE (FC/F)**2 WHERE FC**2=FC2=0.5
C**     CALC FC2*(T*/Nl*/WMAG)/HE**2
C**   NOTE - USE (HE+10**-5)**2 - PREVENTS DIVISION BY ZERO (HE>>10**-5)
      temp(mg)=( fc2*tg(mg)/(bvng(mg)*wmag(mg)) ) /
     & ( (hel(mg)+1.0e-05)*(hel(mg)+1.0e-05) )
  711 continue
C****     CALC MAX(1-DELTA*FC**2/F**2,0) : PUT IN FNI()
      do 712 k=1,nl
      do 712 mg=1,ln2
      xxx      =(  (prf(mg,k)/pg(mg))*temp(mg)/
     & (bvnf(mg,k)*ttg(mg,k)) )
     &  *uu(mg,k)*uu(mg,k)*uu(mg,k)
  712 fni(mg,k)=max(0.0,1.0-delt(k)*xxx)
C**** FORM INTEGRAL OF ABOVE*UU**2 FROM SIG=1 TO SIG=0
      do 716 mg=1,ln2
  716 fnii(mg)=uu(mg,1)*uu(mg,1)*fni(mg,1)*muf(mg,1)*dsk(1)
      do 717 k=2,nl
      do 717 mg=1,ln2
  717 fnii(mg)=fnii(mg)+uu(mg,k)*uu(mg,k)*fni(mg,k)*muf(mg,k)*dsk(k)
C**     IF INTEGRAL=0.0, RESET TO SOME +IVE VALUE
      do 718 mg=1,ln2
  718 if(fnii(mg).eq.0.0)fnii(mg)=1.0
C**** FORM ALAM=G.ALPHA.HE.RHO*Nl*.WMAG/INTEGRAL(ABOVE)
c     alpha=0.01
      alpha=0.0075
      alpha=alpha*grav/rdry
      do 719 mg=1,ln2
      alam(mg)=(alpha*hel(mg)*pg(mg)*bvng(mg)*wmag(mg)
     & /(tg(mg)*fnii(mg))) *unst(mg)
  719 continue
C**** FORM FNI=ALAM*MAX(--,0) AND
C**** SOLVE FOR UU AT T+1 (IMPLICIT SOLUTION)
      do 820 k=1,nl
      do 820 mg=1,ln2
      xxx=sqrt(1.0+4.0*tdt*alam(mg)*fni(mg,k)*uu(mg,k))
      uu(mg,k)=2.0*uu(mg,k)/(1.0+xxx)
  820 continue

C**** DEFINE APUW=ALAM.U1/WMAG , APVW=ALAM.V1/WMAG
      do 821 mg=1,ln2
      apuw(mg)=alam(mg)*u(mg,1)/wmag(mg)
  821 apvw(mg)=alam(mg)*v(mg,1)/wmag(mg)

C**** FORM DV/DT DUE TO GW-DRAG AT EACH LEVEL
C**** = -ALAM.V*/WMAG.UU(T+1)**2.MAX(--,0)
C**** (NOTE : NO P*.COS(LAT)/RAD WEIGHTING YET IF SIGMA)
C**** (NOTE : NO MUF.COS(LAT)/RAD WEIGHTING YET IF HYBRID)
      do 720 k=1,nl
      do 720 mg=1,ln2
      xxx=uu(mg,k)*uu(mg,k)*fni(mg,k)
c     xxx=0.
      ron(mg,k)=-apuw(mg)*xxx
      son(mg,k)=-apvw(mg)*xxx
      dkegw(mg,k)=u(mg,k)*ron(mg,k)+v(mg,k)*son(mg,k)
  720 continue

      if(debug)then
        if(lg.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,i1)')'After gwdrag. IPASS = 1'
          write(25,99)'ron ',(ron(mg,k),k=1,nl)
          write(25,99)'son ',(son(mg,k),k=1,nl)
          write(25,*)
        endif
      endif
 99   format(a,30g10.3)

      return
      end
