c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Numerous corrections to references to arrays WRK1I, WRK1R, WRK2I and WRK2R,
c as these had not been modified to reflect the change in the ordering of the
c indices relative to semii.f
c SJP 2003/04/23
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/04/02
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: semiis.f,v $
c Revision 1.8  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.7  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.6  1997/01/10  06:06:36  ldr
c Replace pr and pi by prr and pri respectively, so that these routines work
c again.
c
c Revision 1.5  1996/10/24  01:03:30  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1996/06/13  02:08:54  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.3  1996/03/21  03:19:15  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.2  1993/12/17  15:34:43  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  93/12/06  16:53:32  ldr
c Initial revision
c
c     Input:   from common/cnsta in CNSTA.f
c                  dsk - sigma thicknesses (1 to nl) 
c
c              from common/cnste in CNSTE.f
c                  flm2 - l(l+1) part of spectral form
c                  flm2i - reciprocal of flm2 
c
c              from common/fewflags in FEWFLAGS.f
c                  chenflag - flag for CHEN version of model
c
c              from common/matzer in this subroutine
c                  ada,ag,amd,jm,uk2di,uk2i -      )semi implicit time 
c                  gmd - matrix for constant tmean )integration matrices
c                  iagji - semi-implicit time integration matrices
c
c              from common/timex in TIMEX.f
c                  nsteps - time step counter
c
c              from common/worka in WORKA.f
c                  tmean - preset isothermal temperature
c
c              from arguments
c                  tdt - 2 X timestep (dt)
c
c     In/Out:  from common/fldmri in FLDMRI.f
c                  pmi - spectral surface pressure, imaginary part
c                  pmr - spectral surface pressure, real part
c                  psimi - stream function, imaginary part
c                  psimr - stream function, real part
c                  temi - spectral temp, imaginary part at t-1
c                  temr - spectral temp, real part at t-1
c                  xhimi - velocity potential at t-1, imaginary part
c                  xhimr - velocity potential at t-1, real part
c
c              from common/fldri in FLDRI.f
c                  pi - spectral pressure, imaginary part
c                  pr - spectral pressure, real part
c                  psii - spectral stream function, imaginary part
c                  psir - spectral stream function, real part
c                  tei - spectral temp, imaginary part
c                  ter - spectral temp, real part
c                  xhii - velocity potental, imaginary part
c                  xhir - velocity potental, real part
c

c This is the parallel SGI ("Seca") version of semii. The common blocks
c should be identical to semii above. The logical structure is largely 
c unchanged, except that three relatively long loops over mm have been
c "factored out" and the array indexing on the work arrays has been altered.
c On a non-parallel scalar machine, the usual semii should really be used,
c but the penalty for using this one will be slight (LDR 12/93).

      subroutine semiis(tdt)

C Global parameters
      include 'PARAMS.f'
      parameter (lmx=lw+mw-2)
      include 'PHYSPARAMS.f'

C Argument list
      real tdt

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'CNSTE.f'
      include 'FEWFLAGS.f'
      include 'FLDRI.f'
      include 'FLDMRI.f'
      include 'TIMEX.f'
      include 'WORKA.f'
      include 'WORKRID.f'
      real uk2di(nl,nl),uk2i(nl,nl),jm(nl,nl),gmd(nl,nl),ag(nl,nl)
     & ,amd(nl,nl),ada(nl,nl),iagji(0:lmx,nl,nl)
      common/matzer/uk2di,uk2i,jm,gmd,ag,amd,ada,iagji

C Local work arrays and variables
      real rwrk(lw,mw),rwrk2(lw,mw),temp(nl)
      save rwrk,rwrk2,temp
      real wrk1r(lw,nl,mw),wrk1i(lw,nl,mw)
     & ,wrk2r(lw,nl,mw),wrk2i(lw,nl,mw)
      real c2r

C Local data, functions etc
      include 'LLMAX.f'

C Start code : ----------------------------------------------------------

CHEN
      t00ec=288.0
c---- t00ec replaces tmean(k)
CHEN
C*    INTEGRATION WITH ASSELIN FILTER
C*    NOTE : INCLUDE ONLY NON T+1 ASSELIN FILTER TERMS FOR TEMP
C**** AND MOISTURE.
      asf=0.05
      asfx=1.0-2.0*asf
      if(nsteps.le.2)then
        do k=1,nl
          temp(k)=2.0*rdry*tmean(k)
          if(chenflag)temp(k)=2.0*rdry*t00ec
        enddo
        g2=tdt*0.5
        g1=g2/eradsq
        do 8 mm=1,mw
        do 8 ll=1,llmax(mm)
        rwrk2(ll,mm)=g1*flm2(ll,mm)
 8      rwrk(ll,mm)=tdt*flm2i(ll,mm)
      endif

C****
C**** TIME INTEGRATE THE MOISTURE EQUATION (NON-IMPLICIT)
C**** Semi Lagrangian version now done in JMCGSLT
C****
C****
C**** TIME INTEGRATE THE VORTICITY EQUATION (NON IMPLICIT)
C****

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (j, k, ll, mm, xxx)
!$OMP& SHARED  (ada, ag, amd, asf, asfx, flm2, ipi, ipr, iti, itr, ixi,
!$OMP&          ixr, jm, ncepstrt, pmi, pmr, psimi, psimr, psii, psir,
!$OMP&          rwrk, rwrk2, tdt, temi, temp, temr, uk2i, wrk1i, wrk1r,
!$OMP&          wrk2i, wrk2r, xhimi, xhimr)

      do 9990 mm=1,mw

      do 10 k=1,nl
      do 10 ll=1,llmax(mm)
      ipr(ll,mm,k)=psimr(ll,mm,k)-rwrk(ll,mm)*ipr(ll,mm,k)
   10 ipi(ll,mm,k)=psimi(ll,mm,k)-rwrk(ll,mm)*ipi(ll,mm,k)
      do 11 k=1,nl
      do 111 ll=1,llmax(mm)
      wrk1r(ll,k,mm)=uk2i(k,1)*ipr(ll,mm,1)
  111 wrk1i(ll,k,mm)=uk2i(k,1)*ipi(ll,mm,1)
      do 12 j=2,nl
      do 112 ll=1,llmax(mm)
      wrk1r(ll,k,mm)=wrk1r(ll,k,mm)+uk2i(k,j)*ipr(ll,mm,j)
  112 wrk1i(ll,k,mm)=wrk1i(ll,k,mm)+uk2i(k,j)*ipi(ll,mm,j)
   12 continue
c....
c.... NOTES ABOUT STARTING FROM NCEP. (ncepagcm=T, ncepstrt used)
c.... When starting from NCEP, there are two half steps
c.... to start the model.
c....
c....   ncepstrt=2 is the first half step
c.... with variables at Tau=-0.5 and Tau=0 (using same data, but Tau=0
c.... will have been adjusted by Physics). The leapfrog time integration
c.... then gives new variables at Tau+0.5
c.... These are updated  in the usual way (eg xhi replaces xhim
c....  with asselin filtering, and xhip replaces xhi  etc.)
c....
c....   ncepstrt=1 is the second half step.
c.... Variables are at Tau=0 and Tau=0.5. The leapfrog time integration
c.... gives new variables at Tau+1. We next want to proceed with the
c.... usual leapfrog steps using full timesteps. I.E. From Tau=0 and
c.... Tau=1 to Tau=2. So do not retain the Tau=0.5 variables, and only
c.... update the Tau+1 variables (eg keep xhim, and xhip replaces xhi etc.
c.... Note that there is no asselin filtering on this step, so be
c.... careful with the Tau+1 part of the Asselin filter for Temp and
c.... Moisture which is done after the next Physical adjustments, in
c.... subroutine assel.f. There was no Asselin filter in semii.f for
c.... ncepstrt=1, so do not apply in assel.f when ncepstrt=0 (ncepstrt
c.... will have been updated by -1 at that stage). There is no Asselin
c.... filter in assel.f when ncepstrt=2 (initial data time).
c....
c....   ncepstrt=0,-1,-2,-3,.... are the usual leapfrog + Asselin filter
c.... timesteps.
c.... IF NOT STARTING FROM NCEP (normal AGCM in climate model mode)
c.... then ncepstrt=-1,-2,-3... and model code reverts to climate model
c.... coding.
c....
      if(ncepstrt.eq.1)then
       do ll=1,llmax(mm)
        psir(ll,mm,k)=wrk1r(ll,k,mm)
        psii(ll,mm,k)=wrk1i(ll,k,mm)
       enddo
      else
       do ll=1,llmax(mm)
        psimr(ll,mm,k)=asfx*psir(ll,mm,k)+
     &   asf*(psimr(ll,mm,k)+wrk1r(ll,k,mm))
        psir(ll,mm,k)=wrk1r(ll,k,mm)
        psimi(ll,mm,k)=asfx*psii(ll,mm,k)+
     &   asf*(psimi(ll,mm,k)+wrk1i(ll,k,mm))
        psii(ll,mm,k)=wrk1i(ll,k,mm)
       enddo
      endif
   11 continue
C****
C**** SEMI-IMPLCIT TIME INTEGRATION
C****
C***  SPECIAL CASE WHEN M=0(MM=1),L=0(LL=1) & => MEAN VORTICITY,
C**** DIVERGENCE UNCHANGING. CHANGES IN MEAN TEMPS AND MOISTURE
C**** ALSO CALCULATED.
C****
C****  DIVERGENCE,TEMPERATURE, AND SURFACE PRESSURE
C****
C**      FORM T-1 DIVERGENCE (PLACE IN XHIM)
      do 16 k=1,nl
      do 16 ll=1,llmax(mm)
      xhimr(ll,mm,k)=-flm2(ll,mm)*xhimr(ll,mm,k)
      xhimi(ll,mm,k)=-flm2(ll,mm)*xhimi(ll,mm,k)
      wrk2r(ll,k,mm)=flm2(ll,mm)*xhimr(ll,mm,k)
   16 wrk2i(ll,k,mm)=flm2(ll,mm)*xhimi(ll,mm,k)
c**   form the t+1 divergence components
c**    form xm()=iagjm()*dvgm()=(um-fllp*(ag+jm))*dvgm - in wrk1
      do 19 k=1,nl
      do 118 ll=1,llmax(mm)
      wrk1r(ll,k,mm)=xhimr(ll,mm,k)
  118 wrk1i(ll,k,mm)=xhimi(ll,mm,k)
      do 18 j=1,nl
      xxx=ag(k,j)+jm(k,j)
      do 119 ll=1,llmax(mm)
      wrk1r(ll,k,mm)=wrk1r(ll,k,mm)-xxx*wrk2r(ll,j,mm)
  119 wrk1i(ll,k,mm)=wrk1i(ll,k,mm)-xxx*wrk2i(ll,j,mm)
   18 continue
   19 continue

c**    form atm()+tdt*ay()+temp()*pmx - hold in wrk2
      do 20 k=1,nl
      do 120 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=temp(k)*pmr(ll,mm)
  120 wrk2i(ll,k,mm)=temp(k)*pmi(ll,mm)
      do 20 j=1,nl
      xxx=tdt*amd(k,j)
      do 121 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=wrk2r(ll,k,mm)+ada(k,j)*temr(ll,mm,j)
     &  +xxx*itr(ll,mm,j)
      wrk2i(ll,k,mm)=wrk2i(ll,k,mm)+ada(k,j)*temi(ll,mm,j)
     &  +xxx*iti(ll,mm,j)
  121 continue
   20 continue
c**   form rhs of si divergence eqn - hold in wrk1
      do 122 k=1,nl
      do 122 ll=1,llmax(mm)
      wrk1r(ll,k,mm)=wrk1r(ll,k,mm)+tdt*ixr(ll,mm,k) 
     & +rwrk2(ll,mm)*wrk2r(ll,k,mm)
      wrk1i(ll,k,mm)=wrk1i(ll,k,mm)+tdt*ixi(ll,mm,k) 
     & +rwrk2(ll,mm)*wrk2i(ll,k,mm)
  122 continue

 9990 continue 

!$OMP END PARALLEL DO

      do 22 k=1,nl
      wrk1r(1,k,1)=0.0
   22 wrk1i(1,k,1)=0.0

C**   SOLVE IAGJ()*DVGP()=DVGX()   (WRK1=DVGX)
C**   THIS IS NOT FULLY VECTORISABLE AT PRESENT
C***      ls=2
C***      do 24 mm=1,mw
C***      if(mm.eq.2)ls=1
C***      do 24 ll=ls,lw
C***      l=mm+ll-2
C***      do 26 k=1,nl
C***      wrk2r(ll,k,mm)=iagji(l,k,1)*wrk1r(ll,1,mm)
C***      wrk2i(ll,k,mm)=iagji(l,k,1)*wrk1i(ll,1,mm)
C***      do 26 j=2,nl
C***      wrk2r(ll,k,mm)=wrk2r(ll,k,mm)+iagji(l,k,j)*wrk1r(ll,j,mm)
C***   26 wrk2i(ll,k,mm)=wrk2i(ll,k,mm)+iagji(l,k,j)*wrk1i(ll,j,mm)
C***   24 continue
C***      do 27 k=1,nl
C***        wrk2r(1,k,1)=0.0
C***        wrk2i(1,k,1)=0.0
C*** 27   continue

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (j, k, ll, mm)
!$OMP& SHARED  (iagji, wrk1i, wrk1r, wrk2i, wrk2r)

      do 9995 mm=1,mw 

      do 26 k=1,nl
      do 126 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=iagji(ll+mm-2,k,1)*wrk1r(ll,1,mm)
  126 wrk2i(ll,k,mm)=iagji(ll+mm-2,k,1)*wrk1i(ll,1,mm)
      do 24 j=2,nl
      do 24 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=wrk2r(ll,k,mm)
     & +iagji(ll+mm-2,k,j)*wrk1r(ll,j,mm)
      wrk2i(ll,k,mm)=wrk2i(ll,k,mm)
     & +iagji(ll+mm-2,k,j)*wrk1i(ll,j,mm)
   24 continue
   26 continue

 9995 continue

!$OMP END PARALLEL DO

      do 27 k=1,nl
      wrk2r(1,k,1)=0.0
   27 wrk2i(1,k,1)=0.0

c**   wrk2 now holds divergence at t+1
c**   form implicit divergence at t ( hold in wrk1)

      c2r=prr(1,1)

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (c1i, c1r, j, k, ll, mm, xxx)
!$OMP& SHARED  (asf, asfx, dsk, flm2i, gmd, iti, itr, ncepstrt, pmi,
!$OMP&          pmr, pri, prr, tdt, tei, temi, temr, ter, uk2di, wrk1i,
!$OMP&          wrk1r, wrk2i, wrk2r, xhii, xhimi, xhimr, xhir)

      do 9999 mm=1,mw

      do 28 k=1,nl
      do 28 ll=1,llmax(mm)
      wrk1r(ll,k,mm)=0.5*(xhimr(ll,mm,k)+wrk2r(ll,k,mm))
   28 wrk1i(ll,k,mm)=0.5*(xhimi(ll,mm,k)+wrk2i(ll,k,mm))
c**   update the xhi fields
      if(ncepstrt.eq.1)then
       do k=1,nl
       do ll=1,llmax(mm)
        xhir(ll,mm,k)=-flm2i(ll,mm)*wrk2r(ll,k,mm)
        xhii(ll,mm,k)=-flm2i(ll,mm)*wrk2i(ll,k,mm)
        xhimr(ll,mm,k)=-flm2i(ll,mm)*xhimr(ll,mm,k) ! T-1 Div to xhi
        xhimi(ll,mm,k)=-flm2i(ll,mm)*xhimi(ll,mm,k) ! T-1 Div to xhi
       enddo
       enddo
      else
       do k=1,nl
       do ll=1,llmax(mm)
        c1r=-flm2i(ll,mm)*wrk2r(ll,k,mm)
        xhimr(ll,mm,k)=asfx*xhir(ll,mm,k)+
     &   asf*(c1r-flm2i(ll,mm)*xhimr(ll,mm,k))
        xhir(ll,mm,k)=c1r
        c1i=-flm2i(ll,mm)*wrk2i(ll,k,mm)
        xhimi(ll,mm,k)=asfx*xhii(ll,mm,k)+
     &   asf*(c1i-flm2i(ll,mm)*xhimi(ll,mm,k))
        xhii(ll,mm,k)=c1i
       enddo
       enddo
      endif
c****
c**** compute the si term for the temp equation
c****
      do 32 k=1,nl
      do 132 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=gmd(k,1)*wrk1r(ll,1,mm)
  132 wrk2i(ll,k,mm)=gmd(k,1)*wrk1i(ll,1,mm)
      do 33 j=2,nl
       xxx=gmd(k,j)
      do 33 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=wrk2r(ll,k,mm)+xxx*wrk1r(ll,j,mm)
   33 wrk2i(ll,k,mm)=wrk2i(ll,k,mm)+xxx*wrk1i(ll,j,mm)
   32 continue
      do 34 k=1,nl
      do 34 ll=1,llmax(mm)
      itr(ll,mm,k)=temr(ll,mm,k)+tdt*(itr(ll,mm,k)-
     &  wrk2r(ll,k,mm))
      iti(ll,mm,k)=temi(ll,mm,k)+tdt*(iti(ll,mm,k)-
     &  wrk2i(ll,k,mm))
   34 continue
      do 35 k=1,nl
      do 136 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=uk2di(k,1)*itr(ll,mm,1)
  136 wrk2i(ll,k,mm)=uk2di(k,1)*iti(ll,mm,1)
      do 35 j=2,nl
      do 35 ll=1,llmax(mm)
      wrk2r(ll,k,mm)=wrk2r(ll,k,mm)+uk2di(k,j)*
     &   itr(ll,mm,j)
      wrk2i(ll,k,mm)=wrk2i(ll,k,mm)+uk2di(k,j)*
     &   iti(ll,mm,j)
   35 continue
c**   no vert diffusion of mean  p*.t'
      if(mm.eq.1)then
        do 36 k=1,nl
        wrk2r(1,k,mm)=itr(1,mm,k)
   36   wrk2i(1,k,mm)=iti(1,mm,k)
      endif
c**   update the temp arrays
      if(ncepstrt.eq.1)then
       do k=1,nl
       do ll=1,llmax(mm)
        ter(ll,mm,k)=wrk2r(ll,k,mm)
        tei(ll,mm,k)=wrk2i(ll,k,mm)
       enddo
       enddo
      else
       do k=1,nl
       do ll=1,llmax(mm)
        temr(ll,mm,k)=asfx*ter(ll,mm,k)+asf*temr(ll,mm,k)
        ter(ll,mm,k)=wrk2r(ll,k,mm)
        temi(ll,mm,k)=asfx*tei(ll,mm,k)+asf*temi(ll,mm,k)
        tei(ll,mm,k)=wrk2i(ll,k,mm)
       enddo
       enddo
      endif
c****
c**** time integrate the surface pressure
c****
      do 139 ll=1,llmax(mm)
      wrk2r(ll,1,mm)=dsk(1)*wrk1r(ll,1,mm)
  139 wrk2i(ll,1,mm)=dsk(1)*wrk1i(ll,1,mm)
      do 40 k=2,nl
      do 40 ll=1,llmax(mm)
      wrk2r(ll,1,mm)=wrk2r(ll,1,mm)+dsk(k)*wrk1r(ll,k,mm)
   40 wrk2i(ll,1,mm)=wrk2i(ll,1,mm)+dsk(k)*wrk1i(ll,k,mm)
c**     update the p* arrays
      if(ncepstrt.eq.1)then
       do ll=1,llmax(mm)
        prr(ll,mm)=pmr(ll,mm)-tdt*wrk2r(ll,1,mm)
        pri(ll,mm)=pmi(ll,mm)-tdt*wrk2i(ll,1,mm)
       enddo
      else
       do ll=1,llmax(mm)
        c1r=pmr(ll,mm)-tdt*wrk2r(ll,1,mm)
        pmr(ll,mm)=asfx*prr(ll,mm)+asf*(pmr(ll,mm)+c1r)
        prr(ll,mm)=c1r
        c1i=pmi(ll,mm)-tdt*wrk2i(ll,1,mm)
        pmi(ll,mm)=asfx*pri(ll,mm)+asf*(pmi(ll,mm)+c1i)
        pri(ll,mm)=c1i
       enddo
      endif

 9999 continue 

!$OMP END PARALLEL DO

      pmr(1,1)=c2r
      pmi(1,1)=0.0
      prr(1,1)=c2r
      pri(1,1)=0.0
      return
      end
