c $Log: initax.f,v $
c Revision 1.19  2001/11/07 06:56:28  rot032
c Some further minor tuning from HBG
c
c Revision 1.18  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.17  2001/02/22 05:56:37  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.16  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.15  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.14  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.13  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.12  1998/12/10  00:55:38  ldr
c HBG changes to V5-1-21
c
c Revision 1.11  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.10  1996/10/24  01:02:58  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.9  1996/03/21  03:18:53  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.8  1995/11/30  04:40:16  ldr
c Final updates to qcloud scheme to bring it to long run h63 used for
c BMRC workshop and write-up.
c
c Revision 1.7  1995/11/23  06:03:30  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.6  1995/08/31  05:40:42  ldr
c Correct HBG's setting of aftsea to be consistent with LDR's value, and
c add diagnostic for jclb.
c
c Revision 1.5  1995/08/31  04:30:45  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.4  1993/12/17  15:32:59  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  92/12/09  14:43:38  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.2  91/03/13  12:58:49  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:34  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT
c     Input:   from common/cnsta in CNSTA.f
c                  sig  - sig - sigma values    sigk - sigma**cappa
c                  dsk  - sigma thicknesses (1 to nl) 
c                  algf - log(sig(k))
c
c              from common/gausl in GAUSL.f
c               coa - Array of sin(gauslat), running from North pole to equator
c
c     Output:  from common/dystb in this subroutine
c                  profile arrays:        wpb, wpc - winds
c                  qpb, qpc - moisture    tpb, tpc - temperature
c
c              from common/levdata in this subroutine
c                  aftsea - neutral drag coefficient over sea
c                  kbgel  - level at which shallow convection begins
c                  klcmc  - convective mapping indicator
c
c     In/Out:  from common/dystb in this subroutine
c                  profile arrays:
c                  qpa - moisture  tpa - temperature  wpa - winds
c
c              from common/levdata in this subroutine
c                  jclb   - convection cloud base
c                  klowr  - low rainfall (up to sig=0.9)
c                  kmidrn - middle rain (up to sig=0.45) 
c                  nlow   - top level for low cloud
c                  nmid   - top level for mid cloud
c
c 
      subroutine initax

C Global parameters
      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,jlp=jl+1,il2=il+2,jl2=jl+2)
      parameter(jl3=jl+3,il4=il+4,jl4=jl+4)
      include 'PHYSPARAMS.f'

C Argument list

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'TIMEX.f'
      common/dystb/tpa(nl),tpb(nl),tpc(nl),wpa(nl),wpb(nl),wpc(nl)
     &,qpa(nl),qpb(nl),qpc(nl),plev(nl),fred(nl)
      common/hverdat/roncp,rong,rlogs1,rlogs2,rlogh1,rlog12
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn
      common/savetrig/rlat(jl4),rlong(il4),cs(jl4),sn(jl4),dlong
     & ,snlo(il4),cslo(il4),ucon(jl4),vcona(jl4),vconb(jl4),vconc(jl4)

C Local work arrays and variables
      real alat1(jl)

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** SET UP INTERPOLATOR ARRAYS (A,B,C) FOR VARIOUS FIELDS
C**** FOR DYNMST
      do 13 k=2,nl

C**** TEMP uses profile A+B*log(sigma)
C**** GPH also uses profile A+B*log(sigma)
      tpa(k)=algf(k-1)-algf(k)
      tpb(k)=algf(k-1)/tpa(k)
      tpc(k)=algf(k)/tpa(k)

C**** VELOCITY uses profile A+B*sigma
      wpa(k)=sig(k-1)-sig(k)
      wpb(k)=sig(k-1)/wpa(k)
      wpc(k)=sig(k)/wpa(k)

C**** MOISTURE uses profile A+B*(sigma)**3
      qpa(k)=sig(k-1)**3-sig(k)**3
      qpb(k)=sig(k-1)**3/qpa(k)
   13 qpc(k)=sig(k)**3/qpa(k)

      algf0=0.0
      tpa(1)=algf0-algf(1)
      tpb(1)=algf0/tpa(1)
      tpc(1)=algf(1)/tpa(1)

      wpa(1)=1.0-sig(1)
      wpb(1)=1.0/wpa(1)
      wpc(1)=sig(1)/wpa(1)

      qpa(1)=1.0-sig(1)**3
      qpb(1)=1.0/qpa(1)
      qpc(1)=sig(1)**3/qpa(1)

      if(.not.statsflag)then
c----  If statsflag=F, set dummy plev values
        do k=1,nl
          plev(k)=1000.0*sig(k)
        enddo
      endif

c**** Set up level from which convection should take place 
c**** (for conv.f)
      IF((.not.ukconv).and.(.not.kuoflag))THEN

c.... Preset convective cloud base minimum at 500m (p/Ps=0.94)
      k=1
      f2=0.94-sig(k)
   12 k=k+1
      f1=f2
      f2=0.94-sig(k)
      if((f1.le.0.0).and.(f2.gt.0.0))go to 20
      go to 12
   20 jclb=k-1
      if(f2.lt.abs(f1))jclb=k
      write(6,'(a,i2)')' initax: Convective cloud base set to k =',jclb

      ENDIF

c.... Set up digital mapping value for convection. Cut off level for 
c.... indicating Low/Middle/Penetrating convection set at p/Ps=0.600 

      k=1
      f2=0.60-sig(k)
   14 k=k+1
      f1=f2
      f2=0.60-sig(k)
      if((f1.le.0.0).and.(f2.gt.0.0))go to 16
      go to 14
   16 klcmc=k-1
      if(f2.lt.abs(f1))klcmc=k

c**** Set up levels for low/mid/high cloud counting (for cloud2.f)

c.... "low" cloud up to half level closest to p/Ps=0.800 (~800mbs)
      k=1
      f2=0.800-sigh(k)
   22 k=k+1
      f1=f2
      f2=0.800-sigh(k)
      if((f1.le.0.0).and.(f2.gt.0.0))go to 30
      go to 22
   30 nlow=k-1
      if(f2.lt.abs(f1))nlow=k
c..   The FULL level is 1 below the half level indicator
      nlow=nlow-1

c.... "mid" cloud up to half level closest to p/Ps=0.400 (~400mbs)
      f2=0.400-sigh(k)
   32 k=k+1
      f1=f2
      f2=0.400-sigh(k)
      if((f1.le.0.0).and.(f2.gt.0.0))go to 40
      go to 32
   40 nmid=k-1
      if(f2.lt.abs(f1))nmid=k
c..   The FULL level is 1 below the half level indicator
      nmid=nmid-1

c**** Set up aftsea (for hsflux.f) for specific number of model levels
c**** This has to be tuned for a given model

      if(nl.lt.9)then
        aftsea=.0010
      elseif(nl.lt.18)then
        aftsea=.0012
      else
        aftsea=.0014
        if(ukconv)then
c         aftsea=0.0012
          aftsea=0.0014
          if(ncarpbl)then
            if(lw.eq.22)then
              aftsea=0.0010
            else
              aftsea=0.0008
            endif
          endif
        endif
      endif

c.. for a specific model use :
c     if(nl.eq.??) aftsea=.001?


c**** Set up digital map values for rainda.f

c.... Low rainfall (up to sig=0.9), Middle Rain (up to sig=0.45), High above
      klowrn=1
      kmidrn=nl-1
      do 24 k=1,nl
      if(sig(k).gt.0.9) klowrn=k
   24 if(sig(k).gt.0.45) kmidrn=k
      klowrn=max(klowrn,1)
      kmidrn=min(kmidrn,nl-1)


c**** Set up values used in hvertmx.f
      roncp=rdry/cp
      rong=rdry/grav
      rlogs1=log(sig(1))
      rlogs2=log(sig(2))
      rlogh1=log(sigh(2))
      rlog12=1./(rlogs1-rlogs2)
c- ksc in FEWFLAGS
      if(ukconv)then
        ksc=0 ! No shallow conv if UKMO convection
      else
        if(qcloud)then  !Use Tiedtke shallow convection scheme with ktop=ksc
          if(nl.gt.9)then
            ksc=nl/3
          else
            ksc=-1
          endif
        else            !Use Geleyn shallow convection scheme
          ksc=-1
        endif
      endif

c**** Set up any shallow convection (SC) values for hvertmx.f

      IF(ksc.gt.1)THEN

c.... For Tiedtke SC :
c Look for the lowest layer s.t. the top of the layer is above the LCL,
c and call that layer cloud base. ktied is level at which search starts.

c..   Use p/Ps=0.95 as the level to start down from.
      k=1
      f2=0.95-sig(k)
   42 k=k+1
      f1=f2
      f2=0.95-sig(k)
      if((f1.le.0.0).and.(f2.gt.0.0))go to 50
      go to 42
   50 ktied=k-1
      if(f2.lt.abs(f1))ktied=k
      write(6,'(a,i2)')' initax: Highest s-conv cloud base is k =',ktied

      ENDIF

      IF(ksc.lt.0)THEN

c.... For Geleyn SC : kbgel is level at which SC starts.

c..   Preset SC base at 500m (p/Ps=0.94)
      k=1
      f2=0.94-sig(k)
   52 k=k+1
      f1=f2
      f2=0.94-sig(k)
      if((f1.le.0.0).and.(f2.gt.0.0))go to 60
      go to 52
   60 kbgel=k-1
      if(f2.lt.abs(f1))kbgel=k

      ENDIF


c---- Set up stuff for SLT scheme (hadvect)
      do 80 j=1,lat
        alat1(j) = asin(coa(j)) ! Set up Gaussian lats
   80   alat1(2*lat+1-j) = - alat1(j)
      do 81 j=3,jl2
        rlat(j)=alat1(lat2p-(j-2)) ! alat1 runs from Nl to S
        cs(j)=cos(rlat(j))
   81   sn(j)=sin(rlat(j))
      rlat(2)=-pi-rlat(3)
      rlat(jl3)=pi-rlat(jl2)
      rlat(1)=-pi-rlat(4)
      rlat(jl4)=pi-rlat(jlp)
      dlong=2.*pi/il
      do 82 i=3,il2
        rlong(i)=(i-3)*dlong -pi
        snlo(i)=sin(rlong(i))
   82   cslo(i)=cos(rlong(i))
      dtx=2.0*dt
      if(ncepstrt.gt.0)then ! ncepstrt =2,1,0,0,0 ....
        dtx=0.5*mstep*60.0
      endif
      do 84 j=3,jl2
        ucon(j)=dtx/(erad*cs(j)*2.*dlong)
        vcona(j)=dtx*(rlat(j)-rlat(j+1))/
     &         (erad*(rlat(j-1)-rlat(j))*(rlat(j-1)-rlat(j+1)))
        vconb(j)=dtx*(2.*rlat(j)-rlat(j-1)-rlat(j+1))/
     &         (erad*(rlat(j)-rlat(j-1))*(rlat(j)-rlat(j+1)))
        vconc(j)=dtx*(rlat(j)-rlat(j-1))/
     &         (erad*(rlat(j+1)-rlat(j-1))*(rlat(j+1)-rlat(j)))
   84 continue

c.... set up height reduction factor for implicit vorticity
c.... (do not multiply ubx and dubx by fred(k) if not required)
        do 122 k=1,nl
c       fred(k)=((sig(1)-sig(k))/(sig(1)-sig(nl)))**4
c 122   if(fred(k).lt.0.05)fred(k)=0.0
        fred(k)=0.0
        if(sig(k).lt.0.2)
     &  fred(k)=((0.2-sig(k))/(0.2-sig(nl)))**2
  122   continue

      return
      end
