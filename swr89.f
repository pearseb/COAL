c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: swr89.f,v $
c Revision 1.26  1998/12/10 00:55:59  ldr
c HBG changes to V5-1-21
c
c Revision 1.25  1997/12/17  23:23:01  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.24  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.22  1996/02/13  23:27:11  mrd
c Create a new index array ibn to get around bounds problems.
c
c Revision 1.21  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.19.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.20  1994/08/08  17:22:52  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.19  94/05/27  16:54:06  ldr
c Make some variables double precision for Seca to avoid underflows.
c 
c Revision 1.18  94/05/03  11:57:26  ldr
c Simplify conversion of fluxes to heating rates (from MRD).
c 
c Revision 1.17  93/12/17  15:34:03  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.16  93/08/19  15:11:03  ldr
c Minor cosmetic changes.
c 
c Revision 1.15  93/07/20  10:08:00  ldr
c Use block data for initializing data in common, to keep the VP happy.
c 
c Revision 1.14  93/06/23  14:30:37  mrd
c Made padding of vtemp common block resolution independent
c 
c Revision 1.13  93/06/16  15:44:33  ldr
c Minor change for energy conservation calculation (HBG).
c 
c Revision 1.12  92/12/09  14:44:38  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.11  92/10/20  16:54:17  ldr
c Changes to V4-0 for SGI (mainly to get nsib and seaice stuff running.)
c 
        subroutine swr89(fsw,hsw,grdflx,ufsw,dfsw,press_in,press2_in,
     &                   coszro_in,taudar_in,rh2o_in,rrco2,ssolar,
     &                   qo3_in,nclds_in,ktopsw_in,kbtmsw_in,cirab_in,
     &                   cirrf_in,cuvrf_in,camt_in)

!$OMP THREADPRIVATE ( /VTEMP/ )

c===>    *********************************************************
c     -sw- radiation code............................
c        inputs:common block radin (except temp);
c        output:fsw,hsw,grdflx(common block swocom).
c===>    *********************************************************
c
c

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'
c          inlte  =  no. levels used for nlte calcs.
c          nnlte  =  index no. of freq. band in nlte calcs.
c          nb,ko2 are shortwave parameters; other quantities are derived
c                    from the above parameters.

C Argument list
      real fsw(imax*lp1)
      real hsw(imax*l)
      real grdflx(imax)
      real ufsw(imax*lp1)
      real dfsw(imax*lp1)
      real press_in(imax*lp1)
      real press2_in(imax*lp1)
      real coszro_in(imax)
      real taudar_in(imax)
      real rh2o_in(imax*l)
      real rrco2
      real ssolar
      real qo3_in(imax*l)
      integer nclds_in(imax)
      integer ktopsw_in(imax*lp1)
      integer kbtmsw_in(imax*lp1)
      real cirab_in(imax*lp1)
      real cirrf_in(imax*lp1)
      real cuvrf_in(imax*lp1)
      real camt_in(imax*lp1)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common /vtemp/ dfn(lp1,imax),ufn(lp1,imax),
     &  ttd(lp1,imax),ttu(lp1,imax),ttdb1(imax*lp1),ttub1(imax*lp1),
     &  tco2(imax*llp2),ud(imax*lp1),ur(imax*lp1),
     &  pp(lp1,imax),dp(imax*lp1),pptop(imax*lp1),dpcld(imax*lp1),
     &  dfntop(imax*nb),cr(imax*lp1),ct(imax*lp1),
     &  alfa(imax*lp1),alfau(imax*lp1),
     &  tdcl1(imax*lp1),tdcl2(imax*lp1),tucl1(imax*lp1),
     &  tucl1i(imax*lp1),tdcl1i(imax*lp1),tdcl2i(imax*lp1),
     &  tclu(imax*lp1),tcld(imax*lp1),
     &  ufnclu(imax*lp1),dfnclu(imax*lp1),
     &  ufntrn(lp1,imax),dfntrn(lp1,imax),
     &  temp1(imax),
     &  refl(imax),secz(imax),rray(imax)
      common /vtemp/ dummy(imax*(2*l*l-13*l-15-nb))
      double precision alfa, alfau, ff, ffco2, cr, tcld, tclu,
     &                 tdcl1, tdcl2, ffo3, pr2, ufnclu, dfnclu, uo3
c---dimensions of variables equivalenced to those in vtemp--
      dimension ff(lp1,imax),ffco2(lp1,imax),ffo3(lp1,imax)
      dimension pr2(imax*lp1)
      dimension du(imax*lp1),duco2(imax*lp1),duo3(imax*lp1)
      dimension uco2(imax*llp2),uo3(imax*llp2)
      dimension ao3(imax*llp2)
c---dimensions of equivalenced variables for l40 model---
c     dimension ao2(imax*llp2)
c     dimension uo2(imax*lp1)
c--dimensions of local data variables---
c---equivalence statements
      equivalence (ff,alfa),(ffco2,alfau)
      equivalence (ffo3,tdcl1),(pr2,tdcl2)
      equivalence (du,tucl1),(duco2,tucl1i)
      equivalence (duo3,tdcl1i)
      equivalence (uo3(1),ufnclu(1))
      equivalence (uco2(1),tclu(1))
      equivalence (ao3(1),ufntrn(1,1))
c
c---equivalence statements for l40 variables---
c     equivalence (uo2,tdcl2i)
c     equivalence (ao2(1),pptop(1))

C Global data blocks

C Local work arrays and variables
      dimension press(imax*lp1),coszro(imax),taudar(imax)
      dimension press2(imax*lp1)
      dimension rh2o(imax*l),qo3(imax*l)
      dimension nclds(imax),ktopsw(imax*lp1),kbtmsw(imax*lp1)
      dimension cirab(imax*lp1),cirrf(imax*lp1),cuvrf(imax*lp1)
      dimension camt(imax*lp1)
c     Introduce the extra array ibn because can't use max(nb,lp1) as dim
      dimension ib(lp1), iq(lp1), ibn(nb)
      dimension iactb(imax),iacte(imax)
      logical lcoszn(imax)

C Local data, functions etc
c     From GFDL set up
      real abcff(nb),pwts(nb),cfco2,cfo3,reflo3,rrayav
      save abcff,    pwts,    cfco2,cfo3,reflo3,rrayav
c---specification of data statements:
c         abcff=absorption coefficients for bands in k-distri-
c     bution. originally given by lacis and hansen, revised by
c     ramaswamy
c         pwts=corresponding weights assigned to bands in the
c     k-distribution
c         reflo3,rrayav= reflection coefficients given by
c     lacis and hansen to account for effects of rayleigh scattering
c     in the visible frequencies (band 1)
c         cfco2,cfo3=conversion factors from gm/cm**2 to cm-atm(stp)
c
c---the following are the coefficients for the 12-band shortwave
c   radiation code, specified by ramaswamy.
      data abcff/2*4.0e-5,.002,.035,.377,1.95,9.40,44.6,190.,
     &           989.,2706.,39011./
      data pwts/.5000,.121416,.0698,.1558,.0631,.0362,.0243,.0158,.0087,
     &          .001467,.002342,.001075/
c---the original 9-band lacis-hansen coefficients are given here; it
c   the user insists on using these values, she must also change
c   the parameter nb from 12 to 9. this parameter is defined in
c   rdparm.h . no other changes are required!
c     data abcff/2*4.0e-5,.002,.035,.377,1.95,9.40,44.6,190./
c     data pwts/.5000,.1470,.0698,.1443,.0584,.0335,.0225,.0158,.0087/
c
      data cfco2,cfo3/508.96,466.64/
      data reflo3/1.9/
      data rrayav/0.144/

C Start code : ----------------------------------------------------------

      ipts=0
      do 100 i=1,imax
        lcoszn(i)=coszro_in(i).gt.zero
        if(lcoszn(i)) ipts=ipts+1
100   continue
      if(ipts.eq.0) then
        do 105 i=1,imax
        grdflx(i)=zero
105     continue
        do 110 i=1,imax*lp1
        ufsw(i)=zero
        dfsw(i)=zero
        fsw(i)=zero
110     continue
        do 115 i=1,imax*l
        hsw(i)=zero
115     continue
        go to 9000
      endif
c     Copy the input data to internal arrays
      do i=1,imax
	  coszro(i) = coszro_in(i)
	  taudar(i) = taudar_in(i)
	  nclds(i) = nclds_in(i)
      end do
      do i=1,imax*l
          rh2o(i) = rh2o_in(i)
	  qo3(i) = qo3_in(i)
      end do
      do i=1,imax*lp1
          press(i) = press_in(i)
          press2(i) = press2_in(i)
	  ktopsw(i) = ktopsw_in(i)
	  kbtmsw(i) = kbtmsw_in(i)
          cirab(i) = cirab_in(i)
	  cirrf(i) = cirrf_in(i)
	  cuvrf(i) = cuvrf_in(i)
          camt(i) = camt_in(i)
      end do
      do k=1,lp1
         ib(k)=ipts*(k-1)
         iq(k)=imax*(k-1)
      end do
      do nx=1,nb
         ibn(nx)=ipts*(nx-1)
      end do
      ipts1=ipts+1
      ipl=ipts*l
      iplp1=ipts*lp1
      ipllp1=ipts*llp1
      if(ipts.lt.imax) then
        iacblk=0
        if(lcoszn(1)) then
          iacblk=1
          iactb(iacblk)=1
        endif
        do 150 i=2,imax
        if(lcoszn(i).neqv.lcoszn(i-1)) then
          if(lcoszn(i-1)) then
            iacte(iacblk)=i-1
          else
            iacblk=iacblk+1
            iactb(iacblk)=i
          endif
        endif
150     continue
        if(lcoszn(imax)) iacte(iacblk)=imax
        if(iactb(1).eq.1) then
          ibl1st=2
          ip=iacte(1)
        else
          ibl1st=1
          ip=0
        endif
        do 347 iblk=ibl1st,iacblk
cdir$ ivdep
*vdir nodep
          do 347 i=iactb(iblk),iacte(iblk)
            ip=ip+1
            nclds(ip)=nclds(i)
            coszro(ip)=coszro(i)
            taudar(ip)=taudar(i)
            rh2o(ip)=rh2o(i)
            qo3(ip)=qo3(i)
            ktopsw(ip)=ktopsw(i)
            kbtmsw(ip)=kbtmsw(i)
            press(ip)=press(i)
            camt(ip)=camt(i)
            cuvrf(ip)=cuvrf(i)
            cirrf(ip)=cirrf(i)
            cirab(ip)=cirab(i)
347     continue
        do 348 k=2,l
          ip=0
          do 348 iblk=1,iacblk
cdir$ ivdep
*vdir nodep
          do 348 i=iactb(iblk),iacte(iblk)
            ip=ip+1
            rh2o(ib(k)+ip)=rh2o(iq(k)+i)
            qo3(ib(k)+ip)=qo3(iq(k)+i)
348     continue
        do 349 k=2,lp1
          ip=0
          do 349 iblk=1,iacblk
cdir$ ivdep
*vdir nodep
          do 349 i=iactb(iblk),iacte(iblk)
            ip=ip+1
            ktopsw(ib(k)+ip)=ktopsw(iq(k)+i)
            kbtmsw(ib(k)+ip)=kbtmsw(iq(k)+i)
            press(ib(k)+ip)=press(iq(k)+i)
            press2(ib(k)+ip)=press2(iq(k)+i)
            camt(ib(k)+ip)=camt(iq(k)+i)
            cuvrf(ib(k)+ip)=cuvrf(iq(k)+i)
            cirrf(ib(k)+ip)=cirrf(iq(k)+i)
            cirab(ib(k)+ip)=cirab(iq(k)+i)
349     continue
      endif
c
c---find the maximum number of clouds in the latitude row
      kclds=nclds(1)
      do 2106 i=2,ipts
      kclds=max(nclds(i),kclds)
2106  continue
c
c   calculate secant of zenith angle (secz),flux pressures(pp),layer
c   width (dp) and pressure scaling factor (pr2).
      do 833 ip=1,ipts
      secz(ip)=h35e1/sqrt(h1224e3*coszro(ip)*coszro(ip)+one)
c     secz(ip)=1./coszro(ip)
      temp1(ip)=1./press(ipl+ip)
833   continue
      do 1 k=1,lp1
      do 1 i=1,ipts
c     pp(k,i)=haf*(press(ib(k)+i)+press(ib(k-1)+i))
      pp(k,i)=press2(ib(k)+i)  ! mrd
1     continue
      do 2 k=1,l
      do 2 i=1,ipts
      dp(ib(k)+i)=pp(k+1,i)-pp(k,i)
      pr2(ib(k)+i)=haf*(pp(k,i)+pp(k+1,i))
2     continue
      do 222 k=1,l
      do 222 ip=1,ipts
      pr2(ib(k)+ip)=pr2(ib(k)+ip)*temp1(ip)
222   continue
c---set up angular factor to be multiplied to the optical factor.
c   above the highest cloud, this is the secant of the zenith
c   angle (modified slightly for refractive effects) (=secz).
c   below the highest cloud, this is diffctr (o3difctr for ozone,
c   in accordance with lacis-hansen parameterization).this factor
c   is used regardless of cloud amount-and for direct and diffuse
c   radiation (this is not a 2 1/2 stream model!)
      do 351 k=1,lp1
      do 351 i=1,ipts
      ff(k,i)=diffctr
      ffco2(k,i)=diffctr
      ffo3(k,i)=o3difctr
 351  continue
      do 131 ip=1,ipts
      jtop=ktopsw(ib(nclds(ip)+1)+ip)
      do 131 k=1,jtop
      ffo3(k,ip)=secz(ip)
      ffco2(k,ip)=secz(ip)
      ff(k,ip)=secz(ip)
131   continue
c     calculate pressure-weighted optical paths for each layer
c     in units of cm-atm. pressure weighting is using pr2.
c     du= value for h2o;duco2 for co2;duo3 for o3.
      do 3 i=1,ipl
      duo3(i)=(ginv*cfo3)*qo3(i)*dp(i)
      duco2(i)=(rrco2*ginv*cfco2)*dp(i)*pr2(i)
      du(i)=ginv*rh2o(i)*dp(i)*pr2(i)
3     continue
c     obtain the optical path from the top of the atmosphere to the
c     flux pressure. angular factors are now included. ud=downward
c     path for h2o,wigth ur the upward path for h2o. corresponding
c     quantities for co2,o3 are udco2/urco2 and udo3/uro3.
      do 834 ip=1,ipts
      ud(ip)=zero
      uco2(ip)=zero
      uo3(ip)=zero
834   continue
      do 835 k=2,lp1
      do 835 i=1,ipts
      ud(ib(k)+i)=ud(ib(k-1)+i)+du(ib(k-1)+i)*ff(k,i)
      uco2(ib(k)+i)=uco2(ib(k-1)+i)+duco2(ib(k-1)+i)* ffco2(k,i)
      uo3(ib(k)+i)=uo3(ib(k-1)+i)+duo3(ib(k-1)+i)*ffo3(k,i)
835   continue
      do 836 ip=1,ipts
      ur(ipl+ip)=ud(ipl+ip)
      uco2(iplp1+ipl+ip)=uco2(ipl+ip)
      uo3(iplp1+ipl+ip)=uo3(ipl+ip)
836   continue
      do 837 k=l,1,-1
      do 837 ip=1,ipts
      ur(ib(k)+ip)=ur(ib(k+1)+ip)+du(ib(k)+ip)*diffctr
      uco2(iplp1+ib(k)+ip)=uco2(iplp1+ib(k+1)+ip)+
     &                     duco2(ib(k)+ip)*diffctr
      uo3(iplp1+ib(k)+ip)=uo3(iplp1+ib(k+1)+ip)+duo3(ib(k)+ip)*reflo3
837   continue
c---for the skyhi model only, obtain the oxygen optical path,using
c   the o3 angular integration.because of the equivalencing,this
c   (and all other) optical path calculations must be complete before
c   transmissivities are computed.
c     do 321 k=2,ko2
c     do 321 i=1,ipts
c     uo2(ib(k)+i)=pp(k,i)*ffo3(k,i)
c321  continue
c---end skyhi model only
c
c     calculate co2 absorptions . they will be used in near infrared
c     bands.since the absorption amount is given (in the formula used
c     below, derived from sasamori) in terms of the total solar flux,
c     and the absorption is only included in the near ir (50 percent
c     of the solar spectrum), the absorptions are multiplied by 2.
c       since code actually requires transmissions, these are the
c     values actually stored in tco2.
      do 2217 i=ipts1,ipllp1
      tco2(i)=1.-two*(h235m3*exp(hp26*log(uco2(i)+h129m2))-h75826m4)
2217  continue
c     now calculate ozone absorptions. these will be used in
c     the visible band.just as in the co2 case, since this band is
c     50 percent of the solar spectrum,the absorptions are multiplied
c     by 2.
      htemp=h1036e2*h1036e2*h1036e2
      do 31 i=ipts1,ipllp1
      ao3(i)=two*uo3(i)*
     &    (h1p082*exp(hmp805*log(one+h1386e2*uo3(i)))+
     &     h658m2/(one+htemp*uo3(i)*uo3(i)*uo3(i))+
     &     h2118m2/(one+uo3(i)*(h42m2+h323m4*uo3(i))))
31    continue
c---for skyhi model only
c     calculate o2 absorptions (in visible band)
c     formula is: abs=1.02e-5*uo2(k)**.3795 for uo2<673.9057
c                 abs=4.51e-6*uo2(k)**.5048 for uo2>673.9057
c     the absorption is constant below 35 km (skyhi level 12)
c     do 844 i=ipts1,ipts*ko2
c     if ((uo2(i)-h67390e2).le.zero) then
c       ao2(i)=1.02e-5*uo2(i)**.3795
c     else
c       ao2(i)=4.51e-6*uo2(i)**.5048
c     endif
c844  continue
c     do 624 k=ko21,llp1
c     do 624 i=1,ipts
c     ao2(ib(k)+i)=ao2(ib(ko2)+i)
c624  continue
c---add o2 absorption to o3 absorption (needed only for skyhi).
c     do 33 i=ipts1,ipllp1
c     ao3(i)=ao3(i)+ao2(i)
c33   continue
c---end skyhi model only
c
c     calculate entering flux at the top for each band(in cgs units)
c
      do 6 nx=1,nb
      do 6 ip=1,ipts
      dfntop(ibn(nx)+ip)=ssolar*h69766e5*coszro(ip)*taudar(ip)*pwts(nx)
6     continue
c
c   start frequency loop (on nx) here
c
c--- band 1 (visible) includes o3,o2 and (negligible) h2o absorption
c
      do 1208 k=1,l
      do 1208 i=1,ipts
      ttdb1(ib(k+1)+i)=exp(hm1ez*min(fifty,abcff(1)*ud(ib(k+1)+i)))
      ttub1(ib(k)+i)=exp(hm1ez*min(fifty,abcff(1)*ur(ib(k)+i)))
      ttd(k+1,i)=ttdb1(ib(k+1)+i)*(1.-ao3(ib(k+1)+i))
      ttu(k,i)=ttub1(ib(k)+i)*(1.-ao3(iplp1+ib(k)+i))
1208  continue
      do 1209 i=1,ipts
      ttd(1,i)=1.
      ttu(lp1,i)=ttd(lp1,i)
1209  continue
c
c
c---execute the lacis-hansen reflectivity parameterization
c---for the visible band
      do 1115 i=1,ipts
        rray(i)=hp219/(one+hp816*coszro(i))
        refl(i)=rray(i)+(one-rray(i))*(one-rrayav)*cuvrf(i)/
     &  (one-cuvrf(i)*rrayav)
1115  continue
c
c***if no clouds at all exist in the row, the calculations are
c   drastically simplified.
      if (kclds.eq.0) then
        do 1117 k=1,lp1
        do 1117 i=1,ipts
          dfn(k,i)=ttd(k,i)
1117    continue
        do 1119 i=1,ipts
         temp1(i)=refl(i)*dfn(lp1,i)/ttu(lp1,i)
1119    continue
        do 1121 k=1,lp1
        do 1121 i=1,ipts
         ufn(k,i)=temp1(i)*ttu(k,i)
1121    continue
      endif
c
c***compute normal case: at least 1 pt has a cloud.
c
      if (kclds.ne.0) then
c
c---the following calculation is done for nl=1 case only
c
c---obtain the pressure at the top,bottom and the thickness of
c   "thick" clouds (those at least 2 layers thick). this is used
c   later is obtaining fluxes inside the thick clouds, for all
c   frequency bands.
      do 4226 kk=1,kclds
      do 4226 i=1,ipts
        if ((kbtmsw(ib(kk+1)+i)-1).gt.ktopsw(ib(kk+1)+i)) then
           pptop(ib(kk)+i)=pp(ktopsw(ib(kk+1)+i),i)
           dpcld(ib(kk)+i)=one/(pptop(ib(kk)+i)-
     &                     pp(kbtmsw(ib(kk+1)+i),i))
        endif
4226  continue
c
      do 832 i=1,ipts
c---the first cloud is the ground; its properties are given
c   by refl (the transmission (0) is irrelevant for now!).
      cr(i)=refl(i)
832   continue
c***obtain cloud reflection and transmission coefficients for
c   remaining clouds (if any) in the visible band
c---the maximum no of clouds in the row (kclds) is used. this creates
c   extra work (may be removed in a subsequent update).
      do 4201 i=ipts1,ipts*(kclds+1)
          cr(i)=cuvrf(i)*camt(i)
          ct(i)=one-cr(i)
4201  continue
c
c
c***for execution of the cloud loop, it is necessary to separate out
c   transmission fctns at the top and bottom of the clouds, for
c   each band n. the required quantities are:
c      ttd(i,ktopsw(i,k),nl)  k runs from 1 to nclds(i)+1:
c      ttd(i,kbtmsw(i,k),nl)  k runs from 1 to nclds(i)+1:
c      ttu(i,ktopsw(i,k),nl)  k runs from 1 to nclds(i)+1:
c      and inverses of the above. the above quantities are stored
c      in tdcl1,tdcl2,tucl1,tdcl1i,tdcl2i,tucli,respectively, as
c      they have multiple use in the pgm.
c
c---for first cloud layer (ground) tdcl1,tucl1 are known:
      do 2321 i=1,ipts
        tdcl1(i)=ttd(lp1,i)
        tucl1(i)=ttu(lp1,i)
        tdcl1i(i)=one/tdcl1(i)
        tucl1i(i)=tdcl1i(i)
        tdcl2(i)=tdcl1(i)
2321  continue
c---for remaining ones,use gathers:
      do 2323 kk=2,kclds+1
      do 2323 i=1,ipts
        tdcl1(ib(kk)+i)=ttd(ktopsw(ib(kk)+i),i)
        tucl1(ib(kk)+i)=ttu(ktopsw(ib(kk)+i),i)
        tdcl2(ib(kk)+i)=ttd(kbtmsw(ib(kk)+i),i)
2323  continue
c---compute inverses
      do 2325 i=ipts1,ipts*(kclds+1)
         tdcl1i(i)=one/tdcl1(i)
         tucl1i(i)=one/tucl1(i)
2325  continue
c---compute the transmissivity from the top of cloud (k+1) to the
c   top of cloud (k). the cloud transmission (ct) is included. this
c   quantity is called tclu (index k). also, obtain the transmissivity
c   from the bottom of cloud (k+1) to the top of cloud (k)(a path
c   entirely outside clouds). this quantity is called tcld (index k).
      do 2327 i=1,ipts*kclds
         tclu(i)=tdcl1(i)*tdcl1i(ipts+i)*ct(ipts+i)
         tcld(i)=tdcl1(i)/tdcl2(ipts+i)
2327  continue
c
c***the following is the recursion relation for alfa: the reflection
c   coefficient for a system including the cloud in question and the
c   flux coming out of the cloud system including all clouds below
c   the cloud in question.
c---alfau is alfa without the reflection of the cloud in question
      do 2331 i=1,ipts
         alfa(i)=cr(i)
         alfau(i)=zero
2331  continue
      do 2332 kk=2,kclds+1
      do 2334 i=1,ipts
c---again,excessive calculations-may be changed later!
        alfau(ib(kk)+i)=tclu(ib(kk-1)+i)*tclu(ib(kk-1)+i)*
     &                  alfa(ib(kk-1)+i)/(1.-tcld(ib(kk-1)+i)*
     &                  tcld(ib(kk-1)+i)*alfa(ib(kk-1)+i)*cr(ib(kk)+i))
        alfa(ib(kk)+i)=alfau(ib(kk)+i)+cr(ib(kk)+i)
2334  continue
2332  continue
c
c     calculate ufn at cloud tops and dfn at cloud bottoms
c---note that ufnclu(i,kclds+1) gives the upward flux at the top
c   of the highest real cloud (if nclds(i)=kclds). it gives the flux
c   at the top of the atmosphere if nclds(i) < kclds. in the first
c   case, tdcl1 equals the transmission fctn to the top of the
c   highest cloud, as we want. in the second case, tdcl1=1, so ufnclu
c   equals alfa. this is also correct.
      do 2370 i=ipts*kclds+1,ipts*(kclds+1)
      ufnclu(i)=alfa(i)*tdcl1(i)
      dfnclu(i)=tdcl1(i)
2370  continue
c---this calculation is the reverse of the recursion relation used
c  above
      do 2371 kk=kclds,1,-1
      do 2371 i=1,ipts
      ufnclu(ib(kk)+i)=ufnclu(ib(kk+1)+i)*alfau(ib(kk+1)+i)/
     &                 (alfa(ib(kk+1)+i)*tclu(ib(kk)+i))
      dfnclu(ib(kk)+i)=ufnclu(ib(kk)+i)/alfa(ib(kk)+i)
2371  continue
c     now obtain dfn and ufn for levels between the clouds
      do 2381 k=1,kclds+1
      do 2381 i=1,ipts
        ufntrn(k,i)=ufnclu(ib(k)+i)*tucl1i(ib(k)+i)
        dfntrn(k,i)=dfnclu(ib(k)+i)*tdcl1i(ib(k)+i)
2381  continue
c---case of kk=1( from the ground to the bottom of the lowest cloud)
      do 2390 i=1,ipts
      j2=kbtmsw(ipts+i)
      ufnsav=ufntrn(1,i)
      dfnsav=dfntrn(1,i)
      do 2392 k=j2,lp1
        ufn(k,i)=ufnsav*ttu(k,i)
        dfn(k,i)=dfnsav*ttd(k,i)
2392  continue
2390  continue
c---remaining levels (if any!)
      do 2391 kk=2,kclds+1
      do 2394 i=1,ipts
      j1=ktopsw(ib(kk)+i)
      j2=kbtmsw(ib(kk+1)+i)
      if (j1.eq.1) go to 2394
      ufnsav=ufntrn(kk,i)
      dfnsav=dfntrn(kk,i)
      do 2393 k=j2,j1
       ufn(k,i)=ufnsav*ttu(k,i)
       dfn(k,i)=dfnsav*ttd(k,i)
2393  continue
c---for the thick clouds, the flux divergence through the cloud
c   layer is assumed to be constant. the flux derivative is given by
c    tempf (for the upward flux) and tempg (for the downward flux).
      j3=kbtmsw(ib(kk)+i)
      if ((j3-j1).gt.1) then
        ufnsav=ufnclu(ib(kk)+i)
        dfnsav=dfnclu(ib(kk)+i)
        ppsav=pptop(ib(kk-1)+i)
        tempf=(ufnsav-ufn(j3,i))*dpcld(ib(kk-1)+i)
        tempg=(dfnsav-dfn(j3,i))*dpcld(ib(kk-1)+i)
        do 2397 k=j1+1,j3-1
         ufn(k,i)=ufnsav+tempf*(pp(k,i)-ppsav)
         dfn(k,i)=dfnsav+tempg*(pp(k,i)-ppsav)
2397    continue
      endif
2394  continue
2391  continue
c
      endif
c
c     scale visible band fluxes by solar flux at the top of the
c     atmosphere (dfntop(i))
c     dfsw/ufsw will be the fluxes, summed over all bands
      do 1615 k=1,lp1
      do 1615 i=1,ipts
      dfsw(ib(k)+i)=dfn(k,i)*dfntop(i)
      ufsw(ib(k)+i)=ufn(k,i)*dfntop(i)
1615  continue
c
c
c---now obtain fluxes for the near ir bands. the methods are the same
c   as for the visible band, except that the reflection and
c   transmission coefficients (obtained below) are different, as
c   rayleigh scattering need not be considered.
      do 4217 i=1,ipts*(kclds+1)
          cr(i)=cirrf(i)*camt(i)
          ct(i)=one-camt(i)*
     &                  (cirrf(i)+cirab(i))
4217  continue
c
      do 3301 nx=2,nb
c
      if (nx.eq.2) then
c   the water vapor transmission function for band 2 is equal to
c   that of band 1 (saved as ttdb1,ttub1)
        do 2207 k=1,l
        do 2207 i=1,ipts
        ttd(k+1,i)=ttdb1(ib(k+1)+i)*tco2(ib(k+1)+i)
        ttu(k,i)=ttub1(ib(k)+i)*tco2(iplp1+ib(k)+i)
2207  continue
        do 2208 i=1,ipts
        ttd(1,i)=1.
        ttu(lp1,i)=ttd(lp1,i)
2208    continue
      else
c   calculate water vapor transmission functions for near infrared
c   bands. include co2 transmission (tdco2/tuco2), which
c   is the same for all infrared bands.
        do 2209 k=1,l
        do 2209 i=1,ipts
        ttd(k+1,i)=exp(hm1ez*min(fifty,abcff(nx)*ud(ib(k+1)+i)))*
     &             tco2(ib(k+1)+i)
        ttu(k,i)=exp(hm1ez*min(fifty,abcff(nx)*ur(ib(k)+i)))*
     &           tco2(iplp1+ib(k)+i)
2209    continue
c---at this point,include ttd(1),ttu(lp1), noting that ttd(1)=1 for
c   all bands, and that ttu(lp1)=ttd(lp1) for all bands.
        do 45 ip=1,ipts
        ttu(lp1,ip)=ttd(lp1,ip)
        ttd(1,ip)=one
45      continue
      endif
c
c***again, if no clouds at all exist in the row, the calculations are
c   drastically simplified.
      if (kclds.eq.0) then
        do 2118 k=1,lp1
        do 2118 i=1,ipts
          dfn(k,i)=ttd(k,i)
2118    continue
        do 2120 i=1,ipts
         temp1(i)=cirrf(i)*dfn(lp1,i)/ttu(lp1,i)
2120    continue
        do 2122 k=1,lp1
        do 2122 i=1,ipts
         ufn(k,i)=temp1(i)*ttu(k,i)
2122    continue
      endif
c
      if (kclds.ne.0) then
c
c***for execution of the cloud loop, it is necessary to separate out
c   transmission fctns at the top and bottom of the clouds, for
c   each band n. the required quantities are:
c      ttd(i,ktopsw(i,k),nl)  k runs from 1 to nclds(i)+1:
c      ttd(i,kbtmsw(i,k),nl)  k runs from 2 to nclds(i)+1:
c      ttu(i,ktopsw(i,k),nl)  k runs from 1 to nclds(i)+1:
c      and inverses of the above. the above quantities are stored
c      in tdcl1,tdcl2,tucl1,tdcl1i,tdcl2i,tucli,respectively, as
c      they have multiple use in the pgm.
c
c---for first cloud layer (ground) tdcl1,tucl1 are known:
      do 3321 i=1,ipts
        tdcl1(i)=ttd(lp1,i)
        tucl1(i)=ttu(lp1,i)
        tdcl1i(i)=one/tdcl1(i)
        tucl1i(i)=tdcl1i(i)
        tdcl2(i)=tdcl1(i)
3321  continue
c---for remaining ones,use gathers:
      do 3323 kk=2,kclds+1
      do 3323 i=1,ipts
        tdcl1(ib(kk)+i)=ttd(ktopsw(ib(kk)+i),i)
        tucl1(ib(kk)+i)=ttu(ktopsw(ib(kk)+i),i)
        tdcl2(ib(kk)+i)=ttd(kbtmsw(ib(kk)+i),i)
3323  continue
c---compute inverses
      do 3325 i=ipts1,ipts*(kclds+1)
         tdcl1i(i)=one/tdcl1(i)
         tucl1i(i)=one/tucl1(i)
3325  continue
      do 3327 i=1,ipts*kclds
         tclu(i)=tdcl1(i)*tdcl1i(ipts+i)*ct(ipts+i)
         tcld(i)=tdcl1(i)/tdcl2(ipts+i)
3327  continue
c
c***the following is the recursion relation for alfa: the reflection
c   coefficient for a system including the cloud in question and the
c   flux coming out of the cloud system including all clouds below
c   the cloud in question.
      do 3331 i=1,ipts
         alfa(i)=cr(i)
         alfau(i)=zero
3331  continue
      do 3332 kk=2,kclds+1
      do 3334 i=1,ipts
c---again,excessive calculations-may be changed later!
        alfau(ib(kk)+i)=tclu(ib(kk-1)+i)*tclu(ib(kk-1)+i)*
     &                  alfa(ib(kk-1)+i)/(1.-tcld(ib(kk-1)+i)*
     &                  tcld(ib(kk-1)+i)*alfa(ib(kk-1)+i)*cr(ib(kk)+i))
        alfa(ib(kk)+i)=alfau(ib(kk)+i)+cr(ib(kk)+i)
3334  continue
3332  continue
c
c     calculate ufn at cloud tops and dfn at cloud bottoms
c---note that ufnclu(i,kclds+1) gives the upward flux at the top
c   of the highest real cloud (if nclds(i)=kclds). it gives the flux
c   at the top of the atmosphere if nclds(i) < kclds. it the first
c   case, tdcl1 equals the transmission fctn to the top of the
c   highest cloud, as we want. in the second case, tdcl1=1, so ufnclu
c   equals alfa. this is also correct.
      do 3370 i=ipts*kclds+1,ipts*(kclds+1)
      ufnclu(i)=alfa(i)*tdcl1(i)
      dfnclu(i)=tdcl1(i)
3370  continue
      do 3371 kk=kclds,1,-1
      do 3371 i=1,ipts
      if ( cr(ib(kk+1)+i).gt.0. ) then  ! test added mrd 3/9/90
          ufnclu(ib(kk)+i)=ufnclu(ib(kk+1)+i)*alfau(ib(kk+1)+i)/
     &                 (alfa(ib(kk+1)+i)*tclu(ib(kk)+i))
      else ! in this case the alfa's cancel
          ufnclu(ib(kk)+i)=ufnclu(ib(kk+1)+i)/tclu(ib(kk)+i)
      end if
      dfnclu(ib(kk)+i)=ufnclu(ib(kk)+i)/alfa(ib(kk)+i)
3371  continue
c     now obtain dfn and ufn for levels between the clouds
      do 3381 k=1,kclds+1
      do 3381 i=1,ipts
        ufntrn(k,i)=ufnclu(ib(k)+i)*tucl1i(ib(k)+i)
        dfntrn(k,i)=dfnclu(ib(k)+i)*tdcl1i(ib(k)+i)
3381  continue
c---case of kk=1
      do 3390 i=1,ipts
      j2=kbtmsw(ipts+i)
      ufnsav=ufntrn(1,i)
      dfnsav=dfntrn(1,i)
      do 3392 k=j2,lp1
        ufn(k,i)=ufnsav*ttu(k,i)
        dfn(k,i)=dfnsav*ttd(k,i)
3392  continue
3390  continue
      do 3391 kk=2,kclds+1
      do 3394 i=1,ipts
      j1=ktopsw(ib(kk)+i)
      j2=kbtmsw(ib(kk+1)+i)
      if (j1.eq.1) go to 3394
      ufnsav=ufntrn(kk,i)
      dfnsav=dfntrn(kk,i)
      do 3393 k=j2,j1
       ufn(k,i)=ufnsav*ttu(k,i)
       dfn(k,i)=dfnsav*ttd(k,i)
3393  continue
      j3=kbtmsw(ib(kk)+i)
      if ((j3-j1).gt.1) then
        ufnsav=ufnclu(ib(kk)+i)
        dfnsav=dfnclu(ib(kk)+i)
        ppsav=pptop(ib(kk-1)+i)
        tempf=(ufnsav-ufn(j3,i))*dpcld(ib(kk-1)+i)
        tempg=(dfnsav-dfn(j3,i))*dpcld(ib(kk-1)+i)
        do 3397 k=j1+1,j3-1
         ufn(k,i)=ufnsav+tempf*(pp(k,i)-ppsav)
         dfn(k,i)=dfnsav+tempg*(pp(k,i)-ppsav)
3397    continue
      endif
3394  continue
3391  continue
c
      endif
c
c     scale the previously computed fluxes by the flux at the top
c     and sum over bands
c
      do 860 k=1,lp1
      do 860 ip=1,ipts
      dfsw(ib(k)+ip)=dfsw(ib(k)+ip)+dfn(k,ip)*dfntop(ibn(nx)+ip)
      ufsw(ib(k)+ip)=ufsw(ib(k)+ip)+ufn(k,ip)*dfntop(ibn(nx)+ip)
860   continue
c
3301  continue
c
c---end of frequency loop (over nx)
c
      if(ipts.lt.imax) then
        do 505 k=lp1,2,-1
          ip=ipts1
          do 505 iblk=iacblk,1,-1
cdir$ ivdep
*vdir nodep
          do 505 i=iacte(iblk),iactb(iblk),-1
            ip=ip-1
            ufsw(iq(k)+i)=ufsw(ib(k)+ip)
            dfsw(iq(k)+i)=dfsw(ib(k)+ip)
            dp(iq(k)+i)=dp(ib(k)+ip)
505     continue
          ip=ipts1
          do 506 iblk=iacblk,ibl1st,-1
cdir$ ivdep
*vdir nodep
          do 506 i=iacte(iblk),iactb(iblk),-1
            ip=ip-1
            ufsw(i)=ufsw(ip)
            dfsw(i)=dfsw(ip)
            dp(i)=dp(ip)
506     continue
c
        if(iactb(1).gt.1) then
          do 509 k=1,lp1
            do 509 i=1,iactb(1)-1
              ufsw(iq(k)+i)=zero
              dfsw(iq(k)+i)=zero
              dp(iq(k)+i)=one
509       continue
        endif
        do 510 k=1,lp1
          do 510 iblk=2,iacblk
          do 510 i=iacte(iblk-1)+1,iactb(iblk)-1
            ufsw(iq(k)+i)=zero
            dfsw(iq(k)+i)=zero
            dp(iq(k)+i)=one
510     continue
        if(iacte(iacblk).lt.imax) then
          do 511 k=1,lp1
            do 511 i=iacte(iacblk)+1,imax
              ufsw(iq(k)+i)=zero
              dfsw(iq(k)+i)=zero
              dp(iq(k)+i)=one
511       continue
        endif
      endif
      do 15 i=1,imax
      grdflx(i)=dfsw(imax*l+i)-ufsw(imax*l+i)
15    continue
c
c---obtain the net flux and the shortwave heating rate
      do 16 i=1,imax*lp1
      fsw(i)=ufsw(i)-dfsw(i)
16    continue
      do 17 i=1,imax*l
c     hsw(i)=radcon*(fsw(imax+i)-fsw(i))/dp(i)
      hsw(i)=(fsw(imax+i)-fsw(i))
17    continue
9000  continue
      return
      end
