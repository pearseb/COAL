c Removing the declaration of the COMMON block /GRADNS/ to a header file.
c SJP 2009/04/25
c
c $Log: tracera.f,v $
c Revision 1.10  1996/10/24 01:03:20  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.9  1996/03/21  03:19:11  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.8  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.7  1993/12/17  15:34:20  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.6  93/03/15  15:45:18  ldr
c HBG's changes to combine uvpgd and sltjmcg common blocks.
c 
c Revision 1.5  92/12/10  09:55:53  ldr
c Replace np's with nlp's for compatibility with ogcm.
c 
c Revision 1.4  92/12/09  14:44:45  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/08/06  16:41:03  ldr
c Change name to tracera (so as not to clash with ogcm).
c 
c Revision 1.2  91/03/13  13:01:13  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:20  ldr
c Initial release V3-0
c 
      subroutine tracera
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      parameter (nps=63*48)
      include 'GAUSL.f'
      include 'UVPGD.f'
      include 'TIMEX.f'
      include 'GRADNS.f'
      integer lonc(nps),latc(nps),ksc(nps)
      real lonpos,latpos
      common/positn/lonpos(nps),latpos(nps),sigpos(nps),npts
      include 'HYBARR.f'
C**** Tracer routine - determines box containing each particle
C****                - determines the U,V,SIGMA-DOT velocity at each
C****                  corner of the box
C****                - determines the weighted mean velocity of the 
C****                  particle, and then adjusts the position
C****
c---- To determine upper,SW corner of box containing particle
c     lonc is longitude of corner
c     latc is latitude of corner
c     ksc is vertical (down) level of corner
c       lonc=1 to 64 (1=Greenwich)
c       latc=1 to 56 (56 Gaussian lats SP to NP)
c       ksc=1 to 9 : 1=top level , sigma=0.021 (using inverted model notation)
c                    9=bottom level , sigma=0.979
c----  (lonc,latc,ksc) = UPPER, SW corner of box = 1.
c
c                       4..........3
c                      / :       /:
c                   1....:.....2  :
c                    :   : *  :   : 
c                    :  8:....:...:7
c                    : /      :  /
c                   5:........:6
c
CZZZZ 1 particle/step
      numsteps=npts
CZZZZ
      tpi=2.0*pi
      do 10 i=1,numsteps
   10 lonc(i)=int((lonpos(i)/tpi)*lon)+1
    
      do 12 i=1,numsteps
      do 14 k=1,lat2-2
      if(latpos(i).le.glats(k+1))go to 16
   14 continue
      k=lat2-1
   16 latc(i)=k
   12 continue

      do 18 i=1,numsteps
      do 20 k=1,nl-1
      if(sigpos(i).le.sigi(k+1))go to 22
   20 continue
   22 ksc(i)=k
   18 continue
c----
c---- get the velocities of the 8 corners for each particle
c----
      do 30 i=1,numsteps
      lonw=lonc(i)
      lone=lonw+1
      if(lone.eq.(lon+1))lone=1
      lats=latc(i)
      latn=lats+1
      ktop=ksc(i)
      kbot=ktop+1
      lgs=lats
      if(lgs.gt.lat)lgs=lat2+1-lgs
      lgn=latn
      if(lgn.gt.lat)lgn=lat2+1-lgn
c---- velocities ugd,vgd have factor (P*).cos(lat)/erad
c---- do not remove 'rad' factor since arc displacements need this
      siais=1.0/sia(lgs)
      siain=1.0/sia(lgn)
      npkt=nlp-ktop
      npkb=nlp-kbot
c---- pressure factor at corners (may be hybrid coordinate)
      pwst=dadnf(npkt)+dbdnf(npkt)*pgd(lonw,lats)
      pest=dadnf(npkt)+dbdnf(npkt)*pgd(lone,lats)
      pwsb=dadnf(npkb)+dbdnf(npkb)*pgd(lonw,lats)
      pesb=dadnf(npkb)+dbdnf(npkb)*pgd(lone,lats)
      pwnt=dadnf(npkt)+dbdnf(npkt)*pgd(lonw,latn)
      pent=dadnf(npkt)+dbdnf(npkt)*pgd(lone,latn)
      pwnb=dadnf(npkb)+dbdnf(npkb)*pgd(lonw,latn)
      penb=dadnf(npkb)+dbdnf(npkb)*pgd(lone,latn)
c---- u velocities
      u1= ugd(lonw,lats,npkt)*siais/pwst
      u2= ugd(lone,lats,npkt)*siais/pest
      u3= ugd(lonw,lats,npkb)*siais/pwsb
      u4= ugd(lone,lats,npkb)*siais/pesb
      u5= ugd(lonw,latn,npkt)*siain/pwnt
      u6= ugd(lone,latn,npkt)*siain/pent
      u7= ugd(lonw,latn,npkb)*siain/pwnb
      u8= ugd(lone,latn,npkb)*siain/penb
c---- v velocities
      v1= vgd(lonw,lats,npkt)*siais/pwst
      v2= vgd(lone,lats,npkt)*siais/pest
      v3= vgd(lonw,lats,npkb)*siais/pwsb
      v4= vgd(lone,lats,npkb)*siais/pesb
      v5= vgd(lonw,latn,npkt)*siain/pwnt
      v6= vgd(lone,latn,npkt)*siain/pent
      v7= vgd(lonw,latn,npkb)*siain/pwnb
      v8= vgd(lone,latn,npkb)*siain/penb
c---- sigma-dot velocities
      sd1=sdot(lonw,npkt,lats)
      sd2=sdot(lone,npkt,lats)
      sd3=sdot(lonw,npkb,lats)
      sd4=sdot(lone,npkb,lats)
      sd5=sdot(lonw,npkt,latn)
      sd6=sdot(lone,npkt,latn)
      sd7=sdot(lonw,npkb,latn)
      sd8=sdot(lone,npkb,latn)
c----
c---- knowing the particle position, and corner
c---- calculate the weighting factors
c---- and then the particle velocity.
c---- Then compute the change in position.
c----
      x=(lonpos(i)/tpi)*lon+1.0
      xw=lonc(i)
      xe=xw+1
      y=min(max(latpos(i),glats(1)),glats(lat2))
      ys=glats(latc(i))
      yn=glats(latc(i)+1)
      z=min(sigpos(i),sigi(nl))
      ztop=sigi(ksc(i))
      zbot=sigi(ksc(i)+1)
      ax1=abs(x-xw)
      ax2=abs(xe-x)
      ay1=abs(y-ys)
      ay2=abs(yn-y)
      az1=abs(z-ztop)
      az2=abs(zbot-z)
      vol=abs((xe-xw)*(yn-ys)*(zbot-ztop))
      f1=ax2*ay2*az2/vol
      f2=ax1*ay2*az2/vol
      f3=ax1*ay1*az2/vol
      f4=ax2*ay1*az2/vol
      f5=ax2*ay2*az1/vol
      f6=ax1*ay2*az1/vol
      f7=ax1*ay1*az1/vol
      f8=ax2*ay1*az1/vol
      sdint=sd1*f1+sd2*f2+sd3*f3+sd4*f4
     & +sd5*f5+sd6*f6+sd7*f7+sd8*f8
c     factor the sdot value for points below lowest level
       factr=(1.0-max(sigi(nl),sigpos(i)))/(1.0-sigi(nl))
       sdint=sdint*factr
      uint=u1*f1+u2*f2+u3*f3+u4*f4+u5*f5+u6*f6+u7*f7+u8*f8
      vint=v1*f1+v2*f2+v3*f3+v4*f4+v5*f5+v6*f6+v7*f7+v8*f8
      lonpos(i)=lonpos(i)+uint*dt/cos(latpos(i))
      latpos(i)=latpos(i)+vint*dt
      sigpos(i)=max(min(sigpos(i)+sdint*dt,0.999999),sigi(1))
   30 continue

c---- check polar boundaries before EW boundaries
      pi2=pi/2.0
c     kns=0
c     knn=0
      do 84 i=1,numsteps
      if(latpos(i).ge.-pi2)go to 82
      latpos(i)=-pi-latpos(i)
      lonpos(i)=pi+lonpos(i)
c     kns=kns+1
   82 if(latpos(i).le.pi2)go to 84
      latpos(i)=pi-latpos(i)
      lonpos(i)=pi+lonpos(i)
c     knn=knn+1
   84 continue
c     if((knn+kns).gt.0)
c    * write(6,184)knn,kns
c 184 format(1x,' Num crossing NP=',i4,' Num crossing SP=',i4)

c---- check E-W boudaries
      do 74 i=1,numsteps
   70 if(lonpos(i).ge.0.0)go to 72
      lonpos(i)=tpi+lonpos(i)
      go to 70
   72 if(lonpos(i).lt.tpi)go to 74
      lonpos(i)=lonpos(i)-tpi
      go to 72 
   74 continue

      return
      end
