c $Log: lgndre.f,v $
c Revision 1.5  1996/06/13 02:07:07  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.4  1993/12/17  15:33:00  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  92/12/09  14:43:47  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.2  92/04/22  10:12:53  ldr
c Generalize for 64/32 bit machines.
c 
c Revision 1.1  91/02/22  16:37:36  ldr
c Initial release V3-0
c 
      subroutine lgndre(p,coas,sias,deltas,ir,irmax2)
C**** COMPUTES LEGENDRE POLYS AT A GIVEN LATITUDE
c
c     in this version of lgndre all internal calculations are done in
c     double precision.
c
      implicit  double precision(a-h,o-z)
      dimension p(*)
      real coas,sias,deltas
      coa  =coas
      sia  =sias
      delta=deltas
      irpp=ir+2
      sqr2=sqrt(2.0)
C
      theta=delta
      c1=sqr2
      p(1)=1.0/c1
      do 20 nl=1,irpp
      fn=nl
      fn2=2.*fn
      fn2sq=fn2*fn2
      c1=c1*sqrt(1.0-1.0/fn2sq)
      c3=c1/sqrt(fn*(fn+1.0))
      ang=fn*theta
      s1=0.0
      s2=0.0
      c4=1.0
      c5=fn
      a=-1.0
      b=0.0
      n1=nl+1
      do 27 kk=1,n1,2
      k=kk-1
      s2=s2+c5*sin(ang)*c4
      if (k.eq.nl) c4=0.5*c4
      s1=s1+c4*cos(ang)
      a=a+2.0
      b=b+1.0
      fk=k
      ang=theta*(fn-fk-2.0)
      c4=(a*(fn-b+1.0)/(b*(fn2-a)))*c4
      c5=c5-2.0
   27 continue
      if(nl-irpp)23,21,20
   23 p(nl+1)=s1*c1
   21 p(nl+irmax2)=s2*c3
   20 continue
C
C     ***** 20 HAS SET UP LEGENDRE POLYNOMIALS FOR M=0 AND M=1
C
      if(ir   .eq.2)go to 100
      do 40 m=2,ir
      fm=m
      fm1=fm-1.0
      fm2=fm-2.0
      fm3=fm-3.0
      mm1=m-1
      m1=m+1
      c6=sqrt((2.0*fm+1.0)/(2.0*fm))
      p(irmax2*m+1)=c6*sia*p(irmax2*mm1+1)
      mpir=m+ir   +1
      do 40 l=m1,mpir
      fn=l
      c7=(fn*2.0+1.0)/(fn*2.0-1.0)
      c8=(fm1+fn)/((fm+fn)*(fm2+fn))
      c=sqrt((fn*2.0+1.0)/(fn*2.0-3.0)*c8*(fm3+fn))
      d=sqrt(c7*c8*(fn-fm1))
      e=sqrt(c7*(fn-fm)/(fn+fm))
      lm=irmax2*m+l-m+1
      lmm2=irmax2*(m-2)+l-m+3
      lm1mm2=lmm2-1
      lm2mm2=lm1mm2-1
      lm1m=lm-1
      if(l-mpir)43,42,40
   43 p(lm)=c*p(lm2mm2)-d*p(lm1mm2)*coa+e*p(lm1m) *coa
      go to 40
C
C     ***** BELOUSOV EQUATION 11
C
   42 a=sqrt((fn*fn-0.25)/(fn*fn-fm*fm))
      b=sqrt((2.0*fn+1.0)*(fn-fm-1.0)*(fn+fm1)/
     &((2.0*fn-3.0)*(fn-fm)*(fn+fm)))
      lm2m=lm1m-1
      p(lm)=2.0*a*coa*p(lm1m)-b*p(lm2m)
   40 continue
  100 continue
      return
      end
