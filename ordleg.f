c $Log: ordleg.f,v $
c Revision 1.3  1992/12/09 14:44:05  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c Revision 1.2  92/04/22  10:12:55  ldr
c Generalize for 64/32 bit machines.
c 
c Revision 1.1  91/02/22  16:37:46  ldr
c Initial release V3-0
c 
      subroutine ordleg(sx,coa,ir)
C
      implicit double precision (a-h,o-z)
      irpp = ir + 1
      irppm = irpp - 1
      delta = acos(coa)
      sqr2= sqrt(2.0  )
C
      theta=delta
      c1=sqr2
      do 20 nl=1,irppm
      fn=nl
      fn2=fn+fn
      fn2sq=fn2*fn2
      c1=c1*sqrt(1.0-1.0/fn2sq)
   20 continue
C
      nl=irppm
      ang=fn*theta
      s1=0.0
      c4=1.0
      a=-1.0
      b=0.0
      n1=nl+1
      do 27 kk=1,n1,2
      k=kk-1
      if (k.eq.nl) c4=0.5*c4
      s1=s1+c4* cos(ang)
      a=a+2.0
      b=b+1.0
      fk=k
      ang=theta*(fn-fk-2.0)
      c4=(a*(fn-b+1.0)/(b*(fn2-a)))*c4
   27 continue
      sx=s1*c1
C
      return
      end
