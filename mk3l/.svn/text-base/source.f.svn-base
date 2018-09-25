c $Log: source.f,v $
c Revision 1.8  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.7  1996/10/24 01:03:15  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/03/21  03:19:04  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.5  1994/08/04  16:56:36  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c
c Revision 1.4  93/08/10  12:13:56  ldr
c Re-introduce (modified) memory saving trick.
c 
c Revision 1.3  93/07/02  12:52:52  ldr
c Remove memory saving trick (compiler problems on Cray).
c 
c Revision 1.1  93/01/26  16:24:55  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:  from common/cnsta in CNSTA.f
c                 dsk - sigma thicknesses (1 to nl) 
c
c             from common/timex in TIMEX.f
c                 dt - time step in seconds
c                 nsteps - time step counter
c
c             from common/traceblk in TRACEBLK.f
c                 source2,source3 - tracer sources
c
c     In/Out: from common/traceblk in TRACEBLK.f
c                 con - concentration of atmospheric tracers at t
c                 nslfirst - if 1, starting from data set with no backward
c                            step (mainly for use by IGW)
c 
c
      subroutine source

c Sample version of routine to set up sources for Semi-Lagrangian
c tracer modelling

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'TRACEBLK.f'
      include 'UVPGD.f'
      include 'CNSTA.f'
      include 'TIMEX.f'

      if(nx.eq.0)then
        print*,'ERROR: Need to set NX=1 in TRACEBLK.f'
        stop
      endif

c  add sources 2 and 3
c  divide source in kg by kg of atmospheric layer

      if(nslfirst.eq.0) rdsk1=2.*dt*grav/dsk(1)/100.
      if(nslfirst.eq.1) rdsk1=dt*grav/dsk(1)/100.
      if(nslfirst.eq.1.or.nsteps.eq.0) rmass=rdsk1/1000.
c     print *,' in source nslfirst,nsteps ',nslfirst,nsteps,pstar(1,1)

      do ins = 1,2 
        do lg=1,lat
          lgns=lg*(ins-1)+(lat2p-lg)*(2-ins)
          do mg=1,lon
            if(nslfirst.eq.0.and.nsteps.ne.0) rmass=rdsk1/pstar(mg,lgns)
            con(mg,lgns,1,2)=con(mg,lgns,1,2)+source2(mg,lgns)*rmass
            con(mg,lgns,1,3)=con(mg,lgns,1,3)+source3(mg,lgns)*rmass
          enddo
        enddo
      enddo
      nslfirst=0

      return
      end
      
