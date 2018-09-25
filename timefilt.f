c $Log: timefilt.f,v $
c Revision 1.5  1996/06/13 02:08:39  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.4  1996/03/21  03:19:10  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.3  1993/12/17  15:34:11  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  92/12/09  14:44:42  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  92/06/16  11:55:06  ldr
c Initial revision
c 
      subroutine timefilt (fuh, fvh, fu, fv)

c Transfers h-grid fields fuh,fvh to the u-grid, and applies this input
c to fu,fv using a running-mean time smoother.

c fuh = eastward  field stress on h-grid (supplied)
c fvh = northward wind stress on h-grid (supplied)
c fu  = eastward  time-smoothed wind stress on u-grid (returned)
c fv  = northward time-smoothed wind stress on u-grid (returned)
c dt  = dynamic seaice time step, seconds
c first = logical flag indicating first call this run

      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'COMDICE.f' !

      dimension
     &  fuh(plon,plat),             fvh(plon,plat),
     &  fu(0:plon+1,plat+1),        fv(0:plon+1,plat+1)



c If first call, transfer straight to u-grid. If not, transfer to 
c temporary work arrays and apply the time smoother by linear weighting 
c with previous values. The weight of the nl-th timestep in the past
c turns out to be wt * (1-wt)**n.


        call ltou (fuh, fu)
        call ltou (fvh, fv)

C comment out for present as the input fields are not being upated
c     else

c       call ltou (fuh, worku)
c       call ltou (fvh, worku2)

c       wt = dt/(tfilt*86400.)
c       do 100 j=latua,latub
c       do 100 i=0,plon+1
c         fu(i,j) = (1.-wt)*fu(i,j) + wt*worku(i,j)
c         fv(i,j) = (1.-wt)*fv(i,j) + wt*worku2(i,j)
c100    continue

c     endif

      return
      end
