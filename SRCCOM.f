c $Log: SRCCOM.f,v $
c Revision 1.3  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/12/17  15:31:31  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.1  92/04/15  11:13:38  mrd
c Initial revision
c 
c     common block srccom contains Planck function values used for
c     the radiative calculations: 
      real sorc(imax,lp1,nbly) ! Planck fctn, at model temperatures, for all
                               ! bands used in CTS calculations
      real csour1(imax,lp1)    ! Planck fctn for 560-670 cm-1 band 
      real csour2(imax,lp1)    ! Planck fctn for 670-800 cm-1 band 
      real csour (imax,lp1)    ! Planck fctn for 560-800 cm-1 bands
      real osour (imax,lp1)    ! Planck fctn for 990-1070 cm-1 band
      real ss1   (imax,lp1)    ! Planck fctn for 800-990,1070-1200 cm-1 bands

      common /srccom/ sorc, csour1, csour2, osour, csour ,ss1

