C     *==========================================================*
C     | bulkf_variables.f
C     |
C     | Pearse J Buchanan - 2018
C     *==========================================================*
C
C Defines the COMMON block /BULKF/, which contains the main variables
C required for bulk formula forcing in the call to bulkf_formula.f

C ------------- DECLARES MAIN ARRAYS ----------------

C ------ Declare variables for bulk formula --------
      real :: spechum, airtem, slp, swdown, lwdown, fclt, swflx, lwup, 
     +        lwflx, senflx, latflx, hflx, df0dT, evpflx, fflx, rain,
     +        snow, roff, uwnd, vwnd, frzflx, fwiflx, fice, dust, nox

      COMMON/BULKF/
     I             spechum(imt,jmt), airtem(imt,jmt), slp(imt,jmt),
     I             swdown(imt,jmt), lwdown(imt,jmt), fclt(imt,jmt),
     O             swflx(imt,jmt), lwup(imt,jmt), lwflx(imt,jmt), 
     O             senflx(imt,jmt), latflx(imt,jmt), hflx(imt,jmt),
     O             df0dT(imt,jmt), evpflx(imt,jmt), fflx(imt,jmt),
     I             rain(imt,jmt), snow(imt,jmt), roff(imt,jmt),
     I             uwnd(imt,jmt), vwnd(imt,jmt), fice(imt,jmt),
     I             dust(imt,jmt), nox(imt,jmt),
     O             frzflx(imt,jmt), fwiflx(imt,jmt)




