c Wind stresses removed, as these do not need to be saved to a restart file.
c REAL arrray ODAM_FKM replaced with INTEGER arrays KMT and KMU.
c SJP 2003/09/04
c
c Arrays containing integer indices converted to type INTEGER. Some variables
c transferred from COMMON block /fullwd/, with the duplicates held within
c /orestart/ being deleted.
c SJP 2003/09/03
c
c Defines COMMON block /orestart/, which contains the data read to/written from
c the OGCM restart file. This data was previously held in memory by the Ocean
c Direct Access Menager (ODAM).
c SJP 2003/09/02

      integer itt, isz(jmt, lseg), iez(jmt, lseg), isis(nisle),
     &        ieis(nisle), jsis(nisle), jeis(nisle),
     &        istf(njtbft, lsegf, km), ietf(njtbft, lsegf, km),
     &        isuf(njtbfu, lsegf, km), ieuf(njtbfu, lsegf, km),
     &        iszf(njtbfu, lsegf), iezf(njtbfu, lsegf),
     &        kmt(imt, jmt), kmu(imt, jmt)
      real ttsec, area, volume, pb(imt, jmt), p(imt, jmt),
     &     hr(imt, jmt), ptd2(imt, jmt, 2),
     &     odam_t(imt, jmt, km, nt, 2),
     &     odam_u(imt, jmt, km, 2), odam_v(imt, jmt, km, 2)

      common /orestart/ itt, ttsec, area, volume, pb, p, hr, ptd2,
     &                  isz, iez, isis, ieis, jsis, jeis, istf, ietf,
     &                  isuf, ieuf, iszf, iezf, odam_t, odam_u, odam_v,
     &                  kmt, kmu
