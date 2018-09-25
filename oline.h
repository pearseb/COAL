* ncdread.h file
	integer s1(4), c1(4)
	integer  c2(4)
	integer s3(4), c3(4)
	integer s4(3), c4(3)
	common /c_oline/ tfix(imt,jmt,km,1),sfix(imt,jmt,km,1),
	1	ufix(imt,jmt,km,1),vfix(imt,jmt,km,1),
	1  tfold(imt,jmt,km,1),sfold(imt,jmt,km,1),
	1	ufold(imt,jmt,km,1),vfold(imt,jmt,km,1)
	2	, tsurf(imt,jmt,1,1),tsurf_old(imt,jmt,1,1),tsurf_u(imt,jmt)
	3	, isice(imtm2,jmtm2,1,1),isvmo(imtm2,jmtm2,1,1),
	4	svmoref(imtm2,jmtm2,1)
	5 ,smfzon(imt,jmt,1),smfmer(imt,jmt,1),stfht(imt,jmt,1),
	6 stfsal(imt,jmt,1),res(imt,jmt,1)

	data s1 /1,1,1,1/, c1 /imt,jmt,km,1/
	data c2 /imt,jmt,1,1/
	data s3 /1,1,1,1/, c3 /imtm2,jmtm2,1,1/
	data s4 /1,1,1/, c4 /imtm2,jmtm2,1/
