
c======================= include file "bio.h" ===========================
c
c	Additonal variable storage required by the biogeochemical module
c
***** iobiobc	- unit numbers for the random access files
* Options variables
* ivredc - allows the c:p  ratio to vary 
*   =0 constant
*   =1 c/p or o/p varies with pco2 in ocean
*   =2 c/p or o/p varies with pco2 in atmos 
* ivredo - allow the o:p:c ratios to vary
*   =0 p/o fixed at namelist value
*   =1 c/o fixed at namelist value
* ivredpic
*   =0 pic/pop 
*   =1 pic/poc  

* ipic - controls the pic to poc ratio
*    =0 constant
*    =1 varies with omega

* pop		- stores new production and poceralization info
* 	note: new production is stored in pop(i,1,j)
* fmin_pop	- power law function for poceralization
* poc           - store POC
* ratio_pic- ratio_pop ratios of POP for the various tracers	
* sediments	- store sediment accumulation at the ocean bottom
* pic		- production and poceralization of PIC
* rain_ratio	- ratio of PIC/POC fraction
* fmin_pic	- function for PIC poceralization
* ratio_pic	- ratio for PIC poceralization
* terr_input - terrestrial input
* var_satst - variable saturation state scaling based on Ridgwell 2007
* pwr_sat       - additional scaling on the poc:pic based Ridgwell 2007
* dizt	=  1/dzt
c omegaca - store carbonate saturation state for calcite
c omegaar - store carbonate satration state for aragonite
c co2sol  - store the solubility of co2 in surface
c fluxgas	- store air-sea fluxes for oxygen and CO2
c pco2o	- store the sea surface pco2
c kwco2	- gas exchange coefficient for CO2
* pco2_atm	- atmospheric pco2
* pco2a	- air pCO2
* fc13a	- C13 fraction in air 
c igas	- an index to specify the gas exchange to use
c		0 - no gas exchange
c		1 - oxygen gas exchange
c		2 - CO2
c		3 - F12 
c		4 - F11
c*		5 - C13
* tglobal 	- stores the global average surface tracer concentration for the virtual fluxes
*  Tracer indexes
* n_pho  	- index for the nutrient limited tracer, if
*		 0 new production is zero
* n_dic	- index for the DIC tracer, if 0 dic=2000
* n_alk	- index for the Alk tracer, if 0 alk=2431
* n_c14	- index for simple c14 run
* n_si  - index for the silca tracer
* n_oxy	- index for oxygen, usually 3 and must be define for obgc
* n_c13	- dic-c13 index  (needs n_dic to be defined)
* n_fe	- index for iron tracer
* n_no3 - index for nitrate tracer 
* n_n15 - index for Nitrogen 15 tracer
* n_tou - index for True Oxygen Utilisation
* n_age - index for age of the water (years)
* n_tomz - index for age of the water in omzs (years)
* n_n2o - index for nitrous oxide (nmol/m3)
  
* obc  - biological surface bc
* index 1  = windspeed
* index 2  = phytoplankton growth
* index 3  = surface pressure
* index 4  = dust deposition
* index 5  = denitrification
* index 6  = N fixation
* index 7  = new remineralisation and stoichiometric calcs
* index 8  = Nitrogen isotopes
* index 9  = Optimal Uptake Kinetics
* index 1  = Altered remineralisation
* index 2  = Nitrous oxide cycling

* ======================================================
        integer ngas,ntmax,n_pho,n_dic,n_alk,n_c14,n_si,n_oxy,n_c13,
     1  n_fe,n_no3,n_n15,n_tou,n_age,n_tomz,n_n2o,
     1  ntdbc,maxbcr,npar,igas, ntm2,nbio2d,
     1 idump,jdump,kdump,
     1 nsum1,itlast
     1 ,ivredc,ivredo,ipic, latj, loni

      logical :: den, fix, atmdep, vary_stoich, NO3diagnostics,
     9           MichaelisMenton, OptimalUptake, 
     1           ReminRelax, ReminTemp_Marsay, ReminTemp_Q10,
     1           ReminPico, ReminDiagnostics, biomol, solarcarb, 
     8           N15diagnostics,
     2           N2Otemp, N2Odepth,
     7           sedfluxes, subgrid, 
     1           vary_frac13C, diff_out
      
      character (len=7) :: Den_type
     
	real poc,pop,pic,sediments,fmin_pop,ratio_pop,pop_fix,fmin_dia,
     1 dizt,tr_off,fmin_opal,opal,ratio_opal,terr_input,ibdt,
     1 rain_ratio,fmin_pic,ratio_pic,pwr_sat,var_satst(imt,jmt),
     1 fluxgas,poc2,kwco2,pco2a,fc13a,pco2_atm,
     1 pN2O, omz_min, omz_mean, omz_max,
     1 Q10, T_rem,
     1 pco2o,pmodel,omegaca,omegaar,co2sol,fco3,co2sw,rhosw,x_co2,
     1 tk,invtk,pH2O,pco2atm,fco2atm,tmp,K0,co2starair,co2star,dpco2,
     1 f_dis,f_spe,c14o,c14a, del14c_atm,
     1 isurf,fixer,
     1 ttflux,stflux,gCtppm,geo_seq,xdump,pdump,
     1 p_k,s_npp,xmld,fpgr_temp,vmax,n_k,fe_k,phs,no3,
     1 obc,sbcbio,silicate,pwl_pop,pwl_pic,pwl_opal,zpop,zpic,zdia,
     1 atm_box,atm_box_tot,atm_box_ann,atm_box_pco2,
     1 atm_box_rc,atm_box_tot_c13,atm_box_ann_c13,
     1 p11,p12,pcfcatm,
     1 zt,dzt,
     1 Ro, dtr, 
     1 ave_tr,ave_flux,tglobal,
     5 oxyneg, denP, denS, sulS, Boh_a, ave_Ssul, no3change, vol,
     5 o2loss, Nl_Pden, Nl_Sden, no3fix, no3dif, budget, Nexp,
     6 Nfix, fixneg, ave_fixneg, ave_Pden, ave_fix, ave_Sden, s_npp_fix,
     6 Fmax, temp, oxy, dFe, Fe_kd, p_kd, R_fix, FixLim, phs_fix,
     6 CP_fix, NP_fix, HP_fix, OP_fix, o2_rem_f, no3_rem_f, ave_exp,
     6 op_tot_f, op_rem_f, op_den_f, op_unrem_f, R_Olim_f, R_Nlim_f,
     7 carb2P, NtoP, HtoP, OtoP, o2_rem, no3_rem, c13toc, 
     7 ave_carb2P, ave_NtoP, ave_CtoP, ave_o2rem, ave_no3rem, 
     7 op_sed, ave_opsed, ic_sed, ave_icsed, op_sed_f, ave_opsed_f,
     7 si_sed, ave_sised, sgb_frac, sedtot, ssulG, ssulF,
     7 o_n, tot_a, tot_b, sdenTg, sdenTf, sdenG, sdenF, sremG, sremF,
     5 R_den, R_Nlim, deNrelax, op_den, Den_o2lim, Den_nlim, Sden_acl,
     8 n15,n14,pon15,pon14,on_rem,on_rem_f,no3_den,no3_den_f,
     8 IUC, NUC, NUCpel, NUCsed, f_ass, f_pel, f_sed,
     8 d_n15,d_n15f,d_n15d,d_n15df, d15Norg, ave_d15N,
     8 Eass, Erem, Efix, Eatm, Ewc, Esed,
     8 ponf15, ponf14, pon_cons, ponf_cons,
     9 fa, Pou, Nou, VoA,
     1 R_Olim, op_unrem, ave_oxyneg, op_rem,
     1 op_tot, OM, tot, dif, 
     1 ave_pop, ave_optot, ave_oprem, ave_opden,
     1 EP, EPmax, dep, Mcurve_b, ave_MCb, Fpico, Teff,
     2 x_N2O, N2Oo, N2Oa, age, tour, d_n2o, tomz, n2o_Pden, n2o_Cden,
     1 Xdif, Ydif, Wdif, Gar84, ave_Udif, ave_Vdif, ave_Wdif, ave_Gar84

**    MOCSY STUFF - CARBON CHEMISTRY     
*      REAL(kind=1.0), DIMENSION(1) :: CCph,CCpco2,CCfco2,CCco2,CChco3,
*     1                    CCco3, CCOmegaA, CCOmegaC, CCBetaD, CCrhoSW,
*     1                    CCp, CCtempis
*      REAL(kind=1.0), DIMENSION(1) :: CCtem, CCsal, CCalk, CCdic, CCsil,
*     1                    CCpo4, CCpatm, CCdep, CClat, CClon
*      character(10) :: optCON = 'mol/m3'
*      character(10) :: optT = 'Tpot'
*      character(10) :: optP = 'm'


	parameter (ngas=nt-2+1)
	parameter (nbio2d = nt+14)
	parameter (ntm2=nt-2+1)  ! add one to ensure it has a nonzero value
	parameter (ntmax=15)
	character *12 trname(nbio2d)
        common /cbio11/ poc(imt,km,jmt),pop(imt,km,jmt),pic(imt,km,jmt),
     1	sediments(imt,jmt), pop_fix(imt,km,jmt),fmin_dia(imt,jmt,km),
     1	fmin_pop(imt,jmt,km), ratio_pop(ntmax),terr_input,pwl_pop,
     1	rain_ratio,fmin_pic(imt,jmt,km),ratio_pic(ntmax),pwr_sat,
     1	dizt(km),omegaca(imt,jmt),omegaar(imt,jmt),co2sol(imt,jmt),
     1  co2sw(imt,jmt), rhosw(imt,jmt), fco3(imt,jmt), x_co2(imt,jmt),
     1	tr_off(ntmax),n_pho,n_dic,n_alk,n_c14,n_oxy,n_c13,f_dis,f_spe,
     1	fmin_opal(imt,jmt,km),opal(imt,km,jmt),ratio_opal(ntmax),
     1	n_si,n_fe,n_no3,n_n15,n_tou,n_age,n_tomz,n_n2o,
     1  ivredc,ivredo,zpop,zdia,ipic,zpic,pwl_pic,pwl_opal,
     5  den,fix,atmdep,vary_stoich,Den_type,NO3diagnostics,
     6  Fe_kd, p_kd, s_npp_fix, CP_fix, NP_fix, o2_rem_f,no3_rem_f,
     8  N15diagnostics, latj, loni,
     9  MichaelisMenton, OptimalUptake, biomol, solarcarb,
     1  ReminRelax, ReminTemp_Marsay, ReminTemp_Q10, 
     1  ReminPico, ReminDiagnostics, 
     1  op_rem(imt,km,jmt), op_den(imt,km,jmt),
     1  op_tot(imt,km,jmt), OM(imt,jmt),
     2  x_N2O(imt,jmt), N2Otemp, N2Odepth,
     7  sedfluxes, subgrid, 
     1  vary_frac13C, diff_out

	common /cgas/ igas(ntmax),fluxgas(imt,jmt,ngas),
     1	pco2o(imt,jmt,ngas),kwco2(imt,jmt),pco2a,fc13a,
     1  pco2_atm(-1:1),
     2	ttflux(nt),stflux(nt),gCtppm,geo_seq(ntmax),
     3  ave_tr(imt,jmt,km,nt),ave_flux(imt,jmt,nbio2d),
     4	tglobal(ntmax),
     5  oxyneg(imt,km,jmt), denP(imt,km,jmt), no3change(imt,km,jmt),
     5  o2loss(imt,km,jmt), Nl_Pden(imt,km,jmt), no3fix(imt,km,jmt), 
     5  no3dif(imt,km,jmt),budget(imt,km,jmt), ave_Pden(imt,jmt,km,1),
     6  Nfix(imt,km,jmt), Nl_Sden(imt,jmt,km), denS(imt,km,jmt),
     6  ave_fix(imt,jmt,km,1), sulS(imt,km,jmt), ave_Ssul(imt,jmt,km,1),
     6  ave_Sden(imt,jmt,km,1), op_tot_f(imt,km,jmt), Nexp(imt,km,jmt),
     6  op_rem_f(imt,km,jmt),op_den_f(imt,km,jmt),ave_exp(imt,jmt,km,1),
     6 op_unrem_f(imt,km,jmt),R_Olim_f(imt,km,jmt),R_Nlim_f(imt,km,jmt),
     7  carb2P(imt,jmt), NtoP(imt,jmt), HtoP(imt,jmt), OtoP(imt,jmt),
     7  o2_rem(imt,jmt), no3_rem(imt,jmt), c13toc(imt,jmt),
     7  ave_carb2P(imt,jmt,1), ave_NtoP(imt,jmt,1),
     7  ave_o2rem(imt,jmt,1), ave_no3rem(imt,jmt,1),
     7  op_sed(imt,km,jmt),ave_opsed(imt,jmt,km,1),sedtot(imt,jmt),
     7  op_sed_f(imt,km,jmt),ave_opsed_f(imt,jmt,km,1),
     7  ic_sed(imt,km,jmt),ave_icsed(imt,jmt,km,1),
     7  si_sed(imt,km,jmt), ave_sised(imt,jmt,km,1),
     7  sdenG(imt,km,jmt), sdenF(imt,km,jmt),
     7  sremG(imt,km,jmt), sremF(imt,km,jmt),
     7  ssulG(imt,km,jmt), ssulF(imt,km,jmt),
     5  R_den(imt,km,jmt), R_fix(imt,jmt),
     5  Boh_a, Sden_acl, Den_o2lim, Den_nlim, R_Nlim(imt,km,jmt),
     8  d_n15(imt,km,jmt), d_n15f(imt,km,jmt), d_n15d(imt,km,jmt),
     8  d_n15df(imt,km,jmt),d15Norg(imt,km,jmt),
     8  pon15(imt,km,jmt),pon14(imt,km,jmt),
     8  ponf15(imt,km,jmt),ponf14(imt,km,jmt),
     8  Eass, Erem, Efix, Eatm, Ewc, Esed, ave_d15N(imt,jmt,km,1),
     8  pon_cons(km), ponf_cons(km), on_rem(km), on_rem_f(km),
     8  no3_den(km), no3_den_f(km),
     1  R_Olim(imt,km,jmt), op_unrem(imt,km,jmt), 
     1  ave_oxyneg(imt,jmt,km,1), ave_pop(imt,jmt,km,1),
     1  ave_optot(imt,jmt,km,1), ave_oprem(imt,jmt,km,1), 
     1  ave_opden(imt,jmt,km,1),
     1  EP(imt,jmt), Mcurve_b(imt,jmt), ave_MCb(imt,jmt,1),
     2  d_n2o(imt,km,jmt), pN2O, omz_min, omz_mean, omz_max,
     1  Q10, t_rem,
     1  Xdif(imt,jmt,km), Ydif(imt,jmt,km), Wdif(imt,jmt,km),
     1  Gar84(imt,jmt,km), ave_Udif(imt,jmt,km,1),
     1  ave_Vdif(imt,jmt,km,1), ave_Wdif(imt,jmt,km,1),
     1  ave_Gar84(imt,jmt,km,1)

        common /ceco/ p_k,s_npp,xmld(imt,jmt),
     1  fpgr_temp(imt,jmt),vmax(imt,jmt),n_k,fe_k,VoA

        common /cdumping/idump,jdump,kdump,xdump,pdump

        common /obgc_names/ trname

c
c	Time dependent S.B.C. data for biogeochemical module 
c     coded by:     rj matear 
        parameter (ntdbc=6,maxbcr=12)
        common /bio_bcr/ obc(imt,jmt,maxbcr,ntdbc)
        common /bio_bcr/ sbcbio(imt,jmt,ntdbc)
        common /bio_bcr/ sgb_frac(imt,jmt,km)

*   Parameters sizes for the ecosystem model.  See the 
* ecosystem directory file param.h for additional parameters.
        parameter (npar=16) 

*  Initialize the Ecosystem model parameters and useful constants
*
*  Store the constants used by the ecosystem model
        logical qrad
        common /bio_const/ qrad, Ro, dtr

* 	Storage of variables required to run the ecoystem model.  See
* ecosystem directory file eco_input.h for additional variables.
        dimension pmodel(npar)
        common /eco_input/ pmodel
	
c...
        dimension zt(kmp1),dzt(km)
        equivalence (zt,zdzz), (dzt,dz)
c...

*  Include file for the cfc simulations
*     coded by:     rj matear 
        common /cfc_1/  p11(68,2), p12(68,2), pcfcatm(imt,jmt,2)

        common /atm_box/ atm_box_pco2,atm_box_tot,atm_box_ann
     1  ,atm_box_tot_c13,atm_box_ann_c13,atm_box_rc

* include information for the ocmip chemistry run
        common /ocmip1/ silicate(imt,jmt)
