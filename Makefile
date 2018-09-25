# Purpose
# -------
# Makefile for the CSIRO Mk3L climate system model.
#
# Usage
# -----
# 'make' builds the model
# 'make clean' tidies up afterwards
#
# History
# -------
# 2006 May 16	Steven Phipps   Original version
# 2006 Jun 22	Steven Phipps   -fpe0 added to the compiler options for
#				ac.apac.edu.au and lc.apac.edu.au
# 2007 May 31	Steven Phipps   (1) Modified to enable the new ocean model grid
#				(2) Removed redundant configuration
#				    information for sc.apac.edu.au
# 2007 Jun 2	Steven Phipps	Removed redundant source file dencal.f
# 2007 Jun 4	Steven Phipps	Modified for the replacement of the OCEAN_DBL
#				preprocessor macro with OCEAN_LOW
# 2007 Jun 20	Steven Phipps	Added occoup1.f; removed tracer1.f
# 2007 Oct 8	Steven Phipps	Updated the libraries on ac.apac.edu.au to
#				version 9.1 of the Intel Fortran compiler
# 2007 Dec 20	Steven Phipps   Added conserve_fw2.f, downscale.f; removed
#				step0.f
# 2008 Feb 3	Steven Phipps   Added read_real_3d.f90
# 2008 Feb 4	Steven Phipps	Removed r2fcfld.f
# 2008 Feb 5	Steven Phipps	Added read_int_2d.f90, read_real_2d.f90
# 2008 Mar 8	Steven Phipps	Added oflux_read.f, oflux_write.f,
#				orest_read.f, orest_write.f; removed ocemon.f,
#				odam.f
# 2008 Nov 10	Steven Phipps	Macro definitions extracted to separate file
# 2008 Nov 20	Steven Phipps	Updated for conversion of hstring.c to Fortran;
#				merged with Makefile_low
# 2009 Jan      AJA             the makefile components for CABLE are added
# 2009 Apr 8	Steven Phipps	Added m2003.f90
# 2009 May 11	Steven Phipps	Added ocemon.f; tidied up list of object files
# 2009 Aug 6	Steven Phipps	Added hose.f
# 2017 Sep		Pearse Buchanan Added new carbon chemistry routines from MOCSY
#########################
#   Macro definitions   #
#########################

include ../../bld/macros

#######################
#   Compile options   #
#######################

ifeq ($(OCEAN_LOW), yes)
  CPPFLAGS := $(CPPFLAGS) -DOCEAN_LOW
endif

ifeq ($(OCEAN_EOCENE), yes)
  CPPFLAGS := $(CPPFLAGS) -DOCEAN_EOCENE
endif

####################
#   Object files   #
####################

# Atmosphere model

ATM = aplota.o  aploti.o  assel.o  ateday.o  atemon.o  atstart.o  atstep.o \
      bvnfcalc.o  c_aero.o  checkl.o  cldblk.o  clddia.o  cldset.o  clo89.o \
      cloud.o  cloud2.o  cloudm.o  co2_read.o  collst.o  collstx.o  comb_l.o \
      conserve.o  conserve_fw.o  conserve_fw2.o  conv.o  convukmo.o  cvmix.o \
      datard.o  diffn.o  downscale.o  dtogcray.o  dynm.o  dynmnl.o  dynmst.o \
      dynmvo.o  e1e288.o  e3v88.o  ecread.o  energy.o  errcheck.o  esbda.o \
      esibda.o  extend.o  filerd.o  filest.o  filewr.o  finterp.o  fst88.o \
      ftospec.o  ftospec2.o  gauleg.o  gaussv.o  gwdrag.o  hadvect.o \
      hconst.o  hinterp.o  hist_acc.o  hist_cld.o  hist_save.o  hist_wlat.o \
      hkuo.o  hmread.o  hsflux.o  hstring.o  hvertmx.o  icecon.o  icefhx.o \
      icestat.o  icefall.o  icetau.o  inital.o  initax.o  initfs.o \
      initice.o  jmcgslt.o  just_fm.o  landrun.o  lgndre.o  linear.o \
      lwr88.o  main.o  matinv.o  matset.o  mfftga.o  mfftma.o  mmtx.o \
      mtxmtx.o  ncinit.o  ncput.o  nestarc.o  newcloud.o  newrain.o \
      ngwdrag.o  o3_read.o  o3blk.o  o3read_amip.o  o3set.o  o3set_amip.o \
      ocforce.o  ocicurr.o  ocntau.o  ocsave1.o  openfl.o  openhist.o \
      orbpar.o  ordleg.o  phys.o  physgm.o  prdaily.o  prmlomap.o  progcld.o \
      prsmap.o  prtcd.o  prtcl.o  prtt.o  przav.o  ptogcray.o  radcoupl.o \
      radfs.o  radh_l.o  radhyb.o  radin.o  radrhs.o  rads_l.o  radsgrg.o \
      radslab.o  radstres.o  radvars.o  rainda.o  read_int_2d.o \
      read_real_2d.o  read_real_3d.o  readnml1.o  reset.o  seaice.o \
      semiis.o  setqcld.o  solargh.o  source.o  spa88.o  specam.o \
      surfdiag.o  surfset.o  surfa.o  surfupa.o  surfb.o  surfupb.o \
      surfupl.o  swr89.o  table.o  timet.o  tmread.o  traceout.o  tracera.o \
      transf.o  trim.o  trim3x.o  ukall.o  uvharm.o  uvreal.o  vadvect.o \
      vinterp.o  xtchemie.o  xtemiss.o  xtsink.o  xtwetdep.o  zenith.o \
      zerogi.o  zerost.o  hose.o

# Sea ice model

ICE = advect.o  cavit.o  cavit2.o  dynice.o  flatme.o  flatset.o  icebound.o \
      icediag.o  icedrive.o  icefree.o  icesetup.o  ltoh.o  ltou.o \
      polefilt.o  timefilt.o  wrapu.o

# carbon chemistry
MOCSY = gsw_mod_kinds.o gsw_mod_teos10_constants.o gsw_mod_toolbox.o \
		gsw_mod_error_functions.o gsw_mod_baltic_data.o gsw_mod_saar_data.o \
		gsw_mod_specvol_coefficients.o gsw_t_from_ct.o gsw_ct_from_t.o \
		gsw_ct_from_pt.o gsw_pt_from_ct.o gsw_pt_from_t.o gsw_pt0_from_t.o \
		gsw_gibbs_pt0_pt0.o gsw_entropy_part.o gsw_entropy_part_zerop.o \
		gsw_gibbs.o gsw_util_xinterp1.o gsw_util_indx.o gsw_add_barrier.o \
		gsw_add_mean.o gsw_rho.o gsw_specvol.o DNAD.o singledouble.o sw_adtg.o \
		sw_ptmp.o sw_temp.o tpot.o tis.o p80.o phsolvers.o rho.o rhoinsitu.o \
		depth2press.o constants.o varsolver.o vars.o derivauto.o derivnum.o \
		errors.o buffesm.o p2fCO2.o f2pCO2.o gasx.o libmocsy.a 

# Ocean model

OCE = clinic.o  filter.o  findex.o  m2003.o  matrix.o  occoup1.o  ocdatra.o \
      ocdatro.o  ocean.o  ocemon.o  ocend.o  ocfinal.o  ocinit.o  \
      oflux_read.o  oflux_write.o  orest_read.o  orest_write.o  relax.o \
      state.o  step.o  tracer.o bulkf_formula.o bulkf_seaice.o \
      csiro_obgc.o chem.o ocmip_chem.o obgc_write.o

# CABLE

CABLE = cable_define_dimensions.o cable_math_constants.o \
        cable_physical_constants.o cable_photosynthetic_constants.o \
        cable_other_constants.o cable_iovars.o cable_define_types.o \
        cable_abort.o casa_variable.o casa_cable.o casa_cnp.o casa_inout.o \
        cable_soilsnow.o cable_carbon.o cable_radiation.o cable_air.o \
        cable_albedo.o cable_roughness.o cable_parameters.o \
        cable_read.o cable_checks.o cable_initialise.o cable_canopy.o \
        cable_cbm.o cable_input.o cable_write.o cable_output.o CABLE_Mk3L.o        

# Coupled model

OBJ = $(ATM) $(ICE) $(MOCSY) $(OCE) $(CABLE)

############################################
#   Source files requiring preprocessing   #
############################################

F_IN    := $(wildcard *.F)
F_OUT   := $(F_IN:.F=.f)
F90_IN  := $(wildcard *.F90)
F90_OUT := $(F90_IN:.F90=.f90)

####################################
#   Rules for building the model   #
####################################

model : PARAMS.f OPARAMS.f $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIB) -lrfftw -lfftw -lnetcdf -lnetcdff -o model
#	$(FC) $(FFLAGS) $(OBJ) $(LIB) -lrfftw -lfftw -lnetcdf -ltvheap_64 -o model

clean :
	/bin/rm $(F_OUT) $(F90_OUT) *.o *.mod

realclean :
	/bin/rm $(F_OUT) $(F90_OUT) *.o *.mod model

.SUFFIXES : .f .F .f90 .F90 .o

.F.o :
	$(CPP) $(CPPFLAGS) $< > $*.f
	$(FC) $(FFLAGS) -c $(INC) $*.f

.f.o :
	$(FC) $(FFLAGS) -c $(INC) $<

.F90.o :
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(F90) $(F90FLAGS) -c $(INC) $*.f90

.f90.o :
	$(F90) $(F90FLAGS) -c $(INC) $<

# pjb - putting in new carbon chem

GSW = mocsy/src/GSW
GSW_MOD_SRCS := \
    gsw_mod_kinds.o \
    gsw_mod_teos10_constants.o \
    gsw_mod_toolbox.o \
    gsw_mod_error_functions.o \
    gsw_mod_baltic_data.o \
    gsw_mod_saar_data.o \
    gsw_mod_specvol_coefficients.o
#GSW_MOD_OBJS := $(GSW_MOD_SRCS:.f90=.o)
GSW_TOOL_SRCS := \
	gsw_t_from_ct.o \
	gsw_ct_from_t.o \
	gsw_ct_from_pt.o \
	gsw_pt_from_ct.o \
	gsw_pt_from_t.o \
	gsw_pt0_from_t.o \
	gsw_gibbs_pt0_pt0.o \
	gsw_entropy_part.o \
	gsw_entropy_part_zerop.o \
	gsw_gibbs.o \
	gsw_util_xinterp1.o \
	gsw_util_indx.o \
	gsw_add_barrier.o \
	gsw_add_mean.o \
	gsw_rho.o \
	gsw_specvol.o
#GSW_TOOL_OBJS := $(GSW_TOOL_SRCS:.f90=.o)

MOC = mocsy/src
mocsy/src/singledouble.f90 : $(MOC)/singledouble.m4
	m4 -DUSE_PRECISION=$(PRECISION) $^ > $@

MOC_SOURCES = singledouble.o \
              sw_adtg.o \
              sw_ptmp.o \
              sw_temp.o \
              tpot.o \
              tis.o \
              p80.o \
              phsolvers.o \
              rho.o \
              rhoinsitu.o \
              depth2press.o \
              constants.o \
              varsolver.o \
              vars.o \
              derivauto.o \
              derivnum.o \
              errors.o \
              buffesm.o \
              p2fCO2.o \
              f2pCO2.o \
              gasx.o 
#MOBJS := $(MOC_SOURCES:.f90=.o)
library = libmocsy.a

#$(library): $(GSW_MOD_OBJS) $(GSW_TOOL_OBJS) DNAD.o $(MOBJS)
#	ar cr $(INCLUDEFLAGS) $@ $^

$(library): $(GSW_MOD_SRCS) $(GSW_TOOL_SRCS) DNAD.o $(MOC_SOURCES)
	ar cr $@ $^

gsw_mod_kinds.o: $(GSW)/gsw_mod_kinds.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_mod_teos10_constants.o: $(GSW)/gsw_mod_teos10_constants.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_mod_toolbox.o: $(GSW)/gsw_mod_toolbox.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_mod_error_functions.o: $(GSW)/gsw_mod_error_functions.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_mod_baltic_data.o: $(GSW)/gsw_mod_baltic_data.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_mod_saar_data.o: $(GSW)/gsw_mod_saar_data.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_mod_specvol_coefficients.o: $(GSW)/gsw_mod_specvol_coefficients.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_t_from_ct.o: $(GSW)/gsw_t_from_ct.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_ct_from_t.o: $(GSW)/gsw_ct_from_t.f90 
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_ct_from_pt.o: $(GSW)/gsw_ct_from_pt.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_pt_from_ct.o: $(GSW)/gsw_pt_from_ct.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_pt_from_t.o: $(GSW)/gsw_pt_from_t.f90 
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_pt0_from_t.o: $(GSW)/gsw_pt0_from_t.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_gibbs_pt0_pt0.o: $(GSW)/gsw_gibbs_pt0_pt0.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_entropy_part.o: $(GSW)/gsw_entropy_part.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_entropy_part_zerop.o: $(GSW)/gsw_entropy_part_zerop.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_gibbs.o: $(GSW)/gsw_gibbs.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_util_xinterp1.o: $(GSW)/gsw_util_xinterp1.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_util_indx.o: $(GSW)/gsw_util_indx.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_add_barrier.o: $(GSW)/gsw_add_barrier.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_add_mean.o: $(GSW)/gsw_add_mean.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_rho.o: $(GSW)/gsw_rho.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gsw_specvol.o: $(GSW)/gsw_specvol.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

DNAD.o: $(MOC)/DNAD.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

singledouble.o: $(MOC)/singledouble.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

sw_adtg.o: $(MOC)/sw_adtg.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

sw_ptmp.o: $(MOC)/sw_ptmp.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

sw_temp.o: $(MOC)/sw_temp.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

tpot.o: $(MOC)/tpot.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

tis.o: $(MOC)/tis.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

p80.o: $(MOC)/p80.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

phsolvers.o: $(MOC)/phsolvers.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

rho.o: $(MOC)/rho.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

rhoinsitu.o: $(MOC)/rhoinsitu.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

depth2press.o: $(MOC)/depth2press.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

constants.o: $(MOC)/constants.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

varsolver.o: $(MOC)/varsolver.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

vars.o: $(MOC)/vars.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

derivauto.o: $(MOC)/derivauto.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

derivnum.o: $(MOC)/derivnum.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

errors.o: $(MOC)/errors.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

buffesm.o: $(MOC)/buffesm.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

p2fCO2.o: $(MOC)/p2fCO2.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

f2pCO2.o: $(MOC)/f2pCO2.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

gasx.o: $(MOC)/gasx.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@


#%: %.o
#	$(FC) $(FCFLAGS) -o $@ $^
#
#%.o: %.f90
#	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@
#
#%.o: %.F90
#	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

# pjb - end implementation of mocsy carbon chem

#rjm
csiro_obgc.f: $(library) csiro_obgc.F  OPARAMS.f 
	$(CPP) $(crdef) csiro_obgc.F csiro_obgc.f

ocean.f:    ocean.F OPARAMS.f
	$(CPP) $(crdef) ocean.F ocean.f

ocemon.f:    ocemon.F OPARAMS.f
	$(CPP) $(crdef) ocemon.F ocemon.f

# rjm end


PARAMS.f : PARAMS.F
	$(CPP) $(CPPFLAGS) PARAMS.F PARAMS.f

OPARAMS.f : OPARAMS.F
	$(CPP) $(CPPFLAGS) OPARAMS.F OPARAMS.f

bvnfcalc.o : bvnfcalc.f90
	$(F90) $(FIXEDFLAGS) -c $(INC) bvnfcalc.f90

ngwdrag.o : ngwdrag.f90
	$(F90) $(FIXEDFLAGS) -c $(INC) ngwdrag.f90

o3read_amip.o : o3read_amip.f90
	$(F90) $(FIXEDFLAGS) -c $(INC) o3read_amip.f90

o3set_amip.o : o3set_amip.f90
	$(F90) $(FIXEDFLAGS) -c $(INC) o3set_amip.f90

clinic.o : clinic.f
	$(AUTO) $(AUTOFLAGS) -c $(INC) clinic.f

state.o : state.f
	$(AUTO) $(AUTOFLAGS) -c $(INC) state.f

step.o : step.f
	$(AUTO) $(AUTOFLAGS) -c $(INC) step.f

tracer.o : tracer.F
	$(CPP) $(CPPFLAGS) tracer.F > tracer.f
	$(AUTO) $(AUTOFLAGS) -c $(INC) tracer.f

bulkf_formula.o : bulkf_formula.f
	    $(AUTO) $(AUTOFLAGS) -c $(INC) bulkf_formula.f

bulkf_seaice.o : bulkf_seaice.f
	    $(AUTO) $(AUTOFLAGS) -c $(INC) bulkf_seaice.f

# rjm- turn off optimization to get the code to run on osx
ocean.o : ocean.f
	 $(F90)  $(AUTOFLAGS) -O0 -c $(INC) ocean.f

inital.o : inital.f
	 $(F90)  $(AUTOFLAGS) -O0 -c $(INC) inital.f

orest_read.o : orest_read.f
	 $(F90)  $(AUTOFLAGS) -O0 -c $(INC) orest_read.f 

oflux_read.o : oflux_read.f
	 $(F90)  $(AUTOFLAGS) -O0 -c $(INC) oflux_read.f 

readnml1.o : readnml1.f
	 $(F90)  $(AUTOFLAGS) -O0 -c $(INC) readnml1.f 

main.o  : main.F
	$(CPP) $(CPPFLAGS) < main.F  > main.f
	$(FC) $(FFLAGS) -c -O0 $(INC) main.f

collst.o : collst.f
	$(F90)  $(AUTOFLAGS) -O0 -c $(INC) collst.f

filest.o : filest.f
	$(F90)  $(AUTOFLAGS) -O0 -c $(INC) filest.f

ncinit.o : ncinit.f
	$(F90)  $(AUTOFLAGS) -O0 -c $(INC) ncinit.f

atemon.o : atemon.f
	$(F90)  $(AUTOFLAGS) -O0 -c $(INC) atemon.f

# end optimization modifications

# start of the CABLE part

cable_define_dimensions.o: CABLE/cable_define_dimensions.f90
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_define_dimensions.f90

cable_math_constants.o: CABLE/cable_math_constants.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_math_constants.f90

cable_physical_constants.o: CABLE/cable_physical_constants.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_physical_constants.f90

cable_photosynthetic_constants.o: CABLE/cable_photosynthetic_constants.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_photosynthetic_constants.f90

cable_other_constants.o: CABLE/cable_other_constants.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_other_constants.f90

cable_iovars.o: CABLE/cable_iovars.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_iovars.f90

cable_define_types.o: CABLE/cable_define_types.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_define_types.f90

cable_abort.o: CABLE/cable_abort.f90
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_abort.f90

cable_soilsnow.o: CABLE/cable_soilsnow.f90 cable_define_types.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_soilsnow.f90

cable_carbon.o: CABLE/cable_carbon.f90 cable_define_types.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_carbon.f90

cable_radiation.o: CABLE/cable_radiation.f90 cable_define_types.o cable_math_constants.o cable_physical_constants.o cable_other_constants.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_radiation.f90

cable_air.o: CABLE/cable_air.f90 cable_physical_constants.o cable_define_types.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_air.f90

cable_albedo.o: CABLE/cable_albedo.f90 cable_define_types.o cable_other_constants.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_albedo.f90

cable_roughness.o: CABLE/cable_roughness.f90 cable_physical_constants.o cable_define_types.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_roughness.f90

cable_parameters.o: CABLE/cable_parameters.f90 cable_define_types.o cable_iovars.o cable_abort.o
	$(FC) $(CABLEFLAGS) -O0 -c $(INC) CABLE/cable_parameters.f90

cable_read.o: CABLE/cable_read.f90 cable_iovars.o cable_abort.o
	$(FC) $(CABLEFLAGS) -O0 -c $(INC) CABLE/cable_read.f90

cable_checks.o: CABLE/cable_checks.f90 cable_radiation.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_checks.f90

cable_initialise.o: CABLE/cable_initialise.f90 cable_read.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_initialise.f90

cable_canopy.o: CABLE/cable_canopy.f90 cable_radiation.o cable_air.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_canopy.f90

cable_cbm.o: CABLE/cable_cbm.f90 cable_canopy.o cable_roughness.o cable_carbon.o cable_soilsnow.o cable_albedo.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_cbm.f90

cable_input.o: CABLE/cable_input.f90 cable_initialise.o cable_checks.o cable_parameters.o cable_iovars.o cable_read.o cable_abort.o cable_radiation.o
	$(FC) $(CABLEFLAGS) -O0 -c $(INC) CABLE/cable_input.f90

cable_write.o: CABLE/cable_write.f90 cable_iovars.o cable_abort.o cable_define_types.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_write.f90

cable_output.o: CABLE/cable_output.f90 cable_checks.o cable_write.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/cable_output.f90

CABLE_Mk3L.o: CABLE/CABLE_Mk3L.f90 cable_write.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/CABLE_Mk3L.f90

casa_variable.o: CABLE/casa_variable.f90 cable_define_dimensions.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/casa_variable.f90

casa_cable.o: CABLE/casa_cable.f90 cable_define_dimensions.o cable_define_types.o casa_variable.o cable_carbon.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/casa_cable.f90

casa_cnp.o: CABLE/casa_cnp.f90 cable_define_dimensions.o cable_define_types.o casa_variable.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/casa_cnp.f90

casa_inout.o: CABLE/casa_inout.f90 cable_define_dimensions.o cable_define_types.o casa_variable.o cable_iovars.o
	$(FC) $(CABLEFLAGS) -c $(INC) CABLE/casa_inout.f90

