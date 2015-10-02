ifneq ($(MAKECMDGOALS), clean)
include $(MAKECMDGOALS).make	
endif
MISC_O=const_module.o timer_module.o  mpi_module.o fftw3_mod.o check_memory_module.o Distributed_FFT.o
#FILTER_WEIGHTS_O=IntegralCodes.o VolumeBesselTransforms.o 
ifdef FFT_QUAD_OURA
	FILTER_WEIGHTS_O=IntegralCodes.o FFT_Quadruple.o  VolumeBesselTransforms.o 
	OPTS=$(FOPTS) -D FFT_QUAD_OOURA
else
	FILTER_WEIGHTS_O=IntegralCodes.o   VolumeBesselTransforms.o  
	OPTS=$(FOPTS)
endif
MODEL_O=data_types_module.o  mpi_saveload_module.o
IMAGE_O= apq_module.o  ie_kernel_hankel_module.o  rc_kernel_hankel_module.o 
IE_O=Integral_Equation_Module.o calc_ie_tensor_module.o ie_solver_module.o   
RC_O=Continuation_Function_Module.o calc_rc_tensor_module.o
SRC_O=sources_module.o
API_O=giem2g_c_api.o

ALL_O=$(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(IMAGE_O) $(IE_O) $(RC_O) $(SRC_O) $(API_O)

LIB_ZFGMRES=-L./ZFGMRES -lzfgmres

LIBS=  $(LIB_ADD)  $(LIB_FFTW) $(LIB_ZFGMRES) $(LIB_BLAS) 
ifneq ($(MAKECMDGOALS), clean)
$(MAKECMDGOALS): zfgmres giem2g 
endif

zfgmres:
	$(MAKE) -C ZFGMRES FC=$(F77)  FOPTS='$(OPTS)' AR=$(AR) FGMRES_PATH='$(FGMRES_PATH)'

giem2g: giem2g_lib giem2g.F90 
	$(FC_Link)   $(OPTS)  giem2g.F90  -L./ -lgiem2g  $(LIBS)  $(INCLUDE) -o giem2g 
	
ifdef INSTALL_PATH
	cp giem2g $(INSTALL_PATH)
endif

giem2g_lib: $(ALL_O)  $(MAKECMDGOALS).make Makefile	
	$(AR) rcs libgiem2g.a $(ALL_O)


%.o:%.f90
	$(FC) $(OPTS) -c $*.f90 -o $*.o $(INCLUDE)

%.o:%.F90
	$(FC) $(OPTS) -c $*.F90 -o $*.o $(INCLUDE)
clean:
	rm $(ALL_O)   *.mod  FFT_Quadruple.o
	$(MAKE) -C ZFGMRES clean
	rm libgiem2g.a
