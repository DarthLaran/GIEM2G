ifneq ($(MAKECMDGOALS), clean)
include make.inc	
else
	FFT_QUAD_OOURA=1
endif
ifdef FFT_QUAD_OOURA
	FILTER_WEIGHTS_O=IntegralCodes.o FFT_Quadruple.o  VolumeBesselTransforms.o 
	OPTS=$(FOPTS) -D FFT_QUAD_OOURA
else
	FILTER_WEIGHTS_O=IntegralCodes.o   VolumeBesselTransforms.o  
	OPTS=$(FOPTS)
endif
MISC_O=const_module.o  mpi_module.o timer_module.o  Logger_Module.o check_memory_module.o
FFT_O=fftw3_mod.o Distributed_FFT.o Local_OMP_FFT_module.o
SAVE_LOAD_O=mpi_saveload_module.o
MODEL_O=data_types_module.o 
IE_IMAGE= apq_module.o  ie_kernel_hankel_module.o 
RC_IMAGE=rc_kernel_hankel_module.o 
IE_O=Integral_Equation_Module.o calc_ie_tensor_module.o apply_ie_operator_module.o
IE_SOLVER=ie_solver_module.o   
RC_O=Continuation_Function_Module.o calc_rc_tensor_module.o
SRC_O=sources_module.o
API_O=giem2g_interfaces.o
ALL_O=$(FSON) $(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(SAVE_LOAD_O) $(FFT_O) $(IE_IMAGE) $(IE_O) $(IE_SOLVER) $(RC_IMAGE)  $(RC_O) $(SRC_O) 
ENGINE_O=$(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(FFT_O)  $(IE_IMAGE) $(IE_O)  $(API_O) 


LIB_ZFGMRES=-L./ZFGMRES -lzfgmres
LIB_GFGMRES=-L./GFGMRES -lgfgmres
LIB_FGMRES=$(LIB_GFGMRES)
LIB_FSON=-L./FSON -lfson

LIBS=  $(LIB_ADD)  $(LIB_FFTW) $(LIB_FGMRES) $(LIB_BLAS)  $(LIB_FSON)


giem2g: zfgmres  giem2g_lib giem2g.F90 
	$(FC_Link)   $(OPTS)  giem2g.F90  -L./ -lgiem2g  $(LIBS)  $(INCLUDE) -o giem2g 
	
ifdef INSTALL_PATH
	cp giem2g $(INSTALL_PATH)
endif

giem2g_lib_shared: $(ENGINE_O)
	$(FC_Link)   $(ENGINE_O) -shared   -L$(SHARED_BLAS) -L$(SHARED_FFTW)  -o $(DST)/libgiem2g.so  

#	$(SHARED_Link) -fPIC  -shared   -L$(SHARED_BLAS) -L$(SHARED_FFTW) $(ENGINE_O)   -o $(DST)/lib_giem2g.so 


zfgmres:
#	$(MAKE) -C ZFGMRES FC=$(F77)  FOPTS='$(OPTS)' AR=$(AR) FGMRES_PATH='$(FGMRES_PATH)'
	$(MAKE) -C GFGMRES FC=$(F77)  FOPTS='$(OPTS)' AR=$(AR) 
fson:
	$(MAKE) -C FSON FC=$(FC)  FOPTS='$(OPTS)' AR=$(AR) 

giem2g_lib: $(ALL_O)   Makefile	
	$(AR) rcs libgiem2g.a $(ALL_O)


%.o:%.f90
	$(FC) $(OPTS) -c $*.f90 -o $*.o $(INCLUDE)

%.o:%.F90
ifneq ($(SHARED_LIB),1)
	$(FC) $(OPTS) -c $*.F90 -o $*.o $(INCLUDE)
else
	$(FC) $(OPTS) -fPIC -D ExtrEMe   -c $*.F90 -o $*.o -I$(SHARED_BLAS_INC)  -I$(SHARED_FFTW_INC) 
endif
clean:
	$(MAKE) -C ZFGMRES clean
	$(MAKE) -C GFGMRES clean 
	rm $(ALL_O)   *.mod  
	rm libgiem2g.a
