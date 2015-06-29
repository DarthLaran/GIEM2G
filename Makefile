ifneq ($(MAKECMDGOALS), clean)
include $(MAKECMDGOALS).make	
endif
MISC_O=const_module.o mpi_module.o fftw3_mod.o
#FILTER_WEIGHTS_O=IntegralCodes.o VolumeBesselTransforms.o  
FILTER_WEIGHTS_O=IntegralCodes.o FFT_Quadruple.o  VolumeBesselTransforms.o  
MODEL_O=data_types_module.o mpi_saveload_module.o
IMAGE_O= apq_module.o  ie_kernel_hankel_module.o  rc_kernel_hankel_module.o 
IE_O=Integral_Equation_Module.o calc_ie_tensor_module.o ie_solver_module.o   
RC_O=Continuation_Function_Module.o calc_rc_tensor_module.o
SRC_O=sources_module.o


ALL_O=$(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(IMAGE_O) $(IE_O) $(RC_O) $(SRC_O)

LIB_ZFGMRES=-L./ZFGMRES -lzfgmres

LIBS=  $(LIB_ADD)  $(LIB_FFTW) $(LIB_ZFGMRES) $(LIB_BLAS) 
ifneq ($(MAKECMDGOALS), clean)
$(MAKECMDGOALS): zfgmres giem2g 
endif

zfgmres:
		$(MAKE) -C ZFGMRES FC=$(F77)  FOPTS='$(FOPTS)' AR=$(AR) FGMRES_PATH='$(FGMRES_PATH)'

giem2g:  $(ALL_O) giem2g.F90 $(MAKECMDGOALS).make Makefile	
	$(FC_Link)   $(FOPTS)   -o giem2g $(ALL_O) giem2g.F90 $(LIBS)  $(INCLUDE) 
	
ifdef INSTALL_PATH
	cp giem2g $(INSTALL_PATH)/giem2g
endif
%.o:%.f90
	$(FC) $(FOPTS) -c $*.f90 -o $*.o $(INCLUDE)

%.o:%.F90
	$(FC) $(FOPTS) -c $*.F90 -o $*.o $(INCLUDE)
clean:
	rm $(ALL_O)   *.mod
	$(MAKE) -C ZFGMRES clean
