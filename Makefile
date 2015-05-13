include $(MAKECMDGOALS).make	

MISC_O=const_module.o mpi_module.o fftw3_mod.o
FILTER_WEIGHTS_O=IntegralCodes.o VolumeBesselTransforms.o  
MODEL_O=data_types_module.o mpi_saveload_module.o
IMAGE_O= apq_module.o  ie_kernel_hankel_module.o  rc_kernel_hankel_module.o 
IE_O=Integral_Equation_Module.o calc_ie_tensor_module.o ie_solver_module.o   
RC_O=Continuation_Function_Module.o calc_rc_tensor_module.o
SRC_O=sources_module.o


ALL_O=$(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(IMAGE_O) $(IE_O) $(RC_O) $(SRC_O)


LIBS=  $(LIB_FGMRES) $(LIB_BLAS) $(LIB_FFTW)  $(LIB_ADD)


$(MAKECMDGOALS): giem2g

giem2g:  $(ALL_O) giem2g.F90 $(MAKECMDGOALS).make Makefile	
	$(FC_Link)   $(DFLAGS)  -o giem2g $(ALL_O) giem2g.F90 $(LIBS)  -I$(INCLUDE) 

%.o:%.f90
	$(FC) $(DFLAGS) -c $*.f90 -I$(INCLUDE)

%.o:%.F90
	$(FC) $(DFLAGS) -c $*.F90 -I$(INCLUDE)
clean:
	rm $(ALL_O)   *.mod
