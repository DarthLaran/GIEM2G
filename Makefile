
FC=mpif77
FC_OMP=mpif77
FC_Link=mpif77
FFTW=  -L$(HOME)/Soft/lib/lib -lfftw3_omp -lfftw3_mpi  -lfftw3  -lfftw3f_omp  -lfftw3f_mpi  -lfftw3f  


DFLAGS= -g  -free  -openmp  -i4 -mcmodel=large  -O3 -traceback    -xHost -prec-div -inline speed -ipo=9
#DFLAGS= -g  -free  -openmp  -i4 -mcmodel=large  -O0


MISC_O=const_module.o mpi_module.o
FILTER_WEIGHTS_O=IntegralCodes.o VolumeBesselTransforms.o  
MODEL_O=data_types_module.o mpi_saveload_module.o
IMAGE_O= apq_module.o  ie_kernel_hankel_module.o  rc_kernel_hankel_module.o 
IE_O=Integral_Equation_Module.o calc_ie_tensor_module.o ie_solver_module.o   
RC_O=Continuation_Function_Module.o calc_rc_tensor_module.o
#SRC_O=impedances_module.o calc_hankel_images_Nlayers.o sources.o
SRC_O=sources_module.o


ALL_O=$(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(IMAGE_O) $(IE_O) $(RC_O) $(SRC_O) 
LIB_FGMRES=-L$(HOME)/Soft/lib/lib/ -lfgmres
LIB_ZSPMV=#-L$(HOME)/lib/lib/lib -lzspmv #FOR BLUEGENE ONLY!
LIB_BLAS=-L$(HOME)/Soft/lib/lib/ -lopenblas_serial_gnu 
LIB_ADD=-lgfortran
LIBS= $(LIB_ZSPMV) $(LIB_FGMRES) $(LIB_BLAS) $(FFTW)  $(LIB_ADD)
INCLUDE= $(HOME)/Soft/lib/include
MOD_PATH=$(HOME)/Soft/lib/modules
giem2g:  $(ALL_O) giem2g.F90 Makefile
	$(FC_Link)   $(DFLAGS)  -o giem2g $(ALL_O) giem2g.F90 $(LIBS)  -I$(INCLUDE) -I$(MOD_PATH)

%.o:%.f90
	$(FC) $(DFLAGS) -c $*.f90 -I$(MOD_PATH)

%.o:%.F90
	$(FC) $(DFLAGS) -c $*.F90 -I$(MOD_PATH)
clean:
	rm $(ALL_O)   *.mod
