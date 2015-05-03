#FC= /bgsys/drivers/ppcfloor/comm/fast/bin/mpixlf2003_r#gfortran
#FC77= /bgsys/drivers/ppcfloor/comm/fast/bin/mpixlf77_r#gfortran
#FC_OMP= /bgsys/drivers/ppcfloor/comm/fast/bin/mpixlf2003_r
#FC_Link= /bgsys/drivers/ppcfloor/comm/fast/bin/mpixlf2003_r

FC= mpixlf2003_r#gfortran
FC77=mpixlf77_r#gfortran
FC_OMP=mpixlf2003_r
FC_Link=mpixlf2003_r
#FFTW= -L$(HOME)/lib/lib/lib -lfftw3_mpi  -lfftw3_omp   -lfftw3
FFTW=-L/gpfs/opt/bgp/fftw/3.3.3-fast/lib -lfftw3_mpi  -lfftw3_omp   -lfftw3




DFLAGS= -qessl -g -O5   -qsmp=omp   -qalias=std -qarch=450d -qhot -qmaxmem=-1  -qipa=partition=large   #-g -free  -r8  -xHost -openmp -cpp -i4
#DFLAGS= -qessl -g -O0 -qcheck  -qsmp=omp   -qalias=std  -qmaxmem=-1    #-g -free  -r8  -xHost -openmp -cpp -i4



MISC_O=const_module.o mpi_module.o
FILTER_WEIGHTS_O=IntegralCodes.o VolumeBesselTransforms.o  
MODEL_O=data_types_module.o mpi_saveload_module.o
IMAGE_O= apq_module.o  ie_kernel_hankel_module.o  rc_kernel_hankel_module.o 
IE_O=Integral_Equation_Module.o calc_ie_tensor_module.o ie_solver_module.o   
RC_O=Continuation_Function_Module.o calc_rc_tensor_module.o
SRC_O=sources_module.o


ALL_O=$(MISC_O) $(FILTER_WEIGHTS_O) $(MODEL_O) $(IMAGE_O) $(IE_O) $(RC_O) $(SRC_O) 
LIB_FGMRES=-L$(HOME)/lib/lib/lib -lfgmres
LIB_ZSPMV=-L$(HOME)/lib/lib/lib -lzspmv #FOR BLUEGENE ONLY!
LIB_BLAS=-L/opt/ibmmath/lib  -lesslsmpbg 
LIB_ADD=-L/opt/ibmcmp/xlmass/bg/4.4/bglib -lmass -lmassv
LIBS= $(LIB_ZSPMV) $(LIB_FGMRES) $(LIB_BLAS) $(FFTW)  $(LIB_ADD)
INCLUDE= -I/opt/ibmmath/include -I/opt/ibmcmp/xlmass/bg/4.4/include -I/gpfs/opt/bgp/fftw/3.3.3-fast/include

MOD_PATH=$(HOME)/lib/include -I/gpfs/opt/bgp/fftw/3.3.3-fast/include
giem2g:  $(ALL_O) giem2g.F90 Makefile
	$(FC_Link)  -qsmp=omp $(DFLAGS)  -o giem2g $(ALL_O) giem2g.F90 $(LIBS)  -I$(INCLUDE) -I$(MOD_PATH)

%.o:%.f90
	$(FC) $(DFLAGS) -c $*.f90 -I$(MOD_PATH)

%.o:%.F90
	$(FC) $(DFLAGS) -c $*.F90 -I$(MOD_PATH)
clean:
	rm $(ALL_O)   *.mod
