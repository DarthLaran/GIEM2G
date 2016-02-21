
ifeq ($(PLATFORM), LAPTOP)


endif

ifeq ($(PLATFORM), DAINT) 

FC=ftn
F77=ftn
FC_Link=ftn
AR=ar

FOPTS= -O 3 -O aggress  -O cache3 -O inlinelib  -O vector3 -O scalar3   -O ipa4 -s integer32 -D FFTW_SIZE_32=1 -D MEMORY_CHECKING # -D performance_test
FFT_QUAD_OOURA=1
LIB_FFTW=   -lfftw3_omp -lfftw3_mpi  -lfftw3 
FGMRES_PATH=~/src/FGMRES/
LIB_BLAS= -L ~/lib/lib/ -lopenblas_omp
LIB_ADD= -lgfortran
INCLUDE=  # -I$(HOME)/Soft/lib/include
INSTALL_PATH=./

endif

ifeq ($(PLATFORM), BGP)


endif
ifeq ($(PLATFORM), LOMONOSOV)


endif
