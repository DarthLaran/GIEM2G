FC=mpif90
F77=mpif77
FC_Link=mpif90
AR=xiar

FOPTS= -g    -openmp  -i4 -mcmodel=large  -O3 -traceback    -xHost -prec-div -inline speed -ipo=15 -mkl=sequential 
#FOPTS= -g    -openmp  -i4 -mcmodel=large  -O0


LIB_FFTW=  -L$(HOME)/Soft/lib2/lib -lfftw3_omp -lfftw3_mpi  -lfftw3    
LIB_ZFGMRES=-L./ZFGMRES -lzfgmres
LIB_BLAS=#-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
LIB_ADD=#-lgfortran
INCLUDE= $(HOME)/Soft/lib/include

