FC=mpif77
F77=mpif77
FC_Link=mpif77


FOPTS= -g  -free  -openmp  -i4 -mcmodel=large  -O3 -traceback    -xHost -prec-div -inline speed -ipo=9
#FOPTS= -g  -free  -openmp  -i4 -mcmodel=large  -O0


LIB_FFTW=  -L$(HOME)/Soft/lib/lib -lfftw3_omp -lfftw3_mpi  -lfftw3  -lfftw3f_omp  -lfftw3f_mpi  -lfftw3f  
LIB_FGMRES=-L$(HOME)/Soft/lib/lib/ -lfgmres
LIB_BLAS=-L$(HOME)/Soft/lib/lib/ -lopenblas_serial_gnu 
LIB_ADD=-lgfortran
INCLUDE= $(HOME)/Soft/lib/include

