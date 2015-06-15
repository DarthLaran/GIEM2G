FC=mpif90
F77=mpif77
FC_Link=mpif90
AR=xiar

FOPTS= -g    -openmp  -i4 -mcmodel=large  -O3    -xHost -prec-div -inline speed -ipo=22
 #FOPTS= -g    -openmp  -i4 -mcmodel=large  -O0


LIB_FFTW=  -L$(HOME)/Soft/lib/lib -lfftw3_omp -lfftw3_mpi  -lfftw3 
LIB_ZFGMRES=-L./ZFGMRES -lzfgmres
LIB_BLAS=-L$(HOME)/Soft/lib/lib/ -lopenblas_serial_gnu 
LIB_ADD=-lgfortran
INCLUDE= $(HOME)/Soft/lib/include
INSTALL_PATH=~/_scratch/bin/GIEM2G/master

