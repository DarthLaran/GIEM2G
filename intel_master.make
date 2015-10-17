FC=mpif90
F77=mpif77
FC_Link=mpif90
AR=xiar

FOPTS= -g    -openmp  -i4 -mcmodel=large  -O3    -xHost -prec-div -inline speed -ipo=25 -DLEGACY_MPI 

 FFT_QUAD_OURA=1

LIB_FFTW= -L$(HOME)/Soft/lib/fftw3_ompi155_icc/lib  -lfftw3_omp -lfftw3_mpi  -lfftw3 
#LIB_FFTW= -L$(HOME)/Soft/lib/fftw3_ompi155_gcc502/lib  -lfftw3_omp -lfftw3_mpi  -lfftw3 
FGMRES_PATH=/mnt/data/users/dm3/vol6/kruglyakov/Soft/src/fgmres
LIB_BLAS=-L$(HOME)/Soft/lib/OpenBLAS_gcc44/lib -lopenblas_omp 
LIB_ADD=-lgfortran
INCLUDE= -I$(HOME)/Soft/lib/fftw3_ompi155_icc/include
INSTALL_PATH=~/_scratch/bin/GIEM2G/giem2g_master

