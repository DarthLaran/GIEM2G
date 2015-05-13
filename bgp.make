#FC= /bgsys/drivers/ppcfloor/comm/fast/bin/mpixlf2003_r#gfortran
#FC_Link= /bgsys/drivers/ppcfloor/comm/fast/bin/mpixlf2003_r

FC= mpixlf2003_r#gfortran
FC_OMP=mpixlf2003_r
FC_Link=mpixlf2003_r




DFLAGS= -qessl -g -O5   -qsmp=omp   -qalias=std -qarch=450d -qhot -qmaxmem=-1  -qipa=partition=large  
#DFLAGS= -qessl -g -O0 -qcheck  -qsmp=omp   -qalias=std  -qmaxmem=-1    #-g -free  -r8  -xHost -openmp -cpp -i4



LIB_FFTW=-L/gpfs/opt/bgp/fftw/3.3.3-fast/lib -lfftw3_mpi  -lfftw3_omp   -lfftw3
LIB_FGMRES=-L$(HOME)/lib/lib/lib -lfgmres
LIB_BLAS=-L/opt/ibmmath/lib  -lesslsmpbg 
LIB_ADD=-L$(HOME)/lib/lib/lib -lzspmv -L/opt/ibmcmp/xlmass/bg/4.4/bglib -lmass -lmassv

INCLUDE= -I/opt/ibmmath/include -I/opt/ibmcmp/xlmass/bg/4.4/include -I/gpfs/opt/bgp/fftw/3.3.3-fast/include
