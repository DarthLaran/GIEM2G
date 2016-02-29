
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

ifeq ($(PLATFORM), BGP_FAST)

FC= mpixlf2003_r
F77= mpixlf77_r
FC_Link=mpixlf2003_r
AR=ar

FFLAGS= -qessl -g -O5   -qsmp=omp   -qalias=std -qarch=450d -qhot -qmaxmem=-1  -qipa=partition=large -qunroll=yes -qipa=level=2
FDEFS=-WF,-DLEGACY_MPI -WF,-DIBM_Bluegene  -WF,-DMPI_TIMER 
FOPTS=$(FFLAGS) $(FDEFS)

LIB_FFTW=-L/gpfs/opt/bgp/fftw/3.3.3-fast/lib -lfftw3_mpi  -lfftw3_omp   -lfftw3
FGMRES_PATH=/gpfs/data/kruglyakov/lib/src/FGMRES_original
LIB_BLAS=-L/opt/ibmmath/lib  -lesslsmpbg 
LIB_ADD=-L$(HOME)/lib/lib/lib -lzspmv -L/opt/ibmcmp/xlmass/bg/4.4/bglib -lmass -lmassv

INCLUDE= -I/opt/ibmmath/include -I/opt/ibmcmp/xlmass/bg/4.4/include -I/gpfs/opt/bgp/fftw/3.3.3-fast/include
ifdef RECALC_ONLY
INSTALL_PATH=~/bin/fast/giem2g_master_rc
else
INSTALL_PATH=~/bin/fast/giem2g_master_mf
endif

endif

ifeq ($(PLATFORM), BGP_TEST)

FC= mpixlf2003_r
F77= mpixlf77_r
FC_Link=mpixlf2003_r
AR=ar

FFLAGS= -qessl -g -O0  -qsmp=omp   -qalias=std  -qmaxmem=-1 -qcheck 
FDEFS=-WF,-DLEGACY_MPI -WF,-DIBM_Bluegene  -WF,-DMPI_TIMER  -WF,-Dinternal_timer  # -WF,-DVOLUME_WEIGHT_QUAD  #  -WF,-DNO_DISPLACEMENT_CURRENTS # 
FOPTS=$(FFLAGS) $(FDEFS)

LIB_FFTW=-L/gpfs/opt/bgp/fftw/3.3.3-stable/lib -lfftw3_mpi  -lfftw3_omp   -lfftw3
FGMRES_PATH=/gpfs/data/kruglyakov/lib/src/FGMRES_original
LIB_BLAS=-L/opt/ibmmath/lib  -lesslsmpbg 
LIB_ADD=-L$(HOME)/lib/lib/lib -lzspmv 

INCLUDE= -I/opt/ibmmath/include  -I/gpfs/opt/bgp/fftw/3.3.3-stable/include

ifdef RECALC_ONLY
INSTALL_PATH=~/bin/test/giem2g_test_rc
else
INSTALL_PATH=~/bin/test/giem2g_test
endif
endif
ifeq ($(PLATFORM), LOMONOSOV)

FC=mpif90
F77=mpif77
FC_Link=mpif90
AR=xiar

FOPTS= -g    -openmp  -i4 -mcmodel=large  -O3    -xHost -prec-div -inline speed -ipo=25 -D LEGACY_MPI 

FFT_QUAD_OURA=1

LIB_FFTW= -L$(HOME)/Soft/lib/fftw3_ompi155_icc/lib  -lfftw3_omp -lfftw3_mpi  -lfftw3 
FGMRES_PATH=/mnt/data/users/dm3/vol6/kruglyakov/Soft/src/fgmres
LIB_BLAS=-L$(HOME)/Soft/lib/OpenBLAS_gcc44/lib -lopenblas_omp 
LIB_ADD=-lgfortran
INCLUDE= -I$(HOME)/Soft/lib/fftw3_ompi155_icc/include
INSTALL_PATH=~/_scratch/bin/GIEM2G/giem2g_master

endif
