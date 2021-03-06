!Copyright (c) 2016 Mikhail Kruglyakov 
!This file is part of GIEM2G.
!
!GIEM2G is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 2 of the License.
!
!GIEM2G is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with GIEM2G.  If not, see <http://www.gnu.org/licenses/>.
!
!
!

PROGRAM GIEMIEMG
USE CONST_MODULE
USE FFTW3
USE MPI_MODULE
USE DATA_TYPES_MODULE
USE MPI_SAVELOAD_MODULE
USE SOURCES_MODULE
USE INTEGRAL_EQUATION_MODULE
USE Calc_IE_Tensor_Module
USE IE_SOLVER_MODULE

USE CONTINUATION_FUNCTION_MODULE
USE Calc_RC_Tensor_Module
USE CHECK_MEMORY

USE LOGGER_MODULE
USE Timer_Module 

IMPLICIT NONE
#define debug		 

CHARACTER(len=*),PARAMETER::default_input="giem2g_input.json"

INTEGER(MPI_CTL_KIND)::PROVIDED,IERROR
INTEGER(MPI_CTL_KIND):: STATUS(MPI_STATUS_SIZE)
INTEGER::Iz,Ifreq
REAL(REALPARM),POINTER::freqs(:)=>NULL()
TYPE (RECEIVER_TYPE),POINTER::recvs(:)
TYPE(IntegralEquationOperator)::ie_op
TYPE(RC_OPERATOR)::rc_op
COMPLEX(REALPARM),POINTER::FX(:,:,:,:),FY(:,:,:,:)
COMPLEX(REALPARM),POINTER::Ea(:,:,:,:),Ha(:,:,:,:)
COMPLEX(REALPARM),POINTER::Et(:,:,:,:),Ht(:,:,:,:)
COMPLEX(REALPARM),POINTER::E_bkg(:,:,:,:)=>NULL()
COMPLEX(REALPARM),POINTER::E_sol(:,:,:,:)=>NULL()

COMPLEX(REALPARM),POINTER::Eprevy(:,:,:,:)=>NULL()

COMPLEX(REALPARM),POINTER::Eprevx(:,:,:,:)=>NULL()

REAL(REALPARM),POINTER::Stations(:,:)=>NULL()

REAL(REALPARM)::mz
INTEGER(MPI_CTL_KIND),PARAMETER::MPI_TWO=2
INTEGER(MPI_CTL_KIND),PARAMETER::MPI_ONE=1
INTEGER::N,Nfreq,Nr
INTEGER::NT
LOGICAL::threads_ok

LOGICAL::SAVE_SOLUTION
LOGICAL::SOLVE_EQUATION
LOGICAL::RECALC_FIELD

LOGICAL::success
LOGICAL::save_all_output
INTEGER(MPI_CTL_KIND)::wcomm,wsize,me,real_comm
TYPE (BKG_DATA_TYPE)::bkg
TYPE (ANOMALY_TYPE)::anomaly
TYPE (FGMRES_CTL_TYPE)::fgmres_ctl

TYPE(FSON_VALUE), POINTER :: input_data
TYPE(FSON_VALUE), POINTER :: test_item
REAL(DOUBLEPARM)::time2
REAL(DOUBLEPARM),ALLOCATABLE::recv_depths(:)

CHARACTER(len=22)::fnum1,fnum2
CHARACTER(len=1024),POINTER::anom_list(:)=>NULL()

INTEGER::Istr,Ir,Nf
CHARACTER(LEN=2048,KIND=C_CHAR)::message 

CHARACTER(LEN=2048)::input_file

INTEGER(C_INTPTR_T)::kernel_len,fft_len
TYPE(C_PTR)::p1,p2,p3
INTEGER::ll
!-------------------MPI INITIALIZATION-------------------------------------!
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, PROVIDED, IERROR)
CALL InitTimer
LOGGER=>NATIVE_LOGGER

wcomm=MPI_COMM_WORLD

!CALL MPI_ERRHANDLER_SET(wcomm,MPI_ERRORS_RETURN, IERROR)
CALL MPI_COMM_SET_ERRHANDLER(wcomm,MPI_ERRORS_RETURN, IERROR)

CALL MPI_COMM_RANK(wcomm, me, IERROR)
CALL MPI_COMM_SIZE(wcomm,wsize,IERROR) 
LOGGER_MASTER=me
#ifndef MPI_TIMER
CALL LOGGER("Timer is based on SYSTEM_CLOCK")
#else
CALL LOGGER("Timer is based on MPI_Wtime()")
#endif
#ifndef VOLUME_WEIGHT_QUAD
CALL LOGGER("Convolution weights are computed using double	 precision")
#else
CALL LOGGER("Convolution weights are computed using quadruple	 precision")
#endif

#ifndef NO_DISPLACEMENT_CURRENTS 
CALL LOGGER("Displacement currents are modelled")
#else
CALL LOGGER("Displacement currents are ignored")
#endif

IF ( .NOT. FFTW_THREADS_DISABLE) THEN
	IF (PROVIDED>=MPI_THREAD_FUNNELED) THEN 
		threads_ok=.TRUE.
		NF= fftw_init_threads(); 
		IF (NF/=0) THEN
			CALL LOGGER('FFTW THREADS ENABLE')
		ELSE
			threads_ok=.FALSE.
			CALL LOGGER("FFTW THREADS DISABLE")
		ENDIF			
		ELSEIF(FFTW_THREADS_FORCE) THEN
		threads_ok=.TRUE.
		NF= fftw_init_threads(); 
		IF (NF/=0) THEN
			CALL LOGGER('FFTW THREADS FORCED ENABLE')
		ELSE
			threads_ok=.FALSE.
			CALL LOGGER("FFTW THREADS DISABLE")
		ENDIF			
	ELSE
		threads_ok=.FALSE.
		CALL LOGGER("FFTW THREADS DISABLE")
	ENDIF
ELSE	
	threads_ok=.FALSE.
	CALL LOGGER ("FFTW FORCED THREADS DISABLE")
ENDIF
CALL PRINT_CALC_NUMBER('Number of processes:',wsize)
NT=OMP_GET_MAX_THREADS()
CALL PRINT_CALC_NUMBER('Number of threads:',NT)


ll=COMMAND_ARGUMENT_COUNT()

if (ll==0) THEN
        input_file=default_input
else
        CALL GET_COMMAND_ARGUMENT(1,input_file)
endif
        input_data => fson_parse(trim(input_file))
CALL Logger(" ")
CALL Logger("INPUT FILE: "//trim(input_file))
CALL Logger(" ")

SAVE_SOLUTION=.TRUE.
SOLVE_EQUATION=.TRUE.
RECALC_FIELD=.TRUE.
save_all_output=.TRUE.

CALL FSON_GET(input_data, "Save Solution", test_item)
IF (ASSOCIATED(test_item)) CALL FSON_GET(test_item, "", SAVE_SOLUTION)

CALL FSON_GET(input_data, "Solve Equation", test_item)
IF (ASSOCIATED(test_item)) CALL FSON_GET(test_item, "", SOLVE_EQUATION)

CALL FSON_GET(input_data, "Recalculate Fields", test_item)
IF (ASSOCIATED(test_item)) CALL FSON_GET(test_item, "", RECALC_FIELD)



#ifdef performance_test
CALL LOGGER('PERFORMACE TEST')
CALL LOGGER('NO CONTINUATION TO RECIEVERS')
CALL LOGGER('NO RESULT STORAGE')
SAVE_SOLUION=.FALSE.
#endif


IF (SAVE_SOLUTION) THEN

	CALL LOGGER("IE solution will be save to disk")
ELSE
	CALL LOGGER("IE solution will NOT be stored")
ENDIF

IF (SOLVE_EQUATION) CALL LOGGER("IE equation will be solved")
IF (RECALC_FIELD) CALL LOGGER("Electrical field will be continued to recievers")



CALL  FFTW_MPI_INIT
freqs=>NULL()
recvs=>NULL()
!CALL LoadBackground(bkg,wcomm,'background.dat')
!CALL LoadAnomalyShape(anomaly,bkg,wcomm,'anomaly_shape.dat',.TRUE.)
!CALL LoadFrequencies(freqs,wcomm,'frequencies.dat')
!CALL LoadFGMRES_Ctl(fgmres_ctl,wcomm,'fgmres_ctl.dat')
!CALL LoadAnomalySigmaList('anomaly_list.dat',wcomm,anom_list,Na)

CALL LoadBackground(bkg,wcomm,input_data)
CALL LoadAnomalyShape(anomaly,bkg,wcomm,input_data)
CALL LoadFrequencies(freqs,wcomm,input_data)
CALL LoadFGMRES_Ctl(fgmres_ctl,wcomm,input_data)
Nr=SIZE(recvs)	
Nfreq=SIZE(freqs)
N=anomaly%Nz
CALL PRINT_BORDER
WRITE (message,'(A, I5, A, I5, A, I5)') 'Nx=',anomaly%Nx, ' Ny=',anomaly%Ny,' Nz=',anomaly%Nz
CALL LOGGER(message)

#ifdef performance_test
SAVE_SOLUION=.FALSE.
#endif



IF (SOLVE_EQUATION) THEN

	kernel_len=CalcLocalKernelSize(anomaly,wsize) 

	CALL PRINT_STORAGE_SIZE('IE matrix needs', kernel_len*16.0_REALPARM/1024/1024/1024)

	fft_len=CalcLocalFFTSize(2*anomaly%Nx,2*anomaly%Ny,3*anomaly%Nz,wsize)
	p1=ALLOCATE_BUFF(kernel_len)
	p2=ALLOCATE_BUFF(fft_len)
	p3=ALLOCATE_BUFF(fft_len)

	CALL PrepareIntegralEquation(ie_op,anomaly,wcomm,p1,p2,p3,kernel_len,fft_len,threads_ok)
	IF (ie_op%real_space) THEN
		     CALL  MPI_COMM_SPLIT(ie_op%ie_comm, MPI_ONE, ie_op%me,    real_comm, IERROR)
	ELSE
		  CALL   MPI_COMM_SPLIT(ie_op%ie_comm, MPI_TWO, ie_op%me,   real_comm, IERROR)
	ENDIF
	IF (ie_op%real_space) THEN
		anomaly%Ny_loc=ie_op%Ny_loc
		CALL AllocateSiga(anomaly)
		ALLOCATE(E_bkg(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz,3))
		ALLOCATE(E_sol(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz,3))
		ALLOCATE(Eprevx(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz,3))
		ALLOCATE(Eprevy(ie_op%Nx,ie_op%Ny_loc,ie_op%Nz,3))
		FX=>NULL()
		FY=>NULL();
		Ea=>NULL()
		Ha=>NULL()
		Et=>NULL()
		Ht=>NULL()
	ELSE
		FX=>NULL()
		FY=>NULL();
		Ea=>NULL()
		Ha=>NULL()
		Et=>NULL()
		Ht=>NULL()
		E_sol=>NULL()
		E_bkg=>NULL()
		Eprevx=>NULL()
		Eprevy=>NULL()
	ENDIF
ENDIF


IF (RECALC_FIELD) THEN

	CALL LoadRecievers(recvs,wcomm,input_data)
	Nr=SIZE(recvs)	
	ALLOCATE(recv_depths(Nr))	
	CALL PrepareRecvs(recvs,anomaly,bkg)
	DO Ir=1,Nr
		recv_depths(Ir)=recvs(Ir)%zrecv
	ENDDO
	CALL PrepareContinuationOperator(rc_op,anomaly,recvs,wcomm,threads_ok);

	IF (rc_op%real_space) THEN
		     CALL  MPI_COMM_SPLIT(rc_op%matrix_comm, MPI_ONE, rc_op%me,    real_comm, IERROR)
	ELSE
		  CALL   MPI_COMM_SPLIT(rc_op%matrix_comm, MPI_TWO, rc_op%me,   real_comm, IERROR)
	ENDIF

	mz=anomaly%dz(1)*1.1;
	IF (rc_op%real_space) THEN
		anomaly%Ny_loc=rc_op%Ny_loc
		CALL AllocateSiga(anomaly)
		ALLOCATE(FX(Nr,1,1,6),FY(Nr,1,1,6))

		ALLOCATE(Ea(Nr,EX:EZ,rc_op%Nx,rc_op%Ny_loc),Ha(Nr,HX:HZ,rc_op%Nx,rc_op%Ny_loc))
		ALLOCATE(Et(Nr,EX:EZ,rc_op%Nx,rc_op%Ny_loc),Ht(Nr,HX:HZ,rc_op%Nx,rc_op%Ny_loc))

		IF (.NOT. ASSOCIATED(E_sol) ) ALLOCATE(E_sol(rc_op%Nx,rc_op%Ny_loc,rc_op%Nz,3))
		IF (.NOT. ASSOCIATED(Eprevx)) ALLOCATE(Eprevx(rc_op%Nx,rc_op%Ny_loc,rc_op%Nz,3))
		IF (.NOT. ASSOCIATED(Eprevy)) ALLOCATE(Eprevy(rc_op%Nx,rc_op%Ny_loc,rc_op%Nz,3))
	ELSE
	    FX=>NULL()
	    FY=>NULL();
	    Ea=>NULL()
	    Ha=>NULL()
	    Et=>NULL()
	    Ht=>NULL()
	    E_sol=>NULL()
	    E_bkg=>NULL()
	    Eprevx=>NULL()
	    Eprevy=>NULL()
	ENDIF

ENDIF

CALL LOGGER("After allocating memory for kernels")
CALL CHECK_MEM(me,0,wcomm)
!----------------------------------------------------------------------------!
DO Ifreq=1,Nfreq
	WRITE (fnum1,'(F22.10)') freqs(Ifreq)

	DO Istr=1,11
		IF (fnum1(Istr:Istr)==' ') fnum1(Istr:Istr)='0'
	ENDDO


	WRITE (message,'(A, F22.10, A)') 'Frequency:', freqs(Ifreq), 'Hz'
	CALL LOGGER(message)

	CALL Set_Freq(bkg,freqs(Ifreq))
	IF (SOLVE_EQUATION) THEN
		CALL CalcIntegralGreenTensor(ie_op,bkg,anomaly)
		CALL CalcFFTofIETensor(ie_op)
       ENDIF

	IF (RECALC_FIELD) THEN
	    CALL CalcRecalculationGreenTensor(rc_op,bkg,anomaly)
	    CALL CalcFFTofRCTensor(rc_op)
       ENDIF
       IF (SOLVE_EQUATION) THEN
                IF (ie_op%real_space) THEN
                        CALL PlaneWaveIntegral(EY,bkg,anomaly,E_bkg)
                        IF (Ifreq==1) THEN
                                CALL LoadAnomalySigma(anomaly,real_comm,input_data)
                        ENDIF
                ENDIF
                CALL SetAnomalySigma(ie_op,anomaly%siga,freqs(Ifreq))
                CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
                time2=GetTime()
                success=.TRUE.
                IF (ie_op%real_space) THEN
                        CALL LoadIESolutionOneFIleBinary(Eprevy,real_comm,'SOL_PY_F'//trim(fnum1),success)
                ENDIF
                time2=GetTime()-time2;
                CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
                CALL MPI_BCAST(success, 1, MPI_LOGICAL, 0,wcomm, IERROR)
                IF (success) THEN
                        CALL PRINT_CALC_TIME("Load IE Solution in",time2)
                ELSE
                    CALL LOGGER("Solution guess not found")
                ENDIF
                IF (success) THEN
                                CALL SolveEquation(ie_op,fgmres_ctl,E_bkg,E_sol,Eprevy)
                ELSEIF (Ifreq==1) THEN
                        CALL SolveEquation(ie_op,fgmres_ctl,E_bkg,E_sol)
                ELSE
                        CALL SolveEquation(ie_op,fgmres_ctl,E_bkg,E_sol,Eprevy)
                ENDIF
                IF (ie_op%real_space)	Eprevy=E_sol

                CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
                IF (SAVE_SOLUTION) THEN
                        time2=GetTime()
                        IF (ie_op%real_space) THEN
                                CALL SaveIESolutionOneFIleBinary(E_sol,real_comm,'SOL_PY_F'//trim(fnum1))
                        ENDIF
                        time2=GetTime()-time2;
                        CALL PRINT_CALC_TIME("Save IE Solution in",time2) 
                ELSE
                        time2=GetTime()-time2;
                ENDIF
                IF (RECALC_FIELD) THEN
                    rc_op%csigb=>ie_op%csigb
                    rc_op%csiga=>ie_op%csiga
                ENDIF
            ELSEIF (RECALC_FIELD) THEN
                IF (rc_op%real_space) THEN
                        IF (Ifreq==1)THEN
                                CALL LoadAnomalySigma(anomaly,real_comm,input_data)
                        ENDIF
                ENDIF
                    CALL FinalizeRCOperator(rc_op,bkg,anomaly,freqs(Ifreq))

                time2=GetTime()
                success=.TRUE.
                IF (rc_op%real_space) THEN
                        CALL LoadIESolutionOneFIleBinary(Eprevy,real_comm,'SOL_PY_F'//trim(fnum1),success)
                ENDIF
                time2=GetTime()-time2;
                CALL MPI_BCAST(success, 1, MPI_LOGICAL, 0,wcomm, IERROR)
                IF (success) THEN
                        CALL PRINT_CALC_TIME("Load IE Solution in",time2)
                ELSE
                    CALL LOGGER("Solution  not found")
                    CYCLE
                ENDIF
            ENDIF
            IF (RECALC_FIELD) THEN
                CALL ReCalculation(rc_op,Eprevy,Ea,Ha)
                time2=GetTime()
                IF (rc_op%real_space) THEN
                        CALL PlaneWave(EY,bkg,recv_depths,FY)
                    DO Ir=1,Nr

                            Et(Ir,EX,:,:)=Ea(Ir,EX,:,:)+FY(Ir,1,1,EX)
                            Et(Ir,EY,:,:)=Ea(Ir,EY,:,:)+FY(Ir,1,1,EY)
                            Et(Ir,EZ,:,:)=Ea(Ir,EZ,:,:)+FY(Ir,1,1,EZ)

                            Ht(Ir,HX,:,:)=Ha(Ir,HX,:,:)+FY(Ir,1,1,HX)
                            Ht(Ir,HY,:,:)=Ha(Ir,HY,:,:)+FY(Ir,1,1,HY)
                            Ht(Ir,HZ,:,:)=Ha(Ir,HZ,:,:)+FY(Ir,1,1,HZ)
                    ENDDO
                        CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),&
                                &real_comm,'PY_F'//trim(fnum1))
                ENDIF
                time2=GetTime()-time2;
                CALL PRINT_CALC_TIME("Save fields in",time2) 
            ENDIF
            IF (SOLVE_EQUATION) THEN
                        IF (ie_op%real_space) THEN
                                CALL PlaneWaveIntegral(EX,bkg,anomaly,E_bkg)
                        ENDIF
                        success=.TRUE.
                        IF (ie_op%real_space) THEN
                                CALL LoadIESolutionOneFIleBinary(Eprevx,real_comm,&
                                        &'SOL_PX_F'//trim(fnum1),success)
                        ENDIF
                        time2=GetTime()-time2;
                        CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
                        CALL MPI_BCAST(success, 1, MPI_LOGICAL, 0,wcomm, IERROR)
                        IF (success) THEN
                                CALL PRINT_CALC_TIME("Load IE Solution in",time2)
                        ELSE
                            CALL LOGGER("Solution guess not found")
                        ENDIF
                        IF (success) THEN
                                CALL SolveEquation(ie_op,fgmres_ctl,E_bkg,E_sol,Eprevx)
                        ELSEIF (Ifreq==1) THEN
                                CALL SolveEquation(ie_op,fgmres_ctl,E_bkg,E_sol)
                        ELSE
                                CALL SolveEquation(ie_op,fgmres_ctl,E_bkg,E_sol,Eprevx)
                        ENDIF
                        IF (ie_op%real_space)	Eprevx=E_sol
                        CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
                        IF (SAVE_SOLUTION) THEN
                                time2=GetTime()
                                IF (ie_op%real_space) THEN
                                        CALL SaveIESolutionOneFIleBinary(E_sol,real_comm,&
                                                &'SOL_PX_F'//trim(fnum1))
                                ENDIF
                                time2=GetTime()-time2;
                                CALL PRINT_CALC_TIME("Save IE Solution in",time2) 
                        ELSE
                                time2=GetTime()-time2;
                        ENDIF
            ELSEIF (RECALC_FIELD) THEN
                time2=GetTime()
                success=.TRUE.
                IF (rc_op%real_space) THEN
                        CALL LoadIESolutionOneFIleBinary(Eprevx,real_comm,&
                                &'SOL_PX_F'//trim(fnum1),success)
                ENDIF
                time2=GetTime()-time2;
                CALL MPI_BCAST(success, 1, MPI_LOGICAL, 0,wcomm, IERROR)
                IF (success) THEN
                        CALL PRINT_CALC_TIME("Load IE Solution in",time2)
                ELSE
                    CALL LOGGER("Solution  not found")
                    CYCLE
                ENDIF
            ENDIF

            IF (RECALC_FIELD) THEN
                CALL ReCalculation(rc_op,Eprevx,Ea,Ha)
                time2=GetTime()
                IF (rc_op%real_space) THEN
                        CALL PlaneWave(EX,bkg,recv_depths,FX)
                        DO Ir=1,Nr
                                Et(Ir,EX,:,:)=Ea(Ir,EX,:,:)+FX(Ir,1,1,EX)
                                Et(Ir,EY,:,:)=Ea(Ir,EY,:,:)+FX(Ir,1,1,EY)
                                Et(Ir,EZ,:,:)=Ea(Ir,EZ,:,:)+FX(Ir,1,1,EZ)

                                Ht(Ir,HX,:,:)=Ha(Ir,HX,:,:)+FX(Ir,1,1,HX)
                                Ht(Ir,HY,:,:)=Ha(Ir,HY,:,:)+FX(Ir,1,1,HY)
                                Ht(Ir,HZ,:,:)=Ha(Ir,HZ,:,:)+FX(Ir,1,1,HZ)
                        ENDDO
                        CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,&
                                &'PX_F'//trim(fnum1))
                ENDIF
                time2=GetTime()-time2;
                CALL PRINT_CALC_TIME("Save fields in",time2) 
            ENDIF
ENDDO

CALL LOGGER("Before exit")
CALL CHECK_MEM(me,0,wcomm)
IF (SOLVE_EQUATION) THEN
         CALL DeleteIE_OP(ie_op)
        CALL fftw_free(p1)
        CALL fftw_free(p2)
        CALL fftw_free(p3)
ENDIF
IF (RECALC_FIELD) CALL DeleteRC_OP(RC_OP)
CALL FSON_DESTROY(input_data)
CALL LOGGER("Finish!")
CALL MPI_FINALIZE(IERROR)
END PROGRAM
