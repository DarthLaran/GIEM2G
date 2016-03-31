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
INTEGER(MPI_CTL_KIND)::PROVIDED,IERROR
INTEGER(MPI_CTL_KIND):: STATUS(MPI_STATUS_SIZE)
INTEGER::Iz,Ifreq
REAL(REALPARM),POINTER::freqs(:)
TYPE (RECEIVER_TYPE),POINTER::recvs(:)
TYPE(IntegralEquationOperator)::ie_op
TYPE(RC_OPERATOR)::rc_op
COMPLEX(REALPARM),POINTER::FX(:,:,:,:),FY(:,:,:,:)
COMPLEX(REALPARM),POINTER::Ea(:,:,:,:),Ha(:,:,:,:)
COMPLEX(REALPARM),POINTER::Et(:,:,:,:),Ht(:,:,:,:)
COMPLEX(REALPARM),POINTER::E_bkg(:,:,:,:)
COMPLEX(REALPARM),POINTER::E_sol(:,:,:,:)
COMPLEX(REALPARM),POINTER::Eprevy(:,:,:,:)
COMPLEX(REALPARM),POINTER::Eprevx(:,:,:,:)
REAL(REALPARM),POINTER::Stations(:,:)
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
INTEGER(MPI_CTL_KIND)::wcomm,wsize,me,real_comm
TYPE (BKG_DATA_TYPE)::bkg
TYPE (ANOMALY_TYPE)::anomaly
TYPE (FGMRES_CTL_TYPE)::fgmres_ctl

REAL(DOUBLEPARM)::time2
REAL(DOUBLEPARM),ALLOCATABLE::recv_depths(:)

CHARACTER(len=22)::fnum1,fnum2
CHARACTER(len=1024),POINTER::anom_list(:)
INTEGER::Istr,Na,Ia,Ir,Nf
CHARACTER(LEN=2048,KIND=C_CHAR)::message 
INTEGER(C_INTPTR_T)::kernel_len,fft_len
TYPE(C_PTR)::p1,p2,p3
INTEGER::ll
!-------------------MPI INITIALIZATION-------------------------------------!
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, PROVIDED, IERROR)
CALL InitTimer
LOGGER=>NATIVE_LOGGER

wcomm=MPI_COMM_WORLD

CALL MPI_COMM_RANK(wcomm, me, IERROR)
CALL MPI_COMM_SIZE(wcomm,wsize,IERROR) 
LOGGER_MASTER=me
#ifndef MPI_TIMER
CALL LOGGER("Timer is based in SYSTEM_CLOCK")
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


SAVE_SOLUTION=.TRUE.

SOLVE_EQUATION=.TRUE.
RECALC_FIELD=.TRUE.

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
CALL LoadBackground(bkg,wcomm,'background.dat')
CALL LoadAnomalyShape(anomaly,bkg,wcomm,'anomaly_shape.dat',.TRUE.)
CALL LoadFrequencies(freqs,wcomm,'frequencies.dat')
CALL LoadRecievers(recvs,wcomm,'recievers.dat')
CALL LoadFGMRES_Ctl(fgmres_ctl,wcomm,'fgmres_ctl.dat')
CALL LoadAnomalySigmaList('anomaly_list.dat',wcomm,anom_list,Na)

Nr=SIZE(recvs)	
Nfreq=SIZE(freqs)
N=anomaly%Nz
ALLOCATE(recv_depths(Nr))	
CALL PrepareRecvs(recvs,anomaly,bkg)
DO Ir=1,Nr
	recv_depths(Ir)=recvs(Ir)%zrecv
	IF (me==0) THEN
	    PRINT*,recvs(Ir)%anom_cell,recvs(Ir)%zrecv
	ENDIF
ENDDO
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
	CALL PrepareContinuationOperator(rc_op,anomaly,recvs,wcomm,threads_ok);

	IF (rc_op%real_space) THEN
		     CALL  MPI_COMM_SPLIT(rc_op%matrix_comm, MPI_ONE, rc_op%me,    real_comm, IERROR)
	ELSE
		  CALL   MPI_COMM_SPLIT(rc_op%matrix_comm, MPI_TWO, rc_op%me,   real_comm, IERROR)
	ENDIF
	IF (me==0) THEN
	    CALL LoadStations('stations.dat',stations)
	ENDIF
	mz=anomaly%dz(1)*0.6;
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
	DO Ia=1,Na
		WRITE (fnum2,'(I5.5)') Ia
		WRITE (message,'(A, I4, 2A )') 'Anomaly ', Ia, ' from ', trim(anom_list(Ia))
		CALL LOGGER(message)
		IF (SOLVE_EQUATION) THEN
			IF (ie_op%real_space) THEN
				CALL PlaneWaveIntegral(EY,bkg,anomaly,E_bkg)
				IF ((Na/=1).OR.(Ifreq==1))THEN
					CALL LoadAnomalySigma(anomaly,real_comm,trim(anom_list(Ia)))
				ENDIF
			ENDIF
			CALL SetAnomalySigma(ie_op,anomaly%siga,freqs(Ifreq))
			CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
			time2=GetTime()
			success=.TRUE.
			IF (ie_op%real_space) THEN
				CALL LoadIESolutionOneFIleBinary(Eprevy,real_comm,'SOL_PY_F'//trim(fnum1)//'T_'//trim(fnum2),success)
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
			Eprevy=E_sol

			CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
			IF (SAVE_SOLUTION) THEN
				time2=GetTime()
				IF (ie_op%real_space) THEN
					CALL SaveIESolutionOneFIleBinary(E_sol,real_comm,'SOL_PY_F'//trim(fnum1)//'T_'//trim(fnum2))
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
				IF ((Na/=1).OR.(Ifreq==1))THEN
					CALL LoadAnomalySigma(anomaly,real_comm,trim(anom_list(Ia)))
				ENDIF
			ENDIF
			    CALL FinalizeRCOperator(rc_op,bkg,anomaly,freqs(Ifreq))

			time2=GetTime()
			success=.TRUE.
			IF (rc_op%real_space) THEN
				CALL LoadIESolutionOneFIleBinary(Eprevy,real_comm,'SOL_PY_F'//trim(fnum1)//'T_'//trim(fnum2),success)
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
			CALL PlaneWave(EY,bkg,recv_depths,FY)
			time2=GetTime()
			IF (rc_op%real_space) THEN
			    DO Ir=1,Nr

				    Et(Ir,EX,:,:)=Ea(Ir,EX,:,:)+FY(Ir,1,1,EX)
				    Et(Ir,EY,:,:)=Ea(Ir,EY,:,:)+FY(Ir,1,1,EY)
				    Et(Ir,EZ,:,:)=Ea(Ir,EZ,:,:)+FY(Ir,1,1,EZ)

				    Ht(Ir,HX,:,:)=Ha(Ir,HX,:,:)+FY(Ir,1,1,HX)
				    Ht(Ir,HY,:,:)=Ha(Ir,HY,:,:)+FY(Ir,1,1,HY)
				    Ht(Ir,HZ,:,:)=Ha(Ir,HZ,:,:)+FY(Ir,1,1,HZ)
			    ENDDO
			    !CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PY_F'//trim(fnum1)//'T_'//trim(fnum2))
			    CALL SaveOutputOneFileStations(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PY_F'//trim(fnum1)//'T_'//trim(fnum2),stations,mz)
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
					CALL LoadIESolutionOneFIleBinary(Eprevx,real_comm,'SOL_PX_F'//trim(fnum1)//'T_'//trim(fnum2),success)
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
				Eprevx=E_sol
				CALL MPI_BARRIER(ie_op%ie_comm,IERROR)
				IF (SAVE_SOLUTION) THEN
					time2=GetTime()
					IF (ie_op%real_space) THEN
						CALL SaveIESolutionOneFIleBinary(E_sol,real_comm,'SOL_PX_F'//trim(fnum1)//'T_'//trim(fnum2))
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
				CALL LoadIESolutionOneFIleBinary(Eprevx,real_comm,'SOL_PX_F'//trim(fnum1)//'T_'//trim(fnum2),success)
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
			CALL PlaneWave(EX,bkg,recv_depths,FX)
			time2=GetTime()
			IF (rc_op%real_space) THEN
				DO Ir=1,Nr
					Et(Ir,EX,:,:)=Ea(Ir,EX,:,:)+FX(Ir,1,1,EX)
					Et(Ir,EY,:,:)=Ea(Ir,EY,:,:)+FX(Ir,1,1,EY)
					Et(Ir,EZ,:,:)=Ea(Ir,EZ,:,:)+FX(Ir,1,1,EZ)

					Ht(Ir,HX,:,:)=Ha(Ir,HX,:,:)+FX(Ir,1,1,HX)
					Ht(Ir,HY,:,:)=Ha(Ir,HY,:,:)+FX(Ir,1,1,HY)
					Ht(Ir,HZ,:,:)=Ha(Ir,HZ,:,:)+FX(Ir,1,1,HZ)
				ENDDO
!			CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PX_F'//trim(fnum1)//'T_'//trim(fnum2))
				CALL SaveOutputOneFileStations(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PX_F'//trim(fnum1)//'T_'//trim(fnum2),stations,mz)
			ENDIF
			time2=GetTime()-time2;
			CALL PRINT_CALC_TIME("Save fields in",time2) 
		    ENDIF
	ENDDO
ENDDO
CALL CHECK_MEM(me,0,wcomm)
IF (SOLVE_EQUATION) CALL DeleteIE_OP(ie_op)
IF (RECALC_FIELD) CALL DeleteRC_OP(RC_OP)
CALL LOGGER("Final!")
CALL MPI_FINALIZE(IERROR)
END PROGRAM
