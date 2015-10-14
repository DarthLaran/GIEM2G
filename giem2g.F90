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

	USE Timer_Module 

	IMPLICIT NONE
#define debug		 
	INTEGER(4)::PROVIDED,IERROR
	INTEGER::shift_y,shift_x
	INTEGER(4)::comm_old,Np
	INTEGER(4)::buff_off(4),I
	INTEGER(8):: loc_shift,Nf
	INTEGER(MPI_CTL_KIND):: STATUS(MPI_STATUS_SIZE)
	INTEGER::Iz,Ifreq,Ix
	REAL(REALPARM):: zrecv(3),dz
	REAL(REALPARM),POINTER::freqs(:)
	TYPE (RECEIVER_TYPE),POINTER::recvs(:)
	TYPE(IntegralEquation)::int_eq
	TYPE(IntegralEquation)::int_eq2
	TYPE(RC_OPERATOR)::rc_op
	COMPLEX(REALPARM),POINTER::FX(:,:,:,:),FY(:,:,:,:)
	COMPLEX(REALPARM),POINTER::Ea(:,:,:,:),Ha(:,:,:,:)
	COMPLEX(REALPARM),POINTER::Et(:,:,:,:),Ht(:,:,:,:)
	COMPLEX(REALPARM),POINTER::Eprevy(:,:,:,:)
	COMPLEX(REALPARM),POINTER::Eprevx(:,:,:,:)
	COMPLEX(REALPARM)::sig1(3),sig2(3),q1,q2
	INTEGER::N,Nfreq,Nr
	INTEGER::NT,OMP_GET_MAX_THREADS
	LOGICAL::threads_ok
	INTEGER(MPI_CTL_KIND)::wcomm,wsize,me,real_comm
	TYPE (BKG_DATA_TYPE)::bkg
	TYPE (ANOMALY_TYPE)::anomaly
	TYPE (FGMRES_CTL_TYPE)::fgmres_ctl
	REAL(REALPARM)::RC_Threshold,IE_Threshold
	CHARACTER(len=11)::fnum1,fnum2
	CHARACTER(len=1024),POINTER::anom_list(:)
	INTEGER::Istr,Na,Ia
	REAL(8)::time2
	REAL(DOUBLEPARM)::tick
!-------------------MPI INITIALIZATION-------------------------------------!
	CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, PROVIDED, IERROR)
	CALL InitTimer
	wcomm=MPI_COMM_WORLD
	
	CALL MPI_COMM_RANK(wcomm, me, IERROR)
	CALL MPI_COMM_SIZE(wcomm,wsize,IERROR) 
	IF (me==0) THEN
#ifndef MPI_TIMER
		 PRINT '(A)', "Timer is based in SYSTEM_CLOCK"
#else
		 PRINT '(A)', "Timer is based on MPI_Wtime()"
#endif
	ENDIF
	IF ( .NOT. FFTW_THREADS_DISABLE) THEN
		IF (PROVIDED>=MPI_THREAD_FUNNELED) THEN 
			threads_ok=.TRUE.
			NF= fftw_init_threads(); 
			IF ((me==0).AND.(NF/=0)) THEN
				PRINT*, 'FFTW THREADS ENABLE'
			ENDIF			 
		ELSEIF(FFTW_THREADS_FORCE) THEN
			threads_ok=.TRUE.
			NF= fftw_init_threads(); 
			IF ((me==0).AND.(NF/=0)) THEN
				PRINT*, 'FFTW THREADS FORCED ENABLE'
			ENDIF			 
		ELSE
			threads_ok=.FALSE.
			IF (me==0) THEN
				PRINT*, 'FFTW THREADS DISABLE'
			ENDIF			 
		ENDIF
	ELSE	
			threads_ok=.FALSE.
			IF (me==0) THEN
				PRINT*, 'FFTW FORCED THREADS DISABLE'
			ENDIF			 
	ENDIF
		IF (me==0)	PRINT*,'Number of processes:',wsize
		NT=OMP_GET_MAX_THREADS()
		IF (me==0)	PRINT*,'Number of threads:',NT

#ifdef performance_test
	IF (me==0) THEN
		PRINT'(A)','PERFORMACE TEST'
		PRINT'(A)','NO RESULT STORAGE'
	ENDIF
#endif

	CALL  FFTW_MPI_INIT
	freqs=>NULL()
        recvs=>NULL()
	CALL LoadBackground(bkg,wcomm,'background.dat')
	CALL LoadAnomalyShape(anomaly,bkg,wcomm,'anomaly_shape.dat',.TRUE.)
	CALL LoadFrequencies(freqs,wcomm,'frequencies.dat')
	CALL LoadRecievers(recvs,wcomm,'recievers.dat')
	CALL LoadFGMRES_Ctl(fgmres_ctl,wcomm,'fgmres_ctl.dat')
	CALL LoadThreshold(IE_Threshold,RC_Threshold,wcomm,'threshold.dat')
	CALL LoadAnomalySigmaList('anomaly_list.dat',wcomm,anom_list,Na)
        
	Nr=SIZE(recvs)	
	Nfreq=SIZE(freqs)
	N=anomaly%Nz
	
	CALL PrepareRecvs(recvs,anomaly,bkg)
	IF (me==0) PRINT'(A80)','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	IF (me==0) PRINT*, 'Nx=',anomaly%Nx, 'Ny=',anomaly%Ny,'Nz=',anomaly%Nz
	CALL PrepareIntegralEquation(int_eq,anomaly,wcomm,threads_ok,1)
	IF (me==0) PRINT*, 'Number of blocks in async FFT', int_eq%DFD%Nb
	real_comm=int_eq%fgmres_comm

	CALL PrepareContinuationOperator(rc_op,anomaly,recvs,wcomm,threads_ok);!,int_eq%DFD)

	CALL  SetSigb(int_eq,anomaly,bkg)
	CALL SetSigbRC(rc_op,anomaly,bkg)

	IF (int_eq%real_space) THEN
		CALL AllocateSiga(anomaly)
		ALLOCATE(FX(1,1,1,6),FY(1,1,1,6))
		ALLOCATE(Ea(Nr,EX:EZ,int_eq%Nx,int_eq%Ny_loc),Ha(Nr,HX:HZ,int_eq%Nx,int_eq%Ny_loc))
		ALLOCATE(Et(Nr,EX:EZ,int_eq%Nx,int_eq%Ny_loc),Ht(Nr,HX:HZ,int_eq%Nx,int_eq%Ny_loc))
		
		int_eq%siga=>anomaly%siga
		rc_op%siga=>anomaly%siga

	ELSE
				FX=>NULL()
				FY=>NULL();
				Ea=>NULL()
				Ha=>NULL()
				Et=>NULL()
				Ht=>NULL()
	ENDIF

!----------------------------------------------------------------------------!



	DO Ifreq=1,Nfreq
		WRITE (fnum1,'(F11.5)') freqs(Ifreq)
		DO Istr=1,11
				IF (fnum1(Istr:Istr)==' ') fnum1(Istr:Istr)='0'
		ENDDO
		IF (me==0) THEN
			PRINT'(A, F11.5, A)', 'Frequency:', freqs(Ifreq), 'Hz'
		ENDIF
		CALL Set_Freq(bkg,freqs(Ifreq))


		CALL CalcIntegralGreenTensor(int_eq,bkg,anomaly,IE_Threshold)
		CALL CalcFFTofIETensor(int_eq)

        IF (me==0)	PRINT'(I7, 8F25.10)' ,me, int_eq%G_symm(1,:,1,1)
        IF (me==0)	PRINT'(I7, 4F25.10)' ,me, int_eq%G_asym(1,1,:,1,1)

		CALL CalcRecalculationGreenTensor(rc_op,bkg,anomaly,RC_Threshold)
		CALL CalcFFTofRCTensor(rc_op)

		DO Ia=1,Na
			WRITE (fnum2,'(I5.5)') Ia
			IF (me==0)	PRINT*, 'Anomaly ', Ia, ' from ', trim(anom_list(Ia))
			IF (int_eq%real_space) THEN
				CALL PlaneWaveIntegral(EY,bkg,anomaly,int_eq%E_n)
				CALL PlaneWave(EY,bkg,(/recvs(1)%zrecv/),FY)
				IF ((Na/=1).OR.(Ifreq==1))THEN
			                CALL LoadAnomalySigma(anomaly,real_comm,trim(anom_list(Ia)))
				ENDIF
			ENDIF

			IF (me==0) PRINT*, 'FY=', FY


			CALL SolveEquation(int_eq,fgmres_ctl)
			time2=GetTime()
#ifndef performance_test
			IF (int_eq%real_space) THEN
!				CALL SaveIESolutionSeparateBinary(int_eq%Esol,int_eq%me,'SOL_PY_F'//trim(fnum1)//'T_'//trim(fnum2))
				CALL SaveIESolutionOneFIleBinary(int_eq%Esol,real_comm,'SOL_PY_F'//trim(fnum1)//'T_'//trim(fnum2))
			ENDIF
			time2=GetTime()-time2;
			IF (int_eq%master) THEN
				PRINT*, "Save IE Solution in",time2, 's'
			ENDIF
#endif
			CALL ReCalculation(rc_op,int_eq%Esol,Ea,Ha)
#ifndef performance_test
			time2=GetTime()
			IF (int_eq%real_space) THEN

				Et(:,EX,:,:)=Ea(:,EX,:,:)+FY(1,1,1,EX)
				Et(:,EY,:,:)=Ea(:,EY,:,:)+FY(1,1,1,EY)
				Et(:,EZ,:,:)=Ea(:,EZ,:,:)+FY(1,1,1,EZ)

				Ht(:,HX,:,:)=Ha(:,HX,:,:)+FY(1,1,1,HX)
				Ht(:,HY,:,:)=Ha(:,HY,:,:)+FY(1,1,1,HY)
				Ht(:,HZ,:,:)=Ha(:,HZ,:,:)+FY(1,1,1,HZ)

				CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PY_F'//trim(fnum1)//'T_'//trim(fnum2))
				!CALL SaveOutputSeparate(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),rc_op%me,rc_op%Ny_offset,'PY_F'//trim(fnum1)//'T_'//trim(fnum2))
			ENDIF
			time2=GetTime()-time2;
			IF (int_eq%master) THEN
				PRINT*, "Save fields in",time2, 's'
			ENDIF

#endif
			IF (int_eq%real_space) THEN
				CALL PlaneWaveIntegral(EX,bkg,anomaly,int_eq%E_n)
				CALL PlaneWave(EX,bkg,(/recvs(1)%zrecv/),FX)
			ENDIF

			CALL SolveEquation(int_eq,fgmres_ctl)
#ifndef performance_test
			time2=GetTime()
			IF (int_eq%real_space) THEN
				!CALL SaveIESolutionSeparateBinary(int_eq%Esol,int_eq%me,'SOL_PX_F'//trim(fnum1)//'T_'//trim(fnum2))
				CALL SaveIESolutionOneFIleBinary(int_eq%Esol,real_comm,'SOL_PX_F'//trim(fnum1)//'T_'//trim(fnum2))
			ENDIF
			time2=GetTime()-time2;
			IF (int_eq%master) THEN
				PRINT*, "Save IE Solution in",time2, 's'
			ENDIF
#endif
			CALL ReCalculation(rc_op,int_eq%Esol,Ea,Ha)

#ifndef performance_test
			time2=GetTime()
			IF (int_eq%real_space) THEN

				Et(:,EX,:,:)=Ea(:,EX,:,:)+FX(1,1,1,EX)
				Et(:,EY,:,:)=Ea(:,EY,:,:)+FX(1,1,1,EY)
				Et(:,EZ,:,:)=Ea(:,EZ,:,:)+FX(1,1,1,EZ)

				Ht(:,HX,:,:)=Ha(:,HX,:,:)+FX(1,1,1,HX)
				Ht(:,HY,:,:)=Ha(:,HY,:,:)+FX(1,1,1,HY)
				Ht(:,HZ,:,:)=Ha(:,HZ,:,:)+FX(1,1,1,HZ)
!				CALL SaveOutputSeparate(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),rc_op%me,rc_op%Ny_offset,'PX_F'//trim(fnum1)//'T_'//trim(fnum2))

				CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PX_F'//trim(fnum1)//'T_'//trim(fnum2))
			ENDIF
			time2=GetTime()-time2;
			IF (int_eq%master) THEN
				PRINT*, "Save fields in",time2, 's'
			ENDIF
#endif
		ENDDO
	ENDDO
	CALL CHECK_MEM(int_eq%me,int_eq%master_proc,int_eq%matrix_comm)
	CALL DeleteIE_OP(int_eq)
	CALL DeleteRC_OP(RC_OP)
	CALL CHECK_MEM(int_eq%me,int_eq%master_proc,int_eq%matrix_comm)
	CALL MPI_FINALIZE(IERROR)
END PROGRAM
