PROGRAM GIEMIEMG
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE DATA_TYPES_MODULE
	USE MPI_MODULE
	USE MPI_SAVELOAD_MODULE
	USE SOURCES_MODULE
	USE FFTW3
   USE INTEGRAL_EQUATION_MODULE
	USE Calc_IE_Tensor_Module
	USE IE_SOLVER_MODULE

	USE CONTINUATION_FUNCTION_MODULE
	USE Calc_RC_Tensor_Module


	IMPLICIT NONE
#define debug		 
#define no_compile
	INTEGER(4)::PROVIDED,IERROR
	INTEGER::shift_y,shift_x
	INTEGER(4)::comm_old,Np
	INTEGER(4)::buff_off(4),I
	INTEGER(8):: loc_shift,Nf
	INTEGER(4):: STATUS(MPI_STATUS_SIZE)
	INTEGER::Iz,Ifreq,Ix
	REAL(REALPARM):: zrecv(3),dz
	REAL(REALPARM),POINTER::freqs(:)
	TYPE (RECEIVER_TYPE),POINTER::recvs(:)
	TYPE(IntegralEquation)::int_eq
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
	INTEGER::wcomm,wsize,me,real_comm
	TYPE (BKG_DATA_TYPE)::bkg
	TYPE (ANOMALY_TYPE)::anomaly
	TYPE (FGMRES_CTL_TYPE)::fgmres_ctl
	REAL(REALPARM)::RC_Threshold,IE_Threshold
	CHARACTER(len=11)::fnum1,fnum2
	CHARACTER(len=1024),POINTER::anom_list(:)
	INTEGER::Istr,Na,Ia
!-------------------MPI INITIALIZATION-------------------------------------!
	CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, PROVIDED, IERROR)
	wcomm=MPI_COMM_WORLD
	
	CALL MPI_COMM_RANK(wcomm, me, IERROR)
	CALL MPI_COMM_SIZE(wcomm,wsize,IERROR) 
	
	
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
		IF (me==0)	PRINT*,'Number of processes:',wsize
		NT=OMP_GET_MAX_THREADS()
		IF (me==0)	PRINT*,'Number of threads:',NT

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
	sig1(1)=1	
	sig1(2)=1*3d0
	sig1(3)=1*10d0
	sig2(1)=.01d0	
	sig2(2)=.01d0/3.d0
	sig2(3)=.01d0/10
	
!	CALL SliceAnomaly(anomaly,bkg,0d0,10000d0,N)
	CALL PrepareRecvs(recvs,anomaly,bkg)
	IF (me==0) PRINT'(A80)','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	IF (me==0) PRINT*, 'Nx=',anomaly%Nx, 'Ny=',anomaly%Ny,'Nz=',anomaly%Nz
	CALL PrepareIntegralEquation(int_eq,anomaly,wcomm,threads_ok)
	real_comm=int_eq%fgmres_comm
	CALL PrepareContinuationOperator(rc_op,anomaly,recvs,wcomm,threads_ok)

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


		CALL CalcRecalculationGreenTensor(rc_op,bkg,anomaly,RC_Threshold)
		CALL CalcFFTofRCTensor(rc_op)

#ifndef not_compile
		DO Ia=1,Na
			WRITE (fnum2,'(I5.5)') Ia
			IF (me==0)	PRINT*, 'Anomaly ', Ia, ' from ', trim(anom_list(Ia))
!#endif			  
			IF (int_eq%real_space) THEN
				CALL PlaneWaveIntegral(EY,bkg,anomaly,int_eq%E_n)
				CALL PlaneWave(EY,bkg,(/recvs(1)%zrecv/),FY)
!				anomaly%siga(:,1:int_eq%Nx/2,:)=sig1(It)
!				anomaly%siga(:,int_eq%Nx/2+1:,:)=sig2(It)
                CALL LoadAnomalySigma(anomaly,real_comm,trim(anom_list(Ia)))

			ENDIF
				IF (me==0) PRINT*, 'FY=', FY
			CALL SolveEquation(int_eq,fgmres_ctl)

			CALL ReCalculation(rc_op,int_eq%Esol,Ea,Ha)
			IF (int_eq%real_space) THEN

				Et(:,EX,:,:)=Ea(:,EX,:,:)+FY(1,1,1,EX)
				Et(:,EY,:,:)=Ea(:,EY,:,:)+FY(1,1,1,EY)
				Et(:,EZ,:,:)=Ea(:,EZ,:,:)+FY(1,1,1,EZ)

				Ht(:,HX,:,:)=Ha(:,HX,:,:)+FY(1,1,1,HX)
				Ht(:,HY,:,:)=Ha(:,HY,:,:)+FY(1,1,1,HY)
				Ht(:,HZ,:,:)=Ha(:,HZ,:,:)+FY(1,1,1,HZ)

				!CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PY_F'//trim(fnum1)//'T_'//trim(fnum2))
				CALL SaveOutputSeparate(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),rc_op%me,rc_op%Ny_offset,'PY_F'//trim(fnum1)//'T_'//trim(fnum2))
			ENDIF

			IF (int_eq%real_space) THEN
				CALL PlaneWaveIntegral(EX,bkg,anomaly,int_eq%E_n)
				CALL PlaneWave(EX,bkg,(/recvs(1)%zrecv/),FX)
			ENDIF

			CALL SolveEquation(int_eq,fgmres_ctl)
			CALL ReCalculation(rc_op,int_eq%Esol,Ea,Ha)

			IF (int_eq%real_space) THEN

				Et(:,EX,:,:)=Ea(:,EX,:,:)+FX(1,1,1,EX)
				Et(:,EY,:,:)=Ea(:,EY,:,:)+FX(1,1,1,EY)
				Et(:,EZ,:,:)=Ea(:,EZ,:,:)+FX(1,1,1,EZ)

				Ht(:,HX,:,:)=Ha(:,HX,:,:)+FX(1,1,1,HX)
				Ht(:,HY,:,:)=Ha(:,HY,:,:)+FX(1,1,1,HY)
				Ht(:,HZ,:,:)=Ha(:,HZ,:,:)+FX(1,1,1,HZ)
				CALL SaveOutputSeparate(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),rc_op%me,rc_op%Ny_offset,'PX_F'//trim(fnum1)//'T_'//trim(fnum2))

!				CALL SaveOutputOneFile(Ea,Et,Ha,Ht,anomaly,recvs,freqs(Ifreq),real_comm,'PX_F'//trim(fnum1)//'T_'//trim(fnum2))
			ENDIF
		ENDDO
#endif
	ENDDO
	 CALL DeleteIE_OP(int_eq)
	 CALL DeleteRC_OP(RC_OP)
	CALL MPI_FINALIZE(IERROR)
END PROGRAM
