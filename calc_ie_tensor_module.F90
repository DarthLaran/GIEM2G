MODULE Calc_IE_Tensor_Module
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE Timer_Module 
	USE DATA_TYPES_MODULE
	USE IE_Kernel_Image_Module
	USE INTEGRAL_EQUATION_MODULE

	USE IntegralCodes 
	USE VolumeBesselTransforms
	USE LOGGER_MODULE
	USE LOCAL_OMP_FFT_MODULE

	IMPLICIT NONE
	PRIVATE
	REAL(REALPARM)::lms(Nfirst:Nlast)
	REAL(REALPARM),PARAMETER::Wt_Threshold=1d-12
	TYPE(LOCAL_OMP_FFT_DATA),POINTER::LFFT
	PUBLIC::CalcIntegralGreenTensor
	INTERFACE
		SUBROUTINE   CalcTensorAlongXSubroutine(xfirst,xlast,ly,G_symm,G_asym,time,anomaly,bkg,dz)
				IMPORT
				INTEGER,INTENT(IN)::xfirst,xlast,ly
				COMPLEX(REALPARM),POINTER,INTENT(IN)::G_symm(:,:,:)
				COMPLEX(REALPARM),POINTER,INTENT(IN)::G_asym(:,:,:,:)
				REAL(DOUBLEPARM),INTENT(INOUT)::time(4)
				TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
				TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
				REAL(REALPARM),POINTER,INTENT(IN)::dz(:)
				INTEGER::Ix
		ENDSUBROUTINE
	ENDINTERFACE
CONTAINS
	SUBROUTINE CalcIntegralGreenTensor(ie_op,bkg,anomaly)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		COMPLEX(REALPARM),POINTER::G_symm(:,:,:,:)
		COMPLEX(REALPARM),POINTER::G_asym(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::pGsend_symm(:,:,:)
		COMPLEX(REALPARM),POINTER::pGrecv_symm(:,:,:)
		COMPLEX(REALPARM),POINTER::pGsend_asym(:,:,:,:)
		COMPLEX(REALPARM),POINTER::pGrecv_asym(:,:,:,:)
		PROCEDURE(CalcTensorAlongXSubroutine),POINTER::CalcTensorAlongX
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER::Ix,Iy,ly,Ic
		INTEGER::Nx,Ny,Ny_loc,Ny_offset,Nz,Nx2
		REAL(REALPARM)::dx,dy,rmin
		INTEGER::Nsend_symm,Nrecv_symm
		INTEGER::Nsend_asym,Nrecv_asym
		INTEGER::recv_proc,send_proc,shift
		INTEGER::xfirst,xlast
		INTEGER::recv_start,send_start
!------------------------- workaround for bug in some realizations of MPI_WAITALL
		INTEGER(MPI_CTL_KIND),TARGET::all_requests(4*ie_op%Ny_loc)
		INTEGER(MPI_CTL_KIND),POINTER::requests(:,:)
!---------------------------------------------------------------------------
		REAL(8)::time(4),time_all!,Wt_Threshold
		REAL(8)::time_max(4),time_min(4)
		REAL(8)::time1,time2,time3,time4
!		CALL MPI_BARRIER(ie_op%ie_op_comm,IERROR)
		time_all=GetTime()
		requests(1:ie_op%Ny_loc,1:4)=>all_requests
		CALL PRINT_BORDER
		CALL LOGGER('Start IE kernel calculation')
		time=0d0
		Nx=ie_op%Nx 
		Nx2=Nx/2
		Ny=ie_op%Ny 
		Ny_loc=ie_op%Ny_loc 
		Ny_offset=ie_op%Ny_offset
		Nz=ie_op%Nz

		G_symm=>ie_op%G_symm
		G_symm=C_ZERO
		G_asym=>ie_op%G_asym
		G_asym=C_ZERO

		IF (ie_op%matrix_kind==GENERAL_MATRIX) THEN
			CalcTensorAlongX=>CalcGeneralTensorAlongX
			LFFT=>NULL()
		ELSEIF (ie_op%matrix_kind==UNIFORM_MATRIX) THEN
			CalcTensorAlongX=>CalcUniformTensorAlongX
			LFFT=>ie_op%LFFT
		ELSE
		ENDIF	



		dx=anomaly%dx
		dy=anomaly%dy
		CALL	PreparePreconditioner(ie_op,bkg,anomaly)
		rmin=SQRT(dx*dx+dy*dy)/2
		CALL VBTransformInit(rmin,lms)

		IF (Ny_offset<Ny) THEN
			Nsend_asym=Nz*Nz*2*Nx2
			Nrecv_asym=Nz*Nz*2*(Nx-Nx2)
			Nsend_symm=Nz*(Nz+1)*2*Nx2
			Nrecv_symm=Nz*(Nz+1)*2*(Nx-Nx2)
			xfirst=0
			xlast=Nx2-1
			recv_start=Nx2
			send_start=0
		ELSE
			Nsend_asym=Nz*Nz*2*(Nx-Nx2)
			Nrecv_asym=Nz*Nz*2*(Nx2)
			Nsend_symm=Nz*(Nz+1)*2*(Nx-Nx2)
			Nrecv_symm=Nz*(Nz+1)*2*(Nx2)
			xfirst=Nx2
			xlast=Nx-1
			recv_start=0
			send_start=Nx2
		ENDIF

		DO Iy=Ny_offset,Ny_offset+Ny_loc-1
			pGsend_symm=>G_symm(:,:,send_start:,Iy)
			pGrecv_symm=>G_symm(:,:,recv_start:,Iy)
			pGsend_asym=>G_asym(:,:,:,send_start:,Iy)
			pGrecv_asym=>G_asym(:,:,:,recv_start:,Iy)
			IF ((Iy==Ny).OR.(Iy==0)) THEN
				shift=Ny
			ELSE
				shift=2*Ny
			ENDIF
			recv_proc=(shift-Iy)/Ny_loc
			send_proc=(shift-Iy)/Ny_loc
			IF (Ny_offset<Ny) THEN
				ly=Iy
			ELSE
				ly=shift-Iy
			ENDIF
			CALL MPI_IRECV(pGrecv_symm,Nrecv_symm,MPI_DOUBLE_COMPLEX,send_proc,Iy,&
				&ie_op%ie_comm,requests(Iy+1-Ny_offset,1),IERROR)
			CALL MPI_IRECV(pGrecv_asym,Nrecv_asym,MPI_DOUBLE_COMPLEX,send_proc,Iy,&
				&ie_op%ie_comm,requests(Iy+1-Ny_offset,3),IERROR)

			CALL	CalcTensorAlongX(xfirst,xlast,ly,G_symm(:,:,:,Iy),G_asym(:,:,:,:,Iy),&
					&time,anomaly,bkg,ie_op%dz)

			CALL MPI_ISEND(pGsend_symm,Nsend_symm,MPI_DOUBLE_COMPLEX,recv_proc,shift-Iy,&
				&ie_op%ie_comm,requests(Iy+1-Ny_offset,2),IERROR)
			CALL MPI_ISEND(pGsend_asym,Nsend_asym,MPI_DOUBLE_COMPLEX,recv_proc,shift-Iy,&
				&ie_op%ie_comm,requests(Iy+1-Ny_offset,4),IERROR)
		ENDDO
		
		CALL MPI_WAITALL(Ny_loc*4, all_requests,  MPI_STATUSES_IGNORE, IERROR)

		IF (Ny_offset>=Ny) THEN
			CALL MirroringAdditionalSpace(Ny_offset,Ny_loc,Ny,Nx, G_symm,G_asym)
!			G_symm(:,S_EXY,:,:)=-G_symm(:,S_EXY,:,:)
!			G_asym(:,:,A_EYZ,:,:)=-G_asym(:,:,A_EYZ,:,:)
		ELSE
!			CALL MirroringRealSpace(Ny_offset,Ny_loc,Ny,Nx, G_symm,G_asym)
		ENDIF

#ifdef internal_timer
		time(3)=time(2)-time(4) 
		CALL MPI_REDUCE(time,time_min,4,MPI_DOUBLE,MPI_MIN,ie_op%master_proc,ie_op%ie_comm,IERROR)
		CALL MPI_REDUCE(time,time_max,4,MPI_DOUBLE,MPI_MAX,ie_op%master_proc,ie_op%ie_comm,IERROR)
#endif
		time_all=GetTime()-time_all

#ifdef internal_timer
		IF (ie_op%master) THEN
			ie_op%counter%tensor_calc(COUNTER_ALL)=time_all	
			ie_op%counter%tensor_calc(COUNTER_WT)=time_max(1)	
			ie_op%counter%tensor_calc(COUNTER_LIN)=time_max(3) 
			ie_op%counter%tensor_calc(COUNTER_SQR)=time_max(4) 
			ie_op%counter%tensor_calc(COUNTER_OTHER)=time_all-time_max(2)-time_max(1)	
			IF (VERBOSE) THEN
				  CALL LOGGER('IE kernel has been computed.')
				 CALL PRINT_CALC_TIME(CALC_IE_FULL_TIME,time_all)
				 CALL PRINT_CALC_TIME(CALC_IE_WEIGHTS, time_max(1))
				 CALL PRINT_CALC_TIME(CALC_IE_LINEAR, time_max(3))
				 CALL PRINT_CALC_TIME(CALC_IE_SQ,    time_max(4))
				 CALL PRINT_CALC_TIME(CALC_IE_OTHER,ie_op%counter%tensor_calc(COUNTER_OTHER))
				CALL PRINT_BORDER
			ENDIF
		ELSE
			ie_op%counter%tensor_calc(COUNTER_ALL)=time_all	
			ie_op%counter%tensor_calc(COUNTER_WT)=time(1)	
			ie_op%counter%tensor_calc(COUNTER_LIN)=time(3) 
			ie_op%counter%tensor_calc(COUNTER_SQR)=time(4) 
			ie_op%counter%tensor_calc(COUNTER_OTHER)=time_all-time(2)-time(1)	
		ENDIF
		ie_op%counter%tensor_calc(COUNTER_ALL)=time_all	
#else
		CALL LOGGER('IE kernel has been computed.')
		CALL PRINT_CALC_TIME(CALC_IE_FULL_TIME,time_all)
		CALL PRINT_BORDER
#endif
	END SUBROUTINE
!---------------------------------------------------------------------------------------------------------------------------------------------!
	SUBROUTINE   CalcGeneralTensorAlongX(xfirst,xlast,ly,G_symm,G_asym,time,anomaly,bkg,dz)
			INTEGER,INTENT(IN)::xfirst,xlast,ly
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_symm(:,:,:)
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_asym(:,:,:,:)
			REAL(DOUBLEPARM),INTENT(INOUT)::time(4)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			REAL(REALPARM),POINTER,INTENT(IN)::dz(:)
			COMPLEX(REALPARM),POINTER::G_symm1(:,:,:)
			COMPLEX(REALPARM),POINTER::G_asym1(:,:,:,:)
			INTEGER::Ix
			G_symm1(1:,S_EXX:,0:)=>G_symm
			G_asym1(1:,1:,A_EXZ:,0:)=>G_asym
			!$OMP PARALLEL PRIVATE(Ix),DEFAULT(SHARED)&
			!$OMP &FIRSTPRIVATE(time)
#ifdef internal_timer
			!$OMP DO SCHEDULE(GUIDED)&
			!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
#else
			!$OMP DO SCHEDULE(GUIDED)
#endif
				DO Ix=xfirst,xlast
					CALL CalcIntegralGreenTensor3dElementGeneral(Ix,ly,G_symm1(:,:,Ix),&
					&G_asym1(:,:,:,Ix),time,anomaly,bkg,dz)
				ENDDO
				!$OMP ENDDO
#ifdef internal_timer
				time1=time(1)
				time2=time(2)
				time3=time(3)
				time4=time(4)
			!$OMP END PARALLEL
			time(1)=time1
			time(2)=time2
			time(3)=time3
			time(4)=time4
#else
			!$OMP END PARALLEL
#endif
	ENDSUBROUTINE

	SUBROUTINE CalcIntegralGreenTensor3dElementGeneral(lx,ly,G_symm,G_asym,calc_time,anomaly,bkg,dz)
		INTEGER, INTENT(IN)::lx,ly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
		COMPLEX(REALPARM),INTENT(INOUT)::G_symm(anomaly%Nz*(anomaly%Nz+1)/2,S_EXX:S_EZZ)
		COMPLEX(REALPARM),INTENT(INOUT)::G_asym(anomaly%Nz,anomaly%Nz,A_EXZ:A_EYZ)
		REAL(REALPARM),INTENT(INOUT)::calc_time(4)
		COMPLEX(REALPARM)::G1(EXX:EZZ,anomaly%Nz,anomaly%Nz)
		REAL(REALPARM2)::W0(1:6),W(1:6),dx,dy
		REAL(REALPARM2)::WT(Nfirst:Nlast,1:6)
		REAL(REALPARM)::lm,Wm(1:6),time1,time2
		INTEGER::K,Iz,Iz0,Ng
		time1=GetTime()
		dx=anomaly%dx
		dy=anomaly%dy
		CALL COMPUTE_WEIGHTS(lx,ly,dx,dy,WT,W0)

#ifdef internal_timer
		time2=GetTime()
		calc_time(1)=calc_time(1)+time2-time1
		time1=time2
#endif
		G1=C_ZERO		
		DO K=Nfirst,Nlast
			Wm=ABS(WT(K,:)/W0(:))
			IF (MAXVAL(Wm)>Wt_Threshold) THEN
				lm=lms(K)!r
				W(IE_DXX)=WT(K,IE_DXX)/lm
				W(IE_DYY)=WT(K,IE_DYY)/lm
				W(IE_DXY)=WT(K,IE_DXY)/lm
				W(IE_DX)=WT(K,IE_DX)*lm
				W(IE_DY)=WT(K,IE_DY)*lm
				W(IE_D0)=WT(K,IE_D0)*lm
				CALL Calc_Double_Integral_General(anomaly,bkg,dz,lm,W,G1,calc_time(3:4))
			ENDIF
		ENDDO
		G_asym(:,:,A_EXZ)=(G1(EXZ,:,:)/PI/4.0_REALPARM)
		G_asym(:,:,A_EYZ)=(G1(EYZ,:,:)/PI/4.0_REALPARM)
				
		DO Iz=1,anomaly%Nz
			DO Iz0=1,Iz
				G_symm(Iz0+Iz*(Iz-1)/2,S_EXX)=G1(EXX,Iz0,Iz)/PI/4.0_REALPARM
				G_symm(Iz0+Iz*(Iz-1)/2,S_EXY)=G1(EXY,Iz0,Iz)/PI/4.0_REALPARM
				G_symm(Iz0+Iz*(Iz-1)/2,S_EYY)=G1(EYY,Iz0,Iz)/PI/4.0_REALPARM
				G_symm(Iz0+Iz*(Iz-1)/2,S_EZZ)=G1(EZZ,Iz0,Iz)/PI/4.0_REALPARM
			ENDDO
		ENDDO
#ifdef internal_timer
		time2=GetTime()
		calc_time(2)=calc_time(2)+time2-time1
#endif
	END SUBROUTINE

	SUBROUTINE   CalcUniformTensorAlongX(xfirst,xlast,ly,G_symm,G_asym,time,anomaly,bkg,dz)
			INTEGER,INTENT(IN)::xfirst,xlast,ly
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_symm(:,:,:)
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_asym(:,:,:,:)
			REAL(DOUBLEPARM),INTENT(INOUT)::time(4)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			REAL(REALPARM),POINTER,INTENT(IN)::dz(:)
			COMPLEX(REALPARM),POINTER::G_symm1(:,:,:)
			COMPLEX(REALPARM),POINTER::G_symm2(:,:,:,:)
			COMPLEX(REALPARM),POINTER::G_asym1(:,:,:,:)
			TYPE(C_PTR)::ptr_c
			COMPLEX(REALPARM),POINTER::ptr(:)
			INTEGER::Ix,N,TN
			G_symm1(1:,S_EXX:,0:)=>G_symm
			G_asym1(1:,1:,A_EXZ:,0:)=>G_asym
			N=anomaly%Nz*4*2*anomaly%Nx
			ptr_c=C_LOC(G_symm1(1,S_EXX,0))
			CALL C_F_POINTER(ptr_c,ptr,(/N/))
			G_symm2(1:2*anomaly%Nz,1:2,S_EXX:S_EZZ,0:2*anomaly%Nx-1)=>ptr
			!$OMP PARALLEL PRIVATE(Ix,TN),DEFAULT(SHARED)&
			!$OMP &FIRSTPRIVATE(time)
			TN=OMP_GET_THREAD_NUM()
#ifdef internal_timer
			!$OMP DO SCHEDULE(GUIDED)&
			!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
#else
			!$OMP DO SCHEDULE(GUIDED)
#endif
				DO Ix=xfirst,xlast
					CALL CalcIntegralGreenTensor3dElementUniform(Ix,ly,G_symm2(:,:,:,Ix),&
					&G_asym1(:,:,:,Ix),time,anomaly,bkg,dz,TN)
				ENDDO
				!$OMP ENDDO
#ifdef internal_timer
				time1=time(1)
				time2=time(2)
				time3=time(3)
				time4=time(4)
			!$OMP END PARALLEL
			time(1)=time1
			time(2)=time2
			time(3)=time3
			time(4)=time4
#else
			!$OMP END PARALLEL
#endif
	ENDSUBROUTINE
	SUBROUTINE CalcIntegralGreenTensor3dElementUniform(lx,ly,G_symm,G_asym,calc_time,anomaly,bkg,dz,TN)
		INTEGER, INTENT(IN)::lx,ly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
		COMPLEX(REALPARM),INTENT(OUT)::G_symm(2*anomaly%Nz,2,S_EXX:S_EZZ)
		COMPLEX(REALPARM),INTENT(OUT)::G_asym(2*anomaly%Nz,2,A_EXZ:A_EYZ)
		REAL(REALPARM),INTENT(INOUT)::calc_time(4)
		INTEGER,INTENT(IN)::TN
		COMPLEX(REALPARM)::G1(EXX:EZZ,anomaly%Nz,4)
		REAL(REALPARM2)::W0(1:6),W(1:6),dx,dy
		REAL(REALPARM2)::WT(Nfirst:Nlast,1:6)
		REAL(REALPARM)::lm,Wm(1:6),time1,time2
		INTEGER::K,Iz,Iz0,Ng
		time1=GetTime()
		dx=anomaly%dx
		dy=anomaly%dy
		CALL COMPUTE_WEIGHTS(lx,ly,dx,dy,WT,W0)


#ifdef internal_timer
		time2=GetTime()
		calc_time(1)=calc_time(1)+time2-time1
		time1=time2
#endif
		G1=C_ZERO		
		DO K=Nfirst,Nlast
			Wm=ABS(WT(K,:)/W0(:))
			IF (MAXVAL(Wm)>Wt_Threshold) THEN
				lm=lms(K)!r
				W(IE_DXX)=WT(K,IE_DXX)/lm
				W(IE_DYY)=WT(K,IE_DYY)/lm
				W(IE_DXY)=WT(K,IE_DXY)/lm
				W(IE_DX)=WT(K,IE_DX)*lm
				W(IE_DY)=WT(K,IE_DY)*lm
				W(IE_D0)=WT(K,IE_D0)*lm
				CALL Calc_Double_Integral_Uniform(anomaly,bkg,dz,lm,W,G1,calc_time(3:4))
			ENDIF
		ENDDO
		CALL SEPARATE_MATRIX(G1,anomaly%Nz)
                CALL EXTRACT_MATRIX(G1,G_symm,G_asym,anomaly%Nz,TN)
#ifdef internal_timer
		time2=GetTime()
		calc_time(2)=calc_time(2)+time2-time1
#endif
	END SUBROUTINE

	
	SUBROUTINE SEPARATE_MATRIX(G,N)
		COMPLEX(REALPARM),INTENT(INOUT)::G(N,4,EXX:EZZ)
		INTEGER,INTENT(IN)::N
		COMPLEX(REALPARM)::T1(2*N),T2(2*N)
		COMPLEX(REALPARM)::a,b
		INTEGER::I,Iz
		DO I=EXX,EZZ
			a=G(1,1,I)
			b=G(1,2,I) 

			DO Iz=1,N-1
				T1(2*Iz-1)=G(Iz,1,I)-a
				T1(2*Iz)=G(Iz,2,I)-b
			ENDDO

			T1(2*N-1)=G(N,1,Iz)-a

			DO Iz=1,N
				T2(Iz)=T1(N-Iz+1)
			ENDDO  

			T2(N+1)=0

			DO Iz=N+2,2*N
				T2(Iz)=T1(3*N+1-Iz)
			ENDDO  

			DO Iz=1,N
				T1(Iz)=G(Iz,3,I)-T2(N-Iz+1)
			ENDDO

			T1(N+1)=C_ZERO
			T1(N+2)=G(1,4,I)-T2(1)

			DO Iz=2,N-1
			    T1(Iz+N+1)=G(Iz,4,I)-T2(2*N-I+2);
			ENDDO

			G(:,1,I)=T1(1:N)
			G(:,2,I)=T1(N+1:2*N)

			G(:,3,I)=T2(1:N)
			G(:,4,I)=T2(N+1:2*N)
		ENDDO
	ENDSUBROUTINE

	SUBROUTINE EXTRACT_MATRIX(G_in,G_symm,G_asym,Nz,TN)
		COMPLEX(REALPARM),INTENT(IN)::G_in(Nz,4,EXX:EZZ)
		COMPLEX(REALPARM),INTENT(OUT)::G_symm(2*Nz,2,S_EXX:S_EZZ)
		COMPLEX(REALPARM),INTENT(OUT)::G_asym(2*Nz,2,A_EXZ:A_EYZ)
		INTEGER,INTENT(IN)::Nz,TN
		COMPLEX(REALPARM),POINTER::data_in(:,:)
		COMPLEX(REALPARM),POINTER::data_out(:,:)

		data_in(1:2*Nz,1:4)=>LFFT%data_in(:,TN)
		data_out(1:2*Nz,1:4)=>LFFT%data_out(:,TN)

                CALL ZCOPY(8*Nz,G_in(:,:,EXX:EXY),ONE,data_in,ONE)
		CALL CALCULATE_FORWARD_AT_THREAD(LFFT,TN)
                G_symm(:,:,S_EXX)=data_out(:,1:2)/PI/4.0_REALPARM
                G_symm(:,:,S_EXY)=data_out(:,3:4)/PI/4.0_REALPARM


                CALL ZCOPY(8*Nz,G_in(:,:,EXZ:EYY),ONE,data_in,ONE)
		CALL CALCULATE_FORWARD_AT_THREAD(LFFT,TN)
                G_asym(:,:,A_EZX)=data_out(:,1:2)/PI/4.0_REALPARM
                G_symm(:,:,S_EXY)=data_out(:,3:4)/PI/4.0_REALPARM

                CALL ZCOPY(8*Nz,G_in(:,:,EYZ:EZZ),ONE,data_in,ONE)
		CALL CALCULATE_FORWARD_AT_THREAD(LFFT,TN)

                G_asym(:,:,A_EYZ)=data_out(:,1:2)/PI/4.0_REALPARM
                G_symm(:,:,S_EZZ)=data_out(:,3:4)/PI/4.0_REALPARM
	ENDSUBROUTINE

	SUBROUTINE PreparePreconditioner(ie_op,bkg,anomaly)
		TYPE(IntegralEquationOperator),INTENT(INOUT)::ie_op
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		INTEGER::Iz
		IF (ie_op%real_space) THEN
			DO Iz=1,ie_op%Nz
					ie_op%csigb(Iz)=bkg%csigma(anomaly%Lnumber(Iz))
					ie_op%sqsigb(Iz)=SQRT(REAL(bkg%csigma(anomaly%Lnumber(Iz))))
			ENDDO
		ENDIF
	END SUBROUTINE

	SUBROUTINE COMPUTE_WEIGHTS(lx,ly,dx,dy,WT,W0)
		INTEGER,INTENT(IN)::lx,ly
		REAL(REALPARM),INTENT(IN)::dx,dy
		REAL(REALPARM),INTENT(INOUT)::WT(Nfirst:Nlast,1:6)
		REAL(REALPARM),INTENT(INOUT)::W0(1:6)
		REAL(REALPARM2)::WT0(Nfirst:Nlast,1:6)
		REAL(REALPARM)::x,y,lm,r
		INTEGER::Iwt(6)
		x=lx*dx+dx/2d0
		y=ly*dy+dy/2d0
		r=SQRT(x*x+y*y)
		CALL VBTransformWeightsAllDoubleInt(x,y,dx,dy,WT0)
		WT(:,IE_D0)=WT0(:,1)
		WT(:,IE_DXX)=WT0(:,2)
		WT(:,IE_DXY)=WT0(:,3)
		WT(:,IE_DX)=WT0(:,4)
		WT(:,IE_DYY)=WT0(:,5)
		WT(:,IE_DY)=WT0(:,6)

		W0(IE_DXX)=maxval(ABS(WT(:,IE_DXX)))
		IF ((lx/=0).AND.((ly/=0))) THEN
			W0(IE_DXY:IE_DY)=maxval(ABS(WT(:,IE_DXY:IE_DY)));
		ELSEIF(lx==0)THEN
			WT(:,IE_DXY)=0d0
			W0(IE_DXY)=-W0(IE_DXX)
			W0(IE_DYY)=maxval(ABS(WT(:,IE_DYY)))
			Iwt(IE_DYY)=MAXLOC(ABS(WT(:,IE_DYY)),1)
			WT(:,IE_DX)=0d0
			W0(IE_DX)=-W0(IE_DXX)
			IF (ly/=0) THEN
				W0(IE_DY)=maxval(ABS(WT(:,IE_DY)))
			ELSE
				WT(:,IE_DY)=0d0
				W0(IE_DY)=-W0(IE_DXX)
			ENDIF
		ELSE
			WT(:,IE_DXY)=0d0
			W0(IE_DXY)=-W0(IE_DXX)
			W0(IE_DYY)=maxval(ABS(WT(:,IE_DYY)))
			W0(IE_DX)=maxval(ABS(WT(:,IE_DX)))
			WT(:,IE_DY)=0d0
			W0(IE_DY)=-W0(IE_DXX)
		ENDIF
		W0(IE_D0)=maxval(ABS(WT(:,IE_D0)))
	ENDSUBROUTINE
	SUBROUTINE MirroringAdditionalSpace(Ny_offset,Ny_loc,Ny,Nx, G_symm,G_asym)
			INTEGER,INTENT(IN)::Nx,Ny,Ny_offset,Ny_loc
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_symm(:,:,:,:)
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_asym(:,:,:,:,:)
			INTEGER::Ix,Iy
			DO Iy=Ny_offset,Ny_offset+Ny_loc-1
				IF (Iy==Ny) THEN
					!$OMP PARALLEL DEFAULT(SHARED)
					!$OMP WORKSHARE
						G_symm(:,:,:,Iy)=C_ZERO
						G_asym(:,:,:,:,Iy)=C_ZERO
					!$OMP END WORKSHARE
					!$OMP END PARALLEL
					CYCLE
				ENDIF
				!$OMP PARALLEL PRIVATE(Ix), DEFAULT(SHARED)
				!$OMP DO SCHEDULE(GUIDED)
				DO Ix=1,Nx-1
!					G_symm(:,S_EXX,2*Nx-Ix,Iy)=G_symm(:,S_EXX,Ix,Iy)
!					G_symm(:,S_EXY,2*Nx-Ix,Iy)=G_symm(:,S_EXY,Ix,Iy)
					G_symm(:,S_EXY,Ix,Iy)=-G_symm(:,S_EXY,Ix,Iy)
!					G_symm(:,S_EYY,2*Nx-Ix,Iy)=G_symm(:,S_EYY,Ix,Iy)
!					G_symm(:,S_EZZ,2*Nx-Ix,Iy)=G_symm(:,S_EZZ,Ix,Iy)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
				G_symm(:,S_EXY,0,Iy)=-G_symm(:,S_EXY,0,Iy)
!				G_symm(:,:,Nx,Iy)=C_ZERO
!				G_asym(:,:,:,Nx,Iy)=C_ZERO
				!$OMP END WORKSHARE
				!$OMP DO SCHEDULE(GUIDED)
				DO Ix=1,Nx-1
!					G_asym(:,:,A_EXZ,2*Nx-Ix,Iy)=-G_asym(:,:,A_EXZ,Ix,Iy)

					G_asym(:,:,A_EYZ,Ix,Iy)=-G_asym(:,:,A_EYZ,Ix,Iy)
!					G_asym(:,:,A_EYZ,2*Nx-Ix,Iy)=G_asym(:,:,A_EYZ,Ix,Iy)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
					G_asym(:,:,A_EYZ,0,Iy)=-G_asym(:,:,A_EYZ,0,Iy)
				!$OMP END WORKSHARE
				!$OMP END PARALLEL
			ENDDO
	ENDSUBROUTINE
	SUBROUTINE MirroringRealSpace(Ny_offset,Ny_loc,Ny,Nx, G_symm,G_asym)
			INTEGER,INTENT(IN)::Nx,Ny,Ny_offset,Ny_loc
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_symm(:,:,:,:)
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_asym(:,:,:,:,:)
			INTEGER::Ix,Iy
			DO Iy=Ny_offset,Ny_offset+Ny_loc-1
				!$OMP PARALLEL PRIVATE(Ix), DEFAULT(SHARED)
				!$OMP DO SCHEDULE(GUIDED)
				DO Ix=1,Nx-1
					G_symm(:,S_EXX,2*Nx-Ix,Iy)=G_symm(:,S_EXX,Ix,Iy)
					G_symm(:,S_EXY,2*Nx-Ix,Iy)=-G_symm(:,S_EXY,Ix,Iy)
					G_symm(:,S_EYY,2*Nx-Ix,Iy)=G_symm(:,S_EYY,Ix,Iy)
					G_symm(:,S_EZZ,2*Nx-Ix,Iy)=G_symm(:,S_EZZ,Ix,Iy)
				ENDDO
				!$OMP ENDDO
				!$OMP DO SCHEDULE(GUIDED)
				DO Ix=1,Nx-1
					G_asym(:,:,A_EXZ,2*Nx-Ix,Iy)=-G_asym(:,:,A_EXZ,Ix,Iy)
					G_asym(:,:,A_EYZ,2*Nx-Ix,Iy)=G_asym(:,:,A_EYZ,Ix,Iy)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
				G_symm(:,:,Nx,Iy)=C_ZERO
				G_asym(:,:,:,Nx,Iy)=C_ZERO
				!$OMP END WORKSHARE
				!$OMP END PARALLEL
			ENDDO
	ENDSUBROUTINE
ENDMODULE
