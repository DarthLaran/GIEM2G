MODULE Calc_IE_Tensor_Module
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE DATA_TYPES_MODULE
	USE IE_Kernel_Image_Module
	USE INTEGRAL_EQUATION_MODULE

	USE IntegralCodes 
	USE VolumeBesselTransforms
	IMPLICIT NONE
	REAL(REALPARM)::lms(Nfirst:Nlast)
	PUBLIC::CalcIntegralGreenTensor
CONTAINS
	SUBROUTINE CalcIntegralGreenTensor(matrix,bkg,anomaly,Wt_Threshold)
		CLASS(IE_MATRIX),INTENT(INOUT)::matrix
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		REAL(REALPARM),INTENT(IN)::Wt_Threshold
		COMPLEX(REALPARM),POINTER::G_symm(:,:,:,:)
		COMPLEX(REALPARM),POINTER::G_asym(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::pGsend_symm(:,:,:)
		COMPLEX(REALPARM),POINTER::pGrecv_symm(:,:,:)
		COMPLEX(REALPARM),POINTER::pGsend_asym(:,:,:,:)
		COMPLEX(REALPARM),POINTER::pGrecv_asym(:,:,:,:)
		INTEGER::Ix,Iy,ly,Ic,IERROR
		INTEGER::Nx,Ny,Ny_loc,Ny_offset,Nz,Nx2
		REAL(REALPARM)::dx,dy
		INTEGER::Nsend_symm,Nrecv_symm
		INTEGER::Nsend_asym,Nrecv_asym
		INTEGER::recv_proc,send_proc,shift
		INTEGER::xfirst,xlast
		INTEGER::recv_start,send_start
		INTEGER::requests(matrix%Ny_loc,4)
		REAL(8)::time(4),time_all!,Wt_Threshold
		REAL(8)::time_max(4),time_min(4)
		REAL(8)::time1,time2,time3,time4
		CALL MPI_BARRIER(matrix%matrix_comm,IERROR)
		time_all=MPI_WTIME()
		IF (VERBOSE) THEN
			IF (matrix%master) THEN
				PRINT'(A80)','***********************************************************************************'
				PRINT*, '  IE kernel with threshold', WT_Threshold, 'start'
			ENDIF
		ENDIF
		time=0d0
		Nx=matrix%Nx 
		Nx2=Nx/2
		Ny=matrix%Ny 
		Ny_loc=matrix%Ny_loc 
		Ny_offset=matrix%Ny_offset
		Nz=matrix%Nz
		G_symm=>matrix%G_symm
		G_asym=>matrix%G_asym
		G_symm=C_ZERO
		G_asym=C_ZERO
		CALL VBTransformInit(lms)
!		CALL Prepare_Hankel(bkg,anomaly)
		dx=anomaly%dx
		dy=anomaly%dy
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
				&matrix%matrix_comm,requests(Iy+1-Ny_offset,1),IERROR)
			CALL MPI_IRECV(pGrecv_asym,Nrecv_asym,MPI_DOUBLE_COMPLEX,send_proc,Iy,&
				&matrix%matrix_comm,requests(Iy+1-Ny_offset,3),IERROR)
			!$OMP PARALLEL PRIVATE(Ix),DEFAULT(SHARED)&
			!$OMP &FIRSTPRIVATE(time)
				!$OMP DO SCHEDULE(GUIDED)&
				!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
				DO Ix=xfirst,xlast
					CALL CalcIntegralGreenTensor3dElement(Ix,ly,G_symm(:,:,Ix,Iy),&
					&G_asym(:,:,:,Ix,Iy),Wt_Threshold,time,anomaly,bkg,matrix%dz)
				ENDDO
				!$OMP ENDDO
				time1=time(1)
				time2=time(2)
				time3=time(3)
				time4=time(4)
			!$OMP END PARALLEL
			time(1)=time1
			time(2)=time2
			time(3)=time3
			time(4)=time4
			CALL MPI_ISEND(pGsend_symm,Nsend_symm,MPI_DOUBLE_COMPLEX,recv_proc,shift-Iy,&
				&matrix%matrix_comm,requests(Iy+1-Ny_offset,2),IERROR)
			CALL MPI_ISEND(pGsend_asym,Nsend_asym,MPI_DOUBLE_COMPLEX,recv_proc,shift-Iy,&
				&matrix%matrix_comm,requests(Iy+1-Ny_offset,4),IERROR)
		ENDDO
		CALL MPI_WAITALL(Ny_loc*4, requests, MPI_STATUSES_IGNORE, IERROR)  
		IF (Ny_offset>=Ny) THEN
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
					G_symm(:,S_EXX,2*Nx-Ix,Iy)=G_symm(:,S_EXX,Ix,Iy)
					G_symm(:,S_EXY,2*Nx-Ix,Iy)=G_symm(:,S_EXY,Ix,Iy)
					G_symm(:,S_EXY,Ix,Iy)=-G_symm(:,S_EXY,Ix,Iy)
					G_symm(:,S_EYY,2*Nx-Ix,Iy)=G_symm(:,S_EYY,Ix,Iy)
					G_symm(:,S_EZZ,2*Nx-Ix,Iy)=G_symm(:,S_EZZ,Ix,Iy)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
				G_symm(:,S_EXY,0,Iy)=-G_symm(:,S_EXY,0,Iy)
				G_symm(:,:,Nx,Iy)=C_ZERO
				G_asym(:,:,:,Nx,Iy)=C_ZERO
				!$OMP END WORKSHARE
				!$OMP DO SCHEDULE(GUIDED)
				DO Ix=1,Nx-1
					G_asym(:,:,A_EXZ,2*Nx-Ix,Iy)=-G_asym(:,:,A_EXZ,Ix,Iy)

					G_asym(:,:,A_EYZ,Ix,Iy)=-G_asym(:,:,A_EYZ,Ix,Iy)
					G_asym(:,:,A_EYZ,2*Nx-Ix,Iy)=G_asym(:,:,A_EYZ,Ix,Iy)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
					G_asym(:,:,A_EYZ,0,Iy)=-G_asym(:,:,A_EYZ,0,Iy)
				!$OMP END WORKSHARE
				!$OMP END PARALLEL
			ENDDO
		ELSE
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
		ENDIF
		time(3)=time(2)-time(4) 
		CALL MPI_REDUCE(time,time_min,4,MPI_DOUBLE,MPI_MIN,matrix%master_proc,matrix%matrix_comm,IERROR)
		CALL MPI_REDUCE(time,time_max,4,MPI_DOUBLE,MPI_MAX,matrix%master_proc,matrix%matrix_comm,IERROR)
		time_all=MPI_WTIME()-time_all
		IF (matrix%master) THEN
			matrix%counter%tensor_calc(COUNTER_ALL)=time_all	
			matrix%counter%tensor_calc(COUNTER_WT)=time_max(1)	
			matrix%counter%tensor_calc(COUNTER_LIN)=time_max(3) 
			matrix%counter%tensor_calc(COUNTER_SQR)=time_max(4) 
			matrix%counter%tensor_calc(COUNTER_OTHER)=time_all-time_max(2)-time_max(1)	
			IF (VERBOSE) THEN
				PRINT *,'Electrical Green Tensor has been computed.'
				PRINT '(A12, 1ES10.2E2, A3)','Full time:  ',time_all,'s'
				PRINT '(A12, 2ES10.2E2, A3)','Weights:	  ',time_min(1), time_max(1),'s'
				PRINT '(A12, 2ES10.2E2, A3)','Linear part:',time_min(3), time_max(3),'s'
				PRINT '(A12, 2ES10.2E2, A3)','Square part:',time_min(4), time_max(4),'s'
				PRINT '(A12, 1ES10.2E2, A3)','Other:	  ',matrix%counter%tensor_calc(COUNTER_OTHER),'s'
				PRINT '(A80)','***********************************************************************************'
			ENDIF
		ELSE
			matrix%counter%tensor_calc(COUNTER_ALL)=time_all	
			matrix%counter%tensor_calc(COUNTER_WT)=time(1)	
			matrix%counter%tensor_calc(COUNTER_LIN)=time(3) 
			matrix%counter%tensor_calc(COUNTER_SQR)=time(4) 
			matrix%counter%tensor_calc(COUNTER_OTHER)=time_all-time(2)-time(1)	
		ENDIF
	END SUBROUTINE
!---------------------------------------------------------------------------------------------------------------------------------------------!
	SUBROUTINE CalcIntegralGreenTensor3dElement(lx,ly,G_symm,G_asym,Wt_Threshold,calc_time,anomaly,bkg,dz)
		IMPLICIT NONE
		INTEGER, INTENT(IN)::lx,ly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
		COMPLEX(REALPARM),INTENT(OUT)::G_symm(anomaly%Nz*(anomaly%Nz+1)/2,S_EXX:S_EZZ)
		COMPLEX(REALPARM),INTENT(OUT)::G_asym(anomaly%Nz,anomaly%Nz,A_EXZ:A_EYZ)
		REAL(REALPARM),INTENT(IN)::Wt_Threshold
		REAL(REALPARM),INTENT(INOUT)::calc_time(4)
		COMPLEX(REALPARM)::G1(EXX:EZZ,anomaly%Nz,anomaly%Nz)
		REAL(REALPARM2)::r,WT(Nfirst:Nlast,1:6),W(1:6),dx,dy
		REAL(REALPARM2)::WT0(Nfirst:Nlast,1:6)
		REAL(REALPARM2)::WT1(Nfirst:Nlast,1:6)
		REAL(REALPARM)::x,y,lm,W0(1:6),Wm(1:6),time1,time2
		INTEGER::K,Iz,IERROR,I,Ik,Iz0,Ic,Iwt(6)
		time1=MPI_Wtime()
		dx=anomaly%dx
		dy=anomaly%dy
		x=lx*dx+dx/2d0
		y=ly*dy+dy/2d0
		r=SQRT(x*x+y*y)
		CALL VBTransformWeightsAllInt42(x,y,dx,dy,WT0)
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
		G1=C_ZERO		
		time2=MPI_Wtime()
		calc_time(1)=calc_time(1)+time2-time1
		time1=time2
		DO K=Nfirst,Nlast
			Wm=ABS(WT(K,:)/W0(:))
			IF (MAXVAL(Wm)>Wt_Threshold) THEN
				lm=lms(K)/r
				W(IE_DXX)=WT(K,IE_DXX)/lm
				W(IE_DYY)=WT(K,IE_DYY)/lm
				W(IE_DXY)=WT(K,IE_DXY)/lm
				W(IE_DX)=WT(K,IE_DX)*lm
				W(IE_DY)=WT(K,IE_DY)*lm
				W(IE_D0)=WT(K,IE_D0)*lm
				CALL Calc_Double_Integral_U(anomaly,bkg,dz,lm,W,G1,calc_time(3:4))

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
		time2=MPI_WTIME()
		calc_time(2)=calc_time(2)+time2-time1
	END SUBROUTINE
ENDMODULE
