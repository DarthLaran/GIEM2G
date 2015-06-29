MODULE Calc_RC_Tensor_Module
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE DATA_TYPES_MODULE
	USE RC_Kernel_Image_Module
	USE CONTINUATION_FUNCTION_MODULE

	USE IntegralCodes 
	USE VolumeBesselTransforms
	IMPLICIT NONE
	REAL(REALPARM)::lms(Nfirst:Nlast)
	PUBLIC::CalcRecalculationGreenTensor
CONTAINS
	SUBROUTINE CalcRecalculationGreenTensor(matrix,bkg,anomaly,Wt_Threshold)
		TYPE(RC_OPERATOR),INTENT(INOUT)::matrix
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		REAL(REALPARM),INTENT(IN)::Wt_Threshold
		COMPLEX(REALPARM),POINTER::G_E(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::G_H(:,:,:,:,:)
		INTEGER::Ix,Iy,ly,Ic,IERROR
		INTEGER::Nx,Ny,Ny_loc,Ny_offset,Nz,Nx2
		REAL(REALPARM)::dx,dy
		INTEGER::xfirst,xlast
		REAL(8)::time(4),time_all
		REAL(8)::time_max(4),time_min(4)
		REAL(8)::time1,time2,time3,time4
		CALL MPI_BARRIER(matrix%matrix_comm,IERROR)
		time_all=MPI_WTIME()
		IF (VERBOSE) THEN
			IF (matrix%master) THEN
				PRINT'(A80)','***********************************************************************************'
				PRINT*, '  RC kernel with threshold', WT_Threshold, 'start'
			ENDIF
		ENDIF
		time=0d0
		Nx=matrix%Nx 
		Nx2=Nx/2
		Ny=matrix%Ny 
		Ny_loc=matrix%Ny_loc 
		Ny_offset=matrix%Ny_offset
		Nz=matrix%Nz
		G_E=>matrix%G_E
		G_H=>matrix%G_H
		G_E=C_ZERO
		G_H=C_ZERO
		CALL VBTransformInit(lms)
		dx=anomaly%dx
		dy=anomaly%dy
		DO Iy=Ny_offset,Ny_offset+Ny_loc-1
			IF (Iy==Ny) THEN

				!$OMP PARALLEL DEFAULT(SHARED)
				!$OMP WORKSHARE
					G_E(:,:,:,:,Iy)=C_ZERO
					G_H(:,:,:,:,Iy)=C_ZERO
				!$OMP END WORKSHARE
				!$OMP END PARALLEL
				CYCLE
			ENDIF
			IF (Ny_offset<Ny) THEN
				ly=Iy
			ELSE
				ly=Iy-2*Ny
			ENDIF

			!$OMP PARALLEL PRIVATE(Ix),DEFAULT(SHARED)&
			!$OMP &FIRSTPRIVATE(time)
				!$OMP DO SCHEDULE(GUIDED)&
				!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
				DO Ix=0,Nx-1
					CALL CalcRecalculationGreenTensor3dElement(Ix,ly,G_H(:,:,:,Ix,Iy),&
					&G_E(:,:,:,Ix,Iy),Wt_Threshold,time,anomaly,bkg,matrix%recvs,matrix%Nr)
				ENDDO
				!$OMP ENDDO
				!$OMP DO SCHEDULE(GUIDED)&
				!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
				DO Ix=Nx+1,2*Nx-1
					CALL CalcRecalculationGreenTensor3dElement(Ix-2*Nx,ly,G_H(:,:,:,Ix,Iy),&
					&G_E(:,:,:,Ix,Iy),Wt_Threshold,time,anomaly,bkg,matrix%recvs,matrix%Nr)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
					G_E(:,:,:,Nx,Iy)=C_ZERO
					G_H(:,:,:,Nx,Iy)=C_ZERO
				!$OMP END WORKSHARE
				time1=time(1)
				time2=time(2)
				time3=time(3)
				time4=time(4)
			!$OMP END PARALLEL
			time(1)=time1
			time(2)=time2
			time(3)=time3
			time(4)=time4

		ENDDO
		CALL MPI_REDUCE(time,time_min,4,MPI_DOUBLE,MPI_MIN,matrix%master_proc,matrix%matrix_comm,IERROR)
		CALL MPI_REDUCE(time,time_max,4,MPI_DOUBLE,MPI_MAX,matrix%master_proc,matrix%matrix_comm,IERROR)
		time_all=MPI_WTIME()-time_all
		IF (matrix%master) THEN
			matrix%counter%tensor_calc(COUNTER_ALL)=time_all	
			matrix%counter%tensor_calc(COUNTER_WT)=time_max(1)	
			matrix%counter%tensor_calc(COUNTER_LIN)=time_max(3) 
			matrix%counter%tensor_calc(COUNTER_SQR)=time_max(4) 
			matrix%counter%tensor_calc(COUNTER_OTHER)=time_all-time_max(2)	
			IF (VERBOSE) THEN
				PRINT *,'Recalculation Green Tensor has been computed.'
				PRINT '(A14, 1ES10.2E2, A3)','Full time:     ',time_all,'s'
				PRINT '(A14, 2ES10.2E2, A3)','Weights:	     ',time_min(1), time_max(1),'s'
				PRINT '(A14, 2ES10.2E2, A3)','Exps and so on:',time_min(3), time_max(3),'s'
				PRINT '(A14, 2ES10.2E2, A3)','Addition:      ',time_min(4), time_max(4),'s'
				PRINT '(A14, 1ES10.2E2, A3)','Other:	     ',matrix%counter%tensor_calc(COUNTER_OTHER),'s'
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
	SUBROUTINE CalcRecalculationGreenTensor3dElement(lx,ly,G_H,G_E,Wt_Threshold,calc_time,anomaly,bkg,recvs,Nr)
		IMPLICIT NONE
		INTEGER, INTENT(IN)::lx,ly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		TYPE (RECEIVER_TYPE),INTENT(IN)::recvs(Nr)
		INTEGER,INTENT(IN)::Nr
		COMPLEX(REALPARM),INTENT(OUT)::G_H(anomaly%Nz,Nr,RHXX:RHZY)
		COMPLEX(REALPARM),INTENT(OUT)::G_E(anomaly%Nz,Nr,REXX:REZZ)
		REAL(REALPARM),INTENT(IN)::Wt_Threshold
		REAL(REALPARM),INTENT(INOUT)::calc_time(4)
		COMPLEX(REALPARM)::G_E1(REXX:REZZ,anomaly%Nz)
		COMPLEX(REALPARM)::G_H1(RHXX:RHZY,anomaly%Nz)
		REAL(REALPARM2)::r,WT(Nfirst:Nlast,1:6),W(1:6),dx,dy,s
		REAL(REALPARM)::x,y,lm,W0(1:6),Wm(1:6),time1,time2,time3
		INTEGER::K,Iz,IERROR,I,Ik,Iz0,Ic,Iwt(6),Irecv
		REAL(REALPARM2)::WT0(Nfirst:Nlast,1:6)
		REAL(REALPARM2)::WT1(Nfirst:Nlast,1:6)
		time3=MPI_Wtime()
		dx=anomaly%dx
		dy=anomaly%dy
		s=1d0/4.0_REALPARM/PI!*dx*dy
		DO Irecv=1,Nr
			time1=MPI_Wtime()
			x=lx*dx-dx/2d0+recvs(Irecv)%x_shift
			y=ly*dy-dy/2d0+recvs(Irecv)%y_shift
			r=SQRT(x*x+y*y)
			CALL VBTransformWeightsAllInt22(x,y,dx,dy,WT0)
			WT(:,RC_D0)=WT0(:,1)
			WT(:,RC_DXX)=WT0(:,2)
			WT(:,RC_DXY)=WT0(:,3)
			WT(:,RC_DX)=WT0(:,4)
			WT(:,RC_DYY)=WT0(:,5)
			WT(:,RC_DY)=WT0(:,6)
			W0=MAXVAL(ABS(WT));
				   
			time2=MPI_Wtime()
			calc_time(1)=calc_time(1)+time2-time1
			G_E1=C_ZERO		
			G_H1=C_ZERO		
			DO K=Nfirst,Nlast
				Wm=ABS(WT(K,:)/W0(:))
!				IF (MAXVAL(Wm)>Wt_Threshold) THEN
					lm=lms(K)/r
					W(RC_DXX)=WT(K,RC_DXX)/lm
					W(RC_DYY)=WT(K,RC_DYY)/lm
					W(RC_DXY)=WT(K,RC_DXY)/lm
					W(RC_DX)=WT(K,RC_DX)*lm
					W(RC_DY)=WT(K,RC_DY)*lm
					W(RC_D0)=WT(K,RC_D0)*lm
					CALL Calc_Integral_U(anomaly,bkg,recvs(Irecv),lm,W,G_E1,G_H1,calc_time(3:4))
!				ENDIF
			ENDDO
			G_E(:,Irecv,:)=TRANSPOSE(G_E1)*s
			G_H(:,Irecv,:)=TRANSPOSE(G_H1)*s
		ENDDO
		time2=MPI_WTIME()
        calc_time(2)=calc_time(2)+time2-time3
	END SUBROUTINE
ENDMODULE
