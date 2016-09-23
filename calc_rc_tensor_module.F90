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
!along with GFMRES.  If not, see <http://www.gnu.org/licenses/>.
!
!
!

MODULE Calc_RC_Tensor_Module
	USE CONST_MODULE
	USE FFTW3
	USE MPI_MODULE

	USE Timer_Module 
	USE DATA_TYPES_MODULE
	USE RC_Kernel_Image_Module
	USE CONTINUATION_FUNCTION_MODULE

	USE IntegralCodes 
	USE VolumeBesselTransforms
        USE LOGGER_MODULE
	IMPLICIT NONE
	REAL(REALPARM)::lms(Nfirst:Nlast)
	REAL(REALPARM),PARAMETER::Wt_Threshold=1d-12
	PUBLIC::CalcRecalculationGreenTensor
CONTAINS
	SUBROUTINE CalcRecalculationGreenTensor(matrix,bkg,anomaly)
		TYPE(RC_OPERATOR),INTENT(INOUT)::matrix
		TYPE (BKG_DATA_TYPE),TARGET,INTENT(INOUT)::bkg
		TYPE (ANOMALY_TYPE),TARGET,INTENT(INOUT)::anomaly
		COMPLEX(REALPARM),POINTER::G_E(:,:,:,:,:)
		COMPLEX(REALPARM),POINTER::G_H(:,:,:,:,:)
		INTEGER::Ix,Iy,ly,Ic
		INTEGER(MPI_CTL_KIND)::IERROR
		INTEGER::Nx,Ny,Ny_loc,Ny_offset,Nz,Nx2
		REAL(REALPARM)::dx,dy,rmin
		INTEGER::xfirst,xlast
		REAL(8)::time(4),time_all
		REAL(8)::time_max(4),time_min(4)
		REAL(8)::time1,time2,time3,time4
		time_all=GetTime()
		CALL PRINT_BORDER
		CALL LOGGER('Start RC kernel calculation')
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
		dx=anomaly%dx
		dy=anomaly%dy
		rmin=SQRT(dx*dx+dy*dy)/2
		CALL VBTransformInit(rmin,lms)
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
#ifdef internal_timer
				!$OMP DO SCHEDULE(GUIDED)&
				!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
#else
				!$OMP DO SCHEDULE(GUIDED)
#endif
				DO Ix=0,Nx-1
					CALL CalcRecalculationGreenTensor3dElement(Ix,ly,G_H(:,:,:,Ix,Iy),&
					&G_E(:,:,:,Ix,Iy),time,anomaly,bkg,matrix%recvs,matrix%Nr)
				ENDDO
				!$OMP ENDDO
#ifdef internal_timer
				!$OMP DO SCHEDULE(GUIDED)&
				!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
#else
				!$OMP DO SCHEDULE(GUIDED)
#endif
				DO Ix=Nx+1,2*Nx-1
					CALL CalcRecalculationGreenTensor3dElement(Ix-2*Nx,ly,G_H(:,:,:,Ix,Iy),&
					&G_E(:,:,:,Ix,Iy),time,anomaly,bkg,matrix%recvs,matrix%Nr)
				ENDDO
				!$OMP ENDDO
				!$OMP WORKSHARE
					G_E(:,:,:,Nx,Iy)=C_ZERO
					G_H(:,:,:,Nx,Iy)=C_ZERO
				!$OMP END WORKSHARE
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
		ENDDO
#ifdef internal_timer
		CALL MPI_REDUCE(time,time_min,4,MPI_DOUBLE,MPI_MIN,matrix%master_proc,matrix%matrix_comm,IERROR)
		CALL MPI_REDUCE(time,time_max,4,MPI_DOUBLE,MPI_MAX,matrix%master_proc,matrix%matrix_comm,IERROR)
#endif
		time_all=GetTime()-time_all
#ifdef internal_timer
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
		matrix%counter%tensor_calc(COUNTER_ALL)=time_all	
#else
	        CALL LOGGER('RC kernel has been computed.')
                CALL PRINT_CALC_TIME(CALC_IE_FULL_TIME,time_all)
                CALL PRINT_BORDER
#endif
	END SUBROUTINE
!---------------------------------------------------------------------------------------------------------------------------------------------!
	SUBROUTINE CalcRecalculationGreenTensor3dElement(lx,ly,G_H,G_E,calc_time,anomaly,bkg,recvs,Nr)
		IMPLICIT NONE
		INTEGER, INTENT(IN)::lx,ly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		INTEGER,INTENT(IN)::Nr
		TYPE (RECEIVER_TYPE),INTENT(IN)::recvs(Nr)
		COMPLEX(REALPARM),INTENT(OUT)::G_H(anomaly%Nz,Nr,RHXX:RHZY)
		COMPLEX(REALPARM),INTENT(OUT)::G_E(anomaly%Nz,Nr,REXX:REZZ)
		REAL(REALPARM),INTENT(INOUT)::calc_time(4)
		COMPLEX(REALPARM)::G_E1(REXX:REZZ,anomaly%Nz)
		COMPLEX(REALPARM)::G_H1(RHXX:RHZY,anomaly%Nz)
		REAL(REALPARM2)::r,WT(Nfirst:Nlast,1:6),W(1:6),dx,dy,s
		REAL(REALPARM)::x,y,lm,W0(1:6),Wm(1:6),time1,time2,time3
		INTEGER::K,Iz,IERROR,I,Ik,Iz0,Ic,Iwt(6),Irecv
		REAL(REALPARM2)::WT0(Nfirst:Nlast,1:6)
		REAL(REALPARM2)::WT1(Nfirst:Nlast,1:6)
#ifdef internal_timer
		time3=GetTime()
#endif
		dx=anomaly%dx
		dy=anomaly%dy
		s=1d0/4.0_REALPARM/PI!*dx*dy
		DO Irecv=1,Nr
#ifdef internal_timer
			time1=GetTime()
#endif
			x=lx*dx-dx/2d0+recvs(Irecv)%x_shift
			y=ly*dy-dy/2d0+recvs(Irecv)%y_shift
			r=SQRT(x*x+y*y)
			CALL VBTransformWeightsAllSingleInt(x,y,dx,dy,WT0)
!		CALL VBTransformWeightsAllDoubleInt(x,y,dx,dy,WT0)
			WT(:,RC_D0)=WT0(:,1)
			WT(:,RC_DXX)=WT0(:,2)
			WT(:,RC_DXY)=WT0(:,3)
			WT(:,RC_DX)=WT0(:,4)
			WT(:,RC_DYY)=WT0(:,5)
			WT(:,RC_DY)=WT0(:,6)
			W0=MAXVAL(ABS(WT));
				   
#ifdef internal_timer
			time2=GetTime()
			calc_time(1)=calc_time(1)+time2-time1
#endif
			G_E1=C_ZERO		
			G_H1=C_ZERO		
			DO K=Nfirst,Nlast
				Wm=ABS(WT(K,:)/W0(:))
				IF (MAXVAL(Wm)>Wt_Threshold) THEN
					lm=lms(K)!/r
					W(RC_DXX)=WT(K,RC_DXX)/lm
					W(RC_DYY)=WT(K,RC_DYY)/lm
					W(RC_DXY)=WT(K,RC_DXY)/lm
					W(RC_DX)=WT(K,RC_DX)*lm
					W(RC_DY)=WT(K,RC_DY)*lm
					W(RC_D0)=WT(K,RC_D0)*lm
					CALL Calc_Integral_U(anomaly,bkg,recvs(Irecv),lm,W,G_E1,G_H1,calc_time(3:4))
				ENDIF
			ENDDO
			G_E(:,Irecv,:)=TRANSPOSE(G_E1)*s
			G_H(:,Irecv,:)=TRANSPOSE(G_H1)*s
		ENDDO
#ifdef internal_timer
		time2=GetTime()
        calc_time(2)=calc_time(2)+time2-time3
#endif
	END SUBROUTINE
ENDMODULE
