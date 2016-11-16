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


			!$OMP PARALLEL PRIVATE(Ix,Iy,ly),DEFAULT(SHARED)&
			!$OMP &FIRSTPRIVATE(time)
#ifdef internal_timer
			!$OMP DO SCHEDULE(GUIDED)&
			!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
#else
			!$OMP DO SCHEDULE(GUIDED)
#endif
				DO Ix=0,Nx2-1
                        		DO Iy=Ny_offset,Ny_offset+Ny_loc-1
                                		ly=MODULO(Iy,Ny)
	        				CALL CalcIntegralGreenTensor3dElementGeneral(Ix+xfirst,ly,G_symm(:,:,Iy,Ix),&
		        			&G_asym(:,:,:,Iy,Ix),time,anomaly,bkg,ie_op%dz)
			        	ENDDO
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
	SUBROUTINE   CalcGeneralTensorAlongXY(xfirst,Nx2,Ny_offset,Ny_loc,G_symm,G_asym,time,anomaly,bkg,dz)
			INTEGER,INTENT(IN)::xfirst,Nx2,Ny_offset,Ny_loc
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_symm(:,:,:,:)
			COMPLEX(REALPARM),POINTER,INTENT(IN)::G_asym(:,:,:,:,:)
			REAL(DOUBLEPARM),INTENT(INOUT)::time(4)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			REAL(REALPARM),POINTER,INTENT(IN)::dz(:)
			COMPLEX(REALPARM),POINTER::G_symm1(:,:,:,:)
			COMPLEX(REALPARM),POINTER::G_asym1(:,:,:,:,:)
			REAL(DOUBLEPARM)::time1,time2,time3,time4
			INTEGER::Ix,Iy,ly,Ny
			G_symm1(1:,S_EXX:,Ny_offset:,0:)=>G_symm
			G_asym1(1:,1:,A_EXZ:,Ny_offset:,0:)=>G_asym
                        Ny=anomaly%Ny
			!$OMP PARALLEL PRIVATE(Ix,Iy,ly),DEFAULT(SHARED)&
			!$OMP &FIRSTPRIVATE(time)
#ifdef internal_timer
			!$OMP DO SCHEDULE(GUIDED)&
			!$OMP& REDUCTION (MAX:time1,time2,time3,time4)
#else
			!$OMP DO SCHEDULE(GUIDED)
#endif
				DO Ix=0,Nx2-1
                        		DO Iy=Ny_offset,Ny_offset+Ny_loc-1
                                		ly=MODULO(Iy,Ny)
	        				CALL CalcIntegralGreenTensor3dElementGeneral(Ix+xfirst,ly,G_symm1(:,:,Iy,Ix),&
		        			&G_asym1(:,:,:,Iy,Ix),time,anomaly,bkg,dz)
			        	ENDDO
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
		COMPLEX(REALPARM),INTENT(OUT)::G_symm(anomaly%Nz*(anomaly%Nz+1)/2,S_EXX:S_EZZ)
		COMPLEX(REALPARM),INTENT(OUT)::G_asym(anomaly%Nz,anomaly%Nz,A_EXZ:A_EYZ)
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
#ifdef internal_timer
		time2=GetTime()
		calc_time(2)=calc_time(2)+time2-time1
#endif
	END SUBROUTINE

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
		REAL(REALPARM),INTENT(OUT)::WT(Nfirst:Nlast,1:6)
		REAL(REALPARM),INTENT(OUT)::W0(1:6)
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
ENDMODULE
