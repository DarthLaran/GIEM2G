MODULE IE_Kernel_Image_Module 
 	USE CONST_MODULE
	USE DATA_TYPES_MODULE
	USE APQ_Module

	USE IntegralCodes 
	IMPLICIT NONE	
	INCLUDE 'mpif.h'
	PRIVATE
	
	INTEGER,PARAMETER::IE_DXX=1
	INTEGER,PARAMETER::IE_DXY=2
	INTEGER,PARAMETER::IE_DYX=IE_DXY
	INTEGER,PARAMETER::IE_DYY=3
	INTEGER,PARAMETER::IE_DX=4
	INTEGER,PARAMETER::IE_DY=5
	INTEGER,PARAMETER::IE_D0=6

	INTEGER,PARAMETER :: W_IND(1:6)=(/INT4DXX,INT4DXY,INT4DYY,INT4DX,INT4DY,INT4/)! codes for integrals
#define no_compile
	PUBLIC::  Calc_Double_Integral_U
	PUBLIC:: IE_DXX,IE_DXY, IE_DYX, IE_DYY, IE_DX, IE_DY, IE_D0, W_IND
	CONTAINS



	SUBROUTINE Calc_Double_Integral_U(anomaly,bkg,dz,lms,WT,G,calc_time)
	    TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
        TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
        REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
		REAL (RealParm), INTENT(IN) ::lms
		REAL (RealParm), INTENT(IN) ::WT(1:6)
		REAL(REALPARM),INTENT(INOUT)::calc_time(2)
		COMPLEX(REALPARM),INTENT(INOUT):: G(EXX:EZZ,1:anomaly%Nz,1:anomaly%Nz)

		COMPLEX(REALPARM)::p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl),Arr(bkg%Nl,2)
		COMPLEX(REALPARM)::expz(anomaly%Nz)

		COMPLEX(REALPARM)::fcont_t(anomaly%Nz,2),t_frame(3,anomaly%Nz),pez(2,anomaly%Nz)
		COMPLEX(REALPARM)::fcont_b(anomaly%Nz),b_frame(anomaly%Nz),qez(2,anomaly%Nz)
		COMPLEX(REALPARM)::Ft(3,anomaly%Nz),Fb(anomaly%Nz)
		COMPLEX(REALPARM)::GII(4,anomaly%Nz)
		COMPLEX(REALPARM) ::f1wdxx,f1wdyy,f1wdxy
		COMPLEX(REALPARM) ::f2wdx,f2wdy,f2wdzz
		COMPLEX(REALPARM) ::f3wdxx,f3wdyy,f3wdxy
		REAL(KIND = RealParm) ::lm2,time1,time2
		INTEGER:: Iz,l
		time1=MPI_WTIME()
		CALL Calc_APQ(bkg,lms,Arr,p,q,eta)
		DO Iz=1,anomaly%Nz
			l=anomaly%Lnumber(Iz)
			expz(Iz)=EXP(-eta(l)*(dz(Iz)))
		ENDDO
		CALL Calc_pexpz(bkg,anomaly,expz,eta,p,pez)

		CALL Calc_tframe_p(bkg,anomaly,expz,eta,pez,fcont_t,t_frame)
		CALL Calc_qexpz(bkg,anomaly,expz,eta,q,qez)
		CALL Calc_bframe_q(bkg,anomaly,expz,eta,q,qez,fcont_b,b_frame)
		CALL Calc_FtFb(bkg,anomaly,dz,expz,pez,qez,eta,Arr,Ft,Fb,GII)
		time2=MPI_WTIME()
		calc_time(1)=calc_time(1)+time2-time1

		lm2=lms*lms
		DO Iz=1,anomaly%Nz
			time1=MPI_WTIME()
			l=anomaly%Lnumber(Iz)
			f1wdxx=GII(1,Iz)*WT(IE_DXX)*lm2
			f1wdyy=GII(1,Iz)*WT(IE_DYY)*lm2
			f1wdxy=GII(1,Iz)*WT(IE_DXY)*lm2

			f2wdx=GII(2,Iz)*WT(IE_DX)
			f2wdy=GII(2,Iz)*WT(IE_DY)
			f2wdzz=(GII(3,Iz))*WT(IE_D0)

			f3wdxx=GII(4,Iz)*WT(IE_DXX)
			f3wdyy=GII(4,Iz)*WT(IE_DYY)
			f3wdxy=GII(4,Iz)*WT(IE_DXY)


			G(EXX,Iz,Iz)=(GII(1,Iz)*WT(IE_D0)*bkg%iwm+(f1wdxx+f3wdxx)/bkg%sigma(l))+G(EXX,Iz,Iz)
			G(EXY,Iz,Iz)=(f1wdxy+f3wdxy)/bkg%sigma(l)+G(EXY,Iz,Iz)
			G(EXZ,Iz,Iz)=f2wdx/bkg%sigma(l)+G(EXZ,Iz,Iz)
			G(EYY,Iz,Iz)=(GII(1,Iz)*WT(IE_D0)*bkg%iwm+(f1wdyy+f3wdyy)/bkg%sigma(l))+G(EYY,Iz,Iz)
			G(EYZ,Iz,Iz)=f2wdy/bkg%sigma(l)+G(EYZ,Iz,Iz)
			G(EZZ,Iz,Iz)=f2wdzz+G(EZZ,Iz,Iz)

			time2=MPI_WTIME()
			calc_time(1)=calc_time(1)+time2-time1
			time1=MPI_WTIME()
			CALL GoByShellsU(Iz,anomaly%Nz,lm2,Ft(:,Iz),fcont_t,t_frame,WT,bkg%iwm,G)
			CALL GoByShellsL(Iz,anomaly%Nz,Fb(Iz),fcont_b,b_frame,WT,G)
			time2=MPI_WTIME()
			calc_time(2)=calc_time(2)+time2-time1
		END DO
	END SUBROUTINE


	SUBROUTINE GoByShellsU(Izs,Nz,lms2,f_base,f_cont,frame,WT,iwm,G)
		INTEGER,INTENT(IN) ::Izs,Nz
		REAL (REALPARM), INTENT(IN) ::lms2
		REAL(REALPARM),INTENT(IN)::WT(6)
		COMPLEX(REALPARM),INTENT(IN)::iwm
		COMPLEX(REALPARM),INTENT(INOUT)::G(EXX:EZZ,Nz,Nz)
		COMPLEX(KIND = RealParm2), INTENT(IN)::f_base(3)
		COMPLEX(KIND = RealParm2), INTENT(IN) ::f_cont(2,1:Nz),frame(3,Nz)
		COMPLEX(KIND = RealParm2) ::f1wdxx,f1wdyy,f1wdxy, f1wd0
		COMPLEX(KIND = RealParm2) ::f2wdx,f2wdy,f2wdzz
		COMPLEX(KIND = RealParm2) ::f3wdxx,f3wdyy,f3wdxy
		INTEGER::Izr
		f1wdxx=f_base(1)*(WT(IE_DXX)+WT(IE_D0))*iwm
		f1wdyy=f_base(1)*(WT(IE_DYY)+WT(IE_D0))*iwm
		f1wdxy=f_base(1)*WT(IE_DXY)*iwm

		f2wdx=f_base(2)*WT(IE_DX)
		f2wdy=f_base(2)*WT(IE_DY)
		f2wdzz=f_base(2)*WT(IE_D0)*lms2

		f3wdxx=f_base(3)*WT(IE_DXX)
		f3wdyy=f_base(3)*WT(IE_DYY)
		f3wdxy=f_base(3)*WT(IE_DXY)


		DO Izr=Izs-1,1,-1

			G(EXX,Izr,Izs)=f1wdxx*frame(1,Izr)+f3wdxx*frame(2,Izr)+G(EXX,Izr,Izs)
			G(EXY,Izr,Izs)=f3wdxy*frame(2,Izr)+f1wdxy*frame(1,Izr)+G(EXY,Izr,Izs)
			G(EXZ,Izr,Izs)=f2wdx*frame(2,Izr)+G(EXZ,Izr,Izs)
			G(EYY,Izr,Izs)=f1wdyy*frame(1,Izr)+f3wdyy*frame(2,Izr)+G(EYY,Izr,Izs)
			G(EYZ,Izr,Izs)=f2wdy*frame(2,Izr)+G(EYZ,Izr,Izs)
			G(EZZ,Izr,Izs)=f2wdzz*frame(3,Izr)+G(EZZ,Izr,Izs)

			f1wdxx=f1wdxx*f_cont(1,Izr)
			f1wdyy=f1wdyy*f_cont(1,Izr)
			f1wdxy=f1wdxy*f_cont(1,Izr)
							  
			f2wdx=f2wdx*f_cont(2,Izr)
			f2wdy=f2wdy*f_cont(2,Izr)
			f2wdzz=f2wdzz*f_cont(2,Izr)
							  
			f3wdxx=f3wdxx*f_cont(2,Izr)
			f3wdyy=f3wdyy*f_cont(2,Izr)
			f3wdxy=f3wdxy*f_cont(2,Izr)


		END DO
	END SUBROUTINE
	
	SUBROUTINE GoByShellsL(Izs,Nz,f_base,f_cont,frame,WT,G)
		INTEGER,INTENT(IN) ::Izs,Nz
		COMPLEX(REALPARM), INTENT(IN)::f_base
		COMPLEX(REALPARM), INTENT(IN) ::f_cont(1:Nz),frame(1:Nz)
		REAL(REALPARM),INTENT(IN)::WT(6)
		COMPLEX(REALPARM),INTENT(INOUT)::G(EXX:EZZ,Nz,Nz)
		COMPLEX(REALPARM) ::f2wdx,f2wdy
		INTEGER::Izr

		f2wdx=f_base*WT(IE_DX)
		f2wdy=f_base*WT(IE_DY)


		DO Izr=Izs+1,Nz

			G(EXZ,Izr,Izs)=f2wdx*frame(Izr)+G(EXZ,Izr,Izs)
			G(EYZ,Izr,Izs)=f2wdy*frame(Izr)+G(EYZ,Izr,Izs)
							  
			f2wdx=f2wdx*f_cont(Izr)
			f2wdy=f2wdy*f_cont(Izr)

		END DO
	END SUBROUTINE
	
	
END MODULE 
