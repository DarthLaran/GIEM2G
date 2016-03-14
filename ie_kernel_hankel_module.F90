MODULE IE_Kernel_Image_Module 
	USE CONST_MODULE
	USE DATA_TYPES_MODULE
	USE APQ_Module
	USE IntegralCodes
	USE Timer_Module 
	IMPLICIT NONE	
	PRIVATE

	INTEGER,PARAMETER::IE_DXX=1
	INTEGER,PARAMETER::IE_DXY=2
	INTEGER,PARAMETER::IE_DYX=IE_DXY
	INTEGER,PARAMETER::IE_DYY=3
	INTEGER,PARAMETER::IE_DX=4
	INTEGER,PARAMETER::IE_DY=5
	INTEGER,PARAMETER::IE_D0=6

	INTEGER,PARAMETER :: W_IND(1:6)=(/INT4DXX,INT4DXY,INT4DYY,INT4DX,INT4DY,INT4/)! codes for integrals
	PUBLIC::  Calc_Double_Integral_General
	PUBLIC::  Calc_Double_Integral_Uniform
	PUBLIC:: IE_DXX,IE_DXY, IE_DYX, IE_DYY, IE_DX, IE_DY, IE_D0, W_IND

CONTAINS



	SUBROUTINE Calc_Double_Integral_General(anomaly,bkg,dz,lms,WT,G,calc_time)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
			REAL (RealParm), INTENT(IN) ::lms
			REAL (RealParm), INTENT(IN) ::WT(1:6)
			REAL(REALPARM),INTENT(INOUT)::calc_time(2)
			COMPLEX(REALPARM),INTENT(INOUT):: G(EXX:EZZ,1:anomaly%Nz,1:anomaly%Nz)

			COMPLEX(REALPARM)::p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl),Arr(bkg%Nl,2)
			COMPLEX(REALPARM)::expz(anomaly%Nz)

			COMPLEX(REALPARM)::fcont_t(2,anomaly%Nz),t_frame(3,anomaly%Nz),pez(2,anomaly%Nz)
			COMPLEX(REALPARM)::fcont_b(anomaly%Nz),b_frame(anomaly%Nz),qez(2,anomaly%Nz)
			COMPLEX(REALPARM)::Ft(3,anomaly%Nz),Fb(3,anomaly%Nz)
			COMPLEX(REALPARM)::GII(4,anomaly%Nz)
			COMPLEX(REALPARM)::Gdiag(EXX:EZZ,anomaly%Nz)
			REAL(KIND = RealParm) ::lm2,time1,time2
			INTEGER:: Iz,l
#ifdef internal_timer
			time1=GetTime()
#endif
			lm2=lms*lms
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

			CALL CalculateMainDiagonal(GII,WT,lm2,anomaly%Nz,bkg,anomaly,Gdiag)
#ifdef internal_timer
			time2=GetTime()
			calc_time(1)=calc_time(1)+time2-time1
#endif
			DO Iz=1,anomaly%Nz
#ifdef internal_timer
				time1=GetTime()
#endif
				G(:,Iz,Iz)=Gdiag(:,Iz)+G(:,Iz,Iz)
#ifdef internal_timer
				time2=GetTime()
				calc_time(1)=calc_time(1)+time2-time1
				time1=GetTime()
#endif
				CALL GoByShellsU(Iz,anomaly%Nz,lm2,Ft(:,Iz),fcont_t,t_frame,WT,bkg%iwm,G(:,:,Iz))
				CALL GoByShellsL(Iz,anomaly%Nz,Fb(2,Iz),fcont_b,b_frame,WT,G(:,:,Iz))
#ifdef internal_timer
				time2=GetTime()	
				calc_time(2)=calc_time(2)+time2-time1
#endif
			END DO
	END SUBROUTINE

	SUBROUTINE Calc_Double_Integral_Uniform&
					&(anomaly,bkg,dz,lms,WT,G,calc_time)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
			REAL (RealParm), INTENT(IN) ::lms
			REAL (RealParm), INTENT(IN) ::WT(1:6)
			REAL(REALPARM),INTENT(INOUT)::calc_time(2)
			COMPLEX(REALPARM),INTENT(INOUT):: G(EXX:EZZ,1:anomaly%Nz,4)

			COMPLEX(REALPARM)::p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl),Arr(bkg%Nl,2)
			COMPLEX(REALPARM)::expz(anomaly%Nz)

			COMPLEX(REALPARM)::fcont_t(2,anomaly%Nz),t_frame(3,anomaly%Nz),pez(2,anomaly%Nz)
			COMPLEX(REALPARM)::fcont_b(2,anomaly%Nz),b_frame(3,anomaly%Nz),qez(2,anomaly%Nz)
			COMPLEX(REALPARM)::Ft(3,anomaly%Nz),Fb(3,anomaly%Nz)
			COMPLEX(REALPARM)::GII(4,anomaly%Nz)
			COMPLEX(REALPARM)::Gdiag(EXX:EZZ,anomaly%Nz)
			REAL(KIND = RealParm) ::lm2,time1,time2
			INTEGER:: Iz,l,Nz
#ifdef internal_timer
			time1=GetTime()
#endif
			lm2=lms*lms
			Nz=anomaly%Nz
			CALL Calc_APQ(bkg,lms,Arr,p,q,eta)
			DO Iz=1,anomaly%Nz
				l=anomaly%Lnumber(Iz)
				expz(Iz)=EXP(-eta(l)*(dz(Iz)))
			ENDDO
			CALL Calc_pexpz(bkg,anomaly,expz,eta,p,pez)

			CALL Calc_tframe_p(bkg,anomaly,expz,eta,pez,fcont_t,t_frame)
			CALL Calc_qexpz(bkg,anomaly,expz,eta,q,qez)
			CALL Calc_bframe_q_all(bkg,anomaly,expz,eta,q,qez,fcont_b,b_frame)
			CALL Calc_FtFb(bkg,anomaly,dz,expz,pez,qez,eta,Arr,Ft,Fb,GII)

			CALL CalculateMainDiagonal(GII,WT,lm2,anomaly%Nz,bkg,anomaly,Gdiag)
			G(:,:,1)=G(:,:,1)+Gdiag
			DO Iz=2,anomaly%Nz
				CALL GoByShellsFull(Iz,Iz-1,-1,anomaly%Nz,lm2,Ft(:,Iz),fcont_t,t_frame,WT,bkg%iwm,G(:,:,2))
			ENDDO
			CALL GoByShellsFull(Nz,1,-1,Nz,lm2,Ft(:,Nz),fcont_t,t_frame,WT,bkg%iwm,G(:,:,3))
			CALL GoByShellsFull(1,Nz,1,Nz,lm2,Fb(:,1),fcont_b,b_frame,WT,bkg%iwm,G(:,:,4))
#ifdef internal_timer
			time2=GetTime()
			calc_time(1)=calc_time(1)+time2-time1
#endif
	END SUBROUTINE

	SUBROUTINE CalculateMainDiagonal(GII,WT,lm2,Nz,bkg,anomaly,Gdiag)
			COMPLEX(REALPARM),INTENT(IN)::GII(4,Nz)
			REAL (RealParm), INTENT(IN) ::WT(1:6),lm2
			INTEGER,INTENT(IN)::Nz
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			COMPLEX(REALPARM),INTENT(INOUT)::Gdiag(EXX:EZZ,Nz)
			COMPLEX(REALPARM) ::f1wdxx,f1wdyy,f1wdxy
			COMPLEX(REALPARM) ::f2wdx,f2wdy,f2wdzz
			COMPLEX(REALPARM) ::f3wdxx,f3wdyy,f3wdxy
			INTEGER::Iz,l

			DO Iz=1,Nz
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

				Gdiag(EXX,Iz)=(GII(1,Iz)*WT(IE_D0)*bkg%iwm+(f1wdxx+f3wdxx)/bkg%csigma(l))
				Gdiag(EXY,Iz)=(f1wdxy+f3wdxy)/bkg%csigma(l)
				Gdiag(EXZ,Iz)=f2wdx/bkg%csigma(l)
				Gdiag(EYY,Iz)=(GII(1,Iz)*WT(IE_D0)*bkg%iwm+(f1wdyy+f3wdyy)/bkg%csigma(l))
				Gdiag(EYZ,Iz)=f2wdy/bkg%csigma(l)
				Gdiag(EZZ,Iz)=f2wdzz
			ENDDO
	ENDSUBROUTINE

	SUBROUTINE GoByShellsU(Izs,Nz,lms2,f_base,f_cont,frame,WT,iwm,G)
			INTEGER,INTENT(IN) ::Izs,Nz
			REAL (REALPARM), INTENT(IN) ::lms2
			REAL(REALPARM),INTENT(IN)::WT(6)
			COMPLEX(REALPARM),INTENT(IN)::iwm
			COMPLEX(REALPARM),INTENT(INOUT)::G(EXX:EZZ,Nz)
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

				G(EXX,Izr)=f1wdxx*frame(1,Izr)+f3wdxx*frame(2,Izr)+G(EXX,Izr)
				G(EXY,Izr)=f3wdxy*frame(2,Izr)+f1wdxy*frame(1,Izr)+G(EXY,Izr)
				G(EXZ,Izr)=f2wdx*frame(2,Izr)+G(EXZ,Izr)
				G(EYY,Izr)=f1wdyy*frame(1,Izr)+f3wdyy*frame(2,Izr)+G(EYY,Izr)
				G(EYZ,Izr)=f2wdy*frame(2,Izr)+G(EYZ,Izr)
				G(EZZ,Izr)=f2wdzz*frame(3,Izr)+G(EZZ,Izr)

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
			COMPLEX(REALPARM),INTENT(INOUT)::G(EXX:EZZ,Nz)
			COMPLEX(REALPARM) ::f2wdx,f2wdy
			INTEGER::Izr

			f2wdx=f_base*WT(IE_DX)
			f2wdy=f_base*WT(IE_DY)


			DO Izr=Izs+1,Nz

				G(EXZ,Izr)=f2wdx*frame(Izr)+G(EXZ,Izr)
				G(EYZ,Izr)=f2wdy*frame(Izr)+G(EYZ,Izr)

				f2wdx=f2wdx*f_cont(Izr)
				f2wdy=f2wdy*f_cont(Izr)

			END DO
	END SUBROUTINE

	SUBROUTINE GoByShellsFull(Izs,Nd,dir,Nz,lms2,f_base,f_cont,frame,WT,iwm,G)
			INTEGER,INTENT(IN) ::Izs,Nz,Nd,dir
			REAL (REALPARM), INTENT(IN) ::lms2
			REAL(REALPARM),INTENT(IN)::WT(6)
			COMPLEX(REALPARM),INTENT(IN)::iwm
			COMPLEX(REALPARM),INTENT(INOUT)::G(EXX:EZZ,Nz)
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


			DO Izr=Izs+dir,Nd,dir

				G(EXX,Izr)=f1wdxx*frame(1,Izr)+f3wdxx*frame(2,Izr)+G(EXX,Izr)
				G(EXY,Izr)=f3wdxy*frame(2,Izr)+f1wdxy*frame(1,Izr)+G(EXY,Izr)
				G(EXZ,Izr)=f2wdx*frame(2,Izr)+G(EXZ,Izr)
				G(EYY,Izr)=f1wdyy*frame(1,Izr)+f3wdyy*frame(2,Izr)+G(EYY,Izr)
				G(EYZ,Izr)=f2wdy*frame(2,Izr)+G(EYZ,Izr)
				G(EZZ,Izr)=f2wdzz*frame(3,Izr)+G(EZZ,Izr)

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

END MODULE 
