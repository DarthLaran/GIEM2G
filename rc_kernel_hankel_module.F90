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

MODULE RC_Kernel_Image_Module 
 	USE CONST_MODULE
	USE DATA_TYPES_MODULE
	USE APQ_Module

	USE IntegralCodes
!	USE MPI_MODULE
	USE Timer_Module 
	IMPLICIT NONE	
	PRIVATE
	
	INTEGER,PARAMETER::RC_DXX=1
	INTEGER,PARAMETER::RC_DXY=2
	INTEGER,PARAMETER::RC_DYX=RC_DXY
	INTEGER,PARAMETER::RC_DYY=3
	INTEGER,PARAMETER::RC_DX=4
	INTEGER,PARAMETER::RC_DY=5
	INTEGER,PARAMETER::RC_D0=6

	INTEGER,PARAMETER :: WR_IND(1:6)=(/INT2DXX,INT2DXY,INT2DYY,INT2DX,INT2DY,INT2/)! codes for integrals


	INTEGER,PARAMETER::G1=1
	INTEGER,PARAMETER::G3DZ=2
	INTEGER,PARAMETER::G2DZ=3
	INTEGER,PARAMETER::G2DZ0=4
	INTEGER,PARAMETER::G2G2DZZ=5
	INTEGER,PARAMETER::G3=6
	INTEGER,PARAMETER::G1DZ=7
	INTEGER,PARAMETER::G2=8
#define no_compile
	PUBLIC::  Calc_Integral_U
	PUBLIC:: RC_DXX,RC_DXY, RC_DYX, RC_DYY, RC_DX, RC_DY, RC_D0, WR_IND
	CONTAINS



	SUBROUTINE Calc_Integral_U(anomaly,bkg,recv,lms,WT,G_E,G_H,calc_time)
	    TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
        TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		TYPE (RECEIVER_TYPE),INTENT(IN)::recv
		REAL (RealParm), INTENT(IN) ::lms
		REAL (RealParm), INTENT(IN) ::WT(1:6)
		COMPLEX(REALPARM),INTENT(INOUT):: G_E(REXX:REZZ,1:anomaly%Nz)
		COMPLEX(REALPARM),INTENT(INOUT):: G_H(RHXX:RHZY,1:anomaly%Nz)
		REAL(REALPARM),INTENT(INOUT)::calc_time(2)

		COMPLEX(REALPARM)::A(bkg%Nl,bkg%Nl,2),p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl),Arr(bkg%Nl,2)
		COMPLEX(REALPARM)::expz(anomaly%Nz)

		COMPLEX(REALPARM)::pez(2,anomaly%Nz)
		COMPLEX(REALPARM)::qez(2,anomaly%Nz)
		COMPLEX(REALPARM)::w1(2,anomaly%Nz)
		COMPLEX(REALPARM)::w2(2,anomaly%Nz)
		COMPLEX(REALPARM)::Ftb(3,anomaly%Nz)
		COMPLEX(REALPARM)::G(G1:G2,anomaly%Nz),G11
		REAL(KIND = RealParm) ::lm2,time1,time2
		INTEGER::N,Nz
		INTEGER:: Iz,l
#ifdef internal_timer
		time1=GetTime()
#endif
		lm2=lms*lms
		CALL Calc_APQ(bkg,lms,Arr,p,q,eta,A)
		DO Iz=1,anomaly%Nz
			l=anomaly%Lnumber(Iz)
			expz(Iz)=EXP(-eta(l)*(anomaly%z(Iz)-anomaly%z(Iz-1)))
		ENDDO
		Nz=anomaly%Nz
		CALL Calc_Recv_Part(bkg,anomaly,p,q,eta,recv,w1,w2)
		Ftb=C_ZERO
		IF (recv%anom_cell<1) THEN
			CALL Calc_qexpz(bkg,anomaly,expz,eta,q,qez)
			CALL Calc_Ftb(bkg,anomaly,expz,eta,qez,1,Nz,Ftb,C_ONE)
			CALL Calc_G_not_match(bkg,anomaly,A,eta,recv,lm2,Ftb,w1,w2,G)
		ELSEIF (recv%anom_cell>anomaly%Nz) THEN
			CALL Calc_pexpz(bkg,anomaly,expz,eta,p,pez)
			CALL Calc_Ftb(bkg,anomaly,expz,eta,pez,1,Nz,Ftb,-C_ONE)
			CALL Calc_G_not_match(bkg,anomaly,A,eta,recv,lm2,Ftb,w1,w2,G)
		ELSE
			N=recv%anom_cell
!			CALL Calc_pexpz(bkg,anomaly,expz,eta,p,1,N,pez)
!			CALL Calc_qexpz(bkg,anomaly,expz,eta,q,N+1,Nz,qez)

			CALL Calc_pexpz(bkg,anomaly,expz,eta,p,pez)
			CALL Calc_qexpz(bkg,anomaly,expz,eta,q,qez)

			CALL Calc_Ftb(bkg,anomaly,expz,eta,pez,1,N,Ftb,-C_ONE)
			CALL Calc_Ftb(bkg,anomaly,expz,eta,qez,N+1,Nz,Ftb,C_ONE)
			CALL Calc_G_not_match(bkg,anomaly,A,eta,recv,lm2,Ftb,w1,w2,G)
			!CALL Calc_G_match(bkg,anomaly,A,p,eta,recv,N,Fb,G)
			
		ENDIF
#ifdef internal_timer
		time2=GetTime()
		calc_time(2)=calc_time(2)+time2-time1
#endif
		DO Iz=1,anomaly%Nz
!--------------------------------- G_E ----------------------------------------------!
			G11=G(G1,Iz)*bkg%iwm

			G_E(REXX,Iz)=G_E(REXX,Iz)+G11*WT(RC_D0)+(G(G3DZ,Iz)+G11)*WT(RC_DXX)
			G_E(REXY,Iz)=G_E(REXY,Iz)+(G(G3DZ,Iz)+G11)*WT(RC_DXY)

			G_E(REXZ,Iz)=G_E(REXZ,Iz)+G(G2DZ,Iz)*WT(RC_DX)

			G_E(REYY,Iz)=G_E(REYY,Iz)+G11*WT(RC_D0)+(G(G3DZ,Iz)+G11)*WT(RC_DYY)

			G_E(REYZ,Iz)=G_E(REYZ,Iz)+G(G2DZ,Iz)*WT(RC_DY)


			G_E(REZX,Iz)=G_E(REZX,Iz)+G(G2DZ0,Iz)*WT(RC_DX)
			G_E(REZY,Iz)=G_E(REZY,Iz)+G(G2DZ0,Iz)*WT(RC_DY)

			G_E(REZZ,Iz)=G_E(REZZ,Iz)+G(G2G2DZZ,Iz)*WT(RC_D0)
!--------------------------------- G_H ----------------------------------------------!
			G_H(RHXX,Iz)=G_H(RHXX,Iz)+G(G3,Iz)*WT(RC_DXY)
			G_H(RHXY,Iz)=G_H(RHXY,Iz)+G(G3,Iz)*WT(RC_DYY)-G(G1DZ,Iz)*WT(RC_D0)
			G_H(RHXZ,Iz)=G_H(RHXZ,Iz)+G(G2,Iz)*WT(RC_DY)

			G_H(RHYX,Iz)=G_H(RHYX,Iz)+G(G1DZ,Iz)*WT(RC_D0)-G(G3,Iz)*WT(RC_DXX)
			G_H(RHYZ,Iz)=G_H(RHYZ,Iz)-G(G2,Iz)*WT(RC_DX)

			G_H(RHZX,Iz)=G_H(RHZX,Iz)-G(G1,Iz)*WT(RC_DY)
			G_H(RHZY,Iz)=G_H(RHZY,Iz)+G(G1,Iz)*WT(RC_DX)

		END DO
#ifdef internal_timer
		time1=GetTime()
		calc_time(1)=calc_time(1)+time1-time2
#endif
	END SUBROUTINE

	SUBROUTINE	 Calc_G_not_match(bkg,anomaly,A,eta,recv,lm2,Ftb,w1,w2,G)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::A(bkg%Nl,bkg%Nl,2),eta(bkg%Nl)
		TYPE (RECEIVER_TYPE),INTENT(IN)::recv
		COMPLEX(REALPARM),INTENT(IN)::w1(2,anomaly%Nz),w2(2,anomaly%Nz)
		REAL(REALPARM),INTENT(IN)::lm2
		COMPLEX(REALPARM),INTENT(IN)::Ftb(3,anomaly%Nz)
		COMPLEX(REALPARM),INTENT(INOUT)::G(8,anomaly%Nz)
		INTEGER::I,l,rl
		COMPLEX(REALPARM)::e1,e2,eta_r,a1,a2
		COMPLEX(REALPARM)::v1(2),v2(2),f(3)
		rl=recv%recv_layer
		DO I=1,anomaly%Nz
			l=anomaly%Lnumber(I)
			a1=A(rl,l,1)
			a2=A(rl,l,2)
			f=Ftb(:,I)
			v1=w1(:,I)
			v2=w2(:,I)
			G(G1,I)=a1*f(1)*v1(1)

			G(G3DZ,I)=(a2*f(3)*v2(2))/bkg%csigma(rl)

			G(G2DZ,I)=(a2*f(2)*v2(2))/bkg%csigma(rl)

			G(G2DZ0,I)=(a2*f(3)*v2(1))/bkg%csigma(rl)

			G(G2G2DZZ,I)=lm2*(a2*f(2)*v2(1))/bkg%csigma(rl)

			G(G3,I)=(a2*f(3)*v2(1)-a1*f(1)*v1(2))
			G(G1DZ,I)=a1*f(1)*v1(2)
			G(G2,I)=a2*f(2)*v2(1)
		ENDDO
	END SUBROUTINE
	SUBROUTINE Calc_Recv_Part(bkg,anomaly,p,q,eta,recv,w1,w2)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl)
		TYPE (RECEIVER_TYPE),INTENT(IN)::recv
		COMPLEX(REALPARM),INTENT(INOUT)::w1(2,anomaly%Nz),w2(2,anomaly%Nz)
		INTEGER::I,l,rl,N,N1,N2
		COMPLEX(REALPARM)::e1,e2,eta_r,a1,a2
		COMPLEX(REALPARM)::eta_zr,eta_dr
		rl=recv%recv_layer
		eta_r=eta(rl)
		eta_zr=eta_r*recv%zrecv
		IF (rl>1) THEN
			eta_dr=eta_r*bkg%depth(rl-1)*C_TWO
		ELSE
			eta_dr=C_ZERO
		ENDIF
		N=recv%anom_cell
		N1=MAX(1,N+1)
		N2=MIN(N,anomaly%Nz)
		DO I=N1,anomaly%Nz
			l=anomaly%Lnumber(I)
			e1=EXP(-anomaly%z(I-1)*eta(l)+eta_zr)
			e2=EXP(-anomaly%z(I-1)*eta(l)-eta_zr+eta_dr)
			w1(1,I)=(e2*p(rl,1)+e1)
			w1(2,I)=(-e2*p(rl,1)+e1)*eta_r
			w2(1,I)=(e2*p(rl,2)+e1)
			w2(2,I)=(-e2*p(rl,2)+e1)*eta_r
		ENDDO
		IF (rl<bkg%Nl) THEN
			eta_dr=eta_r*bkg%depth(rl)*C_TWO
		ELSE
			eta_dr=C_ONE*HUGE(recv%zrecv)
		ENDIF
		DO I=1,N2
			l=anomaly%Lnumber(I)
			e1=EXP(anomaly%z(I)*eta(l)-eta_zr)
			e2=EXP(anomaly%z(I)*eta(l)+eta_zr-eta_dr)
			w1(1,I)=(e2*q(rl,1)+e1)
			w1(2,I)=(e2*q(rl,1)-e1)*eta_r
			w2(1,I)=(e2*q(rl,2)+e1)
			w2(2,I)=(e2*q(rl,2)-e1)*eta_r
		ENDDO
	END SUBROUTINE

	SUBROUTINE Calc_Ftb(bkg,anomaly,expz,eta,pqexpz,N1,N2,Ftb,s)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),pqexpz(2,anomaly%Nz)
		COMPLEX(REALPARM),INTENT(IN)::eta(bkg%Nl)
		COMPLEX(REALPARM),INTENT(IN)::s
		INTEGER,INTENT(IN)::N1,N2
		COMPLEX(REALPARM),INTENT(INOUT)::Ftb(3,anomaly%Nz)
		COMPLEX(REALPARM)::e1(2),e2(2)
		COMPLEX(REALPARM)::a(2),a1(2)
		INTEGER::I,l
			
		DO I=N1,N2
			l=anomaly%Lnumber(I)
			a=(C_ONE-expz(I))
			a1=a/eta(l)	
			e1=C_ONE+pqexpz(:,I)*expz(I)
			e2=(C_ONE-pqexpz(:,I)*expz(I))*s

			Ftb(1:2,I)=a1*e1
			Ftb(3,I)=a(2)*e2(2)
		ENDDO
	ENDSUBROUTINE
	
END MODULE 
