MODULE APQ_Module
	USE Const_module
	USE DATA_TYPES_MODULE
	IMPLICIT NONE
	INTERFACE Calc_pexpz
		MODULE PROCEDURE Calc_pexpz_All,Calc_pexpz_Seg
	END INTERFACE
	INTERFACE Calc_qexpz
		MODULE PROCEDURE Calc_qexpz_All,Calc_qexpz_Seg
	END INTERFACE
CONTAINS

	SUBROUTINE Calc_FtFb(bkg,anomaly,dz,expz,pexpz,qexpz,eta,Arr,Ft,Fb,GII)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		REAL(REALPARM),INTENT(IN)::dz(anomaly%Nz)
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),pexpz(2,anomaly%Nz),qexpz(2,anomaly%Nz)
		COMPLEX(REALPARM),INTENT(IN)::eta(bkg%Nl),Arr(bkg%Nl,2)
		COMPLEX(REALPARM),INTENT(INOUT)::Ft(3,anomaly%Nz),Fb(3,anomaly%Nz),GII(4,anomaly%Nz)
		COMPLEX(REALPARM)::Fb3(3)
		COMPLEX(REALPARM)::p1(2),p2(2),q1(2),q2(2)
		COMPLEX(REALPARM)::pe1(2),pe2(2),qe1(2),qe2(2)
		COMPLEX(REALPARM)::a(2),a1(2),a2(2),eta2,w(2)
		INTEGER::I,l
			
		DO I=1,anomaly%Nz
			l=anomaly%Lnumber(I)
			a=(C_ONE-expz(I))*Arr(l,:)
			a1=a/eta(l)	
			a2=a1/eta(l)	
			eta2=eta(l)*eta(l)
			p1=C_ONE+pexpz(:,I)
			p2=C_ONE-pexpz(:,I)

			pe1=C_ONE+pexpz(:,I)*expz(I)
			pe2=-C_ONE+pexpz(:,I)*expz(I)

			q1=C_ONE+qexpz(:,I)
			q2=C_ONE-qexpz(:,I)

			qe1=C_ONE+qexpz(:,I)*expz(I)
			qe2=C_ONE-qexpz(:,I)*expz(I)

			Ft(1:2,I)=a1*p1*qe1
			Ft(3,I)=a(2)*p1(2)*qe2(2)

			Fb3(1:2)=a1*pe1*q1
			Fb3(3)=a(2)*pe2(2)*q1(2)

			w=-a2*(p2*qe1+pe1*q2)
			GII(1,I)=w(1)+2e0_REALPARM/eta2*dz(I)

			GII(2,I)=Fb3(2)-Ft(2,I)


			GII(3,I)=(w(2)+2e0_REALPARM/eta2*dz(I))*bkg%iwm+eta2/bkg%csigma(l)*w(2)

			GII(4,I)=Fb3(3)-Ft(3,I)-eta2*w(1)
			Fb(:,I)=Fb3(:)
		ENDDO
	ENDSUBROUTINE



	SUBROUTINE Calc_tframe_p(bkg,anomaly,expz,eta,pexpz,fcont,t_frame)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),eta(bkg%Nl),pexpz(2,anomaly%Nz)
		COMPLEX(REALPARM),INTENT(INOUT)::fcont(2,anomaly%Nz),t_frame(3,anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a(2),b(2),c(2),d(2)
		INTEGER::I,l
		
		DO I=1,anomaly%Nz-1
			l=anomaly%Lnumber(I)
			a=pexpz(:,I)
			b=a*expz(I)
			c=b*expz(I)
			
			a=C_ONE+a
			b=C_ONE+b
			c=C_ONE+c

			w=(C_ONE-expz(I))/eta(l)

			fcont(:,I)=expz(I)*a/c
		
			t_frame(1,I)=w*b(1)/c(1)
			t_frame(2,I)=C_ONE-fcont(2,I)
			t_frame(3,I)=w*b(2)/c(2)
			t_frame(2:3,I)=t_frame(2:3,I)/bkg%csigma(l)
		ENDDO
	ENDSUBROUTINE

	SUBROUTINE Calc_pexpz_All(bkg,anomaly,expz,eta,p,pexpz)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),p(bkg%Nl,2),eta(bkg%Nl)
		COMPLEX(REALPARM),INTENT(INOUT)::pexpz(2,anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a(2),b(2),c(2),d(2)
		INTEGER::I,l
		l=anomaly%Lnumber(0)
		IF (l>1) THEN
			e=EXP(-2e0_REALPARM*eta(l)*(anomaly%z(0)-bkg%depth(l-1)))
		ELSEIF(l==1) THEN
			e=EXP(-2e0_REALPARM*eta(1)*(anomaly%z(0)))
		ELSE
			e=C_ONE
		ENDIF
		
		DO I=1,anomaly%Nz
			IF (l/=anomaly%Lnumber(I)) THEN
				e=C_ONE
			ENDIF
			l=anomaly%Lnumber(I)
			a=e*p(l,:)
			pexpz(:,I)=a
			e=e*expz(I)*expz(I)
		ENDDO
	ENDSUBROUTINE
	SUBROUTINE Calc_bframe_q(bkg,anomaly,expz,eta,q,qexpz,fcont,b_frame)
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),q(bkg%Nl,2),eta(bkg%Nl),qexpz(2,anomaly%Nz)
		COMPLEX(REALPARM),INTENT(INOUT)::fcont(anomaly%Nz),b_frame(anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a,b,c,d
		INTEGER::I,l
		DO I=anomaly%Nz,2,-1
			l=anomaly%Lnumber(I)
			a=qexpz(2,I)
			d=C_ONE+a*expz(I)*expz(I)
			a=C_ONE+a
			fcont(I)=expz(I)*a/d
			b_frame(I)=(fcont(I)-C_ONE)/bkg%csigma(l)
		ENDDO
	ENDSUBROUTINE

	SUBROUTINE Calc_bframe_q_all(bkg,anomaly,expz,eta,q,qexpz,fcont,b_frame)
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),q(bkg%Nl,2),eta(bkg%Nl),qexpz(2,anomaly%Nz)
		COMPLEX(REALPARM),INTENT(INOUT)::fcont(2,anomaly%Nz),b_frame(3,anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a(2),b(2),c(2),d(2)
		INTEGER::I,l
		DO I=anomaly%Nz,2,-1
			l=anomaly%Lnumber(I)
			a=qexpz(:,I)
			b=a*expz(I)
			c=b*expz(I)
			
			b=C_ONE+b
			c=C_ONE+c

			w=-(C_ONE-expz(I))/eta(l)
			d=C_ONE+a*expz(I)*expz(I)
			a=C_ONE+a
			fcont(:,I)=expz(I)*a/d

			b_frame(1,I)=w*b(1)/c(1)
                        b_frame(2,I)=(fcont(2,I)-C_ONE)
			b_frame(3,I)=w*b(2)/c(2)

			b_frame(2:3,I)=b_frame(2:3,I)/bkg%csigma(l)
		ENDDO
	ENDSUBROUTINE
	SUBROUTINE Calc_qexpz_All(bkg,anomaly,expz,eta,q,qexpz)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),q(bkg%Nl,2),eta(bkg%Nl)
		COMPLEX(REALPARM),INTENT(INOUT)::qexpz(2,anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a(2),b(2),c(2),d(2)
		INTEGER::I,l
		l=anomaly%Lnumber(anomaly%Nz)
		IF (l/=bkg%Nl) THEN
			e=EXP(-2e0_REALPARM*eta(l)*(bkg%depth(l)-anomaly%z(anomaly%Nz)))
		ELSE
			e=C_ZERO
		ENDIF
		DO I=anomaly%Nz,1,-1
			IF (l/=anomaly%Lnumber(I)) THEN
				e=C_ONE
			ENDIF
			l=anomaly%Lnumber(I)
			a=e*q(l,:)
			qexpz(:,I)=a
			e=e*expz(I)*expz(I)
		ENDDO
	ENDSUBROUTINE
	SUBROUTINE Calc_pexpz_Seg(bkg,anomaly,expz,eta,p,N1,N2,pexpz)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),p(bkg%Nl,2),eta(bkg%Nl)
		INTEGER,INTENT(IN)::N1,N2
		COMPLEX(REALPARM),INTENT(INOUT)::pexpz(2,1:anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a(2),b(2),c(2),d(2)
		INTEGER::I,l
		l=anomaly%Lnumber(N1-1)
		IF (l>1) THEN
			e=EXP(-2e0_REALPARM*eta(l)*(anomaly%z(N1-1)-bkg%depth(l-1)))
		ELSEIF(l==1) THEN
			e=EXP(-2e0_REALPARM*eta(1)*(anomaly%z(N1-1)))
		ELSE
			e=C_ONE
		ENDIF
		
		DO I=1,N2
			IF (l/=anomaly%Lnumber(I)) THEN
				e=C_ONE
			ENDIF
			l=anomaly%Lnumber(I)
			a=e*p(l,:)
			pexpz(:,I)=a
			e=e*expz(I)*expz(I)
		ENDDO
	ENDSUBROUTINE

	SUBROUTINE Calc_qexpz_Seg(bkg,anomaly,expz,eta,q,N1,N2,qexpz)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(IN)::expz(anomaly%Nz),q(bkg%Nl,2),eta(bkg%Nl)
		INTEGER,INTENT(IN)::N1,N2
		COMPLEX(REALPARM),INTENT(INOUT)::qexpz(2,1:anomaly%Nz)
		COMPLEX(REALPARM)::e,e1,e2,w,a(2),b(2),c(2),d(2)
		INTEGER::I,l
		l=anomaly%Lnumber(N2)
		IF (l/=bkg%Nl) THEN
			e=EXP(-2e0_REALPARM*eta(l)*(bkg%depth(l)-anomaly%z(N2)))
		ELSE
			e=C_ZERO
		ENDIF
		DO I=N2,N1,-1
			IF (l/=anomaly%Lnumber(I)) THEN
				e=C_ONE
			ENDIF
			l=anomaly%Lnumber(I)
			a=e*q(l,:)
			qexpz(:,I)=a
			e=e*expz(I)*expz(I)
		ENDDO
	ENDSUBROUTINE
	SUBROUTINE Calc_Apq(bkg,lms,Arr,p,q,eta,A)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		COMPLEX(REALPARM),INTENT(INOUT)::p(bkg%Nl,2),q(bkg%Nl,2)
		COMPLEX(REALPARM),INTENT(INOUT)::eta(bkg%Nl),Arr(bkg%Nl,2)
		COMPLEX(REALPARM),OPTIONAL,INTENT(INOUT)::A(bkg%Nl,bkg%Nl,2)
		REAL(REALPARM),INTENT(IN)::lms
		COMPLEX(REALPARM)::eta0,el(bkg%Nl),w(bkg%Nl,2),tmp
		INTEGER:: I,J
		DO I=1,bkg%Nl
			eta(I)=SQRT(lms*lms-bkg%k2(I))
		ENDDO
		DO I=1,bkg%Nl-1
			el(I)=EXP(-2*eta(I)*bkg%thick(I))
		ENDDO
		el(bkg%Nl)=C_ZERO
		p=Calc_p(bkg,eta,el,lms)
		q=Calc_q(bkg,eta,el)
		DO I=1,bkg%Nl
			Arr(I,:)=C_ONE/eta(I)/(C_ONE-p(I,:)*q(I,:)*el(I))
		ENDDO
		IF (PRESENT(A))THEN
			DO I=1,bkg%Nl-1
				tmp=(bkg%k2(I)-bkg%k2(I+1))
				tmp=bkg%depth(I)*tmp
				tmp=tmp/(eta(I+1)+eta(I))
				w(I,:)=(C_ONE+p(I+1,:))/(C_ONE+p(I,:)*el(I))*EXP(tmp)
			ENDDO
			DO I=1,bkg%Nl
				A(I,I,:)=Arr(I,:)
				DO J=I-1,1,-1
					A(J,I,:)=A(J+1,I,:)*w(J,:)
					A(I,J,1)=A(J,I,1)
					A(I,J,2)=A(J,I,2)*bkg%csigma(I)/bkg%csigma(J)
				ENDDO
			ENDDO
		ENDIF
	ENDSUBROUTINE
	FUNCTION Calc_p(bkg,eta,el,lms) RESULT(p)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			COMPLEX(REALPARM),INTENT(IN)::eta(bkg%Nl),el(bkg%Nl)
			REAL(REALPARM),INTENT(IN)::lms
			COMPLEX(REALPARM)::p(bkg%Nl,2)
			COMPLEX(REALPARM)::eta0,w(2),a(2),cgamma,ceta
			INTEGER:: I
			eta0=sqrt(lms*lms-bkg%k2(0))
			ceta=eta0/eta(1)
			p(1,1)=(C_ONE-ceta)/(C_ONE+ceta)
			cgamma=bkg%csigma(0)/bkg%csigma(1)

			p(1,2)=(cgamma-ceta)/(cgamma+ceta)

			DO I=2,bkg%Nl
				w(1)=el(I-1)*p(I-1,1)
				w(2)=el(I-1)*p(I-1,2)
				w=(w-1)/(w+1)
				a(1)=eta(I-1)/eta(I)
				a(2)=a(1)*bkg%csigma(I)/bkg%csigma(I-1)
				w=w*a
				p(I,:)=(1+w)/(1-w)
			ENDDO
	END FUNCTION
	FUNCTION Calc_q(bkg,eta,el) RESULT(q)
			TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
			COMPLEX(REALPARM),INTENT(IN)::eta(bkg%Nl),el(bkg%Nl)
			COMPLEX(REALPARM)::q(bkg%Nl,2)
			COMPLEX(REALPARM)::eta0,w(2),a(2)
			INTEGER:: I
			q(bkg%Nl,:)=C_ZERO
			IF (bkg%Nl>1) THEN 
				a(1)=eta(bkg%Nl)/eta(bkg%Nl-1)
				a(2)=a(1)*bkg%csigma(bkg%Nl-1)/bkg%csigma(bkg%Nl)
				q(bkg%Nl-1,:)=(1-a)/(1+a)
				DO I=bkg%Nl-2,1,-1
					w(1)=el(I+1)*q(I+1,1)
					w(2)=el(I+1)*q(I+1,2)
					w=(w-1)/(w+1)
					a(1)=eta(I+1)/eta(I)
					a(2)=a(1)*bkg%csigma(I)/bkg%csigma(I+1)
					w=w*a
					q(I,:)=(1+w)/(1-w)
				ENDDO
			ENDIF
	END FUNCTION
END MODULE
