MODULE SOURCES_MODULE
	USE CONST_MODULE
	USE APQ_MODULE
        USE LOGGER_MODULE
	USE DATA_TYPES_MODULE
	IMPLICIT NONE
CONTAINS
	SUBROUTINE PlaneWave(pol,bkg,zrecv,F0)  
		INTEGER,INTENT(IN)::pol
	        TYPE (BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
		REAL(REALPARM),INTENT(IN)::zrecv(:)
		COMPLEX(REALPARM),INTENT(OUT)::F0(:,:,:,:)
		COMPLEX(REALPARM)::G(2)
		COMPLEX(REALPARM)::A(bkg%Nl,bkg%Nl,2),p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl),Arr(bkg%Nl,2)
		COMPLEX(REALPARM)::t1,t2,e1,e2
		INTEGER::Nr,Ir,l
		F0=C_ZERO
		CALL Calc_APQ(bkg,R_ZERO,Arr,p,q,eta,A)
		IF (bkg%Nl>1) THEN
			t1=C_ONE+q(1,1)*EXP(-bkg%depth(1)*eta(1)*C_TWO)
		ELSE
			t1=C_ONE
		ENDIF
		t1=A(1,1,1)*t1
		Nr=SIZE(zrecv)
		DO Ir=1,Nr
			l=GetLayer(zrecv(Ir),bkg)
			l=MAX(l,1)
			t2=A(l,1,1)/t1
			e1=EXP(-eta(l)*zrecv(Ir))
			IF (l<bkg%Nl) THEN
				e2=q(l,1)*(EXP(-eta(l)*(C_TWO*bkg%depth(l)-zrecv(Ir))))
			ELSE
				e2=C_ZERO
			ENDIF
			G(1)=t2*(e1+e2)
			G(2)=-t2*(e1-e2)*eta(l)/bkg%iwm
			IF (pol==EX) THEN
				F0(Ir,:,:,EX)=G(1)
				F0(Ir,:,:,HY)=G(2)
			ELSEIF (pol==EY) THEN
				F0(Ir,:,:,EY)=G(1)
				F0(Ir,:,:,HX)=-G(2)
			ELSE
				CALL LOGGER('INCORRECT POLARIZATION')
			ENDIF
		ENDDO
	END SUBROUTINE
	SUBROUTINE PlaneWaveIntegral(pol,bkg,anomaly,E0)
		INTEGER,INTENT(IN)::pol
        TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
		COMPLEX(REALPARM),INTENT(OUT)::E0(:,:,:,:)
		COMPLEX(REALPARM)::A(bkg%Nl,bkg%Nl,2),p(bkg%Nl,2),q(bkg%Nl,2),eta(bkg%Nl),Arr(bkg%Nl,2)
		COMPLEX(REALPARM)::t1,t2,e1,e2,e3,e4,G
		INTEGER::Nz,Iz,l
		REAL(REALPARM)::dz
		E0=0d0
		CALL Calc_APQ(bkg,R_ZERO,Arr,p,q,eta,A)
		IF (bkg%Nl>1) THEN
			t1=C_ONE+q(1,1)*EXP(-bkg%depth(1)*eta(1)*C_TWO)
		ELSE
			t1=C_ONE
		ENDIF
		t1=A(1,1,1)*t1
		Nz=anomaly%Nz
		DO Iz=1,Nz
			dz=anomaly%z(Iz)-anomaly%z(Iz-1)
			l=anomaly%Lnumber(Iz)
			e1=EXP(-eta(l)*anomaly%z(Iz-1))
			e2=EXP(-eta(l)*anomaly%z(Iz))
			IF (l<bkg%Nl) THEN
				e3=q(l,1)*(EXP(-eta(l)*(C_TWO*bkg%depth(l)-anomaly%z(Iz-1))))
				e4=q(l,1)*(EXP(-eta(l)*(C_TWO*bkg%depth(l)-anomaly%z(Iz))))
			ELSE
				e3=C_ZERO
				e4=C_ZERO
			ENDIF
			G=(e1-e2+e4-e3)*A(l,1,1)/t1
			G=G/dz/eta(l)
			E0(Iz,pol,:,:)=G
		ENDDO
	END SUBROUTINE
END MODULE
