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

MODULE Impedances_Module
	 	USE CONST_MODULE
		USE DATA_TYPES_MODULE
		IMPLICIT NONE	
	
       	REAL (REALPARM),POINTER:: sig(:),lt(:),h(:)
      	COMPLEX (RealParm),POINTER :: k2(:),k(:)
        COMPLEX(REALPARM)::iwm
        REAL(REALPARM),POINTER::z(:)
        INTEGER,POINTER::Ln(:)
        INTEGER::Nz,Nh
        INTERFACE CalcImp_t
	        MODULE PROCEDURE   CalcImp_t_exp,CalcImp_t_len
        END INTERFACE


        INTERFACE CalcImp_b
	        MODULE PROCEDURE   CalcImp_b_exp, CalcImp_b_len
        END INTERFACE


   CONTAINS
	SUBROUTINE Prepare_Hankel(bkg,anomaly)
		TYPE(BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
		TYPE (ANOMALY_TYPE),TARGET,OPTIONAL,INTENT(IN)::anomaly
		h=>bkg%depth
		sig=>bkg%sigma
		lt=>bkg%thick
		Nh=bkg%Nl
		k=>bkg%k
		k2=>bkg%k2
		iwm=bkg%iwm
		IF (PRESENT(anomaly)) THEN
			Nz=anomaly%Nz
			z=>anomaly%z
			Ln=>anomaly%LNumber
		ENDIF
	END SUBROUTINE
	PURE FUNCTION CalcImp_t_len(Z0,eta,x) RESULT(f)
		IMPLICIT NONE
		COMPLEX(KIND = RealParm2),INTENT(IN):: Z0,eta
		COMPLEX(KIND = RealParm2):: f
		REAL (KIND = RealParm),INTENT(IN):: x
		COMPLEX(KIND = RealParm2):: y,Z,Z1
		Z1=Z0/eta
		Z1=(1e0_REALPARM2-Z1)/(1e0_REALPARM2+Z1)
		y=EXP(-2e0_REALPARM2*eta*x);
		Z=-eta*(Z1*y-1e0_REALPARM2)/(Z1*y+1e0_REALPARM2);
		f=Z
	END FUNCTION CalcImp_t_len

	PURE FUNCTION CalcImp_t_exp(Z0,eta,e) RESULT(f)
		IMPLICIT NONE
		COMPLEX(KIND = RealParm2),INTENT(IN):: Z0,eta,e
		COMPLEX(KIND = RealParm2):: f
		COMPLEX(KIND = RealParm2):: y,Z,Z1
		Z1=Z0/eta
		Z1=(1e0_REALPARM2-Z1)/(1e0_REALPARM2+Z1)
		y=e*e!EXP(-2e0_REALPARM2*eta*x);
		Z=-eta*(Z1*y-1e0_REALPARM2)/(Z1*y+1e0_REALPARM2);
		f=Z
	END FUNCTION CalcImp_t_exp
	PURE FUNCTION CalcImp_b_len(Z0,eta,x) RESULT(f)
		IMPLICIT NONE
		COMPLEX(KIND = RealParm2),INTENT(IN):: Z0,eta
		COMPLEX(KIND = RealParm2):: f
		REAL (KIND = RealParm),INTENT(IN):: x
		COMPLEX(KIND = RealParm2):: y,Z,Z1

		Z1=Z0/eta
		Z1=(1e0_REALPARM2+Z1)/(1e0_REALPARM2-Z1)
		y=exp(-2e0_REALPARM2*eta*x);
		Z=eta*(Z1*y-1e0_REALPARM2)/(Z1*y+1e0_REALPARM2);
		f=Z
	END FUNCTION CalcImp_b_len

	PURE FUNCTION CalcImp_b_exp(Z0,eta,e) RESULT(f)
		IMPLICIT NONE
		COMPLEX(KIND = RealParm2),INTENT(IN):: Z0,eta,e
		COMPLEX(KIND = RealParm2):: f
		COMPLEX(KIND = RealParm2):: y,Z,Z1

		Z1=Z0/eta
		Z1=(1e0_REALPARM2+Z1)/(1e0_REALPARM2-Z1)
		y=e*e!exp(-2e0_REALPARM2*eta*x);
		Z=eta*(Z1*y-1e0_REALPARM2)/(Z1*y+1e0_REALPARM2);
		f=Z
	END FUNCTION CalcImp_b_exp
	FUNCTION Calc_Depth_Imp(lms,z0,lnum) RESULT(Z)
		IMPLICIT NONE
	COMPLEX(KIND = RealParm2):: Z
	REAL (KIND = RealParm2),INTENT(IN):: lms
	REAL (KIND = RealParm),INTENT(IN):: z0
	INTEGER,INTENT(IN)::lnum
	COMPLEX(KIND = RealParm2):: eta
	INTEGER::Ih
	Z=-SQRT(lms*lms-k2(Nh))
	IF (lnum.EQ.(Nh)) THEN
		return
	ELSE
		DO Ih=Nh-1,lnum+1,-1
			eta=SQRT(lms*lms-k2(Ih))
			Z=CalcImp_b( Z,eta,lt(Ih))
!			PRINT*,"QQQQ",Ih
		END DO
		IF (lnum/=0) THEN
			eta=SQRT(lms*lms-k2(lnum))
			Z=CalcImp_b( Z,eta,h(lnum)-z0)
		ENDIF
!		PRINT*,"PPP",Nh-1,lnum+1
		
	ENDIF
	END FUNCTION Calc_Depth_Imp
    
	PURE FUNCTION Calc_Depth_Imp2(lms,z0,lnum) RESULT(Z)
		IMPLICIT NONE
	COMPLEX(KIND = RealParm2):: Z
	REAL (KIND = RealParm2),INTENT(IN):: lms
	REAL (KIND = RealParm),INTENT(IN):: z0
	INTEGER,INTENT(IN)::lnum
	COMPLEX(KIND = RealParm2):: eta
	INTEGER:: Ih
	IF (lnum.EQ.(Nh)) THEN
		Z=-SQRT(lms*lms-k2(Nh))
		return
	ELSE
		Z=-SQRT(lms*lms-k2(Nh))/sig(Nh)*sig(Nh-1)
		DO Ih=Nh-1,lnum+1,-1
			eta=SQRT(lms*lms-k2(Ih))
			Z=CalcImp_b( Z,eta,lt(Ih))
			Z=Z/sig(Ih)*sig(Ih-1)
		END DO
		IF (lnum/=0) THEN
			eta=SQRT(lms*lms-k2(lnum))
			Z=CalcImp_b( Z,eta,h(lnum)-z0)
		ENDIF
	ENDIF

	END FUNCTION Calc_Depth_Imp2
    
	PURE FUNCTION Calc_Top_Imp(lms,z0,lnum) RESULT(Z)
		IMPLICIT NONE
	COMPLEX(KIND = RealParm2):: Z
	REAL (KIND = RealParm2),INTENT(IN):: lms
	REAL (KIND = RealParm),INTENT(IN):: z0
	INTEGER,INTENT(IN)::lnum
	COMPLEX(KIND = RealParm2):: eta
		INTEGER::Ih 
	!PRINT*,k2,nh+1
	
	Z=lms
		IF (lnum==0) RETURN
	IF (lnum==1) THEN
		eta=sqrt(lms*lms-k2(1));
		Z=CalcImp_t(Z,eta,z0);
		RETURN
	ELSE
		DO Ih=1,lnum-1,1
			eta=SQRT(lms*lms-k2(Ih))
			Z=CalcImp_t( Z,eta,lt(Ih))
		END DO
		eta=SQRT(lms*lms-k2(lnum))
		Z=CalcImp_t( Z,eta,z0-h(lnum-1))
	ENDIF
	END FUNCTION Calc_Top_Imp

	FUNCTION Calc_Top_Imp2(lms,z0,lnum) RESULT(Z)
	IMPLICIT NONE
    !!ATTENTION!! 
    !IF z0 is not at air-ground the top impedance for  z0 returned
    !	other case -undefined
	COMPLEX(KIND = RealParm2):: Z
	REAL (KIND = RealParm2),INTENT(IN):: lms
	REAL (KIND = RealParm),INTENT(IN):: z0
	INTEGER,INTENT(IN)::lnum
	COMPLEX(KIND = RealParm2):: eta,y
		INTEGER::Ih
		eta=sqrt(lms*lms-k2(1));
		IF (lnum.EQ.(1)) THEN
			IF (z0>0D0) THEN
				y=EXP(-2e0_REALPARM2*sqrt(lms*lms-k2(1))*z0)
				Z=eta*(1e0_REALPARM2+y)/(1e0_REALPARM2-y)
				RETURN
			ELSE	
!			IF (Ln(1).EQ.1) THEN
!				y=EXP(-2e0_REALPARM2*sqrt(lms**2e0_REALPARM2*z(1)**2e0_REALPARM2-k2(1)*z(1)**2e0_REALPARM2))
!				Z=eta*(1e0_REALPARM2+y)/(1e0_REALPARM2-y)
!				RETURN
!			ELSE
!				y=EXP(-2e0_REALPARM2*sqrt(lms**2e0_REALPARM2*h(1)**2e0_REALPARM2-k2(1)*h(1)**2e0_REALPARM2))
!				Z=eta*(1e0_REALPARM2+y)/(1e0_REALPARM2-y)*sig(2)/sig(1)
!				RETURN
!			ENDIF
		!	Z=NaN
		ENDIF
	ELSE
		y=EXP(-2e0_REALPARM2*sqrt(lms*lms-k2(1))*h(1))
		Z=eta*(1e0_REALPARM2+y)/(1e0_REALPARM2-y)/sig(1)*sig(2)
		DO Ih=2,lnum-1,1
			eta=SQRT(lms*lms-k2(Ih))
			Z=CalcImp_t( Z,eta,lt(Ih))
			Z=Z/sig(Ih)*sig(Ih+1)
		END DO
		eta=SQRT(lms*lms-k2(lnum))
		Z=CalcImp_t( Z,eta,z0-h(lnum-1))
		!PRINT*,'src Z',Z,z0,lnum
	ENDIF

	END FUNCTION Calc_Top_Imp2

END MODULE
