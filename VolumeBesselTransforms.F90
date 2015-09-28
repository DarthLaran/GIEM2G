
!ATTENTION REAL(16) QUADRUPLE PRECISION ONLY!!!!!
MODULE VolumeBesselTransforms 
	USE IntegralCodes
	USE Const_module,ONLY: RealParm2,REALPARM
#ifdef  FFT_QUAD_OOURA
	USE FFT_QUAD
#endif
	IMPLICIT NONE
PRIVATE
	INTEGER ,PARAMETER::dp=REALPARM
	INTEGER ,PARAMETER::dp2=dp
	REAL(dp),PARAMETER :: PI=3.1415926535897932384626433832795028841971693993751058209_dp
	REAL(dp),PARAMETER::SPI=SQRT(PI)
	INTEGER,PARAMETER::Log2N=10
	INTEGER,PARAMETER::NI=2**(Log2N/2+1)+2
	INTEGER, PARAMETER::N=2**Log2N
	INTEGER, PARAMETER::N2=-N/2
	REAL(dp),PARAMETER::y_step=1e-1_dp*2e0_dp
	REAL(dp),PARAMETER::xi=-4.82205425726955430357775925664462079e-0002_dp*2e0_dp
	REAL(dp),PARAMETER::q=y_step/2.0_dp
	REAL(dp),PARAMETER::p=q+xi
	COMPLEX(dp2),TARGET::f1(0:N-1)
	INTEGER,PARAMETER::JX=1
	INTEGER,PARAMETER::JY=2

	INTEGER,PARAMETER::D0=1
	INTEGER,PARAMETER::D1=2
	INTEGER,PARAMETER::D2=3

	REAL(dp),TARGET::lambda_shift
	INTEGER::INDS(0:N-1)
#ifdef  FFT_QUAD_OOURA
	INTEGER,PARAMETER::FWD=1
	INTEGER,PARAMETER::BWD=-1
	REAL(dp2),TARGET::WORK(0:N/2-1)
#else
	COMPLEX(dp),TARGET::W_fwd(0:N-1),W_bwd(0:N-1)
#endif
	TYPE EXP_DATA_TYPE
		REAL(dp),DIMENSION(0:N-1)::x1,y1,dx,dy
		REAL(dp),DIMENSION(3,0:N-1,2)::xy
		REAL(dp),DIMENSION(3,0:N-1,2)::xy_exp2
		REAL(dp),DIMENSION(3,0:N-1,2)::xy_erf
		REAL(dp),DIMENSION(3,0:N-1,2)::xy2
        	REAL(dp)::explt(0:N-1)
        	REAL(dp)::explt2(0:N-1)
        	REAL(dp)::explt3(0:N-1)
	END TYPE
INTERFACE
	FUNCTION out_all(edt,I) RESULT(R)
		IMPORT EXP_DATA_TYPE,dp
		INTEGER,INTENT(IN)::I
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		REAL(dp)::R(6)
	END FUNCTION
END INTERFACE



PUBLIC ::VBTransformInit, VBTransformWeightsAllDoubleInt
PUBLIC ::VBTransformWeightsAllSingleInt


CONTAINS
#ifdef  FFT_QUAD_OOURA
SUBROUTINE VBTransformInit(s,lms) !NOT THREAD SAFETY!
	REAL(RealParm2),INTENT(IN)::s
	REAL(RealParm2),INTENT(OUT)::lms(Nfirst:Nlast)
	COMPLEX(dp)::f11(0:N-1),f21(0:N-1),tmp
	REAL(dp)::y1,theta,y2
	INTEGER::I,J,K,I2
	
	CALL makewt(N/2, INDS, WORK)

	lms=xi+(/((J+N2)*y_step,J=Nfirst,Nlast)/)
	lms=EXP(lms)/s
        lambda_shift=LOG(s)	
	DO I=0,N-1,2
		y1=(I+N2)*y_step+q
		f1(I)=inputfunction(y1)
		y1=y1+y_step
		f1(I+1)=-inputfunction(y1)

	ENDDO
	CALL FFT_16(f1,FWD)
	
END SUBROUTINE
#else
SUBROUTINE VBTransformInit(s,lms) !NOT THREAD SAFETY!
	REAL(RealParm2),INTENT(IN)::s
	REAL(RealParm2),INTENT(OUT)::lms(Nfirst:Nlast)
	COMPLEX(dp)::f11(0:N-1),f21(0:N-1),tmp
	REAL(dp)::y1,theta,y2
	INTEGER::I,J,K,I2
	INDS=bit_reverse((/(I,I=0,N-1)/));
	I2=0
	DO I=0,Log2N-1
		K=2**I
		DO J=0,K-1
			theta=-PI/K*J
			W_fwd(I2)=EXP(theta*(0e0_dp,1e0_dp))
			W_bwd(I2)=CONJG(W_fwd(I2))
			I2=I2+1
		ENDDO
	ENDDO
	lms=xi+(/((J+N2)*y_step,J=Nfirst,Nlast)/)
	lms=EXP(lms)/s
        lambda_shift=LOG(s)	
	DO I=0,N-1,2
		y1=(I+N2)*y_step+q!+lambda_shift
		f1(INDS(I))=inputfunction(y1)
		y1=y1+y_step
		f1(INDS(I+1))=-inputfunction(y1)

	ENDDO
	CALL FFT_16(f1,W_fwd)
END SUBROUTINE
#endif

SUBROUTINE CALC_EXP_LT(s,edt)
        REAL(dp),INTENT(IN)::s
	TYPE (EXP_DATA_TYPE),INTENT(INOUT)::edt
	INTEGER::I
	REAL(dp)::y2
	DO I=0,N-1,2
		y2=(I+N2)*y_step+p+s
		edt%explt(I)=EXP(y2)
		edt%explt2(I)=edt%explt(I)*edt%explt(I)
		edt%explt3(I)=edt%explt2(I)*edt%explt(I)
		y2=y2+y_step
		edt%explt(I+1)=EXP(y2)
		edt%explt2(I+1)=edt%explt(I+1)*edt%explt(I+1)
		edt%explt3(I+1)=edt%explt2(I+1)*edt%explt(I+1)
	ENDDO
END SUBROUTINE

SUBROUTINE VBTransformWeightsAllDoubleInt(x,y,hx,hy,WT)
	IMPLICIT NONE
	REAL(RealParm),INTENT(IN)::x,y,hx,hy
	REAL(RealParm),INTENT(OUT)::WT(Nfirst:Nlast,6)
		TYPE (EXP_DATA_TYPE)::edt
	REAL(dp)::rho,phi,alpha,beta
	REAL(dp)::xq,yq,hxq,hyq
	PROCEDURE(out_all),POINTER::out4
	REAL(RealParm)::rp(6)
	
	xq=x;
	yq=y;
	hxq=hx;
	hyq=hy;
	rho=SQRT(xq*xq+yq*yq);

	rp(1)=rho*rho*rho
	rp(2)=rho
	rp(3)=rho
	rp(4)=rho*rho
	rp(5)=rho
	rp(6)=rho*rho

	alpha=hxq/rho;
	beta=hyq/rho;
	phi=ATAN2(yq,xq);
        CALL CALC_EXP_LT(LOG(rho)-lambda_shift,edt)
	CALL  PrepareExps4(phi,alpha,beta,edt)
	out4=>outfunc_int4_all
	CALL CalcWeights22(edt,WT,out4)

	
	WT(:,1)=WT(:,1)/hx/hy*rp(1)
	WT(:,2)=WT(:,2)/hx/hy*rp(2) 
	WT(:,3)=WT(:,3)/hx/hy*rp(3)
	WT(:,4)=WT(:,4)/hx/hy*rp(4)
	WT(:,5)=WT(:,5)/hx/hy*rp(5)
	WT(:,6)=WT(:,6)/hx/hy*rp(6) 
END SUBROUTINE
SUBROUTINE VBTransformWeightsAllSingleInt(x,y,hx,hy,WT)
	IMPLICIT NONE
	REAL(RealParm),INTENT(IN)::x,y,hx,hy
	REAL(RealParm),INTENT(OUT)::WT(Nfirst:Nlast,6)
		TYPE (EXP_DATA_TYPE)::edt
	REAL(dp)::rho,phi,alpha,beta
	REAL(dp)::xq,yq,hxq,hyq
	PROCEDURE(out_all),POINTER::out2
	
	xq=x;
	yq=y;
	hxq=hx;
	hyq=hy;
	rho=SQRT(xq*xq+yq*yq);


	alpha=hxq/rho;
	beta=hyq/rho;
	phi=ATAN2(yq,xq);
        CALL CALC_EXP_LT(LOG(rho)-lambda_shift,edt)
	CALL  PrepareExps2(phi,alpha,beta,edt)
	out2=>outfunc_int2_all
	CALL CalcWeights22(edt,WT,out2)

	
	WT(:,1)=WT(:,1)*rho!/hx/hy
	WT(:,2)=WT(:,2)/rho!/hx/hy
	WT(:,3)=WT(:,3)/rho!/hx/hy
!	WT(:,4)=WT(:,4)/hx/hy
	WT(:,5)=WT(:,5)/rho!/hx/hy
!	WT(:,6)=WT(:,6)/hx/hy
END SUBROUTINE
#ifdef  FFT_QUAD_OOURA
SUBROUTINE CalcWeights22(edt,WT,outfunc)
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		PROCEDURE(out_all),POINTER::outfunc
		REAL(RealParm2),INTENT(OUT)::WT(Nfirst:Nlast,6)
		COMPLEX(dp2)::g11(0:N-1,6)
		COMPLEX(dp2)::h(0:N-1,6),tmp
		INTEGER::I,J
		DO I=0,N-1,2
			g11(I,:)=outfunc(edt,I)
			g11(I+1,:)=-outfunc(edt,I+1)
		ENDDO
		DO I=1,6
			CALL FFT_16(g11(:,I),FWD)
                ENDDO
		DO I=1,6
			DO J=0,N-1,2
				h(J,I)=g11(J,I)/f1(J)
				h(J+1,I)=-g11(J+1,I)/f1(J+1)
			ENDDO
		ENDDO
		DO I=1,6
			CALL FFT_16(h(:,I),BWD)
		ENDDO
		DO I=1,6
			DO J=1,N-1,2
				h(J-1,I)=h(J-1,I)/N
				h(J,I)=-h(J,I)/N
			ENDDO
		ENDDO
		
		WT=REAL(h(Nfirst:Nlast,:),KIND=REALPARM2);
	END SUBROUTINE
#else
	SUBROUTINE CalcWeights22(edt,WT,outfunc)
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		PROCEDURE(out_all),POINTER::outfunc
		REAL(RealParm2),INTENT(OUT)::WT(Nfirst:Nlast,6)
		COMPLEX(dp)::g11(6,0:N-1)
		COMPLEX(dp)::h(6,0:N-1),tmp
		INTEGER::I,J
		g11=1.0_dp
		DO I=0,N-1,2
			g11(:,INDS(I))=outfunc(edt,I)
			g11(:,INDS(I+1))=-outfunc(edt,I+1)
		ENDDO
		CALL FFT_16_6(g11,W_fwd)
		DO J=0,N-1,2
			h(:,INDS(J))=g11(:,J)/f1(J)
			h(:,INDS(J+1))=-g11(:,J+1)/f1(J+1)
		ENDDO
		CALL FFT_16_6(h,W_bwd)
		DO J=1,N-1,2
				h(:,J-1)=h(:,J-1)/N
				h(:,J)=-h(:,J)/N
		ENDDO
		
		WT=TRANSPOSE(REAL(h(:,Nfirst:Nlast),KIND=REALPARM2));
	END SUBROUTINE
#endif
!------------------------ Input Function  -------------------------------------!

FUNCTION inputfunction(t) RESULT (R)
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
	REAL(dp)::y,y2,y4
	y=EXP(-t);
		y2=y*y;
	y4=y2*y2
		R=8.0_dp*EXP(-y2)*(y4-4.0_dp*y2+2.0_dp)*y;
END FUNCTION

!--------------------------------------------------------------------------------!

SUBROUTINE PrepareExps4(phi,alpha,beta,edt)
		REAL(dp),INTENT(IN)::phi,alpha,beta
		TYPE (EXP_DATA_TYPE),INTENT(INOUT)::edt
		REAL(dp)::cphi,sphi,x(3),y(3),x2(3),y2(3)
		INTEGER::I

		cphi=COS(phi)
		sphi=SIN(phi)
		edt%x1=edt%explt*cphi
		edt%y1=edt%explt*sphi

		edt%dx=alpha*edt%explt
		edt%dy=beta*edt%explt

		edt%x1=edt%x1-edt%dx*0.5_dp
		edt%y1=edt%y1-edt%dy*0.5_dp
		DO I=0,N-1
			x(1)=edt%x1(I)
			x(2)=x(1)+edt%dx(I)
			x(3)=x(1)-edt%dx(I)
			x=x*0.5_dp
			x2=x*x
			edt%xy_exp2(:,I,JX)=EXP(-x2)
			edt%xy_erf(:,I,JX)=ERF(x)
			edt%xy(:,I,JX)=x
			edt%xy2(:,I,JX)=x2
			y(1)=edt%y1(I)
			y(2)=y(1)+edt%dy(I)
			y(3)=y(1)-edt%dy(I)
			y=y*0.5_dp
			y2=y*y
			edt%xy_exp2(:,I,JY)=EXP(-y2)
			edt%xy_erf(:,I,JY)=ERF(y)
			edt%xy(:,I,JY)=y
			edt%xy2(:,I,JY)=y2
		ENDDO
ENDSUBROUTINE

SUBROUTINE PrepareExps2(phi,alpha,beta,edt)
		REAL(dp),INTENT(IN)::phi,alpha,beta
		TYPE (EXP_DATA_TYPE),INTENT(INOUT)::edt
		REAL(dp)::cphi,sphi,x(2),y(2),x2(2),y2(2)
		INTEGER::I

		cphi=COS(phi)
		sphi=SIN(phi)
		edt%x1=edt%explt*cphi*0.5_dp
		edt%y1=edt%explt*sphi*0.5_dp
		edt%dx=alpha*edt%explt*0.5_dp
		edt%dy=beta*edt%explt*0.5_dp

		DO I=0,N-1
			x(1)=edt%x1(I)
			x(2)=x(1)+edt%dx(I)
			x2=x*x
			edt%xy_exp2(1:2,I,JX)=EXP(-x2)
			edt%xy_erf(1:2,I,JX)=ERF(x)
			edt%xy(1:2,I,JX)=x
			edt%xy2(1:2,I,JX)=x2
			y(1)=edt%y1(I)
			y(2)=y(1)+edt%dy(I)
			y2=y*y
			edt%xy_exp2(1:2,I,JY)=EXP(-y2)
			edt%xy_erf(1:2,I,JY)=ERF(y)
			edt%xy(1:2,I,JY)=y
			edt%xy2(1:2,I,JY)=y2
		ENDDO
ENDSUBROUTINE
!------------------------ Output Functions	-------------------------------------!
!----------------- INT4 ---------------------------------------------------------!
!lt is distance between bottom left corner of the first rectangle and center of the second. The sizes of rectangles are the same

!------------------------------ ---------------------------------------!
FUNCTION outfunc_int4_all(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R(6)
		REAL(dp)::fx(1:3,3),fy(1:3,3)
		fx(:,D0)=DOUBLE_INT_D0(edt,I,JX)
		fx(:,D1)=DOUBLE_INT_D1(edt,I,JX)
		fx(:,D2)=DOUBLE_INT_D2(edt,I,JX)

		fy(:,D0)=DOUBLE_INT_D0(edt,I,JY)
		fy(:,D1)=DOUBLE_INT_D1(edt,I,JY)
		fy(:,D2)=DOUBLE_INT_D2(edt,I,JY)

		R(1)=fx(3,D0)*fy(1,D0)+fx(1,D0)*fy(3,D0)+2.0_dp*fx(2,D0)*fy(2,D0)
		R(2)=fx(3,D2)*fy(1,D0)+fx(1,D2)*fy(3,D0)+2.0_dp*fx(2,D2)*fy(2,D0)
		R(3)=fx(3,D1)*fy(1,D1)+fx(1,D1)*fy(3,D1)+2.0_dp*fx(2,D1)*fy(2,D1)
		R(4)=fx(3,D1)*fy(1,D0)+fx(1,D1)*fy(3,D0)+2.0_dp*fx(2,D1)*fy(2,D0)
		R(5)=fx(3,D0)*fy(1,D2)+fx(1,D0)*fy(3,D2)+2.0_dp*fx(2,D0)*fy(2,D2)
		R(6)=fx(3,D0)*fy(1,D1)+fx(1,D0)*fy(3,D1)+2.0_dp*fx(2,D0)*fy(2,D1)

		R(1)=R(1)/edt%explt3(I)
		R(2)=R(2)/edt%explt(I)
		R(3)=R(3)/edt%explt(I)
		R(4)=R(4)/edt%explt2(I)
		R(5)=R(5)/edt%explt(I)
		R(6)=R(6)/edt%explt2(I)
END FUNCTION
!-------------------------------- INT2 ------------------------------------!
!lt is distance between bottom left corner of the  rectangle and point.

FUNCTION outfunc_int2_all(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R(6)
		REAL(dp)::fx(1:3,3),fy(1:3,3)
		REAL(dp)::t(4),q(4)
		fx(:,D0)=SINGLE_INT_D0(edt,I,JX)
		fx(:,D1)=SINGLE_INT_D1(edt,I,JX)
		fx(:,D2)=SINGLE_INT_D2(edt,I,JX)

		fy(:,D0)=SINGLE_INT_D0(edt,I,JY)
		fy(:,D1)=SINGLE_INT_D1(edt,I,JY)
		fy(:,D2)=SINGLE_INT_D2(edt,I,JY)

		R(1)=fx(3,D0)*fy(1,D0)+fx(1,D0)*fy(3,D0)+2.0_dp*fx(2,D0)*fy(2,D0)
		R(2)=fx(3,D2)*fy(1,D0)+fx(1,D2)*fy(3,D0)+2.0_dp*fx(2,D2)*fy(2,D0)

		R(3)=fx(3,D1)*fy(1,D1)+fx(1,D1)*fy(3,D1)+2.0_dp*fx(2,D1)*fy(2,D1)

		R(4)=fx(3,D1)*fy(1,D0)+fx(1,D1)*fy(3,D0)+2.0_dp*fx(2,D1)*fy(2,D0)
		R(5)=fx(3,D0)*fy(1,D2)+fx(1,D0)*fy(3,D2)+2.0_dp*fx(2,D0)*fy(2,D2)
		R(6)=fx(3,D0)*fy(1,D1)+fx(1,D0)*fy(3,D1)+2.0_dp*fx(2,D0)*fy(2,D1)

		R(1)=R(1)/edt%explt(I)
		R(2)=R(2)*edt%explt(I)
		R(3)=R(3)*edt%explt(I)
		R(5)=R(5)*edt%explt(I)
END FUNCTION
!-------------------------------------------------------------------------------------!

FUNCTION DOUBLE_INT_D2(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3),x(1:3),x2(1:3),ex(1:3)
	x=edt%xy(:,I,J)
	x2=edt%xy2(:,I,J)
	ex=edt%xy_exp2(:,I,J)/2.0_dp
	f=x2*ex*4.0_dp
	R(2)=f(2)+f(3)-2.0_dp*f(1)
	f=x2*f*4.0_dp
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	f=ex
	R(3)=f(2)+f(3)-2.0_dp*f(1)

END FUNCTION
FUNCTION INT_D2TEXPd4_2d_edt(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3),x(1:3),x2(1:3),ex(1:3)
	x=edt%xy(:,I,J)
	x2=edt%xy2(:,I,J)
	ex=edt%xy_exp2(:,I,J)/2.0_dp
	f=x2*ex*4.0_dp
	R(2)=f(2)+f(3)-2.0_dp*f(1)
	f=x2*f*4.0_dp
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	f=ex
	R(3)=f(2)+f(3)-2.0_dp*f(1)

END FUNCTION

FUNCTION DOUBLE_INT_D1(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3)
	f=INT_1d_edt(edt,I,J)
	R(1)=-8.0_dp*f(3)-12.0_dp*f(2)+12.0_dp*f(1);
	R(2)=-2.0_dp*f(2)+2.0_dp*f(1)
		R(3)=f(1)
!R=R*0.5_dp
END FUNCTION

FUNCTION INT_D1TEXPd4_2d_edt(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3)
	f=INT_1d_edt(edt,I,J)
	R(1)=-8.0_dp*f(3)-12.0_dp*f(2)+12.0_dp*f(1);
	R(2)=-2.0_dp*f(2)+2.0_dp*f(1)
		R(3)=f(1)
!R=R*0.5_dp
END FUNCTION
FUNCTION INT_1d_edt(edt,I,J)RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::x(1:3),f(1:3),x2(1:3)
	x=edt%xy(:,I,J)
	x2=edt%xy2(:,I,J)
	f=edt%xy_erf(:,I,J)
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	R(1)=R(1)*SPI*0.5_dp

	f=edt%xy_exp2(:,I,J)*x
	R(2)=f(2)+f(3)-2.0_dp*f(1)
	f=f*x2
	R(3)=f(2)+f(3)-2.0_dp*f(1)
END FUNCTION
FUNCTION DOUBLE_INT_D0(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3)
	f=INT_2d_edt(edt,I,J)
	R(1)=32.0_dp*(1.5_dp*f(1)+f(2)+f(3))
	R(2)=8.0_dp*f(1)+4.0_dp*f(2)
	R(3)=4.0_dp*f(1)+f(2);
END FUNCTION


FUNCTION INT_TEXPd4_2d_edt(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3)
	f=INT_2d_edt(edt,I,J)
	R(1)=32.0_dp*(1.5_dp*f(1)+f(2)+f(3))
	R(2)=8.0_dp*f(1)+4.0_dp*f(2)
	R(3)=4.0_dp*f(1)+f(2);
END FUNCTION


FUNCTION INT_2d_edt(edt,I,J)RESULT(R)

	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(1:3)
	REAL(dp)::x(1:3),f(1:3),x2(1:3)
	x=edt%xy(:,I,J)
	x2=edt%xy2(:,I,J)
	f=x*edt%xy_erf(:,I,J)
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	R(1)=R(1)*SPI*0.25_dp

	f=edt%xy_exp2(:,I,J)
	R(2)=f(2)+f(3)-2.0_dp*f(1)!h should be lagrger than 1e-5l
   !	 R(2)=EXP(-l*l)*(-2.0_dp+4.0*l*l)*h*h
	f=edt%xy_exp2(:,I,J)*x2/4.0_dp
	R(3)=f(2)+f(3)-2.0_dp*f(1)
END FUNCTION


FUNCTION SINGLE_INT_D0(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(3)
	REAL(dp)::f(2),a,b,c
	f=edt%xy_erf(1:2,I,J)*SPI
	a=f(2)-f(1)
	R(1)=(a)*0.5_dp
	f=edt%xy(1:2,I,J)*edt%xy_exp2(1:2,I,J)
	b=f(2)-f(1)
	R(2)=a-2.0_dp*b
	f=f*edt%xy2(1:2,I,J)
	c=f(2)-f(1)
	R(3)=6.0_dp*a-12.0_dp*b-8.0_dp*c
END FUNCTION

FUNCTION SINGLE_INT_D1(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(3)
	REAL(dp)::f(2)

	f=edt%xy_exp2(1:2,I,J)
	R(1)=f(2)-f(1)
	f=f*edt%xy2(1:2,I,J)
	R(2)=f(2)-f(1)

	f=f*edt%xy2(1:2,I,J)
	R(3)=f(2)-f(1)
	R(1)=R(1)/2.0_dp
	R(2)=R(2)*2.0_dp
	R(3)=R(3)*8.0_dp
END FUNCTION

FUNCTION SINGLE_INT_D2(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(3)
	REAL(dp)::f(2),a,b,c
	f=edt%xy(1:2,I,J)*edt%xy_exp2(1:2,I,J)
	a=f(2)-f(1)
	R(1)=-a*0.5_dp
	f=f*edt%xy2(1:2,I,J)
	b=f(2)-f(1)
	R(2)=2.0_dp*(a-b)
	f=f*edt%xy2(1:2,I,J)
	c=f(2)-f(1)
	R(3)=16.0_dp*b-8.0_dp*c


END FUNCTION


!-------------------------------------------------------------!

#ifdef  FFT_QUAD_OOURA
	SUBROUTINE FFT_16(X,dir)!x MUST be after bit reverse, wp -array of coefficients for forward or backward FT
		USE ISO_C_BINDING
		COMPLEX(dp2),TARGET ,DIMENSION(0:N-1), INTENT(inout) :: x
		INTEGER,INTENT(IN)::dir
		REAL(dp2),POINTER ::pX(:)
		TYPE(C_PTR)::cp
		cp=C_LOC(x)
		CALL C_F_POINTER(cp,pX,(/2*N/))
		CALL cdft(N*2,dir, pX, INDS, WORK)
	END SUBROUTINE 
#else
	ELEMENTAL  FUNCTION bit_reverse(I) RESULT(R)
		INTEGER, INTENT(in) :: I
		INTEGER ::R 
		INTEGER :: J, temp
		temp = I
		DO J = Log2N, 2, -1
			   temp = ISHFTC(temp, -1,J)
		END DO
		R = temp
	END FUNCTION bit_reverse
	SUBROUTINE FFT_16(x,wp)!x MUST be after bit reverse, wp -array of coefficients for forward or backward FT
		COMPLEX(dp), DIMENSION(0:N-1), INTENT(inout) :: x
		COMPLEX(dp),INTENT(IN)::wp(0:N-1)
		COMPLEX(dp) ::	temp
		INTEGER :: I,J,J2,I2,Istep,Ip,Iq,Lstep
			I2=0
		DO I=0,Log2N-1
			J2=2**I
			Lstep = 2 * J2
			DO J = 0, J2-1
				DO Ip=J, N-1, Lstep
					Iq=Ip+J2
					temp = wp(I2)*x(Iq)
					x(Iq) = x(Ip) - temp
					x(Ip) = x(Ip) + temp
				END DO
				I2=I2+1
			END DO
		END DO
	END SUBROUTINE
	SUBROUTINE FFT_16_6(x,wp)!x MUST be after bit reverse, wp -array of coefficients for forward or backward FT
		COMPLEX(dp), DIMENSION(6,0:N-1), INTENT(inout) :: x
		COMPLEX(dp),INTENT(IN)::wp(0:N-1)
		COMPLEX(dp) ::	temp(6)
		INTEGER :: I,J,J2,I2,Istep,Ip,Iq,Lstep
		I2=0
		J2=1
		DO I=0,Log2N-1
			Lstep = 2 * J2
			DO J = 0, J2-1
				DO Ip=J, N-1, Lstep
					Iq=Ip+J2
					temp = wp(I2)*x(:,Iq)
					x(:,Iq) = x(:,Ip) - temp
					x(:,Ip) = x(:,Ip) + temp
				END DO
				I2=I2+1
			END DO
			J2=Lstep
		END DO
	END SUBROUTINE FFT_16_6
#endif
END
