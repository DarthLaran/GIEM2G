
!ATTENTION REAL(16) QUADRUPLE PRECISION ONLY!!!!!
MODULE VolumeBesselTransforms 
	USE IntegralCodes
	USE Const_module,ONLY: RealParm2,REALPARM
PRIVATE
	INTEGER ,PARAMETER::dp=16
		REAL (dp),PARAMETER :: PI=3.1415926535897932384626433832795028841971693993751058209_dp
	REAL(dp),PARAMETER::SPI=SQRT(PI)
	INTEGER,PARAMETER::Log2N=11
	INTEGER, PARAMETER::N=2**Log2N
	INTEGER, PARAMETER::N2=-N/2
	REAL(dp),PARAMETER::y_step=1e-1_dp!*2e0_dp
	REAL(dp),PARAMETER::xi=-4.82205425726955430357775925664462079e-0002_dp!*2e0_dp
	REAL(dp),PARAMETER::q=y_step/2.0_dp
	REAL(dp),PARAMETER::p=q+xi
		COMPLEX(dp),TARGET::f1(0:N-1)
	INTEGER::INDS(0:N-1)
	COMPLEX(dp),TARGET::W_fwd(0:N-1),W_bwd(0:N-1)
	REAL(dp),TARGET::explt(0:N-1)
	REAL(dp),TARGET::explt2(0:N-1)
	REAL(dp),TARGET::explt3(0:N-1)
	INTEGER,PARAMETER::JX=1
	INTEGER,PARAMETER::JY=2

	INTEGER,PARAMETER::D0=1
	INTEGER,PARAMETER::D1=2
	INTEGER,PARAMETER::D2=3
TYPE EXP_DATA_TYPE
	REAL(dp),DIMENSION(0:N-1)::x1,y1,dx,dy
	REAL(dp),DIMENSION(3,0:N-1,2)::xy
	REAL(dp),DIMENSION(3,0:N-1,2)::xy_exp2
	REAL(dp),DIMENSION(3,0:N-1,2)::xy_erf
	REAL(dp),DIMENSION(3,0:N-1,2)::xy2
END TYPE
INTERFACE
	FUNCTION outfunction(lt,cphi,sphi,dxr,dyr) RESULT(R)
		IMPORT dp
		REAL(dp),INTENT(IN)::lt,cphi,sphi,dxr,dyr
		REAL(dp)::R
	END FUNCTION
	FUNCTION outfunction2(edt,I) RESULT(R)
		IMPORT EXP_DATA_TYPE,dp
		INTEGER,INTENT(IN)::I
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		REAL(dp)::R
	END FUNCTION
	FUNCTION out_all(edt,I) RESULT(R)
		IMPORT EXP_DATA_TYPE,dp
		INTEGER,INTENT(IN)::I
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		REAL(dp)::R(6)
	END FUNCTION
END INTERFACE



PUBLIC ::VBTransformWeights,VBTransformInit, VBTransformWeightsAllInt4,VBTransformWeightsAllInt42
PUBLIC ::VBTransformWeightsAllInt22


CONTAINS

SUBROUTINE VBTransformInit(lms) !NOT THREAD SAFETY!
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
	lms=EXP(lms)
	DO I=0,N-1,2
		y1=(I+N2)*y_step+q
		f1(INDS(I))=inputfunction(y1)
		y1=y1+y_step
		f1(INDS(I+1))=-inputfunction(y1)

		y2=(I+N2)*y_step+p
		explt(I)=EXP(y2)
		explt2(I)=explt(I)*explt(I)
		explt3(I)=explt2(I)*explt(I)
		y2=y2+y_step
		explt(I+1)=EXP(y2)
		explt2(I+1)=explt(I+1)*explt(I+1)
		explt3(I+1)=explt2(I+1)*explt(I+1)
	ENDDO
	CALL FFT_16(f1,W_fwd)
	
END SUBROUTINE
SUBROUTINE VBTransformWeights(x,y,hx,hy,WT,c,IERROR)
	IMPLICIT NONE
	REAL(RealParm),INTENT(IN)::x,y,hx,hy
	INTEGER,INTENT(IN)::c
	INTEGER,OPTIONAL,INTENT(OUT)::IERROR
	REAL(RealParm2),INTENT(OUT)::WT(Nfirst:Nlast)
	REAL(dp)::rho,phi,alpha,beta
	REAL(dp)::xq,yq,hxq,hyq
		COMPLEX(dp),POINTER::inputfunc(:)
	REAL(RealParm2)::rp
	PROCEDURE(outfunction),POINTER::outfunc
	
	xq=x;
	yq=y;
	hxq=hx;
	hyq=hy;
	rho=SQRT(xq*xq+yq*yq);
	rp=1e0_dp
		SELECT CASE (c)
				 CASE (INT4DYY:INT4DY,INT2DYY:INT2DY)
					alpha=hyq/rho;
					beta=hxq/rho;
					phi=ATAN2(xq,yq);
				CASE DEFAULT
					alpha=hxq/rho;
					beta=hyq/rho;
					phi=ATAN2(yq,xq);
		END SELECT

		SELECT CASE (c)
!--------------- INT4 ----------------!
				CASE (INT4)
						outfunc=>outfunci4cl
						inputfunc=>f1
						rp=rho*rho*rho
				CASE (INT4DXX)
						outfunc=>outfunci4dxdxcl
						inputfunc=>f1
						rp=rho
				CASE (INT4DXY)
						outfunc=>outfunci4dxdycl
						inputfunc=>f1
						rp=rho
				CASE (INT4DX)
						outfunc=>outfunci4dxcl
						inputfunc=>f1
						rp=rho*rho
				CASE  (INT4DYY)
						outfunc=>outfunci4dxdxcl
						inputfunc=>f1
						rp=rho
				CASE (INT4DY)
						outfunc=>outfunci4dxcl
						inputfunc=>f1
						rp=rho*rho
!--------------- INT2 ----------------!
				CASE (INT2)
						outfunc=>outfunci2cl
						inputfunc=>f1
			rp=rho
				CASE (INT2DXX)
						outfunc=>outfunci2dxdxcl
						inputfunc=>f1
			rp=1e0_dp/rho
				CASE (INT2DXY)
						outfunc=>outfunci2dxdycl
						inputfunc=>f1
			rp=1e0_dp/rho
				CASE (INT2DX)
						outfunc=>outfunci2dxcl
						inputfunc=>f1
				CASE (28)!(INT2DYY)
						outfunc=>outfunci2dxdxcl
						inputfunc=>f1
			rp=1e0_dp/rho
				CASE (30)!(INT2DY)
						outfunc=>outfunci2dxcl
						inputfunc=>f1
			   CASE DEFAULT
				PRINT*, 'ERROR! The wrong type', c ,'. May be in future? '
			IF (PRESENT(IERROR)) THEN
				IF (c/=0) THEN 
					IERROR=c
					PRINT*,x,y,hx,hy
				ELSE
					IERROR=-1
					PRINT*,x,y,hx,hy
				ENDIF
			ENDIF
					return
		END SELECT
	CALL CalcWeights(phi,alpha,beta,WT,inputfunc,outfunc);
	
		WT=WT/hx/hy*rp 
	IF (PRESENT(IERROR)) THEN
		IERROR=0
	ENDIF
END SUBROUTINE

	SUBROUTINE CalcWeights(phi,alpha,beta,WT,inputfunc,outfunc)
		REAL(dp),INTENT(IN)::phi,alpha,beta
		COMPLEX(dp),INTENT(IN)::inputfunc(0:N-1)
		PROCEDURE(outfunction),POINTER,INTENT(IN)::outfunc
		REAL(RealParm2),INTENT(OUT)::WT(Nfirst:Nlast)
		REAL(dp)::y2
		REAL(dp)::cphi,sphi
		COMPLEX(dp)::g11(0:N-1)
		COMPLEX(dp)::h(0:N-1),tmp
		g11=1d0
		cphi=COS(phi)
		sphi=SIN(phi)
		DO I=0,N-1,2
			y2=(I+N2)*y_step+p
			g11(INDS(I))=outfunc(explt(I),cphi,sphi,alpha,beta)
!			g11(INDS(I))=outfunc(y2,phi,alpha,beta)

			y2=y2+y_step
!			g11(INDS(I+1))=-outfunc(y2,phi,alpha,beta)
			g11(INDS(I+1))=-outfunc(explt(I+1),cphi,sphi,alpha,beta)
		ENDDO
		CALL FFT_16(g11,W_fwd)
		h=g11/inputfunc;
		DO J=0,N-1,2
			h(INDS(J))=g11(J)/inputfunc(J)
			h(INDS(J+1))=-g11(J+1)/inputfunc(J+1)
		ENDDO
		CALL FFT_16(h,W_bwd)
		DO J=1,N-1,2
			h(J-1)=h(J-1)/N
			h(J)=-h(J)/N
		ENDDO
		
		WT=REAL(h(Nfirst:Nlast),KIND=REALPARM2);
	END SUBROUTINE
SUBROUTINE VBTransformWeightsAllInt4(x,y,hx,hy,WT)
	IMPLICIT NONE
	REAL(RealParm),INTENT(IN)::x,y,hx,hy
	REAL(RealParm),INTENT(OUT)::WT(Nfirst:Nlast,6)
		TYPE (EXP_DATA_TYPE)::edt
	REAL(dp)::rho,phi,alpha,beta
	REAL(dp)::xq,yq,hxq,hyq
	PROCEDURE(outfunction2),POINTER::outfunc
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
	CALL  PrepareExps4(phi,alpha,beta,edt)
	outfunc=>outfunci4cl2
	CALL CalcWeights2(edt,WT(:,1),outfunc)

	outfunc=>outfunci4dxdxcl2
	CALL CalcWeights2(edt,WT(:,2),outfunc)
	outfunc=>outfunci4dxdycl2
	CALL CalcWeights2(edt,WT(:,3),outfunc)
	outfunc=>outfunci4dxcl2
	CALL CalcWeights2(edt,WT(:,4),outfunc)

	alpha=hyq/rho;
	beta=hxq/rho;
	phi=ATAN2(xq,yq);
!	CALL  PrepareExps4(phi,alpha,beta,edt)

	outfunc=>outfunci4dydycl2
	CALL CalcWeights2(edt,WT(:,5),outfunc)
	outfunc=>outfunci4dycl2
	CALL CalcWeights2(edt,WT(:,6),outfunc)
	
	WT(:,1)=WT(:,1)/hx/hy*rp(1)
	WT(:,2)=WT(:,2)/hx/hy*rp(2) 
	WT(:,3)=WT(:,3)/hx/hy*rp(3)
	WT(:,4)=WT(:,4)/hx/hy*rp(4)
	WT(:,5)=WT(:,5)/hx/hy*rp(5)
	WT(:,6)=WT(:,6)/hx/hy*rp(6) 
END SUBROUTINE
SUBROUTINE VBTransformWeightsAllInt42(x,y,hx,hy,WT)
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
SUBROUTINE VBTransformWeightsAllInt22(x,y,hx,hy,WT)
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
	CALL  PrepareExps2(phi,alpha,beta,edt)
	out2=>outfunc_int2_all
	CALL CalcWeights22(edt,WT,out2)

	
	WT(:,1)=WT(:,1)*rho
	WT(:,2)=WT(:,2)rho 
	WT(:,3)=WT(:,3)/rho
	WT(:,5)=WT(:,5)/rho
END SUBROUTINE
SUBROUTINE CalcWeights2(edt,WT,outfunc)
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		PROCEDURE(outfunction2),POINTER,INTENT(IN)::outfunc
		REAL(RealParm2),INTENT(OUT)::WT(Nfirst:Nlast)
		REAL(dp)::y2
		REAL(dp)::cphi,sphi
		COMPLEX(dp)::g11(0:N-1)
		COMPLEX(dp)::h(0:N-1),tmp
		g11=1.0_dp
		DO I=0,N-1,2
			g11(INDS(I))=outfunc(edt,I)
			g11(INDS(I+1))=-outfunc(edt,I+1)
		ENDDO
		CALL FFT_16(g11,W_fwd)
		h=g11/f1;
		DO J=0,N-1,2
			h(INDS(J))=g11(J)/f1(J)
			h(INDS(J+1))=-g11(J+1)/f1(J+1)
		ENDDO
		CALL FFT_16(h,W_bwd)
		DO J=1,N-1,2
			h(J-1)=h(J-1)/N
			h(J)=-h(J)/N
		ENDDO
		
		WT=REAL(h(Nfirst:Nlast),KIND=REALPARM2);
	END SUBROUTINE
SUBROUTINE CalcWeights22(edt,WT,outfunc)
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		PROCEDURE(out_all),POINTER::outfunc
		REAL(RealParm2),INTENT(OUT)::WT(Nfirst:Nlast,6)
		COMPLEX(dp)::g11(0:N-1,6)
		COMPLEX(dp)::h(0:N-1,6),tmp
		g11=1.0_dp
		DO I=0,N-1,2
			g11(INDS(I),:)=outfunc(edt,I)
			g11(INDS(I+1),:)=-outfunc(edt,I+1)
		ENDDO
		DO I=1,6
			CALL FFT_16(g11(:,I),W_fwd)
			DO J=0,N-1,2
				h(INDS(J),I)=g11(J,I)/f1(J)
				h(INDS(J+1),I)=-g11(J+1,I)/f1(J+1)
			ENDDO
			CALL FFT_16(h(:,I),W_bwd)
			DO J=1,N-1,2
				h(J-1,I)=h(J-1,I)/N
				h(J,I)=-h(J,I)/N
			ENDDO
		ENDDO
		
		WT=REAL(h(Nfirst:Nlast,:),KIND=REALPARM2);
	END SUBROUTINE
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
		edt%x1=explt*cphi
		edt%y1=explt*sphi
		edt%dx=alpha*explt
		edt%dy=beta*explt
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
		edt%x1=explt*cphi*0.5_dp
		edt%y1=explt*sphi*0.5_dp
		edt%dx=alpha*explt*0.5_dp
		edt%dy=beta*explt*0.5_dp

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

FUNCTION outfunci4cl(lt,cphi,sphi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,dxr,dyr
		REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t3,p
!	t=EXP(lt)
	t=lt
	t3=t*t*t
	x1=t*cphi
	y1=t*sphi
	dx=dxr*t
	dy=dyr*t
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d2xd2y(x1,y1,dx,dy)
	R=R/t3
END FUNCTION

FUNCTION outfunci4dxdxcl(lt,cphi,sphi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,dxr,dyr
		REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t3,p
!	t=EXP(lt)
	t=lt
	x1=t*cphi
	y1=t*sphi
	dx=dxr*t
	dy=dyr*t
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d0xd2y(x1,y1,dx,dy)
	R=R/t
END FUNCTION

FUNCTION outfunci4dxdycl(lt,cphi,sphi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,dxr,dyr
		REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t2,p
!	t=EXP(lt)
	t=lt	
	x1=t*cphi
	y1=t*sphi
	dx=dxr*t
	dy=dyr*t
!	PRINT*, x1/dx*2.0_dp-1.0_dp,y1/dy*2.0_dp-1.0_dp
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d1xd1y(x1,y1,dx,dy)
	R=R/t
END FUNCTION
FUNCTION outfunci4dxcl(lt,cphi,sphi,dxr,dyr) RESULT(R)	
	REAL(dp),INTENT(IN)::lt,dxr,dyr
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t2,p
!	t=EXP(lt)
	t=lt
	t2=t*t
	x1=t*cphi
	y1=t*sphi
	dx=dxr*t
	dy=dyr*t
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d1xd2y(x1,y1,dx,dy)
	R=R/t2
END FUNCTION

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

		R(1)=R(1)/explt3(I)
		R(2)=R(2)/explt(I)
		R(3)=R(3)/explt(I)
		R(4)=R(4)/explt2(I)
		R(5)=R(5)/explt(I)
		R(6)=R(6)/explt2(I)
END FUNCTION
FUNCTION outfunci4cl2(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R
		REAL(dp)::fx(1:3),fy(1:3)
		fx=INT_TEXPd4_2d_edt(edt,I,JX)
		fy=INT_TEXPd4_2d_edt(edt,I,JY)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
		R=R/explt3(I)
END FUNCTION

FUNCTION outfunci4dxdxcl2(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R
		REAL(dp)::fx(1:3),fy(1:3)
		fx=INT_D2TEXPd4_2d_edt(edt,I,JX)
		fy=INT_TEXPd4_2d_edt(edt,I,JY)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
		R=R/explt(I)
END FUNCTION

FUNCTION outfunci4dxdycl2(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R
		REAL(dp)::fx(1:3),fy(1:3)
		fx=INT_D1TEXPd4_2d_edt(edt,I,JX)
		fy=INT_D1TEXPd4_2d_edt(edt,I,JY)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
		R=R/explt(I)
END FUNCTION
FUNCTION outfunci4dxcl2(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R
		REAL(dp)::fx(1:3),fy(1:3)
		fx=INT_D1TEXPd4_2d_edt(edt,I,JX)
		fy=INT_TEXPd4_2d_edt(edt,I,JY)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
		R=R/explt2(I)
END FUNCTION
FUNCTION outfunci4dydycl2(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R
		REAL(dp)::fx(1:3),fy(1:3)
		fx=INT_TEXPd4_2d_edt(edt,I,JX)
		fy=INT_D2TEXPd4_2d_edt(edt,I,JY)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
		R=R/explt(I)
END FUNCTION
FUNCTION outfunci4dycl2(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R
		REAL(dp)::fx(1:3),fy(1:3)
		fx=INT_TEXPd4_2d_edt(edt,I,JX)
		fy=INT_D1TEXPd4_2d_edt(edt,I,JY)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
		R=R/explt2(I)
END FUNCTION
!------------------------------ ---------------------------------------!
!lt is distance between centers of two rectangles of the same size dxr*dyr

FUNCTION outfunci4c(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t3,p,q
!	t=EXP(lt)
	t=lt
	t3=t*t*t
		x1=t*COS(phi)
		y1=t*SIN(phi)
		dx=dxr*t
		dy=dyr*t
	R=INT_t4EXP_d2xd2y(x1,y1,dx,dy)
		R=R/t3
END FUNCTION
!-------------------------------- INT2 ------------------------------------!
!lt is distance between bottom left corner of the  rectangle and point.

FUNCTION outfunci2cl(lt,cphi,sphi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,dxr,dyr
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,a,b,f0,f1,f2,g0,g1,g2
	CALL CalcXY2(lt,cphi,sphi,dxr,dyr,x1,x2,y1,y2,t) 
	f0=INT_I0(y1,y2)
	f1=INT_I1(y1,y2)
	f2=INT_I2(y1,y2)
	
	g0=INT_I0(x1,x2)
	g1=INT_I1(x1,x2)
	g2=INT_I2(x1,x2)

	R=f0*g2+2.0_dp*f1*g1+f2*g0

	R=R/t

END FUNCTION
FUNCTION outfunci2dxdxcl(lt,cphi,sphi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,dxr,dyr	
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,f0,f1,f2,g11,g31,g51,g12,g32,g52
	CALL CalcXY2(lt,cphi,sphi,dxr,dyr,x1,x2,y1,y2,t) 
	f0=INT_I0(y1,y2)*2.0_dp
	f1=INT_I1(y1,y2)*2.0_dp
	f2=INT_I2(y1,y2)*2.0_dp

	g11=XEXP2(x1)
	g12=XEXP2(x2)
	g31=x1*x1*g11
	g32=x2*x2*g12
	g51=x1*x1*g31
	g52=x2*x2*g32

	R=8.0_dp*f0*(g32-g31)
	R=R-4.0_dp*f0*(g52-g51)
	R=R+2.0_dp*(g12-g11+g31-g32)*f1
	R=R-0.25_dp*(g12-g11)*f2
	R=R*t

END FUNCTION
FUNCTION outfunci2dxdycl(lt,cphi,sphi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,dxr,dyr
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,f11,f12,f21,f22,r1
	CALL CalcXY2(lt,cphi,sphi,dxr,dyr,x1,x2,y1,y2,t)
	r1=(x1*x1+y1*y1)
	f11=EXP(-r1)*r1*r1*4.0_dp

	r1=(x2*x2+y1*y1)
	f21=EXP(-r1)*r1*r1*4.0_dp

	r1=(x1*x1+y2*y2)
	f12=EXP(-r1)*r1*r1*4.0_dp

	r1=(x2*x2+y2*y2)
	f22=EXP(-r1)*r1*r1*4.0_dp
	R=(f22-f12)-(f21-f11)
	R=R*t
END FUNCTION
FUNCTION outfunci2dxcl(lt,cphi,sphi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,dxr,dyr
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,f0,f1,f2,g01,g21,g41,g02,g22,g42
	CALL CalcXY2(lt,cphi,sphi,dxr,dyr,x1,x2,y1,y2,t)
	f0=INT_I0(y1,y2)
	f1=INT_I1(y1,y2)
	f2=INT_I2(y1,y2)
	g01=EXP(-x1*x1)
	g02=EXP(-x2*x2)
	
	g21=x1*x1*g01
	g22=x2*x2*g02

	g41=x1*x1*g21
	g42=x2*x2*g22
	
	R=8.0_dp*f0*(g42-g41)+4.0_dp*f1*(g22-g21)+f2*(g02-g01)*0.5_dp
END FUNCTION
!---------------------------------------------------------------------!
!lt is distance between centers of rectangle and point
FUNCTION outfunci2c(lt,cphi,sphi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,dxr,dyr
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,a,b,f1,f2,dx,dy,g1,g2
!	t=EXP(lt)
	t=lt
	x1=t*cphi
	y1=t*sphi
	dx=dxr*t
	dy=dyr*t
	x1=x1-dx/2.0_dp
	y1=y1-dy/2.0_dp
	x2=x1+dx
	y2=y1+dy
	f1=ERRFD(x1,x2,1.0_dp)*SPI;
	f2=ERRFD(y1,y2,1.0_dp)*SPI;
	g1=XEXP2(x2/2.0_dp)-XEXP2(x1/2.0_dp)
	g2=XEXP2(y2/2.0_dp)-XEXP2(y1/2.0_dp)
	R=((f1-2.0_dp*g1)*f2+(f2-2.0_dp*g2)*f1)/4.0_dp!/dx/dy/8.0_dp
	R=R/t
END FUNCTION

FUNCTION outfunc_int2_all(edt,I) RESULT(R)	
		TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
		INTEGER,INTENT(IN)::I
		REAL(dp)::R(6)
		REAL(dp)::fx(1:3,3),fy(1:3,3)
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

		R(1)=R(1)/explt(I)
		R(2)=R(2)*explt(I)
		R(3)=R(3)*explt(I)
		R(5)=R(5)*explt(I)
END FUNCTION
!-------------------------------------------------------------------------------------!
SUBROUTINE CalcXY2(lt,cphi,sphi,dxr,dyr,x1,x2,y1,y2,t) 
	REAL(dp),INTENT(IN)::lt,dxr,dyr
	REAL(dp),INTENT(IN)::cphi,sphi
	REAL(dp),INTENT(OUT)::x1,x2,y1,y2,t
	REAL(dp)::dx,dy
!	t=EXP(lt)
t=lt
		x1=t*cphi
		y1=t*sphi
		dx=dxr*t
		dy=dyr*t
	x2=x1+dx
	y2=y1+dy
	x1=x1/2.0_dp
	y1=y1/2.0_dp
	x2=x2/2.0_dp
	y2=y2/2.0_dp
END SUBROUTINE
FUNCTION INT_t4EXP_d2xd2y(lx,ly,hx,hy) RESULT(R)
	REAL(dp),INTENT(IN)::lx,ly,hx,hy
	REAL(dp)::R
	REAL(dp)::fx(1:3),fy(1:3)
	fx=INT_TEXPd4_2d(lx,hx)
	fy=INT_TEXPd4_2d(ly,hy)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
END FUNCTION

FUNCTION INT_t4EXP_d0xd2y(lx,ly,hx,hy) RESULT(R)
	REAL(dp),INTENT(IN)::lx,ly,hx,hy
	REAL(dp)::R
	REAL(dp)::fx(1:3),fy(1:3)
	fx=INT_D2TEXPd4_2d(lx,hx)
	fy=INT_TEXPd4_2d(ly,hy)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
END FUNCTION

FUNCTION INT_t4EXP_d1xd1y(lx,ly,hx,hy) RESULT(R)
	REAL(dp),INTENT(IN)::lx,ly,hx,hy
	REAL(dp)::R
	REAL(dp)::fx(1:3),fy(1:3)
	fx=INT_D1TEXPd4_2d(lx,hx)
	fy=INT_D1TEXPd4_2d(ly,hy)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
END FUNCTION

FUNCTION INT_t4EXP_d1xd2y(lx,ly,hx,hy) RESULT(R)
	REAL(dp),INTENT(IN)::lx,ly,hx,hy
	REAL(dp)::R
	REAL(dp)::fx(1:3),fy(1:3)
	fx=INT_D1TEXPd4_2d(lx,hx)
	fy=INT_TEXPd4_2d(ly,hy)
		R=fx(3)*fy(1)+fx(1)*fy(3)+2.0_dp*fx(2)*fy(2)
END FUNCTION
ELEMENTAL FUNCTION IX4EXP2(t) RESULT(R) !\int e^{-x^2}x^4dx
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
		REAL(dp)::x1,x2
	R=-5.e-1_dp*X3EXP2(t)+1.5_dp*IX2EXP2(t)
END FUNCTION

ELEMENTAL FUNCTION IX2EXP2(t) RESULT(R) !\int e^{-x^2}x^2dx
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
		REAL(dp)::x1,x2
	R=-5.e-1_dp*XEXP2(t)+IEXP2(t)*5e-1_dp
END FUNCTION

ELEMENTAL FUNCTION IEXP2(t) RESULT(R) !\int e^{-x^2}dx=\frac{\sqrt{\pi}}}{2}Erf(x)
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
		REAL(dp)::x1,x2
	R=SPI/2.0_dp*ERF(t)
END FUNCTION

ELEMENTAL FUNCTION X3EXP2(t) RESULT(R)!e^{-x^2}x^3
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
		REAL(dp)::t2
	t2=t*t
	R=t*t2*EXP(-t2)
END FUNCTION
ELEMENTAL FUNCTION XEXP2(t) RESULT(R)!e^{-x^2}x
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
		REAL(dp)::x1,x2
	R=t*EXP(-t*t)
END FUNCTION

FUNCTION ERRFD(t1,t2,a) RESULT(R)
	REAL(dp),INTENT(IN)::t1,t2,a
	REAL(dp)::R
		REAL(dp)::x1,x2
		x1=t1/2.0_dp/sqrt(a)
		x2=t2/2.0_dp/sqrt(a)
!		 PRINT*,x1,x2
		R=ERF(x2)-ERF(x1)
END FUNCTION


FUNCTION INT_D2TEXPd4_2d(l,h) RESULT(R)
	REAL(dp),INTENT(IN)::l,h
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3),l2,h2,x(1:3),x2(1:3),ex(1:3)
	x(1)=l
	x(2)=l+h
	x(3)=l-h
	x2=x*x
	ex=EXP(-x2/4.0_dp)/2.0_dp
	f=x2*ex
	R(2)=f(2)+f(3)-2.0_dp*f(1)
	f=x2*f
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	f=ex
	R(3)=f(2)+f(3)-2.0_dp*f(1)

END FUNCTION

FUNCTION INT_D1TEXPd4_2d(l,h) RESULT(R)
	REAL(dp),INTENT(IN)::l,h
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3),l2,h2
	l2=l*0.5_dp
	h2=h*0.5_dp
	f=INT_1d(l2,h2)
	R(1)=-8.0_dp*f(3)-12.0_dp*f(2)+12.0_dp*f(1);
	R(2)=-2.0_dp*f(2)+2.0_dp*f(1)
		R(3)=f(1)
!R=R*0.5_dp
END FUNCTION
FUNCTION INT_1d(l,h)RESULT(R)
	REAL(dp),INTENT(IN)::l,h
	REAL(dp)::R(1:3)
	REAL(dp)::x(1:3),f(1:3),x2(1:3)
	x(1)=l
	x(2)=l+h
	x(3)=l-h
		x2=x*x
	f=ERF(x)
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	R(1)=R(1)*SPI*0.5_dp

	f=x*EXP(-x2)
	R(2)=f(2)+f(3)-2.0_dp*f(1)
	f=f*x2
	R(3)=f(2)+f(3)-2.0_dp*f(1)
END FUNCTION

FUNCTION INT_TEXPd4_2d(l,h) RESULT(R)
	REAL(dp),INTENT(IN)::l,h
	REAL(dp)::R(1:3)
	REAL(dp)::f(1:3),l2,h2
	l2=l*0.5_dp
	h2=h*0.5_dp
	f=INT_2d(l2,h2)
	R(1)=32.0_dp*(1.5_dp*f(1)+f(2)+f(3))
	R(2)=8.0_dp*f(1)+4.0_dp*f(2)
	R(3)=4.0_dp*f(1)+f(2);
END FUNCTION


FUNCTION INT_2d(l,h)RESULT(R)
	REAL(dp),INTENT(IN)::l,h
	REAL(dp)::R(1:3)
	REAL(dp)::x(1:3),f(1:3)
	x(1)=l
	x(2)=l+h
	x(3)=l-h

	f=x*ERF(x)
	R(1)=f(2)+f(3)-2.0_dp*f(1)
	R(1)=R(1)*SPI*0.25_dp

	f=EXP(-x*x);
	R(2)=f(2)+f(3)-2.0_dp*f(1)!h should be lagrger than 1e-5l
   !	 R(2)=EXP(-l*l)*(-2.0_dp+4.0*l*l)*h*h
	f=EXP(-x*x)*x*x/4.0_dp
	R(3)=f(2)+f(3)-2.0_dp*f(1)
END FUNCTION

FUNCTION INT_I0(y1,y2) RESULT(R)
	REAL(dp),INTENT(IN)::y1,y2
	REAL(dp)::R
	REAL(dp)::f1,f2
	f1=ERF(y1)
	f2=ERF(y2)
	R=SPI*(f2-f1)*0.5_dp
END FUNCTION

FUNCTION INT_I1(y1,y2) RESULT(R)
	REAL(dp),INTENT(IN)::y1,y2
	REAL(dp)::R
	REAL(dp)::f1,f2
	f1=ERF(y1)
	f2=ERF(y2)
	R=SPI*(f2-f1)
	f1=XEXP2(y1)
	f2=XEXP2(y2)
	R=R-2.0_dp*(f2-f1)
END FUNCTION


FUNCTION INT_I2(y1,y2) RESULT(R)
	REAL(dp),INTENT(IN)::y1,y2
	REAL(dp)::R
	REAL(dp)::f1,f2,l2,h2
	f1=ERF(y1)
	f2=ERF(y2)
	R=6.0_dp*SPI*(f2-f1)
	f1=XEXP2(y1)
	f2=XEXP2(y2)
	R=R-12.0_dp*(f2-f1)
	f1=X3EXP2(y1)
	f2=X3EXP2(y2)
	R=R-8.0_dp*(f2-f1)
END FUNCTION
!-----------------------------------------------------------------------------!

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
	R=R/2.0_dp
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
	R=R/2.0_dp
END FUNCTION

FUNCTION SINGLE_INT_D2(edt,I,J) RESULT(R)
	TYPE (EXP_DATA_TYPE),INTENT(IN)::edt
	INTEGER,INTENT(IN)::I,J
	REAL(dp)::R(3)
	REAL(dp)::f(2),a,b,c
	f=edt%xy(1:2,I,J)*edt%xy_exp2(1:2,I,J)
	a=f(2)-f(1)
	R(1)=-a
	f=f*edt%xy2(1:2,I,J)
	b=f(2)-f(1)
	R(2)=-b+a*2.0_dp
	f=f*edt%xy2(1:2,I,J)
	c=f(2)-f(1)
	R(3)=-c+4.0*b
END FUNCTION


!-------------------------------------------------------------!

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
		INTEGER :: I,J,J2,I2,Istep,Ip,Iq
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
	END SUBROUTINE FFT_16
END
