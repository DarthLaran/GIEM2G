
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

INTERFACE
	FUNCTION outfunction(lt,phi,dxr,dyr) RESULT(R)
		INTEGER ,PARAMETER::dp=16
		REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
		REAL(dp)::R
	END FUNCTION
END INTERFACE

PUBLIC ::VBTransformWeights,VBTransformInit 

CONTAINS


SUBROUTINE VBTransformInit(lms) !NOT THREAD SAFETY!
	REAL(RealParm2),INTENT(OUT)::lms(Nfirst:Nlast)
		COMPLEX(dp)::f11(0:N-1),f21(0:N-1),tmp
	REAL(dp)::y1,theta
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
				CASE  (8)  ! (INT4DYY)
						outfunc=>outfunci4dxdxcl
						inputfunc=>f1
			rp=rho
				CASE (10) !(INT4DY)
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
			COMPLEX(dp)::g11(0:N-1)
		COMPLEX(dp)::h(0:N-1),tmp
		g11=1d0
		DO I=0,N-1,2
			y2=(I+N2)*y_step+p
			g11(INDS(I))=outfunc(y2,phi,alpha,beta)

			y2=y2+y_step
			g11(INDS(I+1))=-outfunc(y2,phi,alpha,beta)
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


!------------------------ Output Functions	-------------------------------------!
!----------------- INT4 ---------------------------------------------------------!
!lt is distance between bottom left corner of the first rectangle and center of the second. The sizes of rectangles are the same

FUNCTION outfunci4cl(lt,phi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t3,p
	t=EXP(lt)
	t3=t*t*t
		x1=t*COS(phi)
		y1=t*SIN(phi)
		dx=dxr*t
		dy=dyr*t
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d2xd2y(x1,y1,dx,dy)
	R=R/t3
END FUNCTION

FUNCTION outfunci4dxdxcl(lt,phi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t3,p
	t=EXP(lt)
	t3=t*t*t
		x1=t*COS(phi)
		y1=t*SIN(phi)
		dx=dxr*t
		dy=dyr*t
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d0xd2y(x1,y1,dx,dy)
	R=R/t
END FUNCTION

FUNCTION outfunci4dxdycl(lt,phi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t2,p
	t=EXP(lt)
	t2=t*t
		x1=t*COS(phi)
		y1=t*SIN(phi)
		dx=dxr*t
		dy=dyr*t
!	PRINT*, x1/dx*2.0_dp-1.0_dp,y1/dy*2.0_dp-1.0_dp
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d1xd1y(x1,y1,dx,dy)
	R=R/t
END FUNCTION
FUNCTION outfunci4dxcl(lt,phi,dxr,dyr) RESULT(R)	
		REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t2,p
	t=EXP(lt)
	t2=t*t
		x1=t*COS(phi)
		y1=t*SIN(phi)
		dx=dxr*t
		dy=dyr*t
	x1=x1-dx*0.5_dp
	y1=y1-dy*0.5_dp
	R=INT_t4EXP_d1xd2y(x1,y1,dx,dy)
	R=R/t2
END FUNCTION


FUNCTION outfunci4dxdx2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
		R=outfunci4dxdxcl(lt,phi,dxr,dyr)*EXP(-lt)
END FUNCTION

FUNCTION outfunci4dx2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
		R=outfunci4dxcl(lt,phi,dxr,dyr)*EXP(-lt)
END FUNCTION

FUNCTION outfunci4dxdy2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::phi2
!		phi2=ATAN2(lt*COS(phi),lt*SIN(phi))
		R=outfunci4dxdycl(lt,phi,dxr,dyr)*EXP(-lt)
END FUNCTION
!------------------------------ ---------------------------------------!
!lt is distance between centers of two rectangles of the same size dxr*dyr

FUNCTION outfunci4c(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t3,p,q
	t=EXP(lt)
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

FUNCTION outfunci2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,a,b,f0,f1,f2,g0,g1,g2
	CALL CalcXY2(lt,phi,dxr,dyr,x1,x2,y1,y2,t) 
	f0=INT_I0(y1,y2)
	f1=INT_I1(y1,y2)
	f2=INT_I2(y1,y2)
	
	g0=INT_I0(x1,x2)
	g1=INT_I1(x1,x2)
	g2=INT_I2(x1,x2)

	R=f0*g2+2.0_dp*f1*g1+f2*g0

	R=R/t

END FUNCTION
FUNCTION outfunci2dxdxcl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,f0,f1,f2,g11,g31,g51,g12,g32,g52
	CALL CalcXY2(lt,phi,dxr,dyr,x1,x2,y1,y2,t) 
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
FUNCTION outfunci2dxdycl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,f11,f12,f21,f22,r1
	CALL CalcXY2(lt,phi,dxr,dyr,x1,x2,y1,y2,t)
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
FUNCTION outfunci2dxcl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,f0,f1,f2,g01,g21,g41,g02,g22,g42
	CALL CalcXY2(lt,phi,dxr,dyr,x1,x2,y1,y2,t)
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
FUNCTION outfunci2dxdx2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
		R=outfunci2dxdxcl(lt,phi,dxr,dyr)*EXP(-lt)
END FUNCTION

FUNCTION outfunci2dx2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
		R=outfunci2dxcl(lt,phi,dxr,dyr)*EXP(-lt)
END FUNCTION

FUNCTION outfunci2dxdy2cl(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::phi2
!		phi2=ATAN2(lt*COS(phi),lt*SIN(phi))
		R=outfunci2dxdycl(lt,phi,dxr,dyr)*EXP(-lt)
END FUNCTION
!---------------------------------------------------------------------!
!lt is distance between centers of rectangle and point
FUNCTION outfunci2c(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,x2,y1,y2,a,b,f1,f2,dx,dy,g1,g2
	t=EXP(lt)
		x1=t*COS(phi)
		y1=t*SIN(phi)
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

!-------------------------------------------------------------------------------------!
SUBROUTINE CalcXY2(lt,phi,dxr,dyr,x1,x2,y1,y2,t) 
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp),INTENT(OUT)::x1,x2,y1,y2,t
	REAL(dp)::dx,dy
	t=EXP(lt)
		x1=t*COS(phi)
		y1=t*SIN(phi)
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


!-------------------------------------------------------------!

FUNCTION inpf_a(t) RESULT (R)
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
	REAL(dp)::y,y2
	y=EXP(-t);
	y2=y*y;
	R=EXP(-y2-t)*(1.0_dp-EXP(-y2))
END FUNCTION

FUNCTION outpf_a(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t4
	t=EXP(lt)
	t2=t*t/8.0_dp
	R=EXP(-t2)*(EXP(-t2)-0.5_dp)/2.0_dp*t
END FUNCTION


 FUNCTION inpfd(t) RESULT (R)
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
	REAL(dp)::y,y2
	y=EXP(-t);
		y2=y*y;
		R=y*EXP(-y2)*(1.0_dp-y2);
END FUNCTION
 FUNCTION outpfd(t) RESULT (R)
	REAL(dp),INTENT(IN)::t
	REAL(dp)::R
	REAL(dp)::y,y2
	y=EXP(t);
		y2=y*y;
!		 R=EXP(-y2/4.0_dp)*y2/8.0_dp!*y
	R=y2*XEXP2(y/2.0_dp)/4.0_dp
END FUNCTION

FUNCTION outfunc4(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t4
	t=EXP(lt)
	t4=t*t*t*t
	R=EXP(-t*t/4.0_dp)*t4*t/4.0_dp
END FUNCTION


FUNCTION outfunc2(lt,phi,dxr,dyr) RESULT(R)
	REAL(dp),INTENT(IN)::lt,phi,dxr,dyr
	REAL(dp)::R
	REAL(dp)::t,x1,y1,dx,dy,t2
	t=EXP(lt)
	t2=t*t
	R=t2*XEXP2(t/2.0_dp)/4.0_dp
!		 PRINT*, R
END FUNCTION

!---------------------------   FFT for COMPLEX(16)	 for this module only!-----------------------------------------!

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
	SUBROUTINE fft_16(x,wp)!x MUST be after bit reverse, wp -array of coefficients for forward or backward FT
		COMPLEX(dp), DIMENSION(0:), INTENT(inout) :: x
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
	END SUBROUTINE fft_16
END
