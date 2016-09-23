!Copyright (c) 2016 Mikhail Kruglyakov 
!This file is part of GFGMRES.

!GFGMRES is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 2 of the License.

!GFGMRES is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with GFMRES.  If not, see <http://www.gnu.org/licenses/>.
MODULE FGMRES
      USE UTILITIES
      USE FGMRES_INTERFACES
      IMPLICIT NONE

      PRIVATE
      REAL(REALPARM):: ZERO_THRESHOLD=1d-14
      INTEGER,PARAMETER:: NumberOfIGSIterations=3

      TYPE FGMRES_DATA
	      LOGICAL::FLEX
	      COMPLEX(REALPARM),POINTER::x(:)
	      COMPLEX(REALPARM),POINTER::y(:)
	      COMPLEX(REALPARM),POINTER::r0(:)
	      COMPLEX(REALPARM),POINTER::w(:)
	      COMPLEX(REALPARM),POINTER::dotproducts(:)
	      COMPLEX(REALPARM),POINTER::KrylovBasis(:,:)
	      COMPLEX(REALPARM),POINTER::Hessenberg(:,:)
	      COMPLEX(REALPARM),POINTER::PrecondKrylovBasis(:,:)
	      COMPLEX(REALPARM),POINTER::GivensRotSin(:) 
	      REAL(REALPARM),POINTER::GivensRotCos(:)
	      INTEGER::Nloc
	      INTEGER::BasisSize
	      INTEGER::Iternum
	      TYPE(GMRES_PARAMETERS)::params
	      TYPE(RESULT_INFO)::info
	      PROCEDURE (MatrixVectorMult),POINTER,NOPASS::ApplyOperator  
	      PROCEDURE (MatrixVectorMult),POINTER,NOPASS::PrecondLeft 
	      PROCEDURE (MatrixVectorMult),POINTER,NOPASS::PrecondRight 
	      PROCEDURE (MANYDOTPRODUCTS),POINTER,NOPASS::MANYDP
	      PROCEDURE (InformAboutIteration),POINTER,NOPASS::Inform
	      TYPE(C_PTR)::A
	      TYPE(C_PTR)::MLeft
	      TYPE(C_PTR)::MRight
              TYPE(C_PTR)::dpinstance
      END TYPE

     PUBLIC:: FGMRES_DATA, GMRES_SOLVE, AllocateFGMRES,Set_FGMRES_Operators,DeallocateFGMRES

      CONTAINS


	SUBROUTINE AllocateFGMRES(solver,n,m,flex,l)
        	USE, INTRINSIC :: iso_c_binding
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		INTEGER,INTENT(IN)::n,m
		LOGICAL,OPTIONAL,INTENT(IN)::flex
                INTEGER(C_INTPTR_T),INTENT(OUT),OPTIONAL::l
                INTEGER(C_INTPTR_T)::length
		solver%Nloc=n
		solver%BasisSize=m
		ALLOCATE(solver%dotproducts(m))
		ALLOCATE(solver%Hessenberg(m+1,m+1))
		ALLOCATE(solver%KrylovBasis(n,m+1))
		ALLOCATE(solver%x(n),solver%y(m),solver%w(n),solver%r0(n))
		ALLOCATE(solver%GivensRotCos(m),solver%GivensRotSin(m))
                length=4*m+(m+1)**2+n*(m+1)+3*n
		IF (PRESENT(flex)) THEN
			IF(flex .EQV. FLEXIBLE) THEN
			       ALLOCATE(solver%PrecondKrylovBasis(n,m))
                               length=length+m*n
			ELSE
				NULLIFY(solver%PrecondKrylovBasis)
			ENDIF
			solver%FLEX=flex
		ELSE
		       solver%FLEX=.NOT. FLEXIBLE
		       NULLIFY(solver%PrecondKrylovBasis)
		ENDIF
                IF (PRESENT(l)) l=length
	ENDSUBROUTINE
	
	SUBROUTINE Set_FGMRES_Operators(solver,apply,precond_right,precond_left,dp,inform,A,M1,M2)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::apply
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::precond_right
		PROCEDURE(MatrixVectorMult),POINTER,INTENT(IN)::precond_left
		PROCEDURE(MANYDOTPRODUCTS),POINTER,INTENT(IN)::dp
		PROCEDURE (InformAboutIteration),POINTER,INTENT(IN)::Inform
		TYPE(C_PTR),INTENT(IN)::A,M1,M2
		solver%ApplyOperator=>apply
		solver%PrecondLeft=>precond_left
		solver%PrecondRight=>precond_Right
		IF (.NOT. ASSOCIATED(precond_right)) NULLIFY(solver%PrecondRight)
		IF (.NOT. ASSOCIATED(precond_left)) NULLIFY(solver%PrecondLeft)
		solver%A=A
		solver%MLeft=M1
		solver%MRight=M2
		solver%MANYDP=>dp
		solver%Inform=>Inform
                solver%dpinstance=C_NULL_PTR
		IF (.NOT. ASSOCIATED(Inform)) NULLIFY(solver%Inform)
	ENDSUBROUTINE


        RECURSIVE	SUBROUTINE GMRES_SOLVE(solver,x,x0,b,info)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		COMPLEX(REALPARM),INTENT(INOUT)::x(:)
		COMPLEX(REALPARM),INTENT(IN)::b(:),x0(:)
		TYPE(RESULT_INFO),INTENT(OUT)::info
		INTEGER :: Itnum,Iter
		REAL(REALPARM)::bn,sb,sPb
		REAL(REALPARM)::beta
		REAL(REALPARM)::be,bea
		INTEGER::I,J,K,jH
		COMPLEX(REALPARM),POINTER::r0(:)
		COMPLEX(REALPARM),POINTER::w(:)
		COMPLEX(REALPARM),POINTER::KrylovBasis(:,:)
		COMPLEX(REALPARM),POINTER::Hessenberg(:,:)
		INTEGER::N,M
		r0=>solver%r0
		w=>solver%w
		KrylovBasis=>solver%KrylovBasis
		Hessenberg=>solver%Hessenberg
                Hessenberg=C_ZERO
		N=solver%Nloc
		M=solver%BasisSize
		solver%Iternum=0
		Itnum=0
		Iter=0
		be=-R_ONE
		bea=-R_ONE


		bn=CalculateVectorNorm(solver,b)
		IF (bn< ZERO_THRESHOLD) THEN
			x=C_ZERO
			info%stat=CONVERGED
			info%Iterations=0
			RETURN
		ENDIF
!!!!!!!!!!!!---   Scales for preconditioned/unpreconditioned systems. Not!implemented
	       sb=bn
	       sPb=bn
!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
                CALL ZCOPY(N,x0,ONE,x,ONE) 
		CALL ApplyOperator(solver,x0,r0)
                !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(I)
                !$OMP DO SCHEDULE(GUIDED)
                        DO I=1,N
	        		r0(I)=b(I)-r0(I)
                        ENDDO
               !$OMP ENDDO
               !$OMP END PARALLEL 

		CALL ApplyLeftPreconditioner(solver,r0,w)
		beta=CalculateVectorNorm(solver,w)

		IF (beta/sb< solver%params%Tol) THEN
			info%stat=CONVERGED
			info%be=beta/sb
			info%Iterations=0
			RETURN
		ENDIF
		Iter=0

		DO 
                        Hessenberg(1,M+1)=beta
			!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(I)
			!$OMP DO SCHEDULE(GUIDED)
                        
				DO I=1,N
					KrylovBasis(I,1)=w(I)/beta
				ENDDO
			!$OMP ENDDO
			!$OMP DO SCHEDULE(GUIDED)
				DO I=1,M
					Hessenberg(I+1,M+1)=C_ZERO
				ENDDO
			!$OMP ENDDO
			!$OMP END PARALLEL
			jH=1
			bea=KrylovSubspaceConstruction(solver,bn,jH)
			IF (bea<0) THEN
				info=solver%info
				info%stat=INTERUPTED
				RETURN
			ELSEIF (bea <solver%params%Tol) THEN
				be=FinalActions(solver,jH, sPb, b,x)
				info%stat=solver%Info%stat
				info%Iterations=solver%Iternum
				info%be=be
				info%bea=bea
				RETURN
			ENDIF


			Iter=Iter+1

			IF (Iter >=solver%params%MaxIt) EXIT

			beta=PrepareRestart(solver,b,x,jH)
		 ENDDO
			be=FinalActions(solver,jH, sPb, b,x)
			info%bea=bea
			info%be=be
			info%Iterations=solver%Iternum
			info%stat=NON_CONVERGED
	ENDSUBROUTINE


        RECURSIVE	FUNCTION  KrylovSubspaceConstruction(solver,bn, jH) RESULT (RES)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		REAL(REALPARM),INTENT(IN)::bn
		INTEGER,INTENT(INOUT)::jH
		REAL(REALPARM)::RES

		COMPLEX(REALPARM),POINTER::Hessenberg(:,:)
		COMPLEX(REALPARM),POINTER::w(:),r0(:)
		COMPLEX(REALPARM),POINTER::KrylovBasis(:,:)
		COMPLEX(REALPARM),POINTER::dotproducts(:)
		COMPLEX(REALPARM),POINTER::GivensRotSin(:) 
		REAL(REALPARM),POINTER::GivensRotCos(:)
		COMPLEX(REALPARM),POINTER::tmpptr(:) 
		REAL(REALPARM)::bea=-R_ONE
		INTEGER::Northo=ZERO
		REAL(REALPARM)::dloo,dnormw,dnormres
                COMPLEX(REALPARM)::H1(1),H2(1)
                REAL(REALPARM)::DZNRM2

		INTEGER::N,M
		INTEGER::I,J
		INTEGER::interp
		w=>solver%w
		r0=>solver%r0
		KrylovBasis=>solver%KrylovBasis
		Hessenberg=>solver%Hessenberg
		dotproducts=>solver%dotproducts
		GivensRotSin=>solver%GivensRotSin
		GivensRotCos=>solver%GivensRotCos

		N=solver%Nloc
		M=solver%BasisSize

		DO    
                        tmpptr(1:N)=>KrylovBasis(:,jH)
			CALL ApplyRightPreconditioner(solver,tmpptr,w, jH);

                        CALL ApplyOperator(solver,w,r0)
			CALL ApplyLeftPreconditioner(solver,r0,w);
			!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
			!$OMP DO SCHEDULE(GUIDED)
			DO I=1,jH
			    Hessenberg(I, jH) = C_ZERO
			ENDDO
			!$OMP ENDDO
			!$OMP ENDPARALLEL

			dloo=R_ZERO 
			Northo=ZERO
			DO
			    Northo=Northo+1
			    
                            CALL   CalculateDotProduct(solver,KrylovBasis,&
                                    &w,jH,dotproducts)


			    CALL ZAXPY(jH,C_ONE,dotproducts,ONE,Hessenberg(:,jH),ONE)


			    CALL ZGEMV('N',N,jH,-C_ONE,KrylovBasis,N,dotproducts,ONE,C_ONE,w,ONE)

			    dloo=DZNRM2(jH,dotproducts,ONE)

                            dnormw = CalculateVectorNorm(solver,w)
			    IF ((2.0*dnormw>dloo) .OR. (Northo > NumberOfIGSIterations)) EXIT
			ENDDO  

			Hessenberg(jH+1, jH) = dnormw


			!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
			!$OMP DO SCHEDULE(GUIDED)
			DO I=1,N 
				KrylovBasis(I,jH+1)=w(I)/dnormw
			ENDDO
			!$OMP ENDDO
			!$OMP ENDPARALLEL


			DO I=1,jH-1
                                CALL	ZROT(ONE,Hessenberg(I, jH), ONE,Hessenberg(I + 1, jH),ONE, GivensRotCos(I), GivensRotSin(I))
			ENDDO

                        H1(1)=Hessenberg(jH, jH)
                        H2(1)=Hessenberg(jH + 1, jH)
			CALL ZROTG(H1, H2, GivensRotCos(jH),GivensRotSin(jH))


			CALL ZROT(ONE,Hessenberg(jH, M+1),ONE,Hessenberg(jH + 1, M+1),ONE ,GivensRotCos(jH), GivensRotSin(jH))

			CALL ZROT(ONE,Hessenberg(jH, jH),ONE,Hessenberg(jH + 1, jH),ONE ,GivensRotCos(jH), GivensRotSin(jH))

			Hessenberg(jH + 1, jH) =C_ZERO;
			dnormres=ABS(Hessenberg(jH+1,M+1))

			bea = dnormres / bn;



			solver%IterNum=solver%IterNum+1
			solver%info%bea=bea
			solver%info%iterations=solver%IterNum
			solver%info%stat=INPROCESS
			interp=Inform(solver)
			IF (interp==INTERUPT_SOLVER) THEN
				res=-R_ONE
				RETURN
			ENDIF

			jH=jH+1
			IF ( (bea < solver%params%Tol).OR. (jH>M)) EXIT
		ENDDO
	    res=bea
    ENDFUNCTION
    


	FUNCTION PrepareRestart(solver,rhs,solution,jH) RESULT(beta)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		INTEGER,INTENT(IN)::jH
		COMPLEX(REALPARM),TARGET,INTENT(IN)::rhs(:)
		COMPLEX(REALPARM),TARGET,INTENT(INOUT)::solution(:)
		COMPLEX(REALPARM),POINTER::x(:),w(:),r0(:)
		COMPLEX(REALPARM),POINTER::Hessenberg(:,:)
		COMPLEX(REALPARM),POINTER::KrylovBasis(:,:)
		COMPLEX(REALPARM),POINTER::GivensRotSin(:) 
		REAL(REALPARM),POINTER::GivensRotCos(:)
		REAL(REALPARM)::beta		    
		INTEGER::M,N
		INTEGER::I,J
		N=solver%Nloc
		M=solver%BasisSize
		x(1:N)=>solver%x
		w(1:N)=>solver%w
		r0(1:N)=>solver%r0
		Hessenberg=>solver%Hessenberg
		KrylovBasis=>solver%KrylovBasis
		GivensRotSin=>solver%GivensRotSin
		GivensRotCos=>solver%GivensRotCos
		CALL ConstructCurrentSolution(solver,jH,solution)

		IF (solver%params%RESTART_RESIDUAL .EQV.  TRUE_RESIDUAL_AT_RESTART) THEN
			CALL ZCOPY(N,x,ONE,w,ONE)
			CALL ApplyOperator(solver,w,r0)
			!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
			!$OMP DO SCHEDULE(GUIDED)
				DO I=1,N
					r0(I)=rhs(I)-r0(I)
				ENDDO
			!$OMD ENDDO
			!$OMP END PARALLEL
			CALL ApplyLeftPreconditioner(solver,r0,w)
			beta=CalculateVectorNorm(solver,w)
		ELSE
			beta=ABS(Hessenberg(M+1,M+1))
			DO J=M,1,-1
				Hessenberg(J,M+1)=C_ZERO
				CALL ZROT(ONE,Hessenberg(J,M+1),ONE,Hessenberg(J+1,M+1),&
					&ONE,GivensRotCos(J),-GivensRotSin(J))
			ENDDO
			CALL ZGEMV('N',N,M+1,C_ONE,KrylovBasis,N,Hessenberg(:,M+1),ONE,C_ZERO,w,ONE)
		ENDIF

	ENDFUNCTION


	FUNCTION FinalActions(solver,jH, sPb, rhs, solution) RESULT(res)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		INTEGER,INTENT(IN)::jH
		REAL(REALPARM),INTENT(IN)::sPb
		COMPLEX(REALPARM),TARGET,INTENT(IN)::rhs(:)
                COMPLEX(REALPARM),TARGET,INTENT(INOUT)::solution(:)
		REAL(REALPARM)::res
		REAL(REALPARM)::TrueResNorm
		COMPLEX(REALPARM),POINTER::Hessenberg(:,:)
		COMPLEX(REALPARM),POINTER::x(:),r0(:),w(:)
		INTEGER::N,M
		INTEGER::I,J
		N=solver%Nloc
		M=solver%BasisSize
		x=>solver%x
		r0=>solver%r0
		w=>solver%w
                Hessenberg=>solver%Hessenberg
                CALL ConstructCurrentSolution(solver,jH, solution)
                
		CALL ZCOPY(N,solution,ONE,w,ONE)
		CALL ApplyOperator(solver,w,r0)
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
		!$OMP DO SCHEDULE(GUIDED)
			DO I=1,N
				r0(I)=rhs(I)-r0(I)
                                w(I)=r0(I)
			ENDDO
		!$OMD ENDDO
		!$OMP END PARALLEL

		TrueResNorm=CalculateVectorNorm(solver,w)
		res=TrueResNorm/sPb
		IF (res>=solver%params%Tol) THEN
		       solver%info%stat=NON_TRUE_CONVERGED
	       ELSE
		       solver%info%stat=CONVERGED

		ENDIF
    ENDFUNCTION
	
	SUBROUTINE  ConstructCurrentSolution(solver,jH,solution)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		INTEGER,INTENT(IN)::jH
		COMPLEX(REALPARM),TARGET,INTENT(INOUT)::solution(:)
		COMPLEX(REALPARM),POINTER::Hessenberg(:,:)
		COMPLEX(REALPARM),POINTER::x(:),r0(:),y(:)
		INTEGER::N,M
		INTEGER::I,J
                INTEGER::L
		N=solver%Nloc
		M=solver%BasisSize
		Hessenberg=>solver%Hessenberg

		x(1:N)=>solver%x
		y(1:M)=>solver%y
		r0(1:N)=>solver%r0


		L=jH-1
                        

		CALL ZCOPY(L,Hessenberg(:,M+1),ONE,y,ONE)
		CALL ZTRSV('U','N','N',L,Hessenberg,M+1,y,ONE)

		CALL CalculateSolutionProjectionToKrylovSubspace(solver,L)

		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
		!$OMP DO SCHEDULE(GUIDED)
			DO I=1,N
				x(I)=r0(I)+solution(I)
				solution(I)=x(I)
			ENDDO
		!$OMD ENDDO
		!$OMP END PARALLEL

	ENDSUBROUTINE


	SUBROUTINE CalculateSolutionProjectionToKrylovSubspace(solver,jH)
		TYPE(FGMRES_DATA),INTENT(IN)::solver
		INTEGER,INTENT(IN)::jH
		COMPLEX(REALPARM),POINTER::basis(:,:)
		COMPLEX(REALPARM),POINTER::x(:),y(:),r0(:)
		INTEGER::N,M
		N=solver%Nloc
		M=solver%BasisSize
		x(1:N)=>solver%x
		r0(1:N)=>solver%r0
		y(1:M)=>solver%y
		IF (solver%flex .EQV.   FLEXIBLE) THEN
			basis=>solver%PrecondKrylovBasis
		ELSE
			basis=>solver%KrylovBasis
		ENDIF

		CALL ZGEMV('N',N,jH,C_ONE,basis,N,y,ONE,C_ZERO,x,ONE)
		IF (solver%flex .EQV. FLEXIBLE) THEN
		    CALL ZCOPY(solver%Nloc,x,ONE,r0,ONE)
		ELSE
		    CALL ApplyRightPreconditioner(solver,x,r0,jH)
		ENDIF

	ENDSUBROUTINE





RECURSIVE	SUBROUTINE ApplyRightPreconditioner(solver,v_in,v_out,K)
		TYPE(FGMRES_DATA),TARGET,INTENT(IN)::solver
		COMPLEX(REALPARM),TARGET,INTENT(IN)::v_in(:)
                COMPLEX(REALPARM),TARGET,INTENT(INOUT)::v_out(:)
		INTEGER,INTENT(IN)::K
		COMPLEX(REALPARM),POINTER::pH(:)
		IF (ASSOCIATED(solver%PrecondRight)) THEN
                        CALL solver%PrecondRight(solver%MRight,v_in,v_out)
                ELSE
                        CALL ZCOPY(solver%Nloc,v_in,ONE,v_out,ONE)
                ENDIF
                IF  (solver%flex .EQV. FLEXIBLE) THEN
                        pH=>solver%PrecondKrylovBasis(:,K)
                        CALL   ZCOPY(solver%Nloc,v_out,ONE,pH,ONE)
                ENDIF
        ENDSUBROUTINE

	SUBROUTINE ApplyOperator(solver,v_in,v_out)
		TYPE(FGMRES_DATA),TARGET,INTENT(IN)::solver
		COMPLEX(REALPARM),INTENT(IN)::v_in(:)
		COMPLEX(REALPARM),INTENT(INOUT)::v_out(:)
		CALL solver%ApplyOperator(solver%A,v_in,v_out)
	ENDSUBROUTINE

	SUBROUTINE ApplyLeftPreconditioner(solver,v_in,v_out)
		TYPE(FGMRES_DATA),TARGET,INTENT(IN)::solver
		COMPLEX(REALPARM),INTENT(IN)::v_in(:)
		COMPLEX(REALPARM),INTENT(INOUT)::v_out(:)
		IF (ASSOCIATED(solver%PrecondLeft)) THEN
			CALL solver%PrecondLeft(solver%MLeft,v_in,v_out)
		ELSE
		    CALL ZCOPY(solver%Nloc,v_in,ONE,v_out,ONE)
		ENDIF
	ENDSUBROUTINE

        SUBROUTINE CalculateDotProduct(solver,m,v,K,res)
		TYPE(FGMRES_DATA),TARGET,INTENT(IN)::solver
		COMPLEX(REALPARM),INTENT(IN)::m(:,:)
		COMPLEX(REALPARM),INTENT(IN)::v(:)
		COMPLEX(REALPARM),INTENT(INOUT)::res(:)
                INTEGER,INTENT(IN)::K
		CALL solver%MANYDP(m,v,K,res,solver%dpinstance)
        ENDSUBROUTINE
                 


	FUNCTION CalculateVectorNorm(solver,v) RESULT(res)
		TYPE(FGMRES_DATA),TARGET,INTENT(IN)::solver
		COMPLEX(REALPARM),TARGET,INTENT(IN)::v(:)
		COMPLEX(REALPARM),POINTER::m(:,:),pv(:)
                TYPE(C_PTR)::cptr
		COMPLEX(REALPARM),TARGET::tt(1)
		REAL(REALPARM)::res
                pv=>v
                cptr=C_LOC(pv(1))
                CALL C_F_POINTER(cptr,m,(/SIZE(v),1/))
		CALL solver%MANYDP(m,v,ONE,tt,solver%dpinstance)
		res=SQRT(REAL(tt(1)))
	ENDFUNCTION

	FUNCTION Inform(solver) RESULT(res) 
	     TYPE(FGMRES_DATA),TARGET,INTENT(IN)::solver
	     INTEGER::res
	     res=NOINTERUPT
	     IF (ASSOCIATED(solver%Inform)) THEN
		     res=solver%Inform(solver%info)
	     ENDIF
     ENDFUNCTION




	SUBROUTINE DeallocateFGMRES(solver)
		TYPE(FGMRES_DATA),INTENT(INOUT)::solver
		CALL Free(solver%x)
		CALL Free(solver%y)
		CALL Free(solver%w)
		CALL Free(solver%r0)
		CALL Free(solver%KrylovBasis)
		CALL Free(solver%PrecondKrylovBasis)
		CALL Free(solver%Hessenberg)
		CALL Free(solver%GivensRotCos)
		CALL Free(solver%GivensRotSin)
	ENDSUBROUTINE
!-----------------------------------------------------------------------------------------------------------------------------!


      END MODULE
