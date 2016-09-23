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

MODULE Data_Types_Module
	USE CONST_MODULE
	IMPLICIT NONE
	PRIVATE
	TYPE BKG_DATA_TYPE
		INTEGER :: Nl
		REAL (REALPARM),POINTER:: sigma(:), thick(:), depth(:)
		COMPLEX(REALPARM),POINTER::csigma(:)
		REAL (REALPARM) :: freq,omega
		COMPLEX (RealParm),POINTER :: k(:),k2(:)
		COMPLEX(REALPARM)::iwm
	ENDTYPE

	TYPE ANOMALY_TYPE
		INTEGER :: Nx,Ny,Nz
		INTEGER::Ny_loc
		REAL(REALPARM)::dx,dy
		REAL(REALPARM),POINTER::z(:)=>NULL()
		REAL(REALPARM),POINTER::siga(:,:,:)=>NULL()
		REAL(REALPARM),POINTER::epsa(:,:,:)=>NULL()
		REAL(REALPARM),POINTER::dz(:)=>NULL()
		INTEGER,POINTER::Lnumber(:)=>NULL()
	ENDTYPE

	TYPE RECEIVER_TYPE
		REAL(REALPARM)::zrecv,x_shift,y_shift
		INTEGER::recv_layer
		INTEGER:: anom_cell !anom_cell<1 recv is upper anomaly , anom_cell >Nz recv is lower anomaly
	END TYPE

	TYPE FGMRES_CTL_TYPE
		REAL(REALPARM)::misfit
		INTEGER::fgmres_maxit
		INTEGER::gmres_buf
		INTEGER::fgmres_buf
		INTEGER::ort_type
	END TYPE

	INTERFACE GetLayer
		MODULE PROCEDURE GetLayerSingle,GetLayerArray
	END INTERFACE

	PUBLIC:: BKG_DATA_TYPE,ANOMALY_TYPE, RECEIVER_TYPE,FGMRES_CTL_TYPE

	PUBLIC::GetLayer,SET_FREQ
	PUBLIC::PrepareRecvs, AllocateSiga
	CONTAINS
	FUNCTION GetLayerArray(z,bkg) RESULT(lnum)
		REAL (KIND = RealParm), INTENT(IN) ::z(:)
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		INTEGER ::lnum(SIZE(z))
		INTEGER ::I
		DO I=1,SIZE(z)
			lnum(I)=GetLayerSingle(z(I),bkg)
		ENDDO
	END FUNCTION
	FUNCTION GetLayerSingle(z,bkg) RESULT(lnum)
		REAL (KIND = RealParm), INTENT(IN) ::z
		TYPE (BKG_DATA_TYPE),INTENT(IN)::bkg
		INTEGER ::lnum
		IF (ABS(z)<1d-15) THEN
			lnum=0
			RETURN
		ENDIF
		
		IF (z<0) THEN
			lnum=-1
			RETURN
		ENDIF
		
		IF (bkg%Nl==1) THEN
			lnum=1
			RETURN
		ENDIF
			
		lnum=1
		IF (z>bkg%depth(bkg%Nl-1)) THEN
			lnum=bkg%Nl
			RETURN
		ELSE
			DO WHILE (z>bkg%depth(lnum))
				lnum=lnum+1
			ENDDO
		ENDIF
	END FUNCTION

	SUBROUTINE SET_FREQ(bkg,f)
		TYPE (BKG_DATA_TYPE),INTENT(INOUT)::bkg
		REAL(KIND=RealParm),INTENT(IN)::f
		bkg%freq=f
		bkg%omega=f*2D0*PI
#ifndef NO_DISPLACEMENT_CURRENTS
		bkg%csigma(0)=AIR_CONDUCTIVITY-C_IONE*bkg%omega*EPS0
		bkg%csigma(1:)=bkg%sigma-C_IONE*bkg%omega*EPS0
#else
		bkg%csigma(0)=0!AIR_CONDUCTIVITY
		bkg%csigma(1:)=bkg%sigma
#endif
		bkg%k2=C_IONE*bkg%csigma*MU0*bkg%omega
		bkg%k=SQRT(bkg%k2)
		bkg%iwm=C_IONE*MU0*bkg%omega
	END SUBROUTINE



	SUBROUTINE PrepareRecvs(recvs,anomaly,bkg)
			TYPE (RECEIVER_TYPE),POINTER,INTENT(INOUT)::recvs(:)
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			TYPE (BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
			INTEGER::I,Nr,l,Iz
			REAL(REALPARM)::zrecv
			Nr=SIZE(recvs)
			DO I=1,Nr
				zrecv=recvs(I)%zrecv
				recvs(I)%recv_layer=MAX(GetLayer(zrecv,bkg),1)
				IF (zrecv<=anomaly%z(0)) THEN
					recvs(I)%anom_cell=0
					
				ELSEIF (zrecv>=anomaly%z(anomaly%Nz)) THEN
					recvs(I)%anom_cell=anomaly%Nz+1
				ELSE
				    DO Iz=1,anomaly%Nz
					IF (abs(zrecv/anomaly%z(Iz)-1)<1e-7) THEN
					    recvs(I)%anom_cell=Iz
					    EXIT
					ELSEIF (zrecv>anomaly%z(Iz-1) .AND. zrecv< anomaly%z(Iz)) THEN
						recvs(I)%zrecv=anomaly%z(Iz-1)
						recvs(I)%anom_cell=Iz-1
						EXIT
					ENDIF
				    ENDDO
				ENDIF
			ENDDO
	END SUBROUTINE

	SUBROUTINE PrepareRecvs2(recvs,anomaly,bkg)
			TYPE (RECEIVER_TYPE),POINTER,INTENT(INOUT)::recvs(:)
			TYPE (ANOMALY_TYPE),INTENT(IN)::anomaly
			TYPE (BKG_DATA_TYPE),TARGET,INTENT(IN)::bkg
			INTEGER::I,Nr,l
			REAL(REALPARM)::zrecv
			Nr=SIZE(recvs)
			DO I=1,Nr
				zrecv=recvs(I)%zrecv
				recvs(I)%recv_layer=MAX(GetLayer(zrecv,bkg),1)
				IF (zrecv<=anomaly%z(0)) THEN
					recvs(I)%anom_cell=0
					
				ELSEIF (zrecv>=anomaly%z(anomaly%Nz)) THEN
					recvs(I)%anom_cell=anomaly%Nz+1
				ELSE
					l=1
					DO WHILE (zrecv>anomaly%z(l))
						l=l+1
					ENDDO
					recvs(I)%anom_cell=l
				ENDIF
			ENDDO
	END SUBROUTINE
	SUBROUTINE AllocateSiga(anomaly)
		TYPE (ANOMALY_TYPE),INTENT(INOUT)::anomaly
		IF (ASSOCIATED(anomaly%siga))	DEALLOCATE(anomaly%siga)
			ALLOCATE(anomaly%siga(anomaly%Nz,anomaly%Nx,anomaly%Ny_loc))
	ENDSUBROUTINE

	ENDMODULE

