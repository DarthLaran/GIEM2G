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

MODULE LOCAL_OMP_FFT_MODULE
        USE CONST_MODULE
        USE FFTW3
        USE MPI_MODULE
        USE Timer_Module 
        USE DATA_TYPES_MODULE
        USE LOGGER_MODULE
        IMPLICIT NONE

        TYPE LOCAL_OMP_FFT_DATA
                TYPE(C_PTR)::plan_fwd
                TYPE(C_PTR)::plan_bwd
                COMPLEX(REALPARM),POINTER::data_in(:,:)
                COMPLEX(REALPARM),POINTER::data_out(:,:)
                INTEGER::Nc,NT
        ENDTYPE
	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_PAT=IOR(FFTW_PATIENT,FFTW_DESTROY_INPUT )
	INTEGER(KIND=KIND(FFTW_PATIENT)), PARAMETER  :: FFT_MEA=IOR(FFTW_MEASURE,FFTW_DESTROY_INPUT )
CONTAINS
        FUNCTION CALC_LOCAL_OMP_FFT_SIZE(Nc) RESULT (R)
                INTEGER,INTENT(IN)::Nc
                INTEGER::R
                INTEGER::NT
                NT=OMP_GET_NUM_THREADS()
                R=(Nc*6+4)*NT
        ENDFUNCTION

        SUBROUTINE PREPARE_LOCAL_OMP_FFT(LFFT,Nc,buff_in,buff_out)
                TYPE(LOCAL_OMP_FFT_DATA),INTENT(INOUT)::LFFT
                INTEGER,INTENT(IN)::Nc
                TYPE(C_PTR),INTENT(IN)::buff_in,buff_out
                COMPLEX(REALPARM),POINTER::ptr(:)
                INTEGER::NT,N
                INTEGER::pN(1)
                INTEGER::M4,M6
                NT=OMP_GET_NUM_THREADS()
                LFFT%NT=NT
                N=(Nc*6+4)*NT
                CALL C_F_POINTER(buff_in,ptr,(/N/))
                LFFT%data_in(1:Nc*6+4,0:NT-1)=>ptr
                CALL C_F_POINTER(buff_out,ptr,(/N/))
                LFFT%data_out(1:Nc*6+4,0:NT-1)=>ptr

		CALL FFTW_PLAN_WITH_NTHREADS(1)
                pN=Nc
                M4=4
                M6=6
		LFFT%plan_fwd=fftw_plan_many_dft(ONE,pN,M4,LFFT%data_in(:,0),pN,ONE,Nc,&
                                                        &LFFT%data_out(:,0),pN,ONE,Nc&
                                                        &,FFTW_FORWARD, FFT_MEA)

		LFFT%plan_bwd=fftw_plan_many_dft(ONE,pN,M6,LFFT%data_in(:,0),pN,ONE,Nc,&
                                                        &LFFT%data_out(:,0),pN,ONE,Nc&
                                                        &,FFTW_BACKWARD, FFT_MEA)

		CALL FFTW_PLAN_WITH_NTHREADS(NT)

        ENDSUBROUTINE

        SUBROUTINE CALCULATE_FORWARD_AT_THREAD(LFFT,I)
                TYPE(LOCAL_OMP_FFT_DATA),INTENT(INOUT)::LFFT
                INTEGER,INTENT(IN)::I
		CALL fftw_execute_dft(LFFT%plan_fwd,LFFT%data_in(:,I),LFFT%data_out(:,I))
        ENDSUBROUTINE

        SUBROUTINE CALCULATE_BACKWARD_AT_THREAD(LFFT,I)
                TYPE(LOCAL_OMP_FFT_DATA),INTENT(INOUT)::LFFT
                INTEGER,INTENT(IN)::I
		CALL fftw_execute_dft(LFFT%plan_bwd,LFFT%data_in(:,I),LFFT%data_out(:,I))
        ENDSUBROUTINE
ENDMODULE
