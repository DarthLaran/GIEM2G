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

MODULE LOGGER_MODULE
USE CONST_MODULE
USE MPI_MODULE
USE TIMER_MODULE

IMPLICIT NONE
CHARACTER(LEN=*),PARAMETER  :: CALC_TIME_FMT='(A, 1ES12.2E3, A3)'
CHARACTER(LEN=*),PARAMETER  :: CALC_TIME_RNG_FMT='(A, 1ES10.2E2, A, ES10.2E2 ,A)'
CHARACTER(LEN=*),PARAMETER  :: CALC_NUM_FMT='(A, I10)'
CHARACTER(LEN=*),PARAMETER  :: FFT_INFO_FMT = '(A36, ES10.2E3, A3, F10.3, A3)'


!--------------------------CALC IE TENSOR--------------------!
CHARACTER(LEN=*),PARAMETER  :: CALC_IE_FULL_TIME = 'Full time:   '
CHARACTER(LEN=*),PARAMETER  :: CALC_IE_WEIGHTS   = 'Weights:     '	
CHARACTER(LEN=*),PARAMETER  :: CALC_IE_LINEAR    = 'Linear part: '
CHARACTER(LEN=*),PARAMETER  :: CALC_IE_SQ        = 'Square part: '
CHARACTER(LEN=*),PARAMETER  :: CALC_IE_OTHER     = 'Other:       '
!-------------------------- FFT --------------------!
!-------------------------- IE SOLVER --------------------!

PROCEDURE(LOGGER_SUBROUTINE),POINTER::LOGGER
PROCEDURE(LOGGER_SUBROUTINE),POINTER::C_LOGGER
INTEGER(MPI_CTL_KIND)::LOGGER_MASTER


INTERFACE 
    SUBROUTINE LOGGER_SUBROUTINE(MESSAGE) 
                USE, INTRINSIC :: iso_c_binding
		CHARACTER(LEN=*,KIND=C_CHAR),INTENT(IN):: MESSAGE
    END SUBROUTINE
END INTERFACE

INTERFACE PRINT_CALC_NUMBER
	MODULE PROCEDURE PRINT_CALC_NUMBER_INT32,PRINT_CALC_NUMBER_INT64
END INTERFACE

CONTAINS
SUBROUTINE NATIVE_LOGGER(MESSAGE)
		CHARACTER(LEN=*,KIND=C_CHAR),INTENT(IN):: MESSAGE
		REAL(DOUBLEPARM)::t1,t2
		t1=GetTime()
		IF (LOGGER_MASTER==0) PRINT '(F12.5, A, A)',t1, " s: ",TRIM(MESSAGE)
END SUBROUTINE

SUBROUTINE EXTREME_LOGGER(MESSAGE)
		CHARACTER(LEN=*,KIND=C_CHAR),INTENT(IN):: MESSAGE
                CALL C_LOGGER(TRIM(message)//C_NULL_CHAR)
END SUBROUTINE

SUBROUTINE PRINT_CALC_TIME(text,time)
		CHARACTER(LEN=*),INTENT(IN):: text
                REAL(DOUBLEPARM),INTENT(IN)::time
		CHARACTER(LEN=2048)::message 
		WRITE (message,CALC_TIME_FMT) text, time,' s'
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE

SUBROUTINE PRINT_CALC_TIME_RANGE(text,time1,time2)
		CHARACTER(LEN=*),INTENT(IN):: text
                REAL(DOUBLEPARM),INTENT(IN)::time1,time2
		CHARACTER(LEN=2048)::message 
		WRITE (message,CALC_TIME_RNG_FMT) text, time1,' to ', time2,' s'
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE

SUBROUTINE PRINT_STORAGE_SIZE(text,s)
		CHARACTER(LEN=*),INTENT(IN):: text
                REAL(DOUBLEPARM),INTENT(IN)::s
		CHARACTER(LEN=2048)::message 
		WRITE (message,'(A, 1ES10.2E2, A)')  TRIM(text) , s, ' GB per	process'
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE

SUBROUTINE PRINT_CALC_RANGE(text,time1,time2,text2)
		CHARACTER(LEN=*),INTENT(IN):: text,text2
                REAL(DOUBLEPARM),INTENT(IN)::time1,time2
		CHARACTER(LEN=2048)::message 
		WRITE (message,CALC_TIME_RNG_FMT) text, time1,' to ',	time2,text2
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE
SUBROUTINE PRINT_CALC_NUMBER_INT32(text,n)
		CHARACTER(LEN=*),INTENT(IN):: text
                INTEGER(SHORT_INT),INTENT(IN)::n
		CHARACTER(LEN=2048,KIND=C_CHAR)::message 
		WRITE (message,CALC_NUM_FMT) text, n
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE
SUBROUTINE PRINT_CALC_NUMBER_INT64(text,n)
		CHARACTER(LEN=*),INTENT(IN):: text
                INTEGER(LONG_INT),INTENT(IN)::n
		CHARACTER(LEN=2048,KIND=C_CHAR)::message 
		WRITE (message,CALC_NUM_FMT) text, n
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE
SUBROUTINE PRINT_BORDER
                CHARACTER(LEN=*),PARAMETER  :: BORDER_STR ='******************************&
                &******************************************************************'
                CHARACTER(LEN=*),PARAMETER  :: BORDER_FMT ='(A80)'
		CHARACTER(LEN=80,KIND=C_CHAR)::message 
		WRITE (message,BORDER_FMT) BORDER_STR
               CALL LOGGER(TRIM(message)) 
END SUBROUTINE
END MODULE
