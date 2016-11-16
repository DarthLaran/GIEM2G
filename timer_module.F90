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

MODULE TIMER_MODULE
	USE CONST_MODULE 	
	PRIVATE
	REAL(DOUBLEPARM)::StartTime
	INTEGER,PARAMETER::KIND8=SELECTED_INT_KIND(15)
	INTEGER(KIND8)::StartCount
	INTEGER(KIND8)::CountRate
	PUBLIC:: InitTimer, GetTime,StartTime
	CONTAINS
#ifndef MPI_TIMER
	SUBROUTINE InitTimer
		REAL(DOUBLEPARM)::t1,t2
		CALL SYSTEM_CLOCK(StartCount, CountRate)
		t1=StartCount
		t2=CountRate
		StartTime=t1/t2
	ENDSUBROUTINE
#else
	SUBROUTINE InitTimer
		StartTime=MPI_Wtime()
	ENDSUBROUTINE
#endif
#ifndef MPI_TIMER
	FUNCTION GetTime() RESULT(t)
		INTEGER(KIND8)::Ticks
		REAL(DOUBLEPARM)::t1
		CALL SYSTEM_CLOCK(Ticks)
		t1=(Ticks-StartCount)
		t=t1/CountRate
	ENDFUNCTION
#else
	FUNCTION GetTime() RESULT(t)
		USE MPI_MODULE
		REAL(DOUBLEPARM)::t
		t=MPI_Wtime()
	ENDFUNCTION
#endif
END MODULE 
