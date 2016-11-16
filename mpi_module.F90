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

MODULE MPI_MODULE
#ifdef LEGACY_MPI
	INCLUDE 'mpif.h'
#else
	USE MPI
#endif
INTEGER,PARAMETER::MPI_CTL_KIND=KIND(MPI_COMM_WORLD)
END MODULE
