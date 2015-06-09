MODULE IntegralCodes
	INTEGER, PUBLIC, PARAMETER::NFIRST=424/2!0
	INTEGER, PUBLIC, PARAMETER::NLAST=1424/2!2047
!Integrals by both points quadruple integrals!
        INTEGER, PARAMETER :: INT4=1
        INTEGER, PARAMETER :: INT4DXX=2

        INTEGER, PARAMETER :: INT4DXY=4
        INTEGER, PARAMETER :: INT4DXY2=5

        INTEGER, PARAMETER :: INT4DX=6
        INTEGER, PARAMETER :: INT4DX2=7

        INTEGER, PARAMETER :: INT4DYX=INT4DXY

        INTEGER, PARAMETER :: INT4DYY=8
        INTEGER, PARAMETER :: INT4DYY2=9

        INTEGER, PARAMETER :: INT4DY=10

        INTEGER, PARAMETER :: INT4LAST=INT4DY

!Integrals by one  point double integrals!

        INTEGER, PARAMETER :: INT2=21
        INTEGER, PARAMETER :: INT2DXX=22

        INTEGER, PARAMETER :: INT2DXY=24

        INTEGER, PARAMETER :: INT2DX=26

        INTEGER, PARAMETER :: INT2DYX=INT2DXY

        INTEGER, PARAMETER :: INT2DYY=28

        INTEGER, PARAMETER :: INT2DY=30
        INTEGER, PARAMETER :: INT2LAST=INT2DY
!-----------------------For sources------------------------------------!
!----------------------integrals along axis. --------------------------------------------------!
        INTEGER, PARAMETER :: INT1=41
        INTEGER, PARAMETER :: INT1DXX=42

        INTEGER, PARAMETER :: INT1DXY=44

        INTEGER, PARAMETER :: INT1DX=46

        INTEGER, PARAMETER :: INT1DYX=INT4DXY

        INTEGER, PARAMETER :: INT1DYY=48

        INTEGER, PARAMETER :: INT1DY=50
        INTEGER, PARAMETER :: INT1LAST=INT1DY
!---------------------------- any line---------!!!        
        INTEGER, PARAMETER :: INT1A=61
        INTEGER, PARAMETER :: INT1ADXX=62

        INTEGER, PARAMETER :: INT1ADXY=64

        INTEGER, PARAMETER :: INT1ADX=66

        INTEGER, PARAMETER :: INT1ADYX=INT4DXY

        INTEGER, PARAMETER :: INT1ADYY=68

        INTEGER, PARAMETER :: INT1ADY=70
        INTEGER, PARAMETER :: INT1ALAST=INT1ADY
!---- along axis by one point and at the rectangle by another one ---- !
        INTEGER, PARAMETER :: INT3=81
        INTEGER, PARAMETER :: INT3DXX=82

        INTEGER, PARAMETER :: INT3DXY=84

        INTEGER, PARAMETER :: INT3DX=86

        INTEGER, PARAMETER :: INT3DYX=INT4DXY

        INTEGER, PARAMETER :: INT3DYY=88

        INTEGER, PARAMETER :: INT3DY=90
        INTEGER, PARAMETER :: INT3LAST=INT3DY
!---------------------------- any line and rectangle ---------!!!        
        INTEGER, PARAMETER :: INT3A=101
        INTEGER, PARAMETER :: INT3ADXX=102

        INTEGER, PARAMETER :: INT3ADXY=104

        INTEGER, PARAMETER :: INT3ADX=106

        INTEGER, PARAMETER :: INT3ADYX=INT4DXY

        INTEGER, PARAMETER :: INT3ADYY=108

        INTEGER, PARAMETER :: INT3ADY=110
        INTEGER, PARAMETER :: INT3ALAST=INT3ADY
!-------------------------------  Interfaces to functions for computing integrals    ----------------------!!!!!!
END MODULE

