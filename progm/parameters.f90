MODULE PARAM
	IMPLICIT NONE
	INTEGER,	    PARAMETER :: L = 64						! Lattice size
	INTEGER,	    PARAMETER :: tauQ = 8                   ! Rate of quench
    ! value of betas
    REAL(8),        PARAMETER :: betas(1) = [0.454165D0]
    ! value of deltas
    REAL(8),        PARAMETER :: deltas(1) = [2.3D0]

	INTEGER,	    PARAMETER :: Ntherm=1D4				! Num of swweps for thermalization
	INTEGER,	    PARAMETER :: Nskip=50					! Num of sweeps between each measurmt
	INTEGER,	    PARAMETER :: Nmeas=1D3					! Num of measurements

	INTEGER,	    PARAMETER :: L2=L**2
	INTEGER,	    PARAMETER :: L3=L**3
	REAL(8),	    PARAMETER :: invL=1.0d0/dble(L)
	
	REAL(8),	    PARAMETER :: PI = ACOS(-1.0D0)
	REAL(8),	    PARAMETER :: dPI = 2.D0*ACOS(-1.0D0)
	REAL(8), 	    PARAMETER :: eps = 1.D-50
	REAL(8), 	    PARAMETER :: deltaM = 1.D0

	INTEGER,	    PARAMETER :: rvals=11, wvals=12, ifile=13, ofile=14, wmeas=16, info0=17, info1=18
    INTEGER,	    PARAMETER :: wmeas1=19, wmeas2=20, wmeas3=21

	CHARACTER(*),   PARAMETER :: folder_main_pc = '/home/elias/Git/3dXYO2e/'		! path of code
	CHARACTER(*),   PARAMETER :: folder_main_remote = '/home/eliaspolanco/2dxy/'	! path of program in ICN clustr

	! Choose a path (folder_main_pc or folder_main_remote)
	CHARACTER(*),   PARAMETER :: folder = folder_main_pc 

	SAVE
END MODULE PARAM