PROGRAM generate_scripts
    IMPLICIT NONE
    INTEGER, PARAMETER      :: L = 16
    INTEGER, PARAMETER      :: tau_Q(3) = [8, 32, 128]

    REAL(4), PARAMETER      :: beta_c = 0.454165
    REAL(4), PARAMETER      :: beta_initial = 0.3
    REAL(4), PARAMETER      :: beta_final = 0.5
    REAL(4), PARAMETER      :: step_betas = 0.01
    INTEGER, PARAMETER      :: total_betas = 5

    CHARACTER(*), PARAMETER :: folder = '/home/elias/Git/3dXYO2e/'
    CHARACTER(*), PARAMETER :: equilibrium_proc = "CALL m3DXYO2_standard"
    CHARACTER(*), PARAMETER :: quench_proc = "CALL quench_algorithm"
    CHARACTER(*), PARAMETER :: gen_betas_by_step = 'generate_betas_by_step_value'
    CHARACTER(*), PARAMETER :: gen_betas_by_total = 'generate_betas_by_total_value'

    CHARACTER(*), PARAMETER :: main_subroutine = equilibrium_proc
    CHARACTER(*), PARAMETER :: gen_betas_by = gen_betas_by_total

    INTEGER, PARAMETER      :: iparam=11, oparam=12, ifcomp=13, ofcomp=14, imain=15, omain=16

    ! Variables
    REAL(4), ALLOCATABLE :: betas_lat(:)
    INTEGER              :: Nbetas

    !CALL betas_array(gen_betas_by)
    CALL change_param_equil(2)
    
    CONTAINS
SUBROUTINE betas_array(g_b_b)
    INTEGER :: i, b_tot
    CHARACTER(*), INTENT(IN) :: g_b_b
    REAL(4) :: aux, b_step

    IF (g_b_b == gen_betas_by_step) THEN
        b_tot = (beta_final - beta_initial) / step_betas + 1
        ALLOCATE(betas_lat(b_tot))
        betas_lat(1) = beta_initial
        aux = betas_lat(1)
        DO i=2,b_tot
            aux = aux + step_betas
            betas_lat(i) = aux
        ENDDO
        WRITE(*,*) b_tot, betas_lat
    ENDIF

    IF (g_b_b == gen_betas_by_total) THEN
        ALLOCATE(betas_lat(total_betas))
        b_step = (beta_final - beta_initial) / REAL(total_betas)
        betas_lat(1) = beta_initial
        aux = betas_lat(1)
        DO i=2,total_betas
            aux = aux + b_step
            betas_lat(i) = aux 
        ENDDO
        WRITE(*,*) b_step, betas_lat
    ENDIF
END SUBROUTINE

SUBROUTINE change_param_equil(num_betas)
    INTEGER, INTENT(IN) :: num_betas
    INTEGER :: i,j, stat, pos1,pos2, Nlines, Nnlines
    CHARACTER(100) :: param_path, new_param_path
    CHARACTER(135), DIMENSION(:), ALLOCATABLE :: lineas
    CHARACTER(135), DIMENSION(:), ALLOCATABLE :: new_lines
    
    param_path = folder//'progm/parameters.f90'
    ! Lee todas las líneas del archivo de entrada
    CALL read_file(lineas, param_path, iparam)
    Nlines = SIZE(lineas)
    pos1 = 0; pos2 = 0
    DO i=1,Nlines
        pos1 = INDEX(lineas(i), 'PARAMETER :: L =')
        IF (pos1 /= 0) THEN
            lineas(i) = '	INTEGER,	    PARAMETER :: L = 16'
        ENDIF

        pos2 = INDEX(lineas(i), 'betas(')
        IF (pos2 /= 0) THEN
            lineas(i) = '    REAL(8),        PARAMETER :: betas(1) = [0.5D0&'
            ALLOCATE(new_lines(Nlines+num_betas+1))
            Nnlines = SIZE(new_lines)
            new_lines(1:i) = lineas(1:i)
            DO j=1,num_betas
                new_lines(i+j) = '                      ,D0&'
            ENDDO
            new_lines(i+num_betas+1) = '                        ]'
            new_lines(i+num_betas+2:Nnlines) = lineas(i+1:Nlines)
        ENDIF
    ENDDO

    IF (ALLOCATED(new_lines)) THEN
        new_param_path = folder//'progm/parameters0.f90'
        OPEN(UNIT=oparam, FILE=TRIM(new_param_path))
        DO i=1,Nnlines
            WRITE(oparam,'(A)') new_lines(i)
        ENDDO
        DEALLOCATE(new_lines)
        CLOSE(oparam)
    ENDIF

END SUBROUTINE change_param_equil

! Lee un archivo de caracteres
SUBROUTINE read_file(lines0,fname,unit0)
	CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: lines0(:)
	CHARACTER(*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: unit0
	INTEGER :: i, stat
	CHARACTER(135) :: buf

	OPEN(UNIT=unit0, FILE=TRIM(fname), STATUS='old')
	DO
		READ(unit0,'(A)',IOSTAT=stat) buf
		IF (stat/=0) EXIT
		CALL AddStrToArray(lines0, buf)
	ENDDO
    CLOSE(iparam)
END SUBROUTINE read_file

! Subrrutina para arrays dinámicos con elementos character
SUBROUTINE AddStrToArray(lines,line)
	INTEGER :: i, isize
	CHARACTER(*), INTENT(IN) :: line
	CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: lines(:)
	CHARACTER(135), ALLOCATABLE :: elist(:)

	IF (ALLOCATED(lines)) THEN
		isize = SIZE(lines)
		ALLOCATE(elist(isize+1))
		elist(1:isize) = lines
		elist(isize+1) = line

		DEALLOCATE(lines)
		CALL MOVE_ALLOC(elist, lines)
	
	ELSE
		ALLOCATE(lines(1))
		lines(1) = line
	ENDIF
END SUBROUTINE AddStrToArray

END PROGRAM generate_scripts