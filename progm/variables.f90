MODULE VARS
	USE PARAM
	IMPLICIT NONE

	REAL(8)                 :: S(L,L,L,2) 		 			! Red de espines
	INTEGER                 :: Sv3(L,L,L,3)					! Red de vortices

	REAL(8)			        :: beta, delta      			! Elementos del vector betas/deltas

	! Observables:
	REAL(8) 		    :: E, Mag(2), Corr(L), F, acceptance_rate
	INTEGER 		    :: num_of_clus
	! (defectos topologicos)
    REAL(8) 		    :: l_str_avr
	INTEGER 		    :: num_of_vort, num_of_av, nvx, nvy, nvz
	INTEGER 		    :: num_of_str, size_of_str, max_size_str
    ! Auxiliares para calcular observables
	INTEGER, ALLOCATABLE :: vortex_array(:,:)
	REAL(8), ALLOCATABLE :: acc_rate_array(:,:), energy_array(:,:), mag_array(:,:,:)
    REAL(8)			    :: O(L,2), cos_arr(L,L)
	REAL(8) 		    :: vortex1, vortex2, vortex3, sum_vortex
	REAL(8) 		    :: Dphi12, Dphi23, Dphi34, Dphi41, Dphi21, Dphi15
	REAL(8) 		    :: Dphi56, Dphi62, Dphi14, Dphi47, Dphi75, Dphi51

	! Variables para formatos
	CHARACTER(100)	:: FMT1, FMT2, FMT3

    ! Variables para las rutas
    CHARACTER(100)	:: data_path, fprmpt0, fprmpt1, data_path1, data_path2, data_path3

    ! Variables para medir el tiempo
    REAL(4)         :: tTherm1, tTherm2, tMeas1, tMeas2, t1, t2, days, hours, mins, secs, res_secs
    INTEGER         :: c_year, c_month, c_day, c_hour, c_min, c_sec

	! Otros tipos de variable
	TYPE :: my_list
		INTEGER :: value
		TYPE(my_list), POINTER :: next
	END TYPE	

	TYPE :: node
		INTEGER, DIMENSION(:), POINTER :: p
		TYPE (my_list), POINTER :: out0
		TYPE (my_list), POINTER :: head
		TYPE (my_list), POINTER :: tail
	END TYPE
END MODULE VARS