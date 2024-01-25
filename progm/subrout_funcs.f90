MODULE MC
	USE VARS
	IMPLICIT NONE

	INTEGER			:: tags(L,L,L), parents(L3)
	INTEGER			:: d6(6)

	REAL(8)			:: rW(2), DH, bprob, S_prime(2), cosd
	REAL(8)			:: cos_phi, dot_product_OO

	REAL(8)			:: delt_phi_S(6), delt_phi_P(6), delt_phi_M(6), upper, lower, new_angle

	TYPE(node)		:: row(L3)
	TYPE (my_list), POINTER :: out1, out2
	TYPE (my_list), POINTER :: head1, head2
	TYPE (my_list), POINTER :: tail1, tail2

! ------------------------ SUBROUTINES ---------------------------
	CONTAINS
SUBROUTINE m3DXYO2_standard
    INTEGER :: i,j,k

    CALL init_random_seed
    CALL init_arrays

    CALL FechaYHora
    fprmpt0 = folder//'prompt/L'//TRIM(str(L))//'b_'//TRIM(str(c_year))//'_'&
        //TRIM(str(c_month))//'_'//TRIM(str(c_day))//'_'//TRIM(str(c_hour))//'_'&
        //TRIM(str(c_min))//'_'//TRIM(str(c_sec))//'.txt'
    OPEN(info0, FILE=TRIM(fprmpt0))
    WRITE(info0,"('standard action')")
    WRITE(info0,"('L  = ',I3)") L
    FMT2 = "(A8,"//TRIM(str(SIZE(betas)))//"(F17.6,','),/)"
    WRITE(info0,FMT2)'betas:',betas
    WRITE(info0,"('Ntherm = ', I6, /'Nmeas  = ', I6, /'Nskip  = ', I6/)")Ntherm,Nmeas,Nskip
    CLOSE(info0)

    CALL CPU_TIME(t1)
    DO i=1,SIZE(betas)
        beta = betas(i)
        WRITE(*,*) beta
        fprmpt1 = folder//'prompt/L'//TRIM(str(L))//'b'//TRIM(dble2str(beta))//'.txt'
        OPEN(info1, FILE=TRIM(fprmpt1))
        CALL FechaYHora
        WRITE(info1, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
        c_year, c_month, c_day, c_hour, c_min, c_sec
        WRITE(info1,'("beta: ", F17.6)') beta
        CLOSE(info1)

        CALL CPU_TIME(tTherm1)
        CALL hot_start

        DO j=1,Ntherm
            !CALL multicluster
            CALL metropolis
            !CALL metropolis_random
            !CALL Glauber
            !CALL heat_bath
            !CALL heat_bath_random
        ENDDO
        CALL CPU_TIME(tTherm2)

        CALL CPU_TIME(tMeas1)
        data_path = folder//'Data/L'//TRIM(str(L))//'b'//TRIM(dble2str(beta))//".txt"
        WRITE(*,*) data_path
        OPEN(UNIT=wmeas, FILE=TRIM(data_path))
        WRITE(wmeas,'(2A25,26X,A26,1X,7A11,2A22)') "Energy", "Mag", "F", "N clusters",&
            "N vorts", "nvx", "nvy", "nvz", "N strings", "Max s str", "<l string>", "Correl"

        DO j=1,Nskip*Nmeas
            !CALL multicluster
            CALL metropolis
            !CALL metropolis_random
            !CALL Glauber
            !CALL heat_bath
            !CALL heat_bath_random
            IF ( MOD(j,Nskip) == 0 ) THEN
                CALL take_meas
            ENDIF
        ENDDO

        WRITE(wmeas,'(2A25,26X,A26,1X,7A11,2A22)') "Energy", "Mag", "F", "N clusters",&
            "N vorts", "nvx", "nvy", "nvz", "N strings", "Max s str", "<l string>", "Correl"
        CALL CPU_TIME(tMeas2)
        WRITE(wmeas,*)
        WRITE(wmeas,"(F0.6,' hours of thermalization')") (tTherm2 - tTherm1) / 3600.
        WRITE(wmeas,"(F0.6,' hours of measurements')") (tMeas2 - tMeas1) / 3600.
        CLOSE(wmeas)
        OPEN(info1, FILE=TRIM(fprmpt1), STATUS='unknown',POSITION="append", action = "write")
        WRITE(info1,"(/F0.6,' hours of simulation')") (tMeas2 - tTherm1) / 3600.
        CLOSE(info1)
        
    ENDDO
    CALL CPU_TIME(t2)

    secs = t2 - t1
    days = secs / 86400.; res_secs = secs - FLOOR(days) * 86400.
    hours = res_secs / 3600.; res_secs = res_secs - FLOOR(hours) * 3600.
    mins = res_secs / 60.; res_secs = res_secs - FLOOR(mins) * 60.
    
    OPEN(info0, FILE=TRIM(fprmpt0), STATUS='unknown',POSITION="append", action = "write")
    30 FORMAT('Computing time: ',I4,' days ', I3,' hours ',I3,' minutes ',F8.4,' seconds.')
    WRITE(info0,30) FLOOR(days), FLOOR(hours), FLOOR(mins), res_secs
    WRITE(info0,"('time in seconds: ', F25.4 )") secs
    CLOSE(info0)
END SUBROUTINE m3DXYO2_standard

SUBROUTINE quench_algorithm
    INTEGER :: i,j,k
    CHARACTER(10) :: Qname = 'rm'

    CALL FechaYHora
    fprmpt1 = folder//'prompt/L'//TRIM(str(L))//'q'//TRIM(str(tauQ))//'_'//TRIM(str(c_year))//'_'&
        //TRIM(str(c_month))//'_'//TRIM(str(c_day))//'_'//TRIM(str(c_hour))//'_'&
        //TRIM(str(c_min))//'_'//TRIM(str(c_sec))//'.txt'
    OPEN(info1, FILE=TRIM(fprmpt1))
    WRITE(info1,"('quench algorithm')")
    WRITE(info1,"('L  = ',I3)") L
    WRITE(info1,"('tau_Q  = ',I3)") tauQ
    FMT2 = "(A8,"//TRIM(str(SIZE(betas)))//"(F13.6,','),/)"
    WRITE(info1,FMT2)'betas:',betas
    WRITE(info1,"('Ntherm = ', I6, /'Nmeas  = ', I6, /)")Ntherm,Nmeas
    CLOSE(info1)

    ALLOCATE(vortex_array(Nmeas, SIZE(betas)))
    ALLOCATE(energy_array(Nmeas, SIZE(betas)))
    ALLOCATE(mag_array(Nmeas, SIZE(betas),2))
    ALLOCATE(acc_rate_array(Nmeas, SIZE(betas)))
    vortex_array=0.D0
    energy_array = 0.D0
    mag_array = 0.D0
    acc_rate_array = 0.D0

    CALL CPU_TIME(tMeas1)
    DO i=1,Nmeas
        IF (MOD(i,50)==0) THEN
            WRITE(*,*) i
        ENDIF

        CALL hot_start

        beta = betas(1)
        DO j=1,Ntherm
            !CALL multicluster
            CALL metropolis_random
            !CALL Glauber
            !CALL heat_bath
            !CALL heat_bath_random
        ENDDO

        DO j=1,SIZE(betas)
            beta=betas(j)
            
            !CALL multicluster
            CALL metropolis_random
            !CALL Glauber
            !CALL heat_bath_random
            !acc_rate_array(k,1) = acceptance_rate

            CALL init_vars
            CALL calc_vortices
            CALL take_meas_q
            vortex_array(i,j) = num_of_vort
            acc_rate_array(i,j) = acceptance_rate
            energy_array(i,j) = E
            mag_array(i,j,:) = Mag
        ENDDO

    ENDDO

    data_path = folder//'Data/L'//trim(str(L))//'tauQ'//TRIM(str(tauQ))//'v_'//TRIM(Qname)//'.txt'
    OPEN(UNIT=wmeas, FILE=TRIM(data_path))
    DO i=1,Nmeas
        WRITE(wmeas,*) vortex_array(i,:)
    ENDDO
    CALL CPU_TIME(tMeas2)
    WRITE(wmeas,*)
    WRITE(wmeas,*) (tMeas2 - tMeas1 ) / 3600., "hours" 
    CLOSE(wmeas)

    data_path1 = folder//'Data/L'//trim(str(L))//'tauQ'//TRIM(str(tauQ))//'E_'//TRIM(Qname)//'.txt'
    OPEN(UNIT=wmeas1, FILE=TRIM(data_path1))
    DO i=1,Nmeas
        WRITE(wmeas1,*) energy_array(i,:)
    ENDDO
    CLOSE(wmeas1)

    data_path2 = folder//'Data/L'//trim(str(L))//'tauQ'//TRIM(str(tauQ))//'Mag_'//TRIM(Qname)//'.txt'
    OPEN(UNIT=wmeas2, FILE=TRIM(data_path2))
    DO i=1,Nmeas
        WRITE(wmeas2,*) mag_array(i,:,1)
        WRITE(wmeas2,*) mag_array(i,:,2)
    ENDDO
    CLOSE(wmeas2)

    data_path3 = folder//'Data/L'//trim(str(L))//'tauQ'//TRIM(str(tauQ))//'acc_'//TRIM(Qname)//'txt'
    OPEN(UNIT=wmeas3, FILE=TRIM(data_path3))
    DO i=1,Nmeas
        WRITE(wmeas3,*) acc_rate_array(i,:)
    ENDDO
    CLOSE(wmeas3)

    DEALLOCATE(vortex_array, energy_array, mag_array, acc_rate_array)

    secs = tMeas2 - tTherm1
    days = secs / 86400.; res_secs = secs - FLOOR(days) * 86400.
    hours = res_secs / 3600.; res_secs = res_secs - FLOOR(hours) * 3600.
    mins = res_secs / 60.; res_secs = res_secs - FLOOR(mins) * 60.

    OPEN(info1, FILE=TRIM(fprmpt1), STATUS='unknown',POSITION="append", action = "write")
    23 FORMAT('/Computing time: ',I4,' days ', I3,' hours ',I3,' minutes ',F8.4,' seconds.')
    WRITE(info1,23) FLOOR(days), FLOOR(hours), FLOOR(mins), res_secs
    WRITE(info1,"('time in seconds: ', F25.4 )") secs
    CLOSE(info1)

END SUBROUTINE quench_algorithm

SUBROUTINE m3DXYO2_constrain
    INTEGER :: i,j,k

    CALL init_random_seed
    CALL init_arrays

    CALL FechaYHora
    fprmpt0 = folder//'prompt/L'//TRIM(str(L))//'d_'//TRIM(str(c_year))//'_'&
        //TRIM(str(c_month))//'_'//TRIM(str(c_day))//'_'//TRIM(str(c_hour))//'_'&
        //TRIM(str(c_min))//'_'//TRIM(str(c_sec))//'.txt'
    OPEN(info0, FILE=TRIM(fprmpt0))
    WRITE(info0,"('constrained action')")
    WRITE(info0,"('L  = ',I3)") L
    FMT2 = "(A8,"//TRIM(str(SIZE(betas)))//"(F17.6,','),/)"
    WRITE(info0,FMT2)'deltas:',deltas
    WRITE(info0,"('Ntherm = ', I6, /'Nmeas  = ', I6, /'Nskip  = ', I6/)")Ntherm,Nmeas,Nskip
    CLOSE(info0)

    DO i=1,SIZE(deltas)
        delta = deltas(i)
        WRITE(*,*) delta
        cosd = DCOS(delta)

        fprmpt1 = folder//'prompt/L'//TRIM(str(L))//'d'//TRIM(dble2str(delta))//'.txt'
        OPEN(info1, FILE=TRIM(fprmpt1))
        CALL FechaYHora
        WRITE(info1, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
        c_year, c_month, c_day, c_hour, c_min, c_sec
        WRITE(info1,'("delta: ", F17.6)') delta
        CLOSE(info1)

        CALL CPU_TIME(tTherm1)
        CALL warm_start

        DO j=1,Ntherm
            CALL constraint_algorithm
        ENDDO
        CALL CPU_TIME(tTherm2)

        CALL CPU_TIME(tMeas1)
        data_path = folder//'Data/L'//TRIM(str(L))//'d'//TRIM(dble2str(delta))//".txt"
        WRITE(*,*) data_path

        OPEN(UNIT=wmeas, FILE=TRIM(data_path))
        WRITE(wmeas,'(A25,26X,A26,1X,7A11,2A22)') "Mag", "F", "N clusters",&
            "N vorts", "nvx", "nvy", "nvz", "N strings", "Max s str", "<l string>", "Correl"

        DO j=1,Nskip*Nmeas
            CALL constraint_algorithm
            IF ( MOD(j,Nskip) == 0 ) THEN
                CALL take_meas_c
            ENDIF
        ENDDO

        WRITE(wmeas,'(A25,26X,A26,1X,7A11,2A22)') "Mag", "F", "N clusters",&
            "N vorts", "nvx", "nvy", "nvz", "N strings", "Max s str", "<l string>", "Correl"
        CALL CPU_TIME(tMeas2)
        WRITE(wmeas,*)
        WRITE(wmeas,"(F0.6,' hours of thermalization')") (tTherm2 - tTherm1) / 3600.
        WRITE(wmeas,"(F0.6,' hours of measurements')") (tMeas2 - tMeas1) / 3600.
        CLOSE(wmeas)
        OPEN(info1, FILE=TRIM(fprmpt1), STATUS='unknown',POSITION="append", action = "write")
        WRITE(info1,"(/F0.6,' hours of simulation')") (tMeas2 - tTherm1) / 3600.
        CLOSE(info1)
    ENDDO
    CALL CPU_TIME(t2)

    secs = t2 - t1
    days = secs / 86400.; res_secs = secs - FLOOR(days) * 86400.
    hours = res_secs / 3600.; res_secs = res_secs - FLOOR(hours) * 3600.
    mins = res_secs / 60.; res_secs = res_secs - FLOOR(mins) * 60.
    
    OPEN(info0, FILE=TRIM(fprmpt0), STATUS='unknown',POSITION="append", action = "write")
    31 FORMAT('Computing time: ',I4,' days ', I3,' hours ',I3,' minutes ',F8.4,' seconds.')
    WRITE(info0,31) FLOOR(days), FLOOR(hours), FLOOR(mins), res_secs
    WRITE(info0,"('time in seconds: ', F25.4 )") secs
    CLOSE(info0)
END SUBROUTINE m3DXYO2_constrain

SUBROUTINE take_meas
	INTEGER :: i,j,k
	
	CALL init_vars
	CALL calc_vortices
	!CALL calc_F

	DO i=1,L
		DO j=1,L; DO k=1,L
 
			E = E - DOT_PRODUCT(S(i,j,k,:), S(ip(i),j,k,:)) &
			- DOT_PRODUCT(S(i,j,k,:), S(i,ip(j),k,:)) - DOT_PRODUCT(S(i,j,k,:),S(i,j,ip(k),:))

			Mag = Mag + S(i,j,k,:)
		ENDDO; ENDDO
		O(i,1) = sum(S(i,:,:,1))
		O(i,2) = sum(S(i,:,:,2))

		dot_product_OO =  dot_product(O(1,:),O(i,:))
		Corr(i) = Corr(i) + dot_product_OO
	ENDDO

	FMT1='(4(F25.16,1X),7(I10,1X),F25.16,1X,'//trim(str(L))//'(F0.16,1X))'
	WRITE(wmeas,FMT1) E, Mag, F, num_of_clus, num_of_vort, nvx, nvy, nvz,&
		num_of_str, max_size_str,l_str_avr, Corr
END SUBROUTINE take_meas

SUBROUTINE take_meas_c
	INTEGER :: i,j,k
	
	CALL init_vars
	CALL calc_vortices
	!CALL calc_F

	DO i=1,L
		DO j=1,L; DO k=1,L
			Mag = Mag + S(i,j,k,:)
		ENDDO; ENDDO
		O(i,1) = sum(S(i,:,:,1))
		O(i,2) = sum(S(i,:,:,2))

		dot_product_OO =  dot_product(O(1,:),O(i,:))
		Corr(i) = Corr(i) + dot_product_OO
	ENDDO

	FMT1='(3(F25.16,1X),7(I10,1X),F25.16,1X,'//trim(str(L))//'(F0.16,1X))'
	WRITE(wmeas,FMT1) Mag, F, num_of_clus, num_of_vort, nvx, nvy, nvz,&
		num_of_str, max_size_str,l_str_avr, Corr
END SUBROUTINE take_meas_c

SUBROUTINE take_meas_q
	INTEGER :: i,j,k

	CALL init_vars
	!CALL calc_F

	DO i=1,L
		DO j=1,L; DO k=1,L
            E = E - DOT_PRODUCT(S(i,j,k,:), S(ip(i),j,k,:)) &
			- DOT_PRODUCT(S(i,j,k,:), S(i,ip(j),k,:)) - DOT_PRODUCT(S(i,j,k,:),S(i,j,ip(k),:))

			Mag = Mag + S(i,j,k,:)
		ENDDO; ENDDO
	ENDDO
END SUBROUTINE take_meas_q

! Inicio caliente (totalmente aleatorio)
SUBROUTINE hot_start
	REAL(8)::x, theta
	INTEGER::i,j,k

	CALL init_random_seed

	S=0.D0
	DO i=1,L; DO j=1,L; DO k=1,L
		CALL random_number(x)
		theta = DPI*(x-0.5D0)

		S(i,j,k,:) = [DCOS(theta), DSIN(theta)]
		!WRITE(*,'(3I3,3F25.16)') i,j,k, S(i,j,k,:), DOT_PRODUCT(S(i,j,k,:),S(i,j,k,:))
	ENDDO; ENDDO; ENDDO
	!WRITE(*,*)
END SUBROUTINE hot_start

! Inicio tibio (parcialemte aleatorio: restringido a un intervalo de ángulo)
SUBROUTINE warm_start
	REAL(8)::alpha, theta, de, r
	INTEGER::i,j,k

	CALL init_random_seed
	CALL random_number(r)

	alpha = DPI*(r-0.5D0)

	S=0.D0
	DO i=1,L; DO j=1,L; DO k=1,L

		CALL RANDOM_NUMBER(de)

		theta = alpha + delta*(de - 0.5D0)

		S(i,j,k,:) = [DCOS(theta), DSIN(theta)]
		!WRITE(*,'(3I3,1X,4F20.16)') i,j,k, S(i,j,k,:), theta, alpha!DOT_PRODUCT(S(i,j,k,:),S(i,j,k,:))
	ENDDO; ENDDO; ENDDO

	!WRITE(*,*) alpha, alpha+delta*0.5, (alpha-delta*0.5)
END SUBROUTINE warm_start

! Crea vector aleatorio
SUBROUTINE rand_vec
	REAL(8)::alpha, theta

	CALL random_number(alpha)
	theta = DPI*(alpha-0.5d0)

	rW(:) = [ DCOS(theta), DSIN(theta) ]
END SUBROUTINE rand_vec

SUBROUTINE multicluster
	INTEGER:: i,j,k
	REAL(8):: r1, r2, r3
	INTEGER :: down, right, outside

	DH = 0.d0
	tags = 0

	CALL rand_vec

	DO i=1,L; DO j=1,L; DO k=1,L
		outside=0; down=0; right=0

		! Conexión hacia abajo
		DH = 2*DOT_PRODUCT(rW,S(i,j,k,:))*DOT_PRODUCT(rW,S(i,ip(j),k,:)) 
		IF (DH>0) THEN
			bprob = 1.D0-EXP(-DH*beta)
			CALL RANDOM_NUMBER(r2)
			IF (r2 < bprob) THEN
				down=1
				!WRITE(*,*) "*ABAJO", r2,bprob, vertex, vertex2
			ENDIF
		ENDIF 

		! Conexión hacia la derecha
		DH = 2*DOT_PRODUCT(rW,S(i,j,k,:))*DOT_PRODUCT(rW,S(i,j,ip(k),:))
		IF (DH>0) THEN
			bprob = 1.D0-EXP(-DH*beta)
			CALL RANDOM_NUMBER(r3)
			IF (r3 < bprob) THEN
				right=1
				!WRITE(*,*) "*DERECHA", r3,bprob, vertex, vertex2
			ENDIF
		ENDIF

		! Conexión hacia afuera
		DH = 2*DOT_PRODUCT(rW,S(i,j,k,:))*DOT_PRODUCT(rW,S(ip(i),j,k,:))
		IF (DH>0) THEN
			bprob = 1.D0-EXP(-DH*beta)
			CALL RANDOM_NUMBER(r1)
			IF (r1 < bprob) THEN
				outside=1
				!WRITE(*,*) "*AFUERA", r1,bprob, vertex, vertex2
			ENDIF
		ENDIF

		links: SELECT CASE(outside*4+down*2+right)
		CASE(0); tags(i,j,k)=0	!ninguna
		CASE(1); tags(i,j,k)=1	!derecha
		CASE(2); tags(i,j,k)=2	!abajo
		CASE(3); tags(i,j,k)=3	!derecha y abajo
		CASE(4); tags(i,j,k)=4	!afuera
		CASE(5); tags(i,j,k)=5	!afuera y derecha
		CASE(6); tags(i,j,k)=6	!afuera y abajo
		CASE(7); tags(i,j,k)=7	!todas
		CASE DEFAULT; EXIT
		END SELECT links

	ENDDO; ENDDO; ENDDO

	CALL identify_clusters
END SUBROUTINE multicluster

! Acción del constraint
SUBROUTINE constraint_algorithm
	INTEGER:: i,j,k
	INTEGER :: down, right, outside
	REAL(8):: cos_theta

	tags=0

	cos_theta=0.D0

	CALL rand_vec

	DO i=1,L; DO j=1,L; DO k=1,L
		outside=0; down=0; right=0;		

		! Conexión hacia abajo
		S_prime = S(i,ip(j),k,:)-2*rW*DOT_PRODUCT(rW,S(i,ip(j),k,:))
		cos_theta=DOT_PRODUCT(S(i,j,k,:),S_prime)

		IF (cos_theta<cosd) THEN
			down=1
		ENDIF
		

		! Conexión hacia la derecha
		S_prime = S(i,j,ip(k),:)-2*rW*DOT_PRODUCT(rW,S(i,j,ip(k),:))
		cos_theta=DOT_PRODUCT(S(i,j,k,:),S_prime)
		
		IF (cos_theta<cosd) THEN
			right=1
		ENDIF					

		! Conexión hacia afuera
		S_prime = S(ip(i),j,k,:)-2*rW*DOT_PRODUCT(rW,S(ip(i),j,k,:)) 
		cos_theta=DOT_PRODUCT(S(i,j,k,:),S_prime)
		
		IF (cos_theta<cosd) THEN
			outside=1
		ENDIF

		links: SELECT CASE(outside*4+down*2+right)
		CASE(0); tags(i,j,k)=0	!ninguna
		CASE(1); tags(i,j,k)=1	!derecha
		CASE(2); tags(i,j,k)=2	!abajo
		CASE(3); tags(i,j,k)=3	!derecha y abajo
		CASE(4); tags(i,j,k)=4	!afuera
		CASE(5); tags(i,j,k)=5	!afuera y derecha
		CASE(6); tags(i,j,k)=6	!afuera y abajo
		CASE(7); tags(i,j,k)=7	!todas
		CASE DEFAULT; EXIT
		END SELECT links

	ENDDO; ENDDO; ENDDO

	CALL identify_clusters
END SUBROUTINE constraint_algorithm

INTEGER FUNCTION find(p)
	INTEGER :: p
	INTEGER :: root, next

	!Recursive command to find the root/parent of the component
	root=p
	DO WHILE(root/=parents(root))
		root = parents(root)
	ENDDO

	!Recursive command to assign the minimun value to parent[p]
	DO WHILE(p/=root)
		next = parents(p)
		parents(p) = root
		p = next
	ENDDO

	find = root
	RETURN
END FUNCTION find

SUBROUTINE unify(p,q)
	INTEGER, INTENT(IN) :: p,q
	INTEGER :: root1, root2

	root1 = find(p)
	root2 = find(q)

	IF (root1==root2) RETURN
	
	! Merge two components together
	IF (root1<root2) THEN
		parents(root2) = root1
	ELSE
		parents(root1) = root2
	ENDIF
END SUBROUTINE unify

SUBROUTINE identify_clusters
	INTEGER :: i,j,k, ii,jj,kk, vertex
	INTEGER :: outside, right, down
	INTEGER :: small, medium, large
	REAL(8) :: r05 	! random number

	parents = [ (i, i=1,L**3) ]

	!WRITE(*,'(64I8)') parents

	DO i=1,L; DO j=1,L; DO k=1,L
		vertex = (i-1)*L2 + (j-1)*L + k
		outside = (ip(i)-1)*L2 + (j-1)*L + k
		right = (i-1)*L2 + (j-1)*L + ip(k)
		down = (i-1)*L2 + (ip(j)-1)*L + k

		search: SELECT CASE(tags(i,j,k))
		CASE(1)
			!WRITE(*,*) '----',vertex,right,'CASE 1' 
			CALL unify(vertex,right)
		CASE(2)
			!WRITE(*,*) '----',vertex,down,'CASE 2' 
			CALL unify(vertex,down)
		CASE(3)
			!WRITE(*,*) '----',vertex,right,down,'CASE 3'
			CALL unify(vertex,right)
			CALL unify(vertex,down)
		CASE(4)
			!WRITE(*,*) '----',vertex,outside,'CASE 4' 
			CALL unify(vertex,outside)
		CASE(5)
			!WRITE(*,*) '----',vertex,outside,right,'CASE 5'
			CALL unify(vertex,outside)
			CALL unify(vertex,right)
		CASE(6)
			!WRITE(*,*) '----',vertex,outside,down,'CASE 6'
			CALL unify(vertex,outside)
			CALL unify(vertex,down)
		CASE(7)
			!WRITE(*,*) '----',vertex,outside,down,right,'CASE 7'
			CALL unify(vertex,outside)
			CALL unify(vertex,down)
			CALL unify(vertex,right)
		END SELECT search
		
	ENDDO; ENDDO; ENDDO

	DO i=1,L3
		NULLIFY(row(i)%head,row(i)%tail)
	ENDDO

	DO i=1,L3
		!WRITE(*,*) i, parents(i)
		CALL add_to_linked(i,row(find(parents(i)))%head,row(find(parents(i)))%tail)
	ENDDO

	num_of_clus = 0

	!WRITE(*,*)'-- LISTA DE CLUSTERS UF --'
	DO i=1,L3
		row(i)%out0=>row(i)%head
		IF (ASSOCIATED(row(i)%out0)) num_of_clus = num_of_clus + 1
		!WRITE(*,*) 'node', i, ASSOCIATED(row(i)%out0)

		CALL RANDOM_NUMBER(r05)
		IF (r05<0.5) THEN
			print0: DO
				IF(.NOT. ASSOCIATED(row(i)%out0)) EXIT	! Pointer valid?
				ii=(row(i)%out0%value-1)/L**2+1
				jj=(MOD(row(i)%out0%value-1,L**2))/L+1
				kk=MOD(row(i)%out0%value-1,L)+1
				!WRITE(*,*) row(i)%out0%value, ii,jj,kk ! elementos de la lista

				S(ii,jj,kk,:) = S(ii,jj,kk,:)-2.D0*rW*DOT_PRODUCT(rW,S(ii,jj,kk,:))
				row(i)%out0=>row(i)%out0%next
			END DO print0
		ENDIF
		!WRITE(*,*)
	ENDDO

	DO i=1,L3
		CALL destroyLList(row(i)%head)
	ENDDO
END SUBROUTINE identify_clusters

! Algoritmo metropolis
SUBROUTINE metropolis
	INTEGER :: i,j,k, ll
	REAL(8):: r1, x

	DH=0.0D0
	acceptance_rate = 0.D0

	DO i=1,L; DO j=1,L; DO k=1,L

		CALL RANDOM_NUMBER(x)
		x = dPI*(x-0.5D0)
		S_prime = [DCOS(x), DSIN(x)]

		DH = DOT_PRODUCT(S(i,j,k,:),S(ip(i),j,k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,ip(j),k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,j,ip(k),:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(im(i),j,k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,im(j),k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,j,im(k),:))&
			- DOT_PRODUCT(S_prime,S(ip(i),j,k,:))&
			- DOT_PRODUCT(S_prime,S(i,ip(j),k,:))&
			- DOT_PRODUCT(S_prime,S(i,j,ip(k),:))&
			- DOT_PRODUCT(S_prime,S(im(i),j,k,:))&
			- DOT_PRODUCT(S_prime,S(i,im(j),k,:))&
			- DOT_PRODUCT(S_prime,S(i,j,im(k),:))
			

		IF (DH < 0.0) THEN
			S(i,j,k,:) = S_prime
			!acceptance_rate = acceptance_rate + 1.D0
		ELSE
			bprob = DEXP(-DH*beta)
			CALL RANDOM_NUMBER(r1)
			IF (r1<bprob) THEN
				S(i,j,k,:) = S_prime
				!acceptance_rate = acceptance_rate + 1.D0
			ENDIF
		ENDIF

	ENDDO; ENDDO; ENDDO
END SUBROUTINE metropolis

SUBROUTINE Glauber
	INTEGER :: i,j,k, ll
	REAL(8):: r1, x

	DH=0.0D0

	DO i=1,L; DO j=1,L; DO k=1,L

		CALL RANDOM_NUMBER(x)
		x = dPI*(x-0.5D0)
		S_prime = [DCOS(x), DSIN(x)]

		DH = DOT_PRODUCT(S(i,j,k,:),S(ip(i),j,k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,ip(j),k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,j,ip(k),:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(im(i),j,k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,im(j),k,:))&
			+ DOT_PRODUCT(S(i,j,k,:),S(i,j,im(k),:))&
			- DOT_PRODUCT(S_prime,S(ip(i),j,k,:))&
			- DOT_PRODUCT(S_prime,S(i,ip(j),k,:))&
			- DOT_PRODUCT(S_prime,S(i,j,ip(k),:))&
			- DOT_PRODUCT(S_prime,S(im(i),j,k,:))&
			- DOT_PRODUCT(S_prime,S(i,im(j),k,:))&
			- DOT_PRODUCT(S_prime,S(i,j,im(k),:))
			

		bprob = DEXP(-DH*beta)/(1+DEXP(-DH*beta))
		CALL RANDOM_NUMBER(r1)
		IF (r1<bprob) THEN
			S(i,j,k,:) = S_prime
		ENDIF

	ENDDO; ENDDO; ENDDO
END SUBROUTINE Glauber

! Algoritmo metropolis random
SUBROUTINE metropolis_random
	INTEGER :: i,j,k, ii,jj,kk
	REAL(8):: r1, x, x0

	DH=0.0D0
	!acceptance_rate = 0.D0

	DO i=1,L3

		CALL RANDOM_NUMBER(x0); ii=x0*L+1
		CALL RANDOM_NUMBER(x0); jj=x0*L+1
		CALL RANDOM_NUMBER(x0); kk=x0*L+1

		CALL RANDOM_NUMBER(x)
		x = dPI*(x-0.5D0)
		S_prime = [DCOS(x), DSIN(x)]

		DH = DOT_PRODUCT(S(ii,jj,kk,:),S(ip(ii),jj,kk,:))&
			+ DOT_PRODUCT(S(ii,jj,kk,:),S(ii,ip(jj),kk,:))&
			+ DOT_PRODUCT(S(ii,jj,kk,:),S(ii,jj,ip(kk),:))&
			+ DOT_PRODUCT(S(ii,jj,kk,:),S(im(ii),jj,kk,:))&
			+ DOT_PRODUCT(S(ii,jj,kk,:),S(ii,im(jj),kk,:))&
			+ DOT_PRODUCT(S(ii,jj,kk,:),S(ii,jj,im(kk),:))&
			- DOT_PRODUCT(S_prime,S(ip(ii),jj,kk,:))&
			- DOT_PRODUCT(S_prime,S(ii,ip(jj),kk,:))&
			- DOT_PRODUCT(S_prime,S(ii,jj,ip(kk),:))&
			- DOT_PRODUCT(S_prime,S(im(ii),jj,kk,:))&
			- DOT_PRODUCT(S_prime,S(ii,im(jj),kk,:))&
			- DOT_PRODUCT(S_prime,S(ii,jj,im(kk),:))

		IF (DH < 0.0) THEN
			S(ii,jj,kk,:) = S_prime
			!acceptance_rate = acceptance_rate + 1.D0
		ELSE
			bprob = DEXP(-DH*beta)
			CALL RANDOM_NUMBER(r1)
			IF (r1<bprob) THEN
				S(ii,jj,kk,:) = S_prime
				!acceptance_rate = acceptance_rate + 1.D0
			ENDIF
		ENDIF

	ENDDO
END SUBROUTINE metropolis_random

SUBROUTINE heat_bath_random
    INTEGER :: i,j,k, ii,jj,kk
	REAL(8) :: x0
    REAL(8) :: sigma(2), th_sigma, norm_sigma, cos_th
    REAL(8) :: bH, theta = 0.D0, Pxtheta = 0.D0
    REAL(8) :: rx

    DO i=1,L3
        CALL RANDOM_NUMBER(x0); ii=x0*L+1
        CALL RANDOM_NUMBER(x0); jj=x0*L+1
        CALL RANDOM_NUMBER(x0); kk=x0*L+1

        sigma(:) = S(ip(ii),jj,kk,:) + S(ii,ip(jj),kk,:) +&
            S(im(ii),jj,kk,:) + S(ii,im(jj),kk,:) +&
            S(ii,jj,ip(kk),:) + S(ii,jj,im(kk),:)

        IF (sigma(1) > 0 .AND. sigma(2) > 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1))
        ELSEIF (sigma(1) < 0 .AND. sigma(2) > 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1)) + PI
        ELSEIF (sigma(1) < 0 .AND. sigma(2) < 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1)) - PI
        ELSEIF (sigma(1) > 0 .AND. sigma(2) < 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1)) 
        ENDIF

        norm_sigma = NORM2(sigma)
        bH = beta * norm_sigma
        CALL urand(theta, Pxtheta, bH)
        theta = theta + th_sigma

        S(ii,jj,kk,:) = [DCOS(theta), DSIN(theta)]
    ENDDO
END SUBROUTINE heat_bath_random

SUBROUTINE heat_bath
	INTEGER :: i,j,k
    REAL(8) :: sigma(2), th_sigma, norm_sigma, cos_th
    REAL(8) :: bH, theta = 0.D0, Pxtheta = 0.D0
    REAL(8) :: rx
    
    DO i=1,L; DO j=1,L; DO k=1,L
        sigma(:) = S(ip(i),j,k,:) + S(i,ip(j),k,:) + S(im(i),j,k,:) + S(i,im(j),k,:) + S(i,j,ip(k),:) + S(i,j,im(k),:)

        IF (sigma(1) > 0 .AND. sigma(2) > 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1))
        ELSEIF (sigma(1) < 0 .AND. sigma(2) > 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1)) + PI
        ELSEIF (sigma(1) < 0 .AND. sigma(2) < 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1)) - PI
        ELSEIF (sigma(1) > 0 .AND. sigma(2) < 0) THEN
            th_sigma = ATAN(sigma(2) / sigma(1)) 
        ENDIF

        norm_sigma = NORM2(sigma)
        bH = beta * norm_sigma
        CALL urand(theta, Pxtheta, bH)
        theta = theta + th_sigma

        S(i,j,k,:) = [DCOS(theta), DSIN(theta)]
    ENDDO; ENDDO; ENDDO
END SUBROUTINE heat_bath

SUBROUTINE urand(thet, Pxthet, a)
    REAL(8), INTENT(IN) :: a
    REAL(8), INTENT(INOUT) :: thet, Pxthet
    REAL(8), PARAMETER :: epsilon = 0.001
    REAL(8), PARAMETER :: a_star = 0.798953686083986
    REAL(8) :: omega, omega_p
    REAL(8) :: delta0, alpha, beta0, beta1, beta2, Harg,Hab,factorb, Gab, gax
    REAL(8) :: Ptheta
    INTEGER :: i
    REAL(8) :: inv_beta0, inv_alpha, inv_beta0_1, inv_a, factorb1

    inv_a = 1.D0 / a
    
    delta0 = 0.35D0 * MAX(0.D0, a - a_star) + &
        & 1.03D0 * DSQRT(MAX(0.D0, a - a_star))
    !WRITE(*,*) delta0

    alpha = MIN(DSQRT(a * (2.D0 - epsilon)), &
            MAX(DSQRT(epsilon * a), delta0))
    !WRITE(*,*) alpha
    inv_alpha = 1.D0 / alpha

    beta1 = alpha**2 / a
    beta2 = (DCOSH(PI * alpha) - 1.D0) / (DEXP(2* a) - 1)
    beta0 = MAX(beta1, beta2) - 1.D0
    !WRITE(*,*) beta0
    inv_beta0 = 1.D0 / beta0
    inv_beta0_1 = 1.D0 / (1 + beta0)

    CALL RANDOM_SEED
    omega_p = 1.D0; gax = 0.D0
    DO WHILE(omega_p > gax)
        CALL RANDOM_NUMBER(omega)
        CALL RANDOM_NUMBER(omega_p)
        !WRITE(*,*) omega, omega_p

        factorb = DSQRT((1 + beta0) / (1 - beta0))
        factorb1 = DSQRT((1 - beta0) / (1 + beta0))
        Harg = factorb * DTAN((2 * omega - 1) * ATAN(DTANH(PI * alpha * 0.5D0) * factorb1)) 
        Hab = DLOG((1 + Harg) / (1 - Harg)) * inv_alpha
        !WRITE(*,*) "h", Hab

        Gab = 1 - DCOS(Hab) - DLOG(1 + (DCOSH(alpha * Hab) - 1) * inv_beta0_1) * inv_a
        !WRITE(*,*) Gab

        gax = DEXP(-a * Gab)
        !WRITE(*,*) omega_p, gax

        Ptheta = DEXP(a * DCOS(Hab)) / (2 * PI * BESSEL_J0(a))
    END DO

    thet = Hab
    Pxthet = Ptheta
END SUBROUTINE urand

! Función para calcular el ángulo dirigido ente dos espines
FUNCTION delt_phi(s1, s2)
	REAL(8), INTENT(IN):: s1(2), s2(2)
	REAL(8) :: delt_phi

	cos_phi = DOT_PRODUCT(s1,s2)/(NORM2(s1)*NORM2(s2))
	!cos_phi = DOT_PRODUCT(s1,s2)/(DSQRT(DOT_PRODUCT(s1,s1)*DOT_PRODUCT(s2,s2)))

	IF (ABS(cos_phi) > 1) THEN
		cos_phi = 1
	ENDIF

	IF ( Cross2(s1,s2) < 0 ) THEN
		delt_phi = -ACOS(cos_phi)
	ELSE
		delt_phi = ACOS(cos_phi)
	ENDIF
END FUNCTION delt_phi

! Calcula el producto cruz de vectores de 2 dimensiones
FUNCTION Cross2(s1,s2)
	REAL(8) :: Cross2
	REAL(8), INTENT(IN)::s1(2), s2(2)

	Cross2 = s1(1)*s2(2)-s1(2)*s2(1)
END FUNCTION Cross2

! Calcula vortices
SUBROUTINE calc_vortices
	INTEGER :: i,j,k, v

	nvx=0; nvy=0; nvz=0
	num_of_vort=0; num_of_av=0

	!	Fondo(1)		Piso(2)			Lateral Izq(3)
	!	
	!	x4 <----- x3	x1 <----- x2	x7 <----- x4	
	!	|		   |	|		   |	|		   |
	!	|		   |	|		   |	|		   |
	!	x1 -----> x2	x5 -----> x6	x5 -----> x1
	!

	DO i=1,L; DO j=1,L; DO k=1,L
		Dphi12 = delt_phi(S(i,j,k,:), S(i,ip(j),k,:))
		Dphi23 = delt_phi(S(i,ip(j),k,:), S(i,ip(j),ip(k),:))
		Dphi34 = delt_phi(S(i,ip(j),ip(k),:), S(i,j,ip(k),:))
		Dphi41 = delt_phi(S(i,j,ip(k),:), S(i,j,k,:))

		vortex2=Dphi12+Dphi23+Dphi34+Dphi41 ! Fondo(1)
		IF (NINT(vortex2) > 0) THEN
			num_of_vort = num_of_vort+1
		ELSE IF (NINT(vortex2)<0) THEN
			num_of_av = num_of_av+1
		ENDIF

		Dphi15 = delt_phi(S(i,j,k,:), S(ip(i),j,k,:))
		Dphi56 = delt_phi(S(ip(i),j,k,:), S(ip(i),ip(j),k,:))
		Dphi62 = delt_phi(S(ip(i),ip(j),k,:), S(i,ip(j),k,:))
		Dphi21 = delt_phi(S(i,ip(j),k,:), S(i,j,k,:))

		vortex1=Dphi21+Dphi15+Dphi56+Dphi62 ! Piso(2)
		IF (NINT(vortex1) > 0) THEN
			num_of_vort = num_of_vort+1
		ELSE IF (NINT(vortex1)<0) THEN
			num_of_av = num_of_av+1
		ENDIF			

		Dphi14 = delt_phi(S(i,j,k,:), S(i,j,ip(k),:))
		Dphi47 = delt_phi(S(i,j,ip(k),:), S(ip(i),j,ip(k),:))
		Dphi75 = delt_phi(S(ip(i),j,ip(k),:), S(ip(i),j,k,:))
		Dphi51 = delt_phi(S(ip(i),j,k,:), S(i,j,k,:))

		vortex3=Dphi14+Dphi47+Dphi75+Dphi51 ! Lateral(3)

		IF (NINT(vortex3) > 0) THEN
			num_of_vort = num_of_vort+1
		ELSE IF (NINT(vortex3)<0) THEN
			num_of_av = num_of_av+1
		ENDIF

		sum_vortex=sum_vortex+vortex1+vortex2+vortex3

		vortex1 = vortex1/(2*PI) !-Piso(2)
		vortex2 = vortex2/(2*PI) !-Fondo(1)
		vortex3 = vortex3/(2*PI) !-Lateral Izq(3)

		Sv3(i,j,k,1) = NINT(vortex1)	! Piso(2)
		Sv3(i,j,k,2) = NINT(vortex2)	! Fondo(1)
		Sv3(i,j,k,3) = NINT(vortex3)	! Lateral(3)

		nvx=nvx+ABS(vortex1) !-Piso(2)
		nvy=nvy+ABS(vortex2) !-Fondo(1)
		nvz=nvz+ABS(vortex3) !-Lateral Izq(3)

	!WRITE(unit5,'(3I3,1X,3I3)') i-1, j-1, k-1, NINT(vortex1),NINT(vortex2),NINT(vortex3)!, S(i,j,k,:), S(ip(i),j,k,:), &
	!S(ip(i),ip(j),k,:), S(i,ip(j),k,:)
		
	!v=(i-1)*L2 + (j-1)*L + k
	!WRITE(*,'(I6,1X,3I3)') 3*v-2, NINT(vortex1),NINT(vortex2),NINT(vortex3)!, S(i,j,k,:), S(ip(i),j,k,:)

	ENDDO; ENDDO; ENDDO
END SUBROUTINE calc_vortices

SUBROUTINE calc_F
	INTEGER :: i,j,k, i2, j2, k2

	DO i=1,L; DO j=1,L; DO k=1,L
		do i2 = 1,L; do j2 = 1,L; do k2 = 1,L
			F = F + dot_product(S(i2,j2,k2,:),S(i,j,k,:))*cos_arr(i, i2)
		enddo; enddo; enddo
	ENDDO; ENDDO; ENDDO
END SUBROUTINE calc_F

! Creates a linked list
SUBROUTINE add_to_linked(temp, head_, tail_)

	! Data dictionary: declare variable  types & definitions
	TYPE (my_list), POINTER :: head_		! Ptr to head of list
	TYPE (my_list), POINTER :: tail_		! Ptr to tail list
	!TYPE (my_value), POINTER :: ptr_out	! Temporary pointer
	INTEGER, INTENT(IN) :: temp				! Temporary variable
	INTEGER :: istat

	IF(.NOT. ASSOCIATED(head_)) THEN	! No values in list
		ALLOCATE(head_, STAT=istat)		! Allocate new value
		tail_ => head_					! Tail points to new value
		NULLIFY (tail_%next)			! Nullify 'next' in new value
		tail_%value = temp				! Store new number
	ELSE								! Values already in list
		ALLOCATE(tail_%next,STAT=istat)	! Allocate new value
		tail_=>tail_%next				! Tail points to new value
		NULLIFY(tail_%next)				! Nullify 'next' in new value
		tail_%value=temp				! Store number
	END IF

END SUBROUTINE add_to_linked

! Free memory space
SUBROUTINE destroyLList(list)
	TYPE (my_list), POINTER :: list
	TYPE (my_list), POINTER :: next
	TYPE (my_list), POINTER :: dCurrent=>NULL(), dNext=>NULL()

	IF(.NOT. ASSOCIATED(list)) RETURN

	dCurrent => list

	!-Dealocate all data nodes in list
	DO WHILE(ASSOCIATED(dCurrent))
		dNext => dCurrent%next
		NULLIFY(dCurrent%next)
		DEALLOCATE(dCurrent)
		dCurrent => dNext
	END DO

	!-Deallocate the list itself
	!DEALLOCATE(list)
	NULLIFY(list)
END SUBROUTINE destroyLList

SUBROUTINE init_vars
	INTEGER :: i

	E=0.d0;	Mag=0.d0; Corr=0.0d0; O=0.0d0; F=0.D0

	Dphi12=0d0;	Dphi23=0d0;	Dphi34=0d0;	Dphi41=0d0
	Dphi21=0d0;	Dphi15=0d0;	Dphi56=0d0;	Dphi62=0d0
	Dphi14=0d0;	Dphi47=0d0;	Dphi75=0d0;	Dphi51=0d0

	vortex1=0d0; vortex2=0d0; vortex3=0d0
	sum_vortex=0d0

	d6=(/(i, i=1,6)/)
END SUBROUTINE init_vars

SUBROUTINE init_arrays
	INTEGER :: i, i2

	do i = 1,L; do i2 = 1,L
		cos_arr(i, i2) = DCOS(dPI*(i-i2)*invL)
	enddo; enddo
END SUBROUTINE init_arrays

! Subrrutina para arrays dinámicos con elementos enteros
SUBROUTINE AddToArray(list,element)
	INTEGER :: i, isize
	INTEGER, INTENT(IN) :: element
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: list
	INTEGER, DIMENSION(:), ALLOCATABLE :: elist

	IF (ALLOCATED(list)) THEN
		isize = SIZE(list)
		ALLOCATE(elist(isize+1))
		elist(1:isize) = list
		elist(isize+1) = element

		DEALLOCATE(list)
		CALL MOVE_ALLOC(elist, list)
	ELSE
		ALLOCATE(list(1))
		list(1) = element
	ENDIF
END SUBROUTINE AddToArray

! Subrrutina para arrays dinámicos con elementos reales(8)
SUBROUTINE AddDblToArray(list,element)
	INTEGER :: i, isize
	REAL(8), INTENT(IN) :: element
	REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: list
	REAL(8), DIMENSION(:), ALLOCATABLE :: elist

	IF (ALLOCATED(list)) THEN
		isize = SIZE(list)
		ALLOCATE(elist(isize+1))
		elist(1:isize) = list
		elist(isize+1) = element

		DEALLOCATE(list)
		CALL MOVE_ALLOC(elist, list)
	
	ELSE
		ALLOCATE(list(1))
		list(1) = element
	ENDIF
END SUBROUTINE AddDblToArray

! Imprime fecha y hora
SUBROUTINE FechaYHora
  CHARACTER(20) :: fecha_hora
  INTEGER :: date_arry(8)
  CALL DATE_AND_TIME(VALUES = date_arry)

  ! Obtener los valores de la fecha y hora
  c_year = date_arry(1)
  c_month = date_arry(2)
  c_day = date_arry(3)
  c_hour = date_arry(5)
  c_min = date_arry(6)
  c_sec = date_arry(7)

  ! Formatear la fecha y hora
  WRITE(fecha_hora, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
       c_year, c_month, c_day, c_hour, c_min, c_sec

  !WRITE(ounit,*) "Fecha y hora actual:", fecha_hora
END SUBROUTINE FechaYHora

!Condiciones de frontera periodicas
INTEGER FUNCTION ip(pos) RESULT(ppos)
        INTEGER, INTENT(IN) :: pos
        ppos = MODULO(pos,L)+1
    END FUNCTION ip

    INTEGER FUNCTION im(pos) RESULT(mpos)
        INTEGER, INTENT(IN) :: pos
        mpos = MODULO(pos-2,L)+1
    END FUNCTION im

! Convierte enteros a caracter
CHARACTER(len=20) FUNCTION str(k)
	INTEGER, INTENT(IN) :: k

	WRITE(str,*) k
	str = ADJUSTL(str)
END FUNCTION str

! Convierte reales a caracter
CHARACTER(len=12) FUNCTION real2str(r)
	REAL(4), INTENT(IN) :: r
	CHARACTER(12) :: str_

	WRITE(str_,"(F12.4)") r
	real2str = ADJUSTL(str_)
END FUNCTION real2str

! Convierte double a caracter
CHARACTER(len=12) FUNCTION dble2str(d)
	REAL(8), INTENT(IN) :: d
	CHARACTER(12) :: str_

	WRITE(str_,"(F12.6 )") d
	dble2str = ADJUSTL(str_)
END FUNCTION dble2str

!Programa que llama a la semilla
SUBROUTINE init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 1000000*(/ (i - 1, i = 1, n) /)
  !write (*, *) seed
  call random_seed(put = seed)

  deallocate(seed)
END SUBROUTINE init_random_seed
END MODULE MC
