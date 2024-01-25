PROGRAM ModificarArchivo
  IMPLICIT NONE

  CHARACTER(1000) :: linea
  INTEGER :: iunit_in, iunit_out, i, pos, ios
  LOGICAL :: found_L, found_beta

  ! Abrir el archivo de entrada
  OPEN(unit=iunit_in, file="parameters.f90", status="OLD")

  ! Abrir el archivo de salida
  OPEN(unit=iunit_out, file="archivo_modificado.txt", status="UNKNOWN")

  ! Inicializar banderas para controlar las modificaciones
  found_L = .FALSE.
  found_beta = .FALSE.

  ! Leer y modificar el archivo línea por línea
  DO
     READ(iunit_in, '(A)', iostat=ios) linea
     IF (ios /= 0) EXIT  ! Salir al final del archivo

     ! Buscar y modificar la línea que contiene "L = 8"
     IF (.NOT. found_L) THEN
        pos = INDEX(linea, "L = 8")
        IF (pos /= 0) THEN
           linea = linea(1:pos+1) // "16" // linea(pos+3:)
           found_L = .TRUE.
        END IF
     END IF

     ! Buscar y modificar la línea que contiene "beta(1) = [0.5D]"
     IF (.NOT. found_beta) THEN
        pos = INDEX(linea, "betas(")
        IF (pos /= 0) THEN
           linea = linea(1:pos+7) // "0.4D0, 0.5D0]" // linea(pos+18:)
           found_beta = .TRUE.
        END IF
     END IF

     ! Escribir la línea modificada en el archivo de salida
     WRITE(iunit_out, '(A)') linea
  END DO

  ! Cerrar los archivos
  CLOSE(iunit_in)
  CLOSE(iunit_out)

  IF (.NOT. found_L) PRINT *, "No se encontró la línea con 'L = 8'."
  IF (.NOT. found_beta) PRINT *, "No se encontró la línea con 'beta(1) = [0.5D0]'."

  PRINT *, "Archivo modificado y guardado como 'archivo_modificado.txt'."

END PROGRAM ModificarArchivo
