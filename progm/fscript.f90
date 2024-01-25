program EditarParametros
  implicit none
  character(80) :: linea, chars
  integer :: unidad_entrada=11, unidad_salida=12, ios
  real(8) :: nuevo_valor, numero_real

  ! Abre el archivo de entrada para lectura
  open(unit=unidad_entrada, file='parameters.f90', status='old')

  ! Abre un nuevo archivo de salida
  open(unit=unidad_salida, file='param_modificado.f90', status='unknown')

  ! Lee el archivo de entrada línea por línea
  do
    read(unidad_entrada, '(A)', iostat=ios) linea
    if (ios /= 0) exit ! Finaliza al llegar al final del archivo

    ! Busca la línea que contiene betaL y modifica el valor
    if (index(linea, 'betas(') /= 0) then
        numero_real = real(0.9d0, 8)
        WRITE(chars, '(F12.6)') numero_real
      ! Modifica el valor de betaL según sea necesario
      linea = '  REAL(8), PARAMETER :: betaL = ' // &
              trim(adjustl(chars)) // 'D0'
    end if

    ! Escribe la línea en el archivo de salida
    write(unidad_salida, '(A)') linea
  end do

  ! Cierra los archivos
  close(unidad_entrada)
  close(unidad_salida)

  print *, 'Archivo modificado: param_modificado.f90'

end program EditarParametros