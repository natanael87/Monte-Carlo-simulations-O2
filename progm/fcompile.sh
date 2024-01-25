rm *.o *.mod

gfortran -c parameters.f90
gfortran -c variables.f90
gfortran -c subrout_funcs.f90
gfortran -c mainf.f90

gfortran *.o -o L64_b0.454165.exe -llapack

rm *.o *.mod