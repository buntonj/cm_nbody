nbody : nbody.f90
	gfortran -c constants.f90
	gfortran -c lennardjones.f90
	gfortran -c morse.f90
	gfortran -O3 nbody.f90 constants.o lennardjones.o morse.o
