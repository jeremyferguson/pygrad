density: density.o
	gfortran -o density density.o
	python3 runfortran.py density density-tests/ DENSITY.OUT

density.o: density.f
	gfortran -c density.f
