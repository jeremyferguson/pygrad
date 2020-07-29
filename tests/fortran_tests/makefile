all: density mixerc
density: density.f runfortran.py density-tests/
	gfortran -o density density.f
	python3 runfortran.py density density-tests/ 1

mixerc: mixerc.f runfortran.py mixerc-tests/
	gfortran -o mixerc mixerc.f
	python3 runfortran.py mixerc mixerc-tests/ 7
