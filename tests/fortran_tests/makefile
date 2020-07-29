all: density mixerc
density: density.f runfortran.py density-tests/
	gfortran -o density density.f
	python3 runfortran.py density density-tests/ DENSITY.OUT

mixerc: mixerc.f runfortran.py mixerc-tests/
	gfortran -o mixerc mixerc.f
	python3 runfortran.py mixerc mixerc-tests/ MIXERC1.OUT MIXERC2.OUT MIXERC3.OUT MIXERC4.OUT MIXERC5.OUT MIXERC6.OUT MIXERC7.OUT
