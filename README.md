# PyGrad

> Python port of Degrad: http://degrad.web.cern.ch/degrad/

## Prerequisites

- Python 3.x
- Python libraries:
	- numpy
	- h5py

## Setup
`git clone https://github.com/jeremyferguson/pygrad
./config.sh`

The shell script `setup.sh` sets all the necessary environment variables, so it is necessary to `source` it upon opening a new terminal window.

## Features

- The current functionality includes loading some of the gas data and doing the initial density effect calculations

## Usage

- `python3 pygrad.py [-h] infile [outfile]`
###### Positional arguments:
- infile: Filename containing initial parameters
- outfile: Filename for the output of the file

 
## Documentation

- Currently, no documentation exists

## Tests

- Tests are created using the unittest module and can be run on their own by python
- To run all tests:
`python3 -m unittest discover -s tests`
You can also go into the `tests` directory and run any of the test files individually.

## Contributors
- Jeremy Ferguson: [github](https://github.com/jeremyferguson)  
email:jmfergie@gmail.com

## Possible changes from Degrad
- random number generator
