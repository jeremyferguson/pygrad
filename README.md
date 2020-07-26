# PyGrad

> Python port of Degrad: http://degrad.web.cern.ch/degrad/

## Prerequisites

- Python 3.x
- argparse, numpy, os python libraries
- unittest library is necessary for running tests


## Features

- Currently, the only functionality is to do the initial setup and density effect calculations

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
`cd tests/fortran_tests`
`make`
In order to make all the tests for the old Fortran functions.  Then, 
`cd ../..` to return to the top level and 
`python3 -m unittest discover -s tests` to run all of the tests

## Contributors
- Jeremy Ferguson: [github](https://github.com/jeremyferguson) 
email:jmfergie@gmail.com

## TODO list
- refactor plot arrays

## Possible changes from Degrad
- random number generator
