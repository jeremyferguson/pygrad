
export PYGRAD_HOME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PYTHONPATH="$PYTHONPATH:$PYGRAD_HOME"

echo "Making Fortran test files"
cd $PYGRAD_HOME/tests/fortran_tests
make
cd ../..

echo "Creating gas database"
python3 writeh5.py
