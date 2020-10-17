from Pygrad cimport Pygrad
from libc.math cimport sin, cos, acos, asin, log, sqrt
from PyGasMix.Gasmix cimport Gasmix
from Ang cimport Ang
import sys
import cython

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef Mixer(Pygrad object):
    """Loads the initial gas values for the mixture, particularly the momentum cross sections and other related values."""
    cdef double AttachmentCrossSection[6][4000], ElectronCharge, JHI, JLOW, EnergyHigh, F2, BP, EnergyLow
    cdef int iEnergy, GasIndex, iProcess, p, Sum, J, i, j, iIonization, JJ, IL, I

    ElectronCharge = 1.60217656e-19


