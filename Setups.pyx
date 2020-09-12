from Pygrad.Pygrad cimport Pygrad
from libc.math cimport acos, sqrt
cimport numpy as np
import  numpy as np
import cython

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef Setup(Pygrad object):
    """
    This function sets up the given Boltz object. It fills the values of the main constants. 
    
    The object parameter is the Boltz object to be setup.
    """
    cdef double TwoPi,  ElectronCharge, ElectronMass, AMU, BoltzmannConst_eV, BoltzmannConst_eVJ, MassOverChargeDivTen, ALOSCH,  ZeroCelcius, OneAtmosphere, TotFrac, HBAR
    cdef long long IH,  i
    cdef double FracMol = 0.0

    TwoPi = 2.0 * np.pi
    object.RhydbergConst = <float>(13.60569253)
    object.PIR2 = 8.7973554297e-17
    ElectronCharge = 1.602176565e-19
    ElectronMass = 9.10938291e-31
    AMU = 1.660538921e-27
    BoltzmannConst_eV = 8.6173324e-5
    BoltzmannConst_eVJ = 1.3806488e-23
    MassOverChargeDivTen = 1.758820088e10
    ALOSCH = 2.6867805e19
    ZeroCelcius = 273.15
    OneAtmosphere = 760.0
    HBAR = 6.58211928e-16                                     
    EMS = 510998.928
    VC = 299792458.0
    RE = 2.8179403267e-13    
    ALPH = 137.035999074
    EOVM = np.sqrt(2.0*ElectronCharge/ElectronMass)*100.0 
    object.CONST1 = MassOverChargeDivTen / 2.0 * 1.0e-19
    object.CONST2 = object.CONST1 * 1.0e-02
    object.CONST3 = sqrt(0.2 * MassOverChargeDivTen) * 1.0e-9
    object.CONST4 = object.CONST3 * ALOSCH * 1e-15
    object.CONST5 = object.CONST3 / 2.0
    object.NANISO = 2
    object.PresTempCor = ZeroCelcius * object.Pressure_Torr / (OneAtmosphere * (ZeroCelcius + object.TemperatureCentigrade) * 100.0)

    TotFrac = 0.0

    if object.NumberOfGases == 0 or object.NumberOfGases > 6:
        raise ValueError("Error in Gas Input")
    for J in range(object.NumberOfGases):
        if object.GasIDs[J] == 0 or object.GasFractions[J] <= 0:
            raise ValueError("Error in Gas Input")
        TotFrac += object.GasFractions[J]
    if abs(TotFrac - 100) >= 1e-6:
        raise ValueError("Error in Gas Input")

    object.icount = 0
    if object.imip == 1:
        object.icount = 1
    object.EnergySteps = 20000
    check_parameters(object)
    object.Max_Electron_Energy = object.InitialElectronEnergy * 50.0
    object.FinalElectronEnergy = object.InitialElectronEnergy*1.0001 + 760.0*object.Max_Electron_Energy/object.Pressure_Torr*(object.TemperatureCentigrade+ZeroCelsius)/293.15 + object.EField
    if object.FinalElectronEnergy < 1.01 * object.InitialElectronEnergy:
        object.FinalElectronEnery = 1.01 * object.InitialElectronEnergy
    if object.BeamDirection == 1:
        object.Phi = 0.0                          
        object.Theta = 0.0
    elif object.BeamDirection == -1:
        object.Phi = 0
        object.Theta = np.pi
    elif object.BeamDirection == 0:
        object.Phi = 0.0
        object.Theta = np.pi/2.0
    elif object.BeamDirection == 2:
        R3 = drand48(0)
        Phi = TWOPI * R3
        R4 = drand48(0)
        Theta = np.acos(1-2*R4)

#Check if a variable is in a valid range
cpdef check_var( var, name, low=0, high=1):
    if var < low or var > high:
        raise ValueError('Invalid value for '+name)

cpdef check_parameters(object):
    check_var(object.imip,'IMIP',1,5)
    check_var(object.BeamDirection,'NDVEC',-1,2)
    check_var(object.OutputVerbosity,'IWRITE',0,2)
    check_var(object.DetectorEfficiency,'DETEFF',0,100)
    if (object.imip == 1 and object.nDelta > 1e5) or (object.imip != 1 and object.nDelta > 1e4):
        raise ValueError('Invalid value for NDELTA')
    if object.imip == 3 and object.InitialElectronEnergy > 3e6:
        raise ValueError('X-ray energy is too high')
    if object.estart <= 0 or object.etherm <= 0 or object.ecut < 0:
        raise ValueError('Negative energy')
    if object.TemperatureCentigrade <= -ZeroCelsius:
        raise ValueError('Absolute zero')
    if object.Pressure_Torr <= 0:
        raise ValueError('Negative pressure')
