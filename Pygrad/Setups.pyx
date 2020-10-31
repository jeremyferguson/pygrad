from Pygrad.Pygrad cimport Pygrad
from libc.math cimport acos, sqrt
cimport numpy as np
import  numpy as np
import cython

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)

cdef extern from "C/RM48.h":
    double DRAND48(double dummy)
    void RM48(double lenv)
    void RM48IN(int IJL, int NTOT, intNTOT2)


cdef double drand48(double dummy)
cdef double drand48(double dummy):
    return DRAND48(dummy)
cdef void setSeed(int seed):
    RM48IN(seed, 0, 0)
    return

cpdef Setup(Pygrad object):
    """
    This function sets up the given Pygrad object. It fills the values of the main constants. 
    
    The object parameter is the Pygrad object to be setup.
    """
    cdef double TotFrac
    object.last = 0
    object.tmax = 100.0
    object.nout = 10

    object.TwoPi = 2.0 * np.pi
    object.RhydbergConst = <float>(13.60569253)
    object.PIR2 = 8.7973554297e-17
    object.ElectronCharge = 1.602176565e-19
    object.ElectronMass = 9.10938291e-31
    object.AMU = 1.660538921e-27
    BoltzmannConst_eV = 8.6173324e-5
    BoltzmannConst_eVJ = 1.3806488e-23
    MassOverChargeDivTen = 1.758820088e10
    object.ALOSCH = 2.6867805e19
    object.ZeroCelcius = 273.15
    object.OneAtmosphere = 760.0
    object.HBAR = 6.58211928e-16                                     
    object.EMS = 510998.928
    object.VC = 299792458.0
    object.RE = 2.8179403267e-13    
    object.ALPH = 137.035999074
    object.EOVM = np.sqrt(2.0*object.ElectronCharge/object.ElectronMass)*100.0 

    object.EMASS2 = 1021997.804
    object.CONST = 1.873884e-20
    object.CONST1 = MassOverChargeDivTen / 2.0 * 1.0e-19
    object.CONST2 = object.CONST1 * 1.0e-02
    object.CONST3 = sqrt(0.2 * MassOverChargeDivTen) * 1.0e-9
    object.CONST4 = object.CONST3 * object.ALOSCH * 1e-15
    object.CONST5 = object.CONST3 / 2.0
    object.WhichAngularModel = 2
    object.A0 = 0.52917720859e-8
    object.BBCONST=16.0*np.pi*object.A0*object.A0*object.RhydbergConstant*object.RhydbergConstant/object.EMASS2
    
    #Density Effect Constants
    object.EIAV = np.array([np.nan,115.0,188.0,41.8,41.8,137.0,352.0,482.0,41.7,45.4,47.1,48.3,85.0,0.0,71.6,95.0,82.0,0.0,84.9,0.0,0.0,19.2,19.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,128.0,53.7,0.0,0.0,67.6,62.9,61.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,48.3,0.0,61.1,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    object.nElectrons = np.array([np.nan,42,18,2,2,10,36,54,10,18,26,34,22,0,10,16,14,0,22,0,0,2,2,0,0,0,0,0,0,0,70,10,0,0,18,26,34,0,0,0,0,0,0,0,34,0,34,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    object.X00 = np.array([np.nan,1.70,1.7635,2.2017,2.2017,2.0735,1.7158,1.5630,1.6263,1.5090,1.4339,1.3788,1.6294,0.0,1.7952,1.7541,1.7378,0.0,1.6477,0.0,0.0,1.8639,1.8639,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.6,1.6822,0.0,0.0,0.2529,0.2218,0.2046,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.3788,0.0,0.2046,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    object.AKS = np.array([np.nan,3.00,2.9618,5.8347,5.8347,3.5771,3.4051,2.7414,3.6257,3.6095,3.5920,3.4884,3.3227,0.0,3.5901,3.2913,3.2125,0.0,3.3318,0.0,0.0,5.7273,5.7273,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,3.6464,0.0,0.0,3.5477,3.4834,3.5415,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.4884,0.0,3.5415,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    object.AAA = np.array([np.nan,.18551,.19714,.13443,.13443,.08064,.07446,.23314,.09253,0.09627,0.09916,.10852,.11768,0.0,.08101,.11778,.15349,0.0,.11992,0.0,0.0,.14092,.14092,0.0,0.0,0.0,0.0,0.0,0.0,0.0,.177484,.08315,0.0,0.0,.08970,0.09878,0.09644,0.0,0.0,0.0,0.0,0.0,0.0,0.0,.10852,0.0,0.09644,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    object.PresTempCor = object.ZeroCelcius * object.Pressure_Torr / (object.OneAtmosphere * (object.ZeroCelcius + object.TemperatureCentigrade) * 100.0)
    object.ICFLG=0
    object.IRFLG=0
    object.IPFLG=0
    object.IBFLG=0
    object.LPEFLG=0

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
    object.FinalElectronEnergy = object.InitialElectronEnergy*1.0001 + 760.0*object.Max_Electron_Energy/object.Pressure_Torr*(object.TemperatureCentigrade+object.ZeroCelsius)/293.15 + object.EField
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
        Phi = object.TwoPi * R3
        R4 = drand48(0)
        Theta = np.acos(1-2*R4)

    object.drzinit = np.cos(object.Theta)
    object.drxinit = np.sin(object.Theta) * np.cos(object.Phi)
    object.dryinit = np.sin(object.Theta) * np.sin(object.Phi)
    object.ThermalEnergy = (ZeroCelcius + object.TemperatureCentigrade) * BoltzmannConst_eV
    for i in range(6):
        object.MoleculesPerCm3PerGas[i] = object.GasFractions[i] * object.PresTempCor * object.ALOSCH
    object.AN = 100.0 * object.PresTempCor * object.ALOSCH
    for i in range(6):
        object.VMoleculesPerCm3PerGas[i] = object.GasFractions[i] * object.PresTempCor * object.CONST3 * object.ALOSCH
    object.VAN = 100.0 * object.PresTempCor * object.CONST3 * object.ALOSCH
    

    if object.FinalElectronEnergy <= 20000.0:
        object.ElectronEnergyStep = object.FinalElectronEnergy/object.EnergySteps
        object.EHalf = object.EnergySteps/2.0
        for i in range(20000):
            assignIndex(object,i,i)
    elif object.FinalElectronEnergy <= 140000.0:
        object.ElectronEnergyStep = 1.0
        object.EHalf = 0.5
        for i in range(16000):
            assignIndex(object,i,i)
        EnergySteps1 = object.ElectronEnergyStep
        object.ElectronEnergyStep = (object.FinalElectronEnergy-16000.0)/4000.0
        for i in range(16000,20000):
            assignIndex(object,i,i-15999)
        object.ElectronEnergyStep = EnergySteps1
    else:
        object.ElectronEnergyStep = 1.0
        object.EHalf = 0.5
        for i in range(12000):
            assignIndex(object,i,i)
        EnergySteps1 = object.ElectronEnergyStep
        object.ElectronEnergyStep = 20.0
        for i in range(12000,16000):
            assignIndex(object,i,i-11999)
        object.ElectronEnergyStep = (object.FinalElectronEnergy-92000.0)/4000.0
        for i in range(16000,20000):
            assignIndex(object,i,i-15999)
        object.ElectronEnergyStep = EnergySteps1
    
    object.AngularSpeedOfRotation = MassOverChargeDivTen * object.BField_Mag * 1e-12
    if object.BField_Mag != 0.0:
        object.ElectricOverMag = object.EField * 1e-9 / object.BField_Mag

cpdef assignIndex(object,i,aj):
    object.E[i] = object.EHalf + object.EStep * aj
    object.SqrtEnergy[i] = sqrt(object.E[i])
    object.Gamma[i] = (object.EMS + object.E[i])/object.EMS
    object.Beta[i] = sqrt(1.0-1.0/(object.Gamma[i] ** 2))

#Check if a variable is in a valid range
def check_var( var, name, low=0, high=1):
    if var < low or var > high:
        raise ValueError('Invalid value for '+name)

cpdef check_parameters(Pygrad object):
    check_var(object.imip,'IMIP',1,5)
    check_var(object.BeamDirection,'NDVEC',-1,2)
    check_var(object.OutputVerbosity,'IWRITE',0,2)
    check_var(object.DetectorEfficiency,'DETEFF',0,100)
    if (object.imip == 1 and object.nDelta > 1e5) or (object.imip != 1 and object.nDelta > 1e4):
        raise ValueError('Invalid value for NDELTA')
    if object.imip == 3 and object.InitialElectronEnergy > 3e6:
        raise ValueError('X-ray energy is too high')
    if object.InitialElectronEnergy <= 0 or object.ThermalEnergy <= 0 or object.EnergyCut < 0:
        raise ValueError('Negative energy')
    if object.TemperatureCentigrade <= -ZeroCelsius:
        raise ValueError('Absolute zero')
    if object.Pressure_Torr <= 0:
        raise ValueError('Negative pressure')


#Calculate density effect. DENSITY subroutine in Degrad
cpdef calcDensity(Pygrad object):
    sum1 = 0.0
    sumdnom = 0.0
    hsum = 0.0
    for i in range(6):
        gasn = object.GasIDs[i]
        fracn = object.GasFractions[i]
        sum1 += fracn * object.nElectrons[gasn - 1] * np.log(object.eiav[gasn - 1])
        sumdnom += fracn * object.nElectrons[gasn - 1]
        hsum +=  object.MoleculesPerCm3PerGas[i] * object.nElectrons[gasn-1] 
    eibar = np.e ** (sum1/sumdnom)
    hwp1 = sqrt(4.0 * np.pi * hsum * object.RE ** 3) * object.ALPH * object.EMS
    delden = np.log(eibar/hwp1)
    cbar = 1.0 + 2.0 * delden
    if object.ngas != 1:
        if cbar < 10.0:
            x0 = 1.6
            x1 = 4.0
        elif cbar >= 4.0 and cbar < 10.5:
            x0 = 1.7
            x1 = 4.0
        elif cbar >= 10.5 and cbar < 11.0:
            x0 = 1.8
            x1 = 4.0
        elif cbar >= 11.0 and cbar < 11.5:
            x0 = 1.8
            x1 = 4.0
        elif cbar >= 11.5 and cbar < 12.25:
            x0 = 2.0
            x1 = 4.0
        elif cbar >= 12.25 and cbar < 13.804:
            x0 = 2.0
            x1 = 5.0
        else: 
            x0=0.326*cbar-1.5
            x1=5.0
        akbar=3.0
        abar=(cbar-2.0*np.log(10.0)*x0)/((x1-x0)**3)
    else:
        akbar=object.AKS[object.ngasn[0]]
        x0=object.X00[object.ngasn[0]]
        x1=object.X11[object.ngasn[0]]
        abar=object.AAA[object.ngasn[0]]
    dcor = 0.5*np.log10(object.torr*293.15/(760.0*(object.temp+object.ZeroCelsius)))
    x0 -= dcor
    x1 -= dcor
    afc = 2.0*np.log(10.0)
    for i in range(20000):
        bg=object.bet[i]*object.gam[i]
        x = np.log10(bg)
        if x < x0:
            object.den[i] = 0.0
        elif x >= x0 and x <= x1:
            object.den[i] = abar * np.e ** (akbar * np.log(x1 - x)) + afc * x - cbar
        else:
            object.den[i] = afc * x - cbar
