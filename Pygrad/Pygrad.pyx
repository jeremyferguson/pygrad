
import math

from libc.math cimport sin, cos, acos, asin, log, sqrt, pow
from libc.string cimport memset
import Setups
#import Mixers
#import EnergyLimits
#import PyBoltz.MonteFuncs
#import PyBoltz.Townsend
#cimport PyBoltz.MonteFuncs
#cimport PyBoltz.Townsend
#from PyBoltz.MonteFuncs cimport MONTE,MONTET,MONTEB,MONTEBT,MONTEC,MONTECT
#from PyBoltz.Townsend cimport ALPCALCT
from PyGasMix.Gasmix cimport Gasmix



cdef extern from "C/RM48.h":
    double DRAND48(double dummy)
    void RM48(double lenv)
    void RM48IN(int IJL, int NTOT, intNTOT2)


cdef double drand48(double dummy):
    return DRAND48(dummy)
cdef void setSeed(int seed):
    RM48IN(seed, 0, 0)
    return

cdef class Pygrad:
    """
    This is the main object used to start the simulation, as well as store the information of the simulation.
    It has most of the needed arrays, and variables.
    TODO:Add documentation link for Degrad
    .. note::
    """

    def __init__(self):
        '''
        Fill all the variables needed with zeros.This function uses memset as it is fast.
        '''
        memset(self.icolnn, 0, 60 * sizeof(double))
        memset(self.icoln, 0, 512 * sizeof(double))
        memset(self.GasIDs, 0, 6 * sizeof(double))
        memset(self.GasFractions, 0, 6 * sizeof(double))
        memset(self.E, 0, 20000 * sizeof(double))
        memset(self.ElasticCrossSection, 0, 20000 * sizeof(double))
        memset(self.CrossSectionSum, 0, 20000 * sizeof(double))
        memset(self.IonizationCrossSection, 0, 6 * 20000 * sizeof(double))
        memset(self.InelasticCrossSectionPerGas, 0, 6 * 250 * 20000 * sizeof(double))
        memset(self.AttachmentSectionSum, 0, 20000 * sizeof(double))
        memset(self.TotalCrossSection, 0, 20000 * sizeof(double))
        memset(self.RelativeIonMinusAttachCrossSection, 0, 20000 * sizeof(double))
        memset(self.InelasticCrossSection, 0, 20000 * sizeof(double))
        memset(self.ElasticCrosSection, 0, 20000 * sizeof(double))
        memset(self.SqrtEnergy, 0, 20000 * sizeof(double))
        memset(self.Msum, 0, 10000 * sizeof(int))
        memset(self.Mcomp, 0, 10000 * sizeof(int))
        memset(self.Mrayl, 0, 10000 * sizeof(int))
        memset(self.Mpair, 0, 10000 * sizeof(int))
        memset(self.Mphot, 0, 10000 * sizeof(int))
        memset(self.Mvac, 0, 10000 * sizeof(int))
        memset(self.Time, 0, 300 * sizeof(double))
        memset(self.Icoll, 0, 30 * sizeof(int))
        memset(self.Tcfmax, 0, 10 * sizeof(double))
        memset(self.nbrem, 0, 6 * sizeof(int))
        memset(self.ebrtot, 0, 6 * sizeof(double))
        memset(self.Dx, 0, 100000 * sizeof(double))
        memset(self.Dy, 0, 100000 * sizeof(double))
        memset(self.Dz, 0, 100000 * sizeof(double))
        memset(self.Dt, 0, 100000 * sizeof(double))
        memset(self.Dxy, 0, 100000 * sizeof(double))
        memset(self.Dxyz, 0, 100000 * sizeof(double))
        memset(self.Gamma, 0, 20000 * sizeof(double))
        memset(self.Beta, 0, 20000 * sizeof(double))
        memset(self.Density, 0, 20000 * sizeof(double))
        memset(self.Rmax1, 0, 100000 * sizeof(double))
        memset(self.Tsum, 0, 100000 * sizeof(double))
        memset(self.Xneg, 0, 100000 * sizeof(double))
        memset(self.Yneg, 0, 100000 * sizeof(double))
        memset(self.Zneg, 0, 100000 * sizeof(double))
        memset(self.Edelta, 0, 100000 * sizeof(double))
        memset(self.Edelta2, 0, 100000 * sizeof(double))
        memset(self.nClusters, 0, 100000 * sizeof(double))
        memset(self.nClustersExcitation, 0, 100000 * sizeof(double))
        memset(self.ClusterDistributionBins1, 0, 1000 * sizeof(double))
        memset(self.ClusterDistributionBins3, 0, 1000 * sizeof(double))
        memset(self.ClusterDistributionBins10, 0, 1000 * sizeof(double))
        memset(self.ClusterDistributionBins30, 0, 1000 * sizeof(double))
        memset(self.ClusterDistributionBins100, 0, 1000 * sizeof(double))
        memset(self.ClusterDistributionBins300, 0, 1000 * sizeof(double))
        memset(self.XAverage,0,100000 * sizeof(double))
        memset(self.YAverage,0,100000 * sizeof(double))
        memset(self.ZAverage,0,100000 * sizeof(double))
        memset(self.XYAverage,0,100000 * sizeof(double))
        memset(self.XYZAverage,0,100000 * sizeof(double))
        memset(self.TAverage,0,100000 * sizeof(double))
        memset(self.XMaxRange,0,100000 * sizeof(double))
        memset(self.YMaxRange,0,100000 * sizeof(double))
        memset(self.ZMaxRange,0,100000 * sizeof(double))
        memset(self.XYMaxRange,0,100000 * sizeof(double))
        memset(self.RDistributionBins2,0,31 * sizeof(int))
        memset(self.RDistributionBins10,0,31 * sizeof(int))
        memset(self.RDistributionBins40,0,31 * sizeof(int))
        memset(self.RDistributionBins100,0,31 * sizeof(int))
        memset(self.RDistributionBins400,0,31 * sizeof(int))
        memset(self.RDistributionBins1000,0,31 * sizeof(int))
        memset(self.RDistributionBins4000,0,31 * sizeof(int))
        memset(self.RDistributionBins10000,0,31 * sizeof(int))
        memset(self.RDistributionBins40000,0,31 * sizeof(int))
        memset(self.RDistributionBins100000,0,31 * sizeof(int))
        memset(self.XDistributionBins2,0,31 * sizeof(int))
        memset(self.XDistributionBins10,0,31 * sizeof(int))
        memset(self.XDistributionBins40,0,31 * sizeof(int))
        memset(self.XDistributionBins100,0,31 * sizeof(int))
        memset(self.XDistributionBins400,0,31 * sizeof(int))
        memset(self.XDistributionBins1000,0,31 * sizeof(int))
        memset(self.XDistributionBins4000,0,31 * sizeof(int))
        memset(self.XDistributionBins10000,0,31 * sizeof(int))
        memset(self.XDistributionBins40000,0,31 * sizeof(int))
        memset(self.XDistributionBins100000,0,31 * sizeof(int))
        memset(self.YDistributionBins2,0,31 * sizeof(int))
        memset(self.YDistributionBins10,0,31 * sizeof(int))
        memset(self.YDistributionBins40,0,31 * sizeof(int))
        memset(self.YDistributionBins100,0,31 * sizeof(int))
        memset(self.YDistributionBins400,0,31 * sizeof(int))
        memset(self.YDistributionBins1000,0,31 * sizeof(int))
        memset(self.YDistributionBins4000,0,31 * sizeof(int))
        memset(self.YDistributionBins10000,0,31 * sizeof(int))
        memset(self.YDistributionBins40000,0,31 * sizeof(int))
        memset(self.YDistributionBins100000,0,31 * sizeof(int))
        memset(self.ZDistributionBins2,0,31 * sizeof(int))
        memset(self.ZDistributionBins10,0,31 * sizeof(int))
        memset(self.ZDistributionBins40,0,31 * sizeof(int))
        memset(self.ZDistributionBins100,0,31 * sizeof(int))
        memset(self.ZDistributionBins400,0,31 * sizeof(int))
        memset(self.ZDistributionBins1000,0,31 * sizeof(int))
        memset(self.ZDistributionBins4000,0,31 * sizeof(int))
        memset(self.ZDistributionBins10000,0,31 * sizeof(int))
        memset(self.ZDistributionBins40000,0,31 * sizeof(int))
        memset(self.ZDistributionBins100000,0,31 * sizeof(int))
        memset(self.ElectronDistributionBins1,0,31 * sizeof(int))
        memset(self.ElectronDistributionBins10,0,31 * sizeof(int))
        memset(self.ElectronDistributionBins100,0,31 * sizeof(int))

        # Input parameters / settings
        self.BField_Angle = 0.0
        self.BField_Mag = 0.0
        self.NumberOfGases = 0
        self.TemperatureCentigrade = 0.0
        self.Pressure_Torr = 0.0
        self.Enable_Penning = 0
        self.EField = 0.0

        # Calculated Constants
        self.CONST1 = 0.0
        self.CONST2 = 0.0
        self.CONST3 = 0.0
        self.PIR2 = 0.0
        self.RhydbergConst = 0.0
        self.ThermalEnergy = 0.0
        self.ThermalCut = 0.0
        self.SmallNumber = 1e-20
        self.PresTempCor = 0.0

        # Dynamically set
        self.EnergySteps = 0
        self.AnisotropicDetected = 0
        self.Max_Electron_Energy = 0.0
        self.ElectronEnergyStep = 0
        self.InitialElectronEnergy = 0.0

        self.Random_Seed = 54217137
        self.MixObject = Gasmix()

    cpdef Start(self):
        setSeed(self.RandomSeed)
        Setups.setup(self)
