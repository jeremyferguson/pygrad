cimport numpy as np
import math
from libc.stdlib cimport malloc, free
from libc.string cimport memset
from PyGasMix.Gasmix cimport Gasmix

cdef double drand48(double dummy)

cdef class Pygrad:
    #cpdef GetSimFunctions(self, BFieldMag, BFieldAngle, EnableThermalMotion)
    cpdef Start(self)

    cdef public:
        double BField_Angle
        '''This is the angle between the magnetic field and the electric field.'''
        double BField_Mag
        '''This is the magnitude of the magentic field.'''
        double Max_Electron_Energy
        '''This is the upper limit for the electron energy integration.'''
        double ElectronEnergyStep
        '''PyGrad does the electron energy integration in 20000 steps this variable has the difference in energy between each step.'''
        double ThermalCut
        '''PyGrad will track electrons until they fall below this threshold.'''
        double ThermalEnergy
        '''This indicates the amount of energy in the gas (it is equal to the Boltzman constant * absolute tempreture).'''
        double RhydbergConst
        '''This is Rydberg constant times hc in eV.'''
        double TemperatureCentigrade
        '''This is the tempreture in degrees Centigrade.'''
        double Pressure_Torr
        '''This is the pressure in Torr.'''
        int nDelta
        '''Number of delta electrons'''
        int imip
        '''Type of beam'''
        int BeamDirection
        '''direction of beam'''
        int icount
        '''icount'''
        double EnergyCut
        '''Energy cut for MIPS'''
        double ExcitationWeight
        '''Weight given to excitation in Fano factor calculation'''
        int kgas
        '''Which gas has beta decayed identity numbers'''
        int lgas
        '''Which element in kgas has beta decayed identity numbers'''
        int lcmp
        '''Whether compton scattering is included'''
        int lray
        '''Wether Rayleigh scattering is included'''
        int lpap
        '''Whether pair production is included'''
        int lbrm
        '''Whether Bremstrahlung is included'''
        int IECascade
        '''Whether to use parameterised or exact cascade'''
        int OutputVerbosity 
        '''How much output to write'''
        double DetectorEfficiency
        '''DetectorEfficiency'''
        double SmallNumber
        '''This constant is equal to 1e-20. Used to be a small constant.'''
        double Phi
        '''Angle between electric and magnetic field'''
        double Theta
        '''The lower limit of the electron energy integration.'''
        double InitialElectronEnergy
        '''The final electron energy'''
        double FinalElectronEnergy
        '''Electric field [V/cm].'''
        double EField
        '''Constant that is equal to (electron mass charge ratio) / 2 * 1e-9.'''
        double CONST1
        '''Constant that is equal to CONST1 * 1e-2.'''
        double CONST2
        '''Constant that is equal to sqrt(0.2 * (electron mass charge ratio)) * 1e-9.'''
        double CONST3
        '''Constant that is equal to CONST1 * ALOSCH * 1e-15'''
        double CONST4 
        '''Constant that is equal to CONST3 / 2.0'''
        double CONST5 
        '''Whether gas is isotropic or anistropic'''
        int WhichAngularModel

        '''Arrays used in the mixerc subroutines.'''
        dict mixercArrs
        '''Arrays used in the cascade data subroutines.'''
        dict cascdata
        int last 
        double tmax 
        int nout 

        #Constants defined in Setup
        double TwoPi 
        double PIR2 
        double ElectronCharge 
        double ElectronMass 
        double AMU 
        double MassOverChargeDivTen 
        double ALOSCH 
        double ZeroCelcius 
        double OneAtmosphere 
        double HBAR 
        double EMS 
        double VC 
        double RE 
        double ALPH 
        double EOVM 
        double EMASS2
        double CONST 
        double A0 
        double BBCONST
    
        #Density Effect Constants
        double EIAV[80]
        int nElectrons[80]
        double X00[80]
        double X11[80]
        double AKS[80]
        double AAA[80] 
        int ICFLG
        int IRFLG
        int IPFLG
        int IBFLG
        int LPEFLG
        double PresTempCor
        '''Variable used to calculate the correlation constant between the pressure and tempreture. PresTempCor=ABZERO*Pressure/(ATMOS*(ABZERO+TemperatureC)*100.0D0).'''
        long long Random_Seed
        '''Random number generator seed. Not used at the moment.'''
        long long NumberOfGases
        '''Number of gases in the mixture.'''
        long long EnergySteps
        '''Steps for the electron energy integration.'''
        int Enable_Penning
        '''Variable used to indicate the inclusion of penning effects. '''
        int GasIDs[6]
        '''Array used to store the number of the 6 gases in the mixture.'''
        double GasFractions[6]
        '''Array used to store the percentage of each gas in the mixture.'''
        int icoln[512]
        '''icoln'''
        int icolnn[60]
        '''icolnn'''
        Gasmix MixObject
        '''Gas mixer object'''
      
        #Cross section arrays
        double CrossSectionSum[20000]
        double IonizationCrossSection[6][20000]
        double InelasticCrossSectionPerGas[6][250][20000]
        double AttachmentSectionSum[20000]

        double TotalCrossSection[20000]
        double RelativeIonMinusAttachCrossSection[20000]
        double InelasticCrossSection[20000]
        double ElasticCrossSection[20000]
        
        int numExcitationsPerGas[6]

        double MoleculesPerCm3PerGas[6]
        '''Array used to calculate the number of molecules/cm^3 for each gas.'''
        double VMoleculesPerCm3PerGas[6]
        '''Array used to calculate the VAN for each gas.'''
        
        int Msum[10000]
        '''MSUM'''
        int Mcomp[10000]
        '''MCOMP'''
        int Mrayl[10000]
        '''MRAYL'''
        int Mpair[10000]
        '''MPAIR'''
        int Mphot[10000]
        '''MPHOT'''
        int Mvac[10000]
        '''MVAC'''
        double Time[300]
        '''TIME'''
        int Icoll[30]
        '''ICOLL'''
        int Icoln[512]
        '''ICOLN'''
        int Icolnn[60]
        '''ICOLNN'''
        double Tcfmax[10]
        '''TCFMAX'''
        int nbrem[6]
        '''NBREM'''
        double ebrtot[6]
        '''EBRTOT'''

        double Dx[100000]
        '''DX'''
        double Dy[100000]
        '''DY'''
        double Dz[100000]
        '''DZ'''
        double Dt[100000]
        '''DT'''
        double Dxy[100000]
        '''DXY'''
        double Dxyz[100000]
        '''DXYZ'''
        
        double E[20000]
        '''Energy ar each energy step.'''
        double SqrtEnergy[20000]
        '''The square root of each energy step.'''
        double Gamma[20000]
        '''Gamma for each step'''
        double Beta[20000]
        '''Beta for each step'''
        double Density[20000]
        '''Density for each step'''

        double Rmax1[100000]
        '''RMAX1'''
        double Tsum[100000]
        '''TSUM'''
        double Xneg[100000]
        '''XNEG'''
        double Yneg[100000]
        '''YNEG'''
        double Zneg[100000]
        '''ZNEG'''
        double Edelta[100000]
        '''EDELTA'''
        double Edelta2[100000]
        '''EDELTA2'''
        int nClusters[100000]
        '''NCL'''
        int nClustersExcitation[100000]
        '''NCLEXC'''

        int ClusterDistributionBins1[1000]
        '''Electron cluster distributions in bins of 1 '''
        int ClusterDistributionBins3[1000]
        '''Electron cluster distributions in bins of 3 '''
        int ClusterDistributionBins10[1000]
        '''Electron cluster distributions in bins of 10 '''
        int ClusterDistributionBins30[1000]
        '''Electron cluster distributions in bins of 30 '''
        int ClusterDistributionBins100[1000]
        '''Electron cluster distributions in bins of 100 '''
        int ClusterDistributionBins300[1000]
        '''Electron cluster distributions in bins of 300 '''

        double XAverage[100000]
        '''Average X values.'''
        double YAverage[100000]
        '''Average Y values.'''
        double ZAverage[100000]
        '''Average Z values.'''
        double TAverage[100000]
        '''Average T values.'''
        double XYAverage[100000]
        '''Average XY values.'''
        double XYZAverage[100000]
        '''Average XYZ values.'''

        double XMaxRange[100000]
        '''Max Range of X'''
        double YMaxRange[100000]
        '''Max Range of Y'''
        double ZMaxRange[100000]
        '''Max Range of Z'''
        double XYMaxRange[100000]
        '''Max Range of XY'''
        
        #Plot distribution bins.
        int RDistributionBins2[31]
        '''Distribution of r in 2 micron bins'''
        int RDistributionBins10[31]
        '''Distribution of r in 10 micron bins'''
        int RDistributionBins40[31]
        '''Distribution of r in 40 micron bins'''
        int RDistributionBins100[31]
        '''Distribution of r in 100 micron bins'''
        int RDistributionBins400[31]
        '''Distribution of r in 400 micron bins'''
        int RDistributionBins1000[31]
        '''Distribution of r in 1000 micron bins'''
        int RDistributionBins4000[31]
        '''Distribution of r in 4000 micron bins'''
        int RDistributionBins10000[31]
        '''Distribution of r in 10000 micron bins'''
        int RDistributionBins40000[31]
        '''Distribution of r in 40000 micron bins'''
        int RDistributionBins100000[31]
        '''Distribution of r in 100000 micron bins'''
        int XDistributionBins2[31]
        '''Distribution of x in 2 micron bins'''
        int XDistributionBins10[31]
        '''Distribution of x in 10 micron bins'''
        int XDistributionBins40[31]
        '''Distribution of x in 40 micron bins'''
        int XDistributionBins100[31]
        '''Distribution of x in 100 micron bins'''
        int XDistributionBins400[31]
        '''Distribution of x in 400 micron bins'''
        int XDistributionBins1000[31]
        '''Distribution of x in 1000 micron bins'''
        int XDistributionBins4000[31]
        '''Distribution of x in 4000 micron bins'''
        int XDistributionBins10000[31]
        '''Distribution of x in 10000 micron bins'''
        int XDistributionBins40000[31]
        '''Distribution of x in 40000 micron bins'''
        int XDistributionBins100000[31]
        '''Distribution of x in 100000 micron bins'''
        int YDistributionBins2[31]
        '''Distribution of y in 2 micron bins'''
        int YDistributionBins10[31]
        '''Distribution of y in 10 micron bins'''
        int YDistributionBins40[31]
        '''Distribution of y in 40 micron bins'''
        int YDistributionBins100[31]
        '''Distribution of y in 100 micron bins'''
        int YDistributionBins400[31]
        '''Distribution of y in 400 micron bins'''
        int YDistributionBins1000[31]
        '''Distribution of y in 1000 micron bins'''
        int YDistributionBins4000[31]
        '''Distribution of y in 4000 micron bins'''
        int YDistributionBins10000[31]
        '''Distribution of y in 10000 micron bins'''
        int YDistributionBins40000[31]
        '''Distribution of y in 40000 micron bins'''
        int YDistributionBins100000[31]
        '''Distribution of y in 100000 micron bins'''
        int ZDistributionBins2[31]
        '''Distribution of z in 2 micron bins'''
        int ZDistributionBins10[31]
        '''Distribution of z in 10 micron bins'''
        int ZDistributionBins40[31]
        '''Distribution of z in 40 micron bins'''
        int ZDistributionBins100[31]
        '''Distribution of z in 100 micron bins'''
        int ZDistributionBins400[31]
        '''Distribution of z in 400 micron bins'''
        int ZDistributionBins1000[31]
        '''Distribution of z in 1000 micron bins'''
        int ZDistributionBins4000[31]
        '''Distribution of z in 4000 micron bins'''
        int ZDistributionBins10000[31]
        '''Distribution of z in 10000 micron bins'''
        int ZDistributionBins40000[31]
        '''Distribution of z in 40000 micron bins'''
        int ZDistributionBins100000[31]
        '''Distribution of z in 100000 micron bins'''
        int ElectronDistributionBins1[31] 
        '''Distribution of electrons in 1eV bins'''
        int ElectronDistributionBins10[31] 
        '''Distribution of electrons in 10eV bins'''
        int ElectronDistributionBins100[31] 
        '''Distribution of electrons in 100eV bins'''

