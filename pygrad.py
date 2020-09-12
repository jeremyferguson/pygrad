import numpy as np,math,random
import cascdata,mixerc,utils,mixer

#Dictionary that holds information about all gases available in Pygrad.
gas_dict = {
        1: {'name':'CF4','formula':[{'C':1},{'F':4}]},
        2: {'name':'Argon','formula':[{'Ar':1}]},
        3: {'name':'Helium-4','formula':[{'He-4':1}]},
        4: {'name':'Helium-3','formula':[{'He-3':1}]}, 
        5: {'name':'Neon','formula':[{'Ne':1}]},
        6: {'name':'Krypton','formula':[{'Kr':1}]},
        7: {'name':'Xenon','formula':[{'Xe':1}]},
        8: {'name':'Methane','formula':[{'C':1},{'H':4}]},
        9: {'name':'Ethane','formula':[{'C':2},{'H':6}]}, 
        10: {'name':'Propane','formula':[{'C':3},{'H':8}]},
        11: {'name':'Isobutane','formula':[{'C':4},{'H':10}]},
        12: {'name':'CO2','formula':[{'C':1},{'O':2}]},
        14: {'name':'H2O','formula':[{'H':2},{'O':1}]},
        15: {'name':'Oxygen','formula':[{'O':2}]},
        16: {'name':'Nitrogen','formula':[{'N':2}]},
        18: {'name':'Nitrous Oxide','formula':[{'N':2},{'O':1}]},
        21: {'name':'Hydrogen','formula':[{'H':2}]},
        30: {'name':'SF6','formula':[{'S':1},{'F':6}]},
        31: {'name':'NH3','formula':[{'N':1},{'H':3}]},
        34: {'name':'CH3OH','formula':[{'C':1},{'H':4},{'O':1}]},
        35: {'name':'C2H5OH','formula':[{'C':2},{'H':6},{'O':1}]},
        36: {'name':'Iso-Propanol','formula':[{'C':3},{'H':8},{'O':1}]},
        44: {'name':'TMA','formula':[{'C':3},{'N':1},{'H':9}]},
        46: {'name':'N-Propanol','formula':[{'C':3},{'H':8},{'O':1}]}}
shell_order = ['K','L1','L2','L3','M1','M2','M3','M4','M5','N1','N2','N3','N4','N5','O1','O2','O3']

#Density Effect Constants
EIAV = np.array([115.0,188.0,41.8,41.8,137.0,352.0,482.0,41.7,45.4,47.1,48.3,85.0,0.0,71.6,95.0,82.0,0.0,84.9,0.0,0.0,19.2,19.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,128.0,53.7,0.0,0.0,67.6,62.9,61.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,48.3,0.0,61.1,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
JELEC = np.array([42,18,2,2,10,36,54,10,18,26,34,22,0,10,16,14,0,22,0,0,2,2,0,0,0,0,0,0,0,70,10,0,0,18,26,34,0,0,0,0,0,0,0,34,0,34,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
X00 = np.array([1.70,1.7635,2.2017,2.2017,2.0735,1.7158,1.5630,1.6263,1.5090,1.4339,1.3788,1.6294,0.0,1.7952,1.7541,1.7378,0.0,1.6477,0.0,0.0,1.8639,1.8639,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.6,1.6822,0.0,0.0,0.2529,0.2218,0.2046,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.3788,0.0,0.2046,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
X11 = np.array([4.00,4.4855,3.6122,3.6122,4.6421,5.0748,4.7371,3.9716,3.8726,3.8011,3.7524,4.1825,0.0,4.3437,4.3213,4.1323,0.0,4.1565,0.0,0.0,3.2718,3.2718,0.0,0.0,0.0,0.0,0.0,0.0,0.0,4.0,4.1158,0.0,0.0,2.7639,2.7052,2.6681,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.7524,0.0,2.6681,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
AKS = np.array([3.00,2.9618,5.8347,5.8347,3.5771,3.4051,2.7414,3.6257,3.6095,3.5920,3.4884,3.3227,0.0,3.5901,3.2913,3.2125,0.0,3.3318,0.0,0.0,5.7273,5.7273,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,3.6464,0.0,0.0,3.5477,3.4834,3.5415,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.4884,0.0,3.5415,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
AAA = np.array([.18551,.19714,.13443,.13443,.08064,.07446,.23314,.09253,0.09627,0.09916,.10852,.11768,0.0,.08101,.11778,.15349,0.0,.11992,0.0,0.0,.14092,.14092,0.0,0.0,0.0,0.0,0.0,0.0,0.0,.177484,.08315,0.0,0.0,.08970,0.09878,0.09644,0.0,0.0,0.0,0.0,0.0,0.0,0.0,.10852,0.0,0.09644,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#Mixer constants

CONST=1.873884e-20
EMASS2=1021997.804
A0=0.52917720859e-8
RY=13.60569193
BBCONST=16.0*API*A0*A0*RY*RY/EMASS2

def glob():
    return globals()

class PygradException(Exception):
    pass

class InfileFormatException(PygradException):
    pass

class Main():
    NBREM = np.zeros(6,dtype=int)
    EBRTOT = np.zeros(6,dtype=float)
    ICFLG=0
    IRFLG=0
    IPFLG=0
    IBFLG=0
    LPEFLG=0

    #Initializes the main program object.  Analagous to the SETUP subroutine in Degrad.
    def __init__(self, infile):
        self.last=0
        self.tmax=100.0  
        self.nout=10  
        
        #Set up initial angles for beam and electric field
        self.drzinit = math.cos(self.theta)
        self.drxinit = math.sin(self.theta) * math.cos(self.phi)
        self.dryinit = math.sin(self.theta) * math.sin(self.phi)
        
        self.initializeArrays()
        
        self.akt = (ABZERO+self.temp)*BOLTZ
        self.ans = [i * self.corr * ALOSCH for i in self.nfrac]
        self.an = 100.0*self.corr*ALOSCH                                            
        self.vans = [i * VC for i in self.ans]
        self.van= self.an * VC

        #Used to assign indices in the code blocks below
        def assignIndex(i,aj):
            self.e[i] = self.ehalf + self.estep * aj
            self.eroot[i] = math.sqrt(self.e[i])
            self.gam[i] = (EMS + self.e[i])/EMS
            self.bet[i] = math.sqrt(1.0-1.0/(self.gam[i] ** 2))
        
        if self.efinal <= 20000.0:
            self.estep = self.efinal/self.nstep
            self.ehalf = self.estep/2.0
            for i in range(20000):
                assignIndex(i,i)
        elif self.efinal <= 140000.0:
            self.estep = 1.0
            self.ehalf = 0.5
            for i in range(16000):
                assignIndex(i,i)
            estep1 = self.estep
            self.estep = (self.efinal-16000.0)/4000.0
            for i in range(16000,20000):
                assignIndex(i,i-15999)
            self.estep = estep1
        else:
            self.estep = 1.0
            self.ehalf = 0.5
            for i in range(12000):
                assignIndex(i,i)
            estep1 = self.estep
            self.estep = 20.0
            for i in range(12000,16000):
                assignIndex(i,i-11999)
            self.estep = (self.efinal-92000.0)/4000.0
            for i in range(16000,20000):
                assignIndex(i,i-15999)
            self.estep = estep1
        self.wb = AWB * self.bmag * 1e-12
        if self.bmag != 0.0:
            self.eovb = self.efield * 1e-9 / self.bmag


    #Checks all the parameters are in a valid range

    #Initialize all the arrays used in the program to zeros
    def initializeArrays(self):

    #Calculate density effect. DENSITY subroutine in Degrad
    def calcDensity(self):
        jelecMult = lambda x: x[1] * JELEC[x[0] - 1]
        zipped = zip(self.ngasn[:self.ngas],self.nfrac[:self.ngas])
        zippedAn = zip(self.ngasn[:self.ngas],self.ans[:self.ngas])
        sum1 = sum(map(lambda x: jelecMult(x) * math.log(EIAV[x[0] - 1]), list(zipped) ))
        zipped = zip(self.ngasn[:self.ngas],self.nfrac[:self.ngas])
        sumdnom = sum(map(jelecMult, list(zipped)))
        hsum = sum(map(jelecMult, zippedAn))
        eibar = math.e ** (sum1/sumdnom)
        hwp1 = math.sqrt(4.0 * API * hsum * RE ** 3) * ALPH * EMS
        delden = math.log(eibar/hwp1)
        cbar = 1.0 + 2.0 * delden
        if self.ngas != 1:
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
            abar=(cbar-2.0*math.log(10.0)*x0)/((x1-x0)**3)
        else:
            akbar=AKS[self.ngasn[0]]
            x0=X00[self.ngasn[0]]
            x1=X11[self.ngasn[0]]
            abar=AAA[self.ngasn[0]]
        dcor = 0.5*math.log10(self.torr*293.15/(760.0*(self.temp+ABZERO)))
        x0 -= dcor
        x1 -= dcor
        afc = 2.0*math.log(10.0)
        for i in range(20000):
            bg=self.bet[i]*self.gam[i]
            x = math.log10(bg)
            if x < x0:
                self.den[i] = 0.0
            elif x >= x0 and x <= x1:
                self.den[i] = abar * math.e ** (akbar * math.log(x1 - x)) + afc * x - cbar
            else:
                self.den[i] = afc * x - cbar

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Run Degrad')
    parser.add_argument('infile',help='Filename containing initial parameters')
    parser.add_argument('outfile',nargs='?',help='Filename for the output of the file',default='DEGRAD.OUT')
    args = parser.parse_args()
    main = Main(args.infile)
    main.calcDensity()
    cmix = mixerc.MixerC(main)
    cmix.mixc()
    mix = mixer.Mixer(main)
    mix.mix()
