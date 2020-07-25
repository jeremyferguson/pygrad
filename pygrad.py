import numpy as np,math,random

class PygradException(Exception):
    pass

class InfileFormatException(PygradException):
    pass

class Main():
    gas_dict = {1: 'cf4', 2: 'argon', 3: 'helium4', 4: 'helium3', 5: 'neon', 6: 'krypton', 7: 'xenon', 8: 'methane', 9: 'ethane', 10: 'propane', 11: 'isobutane', 12: 'co2', 14: 'h20', 15: 'oxygen', 16: 'nitrogen', 18: 'n20', 21: 'hydrogen', 30: 'sf6', 31: 'nh3', 34: 'ch3oh', 35: 'c2h5oh', 36: 'c3h7oh'}
    API = math.pi
    ABZERO = 273.15
    ARY = 13.60569253
    PIR2 = 8.7973554297e-17
    ECHARG = 1.602176565e-19 
    EMASS = 9.10938291e-31                     
    EMS = 510998.928
    VC = 299792458.0
    AMU = 1.660538921e-27                                             
    BOLTZ = 8.6173324e-5     
    BOLTZJ = 1.3806488e-23                                              
    AWB = 1.758820088e10                                              
    ALOSCH = 2.6867805e19      
    RE = 2.8179403267e-13    
    ALPH = 137.035999074
    HBAR = 6.58211928e-16                                     
    EOVM = math.sqrt(2.0*ECHARG/EMASS)*100.0                            
    ATMOS = 760.0                                                     
    CONST1 = (AWB/2.0)*1e-19                                          
    CONST2 = CONST1*1e-2                                             
    CONST3 = (0.2*AWB)*1e-9                                   
    CONST4 = CONST3*ALOSCH*1e-15                                      
    CONST5 = CONST3/2.0
    TWOPI = 2.0*API
    NANISO = 2
    NBREM = np.zeros(6,dtype=int)
    EBRTOT = np.zeros(6,dtype=float)
    ICFLG=0
    IRFLG=0
    IPFLG=0
    IBFLG=0
    LPEFLG=0
    def __init__(self, infile):
        with open(infile, 'r') as f:
            text = f.read()
        lines = text.split('\n')[:-1]
        if len(lines) < 5:
            raise InfileFormatException('Not enough lines in input file')
        parameters = [line.split(',') for line in lines[:5]]
        lengths = [8,6,8,5,9]
        for line in zip(lengths,parameters):
            if len(line[1]) != line[0]:
                raise InfileFormatException('Incorrect number of parameters in input file')
        self.ngas, self.ndelta, self.imip, self.ndvec, self.nseed = [int(i) for i in parameters[0][:5]]
        self.icount = 0
        if self.imip == 1:
            self.icount = 1
        self.estart, self.etherm, self.ecut = [float(i.replace('−','-')) for i in parameters[0][5:]]
        self.ngasn = [int(i) for i in parameters[1]]
        self.nfrac = [float(i.replace('−','-')) for i in parameters[2][:6]]
        self.temp,self.torr = [float(i.replace('−','-')) for i in parameters[2][6:]]
        self.efield, self.bmag, self.btheta = [float(i) for i in parameters[3][:3]]
        self.iwrite, self.ipen = [int(i) for i in parameters[3][3:]]
        self.deteff, self.excweight = [float(i.replace('−','-')) for i in parameters[4][:2]]
        self.kgas, self.lgas, self.lcmp, self.lray, self.lpap, self.lbrm, self.iecasc = [int(i) for i in parameters[4][2:]] 
        self.check_parameters()
        if self.nseed == 0:
            self.nseed = 54217137
        random.seed(self.nseed)
        self.ebig=0.05*self.estart/1000.0 
        self.efinal=self.estart*1.0001+760.0*self.ebig/self.torr*(self.temp+ABZERO)/293.15*self.efield
        if self.efinal < 1.01 * self.estart :
            self.efinal = 1.01 * self.estart
        self.last=0
        self.tmax=100.0  
        self.nout=10  
        self.nstep=20000
        
        #Set up initial angles for beam and electric field
        if ndvec == 1:
            self.phi = 0.0                          
            self.theta = 0.0
        elif ndvec == -1:
            self.phi = 0
            self.theta = self.API
        elif ndvec == 0:
            self.phi = 0.0
            self.theta = API/2.0
        elif ndvec == 2:
            self.R3 = random.random()
            self.phi = TWOPI * self.R3
            self.R4 = random.random()
            self.theta = math.acos(1-2*self.R4)
        self.drzinit = math.cos(self.theta)
        self.drxinit = math.sin(self.theta) * math.cos(self.phi)
        self.dryinit = math.sin(self.theta) * math.sin(self.phi)
        
        self.initializeArrays()
        
        self.corr = self.ABZERO*self.torr/(self.ATMOS*(self.ABZERO+self.temp)*100.0)
        self.akt = (self.ABZERO+self.temp)*self.BOLTZ
        self.ans = [i * self.corr * self.ALOSCH for i in self.frac]
        self.an = 100.0*self.corr*self.ALOSCH                                            
        self.vans = [i * self.VC for i in self.ans]
        self.van= self.an * self.VC

        if self.efinal <= 20000.0:
            self.estep = self.efinal/self.nstep
            self.ehalf = self.estep/2.0
            for i in range(20000):
                self.e[i] = self.ehalf + self.estep * i
                self.gam[i] = (self.EMS + self.e[i]
                self.bet[i] = math.sqrt(1.0-1.0/(self.gam[i] ** 2))
        
    #Check if a variable is in a valid range
    def check_var(self, var, name, low=0, high=1):
        if var < low or var > high:
            raise InfileFormatException('Invalid value for '+name)

    #Checks all the parameters are in a valid range
    def check_parameters(self):
        self.check_var(self.ngas,'NGAS',1,6)
        self.check_var(self.imip,'IMIP',1,5)
        self.check_var(self.ndvec,'NDVEC',-1,2)
        self.check_var(self.iwrite,'IWRITE',0,2)
        self.check_var(self.ipen,'IPEN')
        self.check_var(self.deteff,'DETEFF',0,100)
        self.check_var(self.lcmp,'LCMP')
        self.check_var(self.lray,'LRAY')
        self.check_var(self.lpap,'LPAP')
        self.check_var(self.lbrm,'LBRM')
        self.check_var(self.iecasc,'IECASC')
        if (self.imip == 1 and self.ndelta > 1e5) or (self.imip != 1 and self.ndelta > 1e4):
            raise InfileFormatException('Invalid value for NDELTA')
        if self.imip == 3 and self.estart > 3e6:
            raise InfileFormatException('X-ray energy is too high')
        if self.estart <= 0 or self.etherm <= 0 or self.ecut < 0:
            raise InfileFormatException('Negative energy')
        if abs(100.0 - sum(self.nfrac)) > 1e-6:
            raise InfileFormatException('NFRAC does not sum to 100%')
        actualngas = 0
        added = []
        for pair in zip(self.ngasn, self.nfrac):
            if pair[0] != 0 and pair[0] not in self.gas_dict:
                raise InfileFormatException('Invalid gas number: '+str(pair[0]))
            if (pair[0] == 0 and pair[1] != 0) or (pair[0] != 0 and pair[1] == 0):
                raise InfileFormatException('Mismatch between NGAS and NFRAC')
            if pair[1] < 0:
                raise InfileFormatException('Negative NFRAC')
            if (pair[0] != 0 and pair[1] != 0):
                if pair[0] in added:
                    raise InfileFormatException('Duplicate gas number')
                actualngas += 1
                added.append(pair[0])
        if actualngas != self.ngas:
            raise InfileFormatException('Mismatch between NGAS and NFRAC')
        if self.temp <= -self.ABZERO:
            raise InfileFormatException('Absolute zero')
        if self.torr <= 0:
            raise InfileFormatException('Negative pressure')

    #Initialize all the arrays used in the program to zeros
    #TODO: Clean this up
    def initializeArrays():
        self.msum = np.zeros(10000,dtype = int)
        self.mcomp = np.zeros(10000,dtype = int)
        self.mrayl = np.zeros(10000,dtype = int)
        self.mpair = np.zeros(10000,dtype = int)
        self.mphot = np.zeros(10000,dtype = int)
        self.mvac = np.zeros(10000,dtype = int)
        self.time = np.zeros(1300,dtype=float)
        self.icoll = np.zeros(30,dtype=int)
        self.icoln = np.zeros(512,dtype=int)
        self.icolnn = np.zeros(60,dtype=int)
        self.tcfmax = np.zeros(10,dtype=float)
        self.nxpl2 = np.zeros(31,dtype=int)
        self.nypl2 = np.zeros(31,dtype=int)
        self.nzpl2 = np.zeros(31,dtype=int)
        self.nxpl10 = np.zeros(31,dtype=int)
        self.nypl10 = np.zeros(31,dtype=int)
        self.nzpl10 = np.zeros(31,dtype=int)
        self.nxpl40 = np.zeros(31,dtype=int)
        self.nypl40 = np.zeros(31,dtype=int)
        self.nzpl40 = np.zeros(31,dtype=int)
        self.nxpl100 = np.zeros(31,dtype=int)
        self.nypl100 = np.zeros(31,dtype=int)
        self.nzpl100 = np.zeros(31,dtype=int)
        self.nxpl400 = np.zeros(31,dtype=int)
        self.nypl400 = np.zeros(31,dtype=int)
        self.nzpl400 = np.zeros(31,dtype=int)
        self.nxpl1000 = np.zeros(31,dtype=int)
        self.nypl1000 = np.zeros(31,dtype=int)
        self.nzpl1000 = np.zeros(31,dtype=int)
        self.nxpl4000 = np.zeros(31,dtype=int)
        self.nypl4000 = np.zeros(31,dtype=int)
        self.nzpl4000 = np.zeros(31,dtype=int)
        self.nxpl10000 = np.zeros(31,dtype=int)
        self.nypl10000 = np.zeros(31,dtype=int)
        self.nzpl10000 = np.zeros(31,dtype=int)
        self.nxpl40000 = np.zeros(31,dtype=int)
        self.nypl40000 = np.zeros(31,dtype=int)
        self.nzpl40000 = np.zeros(31,dtype=int)
        self.nxpl100000 = np.zeros(31,dtype=int)
        self.nypl100000 = np.zeros(31,dtype=int)
        self.nzpl100000 = np.zeros(31,dtype=int)
        self.nrpl2 = np.zeros(31,dtype=int)
        self.nrpl10 = np.zeros(31,dtype=int)
        self.nrpl40 = np.zeros(31,dtype=int)
        self.nrpl100 = np.zeros(31,dtype=int)
        self.nrpl400 = np.zeros(31,dtype=int)
        self.nrpl1000 = np.zeros(31,dtype=int)
        self.nrpl4000 = np.zeros(31,dtype=int)
        self.nrpl10000 = np.zeros(31,dtype=int)
        self.nrpl40000 = np.zeros(31,dtype=int)
        self.nrpl100000 = np.zeros(31,dtype=int)
        self.npl1 = np.zeros(100,dtype=int)
        self.npl10 = np.zeros(100,dtype=int)
        self.npl100 = np.zeros(100,dtype=int)
        melec = np.zeros(1000,dtype=int)
        melec3 = np.zeros(1000,dtype=int)
        melec10 = np.zeros(1000,dtype=int)
        melec30 = np.zeros(1000,dtype=int)
        melec100 = np.zeros(1000,dtype=int)
        melec300 = np.zeros(1000,dtype=int)
        xav = np.zeros(100000,dtype=float)
        yav = np.zeros(100000,dtype=float)
        zav = np.zeros(100000,dtype=float)
        tav = np.zeros(100000,dtype=float)
        xyav = np.zeros(100000,dtype=float)
        xyzav = np.zeros(100000,dtype=float)
        dx = np.zeros(100000,dtype=float)
        dy = np.zeros(100000,dtype=float)
        dz = np.zeros(100000,dtype=float)
        dt = np.zeros(100000,dtype=float)
        dxy = np.zeros(100000,dtype=float)
        dxyz = np.zeros(100000,dtype=float)
        farx1 = np.zeros(100000,dtype=float)
        fary1 = np.zeros(100000,dtype=float)
        farz1 = np.zeros(100000,dtype=float)
        farxy1 = np.zeros(100000,dtype=float)
        rmax1 = np.zeros(100000,dtype=float)
        tsum = np.zeros(100000,dtype=float)
        xneg = np.zeros(100000,dtype=float)
        yneg = np.zeros(100000,dtype=float)
        zneg = np.zeros(100000,dtype=float)
        edelta = np.zeros(100000,dtype=float)
        edelta2 = np.zeros(100000,dtype=float)
        ncl = np.zeros(100000,dtype=int)
        nclexc = np.zeros(100000,dtype=int)
        self.e = np.zeros(20000,dtype = float)
        self.gam = np.zeros(20000,dtype=float)
        self.bet = np.zeros(20000,dtype=float)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Run Degrad')
    parser.add_argument('infile',help='Filename containing initial parameters')
    parser.add_argument('outfile',nargs='?',help='Filename for the output of the file',default='DEGRAD.OUT')
    args = parser.parse_args()
    main = Main(args.infile)

