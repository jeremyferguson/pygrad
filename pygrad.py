import numpy as np

class PygradException(Exception):
    pass

class InfileFormatException(PygradException):
    pass

class Main():
    gas_dict = {1: 'cf4', 2: 'argon', 3: 'helium4', 4: 'helium3', 5: 'neon', 6: 'krypton', 7: 'xenon', 8: 'methane', 9: 'ethane', 10: 'propane', 11: 'isobutane', 12: 'co2', 14: 'h20', 15: 'oxygen', 16: 'nitrogen', 18: 'n20', 21: 'hydrogen', 30: 'sf6', 31: 'nh3', 34: 'ch3oh', 35: 'c2h5oh', 36: 'c3h7oh'}
    
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
        if self.temp <= -273.15:
            raise InfileFormatException('Absolute zero')
        if self.torr <= 0:
            raise InfileFormatException('Negative pressure')
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Run Degrad')
    parser.add_argument('infile',help='Filename containing initial parameters')
    parser.add_argument('outfile',nargs='?',help='Filename for the output of the file',default='DEGRAD.OUT')
    args = parser.parse_args()
    main = Main(args.infile)

