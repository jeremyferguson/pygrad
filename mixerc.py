import numpy as np,pygrad, json

class MixerC():
    def __init__(self,main):
        self.main = main
        self.arrays = {
                "PRSH":np.zeros((6,3,17,17), dtype = float),
                "ESH":np.zeros((6,3,17), dtype = float), 
                "AUG":np.zeros((6,3,17,17,17), dtype = float), 
                "RAD":np.zeros((6,3,17,17), dtype = float),
                "PRSHBT":np.zeros((6,3,17), dtype = float),
                "IZ":np.zeros((6,3),dtype = int),
                "INIOCC":np.zeros((6,3,17), dtype = int),
                "ISHLMX":np.zeros((6,3), dtype = int),
                "AMZ":np.zeros((6,3),dtype = float),
                "XPE":np.zeros((6,3,17,60), dtype = float),
                "YPE":np.zeros((6,3,17,60), dtype = float),
                "XCP":np.zeros((6,3,54), dtype = float),
                "YRY":np.zeros((6,3,54), dtype = float),
                "YCP":np.zeros((6,3,45), dtype = float),
                "YPP":np.zeros((6,3,54), dtype = float),
                "FRMFR":np.zeros((6,3,45), dtype = float),
                "FRMFC":np.zeros((6,3,45), dtype = float)}

    def mixc(self):
        i = 0
        for gas in self.main.ngasn:
            self.gasmixc(gas,i)
            i += 1

    def gasmixc(self,gas,i):
        if gas == 0:
            return
        if gas not in pygrad.gas_dict:
            raise PygradException('Invalid gas number: '+ str(gas))
        
