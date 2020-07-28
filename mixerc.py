import numpy as np,pygrad, h5py, utils

array_pairs = {
        'PRBSH':'PRSH',
        'INIOC':'INIOCC',
        'ES':'ESH',
        'XCOM':'XCP',
        'XENE':'XCP',
        'YRAY':'YRY',
        'YPAP':'YPP',
        'FFR':'FRMFR',
        'FFC':'FRMFC',
        'PRBSHBT':'PRSHBT',
        'YCOM':'YCP',
        'A':'AUG',
        'R':'RAD'}

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
                "YCP":np.zeros((6,3,54), dtype = float),
                "YPP":np.zeros((6,3,54), dtype = float),
                "FRMFR":np.zeros((6,3,45), dtype = float),
                "FRMFC":np.zeros((6,3,45), dtype = float)}
        self.f = h5py.File('gas_data.hdf5','r')

    def mixc(self):
        i = 0
        for gas in self.main.ngasn:
            self.gasmixc(gas,i)
            i += 1
        self.arrays['PRSH'] /= 100.0
        self.arrays['PRSHBT'] /= 100.0

    def gasmixc(self,gas,i):
        group = self.f['elements']['mixerc']
        if gas == 0:
            return
        if gas not in pygrad.gas_dict:
            raise PygradException('Invalid gas number: '+ str(gas))
        print(gas)
        formula = pygrad.gas_dict[gas]['formula']
        j = 0
        for pair in formula:
            element = utils.getSingle(pair)
            number = pair[element]
            elgroup = group[element]
            for array in elgroup:
                self.assignArray(i,j,array,elgroup[array],number)
            attrs = elgroup.attrs
            self.arrays['IZ'][i][j] = attrs.get('IZ')
            self.arrays['AMZ'][i][j] = attrs.get('AMZ') * number
            self.arrays['PRSH'][i][j] /= 100.0
            self.arrays['PRSHBT'][i][j] /= 100.0
            self.arrays['PRSH'][i][j] = self.arrays['PRSH'][i][j].T
            self.arrays['AUG'][i][j][:4] *= 0.0272105
            self.arrays['AUG'][i][j][4:] *= 0.00272105
            self.arrays['RAD'][i][j][4:][5:] *= 6.582119e-16
            self.arrays['YRY'][i][j] = np.log(self.arrays['YRY'][i][j] * number * 1e-24)
            self.arrays['YCP'][i][j] = np.log(self.arrays['YCP'][i][j] * number * 1e-24)
            self.arrays['YPP'][i][j] = np.log(self.arrays['YPP'][i][j] * number * 1e-24)
            self.arrays['XCP'][i][j] = np.log(self.arrays['XCP'][i][j])
            j += 1

    def assignArray(self,i,j,name,array,number):
        if name in array_pairs:
            key = array_pairs[name]
            self.arrays[key][i][j] = array[()]
        else:
            key = name[:3]
            k = pygrad.shell_order.index(name[3:])
            l = array.shape[0]
            self.arrays[key][i][j][k][:l] = array[:]
            if key[0] == 'Y':
                self.arrays[key][i][j][k][:l] *= (1e-24 * number)
            self.arrays[key][i][j][k][:l] = np.log(self.arrays[key][i][j][k][:l])
