import numpy as np,os,pygrad, h5py, utils

#Pairs each array found in element data to the higher-dimensional array 
#in the mixerC class.
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

#This class takes the place of the MIXERC subroutine.  It loads element
#data on each element in each gas from an hdf5 file.'''
cpdef MixerC(Pygrad object):
    #Collection of all arrays which store data loaded in for each gas.
    #Arrays are all 6x3 to store data from all 6 gases in the mixture, and 
    #up to 3 elements per gas. Arrays with 17 as a dimension are storing data per shell.
    cdef dict arrays = {
            
        #Probability of shell shakeoff from one shell to another.
        "PRSH":np.zeros((6,3,17,17), dtype = float),
        
        #Energy of shakeoff for each shell.
        "ESH":np.zeros((6,3,17), dtype = float),
        
        #Auger and Coster-Kronig transition rates for each shell.
        #For each element, K and L shells are stored in milliatomic
        #units and the rest of the shells are in 10**-4 atomic units 
        #and are converted to eV upon being loaded in.
        "AUG":np.zeros((6,3,17,17,17), dtype = float),
        
        #Radiative transition rates. Data is stored in units of
        #1.519e15/sec and are converted to eV upon being loaded in.
        "RAD":np.zeros((6,3,17,17), dtype = float),
        
        #Probability of shakeoff from beta decay. Stored as a
        #percentage but converted to a probability upon being loaded in.
        "PRSHBT":np.zeros((6,3,17), dtype = float),
        
        #Atomic number.
        "IZ":np.zeros((6,3),dtype = int),
        
        #Level occupancy for ground state.
        "INIOCC":np.zeros((6,3,17), dtype = int),
        
        #Highest occupied shell for each element.
        "ISHLMX":np.zeros((6,3), dtype = int),
        
        #Mass of each element multiplied by the number of atoms in the 
        #molecule.
        "AMZ":np.zeros((6,3),dtype = float),
        
        #Photoelectric absorption cross sections for each shell. 
        #Converted to logarithm upon being loaded in.
        "XPE":np.zeros((6,3,17,60), dtype = float),
        
        #Photoelectric absorption cross sections for each shell.
        #Converted to logarithm and multiplied by number of atoms in the
        #molecule * 1e-24 upon being loaded in.
        "YPE":np.zeros((6,3,17,60), dtype = float),
        
        #Compton cross-sections. Converted to logarithm upon being
        #loaded in.
        "XCP":np.zeros((6,3,54), dtype = float),
        
        #Rayleigh cross-sections. Converted to logarithm upon being 
        #loaded in.
        "YRY":np.zeros((6,3,54), dtype = float),
        
        #Compton cross-sections. Converted to logarithm upon being
        #loaded in.
        "YCP":np.zeros((6,3,54), dtype = float),
        
        #Pair production cross-sections. Converted to logarithm upon
        #being loaded in.
        "YPP":np.zeros((6,3,54), dtype = float),
        
        #Rayleigh form factors.
        "FRMFR":np.zeros((6,3,45), dtype = float),
        
        #Compton form factors.
        "FRMFC":np.zeros((6,3,45), dtype = float)}

    cdef dict cgas_data
    with open(os.getenv('PYGRAD_HOME')+'/cgas_data.npy','r') as f:
        cgas_data = pickle.load(f)

    i = 0
    for gas in object.GasIDs:
        gasmixc(object,gas,i,cgas_data)
        i += 1
        
#Read all the gas data for GAS in at mixture position I and do some conversions
#to the proper units.
cpdef gasmixc(Pygrad object, int gas,int i,dict cgas_data):
    group = self.f['elements']
    if gas == 0:
        return
    if gas not in pygrad.gas_dict:
        raise PygradException('Invalid gas number: '+ str(gas))
    #print(gas)
    formula = pygrad.gas_dict[gas]['formula']
    j = 0
    for pair in formula:
        element = utils.getSingle(pair)
        number = pair[element]
        elgroup = group[element]
        for array in elgroup['mixerc']:
            self.assignArray(i,j,array,elgroup['mixerc'][array],number)
        k = 1 
        for occ in self.arrays['INIOCC'][i][j]:
            if occ > 0.0:
                self.arrays['ISHLMX'][i][j] = k
            k += 1
        attrs = elgroup['mixerc'].attrs
        self.arrays['IZ'][i,j] = attrs.get('IZ')
        self.arrays['AMZ'][i,j] = attrs.get('AMZ') * number
        self.arrays['PRSH'][i,j] /= 100.0
        self.arrays['PRSHBT'][i,j] /= 100.0
        self.arrays['PRSH'][i,j] = self.arrays['PRSH'][i][j].T
        self.arrays['AUG'][i,j,:4] *= 0.0272105
        self.arrays['AUG'][i,j,4:] *= 0.00272105
        self.arrays['RAD'][i,j,4:,5:] *= 6.582119e-16
        self.arrays['YRY'][i,j] = np.log(self.arrays['YRY'][i][j] * number * 1e-24)
        self.arrays['YCP'][i,j] = np.log(self.arrays['YCP'][i][j] * number * 1e-24)
        self.arrays['YPP'][i,j] = np.log(self.arrays['YPP'][i][j] * number * 1e-24)
        self.arrays['XCP'][i,j] = np.log(self.arrays['XCP'][i][j])
        j += 1

#Assign ARRAY with NAME as its name to its proper position in the higher-dimensional
#array at position I,J, with NUMBER as the number of times it appears in the gas molecule.
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
