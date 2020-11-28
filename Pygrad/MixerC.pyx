import numpy as np,os,utils

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
#This class takes the place of the MIXERC subroutine.  It loads element
#data on each element in each gas from an hdf5 file.'''
cpdef MixerC(Pygrad object):
    #Collection of all arrays which store data loaded in for each gas.
    #Arrays are all 6x3 to store data from all 6 gases in the mixture, and 
    #up to 3 elements per gas. Arrays with 17 as a dimension are storing data per shell.
    object.mixercArrs = {
            
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
    print(cgas_data.keys())

    i = 0
    for gas in object.GasIDs:
        gasmixc(object,gas,i,cgas_data)
        i += 1
        
#Read all the gas data for GAS in at mixture position I and do some conversions
#to the proper units.
cpdef gasmixc(Pygrad object, int gas,int i,dict cgas_data):
    if gas == 0:
        return
    if gas not in gas_dict:
        raise PygradException('Invalid gas number: '+ str(gas))
    formula = gas_dict[gas]['formula']
    j = 0
    for pair in formula:
        element = utils.getSingle(pair)
        number = pair[element]
        elementData = cgas_data[element]
        for array:
            self.assignArray(i,j,array,elementData[array],number)
        k = 1 
        for occ in object.mixercArrs['INIOCC'][i][j]:
            if occ > 0.0:
                object.mixercArrs['ISHLMX'][i][j] = k
            k += 1
        attrs = elgroup['mixerc'].attrs
        object.mixercArrs['IZ'][i,j] = attrs.get('IZ')
        object.mixercArrs['AMZ'][i,j] = attrs.get('AMZ') * number
        object.mixercArrs['PRSH'][i,j] /= 100.0
        object.mixercArrs['PRSHBT'][i,j] /= 100.0
        object.mixercArrs['PRSH'][i,j] = object.mixercArrs['PRSH'][i][j].T
        object.mixercArrs['AUG'][i,j,:4] *= 0.0272105
        object.mixercArrs['AUG'][i,j,4:] *= 0.00272105
        object.mixercArrs['RAD'][i,j,4:,5:] *= 6.582119e-16
        object.mixercArrs['YRY'][i,j] = np.log(object.mixercArrs['YRY'][i][j] * number * 1e-24)
        object.mixercArrs['YCP'][i,j] = np.log(object.mixercArrs['YCP'][i][j] * number * 1e-24)
        object.mixercArrs['YPP'][i,j] = np.log(object.mixercArrs['YPP'][i][j] * number * 1e-24)
        object.mixercArrs['XCP'][i,j] = np.log(object.mixercArrs['XCP'][i][j])
        j += 1

#Assign ARRAY with NAME as its name to its proper position in the higher-dimensional
#array at position I,J, with NUMBER as the number of times it appears in the gas molecule.
cpdef assignArray(Pygrad object,int i,int j,string name,np.ndarray array,number):
    if name in array_pairs:
        key = array_pairs[name]
        object.mixercArrs[key][i][j] = array[()]
    else:
        key = name[:3]
        k = shell_order.index(name[3:])
        l = array.shape[0]
        object.mixercArrs[key][i][j][k][:l] = array[:]
        if key[0] == 'Y':
            object.mixercArrs[key][i][j][k][:l] *= (1e-24 * number)
        object.mixercArrs[key][i][j][k][:l] = np.log(object.mixercArrs[key][i][j][k][:l])
