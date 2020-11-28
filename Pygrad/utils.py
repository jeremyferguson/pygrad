import numpy as np
import os

#Checks the implicit Fortran type of S.
def checktype(s):
    s = s.upper()
    if ord(s[0]) > ord('H') and ord(s[0]) < ord('O'):
        return np.int64
    return np.double

#Converts X to either a float or an int.
def convert(x):
    if '.' in x or 'e' in x:
        return float(x)
    return int(x)

#Compare every element in arrays A and B to PLACE decimal places.
def compArr(a,b,place=5):
    compArrFunc(a,b,lambda x,y,place: abs(x - y) < .1 ** -place,[],place)

#Assert that FUNC(x,y,PLACE) is true for all elements in A and B,using INDICES to keep track of the previous dimensions
def compArrFunc(a,b,func,indices,place=5):
    assert(len(a) == len(b)),"{} {}".format(len(a),len(b))
    for i in range(len(a)):
        new_indices = indices + [i]
        if type(a[i]) == np.ndarray:
            compArrFunc(a[i],b[i],func,new_indices,place)
        else:
            assert(func(a[i],b[i],place)),"Failed at comparing {} to {} at position {}".format(a[i],b[i],new_indices)

#Ensure every element in A is greater than every element in B to PLACE decimal places.
def compArrGreater(a,b,place=5):
    compArrFunc(a,b,lambda x,y,place:x > y + .1 ** -place,[],place)

#Ensure every element in A is greater than the single value B to PLACE decimal places.
def compArrGreaterValue(a,b,place=5):
    compArrFunc(a,b,lambda x,y,place:x > b + .1 ** -place,[],place)

#Ensure every element in A is less than every element in B to PLACE decimal places.
def compArrLess(a,b,place=5):
    compArrFunc(a,b,lambda x,y,place:x < y + .1 ** -place,[],place)

#Ensure every element in A is less than the single value B to PLACE decimal places.
def compArrLessValue(a,b,place=5):
    compArrFunc(a,b,lambda x,y,place:x < b + .1 ** -place,[],place)

#Get the only key in a dictionary with one key.
def getSingle(d):
    return list(d.keys())[0]

#Return a Pygrad object set up with the parameters found in FNAME at PATH
def createObject(path,fname):
    obj = Pygrad()
    with open(path+fname,'r') as f:
        text = f.read()
    params = [line.split(',') for line in text.split('\n')]
    obj.NumberOfGases = int(params[0][0])
    obj.nDelta = int(params[0][1])
    obj.imip = int(params[0][2])
    obj.BeamDirection = int(params[0][3])
    obj.Random_Seed = int(params[0][4])
    obj.InitialElectronEnergy = double(params[0][5])
    obj.ThermalCut = double(params[0][6])
    obj.EnergyCut = double(params[0][7])
    obj.GasIds = [int(i) for i in params[1]]
    obj.GasFractions = [double(i) for i in params[2][:6]]
    obj.TemperatureCentigrade = double(params[2][6])
    obj.Pressure_Torr = double(params[2][7])
    obj.EField = double(params[3][0])
    obj.BField_Mag = double(params[3][1])
    obj.BField_Angle = double(params[3][2])
    obj.OutputVerbosity = int(params[3][3])
    obj.Enable_Penning = int(params[3][4])
    obj.DetectorEfficiency = double(params[4][0])
    obj.ExcitationWeight = double(params[4][1])
    obj.kgas = int(params[4][2])
    obj.lgas = int(params[4][3])
    obj.lcmp = int(params[4][4])
    obj.lray = int(params[4][5])
    obj.lpap = int(params[4][6])
    obj.lbrm = int(params[4][7])
    obj.iecasc = int(params[4][8])
    return obj


#Runs PYCODE on a variety of inputs found in PATH and checks all the TESTVARS 
#against the Fortran outputs found in subdirectories of PATH with OUTFNAME 
#as the file name, checking to PLACES decimal places for equality.
def checkFortranTests(path,pycode,outfname,places):
    for fname in os.listdir(path):
        if fname[-3:] == '.in':
            #print(fname)
            prefix = fname[:-3]
            testvars = np.array(pycode(fname))
            outpath = path + '/' + prefix + '/' + outfname + '.out'
            with open(outpath, 'r') as outfile:
                text = outfile.read()
            fortran_out = np.array([[float(i) for i in line.split(',')] for line in text.split('\n')]).T
            i = 0
            new_fortran = np.zeros(testvars.shape)
            for testvar in testvars:
                new_fortran[i] = np.reshape(fortran_out[i],testvar.shape,'C')
                i += 1
            compArr(testvars, new_fortran,places)

