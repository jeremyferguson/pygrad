import numpy as np
import os

j = 0
'''Checks the implicit Fortran type of S.'''
def checktype(s):
    s = s.upper()
    if ord(s[0]) > ord('H') and ord(s[0]) < ord('O'):
        return np.int64
    return np.double

'''Compare every element in an array to PLACE decimal points.'''
def compArr(a,b,place,indices):
    assert(len(a) == len(b)),"{} {}".format(len(a),len(b))
    for i in range(len(a)):
        new_indices = indices + [i]
        #print(new_indices)
        if type(a[i]) == np.ndarray:
            compArr(a[i],b[i],place,new_indices)
        else:
            assert(abs(a[i] - b[i]) < .1 ** -place),"{} at position {} is diferent from {}".format(a[i],new_indices,b[i])

'''Get the only key in a dictionary with one key.'''
def getSingle(d):
    return list(d.keys())[0]

'''Runs PYCODE on a variety of inputs found in PATH and checks all the TESTVARS 
against the Fortran outputs found in subdirectories of PATH with OUTFNAME 
as the file name, using decimal PLACES for checking the values'''
def checkFortranTests(path,pycode,outfname,places):
    for fname in os.listdir(path):
        if fname[-3:] == '.in':
            print(fname)
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
            compArr(testvars, new_fortran,places,[])
            
