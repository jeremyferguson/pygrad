import numpy as np
import os

'''Checks the implicit Fortran type of S.'''
def checktype(s):
    s = s.upper()
    if ord(s[0]) > ord('H') and ord(s[0]) < ord('O'):
        return np.int64
    return np.double

'''Compare every element in an array to PLACE decimal points.'''
def compArr(a,b,place):
    assert(len(a) == len(b))
    for i in range(len(a)):
        if type(a[i]) == np.ndarray:
            compArr(a[i],b[i],place)
        else:
            assert(abs(a[i] - b[i]) < .1 ** -place),str(i)

'''Get the only key in a dictionary with one key.'''
def getSingle(d):
    return list(d.keys())[0]

'''Runs PYCODE on a variety of inputs found in PATH and checks all the TESTVARS 
against the Fortran outputs found in subdirectories of PATH with OUTFNAME 
as the file name, using decimal PLACES for checking the values'''
def checkFortranTests(path,pycode,outfname,places):
    for fname in os.listdir(path):
        if fname[-3:] == '.in':
            prefix = fname[:-3]
            testvars = np.array(pycode(fname))
            outpath = path + '/' + prefix + '/' + outfname + '.out'
            with open(outpath, 'r') as outfile:
                text = outfile.read()
            fortran_out = np.array([[float(i) for i in line.split(',')] for line in text.split('\n')]).T
            i = 0
            for testvar in testvars:
                fortran_out[i].reshape(testvar.shape)
                i += 1
            compArr(testvars, fortran_out,places)
            
