import numpy as np

#Checks the implicit Fortran type of S
def checktype(s):
    s = s.upper()
    if ord(s[0]) > ord('H') and ord(s[0]) < ord('O'):
        return np.int64
    return np.double

#Compare every element in an array to PLACE decimal points
def compArr(a,b,place):
    assert(len(a) == len(b))
    for i in range(len(a)):
        if type(a[i]) == np.ndarray:
            compArr(a[i],b[i])
        else:
            assert(abs(a[i] - b[i]) < .1 ** -place),str(i)

#Get the only key in a dictionary with one key
def getSingle(d):
    return list(d.keys())[0]
