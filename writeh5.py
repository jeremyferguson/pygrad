import os,re,h5py, numpy as np

pygrad_home = os.getenv('PYGRAD_HOME')
filepath = pygrad_home + '/degrad-3.9a.f'
array_keys = ['INIOC', 'ES', 'PRBSH', 'PRBSHBT', 'YPEK', 'XPEK', 'YPEL1', 'XPEL1', 'YPEL2', 'XPEL2', 'YPEL3', 'XPEL3', 'YPEM1', 'XPEM1', 'YPEM2', 'XPEM2', 'YPEM3', 'XPEM3', 'YPEM4', 'XPEM4', 'YPEM5', 'XPEM5', 'YPEN1', 'XPEN1', 'YPEN2', 'XPEN2', 'YPEN3', 'XPEN3', 'YPEN4', 'XPEN4', 'YPEN5', 'XPEN5', 'YPEO1', 'XPEO1', 'YPEO2', 'XPEO2', 'YPEO3', 'XPEO3', 'XENE', 'YRAY', 'YCOM', 'YPAP', 'FFR', 'FFC', 'IZ', 'AMZ','XCOMP']
def checktype(s):
    if ord(s[0]) > ord('H') and ord(s[0]) < ord('O'):
        return np.int64
    return np.double

def compArr(a,b):
    assert(len(a) == len(b))
    for i in range(len(a)):
        if type(a[i]) == np.ndarray:
            compArr(a[i],b[i])
        else:
            assert(abs(a[i] - b[i]) < 1e-9),str(i)

def compareGases(i,j):
    for key in arrays[i]:
        if key in arrays[j]:
            print(key)
            compArr(arrays[i][key],arrays[j][key])
with open(filepath, 'r') as f:
    text = f.read()
individual_arrs = ['(IZ)','(AMZ)']

#Isolates the text for each subroutine
find_subs = re.compile("SUBROUTINE CGAS([0-9]+)([\S\s]*?)END")
subs = re.findall(find_subs,text)

#Gets the elements of the individual arrays
find_arrays = re.compile("DATA ([A-Z0-9]+)((?:(?:[\/0-9][0-9\.,/*ED\-]+\/?)[\s]+)+)")
arrays = {}
find_mult = re.compile("([0-9]+)\*([0-9\.e\-]+)")

#Gets the dimensions of all the arrays 
find_dimensions = re.compile("DIMENSION ([A-Z0-9]+(?:[^\n]+\n[\s]+\/)+[^\n]+)")
dimensions = {}
find_dim = re.compile("([A-Z0-9]+)(\([0-9,]+\))")

#Gets the values for the one-dimensional arrays that get their values assigned individually
find_individual_end = r"\(([0-9])\)[\s]*=([0-9\.DE]+)"
find_individual_arrs = [re.compile(i + find_individual_end) for i in individual_arrs]

#Gets the values for the multi-dimensional array elements that get assigned values individually
find_multi_dim = re.compile("DATA ([A-Z0-9]+\([0-9,]+\)\/(?:[\S]+))\n")
find_multi_statement = re.compile("([A-Z]+)\(([0-9,]+)\)\/([0-9DE\.\+\-]+)")
for groups in subs:
    i = int(groups[0])
    print(i)
    
    #Store the dimensions of all arrays assigned in the subroutine
    dims = re.findall(find_dimensions,groups[1])
    dimensions[i] = {}
    for d in dims:
        dimpairs = re.findall(find_dim,d)
        for t in dimpairs:
            key = t[0]
            val = tuple(int(i) for i in t[1][1:-1].split(','))
            dimensions[i][key] = val

    #Store all the array data assigned in the subroutine     
    data = re.findall(find_arrays,groups[1])
    arrays[i] = {}
    for g in data:
        key = g[0]
        print(key)
        val = g[1].replace('E','e').replace('D','e').replace('\n','').replace(' ','').replace('/','')
        val = re.sub(find_mult, lambda match: str(int(match.groups()[0]) * [match.groups()[1]])[1:-1],val)
        val = [i.strip().replace("'",'') for i in val.split(',')]
        val = np.array([float(i) for i in val if i != ''])
        val = np.reshape(val,dimensions[i][key])
        arrays[i][key] = val

    #Store all the individually assigned one-dimensional array information
    for ind in find_individual_arrs:
        individuals = re.findall(ind,groups[1])
        if individuals:
            arrays[i][individuals[0][0]] = np.zeros(dimensions[i][individuals[0][0]],dtype = checktype(individuals[0][0]))
            for match in individuals:
                arrays[i][match[0]][int(match[1])-1] = checktype(match[0])(match[2])

    #Store all the individually assigned multi-dimensional array information
    multi_dims = re.findall(find_multi_dim,groups[1])
    for multi_dim in multi_dims:
        assignments = re.findall(find_multi_statement, multi_dim)
        for a in assignments:
            name = a[0]
            dtype = checktype(name)
            if name not in arrays[i]:
                arrays[i][name] = np.zeros(dimensions[i][name],dtype=dtype)
            curr = arrays[i][name]
            index_arr = a[1].split(',')
            for j in index_arr[:-1]:
                curr = curr[int(j) - 1]
            curr[int(index_arr[-1]) - 1] = dtype(a[2].replace('D','e').replace('E','e'))

