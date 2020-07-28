import os,re,h5py, numpy as np,copy
import pygrad,utils

pygrad_home = os.getenv('PYGRAD_HOME')
filepath = pygrad_home + '/degrad-3.9a.f'
all_array_keys = ['A','R','INIOC', 'ES', 'PRBSH', 'PRBSHBT', 'YPEK', 'XPEK', 'YPEL1', 'XPEL1', 'YPEL2', 'XPEL2', 'YPEL3', 'XPEL3', 'YPEM1', 'XPEM1', 'YPEM2', 'XPEM2', 'YPEM3', 'XPEM3', 'YPEM4', 'XPEM4', 'YPEM5', 'XPEM5', 'YPEN1', 'XPEN1', 'YPEN2', 'XPEN2', 'YPEN3', 'XPEN3', 'YPEN4', 'XPEN4', 'YPEN5', 'XPEN5', 'YPEO1', 'XPEO1', 'YPEO2', 'XPEO2', 'YPEO3', 'XPEO3', 'XENE', 'YRAY', 'YCOM', 'YPAP', 'FFR', 'FFC', 'IZ', 'AMZ','XCOM']
gas_dict = copy.deepcopy(pygrad.gas_dict)
flat_arrs = ['IZ','AMZ']
find_subs = re.compile("SUBROUTINE CGAS([0-9]+)([\S\s]*?)END")
find_arrays = re.compile("DATA ([A-Z0-9]+)((?:(?:[\/0-9][0-9\.,/*ED\-]+\/?)[\s]+)+)")
find_array_elements = re.compile("([0-9]+)\*([0-9\.e\-]+)")
find_dimensions = re.compile("DIMENSION ([A-Z0-9]+(?:[^\n]+\n[\s]+\/)+[^\n]+)")
find_dimension_elements = re.compile("([A-Z0-9]+)(\([0-9,]+\))")
find_flat_arrs_end = r"\(([0-9])\)[\s]*=([0-9\.DE]+)"
find_flat_arrs = [re.compile("({0})".format(i) + find_flat_arrs_end) for i in flat_arrs]
find_multiarr_elements = re.compile("([A-Z]+)\(([0-9,]+)\)\/(?:\n[\s]+\/)?([0-9DE\.\+\-]+)")

#For testing only
def compareGases(i,j):
    for key in arrays[i]:
        if key in arrays[j]:
            print(key)
            compArr(arrays[i][key],arrays[j][key],8)

#Search text for all Fortran dimension statements and return all the 
#array sizes.
def get_dimensions(text):
    dim_statements = re.findall(find_dimensions,text)
    dimensions = {}
    for d in dim_statements:
        dim_pairs = re.findall(find_dimension_elements,d)
        for t in dim_pairs:
            key = t[0]
            val = t[1][1:-1] #strip off parentheses
            val = val.split(',')
            val = tuple(int(i) for i in val)
            dimensions[key] = val
    return dimensions

#Search TEXT for Fortran data statements and return the array elements in 
#arrays of size DIMENSIONS.
def get_arrays(text, dimensions):
    data = re.findall(find_arrays,text)
    arrays = {}
    for g in data:
        key = g[0]
        val = g[1].replace('E','e')
        val = val.replace('D','e')
        val = val.replace('\n','')
        val = val.replace(' ','')
        val = val.replace('/','')
        def repl(match):
            repetitions = int(match.groups()[0])
            element = [match.groups()[1]]
            repl = repetitions * element
            repl = str(repl)
            repl = repl[1:-1] # strip off parentheses
            return repl
        val = re.sub(find_array_elements, repl,val)
        val = [i.strip().replace("'",'') for i in val.split(',')]
        val = np.array([float(i) for i in val if i != ''])
        val = np.reshape(val,dimensions[key])
        arrays[key] = val
    return arrays

#Search TEXT for one-dimensional array element assignments and add the data
#to ARRAYS, creating new arrays of size DIMENSIONS
def get_flat_arrays(text,arrays,dimensions):
    for assignment in find_flat_arrs:
        flat = re.findall(assignment,text)
        if flat:
            key = flat[0][0]
            dtype = utils.checktype(key)
            arrays[key] = np.zeros(dimensions[key],dtype = dtype)
            for match in flat:
                name = match[0]
                index = int(match[1]) - 1 #Fortran indices start at 1
                dtype = utils.checktype(name)
                val = dtype(match[2])
                arrays[name][index] = val
    return arrays

#Search text for multi-dimensional array element assignments and add the data
#to ARRAYS, creating new arrays out of size DIMENSIONS
def get_multi_arrays(text,arrays,dimensions):
    multi_assignments = re.findall(find_multiarr_elements,text)
    for a in multi_assignments:
        name = a[0]
        dtype = utils.checktype(name)
        if name not in arrays:
            arrays[name] = np.zeros(dimensions[name],dtype=dtype)
        curr = arrays[name]
        index_arr = a[1].split(',')
        for j in index_arr[:-1]:
            curr = curr[int(j) - 1]
        index = int(index_arr[-1]) - 1
        val = a[2].replace('D','e')
        val = val.replace('E','e')
        val = dtype(val)
        curr[index] = val
    return arrays

#Return all arrays associated with ELEMENT in ARRAYS, using POSITION to determine
#association, NUMBER to figure out atomic mass data, and DIMENSIONS to create new 
#arrays if necessary
def get_element_data(element,number,position,arrays,dimensions):
    data = {}
    j = 0
    curr_array_keys = list(arrays.keys())
    while j < len(curr_array_keys):
        key = curr_array_keys[j]
        curr_array = arrays[key]
        if key not in flat_arrs:
            if position == 0 and key in all_array_keys:
                data[key] = curr_array
                del arrays[key]
                del curr_array_keys[j]
            elif (key[:-1] in all_array_keys or key[:-2] in all_array_keys) and key[-1] == element[0]:
                if key[:-1] in all_array_keys:
                    data[key[:-1]] = curr_array
                else:
                    data[key[:-2]] = curr_array
                del arrays[key]
                del curr_array_keys[j]
            else:
                j += 1
        else:
            j += 1
    data['IZ'] = arrays['IZ'][position]
    data['AMZ'] = arrays['AMZ'][position] / number
    ignore = ['XPE','YPE','XCOM','XENE']
    def check_ignore(key):
        for name in ignore:
            if key[:len(name)] == name:
                return False
        return True

    for key in all_array_keys:
        if key not in data and check_ignore(key):
            data[key] = np.zeros(dimensions[key],dtype = utils.checktype(key))
    return data

def write_h5(elements):
    f = h5py.File('gas_data.hdf5','w')
    el = f.create_group('/elements/mixerc')
    for element in elements:
        e = el.create_group(element)
        data = elements[element]
        for array in data:
            if type(data[array]) == np.ndarray:
                e.create_dataset(array,data=data[array])
            else:
                e.attrs.create(array,data = data[array])

if __name__ == "__main__":
    with open(filepath, 'r') as f:
        text = f.read()

    subs = re.findall(find_subs,text)
    arrays = {}
    dimensions = {}
    elements = {}


    for sub in subs:
        i = int(sub[0])
        dimensions[i] = get_dimensions(sub[1])
        arrays[i] = get_arrays(sub[1],dimensions[i])
        arrays[i] = get_flat_arrays(sub[1],arrays[i],dimensions[i])
        arrays[i] = get_multi_arrays(sub[1],arrays[i],dimensions[i])

    for gas in gas_dict:
        formula = gas_dict[gas]['formula']
        i = 0
        for pair in formula:
            element = utils.getSingle(pair)
            number = pair[element]
            if element not in elements:
                elements[element] = get_element_data(element,number,i,arrays[gas],dimensions[gas])
            i += 1

    write_h5(elements)
