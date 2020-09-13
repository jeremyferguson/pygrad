import os,re,h5py, numpy as np,copy
import utils,pygrad,cascdata,argparse


verbosity = 0
pygrad_home = os.getenv('PYGRAD_HOME')
filepath = pygrad_home + '/degrad-3.9a.f'
#Names of all arrays found in MIXERC or its subroutines.
all_array_keys = ['A','R','INIOC', 'ES', 'PRBSH', 'PRBSHBT', 'YPEK', 'XPEK', 'YPEL1', 'XPEL1', 'YPEL2', 'XPEL2', 'YPEL3', 'XPEL3', 'YPEM1', 'XPEM1', 'YPEM2', 'XPEM2', 'YPEM3', 'XPEM3', 'YPEM4', 'XPEM4', 'YPEM5', 'XPEM5', 'YPEN1', 'XPEN1', 'YPEN2', 'XPEN2', 'YPEN3', 'XPEN3', 'YPEN4', 'XPEN4', 'YPEN5', 'XPEN5', 'YPEO1', 'XPEO1', 'YPEO2', 'XPEO2', 'YPEO3', 'XPEO3', 'XENE', 'YRAY', 'YCOM', 'YPAP', 'FFR', 'FFC', 'IZ', 'AMZ','XCOM']
gas_dict = copy.deepcopy(pygrad.gas_dict)
local = ['J','K','F','ELOSS','GAMMA1','APOP1','APOP2','APOP3','APOP4','EN']
rot_arrs = ['RSUM','PJ']
fortran_functions = {'ACOS':'np.arccos','DEXP':'np.exp','DFLOAT':'float'}
mixer_arrs_first = ['E','EION']
mixer_arrs_second = ['BEF','KIN','EOBY','IZBR','NC0','EC0','WKLM','EFL','NG1','EG1','NG2','EG2','LEGAS','ISHELL','SCLN']
mixer_arrs_complex = ['IZBR']
flat_arrs = ['IZ','AMZ']
rot_vars = ['B0','QBQA','QBK','GPARA','GORTHO','DBA','DRAT','DBK','AMPV2','AMPV3','APOL','LMAX','AA','DD','FF','A1','B1','A2','EOBFRAC','B0','EATTTH','EATTWD','AMPATT','EATTTH1','EATTWD1','AMPATT1','ESCOBY','SCLOBY','RAT','AR','BR','AEXT20','AGST20']
arrays = {}
dimensions = {}
elements = {}
variables = {}

def fortran_eval(expr,namespace,arrays,lstart=0,lend=0,lincr=1):
    def repl(match):
        def variablesub(expr,space):
            if '(' in expr:
                start = expr.find('(')
                end = expr.find(')')
                arr = expr[:start]
                index = expr[start+1:end]
                if index.isdigit():
                    index = int(index) - 1
                    return str(space[arr][index])
                elif index == 'J' and 'J' not in space:
                    return repr(space[arr][lstart:lend:lincr])
                else:
                    index = fortran_eval(index) - 1
                    return str(space[arr][index])
            else:
                return str(space[expr])

        expr = match.group()
        if expr in fortran_functions:
            return fortran_functions[expr]
        elif expr.isdigit():
            return expr
        elif '.' in expr:
            return expr.replace('D','e').replace('E','e')
        if '(' in expr:
            arr = expr[:expr.find('(')]
        else:
            arr = ''
        if expr in arrays or arr in arrays:
            return variablesub(expr,arrays)
        elif expr in namespace or arr in namespace:
            return variablesub(expr,namespace)
        elif expr in pygrad.glob() or arr in pygrad.glob():    
            return variablesub(expr,pygrad.glob())
        elif expr in cascdata.glob() or arr in cascdata.glob():
            return variablesub(expr,cascdata.glob())
        else:
            return str(utils.checktype(expr)(0))
    pattern = re.compile('[\w]+(?:(?:\.(?:[\d]*[DE][\-\+]?)?[\d]+)|(?:\([\dJK]+\)))?')
    expr = re.sub(pattern,repl,expr)
    return eval(expr,{'array':np.array})

#Search TEXT for all Fortran dimension statements and return all the 
#array sizes.
def get_dimensions(text):
    find_dimensions = re.compile(r"DIMENSION ([A-Z0-9]+(?:[^\n]+\n[\s]+\/)*[^\n]+)")
    find_dimension_elements = re.compile(r"([A-Z0-9]+)(\([0-9,]+\))")
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
def get_arrays(text, dimensions,arrays):
    find_arrays = re.compile(r"DATA ([A-Z0-9]+)((?:(?:[\/0-9][0-9\.,/*ED\-]+\/?)[\s]+)+)")
    find_array_elements = re.compile(r"([0-9]+)\*([0-9\.e\-]+)")
    data = re.findall(find_arrays,text)
    for g in data:
        key = g[0]
        val = g[1].replace('E+','e+')
        val = val.replace('E-','e-')
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
#to ARRAYS, using DIMENSIONS to determine the sizes of any new arrays. 
def get_flat_arrays(text,arrays,dimensions):
    find_flat_arrs_end = r"\(([0-9])\)[\s]*=([0-9\.DE]+)"
    find_flat_arrs = [re.compile("({0})".format(i) + find_flat_arrs_end) for i in flat_arrs]
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

#Search TEXT for multi-dimensional array element assignments and add the data
#to ARRAYS, using DIMENSIONS to determine the sizes of any new arrays. 
def get_multi_arrays(text,arrays,dimensions):
    find_multiarr_elements = re.compile(r"([A-Z]+)\(([\d],[\d,]+)\)\/(?:\n[\s]+\/)?([0-9DE\.\+\-]+)")
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
#association, NUMBER to calculate atomic mass data, and DIMENSIONS to create new 
#arrays if necessary.
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

def get_vars(sub,arrays,dimensions):
    data = get_block_vars(sub,arrays)
    for var in rot_vars:
        data = get_ind_var(sub,var,data,arrays,dimensions)
    data['fullname'] = get_names(sub)
    data['description'] = get_description(sub,data['fullname'])
    data['nulldescription'] = get_null_description(sub)
    data = get_fragment(sub,'DEG',data,arrays,dimensions)
    return data

#Find all the variable assignments in the block near the beginning of SUB
def get_block_vars(sub,arrays):
    find_vars = re.compile(r"(?:[\w]+ *= *[\w\*/\.DE\+\-\(\)]+[\s]+)(?:[\w]+ *= *[\w\*/\.DE\+\-\(\)]+[\s]+)+")
    find_assignments = re.compile(r"([\w]+) *= *([\w\.\+\-\*/DE\(\)]+)")
    blocks = re.findall(find_vars,sub)
    data = {}
    fullblocks = 0
    i = 0
    while 'NIN' not in data:
        block = blocks[i]
        complete = True
        statements = re.findall(find_assignments,block)
        for match in statements:
            if match[0] in local:
                break
            if match[0] not in pygrad.glob():
                data[match[0]] = fortran_eval(match[1],data,arrays)
        i += 1
    return data

#Find both long form names of the gas in SUB.
def get_names(sub):
    def process_name(groups):
        if groups[0] and int(groups[0]) == pygrad.NANISO:
            return groups[1]
        elif groups[2] and int(groups[2]) == pygrad.NANISO:
            return groups[3]
        elif groups[3]:
            return groups[3]
        else:
            return groups[1]

    find_gas_name = re.compile(r"(?:IF\(NANISO\.EQ\.([\d])\) *(?:THEN)?[\s]*\n?[\s]*)?NAME=\'([\s\S]*?)\'(?:[\s]*(?:(?:ELSE)|(?:IF\(NANISO\.EQ\.([\d])\)))[\s]*NAME=\'([\s\S]*?)\')?")
    find_naniso = re.compile(r"NANISO[ ]*=[ ]*([\d])")
    names = re.search(find_gas_name,sub)
    naniso_match = re.search(find_naniso,sub)
    if naniso_match:
        if names.start() > naniso_match.end():
            pygrad.NANISO = int(naniso_match.group(1))
            name = process_name(names.groups(''))
        else:
            name = process_name(names.groups(''))
            pygrad.NANISO = int(naniso_match.group(1))
    else:
        name = process_name(names.groups(''))
    if verbosity > 1:
        print(name)
    return name

#Returns the description found in SUB, substituting NAME in for the first two lines.
def get_description(sub,name):
    find_scrpt = re.compile(r"SCRPT\(([\d])+\)[ ]*[ ]*=[ ]*\'([\s\S]*?)\'")
    matches = re.findall(find_scrpt,sub)
    description = "\n{}\n".format(name)
    for match in matches:
        if int(match[0]) > 2:
            description += match[1] + '\n'
    return description

#Returns the null description found in SUB.
def get_null_description(sub):
    find_scrptn = re.compile(r"SCRPTN\((?:[\d])+\)[ ]*[ ]*=[ ]*\'([\s\S]*?)\'")
    matches = re.findall(find_scrptn,sub)
    description = ""
    if matches:
        for match in matches:
            description += match + '\n'
    return description

#Returns ARRAYS with added data from the do loop assignments in SUB, using DIMENSIONS
#to create new arrays as necessary, and VAR to find variables.
def get_do_assignments(sub,dimensions,arrays,var):
    find_do_loop = re.compile(r"DO *[\d]+ *[JKLN]+ *= *([\w]+),([\w]+),?([\w]+)?[\s]+((?:(?:(?:[\w]+\((?:[\d],)?[JKLN]+\) *= *[\w\.\+\-\*/\(\)]+)|(?:PENSUM=PENSUM\+PENFRAC\(1,K\)))[\s]+)*(?:[\d]+[\s]+)(?:(?:[\w]+\((?:[\d],)?[JKLN]+\) *= *[\w\.\+\-\*/\(\)]+[\s]+)|(?:CONTINUE)))")
    find_do_line = re.compile(r"([\w]+)\(([\d])?,?[JKLN]+\) *= *([\w/\*\+\-\(\)\.]+)")
    matches = re.findall(find_do_loop,sub)
    for match in matches:
        start = fortran_eval(match[0],var,arrays) - 1
        end = fortran_eval(match[1],var,arrays)
        incr = fortran_eval(match[2],var,arrays) if match[2] else 1
        lines = re.findall(find_do_line,match[3])
        for line in lines:
            if verbosity > 1:
                print(line)
            if line[0] not in rot_arrs:
                val =  fortran_eval(line[2],var,arrays,start,end,incr)
                if line[1]:
                    arrays[line[0]][int(line[1])-1,start:end:incr] = val
                else:
                    arrays[line[0]][start:end:incr] = val
    return arrays

#Get individual array assignments of the name NAME in SUB, using DIMENSIONS, ARRAYS, and VAR 
#for evaluation of expressions.
def get_ind_arrs(name,sub,dimensions,arrays,var):
    pattern = re.compile("{}\(([\d]+)\) *= *([\w/\*\+\-\(\)\.]+)".format(name))
    if name in dimensions:
        if name not in arrays:
            arrays[name] = np.zeros(dimensions[name],dtype=utils.checktype(name))
        matches = re.findall(pattern,sub)
        for match in matches:
            arrays[name][int(match[0])-1] = fortran_eval(match[1],var,arrays)
    return arrays

#Get complex individual array assignments of NAME in SUB, using DIMENSIONS, ARRAYS, and VAR for evaluation of expressions.
def get_complex_ind_arrs(name,sub,dimensions,arrays,var):
    pattern = re.compile("{}\(([\w]+\+[\w]+)\) *= *([\w/\*\+\-\(\)\.]+)".format(name))
    if name not in arrays:
        arrays[name] = np.zeros(dimensions[name],dtype = utils.checktype(name))
    matches = re.findall(pattern,sub)
    for match in matches:
        index = fortran_eval(match[0],var,arrays) - 1
        val = fortran_eval(match[1],var,arrays)
        arrays[name][index] = val
    return arrays

#Get the value of NAME in SUB and assign it to VAR, using ARRAYS and DIMENSIONS for
#evaluation of expressions
def get_ind_var(sub,name,var,arrays,dimensions):
    pattern = re.compile("{} *= *([\w/\*\+\-\(\)\.]+)".format(name))
    matches = re.findall(pattern,sub)
    if matches:
        for match in matches:
            var[name] = fortran_eval(match,var,arrays)
    return var

#Get all variable assignments in SUB that start with the fragment FRAG, using VAR,
#ARRAYS, and DIMENSIONS for evaluation of expressions
def get_fragment(sub,frag,var,arrays,dimensions):
    find_fragment = r"({}[\w]+) *= *([\w/\*\+\-\(\)\.]+)"
    pattern = re.compile(find_fragment.format(frag))
    matches = re.findall(pattern,sub)
    for match in matches:
        var[match[0]] = fortran_eval(match[1],var,arrays)
    return var

#Get all the multi-dimensional individual array assignments in SUB, using VAR and ARRAYS
#to evaluate expressions and DIMENSIONS to make new arrays.
def get_multi_ind(sub,var,arrays,dimensions):
    find_multi_ind = re.compile(r"([\w]+)\(([\d],[\d,]+)\) *= *([\w/\*\+\-\(\)\.]+)")
    matches = re.findall(find_multi_ind,sub)
    for match in matches:
        name = match[0]
        if name not in arrays:
            arrays[name] = np.zeros(dimensions[name],dtype=utils.checktype(name))
        curr = arrays[name]
        index_arr = match[1].split(',')
        for j in index_arr[:-1]:
            curr = curr[int(j) - 1]
        index = int(index_arr[-1]) - 1
        val = fortran_eval(match[2],var,arrays)
        curr[index] = val
    return arrays

#Get all the different penfra assignments in SUB, using VAR and ARRAYS for expression
#evaluation and DIMENSIONS to make new arrays.
def get_penfra_special(sub,var,arrays,dimensions):
    #Pattern to find a special form of PENFRA
    find_penfra_if = re.compile(r"DO [\d]+ [JKLN]+=([\w]+),([\w]+)[\s]+((?:PENFRA\([\d],[JLKN]+\)=[\w\.\+\-\*/\(\)]+[\s]+)+)IF\(IPEN.EQ.0\) GO TO [\d]+[\s]+WRITE\([\d]+,[\d]+\) [,\w\.\+\-\*/\(\)]+[\s]+[\d]+ *FORMAT\(['=,\s\w\.\+\-\*/\(\)]+\)[\s]+[\d]+ CONTINUE")
    if 'PENFRA' not in arrays:
        arrays['PENFRA'] = np.zeros(dimensions['PENFRA'],dtype = float)
    if_matches = re.findall(find_penfra_if,sub)
    penfra_assign = re.compile('PENFRA\(([\d]),([JKLN])\)=([\d\.]+)')
    for match in if_matches:
        start = fortran_eval(match[0],var,arrays) - 1
        end = fortran_eval(match[1],var,arrays)
        assign_matches = re.findall(penfra_assign,match[2])
        for m in assign_matches:
            index = fortran_eval(m[0],var,arrays) - 1
            val = fortran_eval(m[2],var,arrays)
            arrays['PENFRA'][index][start:end] = val
    find_penfra_after = re.compile(r"DO [\d]+ ([JKLN]+)=([\w\+\-\.\*/\(\)]+),([\w\+\-\.\*/\(\)]+)[\s]+DO [\d]+ ([JKLN]+)=([\w\+\-\.\*/\(\)]+),([\w\+\-\.\*/\(\)]+)[\s]+ [\d]+ *PENFRA\([\w\+\-\.\*/\(\)],[\w\+\-\.\*/\(\)]\)=[\w\+\-\.\*/\(\)]+[\s]+((?:PENFRA\([\w\+\-\.\*/\(\)],[\w\+\-\.\*/\(\)]\)=[\w\+\-\.\*/\(\)]+[\s]+)+)")
    after_matches = re.findall(find_penfra_after,sub)
    for match in after_matches:
        local_space = copy.deepcopy(var)
        local_space[match[0]] = fortran_eval(match[2],var,arrays)
        local_space[match[3]] = fortran_eval(match[5],var,arrays)
        assign_matches = re.findall(penfra_assign,match[6])
        for m in assign_matches:
            i = fortran_eval(m[0],local_space,arrays) - 1
            j = fortran_eval(m[1],local_space,arrays) - 1
            val = fortran_eval(m[2],local_space,arrays)
            arrays['PENFRA'][i][j] = val
    return arrays

def get_lambda_frag(frag,sub,var,arrays):
    pattern = re.compile(r" ({}[A-HJ-Z0-9]+)=([\(\)\w\+\-\.\*/]+)".format(frag))
    parameters = ['SUM','AKT','TEMPC','TORR']
    matches = re.findall(pattern,sub)
    result = []
    for match in matches:
        left = match[0]
        if left not in var:
            right = tokenize(match[1])
            expr = symbol_eval(right,var,arrays)
            params = [p for p in parameters if p in expr]
            result.append([left,'lambda {}: {}'.format(",".join(params),expr)])
    arrays[frag] = result
    return arrays

#Turn the rest of EXPR into a list and add to TOKENS, using CURR to track the current token.
def tokenize(expr):
    def t(curr,expr,tokens):
        if not expr:
            return tokens + [curr]
        if expr[0] not in '*/+-()':
            return t(curr+expr[0],expr[1:],tokens)
        elif expr[0] in '*/+-' or curr not in fortran_functions:
            return t('',expr[1:],tokens + [curr,expr[0]])
        else:
            return t('',expr[1:],tokens + [fortran_functions[curr],expr[0]])
    return t('',expr,[])

#Evaluate EXPR symbolically using VAR and ARRAYS.
def symbol_eval(expr,var,arrays):
    k = 0
    while k < len(expr):
        token = expr[k]
        expr[k] = token.replace('D-','e-')
        if token in var:
            expr[k] = str(var[token])
        elif token in pygrad.glob():
            expr[k] = str(pygrad.glob()[token])
        elif token in arrays:
            l = expr[k:].index('(')
            r = expr[k:].index(')')
            expr[k+l] = '['
            expr[k+r] = ']'
        k += 1
    string = "".join(expr)
    pattern = re.compile(r"\[([\w])\]")
    string = re.sub(pattern,r"[\1-1]",string)
    return string 

#Carry out the variable assignment in LINE, using VAR, ARRAYS, START,END,INCR, and ITERATOR 
#for evaluation of expressions.
def symbol_assign(line,var,arrays,start,end,incr,iterator):
    left = line[:line.find('=')]
    right = tokenize(line[line.find('=')+1:])
    expr = symbol_eval(right,var,arrays)
    if '(' not in left:
        var[left] = expr
    else:
        wrapped = "lambda AKT: [{} for {} in range({},{},{})]".format(expr,iterator,start,end,incr)
        arrays['PJ'].append([[start-1,end,incr],wrapped])
    return var,arrays

#Evaluate all the tokens in EXPR using VAR and ARRAYS.
def tokens_eval(expr,var,arrays):
    k = 0
    while k < len(expr):
        token = expr[k]
        expr[k] = token.replace('.D','.e')
        if token in var:
            expr[k] = str(var[token])
        elif token in pygrad.glob():
            expr[k] = str(pygrad.glob()[token])
        elif token in arrays:
            l = expr[k:].index('(')
            r = expr[k:].index(')')
            index = fortran_eval(''.join(expr[k+l+1:k+r]),var,arrays) - 1
            expr[k:k+r+1] = [str(arrays[token][index])]
        k += 1
    expr = eval("".join(expr))
    return expr

#Carry out a normal variable assignment in LINE, using VAR and ARRAYS.
def normal_assign(line,var,arrays):
    left = line[:line.find('=')]
    right = tokenize(line[line.find('=')+1:])
    expr = tokens_eval(right,var,arrays)
    if '(' not in left:
        var[left] = expr
    else:
        l = left.index('(')
        r = left.index(')')
        name = left[:l]
        index = fortran_eval(left[l+1:r],var,arrays) - 1
        arrays[name][index] = expr
    return var,arrays

#Turn an if statement found in MATCH into a python conditional statement using 
#VAR and ARRAYS.
def eval_if(match,var,arrays):
    cond = match[0] + match[1]
    cond = cond.replace('\n','')
    cond = cond.replace(' ','')
    cond = cond.replace('/','')
    cond = cond.replace('.EQ.',' == ')
    cond = cond.replace('.OR.',' or ')
    cond = cond.replace('(','[')
    cond = cond.replace(')',']')
    cond = cond[:-1]
    thenvar = match[2]
    thenval = match[3]
    elseval = var[thenvar]
    var[thenvar] = "({} if ({}) else {})".format(thenval,cond,elseval)
    return var

def eval_do_loop(match,var,arrays,dimensions):
    do_pattern = re.compile(r"DO ([\d]+) ([JKLN])=([\w]+),([\w\+]+),?([\w]+)?([\s\S]*?)\n *\1 *([\S]+)")
    if_pattern = re.compile(r"IF\(([\w\(\)\.]+)[\s]+((?:\/[\w\.\(\)]+[\s]+)*)([\w]+)=([\d\.]+)")
    iterator = match[1]
    start = fortran_eval(match[2],var,arrays)
    end = fortran_eval(match[3],var,arrays) + 1
    incr = fortran_eval(match[4],var,arrays) if match[4] else 1
    for i in range(start,end,incr):
        var[iterator] = i
        lines = match[5].strip().split('\n')
        lines = [j.strip() for j in lines]
        lines += [match[6].strip()]
        j = 0
        while j < len(lines):
            line = lines[j]
            if line[:2] == 'DO':
                loop = re.search(do_pattern, "{}\n {} {}".format(match[5],match[0],match[6]))
                j += len(loop.group().split('\n'))
                arrays = eval_do_loop(loop.groups(),var,arrays,dimensions)
            else:
                var, arrays = normal_assign(line,var,arrays)
                j += 1
    return arrays

def eval_pj_loop(match,var,arrays,dimensions):
    if_pattern = re.compile(r"IF\(([\w\(\)\.]+)[\s]+((?:\/[\w\.\(\)]+[\s]+)*)([\w]+)=([\d\.]+)")
    iterator = match[1]
    start = fortran_eval(match[2],var,arrays)
    end = fortran_eval(match[3],var,arrays) + 1
    incr = fortran_eval(match[4],var,arrays) if match[4] else 1
    lines = match[5].strip().split('\n')
    lines = [j.strip() for j in lines]
    lines += [match[6].strip()]
    j = 0
    while j < len(lines):
        line = lines[j]
        if line[:2] == 'IF':
            cond = re.search(if_pattern,"{}\n {} {}".format(match[5],match[0],match[6]))
            j += len(cond.group().split('\n'))
            var = eval_if(cond.groups(),var,arrays)
        else:
            var, arrays = symbol_assign(line,var,arrays,start,end,incr,iterator)
            j += 1
    return arrays

def get_sumstart(sub):
    pattern = re.compile(r' (?:A|R)?SUM=([\d\.]+)')
    match = re.search(pattern,sub)
    if match:
        return float(match[1])

def do_rot_calcs(sub,var,arrays,dimensions):
    do_pattern = re.compile(r"DO ([\d]+) ([JKLN])=([\w]+),([\w\+]+),?([\w]+)?([\s\S]*?)\n *\1 *([\S]+)")
    matches = re.findall(do_pattern,sub)
    elev_arrs = []
    ein_arrs = []
    pj_arrs = []
    for match in matches:
        elevp = re.compile('ELEV\([A-Z]\)=')
        if re.findall(elevp,match[5] + match[6]):
            elev_arrs.append(match)
        einp = re.compile('EIN\([\w\*\-\+]+\)=')
        if re.findall(einp,match[5] + match[6]):
            ein_arrs.append(match)
        pjp = re.compile('PJ\([A-Z]\)=')
        if re.findall(pjp,match[5] + match[6]):
            if 'SUM' not in match[5] + match[6] and 'RSUM' not in match[5] + match[6] and 'ASUM' not in match[5] + match[6]:
                pj_arrs.append(match)
    for arr in elev_arrs:
        local = copy.deepcopy(var)
        local['L'] = 1
        arrays = eval_do_loop(arr,local,arrays,dimensions)
    arrays = get_ind_arrs('EIN',sub,dimensions,arrays,var)
    for arr in ein_arrs:
        local = copy.deepcopy(var)
        arrays = eval_do_loop(arr,local,arrays,dimensions)
    pj_pattern = re.compile('PJ\(([\d]+)\)=([\d\.]+)')
    pj_assign = re.findall(pj_pattern,sub)
    arrays['PJ'] = []
    for match in pj_assign:
        arrays['PJ'].append([[int(match[0])-1,int(match[0]),1],"lambda AKT: [{}]".format(match[1])])
    for arr in pj_arrs:
        local = copy.deepcopy(var)
        arrays = eval_pj_loop(arr,local,arrays,dimensions)
    var['SUMSTART'] = get_sumstart(sub)
    arrays = get_lambda_frag('FROT',sub,var,arrays)
    arrays = get_ind_arrs('ERLVL',sub,dimensions,arrays,var)
    arrays = get_lambda_frag('AP',sub,var,arrays)
    arrays = get_lambda_frag('FA',sub,var,arrays)
    return arrays

def get_loop_params(sub,var,arrays):
    ABpattern = r"A=\(([\w]+)\(J\)\-\1\(J\-1\)\)/\(([\w]+)\(J\)\-\2\(J\-1\)\)[\s]+B=\(\2\(J\-1\)\*\1\(J\)\-\2\(J\)\*\1\(J\-1\)\)/\(\2\(J\-1\)-\2\(J\)\)[\s]+{}(?:\(([\d]+),I\))?=\(A\*EN\+B\)(\*1\.D\-16)?"
    ENpattern = re.compile(r"(IF\(EN\.GT\.(DABS\()?EIN\(1\)\)?\) THEN[\s]+)?GAMMA1=\(EMASS2")
    gammacond = re.search(ENpattern,sub)
    gammacond = [bool(gammacond.group(1)),bool(gammacond.group(2))]
    var['GAMMACOND'] = gammacond

#Write all the data in ELEMENTS to a new hdf5 file.
def write_h5(elements):
    f = h5py.File('gas_data.hdf5','w')
    el = f.create_group('/elements')
    for element in elements:
        e = el.create_group(element+'/mixerc/')
        data = elements[element]
        for array in data:
            if type(data[array]) == np.ndarray:
                e.create_dataset(array,data=data[array])
            else:
                e.attrs.create(array,data = data[array])

def main():
    global arrays,dimensions,subs,variables
    with open(filepath, 'r') as f:
        text = f.read()
    #Comment pattern
    comments = re.compile(r"\n(?:[cC][\s\S]*?\n)+")
    #Pattern to find cgas subroutines in the Fortran file.
    find_csubs = re.compile(r"SUBROUTINE CGAS([0-9]+)([\S\s]*?)END[\s]*\n")
    #Pattern to find gas subroutines in the Fortran file.
    find_subs = re.compile(r"SUBROUTINE GAS([0-9]+)([\S\s]*?)      END[^I]")
    text = re.sub(comments,"\n",text)
    subs = re.findall(find_csubs,text)
    for sub in subs:
        i = int(sub[0])
        if i in gas_dict:
            dimensions[i] = get_dimensions(sub[1])
            arrays[i] = {}
            arrays[i] = get_arrays(sub[1],dimensions[i],arrays[i])
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
    subs = re.findall(find_subs,text)
    arrays = {}
    dimensions = {}
    loops = {}
    for sub in subs:
        oldnaniso = pygrad.NANISO
        i = int(sub[0])
        if i in gas_dict:
            big_loop = re.compile(r"DO [\d]+ [IJ]=1,NSTEP")
            split = re.search(big_loop,sub[1])
            top = sub[1][:split.start()]
            bot = sub[1][split.start():]
            variables[i] = {}
            if verbosity > 0:
                print(i)
            dimensions[i] = get_dimensions(sub[1])
            arrays[i] = {}
            for arr in dimensions[i]:
                arrays[i][arr] = np.zeros(dimensions[i][arr],dtype=utils.checktype(arr))
            arrays[i] = get_arrays(top,dimensions[i],arrays[i])
            for arr in mixer_arrs_first:
                arrays[i] = get_ind_arrs(arr,top,dimensions[i],arrays[i],variables[i])
            variables[i] = get_vars(top,arrays[i],dimensions[i])
            arrays[i] = get_do_assignments(top,dimensions[i],arrays[i],variables[i])
            for arr in mixer_arrs_second:
                arrays[i] = get_ind_arrs(arr,top,dimensions[i],arrays[i],variables[i])
            for arr in mixer_arrs_complex:
                arrays[i] = get_complex_ind_arrs(arr,top,dimensions[i],arrays[i],variables[i])
            arrays[i] = get_multi_ind(top,variables[i],arrays[i],dimensions[i])
            arrays[i] = get_penfra_special(top,variables[i],arrays[i],dimensions[i])
            arrays[i] = do_rot_calcs(top,variables[i],arrays[i],dimensions[i])
            loops[i] = get_loop_params(bot,variables[i],arrays[i])
        pygrad.NANISO = oldnaniso
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v','--verbosity',type=int,default=1)
    verbosity = parser.parse_args().verbosity
    main()
