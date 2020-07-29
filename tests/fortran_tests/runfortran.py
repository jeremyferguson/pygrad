import argparse, os,re

parser = argparse.ArgumentParser(description='Run a fortran program with several inputs and record the outputs')
parser.add_argument('program',type=str,help='fortran program to run')
parser.add_argument('directory',type=str,help='directory where input files are located')
parser.add_argument('outfiles',type=int,help='Number of the output files from the Fortran program')
args = parser.parse_args()

outfiles = [(args.program +str(i) +  '.out').upper() for i in range(1,args.outfiles+1)]
#Runs all input files in the given directory and store their output
for fname in os.listdir(args.directory):
    if len(fname) > 3 and fname[-3:] == '.in':
        prefix = fname[:-3]
        os.system('./'+args.program + ' < ' + args.directory + fname)
        i = 1
        for out in outfiles:
            with open(out,'r') as f:
                text = f.read()
        #Converts all the Fortran output into python ints or floats, then converts them back to strings
            text = text.replace('D+','e')
            text = text.replace('D-','e-')
            def repl(match):
                number = match.groups()[0]
                sign = match.groups()[1]
                exp = match.groups()[2]
                return "{}e{}{}".format(number,sign,exp)
            pattern = re.compile("([0-9\.]+)([-\+])([0-9]+)")
            text = re.sub(pattern,repl,text)
            data = [line for line in text.split('\n')]
            data = data[:-1]
            def replace(i):
                i = i.strip()
                if '.' in i:
                    i = float(i)
                elif i == 'Infinity':
                    i = float('inf')
                elif i == '-Infinity':
                    i = -float('inf')
                else:
                    i = int(i)
                return str(i)
            data = [[replace(i) for i in line.split(' ') if i != ''] for line in data]
            outpath = '{}{}'.format(args.directory,prefix)
            if not os.path.exists(outpath):
                os.system('mkdir {}{}'.format(args.directory,prefix))
            with open("{}/{:02}.out".format(outpath,i),'w') as f:
                f.write("\n".join([','.join(line) for line in data]))
            i += 1


