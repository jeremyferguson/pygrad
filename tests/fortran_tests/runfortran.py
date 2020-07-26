import argparse, os

parser = argparse.ArgumentParser(description='Run a fortran program with several inputs and record the outputs')
parser.add_argument('program',type=str,help='fortran program to run')
parser.add_argument('directory',type=str,help='directory where input files are located')
parser.add_argument('outfile',type=str,help='File that the fortran program outputs to')
args = parser.parse_args()

for fname in os.listdir(args.directory):
    if len(fname) > 3 and fname[-3:] == '.in':
        prefix = fname[:-3]
        os.system('./'+args.program + ' < ' + args.directory + fname)
        with open(args.outfile,'r') as f:
            text = f.read()
        #Converts all the Fortran output into python ints or floats, then converts them back to strings
        data = [[str(float(i.strip()) if '.' in i else int(i.strip())) for i in line.split(' ') if i != ''] for line in text.replace('D+','e').replace('D-','e-').split('\n')[:-1]]
        with open(args.directory + prefix + '.out','w') as f:
            f.write("\n".join([','.join(line) for line in data]))


