import argparse, numpy

parser = argparse.ArgumentParser(description = 'Run Degrad')
parser.add_args('infile',nargs=1,help='File containing initial parameters')
parser.add_args('outfile',nargs='+',help='Filename for the output of the file',default='DEGRAD.OUT')



