import numpy

class PygradException(Exception):
    pass

class InfileFormatException(PygradException):
    pass

def setup(infile):
    pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Run Degrad')
    parser.add_args('infile',nargs=1,help='File containing initial parameters')
    parser.add_args('outfile',nargs='+',help='Filename for the output of the file',default='DEGRAD.OUT')
    args = parser.parse_args()
    setup(args.infile)

