import numpy as np

class PygradException(Exception):
    pass

class InfileFormatException(PygradException):
    pass

class Main():
    def __init__(self, infile):
        with open(infile, 'r') as f:
            text = f.read()
        lines = text.split('\n')[:-1]
        if len(lines) < 5:
            raise InfileFormatException('Not enough lines in input file')
        parameters = [line.split(',') for line in lines[:5]]

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Run Degrad')
    parser.add_args('infile',nargs=1,help='File containing initial parameters')
    parser.add_args('outfile',nargs='+',help='Filename for the output of the file',default='DEGRAD.OUT')
    args = parser.parse_args()
    main = Main(args.infile)

