import unittest, os
from Pygrad.Pygrad import Pygrad
from Pygrad import utils

#A class to test the setup and density functions.
class TestSetup(unittest.TestCase):
    pygrad_home = os.getenv('PYGRAD_HOME')
    setupErrorPath = pygrad_home + '/tests/errors/setup/'
    densityTestPath = pygrad_home + '/tests/fortran_tests/density-tests/'
    def ErrorFileSetup(self,fname,error):
        with self.assertRaises(error):
            pygrad.Main(fname)

    def testDensityCalc(self):
        
        def densityCalc(fname):
            main = pygrad.Main(self.densityTestPath + fname)
            main.calcDensity()
            return [main.den]
        utils.checkFortranTests(self.densityTestPath,densityCalc,"01",2)

    def ErrorFileSetupFormat(self,fname):
        self.ErrorFileSetup(fname, pygrad.InfileFormatException)

    def testMissingLine(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-01.txt')

    def testMissingFile(self):
        self.ErrorFileSetup(self.setupErrorPath+'error-00.txt', FileNotFoundError)
    
    def testBlankFile(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-02.txt')

    def testNoGas(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-03.txt')

    def testAbsZero(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-04.txt')

    def testZeroPressure(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-05.txt')

    def testTooFewParameters(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-06.txt')

    def testTooManyParameters(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-07.txt')
    
    def testNgasBig(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-08.txt')

    def testNgasMismatch(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-09.txt')

    def testNfracMismatch(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-10.txt')
    
    def testNgasNMismatch(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-11.txt')

    def testDuplicateNgas(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-12.txt')
        
    def testIncompletePercentage(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-13.txt')

    def testMaxEnergy(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-14.txt')

    def testMaxXrayEvents(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-15.txt')
    
    def testMaxMIPSEvents(self):
        self.ErrorFileSetupFormat(self.setupErrorPath+'error-16.txt')

    def testFloatNgas(self):
        self.ErrorFileSetup(self.setupErrorPath+'error-17.txt',ValueError)

    def testCorrectSetup(self):
        main = pygrad.Main(self.pygrad_home + '/tests/correct/success-01.txt')
        

if __name__ == '__main__':
    unittest.main()

