import unittest
import pygrad

class TestSetup(unittest.TestCase):
    setupErrorPath = 'tests/errors/setup'
    def ErrorFileSetup(self,fname,error):
        with self.assertRaises(error):
            pygrad.Main(fname)

    def ErrorFileSetupFormat(self,fname):
        self.ErrorFileSetup(fname, pygrad.InfileFormatException)

    def testMissingLine(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-01.txt')

    def testMissingFile(self):
        self.ErrorFileSetup(setupErrorPath+'error-00.txt', FileNotFoundError)
    
    def testBlankFile(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-02.txt')

    def testNoGas(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-03.txt')

    def testAbsZero(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-04.txt')

    def testZeroPressure(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-05.txt')

    def testTooFewParameters(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-06.txt')

    def testTooManyParameters(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-07.txt')
    
    def testNgasBig(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-08.txt')

    def testNgasMismatch(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-09.txt')

    def testNfracMismatch(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-10.txt')
    
    def testNgasNMismatch(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-11.txt')

    def testDuplicateNgas(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-12.txt')
        
    def testIncompletePercentage(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-13.txt')

    def testMaxEnergy(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-14.txt')

    def testMaxXrayEvents(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-15.txt')
    
    def testMaxMIPSEvents(self):
        self.ErrorFileSetupFormat(setupErrorPath+'error-16.txt')

    def testFloatNgas(self):
        self.ErrorFileSetup(setupErrorPath+'error-17.txt',ValueError)

    def testCorrectSetup(self):
        main = pygrad.Main('tests/success-ins/success-01.txt')
        

if __name__ == '__main__':
    unittest.main()

