import unittest
import pygrad

class TestSetup(unittest.TestCase):

    def ErrorFileSetup(self,fname,error):
        with self.assertRaises(error):
            pygrad.Main(fname)

    def ErrorFileSetupFormat(self,fname):
        self.ErrorFileSetup(fname, pygrad.InfileFormatException)

    def testMissingLine(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-01.txt')

    def testMissingFile(self):
        self.ErrorFileSetup('tests/error-ins/error-00.txt', FileNotFoundError)
    
    def testBlankFile(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-02.txt')

    def testNoGas(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-03.txt')

    def testAbsZero(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-04.txt')

    def testZeroPressure(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-05.txt')

    def testTooFewParameters(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-06.txt')

    def testTooManyParameters(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-07.txt')
    
    def testNgasBig(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-08.txt')

    def testNgasMismatch(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-09.txt')

    def testNfracMismatch(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-10.txt')
    
    def testNgasNMismatch(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-11.txt')

    def testDuplicateNgas(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-12.txt')
        
    def testIncompletePercentage(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-13.txt')

    def testMaxEnergy(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-14.txt')

    def testMaxXrayEvents(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-15.txt')
    
    def testMaxMIPSEvents(self):
        self.ErrorFileSetupFormat('tests/error-ins/error-16.txt')

    def testFloatNgas(self):
        self.ErrorFileSetup('tests/error-ins/error-17.txt',ValueError)
if __name__ == '__main__':
    unittest.main()

