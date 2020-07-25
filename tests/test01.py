import unittest
import pygrad

class TestSetup(unittest.TestCase):
    def testErrorFileSetup(self,fname,error):
        with self.assertRaises(error):
            pygrad.setup(fname)

    def testErrorFileSetupFormat(self,fname):
        self.testErrorFileSetup(fname, InfileFormatException)

    def testMissingLine(self):
        self.testErrorFileSetupFormat('error-ins/error-01.txt')

    def testMissingFile(self):
        self.testErrorFileSetup('error-ins/error-00.txt', FileNotFoundError)
    
    def testBlankFile(self):
        self.testErrorFileSetupFormat('error-ins/error-02.txt')

    def testNoGas(self):
        self.testErrorFileSetupFormat('error-ins/error-03.txt')

    def testAbsZero(self):
        self.testErrorFileSetupFormat('error-ins/error-04.txt')

    def testZeroPressure(self):
        self.testErrorFileSetupFormat('error-ins/error-05.txt')

    def testTooFewParameters(self):
        self.testErrorFileSetupFormat('error-ins/error-06.txt')

    def testTooManyParameters(self):
        self.testErrorFileSetupFormat('error-ins/error-07.txt')
    
    def testNgasBig(self):
        self.testErrorFileSetupFormat('error-ins/error-08.txt')

    def testNgasMismatch(self):
        self.testErrorFileSetupFormat('error-ins/error-09.txt')
if __name__ == '__main__':
    unittest.main()

