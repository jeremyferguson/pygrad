import unittest
import pygrad

class TestSetup(unittest.TestCase):
    def testMissingLine(self):
        with self.assertRaises(InfileFormatException):
            pygrad.setup('error-ins/error-01.txt')

    def testMissingFile(self):
        with self.assertRaises(FileNotFoundError):
            pygrad.setup('error-ins/error-00.txt')
    
    def testBlankFile(self):
        with self.assertRaises(InfileFormatException):
            pygrad.setup('error-ins/error-02.txt')

    def testNoGas(self):
        with self.assertRaises(InfileFormatException):
            pygrad.setup('error-ins/error-03.txt')

    def testAbsZero(self):
        with self.assertRaises(InfileFormatexception):
            pygrad.setup('error-ins/error-04.txt')

if __name__ == '__main__':
    unittest.main()

