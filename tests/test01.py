import unittest
import pygrad

class TestSetup(unittest.TestCase):
    def testMissingLine(self):
        pygrad.setup('error-ins/error-01.txt')

    def testMissingFile(self):
        with self.assertRaises(FileNotFoundError):
            pygrad.setup('error-ins/error-00.txt')
    
    def testBlankFile(self):
        with self.assertRaises(InfileFormatException):
            pygrad.setup('error-ins/error-01.txt')

if __name__ == '__main__':
    unittest.main()

