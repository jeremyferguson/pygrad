import unittest,writeh5,os
import numpy as np

#A class to test the gas functions stored by writeh5.py
class TestGasFuncs(unittest.TestCase):
    pygrad_home = os.getenv('PYGRAD_HOME')

    def checksum(self,gas,key):
        self.assertGreater(np.sum(writeh5.arrays[gas][key]),0)

    def testgas1(self):
        arrays = ['XEN','YELM','YELT','YEPS','XVBV4','YVBV4','XVBV3','YVBV3','XVIB5','YVIB5','XVIB6','YVIB6','XCF3','YCF3','XCF2','YCF2','XCF1','YCF1','XCF32','YCF32','XCF0','YCF0','XC0F','YC0F','XCF22','YCF22','XCF','YCF','XCFF','YCFF','XCF2F','YCF2F','XCF3F','YCF3F','XKSHC','YKSHC','XKSHF','YKSHF','XTR1','YTR1','XTR2','YTR2','XTR3','YTR3','XATT','YATT','Z6T','Z9T','EBRM','IZBR','KEL','KIN','E','EIN','EION','EOBY','LEGAS','ISHELL','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2']
        for arr in arrays:
            self.checksum(1,arr)

    def testgas2(self):
        pass

if __name__ == '__main__':
    unittest.main
