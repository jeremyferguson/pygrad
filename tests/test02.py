import unittest,os
from Pygrad.Pygrad import Pygrad
from Pygrad import utils,MixerC

#A class to test the MixerC class
class TestMixerC(unittest.TestCase):
    pygrad_home = os.getenv('PYGRAD_HOME')
    mixercTestPath = pygrad_home + '/tests/fortran_tests/mixerc-tests/'

    #Returns a generic function for passing into the checkFortranTests function.
    def mixerCStart(self,keys):
        def start(fname):
            obj = utils.createObject(self.mixercTestPath,fname)
            MixerC.MixerC(obj)
            return [obj.mixercArrs[key] for key in keys]
        return start

    def testMixerC1(self):
        MixerCCalc = self.mixerCStart(['RAD','PRSH'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"01", 6)

    def testMixerC2(self):
        MixerCCalc = self.mixerCStart(['ESH','PRSHBT','INIOCC'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"02", 6)

    def testMixerC3(self):
        MixerCCalc = self.mixerCStart(['IZ','ISHLMX','AMZ'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"03", 6)
    
    def testMixerC4(self):
        MixerCCalc = self.mixerCStart(['AUG'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"04", 6)

    def testMixerC5(self):
        MixerCCalc = self.mixerCStart(['XPE','YPE'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"05", 6)
   
    def testMixerC6(self):
        MixerCCalc = self.mixerCStart(['XCP','YRY','YCP','YPP'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"06", 6)
    
    def testMixerC7(self):
        MixerCCalc = self.mixerCStart(['FRMFR','FRMFC'])
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"07", 6)

if __name__ == '__main__':
    unittest.main()
