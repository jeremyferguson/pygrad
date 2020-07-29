import unittest,os
import pygrad,utils,mixerc

class TestMixerC(unittest.TestCase):
    pygrad_home = os.getenv('PYGRAD_HOME')
    mixercTestPath = pygrad_home + '/tests/fortran_tests/mixerc-tests/'

    def testMixerC1(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['RAD'],mix.arrays['PRSH']]

        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"01", 6)

    def testMixerC2(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['ESH'],mix.arrays['PRSHBT'],mix.arrays['INIOCC']]

        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"02", 6)

    def testMixerC3(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['IZ'],mix.arrays['ISHLMX'],mix.arrays['AMZ']]

        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"03", 6)
    
    def testMixerC4(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['AUG']]
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"04", 6)

    def testMixerC5(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['XPE'],mix.arrays['YPE']]
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"05", 6)
   
    def testMixerC6(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['XCP'],mix.arrays['YRY'],mix.arrays['YCP'],mix.arrays['YPP']]
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"06", 6)
    
    def testMixerC7(self):
        def MixerCCalc(fname):
            main = pygrad.Main(self.mixercTestPath + fname)
            mix = mixerc.MixerC(main)
            mix.mixc()
            return [mix.arrays['FRMFR'],mix.arrays['FRMFC']]
        utils.checkFortranTests(self.mixercTestPath, MixerCCalc,"07", 6)

if __name__ == '__main__':
    unittest.main()
