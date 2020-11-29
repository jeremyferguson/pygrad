import unittest,os
from Pygrad.Pygrad import Pygrad
from Pygrad import utils,Mixers

#A class to test the MixerC class
class TestMixer(unittest.TestCase):
    pygrad_home = os.getenv('PYGRAD_HOME')
    mixerTestPath = pygrad_home + '/tests/fortran_tests/mixer-tests/'

    def testMixer(self):
        def function(fname):
            print('init')
            obj = utils.createObject(self.mixerTestPath,fname)
            print(obj.MixObject)
            Mixers.Mixer(obj)
            return [obj.ElasticCrossSection,obj.CrossSectionSum,obj.AttachmentSectionSum,obj.TotalCrossSection,obj.RelativeIonMinusAttachCrossSection,obj.InelasticCrossSection,obj.ElasticCrossSection]
        utils.checkFortranTests(self.mixerTestPath, function,"01",3)

if __name__ == '__main__':
    unittest.main()
