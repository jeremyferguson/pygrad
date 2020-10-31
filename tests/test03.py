import unittest,writeh5,os
import numpy as np

#A class to test the gas functions stored by writeh5.py
class TestGasFuncs(unittest.TestCase):
    pygrad_home = os.getenv('PYGRAD_HOME')

    @classmethod
    def setUpClass(cls):
        writeh5.main()

    def checkpartial(self,gas,key,start,mid,end):
        self.assertLess(np.sum(writeh5.arrays[gas][key][start:mid]),0,'{}[{}:{}]'.format(key,start,mid))
        self.assertGreater(np.sum(writeh5.arrays[gas][key][mid:end]),0,'{}[{}:{}]'.format(key,mid,end))
    
    def checksum(self,gas,key):
        self.assertGreater(np.sum(writeh5.arrays[gas][key]),0,key)

    def checksums(self,gas,arrays):
        for arr in arrays:
            self.checksum(gas,arr)

    def checklen(self,gas,key,l):
        self.assertEqual(len(writeh5.arrays[gas][key]),l,'Incorrect {} assignment'.format(key))

    def testgas1(self):
        arrays = ['XEN','YELM','YELT','YEPS','XVBV4','YVBV4','XVBV3','YVBV3','XVIB5','YVIB5','XVIB6','YVIB6','XCF3','YCF3','XCF2','YCF2','XCF1','YCF1','XCF32','YCF32','XCF0','YCF0','XC0F','YC0F','XCF22','YCF22','XCF','YCF','XCFF','YCFF','XCF2F','YCF2F','XCF3F','YCF3F','XKSHC','YKSHC','XKSHF','YKSHF','XTR1','YTR1','XTR2','YTR2','XTR3','YTR3','XATT','YATT','Z6T','Z9T','EBRM','IZBR','KEL','KIN','E','EIN','EION','EOBY','LEGAS','ISHELL','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(1,arrays)
        self.checklen(1,'AP',12)
        self.checklen(1,'FA',0)
        self.assertEqual(writeh5.variables[1]['GAMMACOND'],[False,False])

    def testgas2(self):
        arrays = ['XEN','YSEC','YEL','XEPS','YEPS','XENI','YENI','YENC','YEN1','XEN2','YEN2','XEN3','YEN3','XKSH','YKSH','XL1S','YL1S','XL2S','YL2S','XL3S','YL3S','X1S5','Y1S5','X1S4','Y1S4','X1S3','Y1S3','X1S2','Y1S2','X2P10','Y2P10','X2P9','Y2P9','X2P8','Y2P8','X2P7','Y2P7','X2P6','Y2P6','X2P5','Y2P5','X2P4','Y2P4','X2P3','Y2P3','X2P2','Y2P2','X2P1','Y2P1','X3D6','Y3D6','X3D5','Y3D5','X3D4P','Y3D4P','X3D4','Y3D4','X3D3','Y3D3','X3D1PP','Y3D1PP','X3D1P','X3D1P','X3S1PPPP','Y3S1PPPP','X3S1PPP','Y3S1PPP','X3S1PP','Y3S1PP','X2S5','Y2S5','X2S3','Y2S3','Z18T','EBRM','IZBR','KEL','KIN','E','EOBY','EION','EIN','LEGAS','ISHELL','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(2,arrays)
        self.checklen(2,'FA',0)
        self.checklen(2,'AP',0)
        self.assertEqual(writeh5.variables[2]['GAMMACOND'],[True,False])

    def testgas3(self):
        arrays = ['XEN','YEM','YEL','YEPS','XION','YION','YINC','X23S','Y23S','X21S','Y21S','X23P','Y23P','X21P','Y21P','X33S','Y33S','X31S','Y31S','X33P','Y33P','X33D','Y33D','X31D','Y31D','X31P','Y31P','X43S','Y43S','X41S','Y41S','X43P','Y43P','X43D','Y43D','X41D','Y41D','X43F','Y43F','X41F','Y41F','X41P','Y41P','Z2T','EBRM','IZBR','KEL','KIN','E','EOBY','EION','EIN','EC0','PENFRA']
        self.checksums(3,arrays)
        self.checklen(3,'FA',0)
        self.checklen(3,'AP',0)
        self.assertEqual(writeh5.variables[3]['GAMMACOND'],[True,False])

    def testgas4(self):
        arrays = ['XEN','YEM','YEL','YEPS','XION','YION','YINC','X23S','Y23S','X21S','Y21S','X23P','Y23P','X21P','Y21P','X33S','Y33S','X31S','Y31S','X33P','Y33P','X33D','Y33D','X31D','Y31D','X31P','Y31P','X43S','Y43S','X41S','Y41S','X43P','Y43P','X43D','Y43D','X41D','Y41D','X43F','Y43F','X41F','Y41F','X41P','Y41P','Z2T','EBRM','IZBR','KEL','KIN','E','EOBY','EION','EIN','EC0','PENFRA']
        self.checksums(4,arrays)
        self.checklen(4,'FA',0)
        self.checklen(4,'AP',0)
        self.assertEqual(writeh5.variables[4]['GAMMACOND'],[True,False])

    def testgas5(self):
        arrays = ['XEN','YXSEC','XEL','YEL','XEPS','YEPS','XION','YION','YINC','YIN1','XIN2','YIN2','XIN3','YIN3','XKSH','YKSH','X1S5','Y1S5','X1S4','Y1S4','X1S3','Y1S3','X1S2','Y1S2','X2P10','Y2P10','X2P9','Y2P9','X2P8','Y2P8','X2P7','Y2P7','X2P6','Y2P6','X2P5','Y2P5','X2P4','Y2P4','X2P3','Y2P3','X2P2','Y2P2','X2P1','Y2P1','X2S5','Y2S5','X2S3','Y2S3','X3D6','Y3D6','X3D4P','Y3D4P','X3D4','Y3D4','X3D3','Y3D3','X3D1PP','Y3D1PP','X3D1P','Y3D1P','X3S1PPPP','Y3S1PPPP','X3S1PPP','Y3S1PPP','X3S1PP','Y3S1PP','X3P106','Y3P106','X3P52','Y3P52','X3P1','Y3P1','Z10T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(5,arrays)
        self.checklen(5,'FA',0)
        self.checklen(5,'AP',0)
        self.assertEqual(writeh5.variables[5]['GAMMACOND'],[True,False])

    def testgas6(self):
        arrays = ['XEN','YXSEC','XEL','YEL','XEPS','YEPS','XION','YION','YINC','YIN1','XIN2','YIN2','XIN3','YIN3','XIN4','YIN4','XKSH','YKSH','XL1S','YL1S','XL2S','YL2S','XL3S','YL3S','XM1S','YM1S','XM2S','YM2S','XM3S','YM3S','XM4S','YM4S','XM5S','YM5S','X1S5','Y1S5','X1S2','Y1S2','X2P10','Y2P10','X2P9','Y2P9','X2P8','Y2P8','X2P7','Y2P7','X2P6','Y2P6','X2P5','Y2P5','X3D6','Y3D6','X3D5','Y3D5','X2P4','Y2P4','X3D3','Y3D3','X3D4P','Y3D4P','X2P3','Y2P3','X2P2','Y2P2','X3D4','Y3D4','X2P1','Y2P1','X3D1PP','Y3D1PP','X3D1P','Y3D1P','Y3D1P','X2S5','Y2S5','X3P10','Y3P10','X3P9','Y3P9','X3P8','X3S1PP','Y3S1PP','X3P7','Y3P7','X3P6','Y3P6','X3S1PPPP','Y3S1PPPP','X3S1PPP','Y3S1PPP','X3P5','Y3P5','X4D6','Y4D6','X4D4P','Y4D4P','X4D4','Y4D4','X4D3','Y4D3','X2S3','Y2S3','X4D1PP','Y4D1PP','X4D1P','Y4D1P','X3S5','Y3S5','X4FS','Y4FS','Z36T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(6,arrays)
        self.checklen(6,'FA',0)
        self.checklen(6,'AP',0)
        self.assertEqual(writeh5.variables[6]['GAMMACOND'],[True,False])

    def testgas7(self):
        arrays = ['XEN','YMOM','XEL','YEL','XEPS','YEPS','XION','YION','YINC','YIN1','XIN2','YIN2','XIN3','YIN3','XIN4','YIN4','XIN5','YIN5','XIN6','YIN6','XKSH','YKSH','XL1S','YL1S','XL2S','YL2S','XL3S','YL3S','XM1S','YM1S','XM2S','YM2S','XM3S','YM3S','XM4S','YM4S','XM5S','YM5S','X1S5','Y1S5','X1S4','Y1S4','X1S3','Y1S3','X1S2','Y1S2','X2P10','Y2P10','X2P9','Y2P9','X2P8','Y2P8','X2P7','Y2P7','X2P6','Y2P6','X3D6','Y3D6','X2P5','Y2P5','X3D4P','Y3D4P','X3D3','Y3D3','X3D4','Y3D4','X3D1PP','Y3D1PP','X3D1P','Y3D1P','X2S5','Y2S5','X3P105','Y3P105','X2P4','Y2P4','X4DSUM','Y4DSUM','X2P3','Y2P3','X2P2','Y2P2','X2P1','Y2P1','Z54T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(7,arrays)
        self.checklen(7,'FA',0)
        self.checklen(7,'AP',0)
        self.assertEqual(writeh5.variables[7]['GAMMACOND'],[True,False])

    def testgas8(self):
        arrays = ['XEN','YELM','YELT','YEPS','XATT','YATT','XVBV4','YVBV4','XVBV2','YVBV2','XVBV1','YVBV1','XVBV3','YVBV3','XVBH1','YVBH1','XVBH2','YVBH2','XION','YION','YINC','XINF','YINF','XINF1','YINF1','XINF2','YINF2','XINF3','YINF3','XINF4','YINF4','XINF5','YINF5','XINF6','YINF6','XINPP','YINPP','XDET','YDET','XTR1','YTR1','XTR2','YTR2','XTR3','YTR3','XCHD','YCHD','XCHB','YCHB','XHAL','YHAL','XHBE','YHBE','XKSH','YKSH','Z1T','Z6T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2']
        self.checksums(8,arrays)
        self.checklen(8,'FA',0)
        self.checklen(8,'AP',8)
        self.assertEqual(np.sum(writeh5.arrays[8]['PENFRA']),0)
        self.assertEqual(writeh5.variables[8]['GAMMACOND'],[False,False])

    def testgas9(self):
        arrays = ['XEN','YMT','YEL','YEPS','XATT1','YATT1','XATT2','YATT2','XVIB1','YVIB1','XVIB2','YVIB2','XVIB3','YVIB3','XVIB4','YVIB4','XVIB5','YVIB5','XTR1','YTR1','XTR2','YTR2','XTR3','YTR3','XNUL1','YNUL1','XNUL2','YNUL2','XNUL3','YNUL3','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION','YIONG','YIONC','Z1T','Z6T','SCLN','ESPLIT','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2']
        self.checksums(9,arrays)
        self.checklen(9,'AP',4)
        self.checklen(9,'FA',0)
        self.assertEqual(np.sum(writeh5.arrays[9]['PENFRA']),0)
        self.assertEqual(writeh5.variables[9]['GAMMACOND'],[False,False])

    def testgas10(self):
        arrays = ['XEN','YMT','YEL','YEPS','XION','YIONG','YIONC','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION17','YION17','XION18','YION18','XION19','YION19','XION20','YION20','XION21','YION21','XION22','YION22','XION23','YION23','XION24','YION24','XATT1','YATT1','XATT2','YATT2','XVIB1','YVIB1','XVIB2','YVIB2','XVIB3','YVIB3','XVIB4','YVIB4','XTR1','YTR1','XTR2','YTR2','XTR3','YTR3','XTR4','YTR4','XNUL1','YNUL1','XNUL2','YNUL2','Z1T','Z6T','ESPLIT','SCLN','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2']
        self.checksums(10,arrays)
        self.checklen(10,'AP',3)
        self.checklen(10,'FA',0)
        self.assertEqual(np.sum(writeh5.arrays[10]['PENFRA']),0)
        self.assertEqual(writeh5.variables[10]['GAMMACOND'],[False,False])

    def testgas11(self):
        arrays = ['XEN','YELM','YELT','YEPS','XION','YION','XKSH','YKSH','XVIB1','YVIB1','XVIB2','YVIB2','XVIB3','YVIB3','XVIB4','YVIB4','XVIB5','YVIB5','XEXC1','YEXC1','XEXC2','YEXC2','Z6T','Z1T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(11,arrays)
        self.checklen(11,'AP',34)
        self.checklen(11,'FA',0)
        self.assertEqual(writeh5.variables[11]['GAMMACOND'],[False,False])

    def testgas12(self):
        arrays = ['XEN','YMOM','YEL','YVBMOM','YVBEL','YEPS','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XATT','YATT','XV2','YV2','X2V2','Y2V2','XV1','YV1','X3V2','Y3V2','XV3','YV3','XVPD3','YVPD3','XV130','YV130','XVPD4','YVPD4','XVPD5','YVPD5','XVPD6','YVPD6','XVPD7','YVPD7','XVPD8','YVPD8','XVPD9','YVPD9','XVPDH','YVPDH','XTRP1','YTRP1','XTRP2','YTRP2','XKSHC','YKSHC','XKSHO','YKSHO','Z6T','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(12,arrays)
        self.assertLess(np.sum(writeh5.arrays[12]['EIN'][1:60:2]),0,'EIN first half')
        self.assertGreater(np.sum(writeh5.arrays[12]['EIN'][:60:2]),0,'EIN first half')
        self.checklen(12,'PJ',2)
        self.checklen(12,'FA',0)
        self.checklen(12,'AP',13)
        self.assertEqual(writeh5.variables[12]['GAMMACOND'],[False,False])

    def testgas14(self):
        arrays = ['ELEV','AJL','SALPHA','EROT','AJIN','IMAP','XEL','YEL','XMT','YMT','XEPS','YEPS','XVIB1','YVIB1','XVIB2','YVIB2','XVIB3','YVIB3','XION','YIONC','YIONG','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XKSH','YKSH','XATT1','YATT1','XATT2','YATT2','XATT3','YATT3','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','XNUL1','YNUL1','XNUL2','YNUL2','XNUL3','YNUL3','XNUL4','YNUL4','ENROT','ENRTS','YEPSR','YMTRT','Z1T','Z8T','SCLN','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(14,arrays)
        self.checklen(14,'FA',0)
        self.assertGreater(np.sum(writeh5.arrays[14]['EIN'][:200:2]),0,'EIN first half')
        self.assertLess(np.sum(writeh5.arrays[14]['EIN'][1:200:2]),0,'EIN first half')
        self.checklen(14,'PJ',2)
        self.checklen(14,'AP',4)
        self.assertEqual(writeh5.variables[14]['GAMMACOND'],[False,False])

    def testgas15(self):
        arrays = ['XELA','YELA','YMOM','YEPS','XROT13','YROT13','XROT35','YROT35','XROT57','YROT57','XROT79','YROT79','XROT911','YROT911','XROT1113','YROT1113','XROT1315','YROT1315','XROT1517','YROT1517','XROT1719','YROT1719','XROT1921','YROT1921','XROT2123','YROT2123','XROT2325','YROT2325','XROT2527','YROT2527','XROT2729','YROT2729','XROT2931','YROT2931','XROT3133','YROT3133','XROT3335','YROT3335','XROT3537','YROT3537','XROT3739','YROT3739','XROT3941','YROT3941','XROT4143','YROT4143','XROT4345','YROT4345','XROT4547','YROT4547','XROT4749','YROT4749','XVIB','YVIB1','YVIB2','YVIB3','YVIB4','YVIB5','YVIB6','YVIB7','YVIB8','YVIB9','YVIB10','YVIB11','YVIB12','YVIB13','YVIB14','YVIB15','YVIB16','YVIB17','YVIB18','YVIB19','YVIB20','YVIB21','X3ATT','Y3ATT','XATT','YATT','XEXC1','YEXC1','XEXC2','YEXC2','XEXC3','YEXC3','XEXC4','YEXC4','XEXC5','YEXC5','XEXC6','YEXC6','XEXC7','YEXC7','XEXC8','YEXC8','XEXC9','YEXC9','XROT','YROT','XIONC','YIONC','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XKSH','YKSH','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(15,arrays)
        self.checklen(15,'PJ',2)
        self.checklen(15,'AP',1)
        self.checklen(15,'FA',1)
        self.assertLess(writeh5.arrays[15]['EIN'][0],0,'EIN[0]')
        self.assertGreater(writeh5.arrays[15]['EIN'][24],0,'EIN[24]')
        self.assertGreater(np.sum(writeh5.arrays[15]['EIN'][25:48]),0,'EIN[24:48]')
        self.assertLess(np.sum(writeh5.arrays[15]['EIN'][1:24]),0,'EIN[1:24]')
        self.assertEqual(writeh5.variables[15]['GAMMACOND'],[False,False])

    def testgas16(self):
        arrays = ['XELA','YELA','YMOM','YEPS','XROT','YROT','XVB1','YVB1','XVB2','YVB2','XVB3','YVB3','XVB4','YVB4','XVB5','YVB5','XVB6','YVB6','XVB7','YVB7','XVB8','YVB8','XVB9','YVB9','XVB10','YVB10','XVB11','YVB11','XVB12','YVB12','XVB13','YVB13','XVB14','YVB14','XVB15','YVB15','XTRP1','YTRP1','YTP1M','XTRP2','YTRP2','YTP2M','XTRP3','YTRP3','YTP3M','XTRP4','YTRP4','YTP4M','XTRP5','YTRP5','YTP5M','XTRP6','YTRP6','YTP6M','XTRP7','YTRP7','YTP7M','XTRP8','YTRP8','YTP8M','XTRP9','YTRP9','YTP9M','XTRP10','YTRP10','YTP10M','XTRP11','YTRP11','YTP11M','XTRP12','YTRP12','YTP12M','XTRP13','YTRP13','YTP13M','XTRP14','YTRP14','YTP14M','XSNG1','YSNG1','YSG1M','XSNG2','YSNG2','YSG2M','XSNG3','YSNG3','YSG3M','XSNG4','YSNG4','YSG4M','XSNG5','YSNG5','YSG5M','XSNG6','YSNG6','YSG6M','XSNG7','YSNG7','YSG7M','XSNG8','YSNG8','YSG8M','XSNG9','YSNG9','YSG9M','XSNG10','YSNG10','YSG10M','XSNG11','YSNG11','YSG11M','XSNG12','YSNG12','YSG12M','XSNG13','YSNG13','YSG13M','XSNG14','YSNG14','YSG14M','XSNG15','YSNG15','YSG15M','XKSH','YKSH','XION','YION','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','Z7T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(16,arrays)
        self.checklen(16,'PJ',2)
        self.checklen(16,'FA',0)
        self.checklen(16,'AP',6)
        self.checkpartial(16,'EIN',0,38,76)
        self.assertEqual(writeh5.variables[16]['GAMMACOND'],[False,False])

    def testgas18(self):
        arrays = ['XEN','YELTT','YELMT','XEPS','YEPS','XATT','YATT','XVIBH','YVIBH','XVIBR','YVIBR','XV020','YV020','XV110','YV110','XV120','YV120','XV200','YV200','XV011','YV011','XV210','YV210','XV101','YV101','XV300','YV300','XV002','YV002','XV201','YV201','XV400','YV400','XV500','YV500','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','ENROT','ENRTS','YEPSR','YMTRT','XION','YIONC','YIONG','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','Z7T','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(18,arrays)
        self.checklen(18,'FA',0)
        self.checklen(18,'AP',8)
        self.checkpartial(18,'EIN',0,60,76)
        self.checklen(18,'PJ',1)
        self.checklen(18,'FROT',1)
        self.assertEqual(writeh5.variables[18]['GAMMACOND'],[False,False])

    def testgas21(self):
        arrays = ['XELM','YELM','YELT','YEPS','XROT0','YROT0','XROT1','YROT1','XROT2','YROT2','XROT3','YROT3','XROT4','YROT4','XROT5','YROT5','XROT6','YROT6','XROT7','YROT7','XVIB1','YVIB1','XVIB2','YVIB2','XVIB3','YVIB3','XVIB4','YVIB4','XB3S1','YB3S1','XB3S2','YB3S2','XB3S3','YB3S3','XB3S4','YB3S4','XC3PI','YC3PI','XA3SG','YA3SG','XE3SG','YE3SG','XEFSG','YEFSG','XATT','YATT','XION','YION','XIOND','YIOND','XIONA','YIONA','XION0','YION0','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION17','YION17','AMP','NIONN','ERLVL','BEF','DISLY','DISWR','DISD1P','DISB1S','Z1T','EBRM','E','EOBY','EION','EIN','IZBR','KEL','KIN','PENFRA']
        self.checksums(21,arrays)
        self.checklen(21,'PJ',2)
        self.checklen(21,'AP',0)
        self.checklen(21,'FA',0)
        self.checklen(21,'FROT',10)
        self.assertGreater(np.sum(writeh5.arrays[21]['ERLVL'][1:8:2]),0,"ERLVL[1:8:2]")
        self.assertEqual(writeh5.variables[21]['GAMMACOND'],[True,True])

    def testgas30(self):
        arrays= ['XEN','YELM','YELT','YEPS','XATT','YAT1','YAT2','YAT3','YAT4','YAT5','YAT6','YAT7','XION','YION','YIN1','YIN2','YIN3','YIN4','YIN5','YIN6','YIN7','XL3SH','YL3SH','XL2SH','YL2SH','XL1SH','YL1SH','XKSHS','YKSHS','XKSHF','YKSHF','XV1V1','YV1V1','XV2V1','YV2V1','XV3V1','YV3V1','XV4V1','YV4V1','XV5V1','YV5V1','XVBV3','YVBV3','XTR1','YTR1','XTR2','YTR2','XTR3','YTR3','Z9T','Z16T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2']
        self.checksums(30,arrays)
        self.checklen(30,'FA',0)
        self.checklen(30,'AP',10)
        self.assertGreater(np.sum(writeh5.arrays[30]['KIN'][:9]),0,'KIN[:9]')
        self.assertEqual(np.sum(writeh5.arrays[30]['PENFRA']),0)
        self.assertEqual(writeh5.variables[30]['GAMMACOND'],[False,False])

    def testgas31(self):
        arrays = ['XEL','YEL','YMT','YEPS','XVIBH','YVIBH','XATT','YATT1','YATT2','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','XTRP5','YTRP5','XTRP6','YTRP6','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XKSH','YKSH','XIONC','YIONC','ELEV','AKL','AJL','ENROT','ENRTS','YEPSR','YMTRT','Z1T','Z7T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(31,arrays)
        self.checklen(31,'FA',0)
        self.checklen(31,'AP',1)
        self.checklen(31,'PJ',2)
        self.assertEqual(writeh5.variables[31]['GAMMACOND'],[False,False])

    def testgas34(self):
        arrays = ['XEL','YEL','YMT','YEPS','XVIBH','YVIBH','XATT1','YATT1','XATT2','YATT2','XATT3','YATT3','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','XTRP5','YTRP5','XTRP6','YTRP6','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION17','YION17','XIONC','YIONC','YIONG','ENROT','ENRTS','YEPSR','YMTRT','Z1T','Z6T','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(34,arrays)
        self.checklen(34,'PJ',1)
        self.checklen(34,'FA',0)
        self.checklen(34,'AP',1)
        self.checklen(34,'FROT',1)
        self.checkpartial(34,'EIN',0,26,52)
        self.assertEqual(writeh5.variables[34]['GAMMACOND'],[False,False])

    def testgas35(self):
        arrays = ['XEL','YEL','YMT','YEPS','XVIBH','YVIBH','XATT1','YATT1','XATT2','YATT2','XATT3','YATT3','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','XTRP5','YTRP5','XTRP6','YTRP6','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION17','YION17','XION18','YION18','XION19','YION19','XION20','YION20','XION21','YION21','XION22','YION22','XION23','YION23','XION24','YION24','XION25','YION25','XION26','YION26','XION27','YION27','XIONC','YIONC','YIONG','ENROT','ENRTS','YEPSR','YMTRT','Z1T','Z6T','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(35,arrays)
        self.checklen(35,'AP',1)
        self.checklen(35,'PJ',1)
        self.checklen(35,'FA',0)
        self.checklen(35,'FROT',1)
        self.checkpartial(35,'EIN',0,26,52)
        self.assertEqual(writeh5.variables[35]['GAMMACOND'],[False,False])

    def testgas36(self):
        arrays = ['XEL','YEL','YMT','YEPS','XVIBH','YVIBH','XATT1','YATT1','XATT2','YATT2','XATT3','YATT3','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','XTRP5','YTRP5','XTRP6','YTRP6','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION17','YION17','XION18','YION18','XION19','YION19','XION20','YION20','XION21','YION21','XION22','YION22','XION23','YION23','XION24','YION24','XION25','YION25','XION26','YION26','XION27','YION27','XION28','YION28','XION29','YION29','XIONC','YIONC','YIONG','ENROT','ENRTS','YEPSR','YMTRT','Z1T','Z6T','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(36,arrays)
        self.checklen(36,'FA',0)
        self.checklen(36,'AP',1)
        self.checklen(36,'PJ',1)
        self.checkpartial(36,'EIN',0,26,52)
        self.assertEqual(writeh5.variables[36]['GAMMACOND'],[False,False])

    def testgas44(self):
        arrays = ['XEN','YELM','YELT','YEPS','XION','YION','XKSHC','YKSHC','XKSHN','YKSHN','XTORS','YTORS','XVIB1','YVIB1','XVIB2','YVIB2','XVIB3','YVIB3','XVHAR','YVHAR','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','Z6T','Z7T','Z1T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2']
        self.checksums(44,arrays)
        self.checklen(44,'FA',0)
        self.checklen(44,'AP',3)
        self.assertEqual(np.sum(writeh5.arrays[44]['PENFRA']),0)
        self.assertEqual(writeh5.variables[44]['GAMMACOND'],[False,False])

    def testgas46(self):
        arrays = ['XEL','YEL','YMT','YEPS','XVIBH','YVIBH','XATT1','YATT1','XATT2','YATT2','XATT3','YATT3','XTRP1','YTRP1','XTRP2','YTRP2','XTRP3','YTRP3','XTRP4','YTRP4','XTRP5','YTRP5','XTRP6','YTRP6','XION1','YION1','XION2','YION2','XION3','YION3','XION4','YION4','XION5','YION5','XION6','YION6','XION7','YION7','XION8','YION8','XION9','YION9','XION10','YION10','XION11','YION11','XION12','YION12','XION13','YION13','XION14','YION14','XION15','YION15','XION16','YION16','XION17','YION17','XION18','YION18','XION19','YION19','XION20','YION20','XION21','YION21','XION22','YION22','XION23','YION23','XION24','YION24','XION25','YION25','XION26','YION26','XION27','YION27','XION28','YION28','XION29','YION29','XIONC','YIONC','YIONG','ENROT','ENRTS','YEPSR','YMTRT','Z1T','Z6T','Z8T','EBRM','E','EOBY','EION','EIN','LEGAS','ISHELL','IZBR','KEL','KIN','NC0','EC0','WKLM','EFL','NG1','NG2','EG1','EG2','PENFRA']
        self.checksums(46,arrays)
        self.checklen(46,'FA',0)
        self.checklen(46,'AP',1)
        self.checklen(46,'PJ',1)
        self.assertEqual(writeh5.variables[46]['GAMMACOND'],[False,False])

if __name__ == '__main__':
    unittest.main()
