import pygrad, os, h5py, numpy as np
class Mixer():
    def __init__(self,main):
        self.main = main
        #common
        self.description = ""
        self.qelm = np.zeros(20000,dtype = float)
        self.qsum = np.zeros(20000,dtype = float)
        self.qin = np.zeros((6,250,20000),dtype = float)
        self.qion = np.zeros((6,20000),dtype = float)
        self.qsatt = np.zeros(20000,dtype = float)
        self.qtot = np.zeros(20000,dtype = float)
        self.qrel = np.zeros(20000,dtype = float)
        self.qinel = np.zeros(20000,dtype = float)
        self.qel = np.zeros(20000,dtype = float)
        self.nin = np.zeros(6,dtype = int)
        self.lion = np.zeros(6,dtype = int)
        self.lin = np.zeros((6,250),dtype = int)
        self.alion = np.zeros(6,dtype = float)
        self.alin = np.zeros((6,250),dtype = float)
        self.cf = np.zeros((20000,512),dtype = float)
        self.ein = np.zeros(512,dtype = float)
        self.tcf = np.zeros(20000,dtype = float)
        self.iarry = np.zeros(512,dtype = int)
        self.rgas = np.zeros(512,dtype = float)
        self.ipn = np.zeros(512,dtype = int)
        self.wpl = np.zeros(512,dtype = float)
        self.izbr = np.zeros(512,dtype = int)
        self.penfra = np.zeros((3,512),dtype = float)
        self.iplast = 0
        self.cfn = np.zeros((20000,60),dtype = float)
        self.tcfn = np.zeros(20000,dtype = float)
        self.sclenul = np.zeros(60,dtype = float)
        self.nplast = 0
        self.psct = np.zeros((20000,512),dtype = float)
        self.angct = np.zeros(20000,512),dtype = float)
        self.index = np.zeros(512, dtype = int)
        self.niso = 0
        self.fcion = np.zeros(20000,dtype = float)
        self.fcatt = np.zeros(20000,dtype = float)
        self.negas = np.zeros(512,dtype = int)
        self.legas = np.zeros(512,dtype = int)
        self.ieshell = np.zeros(512,dtype = int)
        self.iecasc = 0
        self.doubles = np.zeros((6,20000), dtype = float)
        self.cminixsc = np.zeros(6,dtype = float)
        self.cminexsc = np.zeros(6,dtype = float)
        self.ecloss = np.zeros(6,dtype = float)
        self.wpln = np.zeros(6,dtype = float)
        self.avpfrac = np.zeros((3,6),dtype = float)
        self.nc0 = np.zeros(512,dtype = int)
        self.ec0 = np.zeros(512,dtype = float)
        self.ng1 = np.zeros(512,dtype = int)
        self.eg1 = np.zeros(512,dtype = float)
        self.ng2 = np.zeros(512,dtype = int)
        self.wklm = np.zeros(512,dtype = float)
        self.efl = np.zeros(512,dtype = float)
        self.ngexc = np.zeros(6,dtype = int)
        self.idg = np.zeros(6,dtype = int)
        self.esplit = np.zeros((512,20),dtype = float)
        self.ionmodel = np.zeros(512,dtype = int)

        #Only used in mixer program
        self.q = np.zeros((6,6,20000),dtype = float)
        self.e = np.zeros((6,6),dtype = float)
        self.ei = np.zeros((6,250),dtype = float)
        self.qatt = np.zeros((6,20000),dtype = float)
        self.eion = np.zeros(6,dtype = float)
        self.peqel = np.zeros((6,6,20000),dtype = float)
        self.peqin = np.zeros((6,250,20000),dtype = float)
        self.penfrac = np.zeros((6,3,250),dtype = float)  
        self.kin = np.zeros((6,250),dtype = int)
        self.kel = np.zeros((6,6),dtype = int)
        self.eions = np.zeros((6,30),dtype = float)
        self.qions = np.zeros((6,30,20000),dtype = float)
        self.peqions = np.zeros(6,30,20000),dtype = float)
        self.legasn = np.zeros((6,30),dtype = int)
        self.ieshel = np.zeros((6,30),dtype = int)
        self.eb = np.zeros((6,30),dtype = float)
        self.nc = np.zeros((6,30),dtype = int)
        self.ec = np.zeros((6,30),dtype = float)
        self.ng = np.zeros((6,30),dtype = int)
        self.eg = np.zeros((6,30),dtype = float)
        self.ngs = np.zeros((6,30),dtype = int)
        self.egs = np.zeros((6,30),dtype = float)
        self.wk = np.zeros((6,30),dtype = float)
        self.efls = np.zeros((6,30),dtype = float)
        self.izbrs = np.zeros((6,250),dtype = int)
        self.qatts = np.zeros((6,8,20000),dtype = float)
        self.qnul = np.zeros((6,10,20000),dtype = float)
        self.scln = np.zeros((6,10), dtype = float)
        self.esplits = np.zeros((6,5,20),dtype = float)
        self.nion = np.zeros(6,dtype = int)
        self.natt = np.zeros(6,dtype = int)
        self.nul = np.zeros(6,dtype = int)


    def mix(self):
        i = 0
        for gas in self.main.ngasn:
            self.mixgas(gas,i)
        n = 0
        self.idg[0] = n
        self.negas[n] = 1 #index?
        self.cf[:,0] = self.q[0,1,:] * self.main.van[0] * self.main.bet
        self.psct[:,0] = 0.5
        self.angct[:,0] = 1.0
        acut = np.vectorize(utils.angcut)
        if kel[0,1] == 1:
            self.angct[:,0],self.psct[:,0] = acut(peqel[0,1,:])
            self.index[0] = 1 #index?
        elif kel[0,1] == 2:
            self.psct[:,0] = peqel[0,1,:]
            self.index[0] = 2 #index?
        rgas1 = 1.0 + self.e[0,1] / 2.0
        self.rgas[0] = rgas1
        L = 1
        self.iarry[0] = L #index?
        self.cminexsc[0] = self.e[0,3] * self.main.an[0]
        self.cminixsc[0] = self.e[0,4] * self.main.an[0]
        self.ecloss[0] = self.e[0,2]
        self.wpln[0] = self.e[0,5]
        if self.main.efinal >= self.e[0,2]: 
            if self.nion[0] > 1:
                start = n + 1
                n += self.nion[0]
                end = n + 1
                self.idg[0] = n
                self.cf[:,start:end] = self.qion[0,:self.nion[0],:] * self.main.van[0] * self.main.bet
                self.fcion += np.apply_along_axis(np.sum,0,self.cf[:,start:end])
                self.legas[start:end] = self.legasn[:self.nion[0]]
                self.ieshell[start:end] = self.ieshel[0,:self.nion[0]]
                peqarr = self.peqions[0,:self.nion[0],:]
                self.ein[start:end] = self.eion[0,:self.nion[0]]/rgas1
            else:
                n += 1
                end = n + 1
                self.idg[0] = n
                if self.main.icount == 1:
                    i = 4
                    self.doubles[0,:] = self.q[0,2,:] / self.q[0,4,:] - 1.0
                else:
                    i = 2
                self.cf[:,n] = self.q[0,i,:] * self.main.van[0] * self.main.bet
                self.fcion += self.cf[:,n]
                if self.main.icount == 1:
                    i = 4
                else:
                    i = 2
                peqarr = self.peqel[0,i,:]
                self.wpl[n] = self.eb[0,0]
                self.ein[n] = self.e[0,2]/rgas1
            length = end - start
            self.negas[start:end] = 1
            self.psct[:,start:end] = 0.5
            self.angct[:,start:end] = 1.0
            if self.kel[0,i] == 1:
                self.angct[:,start:end],self.psct[:,start:end] = acut(peqarr)
                self.index[start:end] = 1#index?
            elif self.kel[0,i] == 2:
                self.psct[:,start:end] = peqarr
                self.index[start:end] = 2#index?
            self.wpl[start:end] = self.eb[0,:length]
            self.ec0[start:end] = self.ec[0,:length]
            self.nc0[start:end] = self.nc[0,:length]
            self.eg1[start:end] = self.eg[0,:length]
            self.ng1[start:end] = self.ng[0,:length]
            self.eg2[start:end] = self.egs[0,:length]
            self.ng2[start:end] = self.ngs[0,:length]
            self.wklm[start:end] = self.wk[0,:length]
            self.efl[start:end] = self.efls[0,:length]
            self.rgas[start:end] = rgas1
            self.ipn[start:end] = 1 #index?
            L = 2
            self.iarry[start:end] = L
            self.esplit[start:end,:20] = self.esplits[0,self.ionmodel[0],:20]
        if self.main.efinal >= self.e[0,3]:  
            if self.natt[0] <= 1:
                n += 1
                start = n
                qarr = self.q[0,3,:]
            else:
                start = n + 1
                n += self.natt[0]
                qarr = self.qatt[0,:self.natt[0],:]
            end = n + 1
            self.idg[0] = n
            self.cf[:,start:end] = qarr * self.main.van[0] * self.main.bet
            self.fcatt += cf[:,start:end]
            self.psct[:,start:end] = 0.5
            self.angct[:,start:end] = 1.0
            self.negas[start:end] = 1
            self.rgas[start:end] = rgas1
            self.ipn[start:end] = -1
            L = 3
            self.iarry[start:end] = L
        if self.nin[0] != 0:
            start = n + 1
            n += self.nin[0]
            self.idg[0] = n
            end = n + 1
            self.negas[start:end] = 1
            self.cf[:,start:end] = self.qin[0,:self.nin[0],:] self.main.van[0] * self.main.bet
            self.angct[:,start:end] = 1.0
            self.penfra[0,start:end] = self.penfrac[0,:self.nin[0]]
            self.penfra[1,start:end] = self.penfrac[1,:self.nin[0]]*1e-6/num.sqrt(3.0)
            self.penfra[2,start:end] = self.penfrac[2,:self.nin[0]]
            for j in range(self.nin[0]):
                if self.main.lbrm == 0:
                    if self.izbrs[0,j] != 0:
                        self.cf[:,start + j] = 0.0
                if self.kin[0,j] == 1:
                    self.angct[:,start + j], self.psct[:,start + j] = acut(self.peqin[0,j,:])
                    self.index[start + j] = 1
                elif self.kin[0,j] == 2:
                    self.psct[:,start + j] = self.peqin[0,j,:]
                    self.index[start + j] = 2
                else:
                    self.psct[:,start + j] = 0.5
                L = 4
                if self.ei[0,j] < 0.0:
                    L = 5
                self.iarry[start + j] = L
                if self.penfra[0,start + j] > self.avpfrac[0,0]:
                    self.avpfrac[0,0] = self.penfra[0,start + j]
                    self.avpfrac[1,0] = self.penfra[1,start + j]
                    self.avpfrac[2,0] = self.penfra[2,start + j]
            self.cminexsc[0] = self.cminexsc[0] * self.avpfrac[0,0]
            self.rgas[start:end] = rgas1
            self.ein[start:end] = self.ei[0,:self.nin[0]]/rgas1
            self.izbr[start:end] = self.izbrs[:self.nin[0]]
            
            

    def mixgas(self,gas,i):
        if gas == 0:
            return
        if gas not in pygrad.gas_dict:
            raise PygradException('Invalid gas value')
