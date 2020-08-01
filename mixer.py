import pygrad, os, h5py, numpy as np,utils
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
        self.qrel = np.zeros(20000,dtype = float)
        self.qinel = np.zeros(20000,dtype = float)
        self.nin = np.zeros(6,dtype = int)
        self.lion = np.zeros(6,dtype = int)
        self.lin = np.zeros((6,250),dtype = int)
        self.alion = np.zeros(6,dtype = float)
        self.alin = np.zeros((6,250),dtype = float)
        self.iplast = 0
        self.nplast = 0
        self.niso = 0
        self.fcion = np.zeros(20000,dtype = float)
        self.fcatt = np.zeros(20000,dtype = float)
        self.iecasc = 0
        self.doubles = np.zeros((6,20000), dtype = float)
        self.cminixsc = np.zeros(6,dtype = float)
        self.cminexsc = np.zeros(6,dtype = float)
        self.ecloss = np.zeros(6,dtype = float)
        self.wpln = np.zeros(6,dtype = float)
        self.avpfrac = np.zeros((3,6),dtype = float)
        self.ngexc = np.zeros(6,dtype = int)
        self.idg = np.zeros(6,dtype = int)

        #Only used in mixer program
        self.q = np.zeros((6,6,20000),dtype = float)
        self.e = np.zeros((6,6),dtype = float)
        self.ei = np.zeros((6,250),dtype = float)
        self.eion = np.zeros(6,dtype = float)
        self.peqel = np.zeros((6,6,20000),dtype = float)
        self.peqin = np.zeros((6,250,20000),dtype = float)
        self.penfrac = np.zeros((6,3,250),dtype = float)  
        self.kin = np.zeros((6,250),dtype = int)
        self.kel = np.zeros((6,6),dtype = int)
        self.eions = np.zeros((6,30),dtype = float)
        self.qions = np.zeros((6,30,20000),dtype = float)
        self.peqions = np.zeros((6,30,20000),dtype = float)
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
        nioncopy = [i if i >= 1 else 1 for i in self.nion]
        nattcopy = [i if i >= 1 else 1 for i in self.natt]
        nincopy = [i if i >= 1 else 0 for i in self.nin]
        total = np.sum(nioncopy) + np.sum(nattcopy) + np.sum(nincopy)

        self.cf = np.zeros((20000,total),dtype = float)
        self.ein = np.zeros(total,dtype = float)
        self.psct = np.zeros((20000,total),dtype = float)
        self.angct = np.zeros((20000,total),dtype = float)
        self.index = np.zeros(total, dtype = int)
        self.iarry = np.zeros(total,dtype = int)
        self.rgas = np.zeros(total,dtype = float)
        self.ipn = np.zeros(total,dtype = int)
        self.wpl = np.zeros(total,dtype = float)
        self.izbr = np.zeros(total,dtype = int)
        self.penfra = np.zeros((3,total),dtype = float)
        self.negas = np.zeros(total,dtype = int)
        self.legas = np.zeros(total,dtype = int)
        self.ieshell = np.zeros(total,dtype = int)
        self.nc0 = np.zeros(total,dtype = int)
        self.ec0 = np.zeros(total,dtype = float)
        self.ng1 = np.zeros(total,dtype = int)
        self.eg1 = np.zeros(total,dtype = float)
        self.ng2 = np.zeros(total,dtype = int)
        self.eg2 = np.zeros(total,dtype = float)
        self.wklm = np.zeros(total,dtype = float)
        self.efl = np.zeros(total,dtype = float)
        self.esplit = np.zeros((total,20),dtype = float)
        self.ionmodel = np.zeros(total,dtype = int)
        n = 0
        acut = np.vectorize(utils.angcut)
        for gas in range(self.main.ngas):
            self.idg[gas] = n
            self.negas[n] = gas + 1 #index?
            self.cf[:,n] = self.q[gas,1,:] * self.main.vans[gas] * self.main.bet
            self.psct[:,n] = 0.5
            self.angct[:,n] = 1.0
            if self.kel[gas,1] == 1:
                self.angct[:,n],self.psct[:,n] = acut(peqel[gas,1])
                self.index[n] = 1 #index?
            elif self.kel[gas,1] == 2:
                self.psct[:,n] = peqel[gas,1]
                self.index[n] = 2 #index?
            rgas1 = 1.0 + self.e[gas,1] / 2.0
            self.rgas[n] = rgas1
            L = 1 + 5*gas
            self.iarry[n] = L #index?
            self.cminexsc[n] = self.e[gas,3] * self.main.ans[gas]
            self.cminixsc[n] = self.e[gas,4] * self.main.ans[gas]
            self.ecloss[n] = self.e[gas,2]
            self.wpln[n] = self.e[gas,5]
            if self.main.efinal >= self.e[gas,2]: 
                if self.nion[gas] > 1:
                    start = n + 1
                    n += self.nion[gas]
                    end = n + 1
                    self.idg[gas] = n
                    i = 2
                    self.cf[:,start:end] = self.qions[gas,:self.nion[gas]] * self.main.vans[gas] * self.main.bet
                    self.fcion += np.sum(self.cf[:,start:end],0)
                    self.legas[start:end] = self.legasn[gas,:self.nion[gas]]
                    self.ieshell[start:end] = self.ieshel[gas,:self.nion[gas]]
                    peqarr = self.peqions[gas,:self.nion[gas]]
                    self.ein[start:end] = self.eion[gas,:self.nion[gas]]/rgas1
                else:
                    n += 1
                    start = n
                    end = n + 1
                    self.idg[gas] = n
                    if self.main.icount == 1:
                        i = 4
                        self.doubles[gas] = self.q[gas,2] / self.q[gas,4] - 1.0
                    else:
                        i = 2
                    self.cf[:,n] = self.q[gas,i] * self.main.vans[gas] * self.main.bet
                    self.fcion += self.cf[:,n]
                    peqarr = self.peqel[gas,i]
                    self.wpl[n] = self.eb[gas,0]
                    self.ein[n] = self.e[gas,2]/rgas1
                length = end - start
                self.negas[start:end] = gas + 1 #index?
                self.psct[:,start:end] = 0.5
                self.angct[:,start:end] = 1.0
                if self.kel[gas,i] == 1:
                    self.angct[:,start:end],self.psct[:,start:end] = acut(peqarr)
                    self.index[start:end] = 1#index?
                elif self.kel[gas,i] == 2:
                    self.psct[:,start:end] = peqarr
                    self.index[start:end] = 2#index?
                self.wpl[start:end] = self.eb[gas,:length]
                self.ec0[start:end] = self.ec[gas,:length]
                self.nc0[start:end] = self.nc[gas,:length]
                self.eg1[start:end] = self.eg[gas,:length]
                self.ng1[start:end] = self.ng[gas,:length]
                self.eg2[start:end] = self.egs[gas,:length]
                self.ng2[start:end] = self.ngs[gas,:length]
                self.wklm[start:end] = self.wk[gas,:length]
                self.efl[start:end] = self.efls[gas,:length]
                self.rgas[start:end] = rgas1
                self.ipn[start:end] = 1 #index?
                L = 2 + 5 * gas
                self.iarry[start:end] = L
                self.esplit[start:end,:20] = self.esplits[gas,self.ionmodel[gas],:20]
            if self.main.efinal >= self.e[gas,3]:  
                if self.natt[gas] <= 1:
                    n += 1
                    start = n
                    qarr = self.q[gas,3]
                    self.cf[:,n] = qarr * self.main.vans[gas] * self.main.bet
                else:
                    start = n + 1
                    n += self.natt[gas]
                    qarr = self.qatts[gas,:self.natt[gas]]
                    self.cf[:,start:end] = qarr * self.main.vans[gas] * self.main.bet
                end = n + 1
                self.idg[gas] = n
                self.fcatt = self.fcatt +  self.cf[:,start:end]
                self.psct[:,start:end] = 0.5
                self.angct[:,start:end] = 1.0
                self.negas[start:end] = 1 + gas #index?
                self.rgas[start:end] = rgas1
                self.ipn[start:end] = -1
                L = 3 + 5 * gas
                self.iarry[start:end] = L
            if self.nin[gas] != 0:
                start = n + 1
                n += self.nin[gas]
                self.idg[gas] = n
                end = n + 1
                self.negas[start:end] = 1
                self.cf[:,start:end] = self.qin[gas,:self.nin[gas]] * self.main.vans[gas] * self.main.bet
                self.angct[:,start:end] = 1.0
                self.penfra[0,start:end] = self.penfrac[0,:self.nin[gas]]
                self.penfra[1,start:end] = self.penfrac[1,:self.nin[gas]]*1e-6/num.sqrt(3.0)
                self.penfra[2,start:end] = self.penfrac[2,:self.nin[gas]]
                for j in range(self.nin[gas]):
                    if self.main.lbrm == 0:
                        if self.izbrs[gas,j] != 0:
                            self.cf[:,start + j] = 0.0
                    if self.kin[gas,j] == 1:
                        self.angct[:,start + j], self.psct[:,start + j] = acut(self.peqin[gas,j])
                        self.index[start + j] = 1
                    elif self.kin[gas,j] == 2:
                        self.psct[:,start + j] = self.peqin[gas,j]
                        self.index[start + j] = 2
                    else:
                        self.psct[:,start + j] = 0.5
                    L = 4 + 5 * gas
                    if self.ei[0,j] < 0.0:
                        L = 5 + 5 * gas
                    self.iarry[start + j] = L
                    if self.penfra[0,start + j] > self.avpfrac[0,gas]:
                        self.avpfrac[0,gas] = self.penfra[0,start + j]
                        self.avpfrac[1,gas] = self.penfra[1,start + j]
                        self.avpfrac[2,gas] = self.penfra[2,start + j]
                self.cminexsc[gas] *= self.avpfrac[0,gas]
                self.rgas[start:end] = rgas1
                self.ein[start:end] = self.ei[gas,:self.nin[gas]]/rgas1
                self.izbr[start:end] = self.izbrs[gas,:self.nin[gas]]
        self.iplast = n
        if np.any(self.cf < 0):
            print(' Warning: negative collision frequency')
        self.tcf = np.sum(self.cf,1)
        safediv = lambda x,y: x / y if y != 0 else 0
        safediv = np.vectorize(safediv)
        self.cf = safediv(self.cf.T,self.tcf)
        self.cf = np.cumsum(self.cf,axis=1)
        self.cf[:,-1] = 1.0
        eroot = self.main.eroot * 1e-10
        self.fcatt *= eroot
        self.fcion *= eroot
        self.tcf *= eroot
        n = 0
        self.nplast = max(0,np.sum(self.nul))
        if np.sum(self.nul) != 0:
            self.cfn = np.zeros((20000,self.nplast),dtype = float)
            self.sclenul = np.zeros(self.nplast,dtype = float)
            for i in range(6):
                j = self.nul[i]
                if j > 0:
                    start = n + 1
                    n += j
                    end = n + 1
                    self.sclenul[start:end] = self.scln[i,:j]
                    self.cfn[:,start:end] = self.qnul[i,:j] * self.main.vans[i] * self.sclenul[start:end] * self.main.bet
        else:
            self.cfn = np.zeros((20000,1),dtype = float)
            self.sclenul = np.zeros(1,dtype = float)
        self.tcfn = np.sum(self.cfn,1)
        if np.any(self.cfn < 0):
            print(' Warning: negative null collision frequency')
        self.cfn = safediv(self.cfn.T,self.tcfn)
        self.tcfn *= 1e-10
        self.cfn = np.cumsum(self.cfn,axis=1)
        self.cfn[:,-1] = 1.0
        print('\n'.join(["Index={}, J={}".format(x[0],x[1]) for x in zip(self.index,range(self.iplast))]))
        kelsum = np.sum(self.kel) + np.sum(self.kin)
        if kelsum > 0:
            self.niso = 1
            print('Anistropic scattering detected NIso = 1')
        bp = self.main.efield ** 2 * pygrad.CONST1
        f2 = self.main.efield * pygrad.CONST3
        elow = self.main.tmax * (self.main.tmax * bp - f2 * np.sqrt(0.5 * self.main.efinal))/self.main.estep-1.0
        elow = min(elow,pygrad.SMALL)
        ehi = self.main.tmax * (self.main.tmax * bp + f2 * np.sqrt(0.5 * self.main.efinal)) / self.main.estep + 1.0
        ehi = min(ehi,20000.0)
        jone = 0
        jlarge = 20000
        for i in range(10):
            jlow = 20000-2000*(10-i) + int(elow)
            jhi = 20000-2000*(9-i) + int(ehi)
            jlow = max(jlow,jone)
            jhi = min(jhi,jlarge)
            self.main.tcfmax[i] = np.max(self.tcf[jlow:jhi])
        self.tcfmax1 = np.max(self.tcf)
        self.qtot,self.qel = np.split(np.sum(self.main.ans * self.q[:,:2].T,0),2)
        self.qtot = self.qtot.T
        self.qel = self.qel.T
        for i in range(6):
            self.qion[i] = self.q[i,2] * self.main.ans[i]
            if self.nion[i] > 1:
                self.qion[i] = self.qions[i,self.nion[i]] * self.main.ans[i]
        self.qatt = (self.q[:,4].T * self.main.ans).T
        self.qsatt = np.sum(self.qatt,0)
        self.qsum = np.sum(self.qion,0) + self.qsatt
        self.qrel = qsum - 2 * self.qsatt
        for i in range(6):
            if i > 0:
                self.qsum += np.sum(self.qin[i,:self.nin[i]],1) * self.main.ans[i]

    def mixgas(self,gas,i):
        if gas == 0:
            return
        if gas not in pygrad.gas_dict:
            raise PygradException('Invalid gas value')
