from Pygrad.Pygrad import Pygrad
from libc.math cimport sin, cos, acos, asin, log, sqrt, pow
from libc.string cimport memset
from Pygrad.Pygrad cimport drand48
import cython
import numpy as np
cimport numpy as np
import sys

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double random_uniform(double dummy):
        cdef double r = drand48(dummy)
        return r

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef run(Pygrad object):
    xs = np.array(150000,float)
    ys = np.array(150000,float)
    zs = np.array(150000,float)
    ts = np.array(150000,float)
    dcx = np.array(150000,float)
    dcy = np.array(150000,float)
    dcz = np.array(150000,float)
    nflgfc = np.array(150000,int)
    nflgppc = np.array(150000,int)
    nflgbrmc = np.array(150000,int)
    temp = np.array(20000,float)
    etemp = np.array(1000,float)
    if (object.FinalElectronEnergy <= 140000.0):
        estep1 = (object.FinalElectronEnergy - 16000.0) / 4000.0
    else:
        estep1 = 20.0
        estep2 = (object.FinalElectronEnergy - 92000.0) / 4000.0
    nprint = 0
    object.tmax1 = 0.0
    emax = 0.0
    RDummy = randStart
    CONST9 = object.CONST3 * 0.01
    object.numReals = 0
    object.numNull = 0
    object.totalElastic = 0
    object.totalExcited = 0
    object.totalInelastic = 0
    nmxadd = 0
    ntempflag = 0
    bp = object.EField ** 2 * object.CONST1
    f1 = object.EField * object.CONST2
    f2 = object.EField * object.CONST3
    f4 = 2.0 * np.pi
    theta1 = object.Theta
    phi1 = object.Phi
    temp = object.TotalCollisionFrequency + TotalNullCollisionFrequency
    tlim = np.min(temp)
    neovfl = 0
    j1 = 0
    for j in range(object.nDelta):
        j1 += 1
        numPrimaries = j1
        object.numExcitationsPerGas = np.zeros(6,int)
        dcz1 = cos(theta1)
        dcx1 = sin(theta1) * cos(phi1)
        dcy1 = sin(theta1) * sin(phi1)
        nflgff = 0
        nflgppp = 0
        nflgbrmm = 0
        nflghigh = 0
        est1 = object.InitialElectronEnergy
        e1 = object.InitialElectronEnergy
        x = 0.0
        y = 0.0
        z = 0.0
        k1 = 0
        kexc = 0
        nstexc = 0
        nexcnul = 0
        nclus = 0
        nelec = 0
        negion = 0
        tlast = 0.0
        st = 0.0
        tdash = 0.0
        if object.imip != 2:
            if object.imip > 2:
                ibadtot, ibad1 = casres(j11,object)
                if ibad1 == 1:
                    j1 -= 1
                    continue
            elif object.imip == 1:
                casrem(j11)
                est = ecas[0]
                est2 = est1
            x = xcas[0]
            y = ycas[0]
            z = zcas[0]
            st = tt1[0]
            ts[0] = tt1[0]
            e1 = ecas[0]
            dcz1 = drzs[0]
            dcy1 = drys[0]
            dcx1 = drxs[0]
            nflgff = nglgf[0]
            nflgppp = nglfpp[0]
            nflgbrmn = 0
            nflghigh = nglgff
            isdum = 0
            xs[1:object.ievntl] = xcas[1:object.ievntl]
            ys[1:object.ievntl] = ycas[1:object.ievntl]
            zs[1:object.ievntl] = zcas[1:object.ievntl]
            ts[1:object.ievntl] = ttl[1:object.ievntl]
            es[1:object.ievntl] = ecas[1:object.ievntl]
            dcx[1:object.ievntl] = drys[1:object.ievntl]
            dcy[1:object.ievntl] = drzs[1:object.ievntl]
            nflgfc[1:object.ievntl] = nflgf[1:object.ievntl]
            nflgppc[1:object.ievntl] = nflgpp[1:object.ievntl]
            nflghigh = max(nflghigh,np.max(nflgf[1:object.ievntl]))
            nclus = object.ievntl - 1
        R5lt = False
        while True:
            R1 = drand48(RDummy)
            T = -log(r1)/ tlim + tdash
            tdash = t
            ap = dcz1 * f2 * sqrt(e1)
            gam1 = (object.EMS + E1) / object.EMS
            bet1 = sqrt(1.0 - 1.0(gam1 * gam1))
            ap = dcz1 * object.EField * bet1 * vc * 1.0e-10
            bp1 = bp / gam1
            e = e1 + (ap + bp1 * T) * T
            if e < 0.0:
                e = 0.001
            if object.imip == 1:
                ie = int(e/object.ElectronEnergyStep) + 1
            else:
                if object.FinalElectronEnergy <= 20000.0:
                    ie = int(e/object.ElectronEnergyStep) + 1
                elif object.FinalElectronEnergy <= 140000.0:
                    if e <= 16000.0:
                        ie = int(e) + 1
                    else:
                        ie = 16000 + int((e-16000.0)/estep1)
                else:
                    if e <= 12000.0:
                        ie = int(e) + 1
                    elif e <= 92000.0:
                        ie = 12000 + int((e - 12000.0) / estep1)
                    else:
                        ie = 16000 + int((e - 92000.0) / estep2)
            ie = min(ie, 20000)
            R5 = drand48(RDummy)
            test1 = object.TotalCollisionFrequency[ie - 1] / tlim
            if R5 <= test1:
                break
            numNull += 1
            test2 = temp[ie - 1] / tlim
            if R5 < test2:
                if object.NPLAST != 0:
                    R2 = drand48(RDummy)
                    i = 0
                    while object.CFN[ie - 1, i] < R2:
                        i += 1
                    nexcnul += 1
                    object.ICOLNN[i] += 1
                    object.XSTN[nexcnul] = x
                    object.YSTN[nexcnul] = y
                    object.ZSTN[nexcnul] = z
                    object.TSTN[nexcnul] = t
                    object.IDNUL[nexcnul] = i
        T2 = T * T
        if e > emax:
            emax = e
        if t > tmax1:
            tmax1 = t
        tdash = 0.0
        numReals += 1
        CONST6 = sqrt(e1/e)
        gam2 = (object.EMS + e)/object.EMS
        gam12 = (gam1 + gam2) / 2.0
        bet2 = sqrt(1.0 - 1.0 / (gam2 * gam2))
        CONST6 = bet1 / bet2
        dcx2 = dcx1 * CONST6
        dcy2 = dcy1 * CONST6
        dcz2 = dcz1 * CONST6 + object.EField * T * 2.0e10 * object.CONST1/(VC * bet2)
        CONST7 = vc * bet1 * 1.0e-12
        A = T * CONST7
        x = x + dcx1 * a
        y = y + dcy1 * a
        z = z + dcz1 * a + T2 * F1/gam12
        st += t
        it = int(t+1.0)
        it = min(it, 300)
        object.Time[it - 1] += 1.0 
        R2 = drand48(RDummy)
        i = 0
        while object.CF[ie - 1,i] < R2:
            i += 1
        if object.izbr[i] != 0 and object.lbrm:
            nflgbrmm = 1
            ipt = object.Iarry[i]
            object.icoll[ipt] += 1
            object.icoln[i] += 1
            kngs = 1
            while kngs <= object.NumberOfGases and ipt != kngs * 5 - 1:
                kngs += 1
            iatomno = object.izbr[i]
            brems(iatomno,e,dcx2,dcy2,eout,edcx,edcy,edcz,egamma,gdcx,gdcy,gdcz)
            object.nbrem[kngs] += 1
            object.ebrtot[kngs] += 1
            e1 = eout
            dcx1 = edcx
            dcy1 = edcy
            dcz1 = edcz
            bremscasc(J11,egamma,x,y,z,st,gdcx,gdcy,gdcz,ilow)
            if ilow != 1: 
                int flag_190 = 0
                etsum = 0.0
                for kbr in range(1, ievntlb + 1):
                    nclus += 1
                    if nclus > 150000:
                        print(nclus, numReals)
                        return #TODO: this is actually a stop statement
                    es[nclus] = ecasb[kbr]
                    etsum += es[nclus]
                    xs[nclus] = xcasb[kbr]
                    ys[nclus] = ycasb[kbr]
                    zs[nclus] = zcasb[kbr]
                    ts[nclus] = ttb1[kbr]
                    dcx[nclus] = drxb[kbr]
                    dcy[nclus] = dryb[kbr]
                    dcz[nclus] = drzb[kbr]
                    nflgfc[nclus] = nflgfb[kbr] + nflhigh
                    nflgppc[nclus] = nflgppb[kbr]
                    nflgbrmc[nclus] = 2
                    if nflgfc[nclus] > nflghigh:
                        nflghigh = nflgfc[nclus]
                        flag_190 = 1
                        break
                if !flag_190:
                    s1 = 1.0 + gam2 * (rgas[i] - 1.0)
                    ei = object.ein[i]
                    if e < ei:
                        ei = e - 0.0001
                    if ipn[i] !=0:
                        if ipn[i] == -1:
                            netot += 1
                            nitot += 1
                            nelec += 1
                            negion += 1
                            ipt = object.InteractionType[i]
                            object.icoll[ipt] += 1
                            object.icoln[i] += 1
                            it = int(t+1.0)
                            it = min(it,J300)
                            time[it] += 1.0
                        else:
                            eistr = ei
                            if object.ionmodel[i] > 0:
                                ionsplit(i,e,ei,esec)
                            else:
                                R9 = drand48(0)
                                esec = wpl[i] * tan(R9 * atan((e-ei)/(2.0*wpl[i])))
                                esec = wpl[i] * (esec\wpl[i]) **0.9524
                            ei += esec
                            nclus += 1
                            nmxadd = max(nclud,nmxadd)
                            if nclus > 150000:
                                print("program stopped. nclus = {0}, nreal = {1}".format(nclus, numReals))
                                return #TODO: should be a stop
                            xs[nclus] = x
                            ys[nclus] = y
                            zs[nclus] = z
                            ts[nclus] = st
                            es[nclus] = esec
                            nflgfc[nclus] = nflgff
                            nflgppc[nclus] = nflgppp
                            nflgbrmc[nclus] = nflgbrmm
                            ntmpflg = 1
                            ncltmp = nclus
                            flag_666 = 0
                            if iecasc != 0 and legas[i] != 0:
                                kg1 = negas[i]
                                lg1 = legas[i]
                                igshel = ieshell[i]
                                cascadee(J11,kg1,lg1,x,y,z,st,esec,igshel)
                                etsum = 0.0
                                for kbr in range(1,ievente+1):
                                    nclus += 1
                                    if nclus > 150000:
                                        print("program stopped. nclus = {0}, nreal = {1}".format(nclus, numReals))
                                        return #TODO: should be a stop
                                    es[nclus] = ecase[kbr]
                                    etsum += es[nclus]
                                    xs[nclus] = xcase[kbr]
                                    ys[nclus] = ycase[kbr]
                                    zs[nclus] = zcase[kbr]
                                    ts[nclus] = tcase[kbr]
                                    dcx[nclus] = drxce[kbr]
                                    dcy[nclus] = dryce[kbr]
                                    dcz[nclus] = drzce[kbr]
                                    nflgfc[nclus] = nflgfe[kbr] + nflghigh
                                    nflgppc[nclus] = nflgppe[kbr]
                                    nflgbrmc[nclus] = nflgbrmm
                                    if nflgfc[nclus] > nflghigh:
                                        nflghigh = nflgfc[nclus]
                                        flag_666 = 1
                            if !flag_666:
                                if eistr > 30.0:
                                    ifltst = 0
                                    if wklm[i] >  0.0:
                                        R9 = drand48(0)
                                        if R9 < wklm[i]:
                                            ifltst = 1
                                    if ifltst == 0:
                                        naug = nc0[i]
                                        eavaug = ec0[i] / float(naug)
                                        for jfl in range(1, nc0[i] + 1):
                                            nclus += 1
                                            xs[nclus] = x
                                            ys[nclus] = y
                                            zs[nclus] = z
                                            ts[nclus] = st
                                            nflgfc[nclus] = nflgff
                                            nflgppc[nclus] = nflgppp
                                            nflgbrmc[nclus] = nflgbrmm
                                            es[nclus] = eavaug
                                            R3 = drand48(0)
                                            F3 = 1.0 - 2.0*R3
                                            theta0 = acos(F3)
                                            F6 = cos(theta0)
                                            F5 = sin(theta0)
                                            R4 = drand48(0)
                                            phi0 = F4 * R4
                                            F8 = sin(phi0)
                                            F9 = cos(phi0)
                                            dcx[nclus] = F9*F5
                                            dcy[nclus] = F8*F5
                                            dcz[nclus] = F6
                                    else:
                                        if ng2[i] != 0:
                                            naug = ng2[i]
                                            eavaug = eg2[i] / float(naug)
                                            for jfl in range(1, ng2[i] + 1):
                                                nclus += 1
                                                xs[nclus] = x
                                                ys[nclus] = y
                                                zs[nclus] = z
                                                nflgfc[nclus] = nflgff
                                                nflgppc[nclus] = nflgppp
                                                nflgbrmc[nclus] = nflgbrmm
                                                ts[nclus] = st
                                                es[nclus] = eavaug
                                                R3 = drand48(0)
                                                F3 = 1.0-2.0* R3
                                                theta0 = acos(F3)
                                                F6 = cos(theta0)
                                                F5 = sin(theta0)
                                                R4 = drand48(0)
                                                phi0 = F4*R4
                                                F8 = sin(phi0)
                                                F9 = cos(phi0)
                                                dcx[nclus] = F9*F5
                                                dcy[nclus] = F8*F5
                                                dcz[ncus] = F6
                                        if ng1[i] != 0:
                                            naug = ng1[i]
                                            eavaug = eg1[i] / float(naug)
                                            R9 = drand48(0)
                                            dfl = -log(R9) * dstfl[i]
                                            for jfl in range(1, ng1[i] + 1):
                                                nclus += 1
                                                R3 = drand48(0)
                                                thefl = acos(1.0-2.0*R3)
                                                R4 = drand48(0)
                                                phifl = F4*R4
                                                xs[nclus] = x + dfl * sin(thefl) * cos(phifl)
                                                ys[nclus] = y+dfl * sin(thefl) * sin(phifl)
                                                zs[nclus] = z + dfl * cos(thefl)
                                                nflgfc[nclus] = nflghigh + 1
                                                nflgppc[nclus] = nflgppp
                                                nflgbrmc[nclus] = nflgbrmm
                                                ts[nclus] = st
                                                es[nclus] = eavaug
                                                R3 = drand48(0)
                                                F3 = 1.0 - 2.0 * R3
                                                theta0 = acos(F3)
                                                F6 = cos(theta0)
                                                F5 = sin(theta0)
                                                R4 = drand48(0)
                                                phi- = F4 * R4
                                                F8 = sin(phi0)
                                                F9 = cos(phi0)
                                                dcx[nclus] = F9*F5
                                                dcy[nclus] = F8*F5
                                                dcz[nclus] = F6
                                                nflghigh = nflgfc[nclus]
                            ipt = object.InteractionType[i]
                            icoll[ipt] += 1
                            icoln[i] += 1
                            if ipen == 0 or object.NumberOfGases == 1:
                            xs[nclus] = x
                            ys[nclus] = y
                            zs[nclus] = z


      J1=0
C START OF PRIMARY EVENT LOOP 
      DO 210 J11=1,NDELTA
       IF(ILOW.EQ.1) GO TO 190
      IF(IPN(I).EQ.0) GO TO 666
      IF(IPN(I).EQ.-1) THEN
       GO TO 335
      ICOLN(I)=ICOLN(I)+1 
C IF EXCITATION THEN ADD PROBABILITY ,PENFRA(1,I),OF TRANSFER TO GIVE   
C IONISATION OF THE OTHER GASES IN MIXTURE
      IF(IPEN.EQ.0.OR.NGAS.EQ.1) GO TO 5 
      IF(PENFRA(1,I).NE.0.0) THEN      
       RAN=drand48(RDUM)
       IF(RAN.GT.PENFRA(1,I)) GO TO 5
       NCLUS=NCLUS+1  
C ENTER HERE POSSIBLE DELOCALISATION LENGTH FOR PENNING TRANSFER
       IF(PENFRA(2,I).EQ.0.0) THEN
        XS(NCLUS)=X
        YS(NCLUS)=Y      
        ZS(NCLUS)=Z             
        NFLGFC(NCLUS)=NFLGFF
        NFLGPPC(NCLUS)=NFLGPPP
        NFLGBRMC(NCLUS)=NFLGBRMM
        GO TO 667
       ENDIF
       ASIGN=1.0
       RAN=drand48(RDUM)
       RAN1=drand48(RDUM)
       IF(RAN1.LT.0.5) ASIGN=-ASIGN
       XS(NCLUS)=X-DLOG(RAN)*PENFRA(2,I)*ASIGN
       RAN=drand48(RDUM)
       RAN1=drand48(RDUM)
       IF(RAN1.LT.0.5) ASIGN=-ASIGN
       YS(NCLUS)=Y-DLOG(RAN)*PENFRA(2,I)*ASIGN
       RAN=drand48(RDUM)
       RAN1=drand48(RDUM)
       IF(RAN1.LT.0.5) ASIGN=-ASIGN
       ZS(NCLUS)=Z-DLOG(RAN)*PENFRA(2,I)*ASIGN
  667  RAN=drand48(RDUM)
       TS(NCLUS)=ST-DLOG(RAN)*PENFRA(3,I)
C ASSIGN EXCESS ENERGY OF 1EV TO PENNING CREATED ELECTRON
       ES(NCLUS)=1.0
       DCX(NCLUS)=DCX1
       DCY(NCLUS)=DCY1
       DCZ(NCLUS)=DCZ1
       GO TO 6
      ENDIF 
C      GO TO 6 
C CALCULATE SUM OF EXCITATION PER CLUSTER AND STORE EXCITATION X Y Z T
   5  IF(IPN(I).EQ.0) THEN
       IF((RGAS(I)*EIN(I)).GT.4.0) THEN
        KEXC=KEXC+1
        IF(KEXC.GT.150000) THEN
         WRITE(6,548) KEXC
 548     FORMAT(2X,' PROGRAM STOPPED . KEXC=',I7)
         STOP
        ENDIF
C FIND GAS IN WHICH EXCITATION OCCURED AND INCREMENT COUNTER
        IF(I.LE.IDG1) THEN 
         NGEXC1=NGEXC1+1
        ELSE IF(I.LE.IDG2) THEN
         NGEXC2=NGEXC2+1
        ELSE IF(I.LE.IDG3) THEN
         NGEXC3=NGEXC3+1
        ELSE IF(I.LE.IDG4) THEN
         NGEXC4=NGEXC4+1
        ELSE IF(I.LE.IDG5) THEN
         NGEXC5=NGEXC5+1
        ELSE IF(I.LE.IDG6) THEN
         NGEXC6=NGEXC6+1
        ELSE
         WRITE(6,9911) 
 9911    FORMAT(' PROGRAM STOPPED BAD GAS ID IN MONTE')
         STOP
        ENDIF
        NEXCTOT=NEXCTOT+1
        NSTEXC=NSTEXC+1
        XSTEXC(KEXC)=X
        YSTEXC(KEXC)=Y
        ZSTEXC(KEXC)=Z
        TSTEXC(KEXC)=ST
       ENDIF
      ENDIF 
   6  S2=(S1*S1)/(S1-1.0D0) 
C   ANISOTROPIC SCATTERING
      R3=drand48(RDUM)
      IF(INDEX(I).EQ.1) THEN
       R31=drand48(RDUM)
       F3=1.0D0-R3*ANGCT(IE,I)          
       IF(R31.GT.PSCT(IE,I))  F3=-F3
      ELSE IF(INDEX(I).EQ.2) THEN
       EPSI=PSCT(IE,I)
       F3=1.0D0-(2.0D0*R3*(1.0D0-EPSI)/(1.0D0+EPSI*(1.0D0-2.0D0*R3))) 
      ELSE 
C ISOTROPIC SCATTERING                                             
       F3=1.0D0-2.0D0*R3
      ENDIF
      THETA0=DACOS(F3)
      R4=drand48(RDUM)
      PHI0=F4*R4                                                        
      F8=DSIN(PHI0)                                                     
      F9=DCOS(PHI0)                                                     
      IF(E.LT.EI) EI=0.0D0                                              
      ARG1=1.0D0-S1*EI/E                                                
      ARG1=DMAX1(ARG1,SMALL)                                            
      D=1.0D0-F3*DSQRT(ARG1)                                            
      E1=E*(1.0D0-EI/(S1*E)-2.0D0*D/S2) 
      E1=DMAX1(E1,SMALL)                                                
      Q=DSQRT((E/E1)*ARG1)/S1                                           
      Q=DMIN1(Q,1.0D0)                                                  
      THETA=DASIN(Q*DSIN(THETA0))                                       
      F6=DCOS(THETA)                                                    
      U=(S1-1.0D0)*(S1-1.0D0)/ARG1                                      
      CSQD=F3*F3                                                        
      IF(F3.LT.0.0D0.AND.CSQD.GT.U) F6=-1.0D0*F6                        
      F5=DSIN(THETA)                                                    
      DCZ2=DMIN1(DCZ2,1.0D0)                                            
      ARGZ=DSQRT(DCX2*DCX2+DCY2*DCY2)
      IF(ARGZ.EQ.0.0D0) THEN
       DCZ1=F6         
       DCX1=F9*F5                             
       DCY1=F8*F5 
       IF(NTMPFLG.EQ.1) THEN
C USE FREE KINEMATICS FOR IONISATION SECONDARY ANGLES
        F5S=F5*DSQRT(E1/ES(NCLTMP))
        IF(F5S.GT.1.0) F5S=1.0
        THSEC=DASIN(F5S)
        F5S=DSIN(THSEC)
        F6S=DCOS(THSEC)
        IF(F6.LT.0.0) F6S=-F6S
        PHIS=PHI0+API   
        IF(PHIS.GT.F4) PHIS=PHI0-F4
        F8S=DSIN(PHIS)
        F9S=DCOS(PHIS)
        DCZ(NCLTMP)=F6S
        DCX(NCLTMP)=F9S*F5S
        DCY(NCLTMP)=F8S*F5S
        NTMPFLG=0
       ENDIF
       GO TO 190
      ENDIF                                            
      DCZ1=DCZ2*F6+ARGZ*F5*F8                                           
      DCY1=DCY2*F6+(F5/ARGZ)*(DCX2*F9-DCY2*DCZ2*F8)                     
      DCX1=DCX2*F6-(F5/ARGZ)*(DCY2*F9+DCX2*DCZ2*F8)
      IF(NTMPFLG.EQ.1) THEN
C USE FREE KINEMATICS FOR IONISATION SECONDARY ANGLES
       F5S=F5*DSQRT(E1/ES(NCLTMP))
       IF(F5S.GT.1.0) F5S=1.0            
       THSEC=DASIN(F5S)
       F5S=DSIN(THSEC)
       F6S=DCOS(THSEC)
       IF(F6.LT.0.0) F6S=-F6S
       PHIS=PHI0+API   
       IF(PHIS.GT.F4) PHIS=PHI0-F4
       F8S=DSIN(PHIS)
       F9S=DCOS(PHIS)
       DCZ(NCLTMP)=DCZ2*F6S+ARGZ*F5S*F8S                               
       DCY(NCLTMP)=DCY2*F6S+(F5S/ARGZ)*(DCX2*F9S-DCY2*DCZ2*F8S)        
       DCX(NCLTMP)=DCX2*F6S-(F5S/ARGZ)*(DCY2*F9S+DCX2*DCZ2*F8S)
       NTMPFLG=0
      ENDIF 
  190 CONTINUE 
C TEST IF ELECTRON IS THERMALISED
      IF(E1.GT.ETHRM) GO TO 1  
C STORE POSITION AND TIME OF ELECTRON AND COLLISION HISTORY
  191 CONTINUE
      K1=K1+1
      XST(K1)=X
      YST(K1)=Y
      ZST(K1)=Z
      TST(K1)=ST
      NFGF(K1)=NFLGFF
      NFGPP(K1)=NFLGPPP
      NFGBR(K1)=NFLGBRMM
      TTIME(K1)=ST-TLAST
      NELEC=NELEC+1
      NETOT=NETOT+1
 335  IF(K1.EQ.150000) GO TO 889
C CATCH SINGLE ELECTRON CLUSTER THAT WAS ATTACHED.
c     IF(NELEC.EQ.1.AND.NCLUS.EQ.0) GO TO 210 
C
        IF(NELEC.EQ.(NCLUS+1)) THEN
C       WRITE(6,884) NELEC,NCLUS,NEGION,J11
C 884 FORMAT(' NELEC=',I6,' NCLUS=',I6,' NEGION=',I3,' J11=',I6)
C LAST ELECTRON IN CLUSTER DO STATISTICS OVER PRIMARY CLUSTER
        CALL STATS(J11,J1)
        GO TO 210
       ENDIF
      IF(NELEC.LT.(NCLUS+1)) THEN
C GET NEW IONISATION ELECTRON FROM STORE
       X=XS(NELEC)
       Y=YS(NELEC)
       Z=ZS(NELEC)
       ST=TS(NELEC)
       NFLGFF=NFLGFC(NELEC)
       NFLGPPP=NFLGPPC(NELEC)
       NFLGBRMM=NFLGBRMC(NELEC)
       TLAST=TS(NELEC)
       E1=ES(NELEC)
       DCX1=DCX(NELEC)
       DCY1=DCY(NELEC)
       DCZ1=DCZ(NELEC)
C      IF(E1.LT.ETHRM) GO TO 191                       
       GO TO 1   
      ENDIF
C  MAIN LOOP END    
  210 CONTINUE
C RESET NUMBER OF EVENTS FOR BAD EVENTS
      IF(IMIP.GT.2) NDELTA=NDELTA-IBADTOT
C
      WRITE(6,887) EMAX,NEOVFL
  887 FORMAT(' EMAX=',D12.7,' NEOVFL=',I5)  
      IF(EMAX.GT.EFINAL) THEN
      WRITE(6,989) EFINAL,EMAX
  989 FORMAT('INCREASE ENERGY LIMIT FROM',D12.6,' EV TO AT LEAST',D12.6,
     /' EV.')
      STOP
      ENDIF                                         
      RETURN 
  889 NLEFT=NCLUS-NELEC
      WRITE(6,992) NPRIME,NLEFT,NCLUS
  992 FORMAT(3(/),' WARNING STOPPED AFTER NPRIME=',I6,' LAST PRIMARY HAS
     /AT LEAST ',I6,' SECONDARIES LEFT TO TRACK OUT OF ',I6,' ELECTRONS 
     /ALREADY IN CLUSTER') 
      STOP                      
      RETURN                                                            
      END
