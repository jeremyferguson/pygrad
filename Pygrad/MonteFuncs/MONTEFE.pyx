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
cpdef run(Pygrad Object):
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
    if (Object.FinalElectronEnergy <= 140000.0):
        estep1 = (Object.FinalElectronEnergy - 16000.0) / 4000.0
    else:
        estep1 = 20.0
        estep2 = (Object.FinalElectronEnergy - 92000.0) / 4000.0
    nprint = 0
    Object.tmax1 = 0.0
    emax = 0.0
    RDummy = randStart
    CONST9 = Object.CONST3 * 0.01
    Object.numReals = 0
    Object.numNull = 0
    Object.totalElastic = 0
    Object.totalExcited = 0
    Object.totalInelastic = 0
    nmxadd = 0
    ntempflag = 0
    bp = Object.EField ** 2 * Object.CONST1
    f1 = Object.EField * Object.CONST2
    f2 = Object.EField * Object.CONST3
    f4 = 2.0 * np.pi
    theta1 = Object.Theta
    phi1 = Object.Phi
    temp = Object.TotalCollisionFrequency + TotalNullCollisionFrequency
    tlim = np.min(temp)
    neovfl = 0
    j1 = 0
    for j in range(object.nDelta):
        j1 += 1
        numPrimaries = j1
        Object.numExcitationsPerGas = np.zeros(6,int)
        dcz1 = cos(theta1)
        dcx1 = sin(theta1) * cos(phi1)
        dcy1 = sin(theta1) * sin(phi1)
        nflgff = 0
        nflgppp = 0
        nflgbrmm = 0
        nflghigh = 0
        est1 = Object.InitialElectronEnergy
        e1 = Object.InitialElectronEnergy
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
        if Object.imip != 2:
            if Object.imip > 2:
                ibadtot, ibad1 = casres(j11,Object)
                if ibad1 == 1:
                    j1 -= 1
                    continue
            elif Object.imip == 1:
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
            xs[1:Object.ievntl] = xcas[1:Object.ievntl]
            ys[1:Object.ievntl] = ycas[1:Object.ievntl]
            zs[1:Object.ievntl] = zcas[1:Object.ievntl]
            ts[1:Object.ievntl] = ttl[1:Object.ievntl]
            es[1:Object.ievntl] = ecas[1:Object.ievntl]
            dcx[1:Object.ievntl] = drys[1:Object.ievntl]
            dcy[1:Object.ievntl] = drzs[1:Object.ievntl]
            nflgfc[1:Object.ievntl] = nflgf[1:Object.ievntl]
            nflgppc[1:Object.ievntl] = nflgpp[1:Object.ievntl]
            nflghigh = max(nflghigh,np.max(nflgf[1:Object.ievntl]))
            nclus = Object.ievntl - 1
        R5lt = False
        while True:
            R1 = drand48(RDummy)
            T = -log(r1)/ tlim + tdash
            tdash = t
            ap = dcz1 * f2 * sqrt(e1)
            gam1 = (Object.EMS + E1) / Object.EMS
            bet1 = sqrt(1.0 - 1.0(gam1 * gam1))
            ap = dcz1 * Object.EField * bet1 * vc * 1.0e-10
            bp1 = bp / gam1
            e = e1 + (ap + bp1 * T) * T
            if e < 0.0:
                e = 0.001
            if Object.imip == 1:
                ie = int(e/Object.ElectronEnergyStep) + 1
            else:
                if Object.FinalElectronEnergy <= 20000.0:
                    ie = int(e/Object.ElectronEnergyStep) + 1
                elif Object.FinalElectronEnergy <= 140000.0:
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
            test1 = Object.TotalCollisionFrequency[ie - 1] / tlim
            if R5 <= test1:
                break
            numNull += 1
            test2 = temp[ie - 1] / tlim
            if R5 < test2:
                if Object.NPLAST != 0:
                    R2 = drand48(RDummy)
                    i = 0
                    while Object.CFN[ie - 1, i] < R2:
                        i += 1
                    nexcnul += 1
                    Object.ICOLNN[i] += 1
                    Object.XSTN[nexcnul] = x
                    Object.YSTN[nexcnul] = y
                    Object.ZSTN[nexcnul] = z
                    Object.TSTN[nexcnul] = t
                    Object.IDNUL[nexcnul] = i
        T2 = T * T
        if e > emax:
            emax = e
        if t > tmax1:
            tmax1 = t
        tdash = 0.0
        numReals += 1
        CONST6 = sqrt(e1/e)
        gam2 = (Object.EMS + e)/Object.EMS
        gam12 = (gam1 + gam2) / 2.0
        bet2 = sqrt(1.0 - 1.0 / (gam2 * gam2))
        CONST6 = bet1 / bet2
        dcx2 = dcx1 * CONST6
        dcy2 = dcy1 * CONST6
        dcz2 = dcz1 * CONST6 + Object.EField * T * 2.0e10 * Object.CONST1/(VC * bet2)
        CONST7 = vc * bet1 * 1.0e-12
        A = T * CONST7
        x = x + dcx1 * a
        y = y + dcy1 * a
        z = z + dcz1 * a + T2 * F1/gam12
        st += t
        it = int(t+1.0)
        it = min(it, 300)
        Object.Time[it - 1] += 1.0 
        R2 = drand48(RDummy)
        i = 0
        while Object.CF[ie - 1,i] < R2:
            i += 1
        if Object.izbr[i] != 0 and Object.lbrm:
            nflgbrmm = 1
            ipt = Object.Iarry[i]
            Object.Icoll[ipt] += 1
            Object.
      J1=0
C START OF PRIMARY EVENT LOOP 
      DO 210 J11=1,NDELTA
      R2=drand48(RDUM)
      IF(IZBR(I).NE.0.AND.LBRM.EQ.1) THEN
       NFLGBRMM=1
       IPT=IARRY(I)
       ICOLL(IPT)=ICOLL(IPT)+1
       ICOLN(I)=ICOLN(I)+1
       DO 141 KNGS=1,NGAS
       IF(IPT.EQ.(KNGS*5)-1) GO TO 142
  141  CONTINUE
  142  IATOMNO=IZBR(I) 
       CALL BREMS(IATOMNO,E,DCX2,DCY2,DCZ2,EOUT,EDCX,EDCY,EDCZ,
     /EGAMMA,GDCX,GDCY,GDCZ)
       NBREM(KNGS)=NBREM(KNGS)+1
       EBRTOT(KNGS)=EBRTOT(KNGS)+EGAMMA
C GET  NEW DRCOS DRCOSY DRCOSX AND ENERGY OF ELECTRON
       E1=EOUT
       DCX1=EDCX
       DCY1=EDCY
       DCZ1=EDCZ
C RUN BREMSSTRAHLUNG GAMMA THROUGH CASCADE THEN STORE CONVERTED 
C ELECTRONS IN COMMON/CASRSB/
C
       CALL BREMSCASC(J11,EGAMMA,X,Y,Z,ST,GDCX,GDCY,GDCZ,ILOW)
C BREMSSTRAHLUNG ENERGY TOO LOW TO IONISE
       IF(ILOW.EQ.1) GO TO 190
C    
C  STORE BREMSSTRAHLUNG DATA IN CLUSTER STORE
C 
       ETSUM=0.0     
       DO 890 KBR=1,IEVNTLB
       NCLUS=NCLUS+1
       IF(NCLUS.GT.150000) THEN 
        WRITE(6,546) NCLUS,NREAL
        STOP
       ENDIF    
       ES(NCLUS)=ECASB(KBR)
       ETSUM=ETSUM+ES(NCLUS)
       XS(NCLUS)=XCASB(KBR)
       YS(NCLUS)=YCASB(KBR)
       ZS(NCLUS)=ZCASB(KBR)
       TS(NCLUS)=TTB1(KBR)
       DCX(NCLUS)=DRXB(KBR)
       DCY(NCLUS)=DRYB(KBR)
       DCZ(NCLUS)=DRZB(KBR)
       NFLGFC(NCLUS)=NFLGFB(KBR)+NFLGHIGH
       NFLGPPC(NCLUS)=NFLGPPB(KBR)
       NFLGBRMC(NCLUS)=2
  890  CONTINUE 
       IF(NFLGFC(NCLUS).GT.NFLGHIGH) NFLGHIGH=NFLGFC(NCLUS)
       GO TO 190
      ENDIF                
  891 CONTINUE  
C*****************************************************************
C     S1=RGAS(I)   
      S1=1.0D0+GAM2*(RGAS(I)-1.0D0)                                    
      EI=EIN(I)
C     WRITE(6,8890) EIN(I),I
C8890 FORMAT(' EIN=',D12.4,' I=',I3)
      IF(E.LT.EI) THEN
      EI=E-0.0001D0
      ENDIF                                                          
      IF(IPN(I).EQ.0) GO TO 666
C ATTACHMENT       
      IF(IPN(I).EQ.-1) THEN
       NETOT=NETOT+1
       NITOT=NITOT+1
       NELEC=NELEC+1
       NEGION=NEGION+1
       IPT=IARRY(I)
       ICOLL(IPT)=ICOLL(IPT)+1
       ICOLN(I)=ICOLN(I)+1 
       IT=DINT(T+1.0D0)
       IT=DMIN0(IT,J300)
       TIME(IT)=TIME(IT)+1.0D0
       GO TO 335
      ENDIF
      EISTR=EI
      IF(IONMODEL(I).GT.0) THEN 
C CALCULATE SECONDARY ENERGY,ESEC,IN IONISATION COLLISION USING
C FIVE DIFFERENT MODELS
       CALL IONSPLIT(I,E,EI,ESEC)
       GO TO 544
      ENDIF
      R9=drand48(RDUM)
C    USE OPAL PETERSON AND BEATY SPLITTING FACTOR.
      ESEC=WPL(I)*TAN(R9*ATAN((E-EI)/(2.0D0*WPL(I))))
      ESEC=WPL(I)*(ESEC/WPL(I))**0.9524
  544 CONTINUE
      EI=ESEC+EI 
C STORE POSITION ,ENERGY, DIRECTION COSINES AND TIME OF GENERATION
C OF SECONDARY IONISATION ELECTRONS
      NCLUS=NCLUS+1
      NMXADD=MAX(NCLUS,NMXADD)
      IF(NCLUS.GT.150000) THEN 
      WRITE(6,546) NCLUS,NREAL
 546  FORMAT(2X,' PROGRAM STOPPED . NCLUS=',I7,' NREAL=',I10)
      STOP
      ENDIF     
      XS(NCLUS)=X       
      YS(NCLUS)=Y
      ZS(NCLUS)=Z
      TS(NCLUS)=ST
      ES(NCLUS)=ESEC   
      NFLGFC(NCLUS)=NFLGFF
      NFLGPPC(NCLUS)=NFLGPPP
      NFLGBRMC(NCLUS)=NFLGBRMM
      NTMPFLG=1
      NCLTMP=NCLUS
C     ES(NCLUS)=ESEC
C RANDOMISE SECONDARY ELECTRON DIRECTION
C     R3=drand48(RDUM)
C     F3=1.0-2.0D0*R3
C     THETA0=DACOS(F3)
C     F6=DCOS(THETA0)
C     F5=DSIN(THETA0)
C     R4=drand48(RDUM)
C     PHI0=F4*R4
C     F8=DSIN(PHI0)
C     F9=DCOS(PHI0)               
C     DCX(NCLUS)=F9*F5
C     DCY(NCLUS)=F8*F5
C     DCZ(NCLUS)=F6   
C*********************************************************
      IF(IECASC.EQ.0)GO TO 333
      IF(LEGAS(I).EQ.0) GO TO 333
C USE COMPLETE CASCADE FOR ELECTRON IONISATION
      KG1=NEGAS(I)
      LG1=LEGAS(I)
      IGSHEL=IESHELL(I)
      CALL CASCADEE(J11,KG1,LG1,X,Y,Z,ST,ESEC,IGSHEL)
C
C STORE CASCADE IN CLUSTER STORE
C
      ETSUM=0.0
      DO 844 KBR=1,IEVENTE
      NCLUS=NCLUS+1
      IF(NCLUS.GT.150000) THEN
       WRITE(6,546) NCLUS,NREAL
       STOP
      ENDIF
      ES(NCLUS)=ECASE(KBR)
      ETSUM=ETSUM+ES(NCLUS)
      XS(NCLUS)=XCASE(KBR)
      YS(NCLUS)=YCASE(KBR)
      ZS(NCLUS)=ZCASE(KBR)
      TS(NCLUS)=TCASE(KBR)
      DCX(NCLUS)=DRXCE(KBR)
      DCY(NCLUS)=DRYCE(KBR)
      DCZ(NCLUS)=DRZCE(KBR)
      NFLGFC(NCLUS)=NFLGFE(KBR)+NFLGHIGH
      NFLGPPC(NCLUS)=NFLGPPE(KBR)
      NFLGBRMC(NCLUS)=NFLGBRMM
  844 CONTINUE
      IF(NFLGFC(NCLUS).GT.NFLGHIGH) NFLGHIGH=NFLGFC(NCLUS)
      GO TO 666
C*********************************************************
C STORE POSSIBLE SHELL EMISSIONS AUGER OR FLUORESCENCE 
  333 IF(EISTR.GT.30.0) THEN
C      WRITE(6,8891) EISTR
C8891  FORMAT(' EISTR=',D12.4)
C TEST IF FLUORESCENCE EMISSION
       IFLTST=0
       IF(WKLM(I).GT.0.0) THEN
        R9=drand48(RDUM)
        IF(R9.LT.WKLM(I)) IFLTST=1
       ENDIF
       IF(IFLTST.EQ.0) THEN
C AUGER EMISSION WITHOUT FLUORESCENCE
        NAUG=NC0(I)
        EAVAUG=EC0(I)/DFLOAT(NAUG)
        DO 700 JFL=1,NC0(I)
        NCLUS=NCLUS+1
        XS(NCLUS)=X
        YS(NCLUS)=Y
        ZS(NCLUS)=Z
        TS(NCLUS)=ST
        NFLGFC(NCLUS)=NFLGFF
        NFLGPPC(NCLUS)=NFLGPPP
        NFLGBRMC(NCLUS)=NFLGBRMM
        ES(NCLUS)=EAVAUG
        R3=drand48(RDUM)
        F3=1.0-2.0D0*R3
        THETA0=DACOS(F3)
        F6=DCOS(THETA0)
        F5=DSIN(THETA0)
        R4=drand48(RDUM)
        PHI0=F4*R4
        F8=DSIN(PHI0)
        F9=DCOS(PHI0)               
        DCX(NCLUS)=F9*F5
        DCY(NCLUS)=F8*F5
        DCZ(NCLUS)=F6
  700   CONTINUE 
       ELSE
C AUGER EMISSION AND FLUORESENCE 
        IF(NG2(I).EQ.0) GO TO 702
        NAUG=NG2(I)
        EAVAUG=EG2(I)/DFLOAT(NAUG)
        DO 701 JFL=1,NG2(I)
        NCLUS=NCLUS+1
        XS(NCLUS)=X
        YS(NCLUS)=Y
        ZS(NCLUS)=Z
        NFLGFC(NCLUS)=NFLGFF
        NFLGPPC(NCLUS)=NFLGPPP
        NFLGBRMC(NCLUS)=NFLGBRMM
        TS(NCLUS)=ST
        ES(NCLUS)=EAVAUG
        R3=drand48(RDUM)
        F3=1.0-2.0D0*R3
        THETA0=DACOS(F3)
        F6=DCOS(THETA0)
        F5=DSIN(THETA0)
        R4=drand48(RDUM)
        PHI0=F4*R4
        F8=DSIN(PHI0)
        F9=DCOS(PHI0)               
        DCX(NCLUS)=F9*F5
        DCY(NCLUS)=F8*F5
        DCZ(NCLUS)=F6
  701   CONTINUE
  702   IF(NG1(I).EQ.0) GO TO 704
        NAUG=NG1(I)
        EAVAUG=EG1(I)/DFLOAT(NAUG)
        R9=drand48(RDUM)
        DFL=-DLOG(R9)*DSTFL(I)
        DO 703 JFL=1,NG1(I)
        NCLUS=NCLUS+1
        R3=drand48(RDUM)
        THEFL=DACOS(1.0-2.0D0*R3)
        R4=drand48(RDUM)
        PHIFL=F4*R4
        XS(NCLUS)=X+DFL*DSIN(THEFL)*DCOS(PHIFL)
        YS(NCLUS)=Y+DFL*DSIN(THEFL)*DSIN(PHIFL)
        ZS(NCLUS)=Z+DFL*DCOS(THEFL)
        NFLGFC(NCLUS)=NFLGHIGH+1
        NFLGPPC(NCLUS)=NFLGPPP
        NFLGBRMC(NCLUS)=NFLGBRMM
        TS(NCLUS)=ST
        ES(NCLUS)=EAVAUG
        R3=drand48(RDUM)
        F3=1.0-2.0D0*R3
        THETA0=DACOS(F3)
        F6=DCOS(THETA0)
        F5=DSIN(THETA0)
        R4=drand48(RDUM)
        PHI0=F4*R4
        F8=DSIN(PHI0)
        F9=DCOS(PHI0)               
        DCX(NCLUS)=F9*F5
        DCY(NCLUS)=F8*F5
        DCZ(NCLUS)=F6
        NFLGHIGH=NFLGFC(NCLUS)
  703   CONTINUE
  704   CONTINUE
       ENDIF
      ENDIF
C                                                                       
C  GENERATE SCATTERING ANGLES AND UPDATE  LABORATORY COSINES AFTER      
C   COLLISION ALSO UPDATE ENERGY OF ELECTRON.                           
C
  666 IPT=IARRY(I)
      ICOLL(IPT)=ICOLL(IPT)+1 
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
