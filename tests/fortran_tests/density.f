      PROGRAM DENSTEST
      CALL SETUP
      CALL DENSITY
      END
      SUBROUTINE DENSITY
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*8 (I-N)
      COMMON/DENS/DEN(20000)
      COMMON/RATIO/AN1,AN2,AN3,AN4,AN5,AN6,AN,FRAC(6)
      COMMON/INPT/NGAS,NSTEP,NANISO,EFINAL,ESTEP,AKT,ARY,TEMPC,TORR,IPEN
      COMMON/GASN/NGASN(6)
      COMMON/MIX2/E(20000),EROOT(20000),QTOT(20000),QREL(20000),
     /QINEL(20000),QEL(20000)
      COMMON/RLTVY/BET(20000),GAM(20000),VC,EMS
      DIMENSION AND(6),EIAV(80),X00(80),X11(80),AKS(80),AAA(80),
     /JELEC(80)
C DENSITY EFFECT CONSTANTS
C EIAV ENERGY IN EV
C JELEC NUMBER OF ELECTRONS PER ATOM OR MOLECULE
      DATA EIAV/115.0,188.0,41.8,41.8,137.0,352.0,482.0,41.7,45.4,47.1,
     /48.3,85.0,0.0,71.6,95.0,82.0,0.0,84.9,0.0,0.0,
     /19.2,19.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,128.0,
     /53.7,0.0,0.0,67.6,62.9,61.1,0.0,0.0,0.0,0.0,
     /0.0,0.0,0.0,48.3,0.0,61.1,
     /34*0.0/
      DATA JELEC/42,18,2,2,10,36,54,10,18,26,
     /34,22,0,10,16,14,0,22,0,0,
     /2,2,0,0,0,0,0,0,0,70,
     /10,0,0,18,26,34,0,0,0,0,
     /0,0,0,34,0,34,
     /34*0/
      DATA X00/1.70,1.7635,2.2017,2.2017,2.0735,1.7158,1.5630,1.6263,
     /1.5090,1.4339,
     /1.3788,1.6294,0.0,1.7952,1.7541,1.7378,0.0,1.6477,0.0,0.0,
     /1.8639,1.8639,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.6,
     /1.6822,0.0,0.0,0.2529,0.2218,0.2046,0.0,0.0,0.0,0.0,
     /0.0,0.0,0.0,1.3788,0.0,0.2046,
     /34*0.0/
      DATA X11/4.00,4.4855,3.6122,3.6122,4.6421,5.0748,4.7371,3.9716,
     /3.8726,3.8011,
     /3.7524,4.1825,0.0,4.3437,4.3213,4.1323,0.0,4.1565,0.0,0.0,
     /3.2718,3.2718,0.0,0.0,0.0,0.0,0.0,0.0,0.0,4.0,
     /4.1158,0.0,0.0,2.7639,2.7052,2.6681,0.0,0.0,0.0,0.0,
     /0.0,0.0,0.0,3.7524,0.0,2.6681,
     /34*0.0/
      DATA AKS/3.00,2.9618,5.8347,5.8347,3.5771,3.4051,2.7414,3.6257,
     /3.6095,3.5920,
     /3.4884,3.3227,0.0,3.5901,3.2913,3.2125,0.0,3.3318,0.0,0.0,
     /5.7273,5.7273,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,
     /3.6464,0.0,0.0,3.5477,3.4834,3.5415,0.0,0.0,0.0,0.0,
     /0.0,0.0,0.0,3.4884,0.0,3.5415,
     /34*0.0/
      DATA AAA/.18551,.19714,.13443,.13443,.08064,.07446,.23314,.09253,
     /0.09627,0.09916,
     /.10852,.11768,0.0,.08101,.11778,.15349,0.0,.11992,0.0,0.0,
     /.14092,.14092,0.0,0.0,0.0,0.0,0.0,0.0,0.0,.177484,
     /.08315,0.0,0.0,.08970,0.09878,0.09644,0.0,0.0,0.0,0.0,
     /0.0,0.0,0.0,.10852,0.0,0.09644,
     /34*0.0/
C
      API=DACOS(-1.0D0)                                                 
      EMS=510998.928D0
      RE=2.8179403267D-13    
      ALPH=137.035999074
      ABZERO=273.15D0                                                   
      ATMOS=760.0D0                                                     
C                                                                       
C DENSITY EFFECT CALCULATION
      AND(1)=AN1
      AND(2)=AN2
      AND(3)=AN3
      AND(4)=AN4
      AND(5)=AN5
      AND(6)=AN6
      HSUM=0.0
      SUM1=0.0
      SUMDNOM=0.0
      DO 120 L1=1,NGAS
      SUM1=SUM1+FRAC(L1)*DFLOAT(JELEC(NGASN(L1)))*DLOG(EIAV(NGASN(L1))) 
      SUMDNOM=SUMDNOM+FRAC(L1)*DFLOAT(JELEC(NGASN(L1)))
  120 HSUM=HSUM+AND(L1)*DFLOAT(JELEC(NGASN(L1)))
      EIBAR=DEXP(SUM1/SUMDNOM)
C PLASMA ENERGY
      HWP1=DSQRT(4.0*API*HSUM*RE**3)*ALPH*EMS
C
      DELDEN=DLOG(EIBAR/HWP1)
      CBAR=1.0+2.0*DELDEN
      IF(NGAS.EQ.1) GO TO 200
C CALC X0 AND X1
      IF(CBAR.LT.10.0) THEN
      X0=1.6
      X1=4.0
      ELSE IF(CBAR.GE.4.0.AND.CBAR.LT.10.5) THEN
      X0=1.7
      X1=4.0
      ELSE IF(CBAR.GE.10.5.AND.CBAR.LT.11.0) THEN
      X0=1.8
      X1=4.0
      ELSE IF(CBAR.GE.11.0.AND.CBAR.LT.11.5) THEN
      X0=1.9
      X1=4.0
      ELSE IF(CBAR.GE.11.5.AND.CBAR.LT.12.25) THEN
      X0=2.0
      X1=4.0
      ELSE IF(CBAR.GE.12.25.AND.CBAR.LT.13.804) THEN
      X0=2.0
      X1=5.0
      ELSE 
      X0=0.326*CBAR-1.5
      X1=5.0
      ENDIF
      AKBAR=3.0
      ABAR=(CBAR-2.0*DLOG(10.0D0)*X0)/((X1-X0)**3)
      GO TO 201
  200 AKBAR=AKS(NGASN(1))
      X0=X00(NGASN(1))
      X1=X11(NGASN(1))
      ABAR=AAA(NGASN(1))
  201 CONTINUE
C CORRECT X0 AND X1 FOR DENSITY CHANGE FROM 20C AND 760 TORR
C NB CORRECTION TO CBAR ALREADY DONE
      DCOR=0.5*DLOG10(TORR*293.15/(760.0*(TEMPC+ABZERO)))
      X0=X0-DCOR
      X1=X1-DCOR
C CALCULATE DENSITY CORRECTION FACTOR ARRAY DEN(20000)
      AFC=2.0*DLOG(10.0D0)
      DO 236 I=1,20000
      BG=BET(I)*GAM(I)
      X=DLOG10(BG)
      IF(X.LT.X0) THEN   
       DEN(I)=0.0
      ELSE IF(X.GT.X0.AND.X.LT.X1) THEN
       DEN(I)=ABAR*DEXP(AKBAR*DLOG(X1-X))+AFC*X-CBAR
      ELSE 
       DEN(I)=AFC*X-CBAR              
      ENDIF
C      WRITE(50,99) E(I)
C  99  FORMAT(D12.5)
  236 CONTINUE
      RETURN
      END
      SUBROUTINE SETUP                                      
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*8 (I-N) 
      INTEGER*4 NSEED                                       
      COMMON/INPT/NGAS,NSTEP,NANISO,EFINAL,ESTEP,AKT,ARY,TEMPC,TORR,IPEN
      COMMON/CNSTS/ECHARG,EMASS,AMU,PIR2
      COMMON/INPT2/KGAS,LGAS,DETEFF,EXCWGHT
      COMMON/INPT1/NDVEC                                
      COMMON/CNSTS1/CONST1,CONST2,CONST3,CONST4,CONST5                  
      COMMON/RATIO/AN1,AN2,AN3,AN4,AN5,AN6,AN,FRAC(6)               
      COMMON/GASN/NGASN(6)                                 
      COMMON/SETP/TMAX,SMALL,API,ESTART,THETA,PHI,TCFMAX(10),TCFMAX1,
     /RSTART,EFIELD,ETHRM,ECUT,NEVENT,IMIP,IWRITE
      COMMON/SET2/DRXINIT,DRYINIT,DRZINIT
      COMMON/BFLD/EOVB,WB,BTHETA,BMAG 
      COMMON/IONC/DOUBLE(6,20000),CMINIXSC(6),CMINEXSC(6),ECLOSS(6),
     /WPLN(6),ICOUNT,AVPFRAC(3,6)
      COMMON/MRATIO/VAN1,VAN2,VAN3,VAN4,VAN5,VAN6,VAN
      COMMON/OUTPT/ICOLL(30),NETOT,NPRIME,TMAX1,TIME(300),NNULL,
     /NITOT,ICOLN(512),ICOLNN(60),NREAL,NEXCTOT
      COMMON/PRIM3/MSUM(10000),MCOMP(10000),MRAYL(10000),MPAIR(10000),
     /MPHOT(10000),MVAC(10000)
      COMMON/RLTVY/BET(20000),GAM(20000),VC,EMS 
      COMMON/COMP/ICMP,ICFLG,IRAY,IRFLG,IPAP,IPFLG,IBRM,IBFLG,LPEFLG 
      COMMON/MIX2/E(20000),EROOT(20000),QTOT(20000),QREL(20000),
     /QINEL(20000),QEL(20000)
      COMMON/PLOT/NXPL10(31),NYPL10(31),NZPL10(31),NXPL40(31),
     /NYPL40(31),NZPL40(31),NXPL100(31),NYPL100(31),NZPL100(31),
     /NXPL400(31),NYPL400(31),NZPL400(31),NXPL1000(31),NYPL1000(31),
     /NZPL1000(31),NXPL2(31),NYPL2(31),NZPL2(31),NXPL4000(31),
     /NYPL4000(31),NZPL4000(31),NXPL10000(31),NYPL10000(31),
     /NZPL10000(31),NXPL40000(31),NYPL40000(31),NZPL40000(31),
     /NXPL100000(31),NYPL100000(31),NZPL100000(31),NRPL2(31),NRPL10(31),
     /NRPL40(31),NRPL100(31),NRPL400(31),NRPL1000(31),NRPL4000(31),
     /NRPL10000(31),NRPL40000(31),NRPL100000(31),NEPL1(100),
     /NEPL10(100),NEPL100(100),MELEC(1000),MELEC3(1000),MELEC10(1000),
     /MELEC30(1000),MELEC100(1000),MELEC300(1000)
      COMMON/BREMG/EBRGAM(10),BRDCOSX(10),BRDCOSY(10),BRDCOSZ(10),
     /BRX(10),BRY(10),BRZ(10),BRT(10),EBRTOT(6),NBREM(6)
      COMMON/CLUS/XAV(100000),YAV(100000),ZAV(100000),TAV(100000),
     /XYAV(100000),XYZAV(100000),DX(100000),DY(100000),DZ(100000),
     /DT(100000),DXY(100000),DXYZ(100000),NCL(100000),FARX1(100000)
     /,FARY1(100000),FARZ1(100000),FARXY1(100000),RMAX1(100000),
     /TSUM(100000),XNEG(100000), 
     /YNEG(100000),ZNEG(100000),EDELTA(100000),EDELTA2(100000),
     /NCLEXC(100000)
      COMMON/KSEED/NSEED
      COMMON/ECASC/NEGAS(512),LEGAS(512),IESHELL(512),IECASC
C                                                                       
C   NEW UPDATE OF CONSTANTS 2010
C
      API=DACOS(-1.0D0)                                                 
      ARY=13.60569253D0                                              
      PIR2=8.7973554297D-17                                       
      ECHARG=1.602176565D-19                                            
      EMASS=9.10938291D-31                     
      EMS=510998.928D0
      VC=299792458.0D0                       
      AMU=1.660538921D-27                                             
      BOLTZ=8.6173324D-5     
      BOLTZJ=1.3806488D-23                                              
      AWB=1.758820088D10                                              
      ALOSCH=2.6867805D19      
      RE=2.8179403267D-13    
      ALPH=137.035999074
      HBAR=6.58211928D-16                                     
      EOVM=DSQRT(2.0D0*ECHARG/EMASS)*100.0D0                            
      ABZERO=273.15D0                                                   
      ATMOS=760.0D0                                                     
      CONST1=AWB/2.0D0*1.0D-19                                          
      CONST2=CONST1*1.0D-02                                             
      CONST3=DSQRT(0.2D0*AWB)*1.0D-09                                   
      CONST4=CONST3*ALOSCH*1.0D-15                                      
      CONST5=CONST3/2.0D0
      TWOPI=2.0D0*API
      NANISO=2
      DO 55 K=1,6
      NBREM(K)=0
      EBRTOT(K)=0.0
   55 CONTINUE
      ICFLG=0
      IRFLG=0
      IPFLG=0
      IBFLG=0
      LPEFLG=0
C  --------------------------------------------       
C                                                                       
C      READ IN OUTPUT CONTROL AND INTEGRATION DATA                      
C                                                                       
      READ(5,2) NGAS,NEVENT,IMIP,NDVEC,NSEED,ESTART,ETHRM,ECUT    
    2 FORMAT(5I10,3F10.5)  
      ICOUNT=0
      IF(IMIP.EQ.1) ICOUNT=1 
      IF(NGAS.EQ.0) GO TO 99 
      IF(ESTART.GT.3.0D6.AND.IMIP.EQ.3) THEN
      WRITE(6,664) ESTART
  664 FORMAT(' PROGRAM STOPPED: X-RAY ENERGY=',D12.3,'EV. MAXIMUM ENERGY
     / 3.0MEV')
       STOP 
      ENDIF
      IF(IMIP.NE.1.AND.NEVENT.GT.10000) THEN 
       WRITE(6,665) NEVENT
  665  FORMAT(' PROGRAM STOPPED NUMBER OF EVENTS =',I7,' LARGER THAN ARR
     /AY LIMIT OF 10000')
       STOP
      ENDIF
      IF(IMIP.EQ.1.AND.NEVENT.GT.100000) THEN
       WRITE(6,666) NEVENT
  666  FORMAT(' PROGRAM STOPPED NUMBER OF EVENTS =',I7,' LARGER THAN ARR
     /AY LIMIT OF 100000')
       STOP
      ENDIF
C 
C   GAS IDENTIFIERS 
C
      READ(5,3) NGASN(1),NGASN(2),NGASN(3),NGASN(4),NGASN(5),NGASN(6)
    3 FORMAT(6I5)        
C      
C      GAS PARAMETERS
C
      READ(5,4) FRAC(1),FRAC(2),FRAC(3),FRAC(4),FRAC(5),FRAC(6),TEMPC,
     /TORR                        
    4 FORMAT(8F10.5)      
C                                                  
C      FIELD VALUES                                                    
C                                                                       
      READ(5,5) EFIELD,BMAG,BTHETA,IWRITE,IPEN                         
    5 FORMAT(3F10.3,2I5)
      READ(5,6) DETEFF,EXCWGHT,KGAS,LGAS,ICMP,IRAY,IPAP,IBRM,IECASC 
    6 FORMAT(2F10.3,7I5)
C     WRITE(6,656) IWRITE
C 656 FORMAT(' IWRITE=',I3)  
      IF(IWRITE.NE.0) OPEN(UNIT=50,FILE='DENSITY.OUT')
C CALCULATE EFINAL FOR DELTAS OR XRAYS 
C INCREASED EFINAL CAUSED BY ELECTRIC FIELD 
      EBIG=0.05*ESTART/1000. 
      EFINAL=ESTART*1.0001+760.0*EBIG/TORR*(TEMPC+ABZERO)/293.15*EFIELD
      IF(EFINAL.LT.(1.01*ESTART)) EFINAL=1.01*ESTART 
C      WRITE(50,20) EFINAL
C   20 FORMAT(1D12.4)

C   CHECK INPUT
      TOTFRAC=0.0D0
      IF(NGAS.EQ.0.OR.NGAS.GT.6) GO TO 999
      DO 10 J=1,NGAS
      IF(NGASN(J).EQ.0.OR.FRAC(J).EQ.0.0D0) GO TO 999
   10 TOTFRAC=TOTFRAC+FRAC(J)
      IF(DABS(TOTFRAC-100.0D0).GT.1.D-6) GO TO 999
      LAST=0
      TMAX=100.0D0  
      NOUT=10  
      NSTEP=20000
C INITIAL ANGLES
      IF(NDVEC.EQ.1) THEN
       PHI=0.0D0                                
       THETA=0.0D0 
      ELSE IF(NDVEC.EQ.(-1)) THEN
       PHI=0.0D0
       THETA=DACOS(-1.D0)
      ELSE IF(NDVEC.EQ.0) THEN
       PHI=0.0D0
       THETA=API/2.0  
      ELSE 
       WRITE(6,992) NDVEC
  992  FORMAT(/,2X,'DIRECTION OF BEAM NOT DEFINED NDVEC =',I5)
       STOP      
      ENDIF
C INITIAL DIRECTION COSINES FOR CASCADE CALCULATION
      DRZINIT=DCOS(THETA)
      DRXINIT=DSIN(THETA)*DCOS(PHI)
      DRYINIT=DSIN(THETA)*DSIN(PHI)
C ZERO COMMON BLOCKS OF OUTPUT RESULTS
      DO 64 J=1,10000
      MSUM(J)=0
      MCOMP(J)=0
      MRAYL(J)=0
      MPAIR(J)=0
      MPHOT(J)=0
   64 MVAC(J)=0
      DO 65 J=1,300                                                     
   65 TIME(J)=0.0D0                                                     
      DO 70 K=1,30                                                      
   70 ICOLL(K)=0  
      DO 80 K=1,512
   80 ICOLN(K)=0                 
      DO 81 K=1,60
   81 ICOLNN(K)=0                                       
      DO 100 K=1,10                                                     
  100 TCFMAX(K)=0.0D0   
C ZERO PLOT ARRAYS
      DO 110 K=1,31
      NXPL2(K)=0
      NYPL2(K)=0
      NZPL2(K)=0
      NXPL10(K)=0
      NYPL10(K)=0
      NZPL10(K)=0
      NXPL40(K)=0
      NYPL40(K)=0
      NZPL40(K)=0
      NXPL100(K)=0
      NYPL100(K)=0
      NZPL100(K)=0
      NXPL400(K)=0
      NYPL400(K)=0
      NZPL400(K)=0
      NXPL1000(K)=0
      NYPL1000(K)=0
      NZPL1000(K)=0
      NXPL4000(K)=0
      NYPL4000(K)=0
      NZPL4000(K)=0
      NXPL10000(K)=0
      NYPL10000(K)=0
      NZPL10000(K)=0
      NXPL40000(K)=0
      NYPL40000(K)=0
      NZPL40000(K)=0
      NXPL100000(K)=0
      NYPL100000(K)=0
      NZPL100000(K)=0
      NRPL2(K)=0
      NRPL10(K)=0
      NRPL40(K)=0
      NRPL100(K)=0
      NRPL400(K)=0
      NRPL1000(K)=0
      NRPL4000(K)=0
      NRPL10000(K)=0
      NRPL40000(K)=0
  110 NRPL100000(K)=0
      DO 111 K=1,100
      NEPL1(K)=0
      NEPL10(K)=0
  111 NEPL100(K)=0
      DO 112 K=1,1000
      MELEC(K)=0
      MELEC3(K)=0
      MELEC10(K)=0
      MELEC30(K)=0
      MELEC100(K)=0
  112 MELEC300(K)=0
C ZERO ARRAYS
      DO 113 KS=1,100000
      XAV(KS)=0.0
      YAV(KS)=0.0
      ZAV(KS)=0.0
      TAV(KS)=0.0
      XYAV(KS)=0.0
      XYZAV(KS)=0.0
      DX(KS)=0.0
      DY(KS)=0.0
      DZ(KS)=0.0
      DT(KS)=0.0
      DXY(KS)=0.0
      DXYZ(KS)=0.0
      FARX1(KS)=0.0
      FARY1(KS)=0.0
      FARZ1(KS)=0.0
      FARXY1(KS)=0.0
      RMAX1(KS)=0.0
      TSUM(KS)=0.0
      XNEG(KS)=0.0
      YNEG(KS)=0.0
      ZNEG(KS)=0.0
      EDELTA(KS)=0.0
      EDELTA2(KS)=0.0
      NCL(KS)=0
      NCLEXC(KS)=0
  113 CONTINUE
      CORR=ABZERO*TORR/(ATMOS*(ABZERO+TEMPC)*100.0D0)                   
      AKT=(ABZERO+TEMPC)*BOLTZ
      AN1=FRAC(1)*CORR*ALOSCH                                           
      AN2=FRAC(2)*CORR*ALOSCH                                           
      AN3=FRAC(3)*CORR*ALOSCH                                           
      AN4=FRAC(4)*CORR*ALOSCH
      AN5=FRAC(5)*CORR*ALOSCH
      AN6=FRAC(6)*CORR*ALOSCH                                           
      AN=100.0D0*CORR*ALOSCH                                            
C     VAN1=FRAC(1)*CORR*CONST4*1.0D15                                   
C     VAN2=FRAC(2)*CORR*CONST4*1.0D15                                   
C     VAN3=FRAC(3)*CORR*CONST4*1.0D15                                   
C     VAN4=FRAC(4)*CORR*CONST4*1.0D15
C     VAN5=FRAC(5)*CORR*CONST4*1.0D15
C     VAN6=FRAC(6)*CORR*CONST4*1.0D15                                   
C     VAN=100.0D0*CORR*CONST4*1.0D15
      VAN1=FRAC(1)*CORR*ALOSCH*VC                                   
      VAN2=FRAC(2)*CORR*ALOSCH*VC                                   
      VAN3=FRAC(3)*CORR*ALOSCH*VC                                  
      VAN4=FRAC(4)*CORR*ALOSCH*VC
      VAN5=FRAC(5)*CORR*ALOSCH*VC
      VAN6=FRAC(6)*CORR*ALOSCH*VC                                  
      VAN=100.0D0*CORR*ALOSCH*VC
C CALCULATE AND STORE ENERGY GRID FOR XRAYS BETAS OR PARTICLES
      IF(EFINAL.LE.20000.0) THEN
       ESTEP=EFINAL/DFLOAT(NSTEP)
       EHALF=ESTEP/2.0D0
       E(1)=EHALF
       GAM(1)=(EMS+E(1))/EMS
       BET(1)=DSQRT(1.0D0-1.0D0/(GAM(1)*GAM(1)))
       DO 203 I=2,20000
       AJ=DFLOAT(I-1)
       E(I)=EHALF+ESTEP*AJ
       GAM(I)=(EMS+E(I))/EMS
  203  BET(I)=DSQRT(1.0D0-1.0D0/(GAM(I)*GAM(I)))
      ELSE IF(EFINAL.GT.20000.0.AND.EFINAL.LE.140000.) THEN
       ESTEP=1.0
       EHALF=0.5
       E(1)=EHALF
       GAM(1)=(EMS+E(1))/EMS
       BET(1)=DSQRT(1.0D0-1.0D0/(GAM(1)*GAM(1)))
       DO 231 I=2,16000
       AJ=DFLOAT(I-1)
       E(I)=EHALF+ESTEP*AJ
       GAM(I)=(EMS+E(I))/EMS
  231  BET(I)=DSQRT(1.0D0-1.0D0/(GAM(I)*GAM(I)))
       ESTEP1=(EFINAL-16000.0)/DFLOAT(4000)
       DO 232 I=16001,20000
       AJ=DFLOAT(I-16000)
       E(I)=16000.0+AJ*ESTEP1
       GAM(I)=(EMS+E(I))/EMS
  232  BET(I)=DSQRT(1.0D0-1.0D0/(GAM(I)*GAM(I)))
      ELSE
       ESTEP=1.0
       EHALF=0.5
       E(1)=EHALF
       GAM(1)=(EMS+E(1))/EMS
       BET(1)=DSQRT(1.0D0-1.0D0/(GAM(1)*GAM(1)))
       DO 233 I=2,12000
       AJ=DFLOAT(I-1)
       E(I)=EHALF+ESTEP*AJ
       GAM(I)=(EMS+E(I))/EMS
  233  BET(I)=DSQRT(1.0D0-1.0D0/(GAM(I)*GAM(I)))
       ESTEP1=20.0
       DO 234 I=12001,16000
       AJ=DFLOAT(I-12000)
       E(I)=12000.0+AJ*ESTEP1
       GAM(I)=(EMS+E(I))/EMS
  234  BET(I)=DSQRT(1.0D0-1.0D0/(GAM(I)*GAM(I)))
       ESTEP2=(EFINAL-92000.0)/DFLOAT(4000)
       DO 235 I=16001,20000
       AJ=DFLOAT(I-16000)
       E(I)=92000.0+AJ*ESTEP2
       GAM(I)=(EMS+E(I))/EMS
  235  BET(I)=DSQRT(1.0D0-1.0D0/(GAM(I)*GAM(I)))
      ENDIF
C  RADIANS PER PICOSECOND                                        
      WB=AWB*BMAG*1.0D-12 
C   METRES PER PICOSECOND
      IF(BMAG.EQ.0.0D0) RETURN
      EOVB=EFIELD*1.D-9/BMAG
      RETURN
  999 WRITE(6,87) NGAS,(J,NGASN(J),FRAC(J),J=1,6) 
   87 FORMAT(3(/),4X,' ERROR IN GAS INPUT : NGAS=',I5,6(/,2X,' N=',I3,' 
     /NGAS=',I5,' FRAC=',F8.3))                                         
   99 LAST=1                                                            
      RETURN                                                            
      END                                                               
