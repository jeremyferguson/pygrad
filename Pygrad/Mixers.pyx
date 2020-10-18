from Pygrad cimport Pygrad
from libc.math cimport sin, cos, acos, asin, log, sqrt
from PyGasMix.Gasmix cimport Gasmix
from Ang cimport Ang
import sys
import cython

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef Mixer(Pygrad object):
    """Loads the initial gas values for the mixture, particularly the momentum cross sections and other related values."""
    cdef double AttachmentCrossSection[6][4000], ElectronCharge, JHI, JLOW, EnergyHigh, F2, BP, EnergyLow
    cdef int iEnergy, GasIndex, iProcess, p, Sum, J, i, j, iIonization, JJ, IL, I

    ElectronCharge = 1.60217656e-19
    object.GasMix.InitWithInfo(object.GasIDs, object.InelasticCrossSectionPerGas, object.N_Inelastic,
            object.PenningFraction, object.E, object.SqrtEnergy, object.NumberOfGases, object.EnergySteps, object.Which_Angular_Model, object.ElectronEnergyStep, object.Max_Electron_Energy, object.ThermalEnergy, object.TemperatureCentigrade, object.Pressure_Torr, object.PIR2, object.RhydbergConst)
    object.GasMix.Run()

    ang = Ang()
    for iEnergy in range(20000):
        object.IonCollisionFreq[iEnergy] = 0.0
        object.AttCollisionFreq[iEnergy] = 0.0
        nProcess = 1
        for idg in range(object.NumberOfGases):
            object.GasExcitationSlots[idg] = 1
            object.negas[nProcess - 1] = 1
            object.legas[nProcess - 1] = 0
            object.ieshell[nProcess - 1] = 0
            object.collisionFrequency[iEnergy][nProcess - 1] = object.Q[idg][1][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * Object.Beta[iEnergy]
            object.ScatteringParameter[iEnergy][nProcess - 1] = 0.5
            object.AngularCut[iEnergy][nProcess - 1] = 0.0
            if object.GasMix.Gases[idg].AngularModel[1] == 1:
                ang.ScatteringParameter1 = object.GasMix.Gases[idg].PEElasticCrossSection[1][iEnergy]
                ang.CalcAngCut()
                object.AngularCut[iEnergy][nProcess - 1] = ang.AngCut
                object.ScatteringParameter[iEnergy][nProcess - 1] = ang.ScatteringParameter2
                object.AngularModel[nProcess - 1] = 1
            elif object.GasMix.Gases[idg].AngularModel[1] == 2:
                object.ScatteringParameter[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].PEElasticCrossSection[1][iEnergy]
                object.AngularModel[nProcess - 1] = 2
            if iEnergy <= 1:
                rGas1 = 1.0 + object.GasMix.Gases[idg].E[1]/2.0
                object.rGas[nProcess - 1] = rGas1
                object.Ein[nProcess - 1] = 0.0
                object.Ipn[nProcess - 1] = 0
                L = 1
                object.InteractionType[nProcess - 1] = L
                object.izbr[nProcess - 1] = 0
                object.PenningFraction[0][nProcess - 1] = 0.0
                object.PenningFraction[1][nProcess - 1] = 0.0
                object.PenningFraction[2][nProcess - 1] = 0.0
                object.avpfrac[0][nProcess - 1] = 0.0
                object.avpfrac[1][nProcess - 1] = 0.0
                object.avpfrac[2][nProcess - 1] = 0.0
                object.cminexsc[0] = object.GasMix.Gases[idg].E[3] * object.MoleculesPerCm3PerGas[idg]
                object.cminixsc[0] = object.GasMix.Gases[idg].E[4] * object.MoleculesPerCm3PerGas[idg]
                object.ecloss[0] = object.GasMix.Gases[idg].E[2]
                object.WPLN[0] = object.GasMix.Gases[idg].E[5]
            if object.FinalElectronEnergy >= object.GasMix.Gases[idg].E[2]:
                if object.GasMix.Gases[idg].N_Ionization <= 1:
                    nProcess += 1
                    object.GasExcitationSlots[idg] = nProcess
                    if object.icount == 1:
                        QIndex = 4
                        object.DOUBLE[0][iEnergy] = object.GasMix.Gases[idg].Q[2][iEnergy] / object.GasMix.Gases[idg].Q[4][iEnergy] - 1.0
                    else:
                        QIndex = 2
                    object.CollisionFrequency[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].Q[QIndex][iEnergy] * object.VMoleculesPerCm3PerGas * object.Beta[iEnergy]
                    object.IonCollisionFreq[iEnergy] += object.CollisionFrequency[iEnergy][nProcess - 1]
                    object.negas[nProcess - 1] = 1
                    object.legas[nProcess - 1] = 0
                    object.ieshell[nProcess - 1] = 0
                    object.ScatteringParameter[iEnergy][nProcess - 1] = 0.5
                    object.AngleCut[iEnergy][nProcess - 1] = 1.0
                    object.AngularModel[nProcess - 1] = 0
                    if object.icount == 1:
                        pIndex = 4
                    else:
                        pIndex = 2
                    object.AngularModel[nProcess - 1] = 2
                    if object.GasMix.Gases[idg].AngularModel[pIndex] == 1:
                        Ang.ScatteringParameter1 = object.GasMix.Gases[idg].PEElasticCrossSection[pIndex][iEnergy]
                        Ang.CalcAngCut()
                        object.AngleCut[iEnergy][nProcess - 1] = Ang.AngCut
                        object.ScatteringParameter = Ang.ScatteringParameter2
                        object.AngularModel[nProcess - 1] = 1
                    object.wpl[nProcess - 1] = object.GasMix.Gases[idg].EB1[0]
                    object.nc0[nProcess - 1] = object.GasMix.Gases[idg].NC0[0]
                    object.ec0[nProcess - 1] = object.GasMix.Gases[idg].EC0[0]
                    object.ng1[nProcess - 1] = object.GasMix.Gases[idg].NG1[0]
                    object.eg1[nProcess - 1] = object.GasMix.Gases[idg].EG1[0]
                    object.ng2[nProcess - 1] = object.GasMix.Gases[idg].NG2[0]
                    object.eg2[nProcess - 1] = object.GasMix.Gases[idg].EG2[0]
                    object.wklm[nProcess - 1] = object.GasMix.Gases[idg].WK[0]
                    object.efl[nProcess - 1] = object.GasMix.Gases[idg].EFL[0]
                    if iEnergy == 1:
                        object.rGas[nProcess - 1] = rgas1
                        object.ein[nProcess - 1] = object.GasMix.Gases[idg]
                        object.ipn[nProcess - 1] = 1
                        L = 2
                        object.InteractionType[nProcess -1] = L
                        object.izbr[nProcess - 1] = 0
                        object.PenningFraction[nProcess - 1][0] = 0.0
                        object.PenningFraction[nProcess - 1][1] = 0.0
                        object.PenningFraction[nProcess - 1][2] = 0.0
                else:
                    for kion in range(object.GasMix.Gases[idg].N_Ionization):
                        nProcess += 1
                        object.GasExcitationSlots[idg] = nProcess
                        object.CollisionFrequency[iEnergy,nProcess - 1] = object.GasMix.Gases[idg].IonizationCrossSection[kin][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * object.Beta[iEnergy]
                        object.IonCollisionFreq[iEnergy] += object.CollisionFrequency[iEnergy,nProcess - 1]
                        object.ScatteringParameter[iEnergy,nProcess - 1] = 0.5
                        object.AngleCut[iEnergy,nProcess - 1] = 1.0
                        object.AngularModel[nProcess - 1] = 0
                        object.negas[nProcess - 1] = 1
                        #object.legas[nProcess - 1] = object.GasMix.Gases[idg].Legas[kion]
                        #object.ieshell[nProcess - 1] = object.GasMix.Gases[idg].IEshell[kion]
                        if object.GasMix.Gases[idg].AngularModel[2] == 1:
                            Ang.ScatteringParameter1[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].PEInelasticCrossSection[kion][iEnergy] 
                            Ang.CalcAngCut()
                            object.AngleCut[iEnergy][nProcess - 1] = Ang.AngCut
                            object.ParameterScattering[iEnergy][nProcess - 1] = Ang.ScatteringParameter2
                            object.AngularModel[nProcess - 1] = 1
                        if object.GasMix.Gases[idg].AngularModel[2] == 2:
                            object.ParameterScattering[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].PEInelasticCrossSection[kion][iEnergy]
                            object.AngularModel[nProcess - 1] = 2
                        object.wpl[nProcess - 1] = object.GasMix.Gases[idg].EB1[kion]
                        object.nc0[nProcess - 1] = object.GasMix.Gases[idg].NC0[kion]
                        object.ec0[nProcess - 1] = object.GasMix.Gases[idg].EC0[kion]
                        object.ng1[nProcess - 1] = object.GasMix.Gases[idg].NG1[kion]
                        object.eg1[nProcess - 1] = object.GasMix.Gases[idg].EG1[kion]
                        object.ng2[nProcess - 1] = object.GasMix.Gases[idg].NG2[kion]
                        object.eg2[nProcess - 1] = object.GasMix.Gases[idg].EG2[kion]
                        object.wklm[nProcess - 1] = object.GasMix.Gases[idg].WK[kion]
                        object.efl[nProcess - 1] = object.GasMix.Gases[idg].EFL[kion]
                        if iEnergy <= 1:
                            object.rGas[nProcess - 1] = rGas1
                            object.ein[nProcess - 1] = object.GasMix.Gases[idg].IonizationEnergy[kion]/rGas1
                            object.ipn[nProcess - 1] = 1
                            L = 2
                            object.InteractionType[nProcess - 1] = L
                            object.izbr[nProcess - 1] = 0
                            object.PenningFraction[nProcess - 1][0] = 0.0
                            object.PenningFraction[nProcess - 1][1] = 0.0
                            object.PenningFraction[nProcess - 1][2] = 0.0
                            #Object.Ionmodel[nProcess - 1][nProcess - 1] = ionmodl1
                            #for i in range(20):
                                #object.esplit[nProcess - 1][i] = object.GasMix.Gases[idg].Esplit[i]
            if object.FinalElectronEnergy >= object.GasMix.Gases[idg].E[3]:
                if object.GasMix.Gases[idg].N_Attachment <= 1:
                    nProcess += 1
                    object.GasExcitationSlots[idg] = nProcess
                    object.CollisionFrequency[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].Q[3][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * object.Beta[iEnergy]
                    object.AttCollisionFreq[iEnergy] += object.CollisionFrequency[iEnergy][nProcess - 1]
                    object.ScatteringParameter[iEnergy][nProcess - 1] = 0.5
                    object.AngleCut[iEnergy][nProcess - 1] = 1.0
                    if iEnergy <= 1:
                        object.negas[nProcess - 1] = 1
                        object.legas[nProcess - 1] = 0
                        object.ieshell[nProcess - 1] = 0
                        object.AngularModel[nProcess - 1] = 0
                        object.rGas[nProcess - 1] = rGas1
                        object.ein[nProcess - 1] = 0.0
                        object.ipn[nProcess - 1] = -1
                        L = 3
                        object.InteractionType[nProcess - 1] = L
                        object.izbr[nProcess - 1] = 0
                        object.PenningFraction[nProcess - 1][0] = 0.0
                        object.PenningFraction[nProcess - 1][1] = 0.0
                        object.PenningFraction[nProcess - 1][2] = 0.0
                    else:
                        for jj in range(object.GasMix.Gases[idg].N_Attachment):
                            nProcess += 1
                            object.GasExcitationSlots[idg] = nProcess
                            object.CollisionFrequency[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].Q[3][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * object.Beta[iEnergy]
                            object.AttCollisionFreq[iEnergy] += object.CollisionFrequency[iEnergy][nProcess - 1]
                            object.ScatteringParameter[iEnergy][nProcess - 1] = 0.5
                            object.AngleCut[iEnergy][nProcess - 1] = 1.0
                            if iEnergy <= 1:
                                object.negas[nProcess - 1] = 1
                                object.legas[nProcess - 1] = 0
                                object.ieshell[nProcess - 1] = 0
                                object.AngularModel[nProcess - 1] = 0
                                object.rGas[nProcess - 1] = rGas1
                                object.ein[nProcess - 1] = 0.0
                                object.ipn[nProcess - 1] = -1
                                L = 3
                                object.InteractionType[nProcess - 1] = L
                                object.izbr[nProcess - 1] = 0
                                object.PenningFraction[nProcess - 1][0] = 0.0
                                object.PenningFraction[nProcess - 1][1] = 0.0
                                object.PenningFraction[nProcess - 1][2] = 0.0
                        if object.GasMix.Gases[idg].N_Inelastic != 0:
                            for j in range(object.GasMix.Gases[idg].N_Inelastic):
                                nProcess += 1
                                object.GasExcitationSlots[idg] = nProcess
                                object.negas[nProcess - 1] = 1
                                object.legas[nProcess - 1] = 0
                                object.ieshell[nProcess - 1] = 0
                                object.CollisionFrequency[iEnergy][nProcess - 1] = object.InelasticCrossSection[j][iEnergy] * object.VMoleculesPerCm3PerGas[0] * object.Beta[iEnergy]
                                #if object.GasMix.Gases[idg].izbr[j] != 0 and object.lbrm == 0:
                                    #object.CollisionFrequency[iEnergy][nProcess - 1] = 0.0
                                object.ScatteringParameter[iEnergy][nProcess] = 0.5
                                object.AngleCut[iEnergy][nProcess] = 1.0
                                object.AngularModel[nProcess - 1] = 0
                                if object.GasMix.Gases[idg].KIN[j] == 1:
                                    Ang.ScatteringParameter1 = object.GasMix.Gases[idg].PEInelasticCrossSection[j][iEnergy]
                                    Ang.CalcAngCut()
                                    object.ScatteringParameter[iEnergy][nProcess - 1] = Ang.ScatteringParameter2
                                    object.AngleCut[iEnergy][nProcess - 1] = Ang.AngleCut
                                    object.AngularModel[nProcess - 1] = 1
                                if object.GasMix.Gases[idg].KIN[j] == 2:
                                    object.ScatteringParameter[iEnergy][nProcess - 1] = object.GasMix.Gases[idg].PEInelasticCrossSection[j][iEnergy]
                                    object.AngularModel[nProcess - 1] = 2
                                if iEnergy <= 1:
                                    object.rGas[nProcess - 1] = rGas1
                                    object.ein[nProcess - 1] = object.GasMix.Gases[idg].EIN[j]/rGas1
                                    L = 4
                                    if object.GasMix.Gases[idg].EIN[j] < 0.0:
                                        L = 5
                                    object.ipn[nProcess - 1] = 0
                                    object.InteractionType[nProcess - 1] = L
                                    object.izbr[nProcess - 1] = object.GasMix.Gases[idg].izbr[j]
                                    object.PenningFraction[nProcess - 1][0] = object.GasMix.Gases[idg].PenningFraction[0][j]
                                    object.PenningFraction[nProcess - 1][1] = object.GasMix.Gases[idg].PenningFraction[1][j] * 1e-6 / sqrt(3.0)
                                    object.PenningFraction[nProcess - 1][2] = object.GasMix.Gases[idg].PenningFraction[2][j]
                                    if object.PenningFraction[nProcess - 1][0] > object.avpfrac[0][0]:
                                        object.avpfrac[0][0] = object.PenningFraction[0][nProcess - 1]
                                        object.avpfrac[1][0] = object.PenningFraction[1][nProcess - 1]
                                        object.avpfrac[2][0] = object.PenningFraction[2][nProcess - 1]
                                    if j == object.GasMix.Gases[idg].N_Inelastic:
                                        object.cminexsc[0] *= object.avpfrac[0][0]





