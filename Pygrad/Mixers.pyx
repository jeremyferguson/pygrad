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
            gasData = object.GasMix.Gases[idg]
            object.GasExcitationSlots[idg] = 1
            object.negas[nProcess - 1] = 1
            object.legas[nProcess - 1] = 0
            object.ieshell[nProcess - 1] = 0
            object.collisionFrequency[iEnergy][nProcess - 1] = object.Q[idg][1][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * Object.Beta[iEnergy]
            object.ScatteringParameter[iEnergy][nProcess - 1] = 0.5
            object.AngularCut[iEnergy][nProcess - 1] = 0.0
            if gasData.AngularModel[1] == 1:
                ang.ScatteringParameter1 = gasData.PEElasticCrossSection[1][iEnergy]
                ang.CalcAngCut()
                object.AngularCut[iEnergy][nProcess - 1] = ang.AngCut
                object.ScatteringParameter[iEnergy][nProcess - 1] = ang.ScatteringParameter2
                object.AngularModel[nProcess - 1] = 1
            elif gasData.AngularModel[1] == 2:
                object.ScatteringParameter[iEnergy][nProcess - 1] = gasData.PEElasticCrossSection[1][iEnergy]
                object.AngularModel[nProcess - 1] = 2
            if iEnergy <= 1:
                rGas1 = 1.0 + gasData.E[1]/2.0
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
                object.cminexsc[0] = gasData.E[3] * object.MoleculesPerCm3PerGas[idg]
                object.cminixsc[0] = gasData.E[4] * object.MoleculesPerCm3PerGas[idg]
                object.ecloss[0] = gasData.E[2]
                object.WPLN[0] = gasData.E[5]
            if object.FinalElectronEnergy >= gasData.E[2]:
                if gasData.N_Ionization <= 1:
                    nProcess += 1
                    object.GasExcitationSlots[idg] = nProcess
                    if object.icount == 1:
                        QIndex = 4
                        object.DOUBLE[0][iEnergy] = gasData.Q[2][iEnergy] / gasData.Q[4][iEnergy] - 1.0
                    else:
                        QIndex = 2
                    object.CollisionFrequency[iEnergy][nProcess - 1] = gasData.Q[QIndex][iEnergy] * object.VMoleculesPerCm3PerGas * object.Beta[iEnergy]
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
                    if gasData.AngularModel[pIndex] == 1:
                        Ang.ScatteringParameter1 = gasData.PEElasticCrossSection[pIndex][iEnergy]
                        Ang.CalcAngCut()
                        object.AngleCut[iEnergy][nProcess - 1] = Ang.AngCut
                        object.ScatteringParameter = Ang.ScatteringParameter2
                        object.AngularModel[nProcess - 1] = 1
                    assignWPL(nProcess,0)
                    if iEnergy == 1:
                        object.rGas[nProcess - 1] = rgas1
                        object.ein[nProcess - 1] = gasData
                        object.ipn[nProcess - 1] = 1
                        L = 2
                        object.InteractionType[nProcess -1] = L
                        object.izbr[nProcess - 1] = 0
                        object.PenningFraction[nProcess - 1][0] = 0.0
                        object.PenningFraction[nProcess - 1][1] = 0.0
                        object.PenningFraction[nProcess - 1][2] = 0.0
                else:
                    for kion in range(gasData.N_Ionization):
                        nProcess += 1
                        object.GasExcitationSlots[idg] = nProcess
                        object.CollisionFrequency[iEnergy,nProcess - 1] = gasData.IonizationCrossSection[kin][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * object.Beta[iEnergy]
                        object.IonCollisionFreq[iEnergy] += object.CollisionFrequency[iEnergy,nProcess - 1]
                        object.ScatteringParameter[iEnergy,nProcess - 1] = 0.5
                        object.AngleCut[iEnergy,nProcess - 1] = 1.0
                        object.AngularModel[nProcess - 1] = 0
                        object.negas[nProcess - 1] = 1
                        if gasData.AngularModel[2] == 1:
                            Ang.ScatteringParameter1[iEnergy][nProcess - 1] = gasData.PEInelasticCrossSection[kion][iEnergy] 
                            Ang.CalcAngCut()
                            object.AngleCut[iEnergy][nProcess - 1] = Ang.AngCut
                            object.ParameterScattering[iEnergy][nProcess - 1] = Ang.ScatteringParameter2
                            object.AngularModel[nProcess - 1] = 1
                        if gasData.AngularModel[2] == 2:
                            object.ParameterScattering[iEnergy][nProcess - 1] = gasData.PEInelasticCrossSection[kion][iEnergy]
                            object.AngularModel[nProcess - 1] = 2
                        assignWPL(nProcess, kion)
                        if iEnergy <= 1:
                            object.rGas[nProcess - 1] = rGas1
                            object.ein[nProcess - 1] = gasData.IonizationEnergy[kion]/rGas1
                            object.ipn[nProcess - 1] = 1
                            L = 2
                            object.InteractionType[nProcess - 1] = L
                            object.izbr[nProcess - 1] = 0
                            object.PenningFraction[nProcess - 1][0] = 0.0
                            object.PenningFraction[nProcess - 1][1] = 0.0
                            object.PenningFraction[nProcess - 1][2] = 0.0
            loadAttachmentData(gasData,idg,nProcess,iEnergy)
            loadInelasticData(gasData,idg,nProcess,iEnergy)

cpdef loadAttachmentData(gasData, idg, nProcess, iEnergy):
    if object.FinalElectronEnergy >= gasData.E[3]:
        if gasData.N_Attachment <= 1:
            nProcess += 1
            object.GasExcitationSlots[idg] = nProcess
            object.CollisionFrequency[iEnergy][nProcess - 1] = gasData.Q[3][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * object.Beta[iEnergy]
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
            for jj in range(gasData.N_Attachment):
                nProcess += 1
                object.GasExcitationSlots[idg] = nProcess
                object.CollisionFrequency[iEnergy][nProcess - 1] = gasData.Q[3][iEnergy] * object.VMoleculesPerCm3PerGas[idg] * object.Beta[iEnergy]
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

cpdef loadInelastic(gasData, idg, nProcess, iEnergy):
    if gasData.N_Inelastic != 0:
        for j in range(gasData.N_Inelastic):
            nProcess += 1
            object.GasExcitationSlots[idg] = nProcess
            object.negas[nProcess - 1] = 1
            object.legas[nProcess - 1] = 0
            object.ieshell[nProcess - 1] = 0
            object.CollisionFrequency[iEnergy][nProcess - 1] = object.InelasticCrossSection[j][iEnergy] * object.VMoleculesPerCm3PerGas[0] * object.Beta[iEnergy]
            object.ScatteringParameter[iEnergy][nProcess] = 0.5
            object.AngleCut[iEnergy][nProcess] = 1.0
            object.AngularModel[nProcess - 1] = 0
            if gasData.KIN[j] == 1:
                Ang.ScatteringParameter1 = gasData.PEInelasticCrossSection[j][iEnergy]
                Ang.CalcAngCut()
                object.ScatteringParameter[iEnergy][nProcess - 1] = Ang.ScatteringParameter2
                object.AngleCut[iEnergy][nProcess - 1] = Ang.AngleCut
                object.AngularModel[nProcess - 1] = 1
            if gasData.KIN[j] == 2:
                object.ScatteringParameter[iEnergy][nProcess - 1] = gasData.PEInelasticCrossSection[j][iEnergy]
                object.AngularModel[nProcess - 1] = 2
            if iEnergy <= 1:
                object.rGas[nProcess - 1] = rGas1
                object.ein[nProcess - 1] = gasData.EIN[j]/rGas1
                L = 4
                if gasData.EIN[j] < 0.0:
                    L = 5
                object.ipn[nProcess - 1] = 0
                object.InteractionType[nProcess - 1] = L
                object.izbr[nProcess - 1] = gasData.izbr[j]
                object.PenningFraction[nProcess - 1][0] = gasData.PenningFraction[0][j]
                object.PenningFraction[nProcess - 1][1] = gasData.PenningFraction[1][j] * 1e-6 / sqrt(3.0)
                object.PenningFraction[nProcess - 1][2] = gasData.PenningFraction[2][j]
                if object.PenningFraction[nProcess - 1][0] > object.avpfrac[0][0]:
                    object.avpfrac[0][0] = object.PenningFraction[0][nProcess - 1]
                    object.avpfrac[1][0] = object.PenningFraction[1][nProcess - 1]
                    object.avpfrac[2][0] = object.PenningFraction[2][nProcess - 1]
                if j == gasData.N_Inelastic:
                    object.cminexsc[0] *= object.avpfrac[0][0]


cpdef assignWPL(index, gIndex):
    object.wpl[index - 1] = gasData.EB1[gIndex]
    object.nc0[index - 1] = gasData.NC0[gIndex]
    object.ec0[index - 1] = gasData.EC0[gIndex]
    object.ng1[index - 1] = gasData.NG1[gIndex]
    object.eg1[index - 1] = gasData.EG1[gIndex]
    object.ng2[index - 1] = gasData.NG2[gIndex]
    object.eg2[index - 1] = gasData.EG2[gIndex]
    object.wklm[index - 1] = gasData.WK[gIndex]
    object.efl[index - 1] = gasData.EFL[gIndex]