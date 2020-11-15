from Pygrad.Pygrad cimport Pygrad
from Pygrad cimport PygradException
from libc.math cimport acos, sqrt
cimport numpy as np
import  numpy as np
import cython

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef Fldist(Pygrad object):
    cdef float eph
    for i in range(512):
        eph = object.efl[i]
        if eph != 0.0:
            idum,kdum,ldum,dist = calcAbsorption(object,3,eph)
            object.efl[i] = dist

cpdef calcAbsorption(Pygrad object, int jf, float eph):
    cdef float xsec[306], xsecc[18], xsecr[18], xsecp[18], absl[306], abslc[18], abslr[18], abslp[18], xsum[360],dummy,dist
    cdef int ishell, kgas, lgas, lcflag, lrflag, lpflag, lpeflag, ephlg, ipt,ipet, ID
    ishell = 0
    kgas = 0
    lgas = 0
    dist = 0.0
    lcflag = 0
    lrflag = 0
    lpflag = 0
    lpeflag = 0
    ephlg = np.log(eph)
    ipt = -1
    dummy = 0.0
    #calculate photoelectric cross section
    for i in range(object.NumberOfGases):
        for j1 in range(3):
            for j in range(17):
                ipt += 1
                xsec[ipt] = 0.0
                absl[ipt] = 0.0
                if j > object.ishlmx[i,j1]:
                    continue
                if ephlg <= object.xpe[i][j1][j][0]:
                    continue
                for k in range(1,60):
                    if ephlg <= object.xpe[i][j1][j][k]:
                        A = (object.ype[i][j1][j][k] - object.ype[i][j1][j][k-1])/(object.xpe[i][j1][j][k-1] - object.xpe[i][j1][j][k-1])
                        B = (object.xpe[i][j1][j][k-1] * object.ype[i][j1][j][k] - object.xpe[i][j1][j][k] * object.ype[i][j1][j][k-1])/(object.xpe[i][j1][j][k-1] - object.xpe[i][j1][j][k])
                        xsec[ipt] = np.exp(A*ephlg+B)
                        absl[ipt] = xsec[ipt]*object.MoleculesPerCm3PerGas[i]
                        break
    ipt = -1
    #calculate Compton cross section and find absorption length
    for i in range(object.NumberOfGases):
        for j1 in range(3):
            ipt += 1
            xsecc[ipt] = 0.0
            abslc[ipt] = 0.0
            if jf == 3 or jf == 2:
                continue
            if object.lcmp != 1:
                continue
            if ephlg < object.xen[i,j1,0]:
                continue
            for k in range(1,54):
                if ephlg <= object.xen[i][j1][k]:
                    A = (object.ycp[i][j1][k] - object.ycp[i][j1][k])/(object.xen[i][j1][k] * object.ycp[i][j1][k-1])
                    B = (object.xen[i][j1][k-1] * object.ycp[i][j1][k] - object.xen[i][j1][k]*object.ycp[i][j1][k-1])/(object.xen[i][j1][k-1] - object.xen[i][j1][k])
                    xsecc[ipt] = np.exp(A*ephlg+B)
                    abslc[ipt] = xsecc[ipt]*object.MoleculesPerCm3PerGas[i]
                    break

    ipt = -1
    #calculate Rayleigh cross section and find absorption length
    for i in range(object.NumberOfGases):
        for j1 in range(3):
            ipt += 1
            xsecr[ipt] = 0.0
            abslr[ipt] = 0.0
            if jf == 3 or jf == 2:
                continue
            if object.lray != 1:
                continue
            if ephlg < object.xen[i][j1][0]:
                continue
            for k in range(1,54):
                if ephlg <= object.xen[i][j1][k]:
                    A = (object.yry[i][j1][k]-object.yry[i][j1][k-1])/(object.xen[i][j1][k]-object.xen[i][j1][k-1])
                    B = (object.xen[i][j1][k-1]*object.yry[i][j1][k]-object.xen[i][j1][k]*object.yry[i][j1][k-1])/(object.xen[i][j1][k-1]-object.xen[i][j1][k])
                    xsecr[ipt] = np.exp(A*ephlg+B)
                    abslr[ipt] = xsecr[ipt]*object.MoleculesPerCm3PerGas[i]
                    break
    ipt = -1
    #calculate pair production cross section and find absorption length
    for i in range(object.NumberOfGases):
        for j1 in range(3):
            ipt += 1
            xsecp[ipt] = 0.0
            abslp[ipt] = 0.0
            if jf == 3 or jf == 2:
                continue
            if object.lpap != 1:
                continue
            if ephlg < object.xen[i][j1][k]:
                A = (object.ypp[i][j1][k]-object.ypp[i][j1][k-1])/(object.xen[i][j1][k]-object.xen[i][j1][k-1])
                B = (object.xen[i][j1][k-1]*object.ypp[i][j1][k]-object.xen[i][j1][k]*object.ypp[i][j1][k-1])/(object.xen[i][j1][k-1]-object.xen[i][j1][k])
                xsecp[ipt] = np.exp(A*ephlg+B)
                abslp[ipt] = xsecp[ipt]*object.MoleculesPerCm3PerGas[i]
                break
    ifin = object.NumberOfGases * 17 * 3
    for j in range(1,ifin):
        xsec[j] += xsec[j-1]
        absl[j] += absl[j-1]
    ifinr = object.NumberOfGases * 3
    for j in range(1,ifinr):
        xsecc[j] += xsecc[j-1]
        abslc[j] += abslc[j-1]
        xsecr[j] += xsecr[j-1]
        abslr[j] += abslr[j-1]
        xsecp[j] += xsecp[j-1]
        abslp[j] += abslp[j-1]
    totalXSec = xsec[ifin-1] + xsecc[ifinr-1]+xsecr[ifinr-1]+xsecp[ifinr-1]
    totalAbsLength = absl[ifin-1]+abslr[ifinr-1]+abslc[ifinr-1]+abslp[ifinr-1]
    if jf == 3:
        dist = 1.0/(totalAbsLength*100.0)
        return (ishell,kgas,lgas,dist)
    if jf == -1:
        if totalAbsLength > 0.0:
            object.absXray = 1.0e4 / totalAbsLength
        elif totalAbsLength == 0.0:
            object.absXray = 1.0e15
        return (ishell,kgas,lgas,dist) 
    if totalAbsLength == 0.0:
        ishell = -1
        return (ishell,kgas,lgas,dist)
    #TODO:normalizing all these values, probably a numpy way to do this
    for j in range(ifin):
        xsec[j] /= totalXSec
    for j in range(ifinr):
        xsecc[j] /= totalXSec
        xsecr[j] /= totalXSec
        xsecp[j] /= totalXSec
    for j in range(ifin):
        xsum[j] = xsec[j]
    iend = ifin
    if object.lcmp == 1:
        istart = ifin
        iend = ifin + ifinr
        for j in range(istart,iend):
            xsum[j] = xsum[istart-1] + xsecc[j-istart+1]
    if object.lray == 1:
        if object.lcmp == 0:
            istart = ifin
            iend = ifin+ifinr
        elif object.lcmp == 1:
            istart = ifin + ifinr
            iend = ifin + 2 * ifinr
        for j in range(istart,iend):
            xsum[j] = xsum[istart-1]+xsecr[j-istart+1]
    if object.lpap == 1:
        flagsum = object.lcmp + object.lray
        if flagsum <= 2 and flagsum >= 0:
            istart = ifin + flagsum * ifinr
            iend = ifin + (flagsum + 1) * ifinr
        else:
            raise PygradException('Invalid flag value')
        for j in range(istart,iend):
            xsum[j] = xsum[istart-1] + xsecp[j-istart+1]
    R1 = Pygrad.drand48(dummy)
    for j in range(iend):
        if xsum[j] < R1:
            ID = j
            break
    ipet = object.NumberOfGases*3*17
    #Photoelectric flag setting
    if ID <= ipet:
        lpeflag = 1
        if ID <= 51:
            kgas = 1
            lgas = int(np.ceil(ID/17.0))
            ishell = ID % 17
            if not ishell and ID:
                ishell = 17
        elif ID <= 102:
            kgas = 2
            lgas = int(np.ceil((ID-51)/17.0))
            ishell = ID % 17
            if not ishell:
                ishell = 17
        elif ID <= 153:
            kgas = 3
            lgas = int(np.ceil((ID-102)/17.0))
            ishell = ID % 17
            if not ishell:
                ishell = 17
        elif ID <= 204:
            kgas = 4
            lgas = int(np.ceil((ID-153)/17.0))
            ishell = ID % 17
            if not ishell:
                ishell = 17
        elif ID <= 255:
            kgas = 5
            lgas = int(np.ceil((ID-204)/17.0))
            ishell = ID % 17
            if not ishell:
                ishell = 17
        else:
            kgas = 6
            lgas = int(np.ceil((ID-255)/17.0))
            ishell = ID % 17
            if not ishell:
                ishell = 17
    #Compton, Rayleigh, or Pair Production
    else:
        if ID <= ipet + ifinr:
            if object.lcmp == 1:
                lcflag = 1
            if object.lcmp == 0 and object.lray == 1:
                lrflag = 1
            if object.lcmp == 0 and object.lray == 0:
                lpflag = 1
            if ID <= ipet + 3:
                kgas = 1
                lgas = ID -ipet
            elif ID <= ipet+6:
                kgas = 2
                lgas = ID - ipet - 3
            elif ID <= ipet + 9:
                kgas = 3
                lgas = ID - ipet - 6
            elif ID <= ipet+12:
                kgas = 4
                lgas = ID - ipet - 9
            elif ID <= ipet+15:
                kgas = 5
                lgas = ID - ipet - 12
            else:
                kgas = 6
                lgas = ID - ipet - 15
        elif ID <= ipet + 2*ifinr:
            if object.lray == 1:
                lrflag = 1
            if object.lray == 0 and object.lpap == 1:
                lpflag = 1
            if ID <= ipet + ifinr + 3:
                kgas = 1
                lgas = ID - ipet - ifinr
            elif ID <= ipet + ifinr + 6:
                kgas = 2
                lgas = ID - ipet - ifinr - 3
            elif ID <= ipet + ifinr + 9:
                kgas = 3
                lgas = ID - ipet - ifinr - 6
            elif ID <= ipet + ifinr + 12:
                kgas = 4
                lgas = ID - ipet - ifinr - 9
            elif ID <= ipet + ifinr + 15:
                kgas = 5
                lgas = ID - ipet - ifinr - 12
            else:
                kgas = 6
                lgas = ID - ipet - ifinr - 15
        else:
            lpflag = 1
            if ID <= ipet + 3*ifinr:
                kgas = 1
                lgas = ID - ipet - ifinr - ifinr
            elif ID <= ipet + ifinr + ifinr + 6:
                kgas = 2
                lgas = ID - ipet - ifinr - ifinr - 3
            elif ID <= ipet + ifinr + ifinr + 9:
                kgas = 3
                lgas = ID - ipet - ifinr - ifinr - 6
            elif ID <= ipet + ifinr + ifinr + 12:
                kgas = 4
                lgas = ID - ipet - ifinr - ifinr - 9
            elif ID <= ipet + ifinr + ifinr + 15:
                kgas = 5
                lgas = ID - ipet - ifinr - ifinr - 12
            else:
                kgas = 6
                lgas = ID - ipet - ifinr - ifinr - 15
    if ID > ipet + 54:
        raise PygradException('identifier in absorption calculation is greater than the limit: {}\nProgram stopped.')
    R1 = Pygrad.drand48(dummy)
    dist = -np.log(R1)/(totalAbsLength * 100.0)
    return (ishell, kgas, lgas, dist)

