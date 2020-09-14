#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import numpy as np
import matplotlib.pyplot as plt

import Analysis_Code.utility.snapClass as sC
import Analysis_Code.plotting as iPlt
import Analysis_Code.utility.miscellaneous.MWtest as MW
from Analysis_Code.utility.miscellaneous.vectors import *
import Analysis_Code.HaloFinding as HF
import Analysis_Code.utility.miscellaneous.LocalValues as LV
import Analysis_Code.utility.miscellaneous.gasDensityDistribution as gDD
import Analysis_Code.utility.miscellaneous.SmoothingLength as SL

def main(snap, \
         fromGalacticSnap = False, fromHaloCatalogue = True, \
         haloCatalogue = None, haloCoords = (np.array([0,0,0]), 200.), \
         hPar = ((1., 10**3), 5., np.array([0.,0.,0.]), \
                 np.array([50.,50.,50.]), 'Mpc'), \
         npArray = False, PTCAM = ['gas', 'star'] , dWAM = False, \
         PTCCM = ['star'], onlyDisk = False, diskHeight = np.inf, \
         reduceToColdGas = False, Tmax = np.inf, \
         doDist = False , distParameters = (10., list()), \
         doDMDensityProfile = False, dpParameters = ([0., 2.], 0.02, 'Msun'), \
         fitDP = False, DPShapes = ["NFW"], plotLocalDensity = False, \
         localDensity = (0., 0.), plotOverDensity = False, \
         doRotCurve = False, rcPar = (100., 0.8),\
         fromMassDist = True, RCparticleTypes = list(), massAverage = False, \
         addObsData = False, url = '', \
         compareToMassDist = False, testMW = False, \
         doSmoothingLengthPlot = False, doProjectionPlots = False,  \
         ppPar = (100, 15., "inferno"), \
         getLocalDens = False, getLocalVrot = False, \
         lvPar = (8.122, 0.05, 20, 0.05), saveImages = False, \
         saveText = False, \
         testRoutine = False \
         ):
    "Does whatever is specified"
    if fromGalacticSnap == True:
        if fromHaloCatalogue == True:
            data = HF.getHaloFromGalacticSize(snap, haloCatalogue, hPar, \
                                              npArray = npArray, \
                                              PTCAM = PTCAM, PTCCM = PTCCM,\
                                              onlyDisk = onlyDisk, \
                                              diskHeight = diskHeight, \
                                              reduceToColdGas = reduceToColdGas, \
                                              Tmax = Tmax, dWAM = dWAM)
        else:
            data = sC.Snap(snap, npArray = npArray)
            print("finished loading the data!")
            data = HF.zoomIn(data, [0,0, haloCoords[1], haloCoords[0]], \
                             npArray = npArray, PTCAM = PTCAM, PTCCM = PTCCM,\
                             onlyDisk = onlyDisk, diskHeight = diskHeight, \
                             lengthUnit = hPar[4], \
                             reduceToColdGas = reduceToColdGas, Tmax = Tmax, \
                             dWAM = dWAM)
            
    elif testRoutine == True:
        data = HF.testRoutine(snap, PTCAM = PTCAM, PTCCM = PTCCM, \
                              onlyDisk = onlyDisk, diskHeight = diskHeight, \
                              lengthUnit = hPar[4], dWAM = dWAM)
    else:
        data = sC.Snap(snap, npArray = npArray)
        sC.alignToHighestDensityGas(data, npArray = npArray, PTCCM = PTCCM)
        sC.subtractCMWeighted(data, 'Velocities', npArray = npArray)
        sC.alignToNewCS(data, sC.calculateAngularMomentum\
                        (data, PTCAM, densityWeighted = dWAM), \
                            npArray = npArray)
    try:
        if doDist == True:
            Isolated = not fromGalacticSnap
            iPlt.main_dist(data.dmPositions, data.dmVelocities, \
                           *distParameters, saveImages = saveImages, \
                            saveText = saveText, Isolated = Isolated)
    except AttributeError:
        print("There are no DM-particles!")
    
    if doDMDensityProfile == True:
        iPlt.main_dP(data, *dpParameters, fitDP, DPShapes, \
                     plotLocalDensity = plotLocalDensity, \
                     localDensity = localDensity, \
                     plotOverDensity = plotOverDensity, \
                     saveImages = saveImages, saveText = saveText)
            
    if doRotCurve == True:
        if fromMassDist == True:
            iPlt.main_rC(data, *rcPar, saveImages = saveImages)
        else:
            for i in range(len(RCparticleTypes)):
                iPlt.rotCurve(xyArray(eval("data." + RCparticleTypes[i] + 'Positions')), \
                              xyArray(eval("data." + RCparticleTypes[i] + 'Velocities')), \
                              *rcPar, RCparticleTypes[i], i, \
                              massAverage = massAverage, \
                              mass = eval("data." + RCparticleTypes[i] + 'Masses'), \
                              saveText = saveText)
            
            if addObsData == True:
                iPlt.plotDataPointsFromURL(url, label = "Eilers et al.")
                
            if compareToMassDist == True:
                iPlt.plotOnlyEffectiveRotationCurve(data, *rcPar, saveText = saveText)
            plt.legend()
            plt.ylabel(r'$v_{rot}$ [km/s]')
            plt.xlabel('R [kpc]')
            if saveImages == True:
                filename = "Images/RotationCurve.png"
                plt.savefig(filename)
            plt.show()
            plt.close()
        
    if testMW == True:
        MW.MWlike_arrayInput(data, saveText = saveText)
        
    if doSmoothingLengthPlot == True:
        SL.PlotSmoothingLength(data, fromGalacticSnap, saveImages)
    
    if doProjectionPlots == True:
        gDD.smoothColumnDensityPlot(data, *ppPar, saveImages = saveImages)
    
    if getLocalDens == True:
        radii = normarray(data.dmPositions)
        reducedRadii, reducedMasses = LV.cut([radii, data.dmMasses], *lvPar)
        LV.averageOverDifferentBinSizes(reducedRadii, reducedMasses, *lvPar, \
                                        "Density", saveText = saveText)
    if getLocalVrot == True:
        for pT in ['gas', 'star', 'disk']:
            print(pT + ':')
            x = xyArray(eval("data." + pT + "Positions"))
            v = azimuthalArray(xyArray \
                                (eval("data." + pT + "Velocities")), \
                                azimuthalUnitVectors(radialUnitVectors(x)))
            m = eval("data." + pT + "Masses")
            radii = normarray(x)
            redR, redV, masses = LV.cut([radii, v, m], *lvPar)
            LV.averageOverDifferentBinSizes(redR, redV, *lvPar, "Vrot", \
                                            massAverage = massAverage, \
                                            masses = masses, \
                                                saveText = saveText, \
                                                    particleType = pT)
            
