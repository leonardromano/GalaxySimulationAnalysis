#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import snapClass as clap
import plotting as iPlt
import MWtest as MW
from utility.vectors import *
import matplotlib.pyplot as plt
import HaloFinding as HF
import numpy as np
import utility.LocalValues as LV

def main(snap, \
         fromGalacticSnap = False, fromHaloCatalogue = True, \
         haloCatalogue = None, haloCoords = (np.array([0,0,0]), 200.), \
         hPar = ((1., 10**3), 5., np.array([0.,0.,0.]), \
                 np.array([50.,50.,50.]), 'Mpc'), \
         npArray = False, PTCAM = ['gas', 'star'] , PTCCM = ['star'], \
         onlyDisk = False, diskHeight = np.inf, \
         reduceToColdGas = False, Tmax = np.inf, \
         doDist = False , distParameters = (10., list()), \
         doRotCurve = False, rcPar = (100., 0.8, list()),\
         fromMassDist = True, addObsData = False, url = '', \
         compareToNFW = False, rNFWPar = (20.55, 0.008), \
         compareToMassDist = False, testMW = False, \
         doDMDensityProfile = False, dpParameters = ([0., 2.], 0.02, 'Msun'), \
         fitNFW = False, NFWParameters = (0.6, -1.2), plotLocalDensity = False, \
         localDensity = (0., 0.), plotOverDensity = False, \
         fitDMDensityProfile = False, DMDensityFits = list(), \
         DMDensityFitDomains = [[0.,1.], [1.,2.]], \
         doProjectionPlots = False,  \
         ppPar = (500, [[-50, 50], [-50, 50]], [[-50, 50], [-15, 15]]), \
         getLocalDens = False, getLocalVrot = False, \
         lvPar = (8.122, 0.05, 20, 0.05), testRoutine = False\
         ):
    "I do what you want me to do if specified"
    if fromGalacticSnap == True:
        if fromHaloCatalogue == True:
            data = HF.getHaloFromGalacticSize(snap, haloCatalogue, hPar, \
                                              npArray = npArray, \
                                              PTCAM = PTCAM, PTCCM = PTCCM,\
                                              onlyDisk = onlyDisk, \
                                              diskHeight = diskHeight, \
                                              reduceToColdGas = reduceToColdGas, \
                                              Tmax = Tmax)
        else:
            data = clap.Snap(snap, npArray = npArray)
            print "finished loading the data!"
            data = HF.zoomIn(data, [0,0, haloCoords[1], haloCoords[0]], \
                             npArray = npArray, PTCAM = PTCAM, PTCCM = PTCCM,\
                             onlyDisk = onlyDisk, diskHeight = diskHeight, \
                             lengthUnit = hPar[4], \
                             reduceToColdGas = reduceToColdGas, Tmax = Tmax)
            
    elif testRoutine == True:
        data = HF.testRoutine(snap, PTCAM = PTCAM, PTCCM = PTCCM, \
                              onlyDisk = onlyDisk, diskHeight = diskHeight, \
                              lengthUnit = hPar[4])
    else:
        data = clap.Snap(snap, npArray = npArray)
        clap.alignToHighestDensityGas(data, npArray = npArray, PTCCM = PTCCM)
        clap.subtractCMWeighted(data, 'Velocities', npArray = npArray)
        clap.alignToNewCS(data, clap.calculateAngularMomentum(data, PTCAM), \
                        npArray = npArray)
    try:
        if doDist == True:
            iPlt.main_dist(data.dmPositions, data.dmVelocities, \
                           distParameters[0], distParameters[1])
    except AttributeError:
        print "There are no DM-particles!"
            
    if doRotCurve == True:
        if fromMassDist == True:
            iPlt.main_rC(data, rcPar[0], rcPar[1])
        else:
            i=0
            for pT in rcPar[2]:
                iPlt.rotCurve(xyArray(eval("data." + pT + 'Positions')), \
                              xyArray(eval("data." + pT + 'Velocities')), \
                              rcPar[0], rcPar[1], pT, i)
                i += 1
            
            if addObsData == True:
                iPlt.plotDataPointsFromURL(url, label = "Eilers et al.")
            
            if compareToNFW == True:
                iPlt.plotNFWRotationCurve(rcPar[0], rcPar[1], \
                                          rNFWPar[0], rNFWPar[1])                
            if compareToMassDist == True:
                iPlt.plotOnlyEffectiveRotationCurve(data, rcPar[0], rcPar[1])
            plt.legend()
            plt.ylabel(r'$v_{rot}$ [km/s]')
            plt.xlabel('R [kpc]')
            plt.show()            
        
    if testMW == True:
        MW.MWlike_arrayInput(data)
        
    if doDMDensityProfile == True:
        iPlt.main_dP(data, dpParameters[0], dpParameters[1], \
                     fitNFW, NFWParameters, \
                     fitDMDensityProfile, DMDensityFits, \
                     DMDensityFitDomains, dpParameters[2], \
                     plotLocalDensity = plotLocalDensity, \
                     localDensity = localDensity, \
                     plotOverDensity = plotOverDensity)
    
    if doProjectionPlots == True:
        iPlt.histogram2D \
            (data.gasPositions[ : , 0 ], data.gasPositions[ : , 1 ], \
             [r'x $\left[kpc \right]$', r'y $\left[kpc \right]$'], \
             ppPar[1], 'gas', data.gasDensity, \
             colorbarLabel = 'logdensity', numberOfBins = ppPar[0])
        iPlt.histogram2D \
            (data.gasPositions[ : , 0 ], data.gasPositions[ : , 2 ], \
             [r'x $\left[kpc \right]$', r'z $\left[kpc \right]$'], \
             ppPar[2], 'gas', data.gasDensity, \
             colorbarLabel = 'logdensity', numberOfBins = ppPar[0])
    
    if getLocalDens == True:
        radii = normarray(data.dmPositions)
        reducedRadii = list()
        reducedMasses = list()
        i=0
        for r in radii:
            if lvPar[0]-(lvPar[1]+(lvPar[2]-1)*lvPar[3])/2 < r \
            < lvPar[0]+(lvPar[1]+(lvPar[2]-1)*lvPar[3])/2:
                reducedRadii.append(r)
                reducedMasses.append(data.dmMasses[i])
            i+=1
        LV.averageOverDifferentBinSizes(np.asarray(reducedRadii), \
                                        np.asarray(reducedMasses), lvPar, "Density")
    if getLocalVrot == True:
        for pT in ['gas', 'star', 'disk']:
            print pT + ':'
            x = xyArray(eval("data." + pT + 'Positions'))
            v = azimuthalArray(xyArray \
                                (eval("data." + pT + 'Velocities')), \
                                azimuthalUnitVectors(radialUnitVectors(x)))
            radii = normarray(x)
            redR = list()
            redV = list()
            i=0
            for r in radii:
                if lvPar[0]-(lvPar[1]+(lvPar[2]-1)*lvPar[3])/2 < r \
            < lvPar[0]+(lvPar[1]+(lvPar[2]-1)*lvPar[3])/2:
                    redR.append(r)
                    redV.append(v[i])
                i+=1
            LV.averageOverDifferentBinSizes(np.asarray(redR), \
                                        np.asarray(redV), lvPar, "Vrot")
            