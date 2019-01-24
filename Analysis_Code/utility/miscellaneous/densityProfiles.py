#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import numpy as np

cvak = 299792458.0
qe = 1.6021766208*10**(-19)
pc = 30856775814913700.0
Msun = 1.98892*10**30

def getXYArray(XArray, YArray):
    "Turns two (n,1) Arrays of into one combined (n,2) array"
    i = 0
    liste = list()
    for X in XArray:
        liste.append([X, YArray[i]])
        i+=1
    return np.asarray(liste)

def createDensityProfileHistogramm(masses, radii, logrmin, logrmax, dlogr, \
                                   unit, saveText = False):
    "Returns density profile from mass and position data"
    logR = np.arange(logrmin, logrmax, dlogr)
    Radii = np.exp(np.log(10)*logR)
    length = Radii.shape[0]
    density = list()
    xValues = list()
    for i in range(length-1):
        j=0
        M=0
        for radius in radii:
            if Radii[i] <= radius < Radii[i+1]:
                M+=masses[j]
            j+=1
        rho = getDensity(M, Radii[i], Radii[i+1], unit)
        density.append(rho)
        xValues.append(np.log10(Radii[i]))
        density.append(rho)
        xValues.append(np.log10(Radii[i+1]))
    densityArray = np.asarray(density)
    radiusArray = np.asarray(xValues)
    logrho = np.log10(densityArray)
    if saveText == True:
        np.savetxt('PlotData/densityProfile', \
                   getXYArray(np.exp(radiusArray*np.log(10)), densityArray), delimiter = ' ')
    return getXYArray(radiusArray, logrho)
    

def getDensity(mass, rmin, rmax, unit):
    "Calculates the local mass density in M_sun pc^(-3)"
    alpha = getConversionFactor(unit)
    return (mass/(4*np.pi*(rmax**3 - rmin**3)/3))*alpha

def getConversionFactor(unit):
    "Returns the conversion factor of the units used"
    if unit == 'Msun':
        return 10.
    elif unit == 'GeV':
        return (10**-14)*Msun*(cvak**2)/(qe*pc**3)
    else:
        return 1.
    
def convertToUnits(logData, unit):
    "Converts the LogData to different unit"
    if unit == 'Msun':
        return logData
    else:
        alpha = getConversionFactor(unit)/10.
        temp = list()
        for logn in logData:
            n = alpha*10**(logn)
            temp.append(np.log10(n))
        return np.asarray(temp)

def getNFWMasses(radii, mS):
    "calculates the theoretical masses, predicted by the \
    NFW profile at radii in units of rS"
    masses = list()
    for r in radii:
        masses.append([12*mS*(np.log(1+r)-r/(1+r))])
    return np.asarray(masses)

def getYLabel(unit):
    "returns the ylabel of the units used"
    if unit == "Msun":
        return r'$log_{10}\left(\rho \left[M_{\odot}pc^{-3}\right]\right)$'
    elif unit == 'GeV':
        return r'$log_{10}\left(\rho \left[GeVcm^{-3}\right]\right)$'
    else:
        return 'A.U.'
    