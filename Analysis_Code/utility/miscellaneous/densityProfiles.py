#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import numpy as np
from Analysis_Code.utility.miscellaneous.vectors import *

cvak = 299792458.0
qe = 1.6021766208*10**(-19)
pc = 30856775814913700.0
Msun = 1.98892*10**30

def createDensityProfileHistogramm(masses, radii, logrmin, logrmax, dlogr, \
                                   unit, fitDP = True, saveText = False):
    "Returns density profile from mass and position data"
    logR = np.arange(logrmin, logrmax, dlogr)
    Radii = 10**logR
    density = list()
    if fitDP == True:
        errors = list()
    xValues = list()
    for i in range(Radii.shape[0]-1):
        M=0
        for j in range(radii.size):
            if Radii[i] <= radii[j] < Radii[i+1]:
                M+=masses[j]
        rho = getDensity(M, Radii[i], Radii[i+1], unit)
        density.append(rho)
        xValues.append(Radii[i])
        density.append(rho)
        xValues.append(Radii[i+1])
        if fitDP == True:
            try:
                #implicitly assumes that all DM particles have the same mass
                srho = rho*np.sqrt(masses[0]/M)
            except ZeroDivisionError:
                srho = 10**10
            errors.append(srho)
            errors.append(srho)
    densityArray = np.asarray(density)
    radiusArray = np.asarray(xValues)
    if fitDP == True:
        errorArray = np.asarray(errors)
        if saveText == True:
            np.savetxt('PlotData/densityProfile', \
                       getNDArray([radiusArray, densityArray, errorArray]), \
                           delimiter = ' ')
        return radiusArray, densityArray, errorArray
    else:
        if saveText == True:
            np.savetxt('PlotData/densityProfile', \
                       getNDArray([radiusArray, densityArray]), \
                           delimiter = ' ')
        return radiusArray, densityArray
    

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
        return r'$\rho \left[M_{\odot}pc^{-3}\right]$'
    elif unit == 'GeV':
        return r'$\rho \left[GeVcm^{-3}\right]$'
    else:
        return 'A.U.'
    