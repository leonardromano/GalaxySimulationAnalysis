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

def averageOverDifferentBinSizes(radii, data, parameters, \
                                 dataType, massAverage = False, \
                                 masses = None, saveText = False, \
                                 particleType = ""):
    "returns the average of the data in the specified region"
    values = list()
    uncertainties =list()
    for i in range(parameters[2]):
        temp = list()
        if massAverage == True:
            temp2 = list()
        dx=parameters[1]+i*parameters[3]
        j=0
        for r in radii:
            if parameters[0]-dx/2.<r<parameters[0]+dx/2:
                temp.append(data[j])
                if massAverage == True:
                    temp2.append(masses[j])
            j+=1
        if dataType == "Density":
            values.append(calculateDensity(temp, dx/2, parameters[0]))
            uncertainties.append(1.)
        elif dataType == 'Vrot':
            values.append(mean(temp, dx, \
                               massAverage = massAverage, masses = temp2))
            uncertainties.append(uncertainty(temp, values[i], \
                                             massAverage = massAverage, \
                                             masses = temp2))
    w = weights(uncertainties)
    if dataType == "Density":
        m = weightedMean(values, w)
        dm = uncertainty(values, m)
        print 'rho(' + str(parameters[0]) + ' kpc) = ' + str(m) + ' +- ' + \
        str(dm) + ' GeV/cm^3'
        if saveText == True:
            np.savetxt("PlotData/LocalDensity", np.asarray([m, dm]), \
                       delimiter = ' +- ')
        return (m, dm)
    elif dataType == 'Vrot':
        m = weightedMean(values, w)
        dm = weightedUncertainty(w)
        print 'vr(' + str(parameters[0]) + ' kpc) = ' + str(m) + ' +- ' + \
        str(dm) + ' km/s'
        if saveText == True:
            np.savetxt("PlotData/LocalRotationVelocity_" + particleType, \
                       np.asarray([m, dm]), delimiter = ' +- ')
        return (m, dm)
        

def calculateDensity(Masses, dx, x):
    "calculates the local Density from the mass data"
    M = 0
    for m in Masses:
        M+=m
    if M == 0:
        print "dx>! " + str(dx)
    return (3*M/(4*np.pi*((x+dx)**3-(x-dx)**3)))*(10**-14)*Msun*(cvak**2)/(qe*pc**3)

def mean(data, dx, massAverage = False, masses = None):
    "calculates the mean of a list of data points"
    N=0.
    M=0
    i=0
    for point in data:
        if massAverage == False:
            N+=1.
            M+=point
        else:
            N+=masses[i]
            M+=masses[i]*point
            i+=1
    if N!=0.:
        return M/N
    else:
        print "dx>! " + str(dx)
        return 0
    
def weights(dx):
    "calculates the weights"
    w = list()
    for dxi in dx:
        w.append(1./dxi**2)
    return w

def weightedMean(x, w):
    "calculates the weighted mean"
    N = 0
    Z = 0
    i = 0
    for wi in w:
        Z+= wi*x[i]
        N+= wi
        i+=1
    if N != 0:
        return Z/N
    else:
        return 0

def uncertainty(data, mean, massAverage = False, masses = None):
    "Calculates the standard deviation of the mean"
    s = 0.
    N = 0.
    i=0
    for point in data:
        if massAverage == False:
            s+= (point-mean)**2
            N+=1.
        else:
            s+=masses[i]*(point-mean)**2
            N+=masses[i]
            i+=1
    if N!=0 and N!=1:
        if massAverage == True:
            return np.sqrt(s/(N*i))
        else:
            return np.sqrt(s/(N*(N-1)))
    else:
        return 1.

def weightedUncertainty(weights):
    "calculates the weighted uncertainty"
    s = 0
    for wi in weights:
        s+=wi
    return 1./np.sqrt(s)