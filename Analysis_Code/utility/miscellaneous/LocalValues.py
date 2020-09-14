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

def averageOverDifferentBinSizes(radii, data, Rs, dr, N, ddr, \
                                 dataType, massAverage = False, \
                                 masses = None, saveText = False, \
                                 particleType = ""):
    "returns the average of the data in the specified region"
    values = np.zeros(N)
    uncertainties =np.zeros(N)
    for i in range(N):
        temp = list()
        if massAverage == True:
            temp2 = list()
        dx = dr+i*ddr
        for j in range(radii.size):
            if Rs-dx/2.<radii[j]<Rs+dx/2:
                temp.append(data[j])
                if massAverage == True:
                    temp2.append(masses[j])
        if dataType == "Density":
            values[i] += calculateDensity(temp, dx/2, Rs)
            uncertainties[i] += 1.
        elif dataType == 'Vrot':
            vdata = np.asarray(temp)
            mdata = np.asarray(temp2)
            values[i]+=mean(vdata, massAverage = massAverage, masses = mdata)
            uncertainties[i]+=uncertainty(vdata, values[i], \
                                       massAverage = massAverage, \
                                           masses = mdata)
    w = weights(uncertainties)
    if dataType == "Density":
        m = weightedMean(values, w)
        dm = uncertainty(values, m)
        print('rho(' + str(Rs) + ' kpc) = ' + str(m) + ' +- ' + \
        str(dm) + ' GeV/cm^3')
        if saveText == True:
            np.savetxt("PlotData/LocalDensity", np.asarray([m, dm]), \
                       delimiter = ' +- ')
        return (m, dm)
    elif dataType == 'Vrot':
        m = weightedMean(values, w)
        dm = weightedUncertainty(w)
        print('vr(' + str(Rs) + ' kpc) = ' + str(m) + ' +- ' + \
        str(dm) + ' km/s')
        if saveText == True:
            np.savetxt("PlotData/LocalRotationVelocity_" + particleType, \
                       np.asarray([m, dm]), delimiter = ' +- ')
        return (m, dm)

def cut(arrays, Rs, drmin, N, ddr):
    reducedList = list()
    for array in arrays:
        reducedList.append([])
    for i in range(arrays[0].size):
        if Rs-(drmin+(N-1)*ddr)/2 < arrays[0][i] < Rs+(drmin+(N-1)*ddr)/2:
            for j in range(len(arrays)):
                reducedList[j].append(arrays[j][i])
    return np.array(reducedList)

def calculateDensity(Masses, dx, x):
    "calculates the local Density from the mass data"
    M = 0
    for m in Masses:
        M+=m
    if M == 0:
        print("dx>! " + str(dx))
    return (3*M/(4*np.pi*((x+dx)**3-(x-dx)**3)))*(10**-14)*Msun*(cvak**2)/(qe*pc**3)

def mean(data, massAverage = False, masses = None):
    "calculates the mean of a list of data points"
    if massAverage == False:
        return weightedMean(data, np.ones(data.size))
    else:
        return weightedMean(data, masses)
    
def weights(dx):
    "calculates the weights"
    w = np.zeros(dx.size)
    for i in range(dx.size):
        w[i] += 1/dx[i]**2
    return w

def weightedMean(x, w):
    "calculates the weighted mean"
    N = 0
    Z = 0
    for i in range(w.size):
        Z+= w[i]*x[i]
        N+= w[i]
    if N != 0:
        return Z/N
    else:
        return 0

def uncertainty(data, mean, massAverage = False, masses = None):
    "Calculates the standard deviation of the mean"
    s = 0.
    N = 0.
    for i in range(data.size):
        if massAverage == False:
            s+= (data[i]-mean)**2
            N+=1.
        else:
            s+=masses[i]*(data[i]-mean)**2
            N+=masses[i]
    if N!=0 and N!=1:
        if massAverage == True:
            return np.sqrt(s/(N*data.size))
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