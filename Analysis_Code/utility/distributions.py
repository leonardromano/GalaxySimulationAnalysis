#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

from Analysis_Code.utility.miscellaneous.vectors import *
import numpy as np

def bins(numberOfBins, symmetric = False):
    "creates an array of bins"
    a = list()
    if symmetric == False:
        for i in range(numberOfBins-1):
            a.append([i,0.,0.])
        return np.asarray(a)
    else:
        for i in range(2*numberOfBins-1):
            a.append([i,0.,0.])
        return np.asarray(a)

def fillBin(bins, data_array, interval, symmetric = False):
    "sorts data by bins"
    if symmetric == False:
        for dataPoint in data_array:
            bins[int(round(dataPoint/interval))] += [0,1.,1.]
        for e in bins:
            if e[2] == 0:
                e[2] = 1.
            e[2] = np.sqrt(e[2])
        return bins
    else:
        for dataPoint in data_array:
            bins[int(round(dataPoint/interval))+bins.shape[0]//2] += [0,1. ,1.]
        for e in bins:
            if e[2] == 0:
                e[2] = 1.
            e[2] = np.sqrt(e[2])
        return bins

def getNormalizedDistribution(filledBinArray, interval, split, \
                              NumberOfParticles, symmetric = False):
    "normalizes the filledBinArray"
    if NumberOfParticles == 0:
        return filledBinArray
    else:
        if symmetric == False:
            for entry in filledBinArray:
                entry[0] *= interval
                entry[1] *= 1./float(interval*NumberOfParticles)
                entry[2] *= 1./float(interval*NumberOfParticles)
            return filledBinArray
        else:
            for entry in filledBinArray:
                entry[0] = interval*(entry[0]-split)
                entry[1] *= 1./float(NumberOfParticles*interval)
                entry[2] *= 1./float(interval*NumberOfParticles)
            return filledBinArray

def distribute(interval, speed_array, NumberOfParticles, symmetric = False):
    "creates a distribution of the speed_array"
    split = minimalSplit(speed_array, interval)
    a = fillBin(bins(split+2, symmetric), speed_array, interval, symmetric)
    return getNormalizedDistribution(a, interval, split , NumberOfParticles, \
                                     symmetric)

def minimalSplit(array, interval):
    "returns minimal Split to plot given data array"
    try:
        return int(round(np.linalg.norm(array , np.inf)/interval))
    except ValueError:
        return 1

def readSnap(positions, velocities, particleType):
    "reads snap data and distributes it into arrays of interest"
    xyVel = xyArray(velocities)
    xyPos = xyArray(positions)
    rads = radialUnitVectors(xyPos)
    azimuthals = azimuthalUnitVectors(rads)
    N = velocities.shape[0]
    return (N, normarray(velocities), radialArray(xyVel, rads), \
            azimuthalArray(xyVel, azimuthals), \
            jArray(velocities, 0), jArray(velocities, 1), jArray(velocities, 2))

def createVelocityDistributions(positions, velocities, \
                                interval, particleType = 'dm'):
    "creates a speed, a radial, a z-component and an azimuthal \
    velocity distribution from snap"
    rS = readSnap(positions, velocities, particleType)
    distributions = [distribute(interval, rS[1], rS[0])]
    for i in range(5):
        distributions.append(distribute(interval,  rS[i+2], rS[0], True))
    return distributions
    
def distToHist(distribution, interval):
    "creates a histogramlike distribution"
    a = list()
    for [x,y,e] in distribution:
        a.append([x,y,e])
        a.append([x+interval,y,e])
    return np.asarray(a)

def createVelocityHistograms(positions, velocities, interval, \
                             particleType = 'dm', saveText = False):
    "creates a speed, a radial, a z-component and an azimuthal velocity \
    histogram from snap"
    distributions = createVelocityDistributions(positions, velocities, \
                                                interval, particleType)
    data = list()
    for distribution in distributions:
        data.append(distToHist(distribution, interval))
    if saveText == True:
        for i in range(len(data)):
            np.savetxt('PlotData/velDistro'+str(i), data[i], delimiter = ' ')
    return data
    
    
    
