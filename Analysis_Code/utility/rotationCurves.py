#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import numpy as np

from Analysis_Code.utility.miscellaneous.vectors import xyArray, normarray
from Analysis_Code.utility.miscellaneous.vectors import getNDArray

G = 6.674*10**(-11)
kpc = 30.9*10**18
M = 2*10**40

def getCombinedArray(massesList, positionsList):
    "reads all sorts of stellar mass and writes them in a single array"
    masses = concat(massesList)
    positions = concat(positionsList)
    return masses, positions

def createRotationCurve(masses, positions, dr, rMax):
    "Creates rotation curve from snap"
    radial = normarray(xyArray(positions))
    radii = radiusRange(0.1 , 8., rMax, dr)
    vrot = np.zeros(radii.size)
    for i in range(radii.size):
        vrot[i] += np.sqrt(M*G*getInsideMass(masses, radial, radii[i])/ \
                               (radii[i]*kpc))*10**(-3)
    return radii, vrot

def radiusRange(rmin, r0, rmax, dr):
    "divides radius Range in inner and outer region"
    if rmin > rmax:
        print("rmax < rmin")
        return radiusRange(rmax,rmin,dr)
    elif rmax < r0:
        return np.arange(rmin, rmax, dr/4)
    elif rmin > r0:
        return np.arange(rmin,rmax, dr)
    else:
        return concat((np.arange(rmin, r0, dr/4), np.arange(r0, rmax, dr)))
    
def concat(arrays):
    "concatinates two arrays"
    temp = list()
    for array in arrays:
        for x in array:
            temp.append(x)
    return np.asarray(temp)

def getInsideMass(masses, radial, r):
    "returns total Mass insinside cylindre with radius r"
    m=0
    i=0
    for radius in radial:
        if radius<=r:
            m += masses[i]
        i+=1
    return m

def distToHist(radii, velocities, dr):
    "Turns distribution into histogram"
    rs = list()
    vs = list()
    for i in range(radii.size):
        if radii[i]<8.0:
            rs.append(radii[i])
            vs.append(velocities[i])
            rs.append(radii[i]+dr/4)
            vs.append(velocities[i])
        else:
            rs.append(radii[i])
            vs.append(velocities[i])
            rs.append(radii[i]+dr)
            vs.append(velocities[i])
    return np.asarray(rs), np.asarray(vs)

def createCombinedRotationCurve(snapObject, rMax, dr):
    "Creates combined rotation curve"
    C = list()
    for particleType in ["dm", 'gas', 'star', 'disk', 'bulge', 'bndry']:
        try:
            C.append(createRotationCurve \
                     (eval('snapObject.' + particleType + 'Masses'), \
                      eval('snapObject.' + particleType + 'Positions'), \
                      dr, rMax))
        except AttributeError:
            print("there are no " + particleType + "-particles!")
    try:
        (stellarMasses, stellarPositions) = \
        getCombinedArray \
        ((snapObject.starMasses, snapObject.diskMasses, snapObject.bulgeMasses), \
         (snapObject.starPositions, snapObject.diskPositions, \
          snapObject.bulgePositions))
        C.append(createRotationCurve(stellarMasses, stellarPositions, dr, rMax))
        (combinedMasses, combinedPositions) = \
        getCombinedArray((stellarMasses, snapObject.dmMasses, \
                          snapObject.gasMasses), \
        (stellarPositions, snapObject.dmPositions, snapObject.gasPositions))
        C.append(createRotationCurve(combinedMasses, \
                                     combinedPositions, dr, rMax))
    except AttributeError:
        print("there are no " + particleType + "-particles!")
    H = list()
    for dist in C:
        H.append(distToHist(*dist, dr))
    return H
    
def createOnlyEffectiveCurve(snapObject, rMax, dr, saveText = False):
    "Creates only effective rotation curve"
    (Masses, Positions) = \
    getCombinedArray \
    ((snapObject.starMasses, snapObject.diskMasses, snapObject.bulgeMasses, \
      snapObject.dmMasses, snapObject.gasMasses), \
        (snapObject.starPositions, snapObject.diskPositions, \
         snapObject.bulgePositions, snapObject.dmPositions, \
         snapObject.gasPositions))
    dist = createRotationCurve(Masses, Positions, dr, rMax)
    if saveText == True:
        np.savetxt("PlotData/effectiveRotationCurve", \
                   getNDArray(dist), delimiter = ' ')
    return distToHist(*dist, dr)
    
