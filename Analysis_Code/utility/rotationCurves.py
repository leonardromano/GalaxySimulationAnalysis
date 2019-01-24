#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import miscellaneous.vectors as vec
import numpy as np
import miscellaneous.Txts as Txt
import miscellaneous.densityProfiles as dP

G = 6.674*10**(-11)
kpc = 30.9*10**18
M = 2*10**40

def getCombinedArray(massesList, positionsList):
    "reads all sorts of stellar mass and writes them in a single array"
    masses = concat(massesList)
    positions = concat(positionsList)
    return (masses, positions)

def createRotationCurve(masses, positions, dr, rMax):
    "Creates rotation curve from snap"
    radial = vec.normarray(vec.xyArray(positions))
    radii = radiusRange(0.1 , 8., rMax, dr)
    vrot = list()
    for r in radii:
        vrot.append([r, np.sqrt(M*G*getInsideMass(masses, radial, r)/ \
                               (r*kpc))*10**(-3)])
    return np.asarray(vrot)

def radiusRange(rmin, r0, rmax, dr):
    "divides radius Range in inner and outer region"
    if rmin > rmax:
        print "rmax < rmin"
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

def distToHist(distribution, dr):
    "Turns distribution into histogram"
    a = list()
    for [r, v] in distribution:
        if r<8.0:
            a.append([r, v])
            a.append([r+dr/4, v])
        else:
            a.append([r, v])
            a.append([r+dr, v])
    return np.asarray(a)

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
            print "there are no " + particleType + "-particles!"
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
        print "there are no " + particleType + "-particles!"
    H = list()
    for dist in C:
        H.append(distToHist(dist, dr))
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
    effective = createRotationCurve(Masses, Positions, dr, rMax)
    if saveText == True:
        Txt.createTxt(effective[:,0], effective[:,1], \
                      "PlotData/effectiveRotationCurve")
    return distToHist(effective, dr)

def createNFWRotationCurve(rMax, dr, rS, rhoS):
    "creates NFWRotationCUrve"
    radii = radiusRange(0.1 , 8., rMax, dr)
    masses = dP.getNFWMasses(radii*(1./rS), (4.*np.pi/30.)*rhoS*rS**3)
    velocities = list()
    i=0
    for m in masses:
        velocities.append([np.sqrt(M*G*m/(radii[i]*kpc))*10**(-3)])
        i+=1
    return Txt.ArrayFrom1DArrays(radii, np.asarray(velocities))
    
