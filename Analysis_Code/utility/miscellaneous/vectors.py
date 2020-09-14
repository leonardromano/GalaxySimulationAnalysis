#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import numpy as np
import yt

def radius(vector):
    "returns the cylindrical radius of a vector"
    return np.sqrt(vector[0]**2 + vector[1]**2)

def maximum(vector):
    "returns the maximum of a 1D-vector"
    return np.linalg.norm(vector, np.inf)

def cart2Pol(vectors):
    "converts [x,y] to [phi, r]"
    coords = list()
    for vector in vectors:
        coords.append(np.array([np.arctan2(vector[1], vector[0])*180/np.pi, \
                          np.linalg.norm(vector)]))
    return np.asarray(coords)

def normalizeVector(vector):
    "divides vector by its norm"
    n = np.linalg.norm(vector)
    if n != 0:
        return vector/n
    else:
        return vector

def azimuthalUnitVectors(radialUnitVectors):
    "generates array of azimuthal unit vectors"
    a = list()
    for vector in radialUnitVectors:
        a.append([-vector[1],vector[0]])
    return np.asarray(a)

def radialUnitVectors(xyVectors):
    "generates array of polar-radial unit vectors"
    a = list()
    for vector in xyVectors:
        a.append(normalizeVector(vector))
    return np.asarray(a)

def normarray(vectors):
    "takes array of vectors and turns it into an array of norms"
    a = list()
    for vector in vectors:
        a.append(np.linalg.norm(vector))
    return np.asarray(a)

def jArray(vectors, j):
    "takes array of vectors and reduces it into an array of its jth-components"
    a = list()
    for vector in vectors:
        a.append(vector[j])
    return np.asarray(a)

def xyArray(vectors):
    "takes array of vectors and reduces it into an array of x- and y-components"
    a = list()
    for vector in vectors:
        a.append([vector[0],vector[1]])
    return np.asarray(a)

def radialArray(xyVelocities, radialUnitVectors):
    "takes array of vectors and reduces it into an array of (ploar)-radial-components"
    a = np.zeros(xyVelocities.shape[0])
    for i in range(xyVelocities.shape[0]):
        a[i] += np.dot(xyVelocities[i], radialUnitVectors[i])
    return a

def azimuthalArray(xyVelocities, azimuthalUnitVectors):
    "takes array of vectors and reduces it into an array of (ploar)-azimuthal-components"
    a = np.zeros(xyVelocities.shape[0])
    for i in range(xyVelocities.shape[0]):
        a[i] += np.dot(xyVelocities[i], azimuthalUnitVectors[i])
    return a

def getRighthandedSystem(zVector):
    "returns righthanded system of unit vectors with ez in direction of zVector"
    if (zVector[0] == 0. and zVector[1] == 0. and zVector[2] == 0.) \
    or zVector.shape != (3,):
        print("Input was the zero-vector!")
        return [np.array([1., 0., 0.]), np.array([0., 1., 0.]), np.array([0., 0., 1.])]
    else:
        z = normalizeVector(zVector)
        if z[0] == 0. and z[1] == 0.:
            return [np.array([1., 0., 0.]), np.array([0., 1., 0.]), z]
        elif z[0] == 0. and z[2] == 0.:
            return [np.array([1., 0., 0.]), z, np.array([0., 0., 1.])]
        elif z[1] == 0. and z[2] == 0.:
            return [z, np.array([0., 1., 0.]), np.array([0., 0., 1.])]
        elif z[0] == 0:
            x = normalizeVector(np.array([0., z[2], -z[1]]))
            y = np.cross(z, x)
            return [x, y, z]
        elif z[1] == 0:
            x = normalizeVector(np.array([z[2], 0, -z[0]]))
            y = np.cross(z, x)
            return [x, y, z]
        elif z[2] == 0:
            x = normalizeVector(np.array([z[1], -z[0], 0.]))
            y = np.cross(z, x)
            return [x, y, z]
        else:
            x = normalizeVector(np.array([z[1], -z[0], 0.]))
            y = np.cross(z, x)
            return [x, y, z]

def getNDList(D1Arrays):
    "creates ND List from list of N 1D Arrays of same length"
    a = list()
    for i in range(D1Arrays[0].size):
        b = list()
        for j in range(len(D1Arrays)):
            b.append(D1Arrays[j][i])
        a.append(b)
    return a

def getNDArray(D1Arrays):
    "creates ND Array from list of N 1D Arrays of same length"
    return np.asarray(getNDList(D1Arrays))

def seperateList3D(listObject, npArray = False):
    "seperates 3D-list into 3 1D-arrays"
    array1 = list()
    array2 = list()
    array3 = list()
    for element in listObject:
        array1.append(element[0])
        array2.append(element[1])
        array3.append(element[2])
    if npArray == False:
        return [yt.YTArray(array1, 'Msun'), \
                yt.YTArray(array2, 'kpc'), \
                yt.YTArray(array3, "km/s")]
    else:
        return [np.asarray(array1), np.asarray(array2), \
                np.asarray(array3)]

def seperateList4D(listObject, npArray = False):
    "seperates 4D-list into 4 1D-arrays"
    array1 = list()
    array2 = list()
    array3 = list()
    array4 = list()
    for element in listObject:
        array1.append(element[0])
        array2.append(element[1])
        array3.append(element[2])
        array4.append(element[3])
    if npArray == False:
        return [yt.YTArray(array1, 'Msun'), \
                yt.YTArray(array2, 'kpc'), \
                yt.YTArray(array3, 'km/s'), \
                yt.YTArray(array4, 'K')]
    else:
        return [np.asarray(array1), np.asarray(array2), \
                np.asarray(array3), np.asarray(array4)]
    
def seperateList5D(listObject, npArray = False):
    "seperates 5D-list into 5 1D-arrays"
    array1 = list()
    array2 = list()
    array3 = list()
    array4 = list()
    array5 = list()
    for element in listObject:
        array1.append(element[0])
        array2.append(element[1])
        array3.append(element[2])
        array4.append(element[3])
        array5.append(element[4])
    if npArray == False:
        return [yt.YTArray(array1, 'Msun'), \
                yt.YTArray(array2, 'kpc'), \
                yt.YTArray(array3, 'km/s'), \
                yt.YTArray(array4, 'g/cm^3'), \
                yt.YTArray(array5, 'K')]
    else:
        return [np.asarray(array1), np.asarray(array2), \
                np.asarray(array3), np.asarray(array4), np.asarray(array5)]
        
def seperateList6D(listObject, npArray = False):
    "seperates 5D-list into 5 1D-arrays"
    array1 = list()
    array2 = list()
    array3 = list()
    array4 = list()
    array5 = list()
    array6 = list()
    for element in listObject:
        array1.append(element[0])
        array2.append(element[1])
        array3.append(element[2])
        array4.append(element[3])
        array5.append(element[4])
        array6.append(element[5])
    if npArray == False:
        return [yt.YTArray(array1, 'Msun'), \
                yt.YTArray(array2, 'kpc'), \
                yt.YTArray(array3, 'km/s'), \
                yt.YTArray(array4, 'g/cm^3'), \
                yt.YTArray(array5, 'kpc'), \
                yt.YTArray(array6, 'K')]
                
    else:
        return [np.asarray(array1), np.asarray(array2), \
                np.asarray(array3), np.asarray(array4), \
                np.asarray(array5), np.asarray(array6)]
