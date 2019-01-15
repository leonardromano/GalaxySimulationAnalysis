#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import numpy as np
import yt

def radius(vector):
    "returns the cylindrical Radius of a vector"
    return np.sqrt(vector[0]**2 + vector[1]**2)

def maximum(vector):
    "returns the maximum of a 1D-vector"
    return np.linalg.norm(vector, np.inf)

def cart2Pol(vectors):
    "converts [x,y] to [phi, r]"
    liste = list()
    for vector in vectors:
        phi = np.arctan2(vector[1], vector[0])*180/np.pi
        r = np.linalg.norm(vector)
        liste.append([phi, r])
    return np.asarray(liste)

def normalizeVector(vector):
    "divides vector by its norm"
    n = np.linalg.norm(vector)
    if n != 0:
        return vector/np.linalg.norm(vector)
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
    for i in vectors:
        a.append(np.linalg.norm(i))
    return np.asarray(a)

def jArray(vectors, j):
    "takes array of vectors and reduces it into an array of its jth-components"
    a = list()
    for i in vectors:
        a.append(i[j])
    return np.asarray(a)

def xyArray(vectors):
    "takes array of vectors and reduces it into an array of x- and y-components"
    a = list()
    for i in vectors:
        a.append([i[0],i[1]])
    return np.asarray(a)

def radialArray(xyVelocities, radialUnitVectors):
    "takes array of vectors and reduces it into an array of (ploar)-radial-components"
    a = list()
    i = 0
    for vel in xyVelocities:
        a.append(np.dot(vel, radialUnitVectors[i]))
        i+=1
    return np.asarray(a)

def azimuthalArray(xyVelocities, azimuthalUnitVectors):
    "takes array of vectors and reduces it into an array of (ploar)-azimuthal-components"
    a = list()
    i = 0
    for vel in xyVelocities:
        a.append(np.dot(vel, azimuthalUnitVectors[i]))
        i+=1
    return np.asarray(a)

def getRighthandedSystem(zVector):
    "returns righthanded system of unit vectors with ez in direction of zVector"
    if (zVector[0] == 0. and zVector[1] == 0. and zVector[2] == 0.) \
    or zVector.shape != (3,):
        print "Input was the zero-vector!"
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
        
def get3DArray(D1Arrays):
    "creates 3D Array from 3 1D Arrays of same length"
    array = list()
    i = 0
    for entry in D1Arrays[0]:
        array.append([entry, D1Arrays[1][i], D1Arrays[2][i]])
        i+=1
    return np.asarray(array)

def get3DList(D1Arrays):
    "creates 3D List from 3 1D Arrays of same length"
    array = list()
    i = 0
    for entry in D1Arrays[0]:
        array.append([ entry, D1Arrays[1][i], D1Arrays[2][i] ])
        i+=1
    return array

def get4DList(D1Arrays):
    "creates 4D List from 4 1D Arrays of same length"
    array = list()
    i = 0
    for entry in D1Arrays[0]:
        array.append([ entry, D1Arrays[1][i], \
                      D1Arrays[2][i], D1Arrays[3][i]])
        i+=1
    return array

def get5DList(D1Arrays):
    "creates 5D List from 5 1D Arrays of same length"
    array = list()
    i = 0
    for entry in D1Arrays[0]:
        array.append([ entry, D1Arrays[1][i], \
                      D1Arrays[2][i], D1Arrays[3][i] ,D1Arrays[4][i]])
        i+=1
    return array

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