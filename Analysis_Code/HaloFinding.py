#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import yt
import numpy as np
from pygadgetreader import *
import snapClass as snappy
from utility.vectors import *


def initialiseHalos(haloCatalogue, z, h):
    "creates a list of all data aspects of the catalogue"
    ds = yt.load(haloCatalogue)
    ad = ds.all_data()
    ds.define_unit("M_halo", (1.0e10, "Msun"))
    masses = ad["halos", "particle_mass"].to("M_halo")
    IDs = ad["halos", "particle_identifier"]
    radii = ad["halos", "virial_radius"].in_units('kpc')*(1.+z)*h
    positions_x = ad["halos", "particle_position_x"].in_units('kpc')*(1.+z)*h
    positions_y = ad["halos", "particle_position_y"].in_units('kpc')*(1.+z)*h
    positions_z = ad["halos", "particle_position_z"].in_units('kpc')*(1.+z)*h
    positions = get3DArray([positions_x, positions_y, positions_z])
    return [np.asarray(IDs), np.asarray(masses), np.asarray(radii), positions]

def distanceFromPoint(array, point):
    "calculates the distance from a point for all points in an array"
    temp = list()
    for element in array:
        temp.append(np.linalg.norm(point-element))
    return np.asarray(temp)

def transformToCMSystem(positions, CM):
    "translates the position vectors in the CM system"
    temp = list()
    for point in positions:
        temp.append(point-CM)
    return np.asarray(temp)

def getHaloDataList(haloData, lowerLeft, upperRight):
    "takes all the halo Data and writes it into a big List"
    center = 0.5*(lowerLeft + upperRight)
    positions = transformToCMSystem(haloData[3], center)
    distances = distanceFromPoint(positions, np.array([0., 0., 0.]))
    i=0
    DataArray = list()
    for ID in haloData[0]:
        DataArray.append([ID, haloData[1][i], haloData[2][i], \
                      positions[i], distances[i]])
        i+=1
    return DataArray

def reduceToInnerRadius(haloDataList, rMax):
    "reduces the data to the data in the inner sphere of radius rMax"
    reduced = list()
    for halo in haloDataList:
        if halo[4] <= rMax:
            reduced.append(halo)
    return reduced

def setMassConstraints(array, bounds = (0, 10**3)):
    "reduces the array to a certain mass region"
    liste = list()
    for halo in array:
        if bounds[0] < halo[1] <bounds[1]:
            liste.append(halo)
    return liste

def getLargestGalaxy(HaloDataList):
    "searches the data list for the largest Galaxy"
    largest = [0., 0., 0., 0., 0.]
    for halo in HaloDataList:
        if halo[1] > largest[1]:
            largest = halo
    print "The most massive halo has a mass of " + str(largest[1]) + r" [10^10 Msun]"
    return largest

def zoomIn(data, halo, ignoreH = False, npArray = False, \
           PTCAM = ['gas', 'star'], PTCCM = ['star'], \
           onlyDisk = False, diskHeight = np.inf, lengthUnit = 'Mpc', \
           reduceToColdGas = False, Tmax = np.inf):
    "Zooms in to halo"
    radius = halo[2]
    if npArray == False:
        CM = yt.YTArray(halo[3], 'kpc')
    else:
        CM = halo[3]
    print "Center of Galaxy: " + str(CM) + " [Mpc]"
    print "Virial Radius: " + str(radius) + " [Mpc]"
    snappy.reduceSnapToGalaxy(data, CM, radius, npArray = npArray)
    print "There are " + str(data.starPositions.size) + " star particles left!"
    print "There are " + str(data.dmPositions.size) + " dm particles left!"
    print "There are " + str(data.gasPositions.size) + " gas particles left!"
    if reduceToColdGas == True:
        snappy.reduceToColdGas(data, Tmax, npArray = npArray)
    if lengthUnit == 'Mpc':
        snappy.MpcTokpc(data)
    snappy.subtractCMWeighted(data, 'Velocities', npArray = npArray)
    snappy.alignToHighestDensityGas(data, npArray = npArray, PTCCM = PTCCM)
    snappy.alignToNewCS(data, snappy.calculateAngularMomentum(data, PTCAM), \
                        npArray = npArray)
    if onlyDisk == True:
        snappy.diskCut(data, diskHeight, npArray = npArray)
        print "successfully reduced gas to disk! There are " \
        + str(data.gasPositions.size) + ' gas particles left!'
        print "successfully reduced stars to disk! There are " \
        + str(data.starPositions.size) + ' star particles left!'
    return data

def testRoutine(snap, npArray = True, PTCAM = ['gas', 'star'], \
                PTCCM = ['star'], onlyDisk = False, diskHeight = np.inf, \
                lengthUnit = 'Mpc'):
    "Test getHaloFromGalactic Size"
    data = snappy.Snap(snap, npArray = npArray)
    print "finished loading the data!"
    z = np.random.rand(3,)
    print z
    snappy.alignToNewCS(data, z, npArray = npArray)
    snappy.kpcToMpc(data)
    snappy.transformToCenterSystem(data, np.asarray([-25.,-25., -25.]))
    print "finished aligning to new CS!"
    halo = [0, 130.169979174, 0.5, np.array([25., 25., 25.]), 0.]
    return zoomIn(data, halo, npArray = npArray, \
                  PTCAM = PTCAM, PTCCM = PTCCM, onlyDisk = onlyDisk, \
                  diskHeight = diskHeight, lengthUnit = lengthUnit)

def getHaloFromGalacticSize(snap, haloCatalogue, hPar, npArray = False, \
                            PTCAM = ['gas', 'star'], PTCCM = ['star'], \
                            onlyDisk = False, diskHeight = np.inf, \
                            reduceToColdGas = False, Tmax = np.inf):
    "gets Data of the most massive Halo in the halo Catalogue"
    data = snappy.Snap(snap, npArray = npArray)
    print "finished loading the data!"
    if npArray == False:
        center = yt.YTArray(0.5*(hPar[2]+hPar[3]), 'kpc')
    else:
        center = 0.5*(hPar[2]+hPar[3])
    snappy.transformToCenterSystem(data, center)
    halo = getLargestGalaxy \
    (setMassConstraints \
     (reduceToInnerRadius \
      (getHaloDataList \
       (initialiseHalos\
        (haloCatalogue, data.z, data.h), hPar[2], hPar[3]), hPar[1]), \
       hPar[0]))
    return zoomIn(data, halo, npArray = npArray, \
                  PTCAM = PTCAM, PTCCM = PTCCM, onlyDisk = onlyDisk, \
                  diskHeight = diskHeight, lengthUnit = hPar[4], \
                  reduceToColdGas = reduceToColdGas, Tmax = Tmax)