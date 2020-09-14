#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import yt
import numpy as np

import Analysis_Code.utility.snapClass as sC
from Analysis_Code.utility.miscellaneous.vectors import getNDArray


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
    positions = getNDArray([positions_x, positions_y, positions_z])
    return [np.asarray(IDs), np.asarray(masses), np.asarray(radii), positions]

def distanceFromPoint(array, point):
    "calculates the distance from a point for all points in an array"
    distances = np.zeros(array.shape[0])
    for i in range(array.shape[0]):
        distances[i] += np.linalg.norm(point-array[i])
    return distances

def transformToCMSystem(positions, CM):
    "translates the position vectors in the CM system"
    temp = list()
    for point in positions:
        temp.append(point-CM)
    return np.asarray(temp)

def getHaloDataList(haloData, center):
    "takes all the halo Data and writes it into a big List"
    positions = transformToCMSystem(haloData[3], center)
    distances = distanceFromPoint(positions, np.array([0., 0., 0.]))
    DataArray = list()
    for i in range(haloData[0].size):
        DataArray.append([haloData[0][i], haloData[1][i], haloData[2][i], \
                      positions[i], distances[i]])
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
    print("The most massive halo has a mass of " + str(largest[1]) \
          + r" [10^10 Msun]")
    return largest

def zoomIn(data, halo, ignoreH = False, npArray = False, \
           PTCAM = ['gas', 'star'], PTCCM = ['star'], \
           onlyDisk = False, diskHeight = np.inf, lengthUnit = 'Mpc', \
           reduceToColdGas = False, Tmax = np.inf, dWAM = False):
    "Zooms in to halo"
    radius = halo[2]
    if npArray == False:
        CM = yt.YTArray(halo[3], 'kpc')
    else:
        CM = halo[3]
    print("Center of Galaxy: " + str(CM) + " [Mpc]")
    print("Virial Radius: " + str(radius) + " [Mpc]")
    sC.reduceSnapToGalaxy(data, CM, radius, npArray = npArray)
    print("There are " + str(data.starPositions.shape[0]) + " star particles left!")
    print("There are " + str(data.dmPositions.shape[0]) + " dm particles left!")
    print("There are " + str(data.gasPositions.shape[0]) + " gas particles left!")
    if reduceToColdGas == True:
        sC.reduceToColdGas(data, Tmax, npArray = npArray)
    if lengthUnit == 'Mpc':
        sC.MpcTokpc(data)
    sC.subtractCMWeighted(data, 'Velocities', npArray = npArray)
    sC.alignToHighestDensityGas(data, npArray = npArray, PTCCM = PTCCM)
    sC.alignToNewCS(data,  sC.calculateAngularMomentum \
                    (data, PTCAM, densityWeighted = dWAM), \
                        npArray = npArray)
    if onlyDisk == True:
        sC.diskCut(data, diskHeight, npArray = npArray)
        print("successfully reduced gas to disk! There are " \
        + str(data.gasPositions.shape[0]) + ' gas particles left!')
        print("successfully reduced stars to disk! There are " \
        + str(data.starPositions.shape[0]) + ' star particles left!')
    return data

def testRoutine(snap, npArray = True, PTCAM = ['gas', 'star'], \
                PTCCM = ['star'], onlyDisk = False, diskHeight = np.inf, \
                lengthUnit = 'Mpc', dWAM = False):
    "Test getHaloFromGalactic Size"
    data = sC.Snap(snap, npArray = npArray)
    print("finished loading the data!")
    z = np.random.rand(3,)
    print(z)
    sC.alignToNewCS(data, z, npArray = npArray)
    sC.kpcToMpc(data)
    sC.transformToCenterSystem(data, np.asarray([-25.,-25., -25.]))
    print("finished aligning to new CS!")
    halo = [0, 130.169979174, 0.5, np.array([25., 25., 25.]), 0.]
    return zoomIn(data, halo, npArray = npArray, \
                  PTCAM = PTCAM, PTCCM = PTCCM, onlyDisk = onlyDisk, \
                  diskHeight = diskHeight, lengthUnit = lengthUnit, dWAM = dWAM)

def getHaloFromGalacticSize(snap, haloCatalogue, hPar, npArray = False, \
                            PTCAM = ['gas', 'star'], PTCCM = ['star'], \
                            onlyDisk = False, diskHeight = np.inf, \
                            reduceToColdGas = False, Tmax = np.inf, \
                            dWAM = False):
    "gets Data of the most massive Halo in the halo Catalogue"
    data = sC.Snap(snap, npArray = npArray)
    print("finished loading the data!")
    if npArray == False:
        center = yt.YTArray(0.5*(hPar[2]+hPar[3]), 'kpc')
    else:
        center = 0.5*(hPar[2]+hPar[3])
    sC.transformToCenterSystem(data, center)
    halo = getLargestGalaxy \
    (setMassConstraints \
     (reduceToInnerRadius \
      (getHaloDataList \
       (initialiseHalos\
        (haloCatalogue, data.z, data.h), center), hPar[1]), \
       hPar[0]))
    return zoomIn(data, halo, npArray = npArray, \
                  PTCAM = PTCAM, PTCCM = PTCCM, onlyDisk = onlyDisk, \
                  diskHeight = diskHeight, lengthUnit = hPar[4], \
                  reduceToColdGas = reduceToColdGas, Tmax = Tmax, dWAM = dWAM)
