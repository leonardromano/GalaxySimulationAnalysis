#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

from pygadgetreader import *
from utility.vectors import *
import numpy as np
import yt

class Snap:
    "A class to contain all useful information contained in the snap"
    def __init__(self, snap, includeBndry = False, npArray = False):
        "initializes the Snap object"
        
        ds = yt.load(snap)
        ad = ds.all_data()
        ds.define_unit("M_halo", (1.0e10, "Msun"))
        
        self.dmMasses = yt.YTArray([])
        self.dmPositions = yt.YTArray([])
        self.dmVelocities = yt.YTArray([])
        self.dm = yt.YTArray([])
        
        self.starMasses = yt.YTArray([])
        self.starPositions = yt.YTArray([])
        self.starVelocities = yt.YTArray([])
        self.star = yt.YTArray([])
        
        self.diskMasses = yt.YTArray([])
        self.diskPositions = yt.YTArray([])
        self.diskVelocities = yt.YTArray([])
        self.disk = yt.YTArray([])
        
        self.bulgeMasses = yt.YTArray([])
        self.bulgePositions = yt.YTArray([])
        self.bulgeVelocities = yt.YTArray([])
        self.bulge = yt.YTArray([])
        
        self.gasMasses = yt.YTArray([])
        self.gasPositions = yt.YTArray([])
        self.gasVelocities = yt.YTArray([])
        self.gasTemperatures = yt.YTArray([])
        self.gasDensity = yt.YTArray([])
        self.gas = yt.YTArray([])
        
        self.bndryMasses = yt.YTArray([])
        self.bndryPositions = yt.YTArray([])
        self.bndryVelocities = yt.YTArray([])
        self.bndry = yt.YTArray([])
        
        self.h = readheader(snap, 'h')
        print "h = " + str(self.h)
        self.z = readheader(snap, 'redshift')
        print "z = " + str(self.z)
        
        if npArray == False:
            typeList = ['PartType0', 'PartType1', 'PartType2', \
                        'PartType3', 'PartType4']
            if includeBndry == True:
                typeList.append("PartType5")
            
            for particleType in typeList:
                try:
                    for variable in ['Masses', 'Coordinates', 'Velocities']:
                        self.setVariable \
                        (ad[particleType, variable], \
                         particleType, variable, \
                         h = self.h, z = self.z, load = True)
                    pT = vanillaParticleTypeFormat(particleType)
                    if particleType == 'PartType0':
                        try:
                            self.setVariable \
                            (np.asarray(ad['PartType0', 'Temperature']), \
                             particleType, 'Temperature')
                        except ZeroDivisionError:
                            print "There is no Temperature information"
                        try:
                            self.setVariable \
                            (np.asarray(ad['PartType0', 'Density']), \
                             particleType, 'Density')
                        except ZeroDivisionError:
                            print "There is no Density information"
                    self.updateCombinedData(pT)
                    print "successfully loading " + \
                    str(eval('self.' + pT + 'Positions.size')) + \
                    ' particles from '+ particleType
                except:
                    print "There are no " + particleType + "-particles!"
        else:
            typeList = ['gas', 'dm', 'disk', \
                        'bulge', 'star']
            if includeBndry == True:
                typeList.append("bndry")
            
            for particleType in typeList:
                try:
                    for variable in ['mass', 'pos', 'vel']:
                        self.setVariable \
                        (readsnap(snap, variable, particleType), \
                         particleType, variable)
                    if particleType == 'gas':
                        try:
                            self.setVariable \
                            (np.asarray(ad['PartType0', 'Temperature']), \
                             particleType, 'Temperature')
                        except ZeroDivisionError:
                            print "There is no Temperature information"
                        try:
                            self.setVariable \
                            (np.asarray(readsnap(snap, 'rho', particleType)), \
                             particleType, 'rho')
                        except ZeroDivisionError:
                            print "There is no density information"
                    self.updateCombinedData(particleType)
                    print "successfully loading " + \
                    str(eval('self.' + particleType + 'Positions.size')) + \
                    ' particles from '+ particleType
                except SystemExit:
                    print "There are no " + particleType + "-particles!"
            
            
    def setVariable(self, value, particleType, variableName, \
                    h = 1., z = 0., load = False):
        "sets a variables Value"
        
        if particleType == 'PartType1' or particleType == 'dm':
            if variableName == 'Masses' or variableName == 'mass':
                if load == False:
                    self.dmMasses = value
                else:
                    self.dmMasses = value.to("M_halo")
            if variableName == 'Coordinates' or variableName == 'pos':
                if load == False:
                    self.dmPositions = value
                else:
                    self.dmPositions = value.in_units('kpc')*h*(1+z)
            if variableName == 'Velocities' or variableName == 'vel':
                if load == False:
                    self.dmVelocities = value
                else:
                    self.dmVelocities = value.in_units('km/s')
            if variableName == particleType:
                self.dm = value
                
        if particleType == 'PartType4' or particleType == 'star':
            if variableName == 'Masses' or variableName == 'mass':
                if load == False:
                    self.starMasses = value
                else:
                    self.starMasses = value.to("M_halo")
            if variableName == 'Coordinates' or variableName == 'pos':
                if load == False:
                    self.starPositions = value
                else:
                    self.starPositions = value.in_units('kpc')*h*(1+z)
            if variableName == 'Velocities' or variableName == 'vel':
                if load == False:
                    self.starVelocities = value
                else:
                    self.starVelocities = value.in_units('km/s')
            if variableName == particleType:
                self.star = value
                
        if particleType == 'PartType2' or particleType == 'disk':
            if variableName == 'Masses' or variableName == 'mass':
                if load == False:
                    self.diskMasses = value
                else:
                    self.diskMasses = value.to("M_halo")
            if variableName == 'Coordinates' or variableName == 'pos':
                if load == False:
                    self.diskPositions = value
                else:
                    self.diskPositions = value.in_units('kpc')*h*(1+z)
            if variableName == 'Velocities' or variableName == 'vel':
                if load == False:
                    self.diskVelocities = value
                else:
                    self.diskVelocities = value.in_units('km/s')
            if variableName == particleType:
                self.disk = value
        
        if particleType == 'PartType3' or particleType == 'bulge':
            if variableName == 'Masses' or variableName == 'mass':
                if load == False:
                    self.bulgeMasses = value
                else:
                    self.bulgeMasses = value.to("M_halo")
            if variableName == 'Coordinates' or variableName == 'pos':
                if load == False:
                    self.bulgePositions = value
                else:
                    self.bulgePositions = value.in_units('kpc')*h*(1+z)
            if variableName == 'Velocities' or variableName == 'vel':
                if load == False:
                    self.bulgeVelocities = value
                else:
                    self.bulgeVelocities = value.in_units('km/s')
            if variableName == particleType:
                self.bulge = value
                
        if particleType == 'PartType0' or particleType == 'gas':
            if variableName == 'Masses' or variableName == 'mass':
                if load == False:
                    self.gasMasses = value
                else:
                    self.gasMasses = value.to("M_halo")
            if variableName == 'Coordinates' or variableName == 'pos':
                if load == False:
                    self.gasPositions = value
                else:
                    self.gasPositions = value.in_units('kpc')*h*(1+z)
            if variableName == 'Velocities' or variableName == 'vel':
                if load == False:
                    self.gasVelocities = value
                else:
                    self.gasVelocities = value.in_units('km/s')
            if variableName == "Temperature":
                self.gasTemperatures = value
            if variableName == "Density" or variableName == "rho":
                self.gasDensity = value
            if variableName == particleType:
                self.gas = value
                
        if particleType == 'PartType5' or particleType == 'bndry':
            if variableName == 'Masses' or variableName == 'mass':
                if load == False:
                    self.bndryMasses = value
                else:
                    self.bndryMasses = value.to("M_halo")
            if variableName == 'Coordinates' or variableName == 'pos':
                if load == False:
                    self.bndryPositions = value
                else:
                    self.bndryPositions = value.in_units('kpc')*h*(1+z)
            if variableName == 'Velocities' or variableName == 'vel':
                if load == False:
                    self.bndryVelocities = value
                else:
                    self.bndryVelocities = value.in_units('km/s')
            if variableName == particleType:
                self.bndry = value
        
            
    def updateCombinedData(self, particleType):
        "updates the combined data, list"
        if particleType == 'gas' and self.gasTemperatures.size != 0:
            self.setVariable \
            (get5DList([self.gasMasses, \
                        self.gasPositions, \
                        self.gasVelocities, \
                        self.gasDensity, \
                        self.gasTemperatures,] \
        ), particleType, particleType)
        elif particleType == 'gas':
            self.setVariable \
            (get4DList([self.gasMasses, \
                        self.gasPositions, \
                        self.gasVelocities, \
                        self.gasDensity] \
        ), particleType, particleType)
        else:
            self.setVariable \
            (get3DList([eval('self.' + particleType + 'Masses'), \
                        eval('self.' + particleType + 'Positions'), \
                        eval('self.' + particleType + 'Velocities')]), \
                        particleType, particleType)

def vanillaParticleTypeFormat(particleType):
    "converts the numbered particle Type Format to Vanilla"
    if particleType == 'PartType0':
        return 'gas'
    if particleType == 'PartType1':
        return 'dm'
    if particleType == 'PartType2':
        return 'disk'
    if particleType == 'PartType3':
        return 'bulge'
    if particleType == 'PartType4':
        return 'star'
    if particleType == 'PartType5':
        return 'bndry'

def reduceSnapToGalaxy(data, CM, radius, npArray = False):
    "reduces the data of a whole Zoom-in snap to only one galaxy"
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            temp = list()
            for particle in eval('data.' + particleType):
                if npArray == False:
                    pos = yt.YTArray(particle[1], 'kpc')
                else:
                    pos = particle[1]
                r = np.linalg.norm( CM - pos )
                if r <= radius:
                    pos -= CM
                    particle[1] = pos
                    temp.append(particle)
            variableList = ['Masses', 'Coordinates', 'Velocities']
            if particleType == 'gas' and data.gasTemperatures.size != 0:
                tempArrays = seperateList5D(temp, npArray)
                variableList.append('Density')
                variableList.append('Temperature')
            elif particleType == 'gas':
                tempArrays = seperateList4D(temp, npArray)
                variableList.append('Density')
            else:
                tempArrays = seperateList3D(temp, npArray)
            i=0
            for variableType in variableList:
                data.setVariable(tempArrays[i], particleType, variableType)
                i += 1
            data.updateCombinedData(particleType)
        except AttributeError:
            print "there are no " + particleType + "-particles!"
    print "finished reducing the snap to Galaxy!"
    
def alignToNewCS(snap, zDirection, npArray = False):
    "aligns the data along a CS with z-axis parallel to angular momentum"
    [x, y, z] = getRighthandedSystem(zDirection)
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            [Pos, Vel] = [list(), list()]
            for particle in eval('snap.' + particleType):
                [pos, vel] = alignToAxis(particle, x, y, z)
                [Pos.append(pos), Vel.append(vel)]
            if npArray == False:
                snap.setVariable(yt.YTArray(Pos, 'kpc'), particleType, 'Coordinates')
                snap.setVariable(yt.YTArray(Vel, 'km/s'), particleType, 'Velocities')
            else:
                snap.setVariable(np.asarray(Pos), particleType, 'Coordinates')
                snap.setVariable(np.asarray(Vel), particleType, 'Velocities')
            snap.updateCombinedData(particleType)
            print "finished rotating all the " + particleType + '-particles!'
        except AttributeError:
            print "there are no " + particleType + "-particles!"
            
def MpcTokpc(snapObject):
    "Converts the position data in snapObject from Mpc to kpc"
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            snapObject.setVariable \
            (eval('snapObject.' + particleType + 'Positions')*10**3, \
             particleType, 'Coordinates')
            snapObject.updateCombinedData(particleType)
        except AttributeError:
            print "there are no " + particleType + "-particles!"
    
def kpcToMpc(snapObject):
    "Converts the position data in snapObject from Mpc to kpc"
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            snapObject.setVariable \
            (eval('snapObject.' + particleType + 'Positions')*10**(-3), \
             particleType, 'Coordinates')
            snapObject.updateCombinedData(particleType)
        except AttributeError:
            print "there are no " + particleType + "-particles!"

def calculateAngularMomentum(snap, pTufCoAM = ['gas', 'star']):
    "calculates the mean angular Momentum of disk Particles"
    ls = list()
    typeList = pTufCoAM
    for pT in typeList:
        try:
            i=0
            for pos in eval('snap.' + pT + 'Positions'):
                ls.append \
                (eval('snap.' + pT + 'Masses[' + str(i) +' ]') * \
                 np.cross(np.asarray(pos), np.asarray(eval('snap.' + pT + \
                          'Velocities[' + str(i) +' ]'))))
                i += 1
        except AttributeError:
            print "There are no " + pT + "-particles!"
    lArray = np.asarray(ls)
    if lArray.size == 0:
        print "there are no particles in this region!"
        return np.asarray([0., 0., 0.])
    else:
        L = mean(lArray)
        print "Angular Momentum: " + str(L)
        return L

def mean(array):
    "calculates the mean of an array of vectors"
    M = np.zeros(array[0].shape)
    i = 0.
    for vector in array:
        M += vector
        i += 1.
    if i!=0.:
        M *= 1./i
    return M
        

def alignToAxis(particle, x, y, z):
    "aligns the coordinates of a particle with x, y and z"
    pos = np.asarray(particle[1])
    vel = np.asarray(particle[2])
    r = [np.dot(x, pos), np.dot(y, pos), np.dot(z, pos)]
    v = [np.dot(x, vel), np.dot(y, vel), np.dot(z, vel)]
    return [yt.YTArray(r, 'kpc'), yt.YTArray(v, 'km/s')]

def Attribute(attribute):
    "Returns the proper attribute"
    if attribute == 'Coordinates':
        return 'Positions'
    else:
        return attribute

def PositionOfDensestGasParticle(snapObject, particleTypes = ['star']):
    "Returns the position of the densest gas particle"
    R = 0
    M = 0
    for pT in particleTypes:
        try:
            i=0
            for r in eval('snapObject.' + pT + 'Positions'):
                m = eval('snapObject.' + pT + 'Masses[' + str(i) +' ]')
                R+=m*r
                M+=m
                i += 1
        except AttributeError:
            print "There are no " + pT + "-particles!"
    gas_positions = snapObject.gasPositions
    density = snapObject.gasDensity
    if M!= 0:
        R*=1./M
    else:
        R = np.array([0., 0., 0.])
    i=0
    maxdens = 0
    K = 0
    KIsZero = True
    for rho in density:
        r = gas_positions[i]
        if np.linalg.norm(r-R)<=1.:
            if rho>maxdens:
                maxdens = rho
                K = r
                KIsZero = False
        i+=1
    if KIsZero == True:
        return R
    else:
        return K
            

def alignToHighestDensityGas(snapObject, npArray = False, PTCCM = ['star']):
    "Shifts everything to the coordinate System centered around the densest \
    central gas particle"
    R = PositionOfDensestGasParticle(snapObject, particleTypes = PTCCM)
    print "position of densest gas particle: " + str(R)
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            vectors = list()
            for v in eval('snapObject.' + particleType + "Positions"):
                vectors.append(v-R)
            if npArray == False:
                snapObject.setVariable \
                (yt.YTArray(vectors, 'kpc'), \
                 particleType, 'Coordinates')
            else:
                snapObject.setVariable(np.asarray(vectors), \
                                       particleType, 'Coordinates')
            snapObject.updateCombinedData(particleType)
        except AttributeError:
            print "There are no " + particleType + "-particles!"

def subtractCMWeighted(snapObject, attribute, npArray = False, onlyStar = False):
    "subtracts the velocity/position of the center of mass from \
    the particles velocity or position"
    V = getCMWeighted(snapObject, Attribute(attribute), \
                      npArray = npArray, onlyStar = onlyStar)
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            vectors = list()
            for v in eval('snapObject.' + particleType + Attribute(attribute)):
                vectors.append(v-V)
            if npArray == False:
                snapObject.setVariable \
                (yt.YTArray(vectors, getUnit(attribute)), \
                 particleType, attribute)
            else:
                snapObject.setVariable(np.asarray(vectors), \
                                       particleType, attribute)
            snapObject.updateCombinedData(particleType)
        except AttributeError:
            print "There are no " + particleType + "-particles!"

def getCMWeighted(snapObject, attribute, npArray = False, onlyStar = False):
    "calculates the center of mass velocity of the snap"
    V = np.array([0., 0., 0.,])
    M = 0
    if onlyStar == False:
        for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
            try:
                i=0
                for v in eval('snapObject.' + particleType + attribute):
                    m = eval('snapObject.' + particleType + \
                             'Masses' + '[' + str(i) + ']')
                    V += m*np.asarray(v)
                    M += m
                    i += 1
            except AttributeError:
                print "There are no " + particleType + "-particles!"
    else:
        try:
            i=0
            for v in eval('snapObject.star' + attribute):
                m = snapObject.starMasses[i]
                V += m*np.asarray(v)
                M += m
                i += 1
        except AttributeError:
            print "There are no " + particleType + "-particles!"
    if npArray == False:
        unit = getUnit(attribute)
    
        if M!=0:
            result = yt.YTArray(V*(1./M), unit)
            print result
            return result
        else:
            return yt.YTArray([0., 0., 0.,], unit)
    else:
        if M!=0:
            result = V*(1./M)
            print result
            return result
        else:
            return np.array([0., 0., 0.,])

def getUnit(attribute):
    "returns the default unit of a given attribute"
    if attribute == 'Velocities':
        return 'km/s'
    if attribute == 'Masses':
        return 'Msun'
    if attribute == 'Coordinates' or 'Positions':
        return 'kpc'
    if attribute == 'Temperature':
        return 'K'

def transformToCenterSystem(snapObject, center):
    "translates the position vectors in the center system"
    for particleType in ['dm', 'star', 'disk', 'bulge', 'gas', 'bndry']:
        try:
            temp = list()
            for point in eval('snapObject.' + particleType + "Positions"):
                temp.append(point - center)
            snapObject.setVariable(np.asarray(temp), particleType, 'Coordinates')
            snapObject.updateCombinedData(particleType)
        except AttributeError:
            print "There are no " + particleType + "-particles!"
    print "finished transforming to the center System!"
    
def diskCut(snapObject, diskHeight, npArray = False):
    "reduces the data to disk"
    for particleType in ['gas', 'star']:
        try:
            temp = list()
            for particle in eval('snapObject.' + particleType):
                if abs(particle[1][2]) < diskHeight:
                    temp.append(particle)
            variableList = ['Masses', 'Coordinates', 'Velocities']
            if particleType == 'gas' and snapObject.gasTemperatures.size != 0:
                tempArrays = seperateList5D(temp, npArray)
                variableList.append('Density')
                variableList.append('Temperature')
            elif particleType == 'gas':
                tempArrays = seperateList4D(temp, npArray)
                variableList.append('Density')
            else:
                tempArrays = seperateList3D(temp, npArray)
            i=0
            for variableType in variableList:
                snapObject.setVariable(tempArrays[i], particleType, variableType)
                i += 1
            snapObject.updateCombinedData(particleType)
        except AttributeError:
            print "there are no " + particleType + "-particles!"

def reduceToColdGas(snapObject, Tmax, npArray = False):
    "reduces the data to cold Gas only"
    try:
        if snapObject.gasTemperatures.size!=0:
            print 'reducing the gas to particles below ' + str(Tmax) + ' K...'
            temp = list()
            for particle in snapObject.gas:
                if particle[4] < Tmax:
                    temp.append(particle)
            tempArrays = seperateList5D(temp, npArray)
            i=0
            for variableType in ['Masses', 'Coordinates', \
                                 'Velocities', 'Density', \
                                 'Temperature']:
                snapObject.setVariable(tempArrays[i], 'gas', variableType)
                i += 1
            snapObject.updateCombinedData('gas')
            print 'finished reducing to Cold Gas!'
            print 'There are ' + str(snapObject.gasPositions.size) + \
            ' gas-particles left'
        else:
            print 'There is no Temperature information!'
    except AttributeError:
        print "there are no gas-particles!"