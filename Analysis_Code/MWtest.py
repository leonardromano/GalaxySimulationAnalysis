#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

from pygadgetreader import *
import snapClass as snappy

def totalMass(masses):
    "calculates the total Mass of an array of masses"
    m = 0
    try:
        for mass in masses:
            m += mass
        return m
    except:
        print "An Error occured in MWtest.py -> totalMass"
        return 0

def MWlike_arrayInput(snapObject):
    "Tests if galaxy is MW like"
    Mhalo = totalMass(snapObject.dmMasses)
    Mgas = totalMass(snapObject.gasMasses)
    Mstellar = totalMass(snapObject.starMasses) + \
    totalMass(snapObject.diskMasses) + \
    totalMass(snapObject.bulgeMasses)
    
    Mtot = Mhalo + Mgas + Mstellar
    print "Mgas = " + str(Mgas)
    if 50<Mhalo<10**4 and 4.5<Mstellar<8.3:
        print "Mhalo = " + str(Mhalo) + " and Mstellar = " + \
        str(Mstellar) + "\n so the galaxy is MW-like!"
        print "total Mass = " + str(Mtot)
        return (True, Mhalo, Mstellar, Mgas)
    else:
        print "Mhalo = " + str(Mhalo) + " and Mstellar = " + \
        str(Mstellar) + "\n so the galaxy is NOT MW-like!"
        print "total Mass = " + str(Mtot)
        return (False, Mhalo, Mstellar, Mgas)