#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 22:57:06 2020

@author: leonard
"""
import numpy as np
import math

def PowerLaw(x, A, n):
    "returns a power law with powerlaw index n"
    return A*x**n

def redChiQuadrat(y, yerr, ymodel, Nparams):
    "returns the chisquared of a fit"
    result = 0
    i=0
    for val in y:
        result+=((val-ymodel[i])/yerr[i])**2
        i+=1
    return result/(y.size-1-Nparams)

def NFW(x, rhos, rs):
    "returns a NFW density profile"
    return rhos/(x/rs)/(1+x/rs)**2

def Burkert(x, rhos, rs):
    "returns a Burkert density profile"
    return rhos/(1+x/rs)/(1+(x/rs)**2)

def Maxwellian(x, vc):
    "Returns a Maxwellian distribution"
    return (4./(vc*np.sqrt(np.pi))) * (x/vc)**2 *np.exp(-(x/vc)**2)

def GeneralizedMaxwellian(x, vc, alpha):
    "Returns a generalized Maxwellian distribution"
    return ((2*alpha*x**2)/(vc**3 * math.gamma(3/(2*alpha))))*np.exp(-(x/vc)**(2*alpha))

def Gaussian(x, mu, s2):
    "Returns a Gaussian distribution"
    return 1/np.sqrt(2*np.pi*s2)*np.exp(-(x-mu)**2 /2/s2)