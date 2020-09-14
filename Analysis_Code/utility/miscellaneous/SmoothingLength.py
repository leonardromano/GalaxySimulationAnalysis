#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 22:50:53 2020

@author: leonard
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as chisq
import Analysis_Code.utility.miscellaneous.functions as func


def PlotSmoothingLength(data, fromGalacticSnap=False, saveImages=False):
    
    if fromGalacticSnap == True:
        x = data.gasDensity*10**(-8)
    else:
        x = 10*data.gasDensity
    y = data.gasSmoothingLength
    popt, pcov = chisq(func.PowerLaw, x, y)
    plt.plot(x, func.PowerLaw(x, *popt), 'r-', label=\
                 r'$h_{sml}\propto \rho^{\gamma}$'+ \
                 '\n' + \
                 r'$\gamma$=%5.4f+-%5.4f'%(popt[1], np.sqrt(abs(pcov[1,1]))))
    plt.scatter(x,y, label = "gas particles")
    plt.xlabel(r"density$[M_\odot\,pc^{-3}]$")
    plt.ylabel(r"$h_{sml}[kpc]$")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    if saveImages == True:
        filename = "Images/SmoothingLengthVSDensity"
        plt.savefig(filename)
    plt.show()
    plt.close()