#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import distributions as dist
import rotationCurves as rC
import utility.densityProfiles as dP
import utility.Txts as Txt
from utility.vectors import *

def plotDataPointsFromURL(url, label):
    "plots the data from the URL"
    data = np.loadtxt(url)
    x = data[:,0]
    y = data[:,1]
    error = np.transpose(data[:, 2:4])
    plt.errorbar(x, y, yerr = error, fmt = 'rx', \
                 markersize = 2., elinewidth = 0.5, \
                 ecolor = 'k', capsize = 2.5, capthick=0.5, label = label)

def plotData(data_array, label, error = False, color=None):
    "plots the data if requested with uncertainty region"
    plt.plot(data_array[:,[0]], data_array[:,[1]], color = color, label = label)
    if error == True:
        x = np.asarray(data_array[:,[0]])
        y = np.asarray(data_array[:,[1]])
        error = np.asarray(data_array[:,[2]])
        plt.fill_between(x[:,0], y[:,0]-error[:,0], y[:,0]+error[:,0])

def histogram2D(x, y, axislabels, axislimits, particleType, colormap = 'r', \
              colorbarLabel = 'color', numberOfBins = 500):
    "creates 2D histogram"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    H = ax.hist2d(x, y, bins = numberOfBins, weights = colormap , \
              norm = LogNorm(), cmap="inferno")
    fig.colorbar(H[3], ax=ax, shrink=0.8, pad=0.01, orientation="vertical", \
                 label=r'$ \rho\ [M_{\odot}\ \mathrm{kpc}^{-3}]$')
    ax.set_title("projection plot of " + particleType + "-particles", \
                 va='bottom')
    plt.xlabel(axislabels[0])
    plt.ylabel(axislabels[1])
    plt.xlim(axislimits[0])
    plt.ylim(axislimits[1])
    plt.show()


def createFit(fitInformation, domain, inc):
    "create fit function from fitInformation = ([fitParameters], fitFunction)"
    fitParameters = fitInformation[0]
    fitFunction = fitInformation[1]
    x = np.arange(domain[0], domain[1], inc)
    if fitFunction == 'genMax':
        alpha = fitParameters[1]
        v0 = fitParameters[0]
        y = ((2*alpha*x**2)/(v0**3 * math.gamma(3/(2*alpha)))) * \
        np.exp(-(x/v0)**(2*alpha)) 
        return [x,y]
    if fitFunction == 'gauss':
        mu = fitParameters[0]
        hwhm = fitParameters[1]
        y = (math.sqrt(math.log(2)/math.pi)/hwhm) * \
        np.exp(-math.log(2)*((x-mu)/hwhm)**2)
        return [x,y]
    if fitFunction == 'constant':
        c = fitParameters[0]
        x = [domain[0], domain[-1]]
        y = [c, c]
        return [x,y]
    if fitFunction == 'linear':
        m = fitParameters[0]
        t = fitParameters[1]
        x = [domain[0], domain[-1]]
        y = [m*x[0]+t, m*x[1]+t]
        return [x,y]
    if fitFunction == 'NFW':
        r1 = np.exp(x*np.log(10))
        rs = np.exp(fitParameters[0]*np.log(10))
        rho = np.exp(fitParameters[1]*np.log(10))
        nfw = 4*rho/((r1/rs)*(1+(r1/rs))**2)
        y = np.log10(nfw)
        return [x, y]
        
def rotCurve(xyPositions, xyVelocities, rMax, dr, particleType, color):
    "creates rotation curve from particles velocity Data"
    radialDomain = np.arange(0., rMax, dr)
    vphi = azimuthalArray(xyVelocities, \
                          azimuthalUnitVectors(radialUnitVectors(xyPositions)))
    velocities = np.zeros(radialDomain.shape[0])
    
    N = np.zeros(radialDomain.shape[0])
    
    i=0
    for particle in xyPositions:
        n = int(round(np.linalg.norm(particle)/dr)) 
        if n < int(rMax/dr):
            velocities[n] += vphi[i]
            N[n] += 1.
        i+=1
        
    i=0
    for n in N:
        if n!=0:
            velocities[i]*=1./n
        i+=1
    nsun = int(round(8.122/dr))
    print "v(" + str((nsun-1)*dr) +" kpc) = " + str(velocities[nsun-1]) + "km/s \n" + \
    "v(" + str(nsun*dr) +" kpc) = " + str(velocities[nsun]) + "km/s \n" + \
    "v(" + str((nsun+1)*dr) +" kpc) = " + str(velocities[nsun+1]) + "km/s \n"
    Txt.createTxt(radialDomain, velocities, "rotCurve_" + particleType)
    (xArray, velArray) = getHistogramm(radialDomain, velocities, dr)
    plt.plot(xArray, velArray, color = 'C'+ str(color), label = particleType)
    
def getHistogramm(domain, image, inc):
    "takes domain and image as input and outputs it in histgramm form"
    x = list()
    y = list()
    for point in domain:
        x.append(point)
        x.append(point + inc)
    for point in image:
        y.append(point)
        y.append(point)
    return (np.asarray(x), np.asarray(y))

def domain(distribution):
    "finds the domain on which the distribution is defined"
    xmin = distribution[[0],[0]]
    xmax = distribution[[-1],[0]]
    return [xmin,xmax]

def plotVelocities(velocityDistributions, labels, fits = list()):
    "plots velocity distributions"
    v = np.arange(0., 650., 1.)
    mw = (4./(220.*math.sqrt(math.pi))) * (v/220.)**2 *np.exp(-(v/220.)**2)
    for i in range(6):
        if i == 0:
            plotData(velocityDistributions[i], labels[i], color = 'C0')
            try:
                fit = createFit(fits[i], domain(velocityDistributions[i]), 1.)
                plt.plot(fit[0],fit[1], label = 'fit', color = 'C1')
            except IndexError:
                print "There is no fit available."
            plt.plot(v, mw, label = 'SHM', linestyle = 'dashed', color = 'black')
            plt.legend()
            plt.xlabel('v [km/s]')
            plt.ylabel(r'f(v) [$(km/s)^{-1}$]')
            plt.show()
        else:
            plotData(velocityDistributions[i], labels[i], color = 'C0')
            try:
                fit = createFit(fits[i], domain(velocityDistributions[i]), 2.)
                plt.plot(fit[0],fit[1], label = 'fit', color = 'C1')
            except IndexError:
                print "There is no fit available."
            plt.legend()
            plt.xlabel('v [km/s]')
            plt.ylabel(r'f(v) [$(km/s)^{-1}$]')
            plt.show()

def main_dist(positions, velocities, interval, fits):
    "Does the thing"
    distros = dist.createVelocityHistograms(positions, velocities, interval)
    plotVelocities(distros, ['speed', 'radial', 'azimuthal', 'x', 'y', 'z'], fits)

def main_rC(snapObject, rMax, dr):
    "Does the thing with rotation Curves"
    rC_data = rC.createCombinedRotationCurve(snapObject, rMax, dr)
    i = -1
    for group in ["Combined mass", "DM", "Gas", "Star", \
                  "Disk", "Bulge", "Boundary Particles", \
                  "Stellar Particles"]:
        try:
            if group == "Combined mass":
                plt.plot(rC_data[i][:,0], rC_data[i][:,1], label = group)
                i+=1
            else:
                plt.plot(rC_data[i][:,0], rC_data[i][:,1], color = 'C' + str(i), label = group, \
                         linestyle = 'dashed')
                i+=1
        except:
            print "There is no " + group + " data to plot!"
            i+=1
    plt.legend()
    plt.ylabel(r'$v_{rot}$ [km/s]')
    plt.xlabel('R [kpc]')
    plt.show()
    
def plotOnlyEffectiveRotationCurve(snapObject, rMax, dr):
    "Does the thing with rotation Curves"
    eff = rC.createOnlyEffectiveCurve(snapObject, rMax, dr)
    plt.plot(eff[:,0], eff[:,1], color = 'black', label = 'Combined mass')

def plotNFWRotationCurve(rMax, dr, rS, rhoS):
    "plots NFW rotation curve"
    NFW = rC.createNFWRotationCurve(rMax, dr, rS, rhoS)
    plt.plot(NFW[:,0], NFW[:,1], color = 'red', label = 'NFW')

def main_dP(snapObject, bounds, dlogr, \
            fitNFW, NFWParameters, \
            fitDMDensityProfile, DMDensityFits, DMDensityFitDomains, unit, \
            plotLocalDensity = False, localDensity = (0., 0.), \
            plotOverDensity = False):
    "Plots the the log(rho)-log(R) density Profile"
    densityArray = dP.createDensityProfileHistogramm \
    (snapObject.dmMasses, normarray(snapObject.dmPositions), \
     bounds[0], bounds[1], dlogr, unit)
    plt.plot(densityArray[:,0], densityArray[:,1], label = 'density Profile')
    if fitNFW == True:
        NFWFit = createFit([[NFWParameters[0], NFWParameters[1]],'NFW'], \
                           bounds,dlogr)
        plt.plot(NFWFit[0], dP.convertToUnits(NFWFit[1], unit), \
                 color = 'C1', label = 'NFW Profile')
    if fitDMDensityProfile == True:
        constantFit = createFit(DMDensityFits[0], DMDensityFitDomains[0], \
                        DMDensityFitDomains[0][1]-DMDensityFitDomains[0][0])
        plt.plot(constantFit[0], dP.convertToUnits(constantFit[1], unit), \
                 label = 'linear fit', color = 'C3')
        linearFit = createFit(DMDensityFits[1], DMDensityFitDomains[1], \
                        DMDensityFitDomains[1][1]-DMDensityFitDomains[1][0])
        plt.plot(linearFit[0], dP.convertToUnits(linearFit[1], unit), \
                 color = 'C3')
    if plotLocalDensity == True:
        y = np.log10(localDensity[0])
        ym = np.log10(localDensity[0]-localDensity[1])
        yp = np.log10(localDensity[0]+localDensity[1])
        plt.plot(bounds, [y, y], '--', color = 'k', label = 'Eilers et al.')
        plt.fill_between(bounds, [yp, yp], [ym, ym], \
                         alpha = 0.2, facecolor = 'k')
    if plotOverDensity == True:
        h = snapObject.h
        y = -5+np.log10(200*1.05352*h**2)
        plt.plot(bounds, [y, y], '--', color = 'r', label = 'Overdensity')
    plt.legend()
    plt.xlabel(r'$log_{10}\left(R \left[kpc\right]\right)$')
    plt.ylabel(dP.getYLabel(unit))
    plt.show