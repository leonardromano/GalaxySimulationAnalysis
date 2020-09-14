#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit as chisq

import Analysis_Code.utility.distributions as dist
import Analysis_Code.utility.rotationCurves as rC
import Analysis_Code.utility.miscellaneous.densityProfiles as dP
from Analysis_Code.utility.miscellaneous.vectors import *
import Analysis_Code.utility.miscellaneous.functions as func

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
    plt.plot(data_array[:,0], data_array[:,1], color = color, label = label)
    if error == True:
        x = np.asarray(data_array[:,0])
        y = np.asarray(data_array[:,1])
        error = np.asarray(data_array[:,2])
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
        
def rotCurve(xyPositions, xyVelocities, rMax, dr, particleType, color, \
             massAverage = False, mass = None, saveText = False):
    "creates rotation curve from particles velocity Data by averaging over radial bins"
    radialDomain = np.arange(0., rMax, dr)
    vphi = azimuthalArray(xyVelocities, \
                          azimuthalUnitVectors(radialUnitVectors(xyPositions)))
    velocities = np.zeros(radialDomain.shape[0])
    M = np.zeros(radialDomain.shape[0])
    
    for i in range(xyPositions.shape[0]):
        n = int(round(np.linalg.norm(xyPositions[i])/dr))
        if n < int(rMax/dr):
            if massAverage == False:
                velocities[n] += vphi[i]
                M[n] +=1.
            else:
                m = mass[i]
                velocities[n] += m*vphi[i]
                M[n] += m
        
    for i in range(M.size):
        if M[i]!=0:
            velocities[i]*=1./M[i]
    
    if saveText == True:
        np.savetxt("PlotData/rotCurve_" + particleType, \
                   getNDArray([radialDomain, velocities]), delimiter = ' ')
        
    histogram = getHistogramm(radialDomain, velocities, dr)
    plt.plot(*histogram, color = 'C'+ str(color), label = particleType)
    
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
    return np.asarray(x), np.asarray(y)

def domain(distribution):
    "finds the domain on which the distribution is defined"
    xmin = distribution[[0],[0]]
    xmax = distribution[[-1],[0]]
    return [xmin,xmax]

def plotVelocities(velocityDistributions, labels, fits = list(), \
                   saveImages = False, Isolated = False):
    "plots velocity distributions"
    for i in range(6):
        if i == 0:
            plotData(velocityDistributions[i], labels[i], color = 'C0')
            x = velocityDistributions[i][:,0]
            y = velocityDistributions[i][:,1]
            yerr = velocityDistributions[i][:,2]
            try:
                popt, pcov = chisq(eval("func." +fits[i][0]), x, y, sigma = yerr, \
                                   bounds = fits[i][1])
                plt.plot(x, eval("func." +fits[i][0]+"(x, *popt)"), \
                         label = fits[i][0] + '\n' + \
                         r'$v_{c}$=%5.4f+-%5.4f'%(popt[0], np.sqrt(abs(pcov[0,0])))+ \
                         '\n' + \
                         r'$\alpha$=%5.4f+-%5.4f'%(popt[1], np.sqrt(abs(pcov[1,1])))+ \
                         '\n' + r'$\chi^2_{\nu}$ = %5.3f'\
                         %(func.redChiQuadrat(y, yerr, \
                         eval("func." + fits[i][0] +"(x, *popt)"), popt.size)), \
                             color = 'C1')
            except IndexError:
                print("There is no fit available.")
            popt, pcov = chisq(func.Maxwellian, x, y, sigma = yerr, \
                               bounds = (50, 600))
            plt.plot(x, func.Maxwellian(x, *popt), \
                         label = 'Maxwellian \n' + \
                         r'$v_{c}$=%5.4f+-%5.4f'%(popt[0], np.sqrt(abs(pcov[0,0])))+ \
                         '\n' + \
                         r'$\chi^2_{\nu}$ = %5.3f'\
                         %(func.redChiQuadrat(y, yerr, \
                           func.Maxwellian(x, *popt), popt.size)), \
                             linestyle = 'dashed', color = 'black')
            plt.legend(bbox_to_anchor=(1, 1))
            plt.xlabel('v [km/s]')
            plt.ylabel(r'f(v) [$(km/s)^{-1}$]')
            if saveImages == True:
                filename = "Images/DM_speedDistribution.png"
                plt.savefig(filename)
            plt.show()
        else:
            plotData(velocityDistributions[i], labels[i], color = 'C0')
            x = velocityDistributions[i][:,0]
            y = velocityDistributions[i][:,1]
            yerr = velocityDistributions[i][:,2]
            try:
                popt, pcov = chisq(eval("func." +fits[i][0]), x, y, sigma = yerr, \
                                   bounds = fits[i][1])
                plt.plot(x, eval("func." +fits[i][0]+"(x, *popt)"), \
                         label = fits[i][0] + '\n' + \
                         r'$\mu$=%5.4f+-%5.4f'%(popt[0], np.sqrt(abs(pcov[0,0])))+ \
                         '\n' + \
                         r'$\sigma^2$=%5.4f+-%5.4f'%(popt[1], np.sqrt(abs(pcov[1,1])))+ \
                         '\n' + r'$\chi^2_{\nu}$ = %5.3f'\
                         %(func.redChiQuadrat(y, yerr, \
                         eval("func." + fits[i][0] +"(x, *popt)"), popt.size)), \
                             color = 'C1')
            except IndexError:
                print("There is no fit available.")
            plt.legend()
            plt.xlabel('v [km/s]')
            plt.ylabel(r'f(v) [$(km/s)^{-1}$]')
            if saveImages == True:
                filename = "Images/DM_velocityDistribution_" + str(i) + ".png"
                plt.savefig(filename)
            plt.show()

def main_dist(positions, velocities, interval, fits, \
              saveImages = False, saveText = False, Isolated = False):
    "creates histograms and plots them"
    histograms = dist.createVelocityHistograms(positions, velocities, \
                                            interval, saveText = saveText)
    plotVelocities(histograms, ['speed', 'radial', 'azimuthal', 'x', 'y', 'z'], \
                   fits, saveImages, Isolated = Isolated)
    
def main_rC(snapObject, rMax, dr, saveImages = False):
    "Makes rotation Curve plots"
    histograms = rC.createCombinedRotationCurve(snapObject, rMax, dr)
    i = -1
    for group in ["Combined mass", "DM", "Gas", "Star", \
                  "Disk", "Bulge", "Boundary Particles", \
                  "Stellar Particles"]:
        try:
            if group == "Combined mass":
                plt.plot(*histograms[i], label = group)
                i+=1
            else:
                plt.plot(*histograms[i], color = 'C' + str(i), label = group, \
                         linestyle = 'dashed')
                i+=1
        except:
            print("There is no " + group + " data to plot!")
            i+=1
    plt.legend()
    plt.ylabel(r'$v_{rot}$ [km/s]')
    plt.xlabel('R [kpc]')
    if saveImages == True:
        filename = "Images/RotationCurve_fromMassDistribution.png"
        plt.savefig(filename)
    plt.show()
    
def plotOnlyEffectiveRotationCurve(snapObject, rMax, dr, saveText = False):
    "Makes effective rotation Curve plots"
    histogram = rC.createOnlyEffectiveCurve(snapObject, rMax, dr, saveText = saveText)
    plt.plot(*histogram, color = 'black', label = 'Combined mass')

def main_dP(snapObject, bounds, dlogr, unit,\
            fitDP, DPShapes, \
            plotLocalDensity = False, localDensity = (0., 0.), \
            plotOverDensity = False, saveImages = False, saveText = False):
    "Plots the the log(rho)-log(R) density Profile"
    densityArrays = dP.createDensityProfileHistogramm \
    (snapObject.dmMasses, normarray(snapObject.dmPositions), \
     *bounds, dlogr, unit, fitDP = fitDP, saveText = saveText)
    x = densityArrays[0]
    y = densityArrays[1]
    rbounds = [10**bounds[0], 10**bounds[1]]
    plt.plot(x, y, label = 'density Profile')
    if fitDP == True:
        yerr = densityArrays[2]
        for [shape, bounds] in DPShapes:
            popt, pcov = chisq(eval("func." + shape), x, y, sigma = yerr, \
                               bounds = bounds)
            plt.plot(x, eval("func."+shape+"(x, *popt)"), '-',label=\
                     shape + "-profile" + '\n' + \
                         r'$\rho_{s}$=%5.4f+-%5.4f'%(popt[0], np.sqrt(abs(pcov[0,0])))+ \
                         '\n' + \
                         r'$r_{s}$=%5.4f+-%5.4f'%(popt[1], np.sqrt(abs(pcov[1,1])))+ \
                         '\n' + r'$\chi^2_{\nu}$ = %5.3f'\
                         %(func.redChiQuadrat(y, yerr, \
                         eval("func." + shape +"(x, *popt)"), popt.size)))
    if plotLocalDensity == True:
        y = localDensity[0]
        ym = localDensity[0]-localDensity[1]
        yp = localDensity[0]+localDensity[1]
        plt.plot(rbounds, [y, y], '--', color = 'k', label = 'Eilers et al.')
        plt.fill_between(rbounds, [yp, yp], [ym, ym], \
                         alpha = 0.2, facecolor = 'k')
    if plotOverDensity == True:
        h = snapObject.h
        y = 200*1.05352*h**2*10**(-5)
        plt.plot(rbounds, [y, y], '--', color = 'r', label = 'Overdensity')
    plt.legend(bbox_to_anchor=(1, 1))
    plt.xlabel(r'$r \left[kpc\right]$')
    plt.ylabel(dP.getYLabel(unit))
    plt.xscale("log")
    plt.yscale("log")
    if saveImages == True:
        filename = "Images/DM_DensityProfile.png"
        plt.savefig(filename)
    plt.show()
    plt.close()
