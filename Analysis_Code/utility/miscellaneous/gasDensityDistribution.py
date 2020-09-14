#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Kazunori Okamoto, translated into python by Leonard Romano
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

mp = 1.672621637*10**(-24)

def W(r, h):
    "calculates the kernel function for the distance r and the smoothing length h"
    a = r/h
    if 0. <=a<=0.5:
        return (1.-6.*a**2+6.*a**3)/(np.pi*h**3)
    elif 0.5 <= a <= 1.:
        return (2.*(1.-a)**3)/(np.pi*h**3)
    elif 1.<a:
        return 0.

def densityDistribution(data, binNumber = 100, rmax = 15.):
    """calculates the column density including the smoothing length
    for gas particles"""
    rho_3d = np.zeros((binNumber,binNumber,binNumber))
    w_3d = np.zeros((binNumber,binNumber,binNumber))
    rho_xz = np.zeros((binNumber,binNumber))
    rho_xy = np.zeros((binNumber,binNumber))
    pos = data.gasPositions
    hsml = data.gasSmoothingLength
    mass = data.gasMasses
    density = data.gasDensity
    
    print("calculating spatial density...\n")
    
    dr = 2.0*rmax/binNumber
    l = 0
    for pt in pos:
        boundsi = [int((rmax + pt[0] - hsml[l])/dr -0.5 - 2.), \
                   int((rmax + pt[0] + hsml[l])/dr -0.5 + 2.)]
        boundsj = [int((rmax + pt[1] - hsml[l])/dr -0.5 - 2.), \
                   int((rmax + pt[1] + hsml[l])/dr -0.5 + 2.)]
        boundsk = [int((rmax + pt[2] - hsml[l])/dr -0.5 - 2.), \
                   int((rmax + pt[2] + hsml[l])/dr -0.5 + 2.)]
        if boundsi[0]<0:
            boundsi[0] = 0
        if boundsj[0]<0:
            boundsj[0] = 0
        if boundsk[0]<0:
            boundsk[0] = 0
        if boundsi[0]>=binNumber:
            l+=1
            continue
        if boundsj[0]>=binNumber:
            l+=1
            continue
        if boundsk[0]>=binNumber:
            l+=1
            continue
        if boundsi[1]>binNumber:
            boundsi[1] = binNumber
        if boundsj[1]>binNumber:
            boundsj[1] = binNumber
        if boundsk[1]>binNumber:
            boundsk[1] = binNumber
        if boundsi[1]<=0:
            l+=1
            continue
        if boundsj[1]<=0:
            l+=1
            continue
        if boundsk[1]<=0:
            l+=1
            continue
        for i in range(boundsi[0], boundsi[1]):
            for j in range(boundsj[0], boundsj[1]):
                for k in range(boundsk[0], boundsk[1]):
                    ds = np.linalg.norm(np.array([-rmax + dr*(i+0.5) - pt[0], \
                                                  -rmax + dr*(j+0.5) - pt[1], \
                                                  -rmax + dr*(k+0.5) - pt[2]]))
                    rho_3d[i,j,k] += mass[l]*W(ds, hsml[l])/mp
                    w_3d[i,j,k] += (mass[l]/density[l])*W(ds, hsml[l])
        l+=1
    for i in range(binNumber):
        for j in range(binNumber):
            for k in range(binNumber):
                if w_3d[i,j,k] != 0:
                    rho_3d[i,j,k] /= w_3d[i,j,k]
    for i in range(binNumber):
        for j in range(binNumber):
            for k in range(binNumber):
                rho_xz[k,i] += rho_3d[i,j,k] * dr
    for i in range(binNumber):
        for j in range(binNumber):
            for k in range(binNumber):
                rho_xy[j,i] += rho_3d[i,j,k] * dr
    print("finished calculating the spatial density!\n")
    return rho_xz, rho_xy

def makePlot(X, Y, data, cmap, mode, saveImages):
    "Makes a colormesh plot and optionally saves it"
    plt.pcolormesh(X, Y, data, cmap = cmap, \
                   norm = LogNorm(), shading = 'gouraud')
    plt.colorbar(label='column gas mass density')
    plt.title("Density distribution of gas (" +mode+")")
    plt.xlabel("X [kpc]")
    if mode == "edge on":
        plt.ylabel("Z [kpc]")
    else:
        plt.ylabel("Y [kpc]")
    if saveImages == True:
        filename = "Images/columndensity_"+ mode +".png"
        plt.savefig(filename)
    plt.show()
    plt.close()

def smoothColumnDensityPlot(data, binNumber = 100, rmax = 15., \
                            cmap = "inferno", saveImages = False):
    "Makes a Column Density plot including smoothing length"
    density_xz, density_xy = densityDistribution(data, binNumber, rmax)
    XX = np.linspace(-rmax, rmax, binNumber)
    YY = np.linspace(-rmax, rmax, binNumber)
    X,Y = np.meshgrid(XX,YY)
    makePlot(X,Y, density_xy, cmap, "face on", saveImages)
    makePlot(X,Y, density_xz, cmap, "edge on", saveImages)
    