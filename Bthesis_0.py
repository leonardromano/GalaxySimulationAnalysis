#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

from Analysis_Code.UserInterface import main
from numpy import array

""" Just get all the information you need from the snap using main."""

"""
Put here the path of the hdf5 file with the data:
Also select the Particle Types used for the Calculation of Angular Momentum (PTCAM)
(to determine the z-axis)
OR densityWeightedAngularMomentum = True which then calculates the density 
weighted mean angular momentum of the gas (regardless of PTCAM)
And the Particle Types used for the Calculation of the Center of Mass (PTCCM)
(to determine the origin (highest density gas particle within 1kpc of the 
                          center of mass))
If saveImages = True, the images are saved in the folder "Images"
If saveText = True, the data for most of the plots is saved in the folder "PlotData"
"""
snap = 'path-of-file.hdf5'
PTCAM = ['star']
densityWeightedAngularMomentum = True
PTCCM = ['star']
saveImages = True
saveText = True
"""
If you are loading the data from galactic zoom-in simulation data set
fromGalacticSimulationData to True
The parameters are the following:
(halo mass limits[10¹⁰Msun], maximum distance from center[lengthUnits], 
coordinates of lower left corner, coordinates of upper right corner, 
lengthUnits)
There are some options:
    - useHaloCatalog = True will result in the halo being extracted from 
        the haloCatalog = 'path-of-file.0.h5'
        (you can create one with HaloCatalogCreation.py)
    - If useHaloCatalog is False, the coordinates of the halo and the
        virial radius will have to be inserted manually in the
        haloCoordinates = (coordinates of the halo, virial radius)
    - useNumpyArrays = True makes sure that all arrays used in the calculations
        are numpy arrays (this is recommended as yt arrays use more memory and
        are slower)
    - alignOnlyToGasAngularMomentum = True will set the z-axis to only the gas
        angular momentum as opposed to gas and stars
    - onlyDiskParticles = True reduces the gas and star particles to those with
        |z|< diskHeight[kpc]
    - onlyColdGas = True reduces the gas particles to those with T < Tmax[K]
"""
fromGalacticSimulationData = True
useHaloCatalog = True
haloCoordinates = (array([0., 0., 0.]), 0.200)
haloCatalog = 'path-of-file.0.h5'
haloParameters = ((1., 10**3), 5., array([0., 0., 0.]), \
                  array([50., 50., 50.]), 'Mpc')
useNumpyArrays = True
onlyDiskParticles = True
diskHeight = 1.
onlyColdGas = True
Tmax = 10**(4.0)

"""
If you want to create DM velocity distributions set 
makeDMVelocityDistributions to True
The parameters are the following: (binsize [km/s], fits)
fits is a list of the fit parameters for each plot (|v|, vr, vphi, vx, vy, vz)
the fit parameters are of the form [[parameters], 'type of fit']
So far there are 5 types of fits:
 1. 'genMax' (generalised Maxwellian): [v0, alpha]
 2. 'gauss': [m, hwhm]
 3. 'constant': [c]
 4. 'linear' : [m, y0]
 5. 'NFW' : [log10(rs), log10(rhos)]
"""
makeDMVelocityDistributions = False
fits = list()
velocityDistributionParameters = (10., fits)



"""
If you want to create a rotation curve set makeRotationCurve to True
The parameters are the following: (rmax[kpc], binsize[kpc], particleTypes)
particlesTypes has the following form ['gas', 'star', 'disk', 'bulge', 'dm']
(of course the exact types can vary)
There are also some options:
    - fromMassDistribution = True makes a plot of only the rotation curve 
        as predicted by the mass distribution
    - massAverage = True calculates the rotational velocity as mass weighted 
        average for each bin (as opposed to just numerical average)
    - addObservationalData = True adds the observational data points stored 
        in a txt file which is pointed at with 
    - URL = 'path-of-file.txt'
    - plotNFW = True plots the rotation curve predicted by only a NFW profile
        in the same picture
    - rotNFWparameters = (rS, rhoS)
        are the parameters for the NFW plot
    - compareToMassDistribution = True plots the rotation curve as predicted 
        by the mass distribution in the same picture
"""
makeRotationCurve = False
rotationCurveParameters = (25., 0.2, ['gas', 'star'])
fromMassDistribution = False
massAverage = True
addObservationalData = True
URL = 'path-of-file.txt'
plotNFW = False
rotNFWParameters = (20., 0.6)
compareToMassDistribution = False



"""
If you want to compute the total mass components of DM, gas and star particles
and want to test if the halo is similar to the MW set 
testMilkyWaylikeness to True
"""
testMilkyWaylikeness = False


"""
If you want to plot the DM density profile set makeDMDensityProfile to True
The parameters are (bounds(logr[kpc]), binsize(log[kpc]), density unit)
density unit can be 'Msun'(Msun/kpc³) or 'GeV' (GeV/cm³)
There are also some options:
    - includeNFWFit = True plots a NFW profile with the parameters
    - NFWParameters = (log10(rs), log10(rhos))
    - includeLocalDensity = True plots a horizontal bar at
    - localDensity = (localDensity, error)
    - includeOverdensity = True plots a horizontal line at the height 
        of the overdensity
    - includeDensityFit = True plots a fitted profile with the parameters
    - DMDensityFitParameters = [[parameters], type of fit] 
        (see velocity Distributions for more information)
    - DMDensityFitDomains = [bounds(fit0), bound(fit1),...]
"""
makeDMDensityProfile = False
DMDensityProfileParameters = ([0, 2.0], 0.01, 'GeV')
includeNFWFit = True
NFWParameters = (1.3, -2.5)
includeLocalDensity = True
localDensity = (0.30, 0.03)
includeOverdensity = True
includeDensityFit = False
DMDensityFitParameters = [[[0.], 'constant'], \
                          [[0., 0.], 'linear']]
DMDensityFitDomains = [[0, 0.6], [0.6, 2.0]]

"""
If you want to plot the SPH smoothing length against the gas density for \
each particle set smoothingLengthVSDensity = True
"""
smoothingLengthVSDensity = False

"""
If you want to make log column density projection plots of the gas particles
set makeProjectionPlots to True
The parameters are (number of bins (100 would correspond to 100x100x100 bins), 
                    |xmax|, colorscheme) 
"""
makeProjectionPlots = False
projectionPlotParameters = (100, 15., "inferno")



"""
If you want to calculate the local DM density or the local rotational velocity
set getLocalDMDensity or getLocalRotationalVelocity to True
The parameters are the following: 
(solar radius[kpc], minimal binsize, number of steps, binsize increment)
The algorithm calculates the average value for each binsize and then averages 
over all obtained values
"""
getLocalDMDensity = True
getLocalRotationalVelocity = True
localValueParameters = (8.122, 0.05, 20, 0.05)



main(snap, \
     fromGalacticSnap = fromGalacticSimulationData, 
     fromHaloCatalogue = useHaloCatalog, haloCatalogue = haloCatalog, \
     haloCoords = haloCoordinates, hPar = haloParameters, \
     npArray = useNumpyArrays, PTCAM = PTCAM, \
     dWAM = densityWeightedAngularMomentum, PTCCM = PTCCM,\
     onlyDisk = onlyDiskParticles, diskHeight = diskHeight, \
     reduceToColdGas = onlyColdGas, Tmax = Tmax, \
     doDist = makeDMVelocityDistributions, \
     distParameters = velocityDistributionParameters, \
     doRotCurve = makeRotationCurve, rcPar = rotationCurveParameters, \
     fromMassDist = fromMassDistribution, massAverage = massAverage, \
     addObsData = addObservationalData, url = URL, \
     compareToNFW = plotNFW, rNFWPar = rotNFWParameters, \
     compareToMassDist = compareToMassDistribution, \
     testMW = testMilkyWaylikeness, \
     doDMDensityProfile = makeDMDensityProfile, \
     dpParameters = DMDensityProfileParameters, \
     fitNFW = includeNFWFit, NFWParameters = NFWParameters , \
     plotLocalDensity = includeLocalDensity, localDensity = localDensity, \
     plotOverDensity = includeOverdensity, \
     fitDMDensityProfile = includeDensityFit, \
     DMDensityFits = DMDensityFitParameters, \
     DMDensityFitDomains = DMDensityFitDomains, \
     doSmoothingLengthPlot = smoothingLengthVSDensity, \
     doProjectionPlots = makeProjectionPlots, \
     ppPar = projectionPlotParameters, \
     getLocalDens = getLocalDMDensity, \
     getLocalVrot = getLocalRotationalVelocity,\
     lvPar = localValueParameters , saveImages = saveImages, \
     saveText = saveText \
     )
