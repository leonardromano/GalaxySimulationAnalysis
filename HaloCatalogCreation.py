#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog

def createHaloCatalogueFromSnap(snap):
    "creates Halo Catalogue from a given snap"
    data_ds = yt.load(snap)
    hc = HaloCatalog(data_ds=data_ds, finder_method='hop')
    hc.create()

"""
Just insert the name of the galactic scale simulation snap inside 
here and run the code to create the corresponding halo catalog.
"""

snap = 'Snapshots/Zoom-in-runs/Osaka_HALO_01_173_10_z0.0.hdf5'
createHaloCatalogueFromSnap(snap)