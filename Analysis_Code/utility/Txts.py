#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analysis code by Leonard Romano
"""

import numpy as np

def createTxt(x_data, y_data, filename):
    "create txt file containing the data"
    np.savetxt(filename, ArrayFrom1DArrays(x_data, y_data), delimiter = ' ')
    
def ArrayFrom1DArrays(x_data,y_data):
    "creates 2D Array from 2 1D Arrays"
    liste = list()
    i=0
    for x in x_data:
        liste.append([x, y_data[i]])
        i+=1
    return np.asarray(liste)