#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 18:27:06 2017

@author: Shane Nichols

As an exercise to learn NumPy, I coded the Matlab function mmPartialWave and
all of its dependencies in Python. This code is about 5 times faster than the
Matlab version because the delta matrix was constructed in a more efficient way
otherwise it would be about half the speed of the Matlab versions.
All dependencies, including the material libraries, are in mm.py.
I only put a few materials in the library to test the code.

run this script in an IPython console by first navigating to the folder mmCalc,
and then inputting:
    %run -t partialWave_testing.py

I'll point out one key difference between the two versions. Whereas in Matlab,
I arrange arrays of matrices as (i, j, N, M, O...), in Python it is better to
do (N, M, O..., i, j). i.e., matrix dimenions at the end.
"""
import numpy as np
import mm


layer1 = ['air', 0, np.array([0, 0, 0]), False, True]
layer2 = ['TiO2', 20, np.array([0, 90, 0]), False, False]
layer3 = ['+quartz', 100000, np.array([24, 0.3, 0]), True, False]
layer4 = ['BK7', 20, np.array([45, 90, 0]), False, True]
layer5 = ['TiO2', 0, np.array([32, 1.2, 0]), False, False]
layerArray = [layer1, layer2, layer3, layer4, layer5]


MM = mm.mmPartialWave(layerArray, np.arange(300, 700, 0.1), 20, True, False)
print(MM[0])
