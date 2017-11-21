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
    %run -plt partialWave_testing.py

I'll point out one key difference between the two versions. Whereas in Matlab,
I arrange arrays of matrices as (i, j, N, M, O...), in Python it is better to
do (N, M, O..., i, j). i.e., matrix dimenions at the end.
"""
import numpy as np
import matplotlib.pyplot as plt
import mm


layer1 = ['air', 0, np.array([0, 0, 0]), False, True]
layer2 = ['TiO2', 20, np.array([0, 90, 0]), False, False]
layer3 = ['+quartz', 100000, np.array([24, 0.3, 0]), True, False]
layer4 = ['BK7', 20, np.array([45, 90, 0]), False, True]
layer5 = ['TiO2', 0, np.array([32, 1.2, 0]), False, False]
layerArray = [layer1, layer2, layer3, layer4, layer5]

wavelengths = np.arange(300, 700, 0.1)
MM = mm.mmPartialWave(layerArray, wavelengths, 20, True, False)

# plot the result: not sure if this will display in all editors. I am using 
# Spyder with the IPython graphics backend set to Automatic.
figW = 11
figH = 6.5
plt.rcParams.update({'font.size': 8,
                     'ytick.direction': 'in',
                     'xtick.direction': 'in'})
fig, ax_lst = plt.subplots(4, 4, sharex='col', figsize=(figW, figH))
for x in range(4):
    for y in range(4):
        ax_lst[x, y].plot(wavelengths, MM[:, x, y])
fig.tight_layout()
fig.subplots_adjust(wspace=0.35, hspace=0.05)
