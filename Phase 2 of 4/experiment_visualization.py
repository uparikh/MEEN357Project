#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy.interpolate import *
import numpy as np
from define_experiment import *
import matplotlib.pyplot as plt

experiment, end_event = experiment1()
alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value='extrapolate')
terrain_dist = np.linspace(np.amin(experiment['alpha_dist']),np.amax(experiment['alpha_dist']),100)
position = alpha_fun(terrain_dist)

plt.plot(terrain_angle,position,color='g',label="interpolated values with cubic fit")
plt.plot(experiment['alpha_dist'],experiment['alpha_deg'],'*',color='red', label='known data values')
plt.title("Terrain Angle vs. Position with cubic approximation")
plt.xlabel("Postion (m)")
plt.ylabel("Terrain Angle (deg)")
plt.legend(loc="best")
plt.show()

'''
Explanation:
    We are cubucally fitting the data points that we are given. Using this fit, we are assuming that
    the terrain in between the known points is changing as a cubic function. Because we are interpolating 
    to find the terrain angles in between all of the known position values, we cannot be certain that the
    terrain does in fact flow in the manner which is found in the figure.
'''

