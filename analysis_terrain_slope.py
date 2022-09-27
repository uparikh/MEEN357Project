#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from define_rover import *
from subfunctions import F_net
from numpy import linspace, zeros, array,ndim
from scipy.optimize import root_scalar
from matplotlib.pyplot import plot,xlabel,ylabel
from random import uniform

rover, planet = rover1()
Crr = 0.2
slope_array_deg = linspace(-10,35,25)
v_max = zeros(len(slope_array_deg), dtype = float)
omega_nl = rover['wheel_assembly']['motor']['speed_noload'] #* (rover['wheel_assembly']['wheel']['radius'])
x0 = uniform(0,omega_nl/10)
x1 = uniform(omega_nl/10,omega_nl)

for a in range(len(slope_array_deg)):
    func_find_root = lambda omega: F_net(omega, float(slope_array_deg[a]), rover, planet, Crr)
    # v_max[a] = (root_scalar(func_find_root, method='secant',x0=0.012,x1=0.03).root) 
    v_max[a] = (root_scalar(func_find_root, method='bisect',bracket=[0,3.8])).root 
    print(v_max[a])
plot(slope_array_deg,v_max)
xlabel('Terrain Angle [deg]')
ylabel('Max Rover Speed [m/s]')


