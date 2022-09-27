#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from define_rover import *
from subfunctions import F_net
from numpy import linspace, zeros, array,ndim,append,NaN
from scipy.optimize import root_scalar
from matplotlib.pyplot import plot,xlabel,ylabel
from random import uniform

rover, planet = rover1()
Crr = 0.2
slope_array_deg = linspace(-10,35,25)
w_max = zeros(len(slope_array_deg), dtype = float)
v_max = zeros(len(slope_array_deg), dtype = float)
omega_nl = rover['wheel_assembly']['motor']['speed_noload']
x0 = uniform(0,(omega_nl-2))
x1 = uniform((omega_nl-2),(omega_nl+2))

# for a in slope_array_deg:
#     # func_find_root = lambda omega: F_net(omega, float(slope_array_deg[a]), rover, planet, Crr)
#     w_max[a] = root_scalar(F_net,(a, rover, planet, Crr),method='bisect',bracket=[0,omega_nl]).root
#     # v_max[a] = (root_scalar(func_find_root, method='bisect',bracket=[0,omega_nl]).root)


for angle in slope_array_deg:
    try:
        omega = root_scalar(F_net, (angle, rover, planet, Crr), method='bisect', bracket=[0, omega_nl]).root
    except ValueError:
        omega = NaN
    w_max = append(w_max,omega)
    # print(w_max)
    
plot(slope_array_deg,w_max)


xlabel('Terrain Angle [deg]')
ylabel('Max Rover Speed [m/s]')


