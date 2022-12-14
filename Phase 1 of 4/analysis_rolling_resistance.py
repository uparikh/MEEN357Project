#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from define_rover import *
from subfunctions import *
from numpy import linspace, zeros, array
from scipy.optimize import root_scalar
from matplotlib.pyplot import plot,xlabel,ylabel
from random import uniform

rover, planet = rover1()
terrain_angle = 0
Crr_array = linspace(0.01,0.4,25)
w_max = zeros(len(Crr_array))
v_max = zeros(len(Crr_array))
omega_nl = rover['wheel_assembly']['motor']['speed_noload']
x0 = uniform(0,omega_nl/2)
x1 = uniform(omega_nl/2,omega_nl)

for a in range(len(Crr_array)):
    func_find_root = lambda omega: F_net(omega,terrain_angle, rover, planet, Crr_array[a])
    root = root_scalar(func_find_root, method='secant',x0=x0,x1=x1)
    w_max[a] = root.root
    
v_max = rover['wheel_assembly']['wheel']['radius'] * w_max / get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

plot(Crr_array,v_max)
xlabel('Coefficient Rolling Resistance (CRR)')
ylabel('Max Rover Speed [m/s]')

