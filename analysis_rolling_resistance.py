from define_rover import *
from subfunctions import *
from numpy import linspace, zeros, array
from scipy.optimize import root_scalar
from matplotlib.pyplot import plot,xlabel,ylabel
from random import uniform

rover, planet = rover1()
Crr = 0.2
slope_array_deg = linspace(-10,35,25)
v_max = zeros(len(slope_array_deg))
omega_nl = rover['wheel_assembly']['motor']['speed_noload']
x0 = uniform(0,omega_nl/2)
x1 = uniform(omega_nl/2,omega_nl)

for a in range(len(slope_array_deg)):
    func_find_root = lambda omega: F_net(omega, float(slope_array_deg[a]), rover, planet, Crr)
    root = root_scalar(func_find_root, method='secant',x0=x0,x1=x1)
    v_max[a] = root.root
    
plot(slope_array_deg,v_max)
xlabel('Terrain Angle [deg]')
ylabel('Max Rover Speed [m/s]')