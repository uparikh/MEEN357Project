'''
output: determine the speed of the rover at various values for the 
coefficient of rolling resistance AND terrain slope.


'''
from define_rover import *
from subfunctions import *
from numpy import linspace, zeros, array
from scipy.optimize import root_scalar
from matplotlib.pyplot import plot,xlabel,ylabel, figure,title
from random import uniform
from numpy import linspace, meshgrid, zeros, shape
from mpl_toolkits.mplot3d import Axes3D


slope_array_deg = linspace(-10,35,25)
Crr_array = linspace(0.01,0.4,25)
omega_nl = rover['wheel_assembly']['motor']['speed_noload']
x0 = uniform(0,omega_nl/2)
x1 = uniform(omega_nl/2,omega_nl)

CRR, SLOPE = meshgrid(Crr_array, slope_array_deg)
VMAX = zeros(shape(CRR), dtype = float)
N = shape(CRR)[0]

for i in range(N):
    for j in range(N):
         #if (SLOPE[i,j] < 0) and (abs(CRR[i,j]) < 1e-5):
             #raise Exception('NAN')
         
         #else:
            Crr_sample = float(CRR[i,j])
            slope_sample = float(SLOPE[i,j])
            
            # here you put code to find the max speed at Crr_sample and slope_sample
            func_find_root = lambda omega: F_net(omega, SLOPE[i,j], rover, planet, CRR[i,j])
            root = root_scalar(func_find_root, method='secant',x0=x0,x1=x1)
    
            VMAX[i,j] = root.root
        
figure = figure()
ax = Axes3D(figure, auto_add_to_figure=False)
ax.plot_surface(CRR, SLOPE, VMAX)


figure.add_axes(ax)
xlabel('Coefficient Rolling Resistance (CRR)')
ylabel('Terrain Angle [deg]')
ax.set_zlabel('Max Rover Speed [m/s]')

title('Max Velocity vs. CRR and Terrain Slope')