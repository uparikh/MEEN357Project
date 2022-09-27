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
from numpy import linspace, meshgrid, zeros, shape, NaN
from mpl_toolkits.mplot3d import Axes3D


slope_array_deg = linspace(-15,45,25)
Crr_array = linspace(0.01,0.4,25)
omega_nl = rover['wheel_assembly']['motor']['speed_noload']
x0 = uniform(0,omega_nl/2)
x1 = uniform(omega_nl/2,omega_nl)

CRR, SLOPE = meshgrid(Crr_array, slope_array_deg)
VMAX = zeros(shape(CRR), dtype = float)
N = shape(CRR)[0]

for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i,j])
        slope_sample = float(SLOPE[i,j])
        try:
            VMAX[i,j] = root_scalar(F_net,(slope_sample,rover,planet,Crr_sample), method='secant', x0=x0, x1=x1).root
        except ValueError:
            VMAX[i,j] = NaN
        
figure = figure()
ax = Axes3D(figure, auto_add_to_figure=False)
ax.plot_surface(CRR, SLOPE, VMAX)


figure.add_axes(ax)
ax.set_xlabel('Coefficient Rolling Resistance (CRR)')
ax.set_ylabel('Terrain Angle [deg]')
ax.set_zlabel('Max Rover Speed [m/s]')

title('Max Velocity vs. CRR and Terrain Slope')

