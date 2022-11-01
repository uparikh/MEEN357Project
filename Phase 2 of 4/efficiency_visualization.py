# Task 5

from scipy.interpolate import *
import numpy as np
from define_experiment import *
import matplotlib.pyplot as plt
from define_rover_phase2 import *

rover, planet = rover1()
effcy_fun = interp1d(rover['wheel_assembly']['motor']['effcy_tau'], rover['wheel_assembly']['motor']['effcy'], kind = 'cubic')
rover_effcy = np.linspace(np.amin(rover['wheel_assembly']['motor']['effcy_tau']),np.amax(rover['wheel_assembly']['motor']['effcy_tau']),100)
efficiency = effcy_fun(rover_effcy)

plt.plot(rover_effcy,efficiency,color='g',label="interpolated values with cubic fit")
plt.plot(rover['wheel_assembly']['motor']['effcy_tau'],rover['wheel_assembly']['motor']['effcy'],'r*', label='known data values')
plt.title("Motor Torque vs. Efficiency")
plt.xlabel("Torque")
plt.ylabel("Efficiency")
plt.legend(loc="best")
plt.show()

'''
Explanation:
    We are cubucally fitting the data points that we are given. Using this fit, we are assuming that
    the efficiency in between the known points is changing as a cubic function. Because we are interpolating 
    to find the efficiencies in between all of the known torque values, we cannot be certain that the
    efficiency does in fact flow in the manner which is found in the figure.
'''