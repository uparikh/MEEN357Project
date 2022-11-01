from scipy.interpolate import *
from subfunctions import *
from define_experiment import *
import matplotlib.pyplot as plt
from define_rover_phase2 import *
import numpy as np



experiment, end_event = experiment1()


end_events = end_event
simulate_rover(rover, planet, experiment, end_events)

# Define variables for graphing
times = rover['telemetry']['Time']
positions = rover['telemetry']['position']
velos = rover['telemetry']['velocity']
pows = rover['telemetry']['power']


fig, ax = plt.subplots(3, 1)
fig.tight_layout(h_pad=4)

# # graph the first equation on the top plot
ax[0].plot(times, positions, 'bo', times, positions, 'b')
ax[0].set(xlabel='Time (s)', ylabel='Position (m)',title='Position vs. Time')

# # graph the second equation on the bottom plot
ax[1].plot(times, velos, 'ro', times, velos, 'r')
ax[1].set(xlabel='Time (s)', ylabel='Velocity (m/s)',title='Velocity vs. Time')

# # graph the third equation on the bottom plot
ax[2].plot(times, pows, 'go', times, pows, 'g')
ax[2].set(xlabel='Time (s)', ylabel='Power (W)',title='Power vs. Time')

# rover, planet = rover1()
# effcy_fun = interp1d(rover['wheel_assembly']['motor']['effcy_tau'], rover['wheel_assembly']['motor']['effcy'], kind = 'cubic')
# rover_effcy = np.linspace(np.amin(rover['wheel_assembly']['motor']['effcy_tau']),np.amax(rover['wheel_assembly']['motor']['effcy_tau']),100)
# efficiency = effcy_fun(rover_effcy)

# plt.plot(rover_effcy,efficiency,color='g',label="interpolated values with cubic fit")
# plt.plot(rover['wheel_assembly']['motor']['effcy_tau'],rover['wheel_assembly']['motor']['effcy'],'r*', label='known data values')
# plt.title("Motor Torque vs. Efficiency")
# plt.xlabel("Torque (N*m)")
# plt.ylabel("Efficiency")
# plt.legend(loc="best")
# plt.show()



# X1 = np.linspace(-5, 5, 20)
# X2 = np.linspace(0, 100, 20)
# X3 = np.linspace(1, 20, 20)

# # create 2 plots that are well spaced apart
# fig, ax = plt.subplots(3, 1)
# fig.tight_layout(h_pad=4)
    
# # graph the first equation on the top plot
# ax[0].plot(X1, 2*(X1), 'bo', X1, 2*(X1), 'b')
# ax[0].set(xlabel='x', ylabel='y',title='Exhibit 1')

# # graph the second equation on the bottom plot
# ax[1].plot(X2, 0.5*(X2), 'ro', X2, 0.5*(X2), 'r')
# ax[1].set(xlabel='work', ylabel='play',title='Exhibit 2')

# # graph the third equation on the bottom plot
# ax[2].plot(X3, 0.5*(X3), 'ro', X3, 0.5*(X3), 'g')
# ax[2].set(xlabel='work', ylabel='play',title='Exhibit 2')


# plt.show()