import matplotlib.pyplot as plt
import numpy as np
from define_rovers import *
from define_edl_system import *
from define_planet import *
from define_mission_events import *
from subfunctions_EDL import *
 
edl_system = define_edl_system_1()
rocket = define_edl_system_1()
heat_shield = define_edl_system_1()
sky_crane = define_edl_system_1()
speed_controller = define_edl_system_1()
position_controller = define_edl_system_1()
planet = define_planet()
mission_events = define_mission_events()
tmax = 2000

diameters = np.arange(14,19.5,0.5)
times = np.array([])
v = np.array([])
fail = np.array([])
for d in diameters:
    edl_system = define_edl_system_1()
    # initial conditions:
    edl_system['altitude'] = 11000    # [m] initial altitude
    edl_system['velocity'] = -578     # [m/s] initial velocity
    rocket['on'] = False              # rocket is off
    edl_system['parachute']['deployed'] = True   # our parachute is open
    edl_system['parachute']['ejected'] = False   # and still attached
    heat_shield['ejected'] = False               # heat shield is not ejected
    sky_crane['on'] = False                      # sky crane is off
    speed_controller['on'] = False               # speed controller is off
    position_controller['on'] = False            # position controller is off
    edl_system['rover']['on_ground'] = False # the rover has not yet landed
    edl_system['parachute']['diameter'] = d
    [t, Y, edl_system] = simulate_edl(edl_system, planet, mission_events, tmax, True)
    times = np.append(times,t[len(t)-1]) 
    v = np.append(np.abs(v),Y[0][len(Y[0])-1])
    if Y[0][len(Y[0])-1] <= edl_system['sky_crane']['danger_speed']:
        fail = np.append(fail,1)
    else:
        fail = np.append(fail,0)
plt.subplot(3,1,1)
plt.plot(diameters, times)
plt.xlabel('Diameter (m)')
plt.ylabel('Time (s)')

plt.subplot(3,1,2)
plt.plot(diameters, v)
plt.xlabel('Diameter (m)')
plt.ylabel('Velocity (m/s)')

plt.subplot(3,1,3)
plt.plot(diameters,fail)
plt.xlabel('Diameter (m)')
plt.ylabel('Success/Failure')

plt.subplots_adjust(bottom=-1.2)
plt.show()








































