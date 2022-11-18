# Task 5: Evaluate the Impact of Changing the Parachute Size 
# i need help determining the time at termination of simulation
from matplotlib.pyplot import plot, xlabel, ylabel, subplots, subplots_adjust, figure
from numpy import linspace, zeros, multiply, max, where
from define_edl_system import *
from subfunctions_EDL import *
from define_planet import *
from define_mission_events import *

edl_system = define_edl_system_1()
rocket = define_edl_system_1()
heat_shield = define_edl_system_1()
sky_crane = define_edl_system_1()
speed_controller = define_edl_system_1()
position_controller = define_edl_system_1()



mars = define_planet()
mission_events = define_mission_events()

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
# edl_system['parachute']['diameter'] = ?


diameters = linspace(14, 19, 11)
tmax = 2000   # [s] maximum simulated time


[t, Y, edl_system] = simulate_edl(edl_system, mars, mission_events, tmax, True)

plot1 = figure(0)
fig, axs = subplots(3) 
axs[0].plot(diameters, t)
axs[0].set_title('time vs. parachute diameter')




# subplot(3,1,1)
# plot(x,y)






'''
from matplotlib.pyplot import plot, xlabel, ylabel, subplot, subplots_adjust
from define_rover import *
from subfunctions import tau_dcmotor
from numpy import linspace, zeros, multiply, max, where

rover, planet = rover1()
omega = linspace(0,3.8,20) #increase last input to make more smooth

subplot(3,1,1)
tau = zeros(len(omega))
for w in range(len(omega)):   
    tau[w] = tau_dcmotor(omega[w], rover['wheel_assembly']['motor'])
plot(tau,omega)
xlabel('Motor Shaft Torque [Nm]')
ylabel('Motor Shaft Speed [rad/s]')

subplot(3,1,2)
power = multiply(omega, tau)
plot(tau,power)
xlabel('Motor Shaft Torque [Nm]')
ylabel('Motor Power [W]')

subplot(3,1,3)
plot(omega,power)
xlabel('Motor Shaft Speed [rad/s]')
ylabel('Motor Power [W]')

subplots_adjust(bottom=-1.0)
'''