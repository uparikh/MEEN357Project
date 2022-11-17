"""###########################################################################
#   This file initializes a rover structure for testing/grading
#
#   Created by: Jonathan Weaver-Rosen
#   Last Modified: 25 June 2021
###########################################################################"""

import numpy as np

def define_rover_1():
    # Initialize Rover dict for testing
    wheel = {'radius':0.30,                 # [m] radius of wheel
             'mass':1}                      # [kg] mass of wheel
    speed_reducer = {'type':'reverted',     # [string] input and output gears are co-axial
                     'diam_pinion':0.04,    # [m] diameter of pinion gear
                     'diam_gear':0.07,      # [m] diameter of main gear
                     'mass':1.5}            # [kg] mass of speed reducer
    motor = {'torque_stall':170,            # [Nm] torque needed to stall motor
             'torque_noload':0,             # [Nm] torque applied by motor with zero load
             'speed_noload':3.80,           # [rad/s] speed of motor with zero load
             'mass':5.0}                    # [kg] mass of motor
    
    
    # phase 2 add ##############################
    motor['effcy_tau'] = np.array([0, 10, 20, 40, 75, 165]) # [Nm] torque data points
    motor['effcy']     = np.array([0,.60,.75,.73,.55, .05]) # [percent] efficiency measurements corresponding to torque data points
    #############################################
    
    
    chassis = {'mass':659}                              # [kg] mass of chassis
    science_payload = {'mass':75}                       # [kg] mass of science payload
    power_subsys = {'mass':90}                          # [kg] mass of power sub-system
    
    wheel_assembly = {'wheel':wheel,                    # [string] wheel assembly has wheel
                      'speed_reducer':speed_reducer,    # [string] wheel assembly has speed reducer
                      'motor':motor}                    # [string] wheel assembly has motor
    
    rover = {'wheel_assembly':wheel_assembly,           # [string] rover has wheel assembly
             'chassis':chassis,                         # [string] rover has chassis
             'science_payload':science_payload,         # [string] rover has science payload
             'power_subsys':power_subsys}               # [string] rover has power sub-system
    
    planet = {'g':3.72}           # [m/s^2] gravity of mars
    
    # end_event = {}
    # end_event['max_distance'] = 50
    # end_event['max_time'] = 2000
    # end_event['min_velocity'] = 0.01
    # end_event['on_surface'] = 'False'
    
    # experiment = {}
    # experiment['crr'] = 0.1
    # experiment['time_range'] = np.array([0, 2000])
    # experiment['initial_conditions'] = np.array([[0.2],[0]])
    # experiment['alpha_deg']  = np.array([0, -3, 3, 10, 5])
    # experiment['alpha_dist'] = np.array([0, 5, 15, 30, 50])
    
    # return everything we need
    return rover, planet

def define_rover_2():
    # Initialize Rover dict for testing
    wheel = {'radius':0.30,                 # [m] radius of wheel
             'mass':2}                      # [kg] mass of wheel
    speed_reducer = {'type':'reverted',     # [string] input and output gears are co-axial
                     'diam_pinion':0.04,    # [m] diameter of pinion gear
                     'diam_gear':0.06,      # [m] diameter of main gear
                     'mass':1.5}            # [kg] mass of speed reducer
    motor = {'torque_stall':180,            # [Nm] torque needed to stall motor
             'torque_noload':0,             # [Nm] torque applied by motor with zero load
             'speed_noload':3.70,           # [rad/s] speed of motor with zero load
             'mass':5.0}                    # [kg] mass of motor
    

    
    # # phase 2 add ##############################
    # motor['effcy_tau'] = [0, 10, 20, 40, 75, 165]     # [Nm] torque data points
    # motor['effcy']     = [0,.60,.75,.73,.55, .05]     # [percent] efficiency measurements corresponding to torque data points
    # #############################################
    

    chassis = {'mass':659}                              # [kg] mass of chassis
    science_payload = {'mass':75}                       # [kg] mass of science payload
    power_subsys = {'mass':90}                          # [kg] mass of power sub-system
    
    wheel_assembly = {'wheel':wheel,                    # [string] wheel assembly has wheel
                      'speed_reducer':speed_reducer,    # [string] wheel assembly has speed reducer
                      'motor':motor}                    # [string] wheel assembly has motor
    

    rover = {'wheel_assembly':wheel_assembly,           # [string] rover has wheel assembly
             'chassis':chassis,                         # [string] rover has chassis
             'science_payload':science_payload,         # [string] rover has science payload
             'power_subsys':power_subsys}               # [string] rover has power sub-system
    
    planet = {'g':3.72}           # [m/s^2] gravity of mars
    
    # end_event = {}
    # end_event['max_distance'] = 50
    # end_event['max_time'] = 2000
    # end_event['min_velocity'] = 0.01
    # end_event['on_surface'] = 'False'
    
    # experiment = {}
    # experiment['crr'] = 0.1
    # experiment['time_range'] = [0, 2000]
    # experiment['initial_conditions'] = [[0.2][0]]
    # experiment['alpha_deg']  = [0, -3, 3, 10, 5]
    # experiment['alpha_dist'] = [0, 5, 15, 30, 50]
    
    # return everything we need
    return rover, planet

def define_rover_3():
    # Initialize Rover dict for testing
    
    wheel = {'radius':0.30,                 # [m] radius of wheel
             'mass':2}                      # [kg] mass of wheel
    speed_reducer = {'type':'standard',     # [string] input and output gears are not co-axial
                     'diam_pinion':0.04,    # [m] diameter of pinion gear
                     'diam_gear':0.06,      # [m] diameter of main gear
                     'mass':1.5}            # [kg] mass of speed reducer
    motor = {'torque_stall':180,            # [Nm] torque needed to stall motor
             'torque_noload':0,             # [Nm] torque applied by motor with zero load
             'speed_noload':3.70,           # [rad/s] speed of motor with zero load
             'mass':5.0}                    # [kg] mass of motor
    
    # # phase 2 add ##############################
    # motor['effcy_tau'] = [0, 10, 20, 40, 75, 165]
    # motor['effcy']     = [0,.60,.75,.73,.55, .05]
    # #############################################
    
    
    chassis = {'mass':659}                              # [kg] mass of chassis
    science_payload = {'mass':75}                       # [kg] mass of science payload
    power_subsys = {'mass':90}                          # [kg] mass of power sub-system
    
    wheel_assembly = {'wheel':wheel,                    # [string] wheel assembly has wheel
                      'speed_reducer':speed_reducer,    # [string] wheel assembly has speed reducer
                      'motor':motor}                    # [string] wheel assembly has motor
    
    rover = {'wheel_assembly':wheel_assembly,           # [string] rover has wheel assembly
             'chassis':chassis,                         # [string] rover has chassis
             'science_payload':science_payload,         # [string] rover has science payload
             'power_subsys':power_subsys}               # [string] rover has power sub-system
    
    planet = {'g':3.72}           # [m/s^2] gravity of mars
    
    # return everything we need
    return rover, planet


def define_rover_4():
    # Initialize Rover dict for testing  
    wheel = {'radius':0.20,                 # [m] radius of wheel
             'mass':2}                      # [kg] mass of wheel
    speed_reducer = {'type':'reverted',     # [string] input and output gears are co-axial
                     'diam_pinion':0.04,    # [m] diameter of pinion gear
                     'diam_gear':0.06,      # [m] diameter of main gear
                     'mass':1.5}            # [kg] mass of speed reducer
    motor = {'torque_stall':165,            # [Nm] torque needed to stall motor
             'torque_noload':0,             # [Nm] torque applied by motor with zero load
             'speed_noload':3.85,           # [rad/s] speed of motor with zero load
             'mass':5.0}                    # [kg] mass of motor
    
    # phase 2 add ##############################
    motor['effcy_tau'] = np.array([0, 10, 20, 40, 75, 170]) # [Nm] torque data points
    motor['effcy']     = np.array([0,.60,.75,.73,.55, .05]) # [percent] efficiency measurements corresponding to torque data points
    ############################################# 
    
    chassis = {'mass':674}                              # [kg] mass of chassis
    science_payload = {'mass':80}                       # [kg] mass of science payload
    power_subsys = {'mass':100}                          # [kg] mass of power sub-system
    
    wheel_assembly = {'wheel':wheel,                    # [string] wheel assembly has wheel
                      'speed_reducer':speed_reducer,    # [string] wheel assembly has speed reducer
                      'motor':motor}                    # [string] wheel assembly has motor
    
    rover = {'wheel_assembly':wheel_assembly,           # [string] rover has wheel assembly
             'chassis':chassis,                         # [string] rover has chassis
             'science_payload':science_payload,         # [string] rover has science payload
             'power_subsys':power_subsys}               # [string] rover has power sub-system
    
    planet = {'g':3.72}           # [m/s^2] gravity of mars
    
    end_event = {'max_distance' : 10,       # [m?] maximum distance before stopping
                 'max_time' : 10000,        # [s] maximum time before stopping
                 'min_velocity' : 0.01}     # [m/s] minimum velocity rover moves
    
    # return everything we need
    return rover, planet #, end_event