"""###########################################################################
#   This file contains subfunctions for Phase 1 of the TAMU MEEN 357 project
#
#   Created by: MEEN 357 Rover Physics Team
#   Last Modified: 19 September 2022
###########################################################################"""

import math
import numpy as np
from define_rover_phase2 import *
from define_experiment import *
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import matplotlib.pyplot as plt

rover,planet = rover1()
experiment, end_event = experiment1()

def get_mass(rover):
    """
    Inputs:  rover:  dict      Data structure containing rover parameters
    
    Outputs:     m:  scalar    Rover mass [kg].
    """
    
    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Input must be a dict')
    
    # add up mass of chassis, power subsystem, science payload, 
    # and components from all six wheel assemblies
    m = rover['chassis']['mass'] \
        + rover['power_subsys']['mass'] \
        + rover['science_payload']['mass'] \
        + 6*rover['wheel_assembly']['motor']['mass'] \
        + 6*rover['wheel_assembly']['speed_reducer']['mass'] \
        + 6*rover['wheel_assembly']['wheel']['mass'] \
    
    return m


def get_gear_ratio(speed_reducer):
    """
    Inputs:  speed_reducer:  dict      Data dictionary specifying speed
                                        reducer parameters
    Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
                                        to output gear shaft. Unitless.
    """
    
    
    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Input must be a dict')
    
    # Check 'type' field (not case sensitive)
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The speed reducer type is not recognized.')
    
    # Main code
    d1 = speed_reducer['diam_pinion']
    d2 = speed_reducer['diam_gear']
    
    # Computes gear ratio using pinion and gear diameters
    Ng = (d2/d1)**2
    
    return Ng


def tau_dcmotor(omega, motor):
    """
    Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
             motor:  dict             Data dictionary specifying motor parameters
    Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
                                      is same size as first input argument.
    """

    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(motor) != dict:
        raise Exception('Second input must be a dict')
        
    # Main code
    tau_s    = motor['torque_stall']
    tau_nl   = motor['torque_noload']
    omega_nl = motor['speed_noload']
    
    # initialize
    tau = np.zeros(len(omega),dtype = float)
    for ii in range(len(omega)):
        if omega[ii] >= 0 and omega[ii] <= omega_nl:
            tau[ii] = tau_s - (tau_s-tau_nl)/omega_nl *omega[ii]
        elif omega[ii] < 0:
            tau[ii] = tau_s
        elif omega[ii] > omega_nl:
            tau[ii] = 0
        
    return tau
    
    


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                    planet:  dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance coefficient
                                              [-]
    
    Outputs:           Frr:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    v_rover = r*omega/Ng
    
    Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
    Frr_simple = -Crr*Fn # simple rolling resistance
    
    Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    """
    Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
                     rover:  dict          Data structure specifying rover 
                                            parameters
                    planet:  dict          Data dictionary specifying planetary 
                                            parameters
    
    Outputs:           Fgt:  numpy array   Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that values of the first input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Check that the third input is a dict
    if type(planet) != dict:
        raise Exception('Third input must be a dict')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    
    Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
    return Fgt


def F_drive(omega, rover):
    """
    Inputs:  omega:  numpy array   Array of motor shaft speeds [rad/s]
             rover:  dict          Data dictionary specifying rover parameters
    
    Outputs:    Fd:  numpy array   Array of drive forces [N]
    """
    
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Main code
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    tau_out = tau*Ng
    
    r = rover['wheel_assembly']['wheel']['radius']
    
    # Drive force for one wheel
    Fd_wheel = tau_out/r 
    
    # Drive force for all six wheels
    Fd = 6*Fd_wheel
    
    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict     Data structure specifying rover 
                                      parameters
                     planet:  dict     Data dictionary specifying planetary 
                                      parameters
                        Crr:  scalar   Value of rolling resistance coefficient
                                      [-]
    
    Outputs:           Fnet:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray,)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
    
    # Main Code
    Fd = F_drive(omega, rover)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fg = F_gravity(terrain_angle, rover, planet)
    
    Fnet = Fd + Frr + Fg # signs are handled in individual functions
    
    return Fnet

def motorW(v, rover): 
    '''
    Compute the rotational speed of the motor shaft [rad/s] given the 
    translational velocity of the rover and the rover dictionary.
    
    inputs: Rover translational velocity (as a scalar or a numpy array), Data structure
    containing rover parameters (as a dictionary)
    
    output: Motor speed in rad/s. 
    
    omega_motor = (V_wheel * d_pinion) / (r_wheel * d_gear)
    
    '''
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, (np.ndarray,np.floating,np.integer))):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(v, np.ndarray):
        v = np.array([v],dtype=float) # make the scalar a numpy array
    elif len(np.shape(v)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    r_wheel = rover['wheel_assembly']['wheel']['radius']
    w = np.zeros(len(v))
    for i in range(len(v)):
        # compute rotational speed of motor using velcoity of rover, gear ratio, and radius of wheel
        w[i] = ((v[i] * gear_ratio) / (r_wheel))
    return w

##motorW TEST##
# v = np.array([0.1,0.3])
# v = 0.3
# print(motorW(v,rover))
    

def rover_dynamics(t, y, rover, planet, experiment): 
    '''
    computes the derivative of the state vector: [velocity, position] for 
    the rover given its current state. It requires rover and experiment 
    dictionary input parameters. It is intended to be passed to an ODE 
    solver.
    
    Inputs: time sample in seconds(scalar)Two element state vector in the form [rover velocity, rover position],
    Data structure containing rover definition, Data structure containing planet definition,
    Data structure containing experiment definition
    
    Output: Two element array first derivative of state vector in the form [rover acceleration, rover velocity]
    '''
    if (type(t) != np.int64 and type(t) != np.float64 and type(t) != float and type(t) != int):
        raise Exception("first input must be a scalar")
    elif (not isinstance(y,np.ndarray)):
        raise Exception("second input must be an array of length 2")
    try:
            (len(y) != 2)
    except TypeError:
            raise Exception("second input must be an array of length 2")
    if len(y) != 2:
            raise Exception("second input must be an array of length 2")
    elif type(rover) != dict:
        raise Exception('Third input must be a dict')
    elif type(planet) != dict:
        raise Exception('Third input must be a dict')
    elif type(experiment) != dict:
        raise Exception('Third input must be a dict')
    if (type(y) == np.ndarray):   
        y = y.tolist() # Added to convert number from <class 'numpy.float64'> to <class 'float'>
    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value='extrapolate')
    terrain_angle = (float(alpha_fun(y[1])))
    Crr = experiment['Crr']
    mass = get_mass(rover)
    mW = motorW(y[0], rover)   
    
    Fn = F_net(mW, terrain_angle, rover, planet, Crr)  
    dydt = np.zeros(2)
    dydt[1] =y[0]
    dydt[0] = Fn/mass
    return dydt

##Rover_Dynamics TEST##
# t = 0
# y = np.array([0.33,0])
# print(rover_dynamics(t, y, rover, planet, experiment))


def mechpower(v, rover):
    '''
    Input: translational speed of the rover (as a numpy array) and the 
    characteristics of a specific rover (as a dictionary)
    
    Output: computed instantaneous mechanical power output by a single DC motor 
    at each point in a given velocity profile in Watts (numpy array)
    P(t) = tau_motor(t) * omega_motor(t)

    '''
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(v, np.ndarray):
        v = np.array([v],dtype=float) # make the scalar a numpy array
    elif len(np.shape(v)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
    if (type(rover) != dict):
        raise Exception('Second input must be a dictionary')
    
    radius = rover['wheel_assembly']['wheel']['radius']
    P = np.zeros(len(v))
    for i in range(len(v)):
        mW = motorW(v[i], rover)
        P[i] = tau_dcmotor(mW, rover['wheel_assembly']['motor']) *  mW
    return P

##mechower TEST##
# print(mechpower(np.array([0.05,0.25]),rover))

def battenergy(t,v,rover): #### NOT DONE ####
    '''
    Inputs: t | 1D numpy array | N-element array of time samples from a rover simulation [s]
            v | 1D numpy array | N-element array of rover velocity data from a simulation [m/s]
            rover | dictionary | Data structure containing rover definition
    Outputs: E | scalar | Total electrical energy consumed from the rover battery pack over the input simulation profile [J]
    '''
    if (not isinstance(v, np.ndarray)):
        raise Exception('Second input must be a vector. If input is a vector, it should be defined as a numpy array.')
    if len(np.shape(v)) != 1:
        raise Exception('Second input must be vector. Matrices are not allowed.')
    if (not isinstance(t, np.ndarray)):
        raise Exception('First input must be a vector. If input is a vector, it should be defined as a numpy array.')
    if len(np.shape(t)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
    if (len(t) != len(v)):
        raise Exception('First input and Second input must be arrays of same length')
    if (type(rover) != dict):
        raise Exception('Second input must be a dictionary')
    effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
    effcy = rover['wheel_assembly']['motor']['effcy']
    mW = motorW(v, rover)
    tauDCM = tau_dcmotor(mW, rover['wheel_assembly']['motor'])
    mP = mechpower(v, rover)
    effcy_fun = interp1d(effcy_tau, effcy,kind='cubic')
    effcyFuncTau = effcy_fun(tauDCM)
    power = np.zeros(len(mP))
    for j in range(len(mP)):
        power[j] = (mP[j] * 6)/effcyFuncTau[j]
    E = integrate.trapz(power, t)
    return E

# t = np.array([0,1,2,3,4,5,6])
# v = np.array([0.33,0.32,0.33,0.2,0.2,0.25,0.28])
# print(battenergy(t, v, rover))

def end_of_mission_event(end_event):
    """
    Defines an event that terminates the mission simulation. Mission is over
    when rover reaches a certain distance, has moved for a maximum simulation 
    time or has reached a minimum velocity.            
    """
    
    mission_distance = end_event['max_distance']
    mission_max_time = end_event['max_time']
    mission_min_velocity = end_event['min_velocity']
    
    # Assume that y[1] is the distance traveled
    distance_left = lambda t,y: mission_distance - y[1]
    distance_left.terminal = True
    
    time_left = lambda t,y: mission_max_time - t
    time_left.terminal = True
    
    velocity_threshold = lambda t,y: y[0] - mission_min_velocity;
    velocity_threshold.terminal = True
    velocity_threshold.direction = -1
    
    # terminal indicates whether any of the conditions can lead to the
    # termination of the ODE solver. In this case all conditions can terminate
    # the simulation independently.
    
    # direction indicates whether the direction along which the different
    # conditions is reached matter or does not matter. In this case, only
    # the direction in which the velocity treshold is arrived at matters
    # (negative)
    
    events = [distance_left, time_left, velocity_threshold]
    
    return events

def simulate_rover(rover, planet,experiment,end_events):
    '''
    Inputs: 
        rover | dictionary | Data structure containing rover definition
        planet | dictionary | Data structure containing planet definition
        experiment | dictionary | Data structure containing parameters of the
        trajectory to be followed by the simulation of rover dynamics
        end_event | dictionary | Data structure containing the conditions necessary
        and sufficient to terminate simulation of rover dynamics
        
    Output:
        rover | dictionary | Data structure containing the parameters of the rover, 
        including updated telemetry information
    '''
    if (type(rover) != dict):
        raise Exception('First input must be a dictionary')
    elif (type(planet) != dict):
        raise Exception('Second input must be a dictionary')
    elif (type(experiment) != dict):
        raise Exception('Third input must be a dictionary')
    elif (type(end_events) != dict):
        raise Exception('Fourth input must be a dictionary')
 
    '''
    Rover{'telemetry'}
        Time | 1D array | n-element array containing the time history of the rover [s]
        completion_time | scalar | time to complete a mission [s]
        velocity | 1D array | array containing the velocty of the rover as it follows a trajectory
        position | 1D array | array containing the position of the rover as it follows a trajectory
        distance_traveled | scalar | total distace traveled of the rover [m]
        max_velocty | scalar | maximum velocity of rover along the given trajectory [m/s]
        average_velocity | scalar | average velocity of rover along a given trajectory [m/s]
        power |  1D numpy array | array containing the instantaneous power outputted by the motor along the trajectory [W]
        battery_energy | scalar | Total energy to be extracted from the battery to complete trajectory [J]
        energy_per_distance | scalar | Total energy spent from battery per meter traveled [J/m]
    '''
    terrainFun = lambda t,y : rover_dynamics(t, y, rover, planet, experiment)
    sol = integrate.solve_ivp(terrainFun, np.array([experiment['time_range'][0],end_event['max_time']]),experiment['initial_conditions'],method='BDF',events=end_of_mission_event(end_event))
    avgV = np.average(sol.y[0])
    por = mechpower(sol.y[0], rover)
    distanceTrav = sol.y[1][len(sol.y[1])-1]
    batenergy = battenergy(sol.t, sol.y[0], rover)
    energyDistance = batenergy/distanceTrav
    rover['telemetry'] = {'Time':sol.t, 
                          'completion_time':sol.t[len(sol.t)-1],
                          'velocity':sol.y[0],
                          'position':sol.y[1],
                          'distance_traveled':sol.y[1][len(sol.y[1])-1],
                          'max_velocity':np.max(sol.y[0]),
                          'average_velocity':avgV,
                          'power':por,
                          'battery_energy':batenergy,
                          'energy_per_distance':energyDistance}              
                          
    return rover
# end_events = end_event
# f = simulate_rover(rover, planet, experiment, end_events)
# print(f)








