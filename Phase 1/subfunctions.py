#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Edited: 09/19/2022
@authors: Udit Parikh, Yajat Jotwani
"""
from define_rover import *
from numpy import sin, radians, cos, tan, array, ndarray, zeros, ndim
from math import erf

rover, planet = rover1()
def get_mass(rover):
    if type(rover) != dict:
        raise Exception("Input to get_mass must be a dictionary")
    return((6*rover['wheel_assembly']['wheel']['mass']) + 6*rover['wheel_assembly']['speed_reducer']['mass'] + 6*rover['wheel_assembly']['motor']['mass'] + rover['chassis']['mass'] + rover['science_payload']['mass'] + rover['power_subsys']['mass'])

def get_gear_ratio(speed_reducer):
    '''
    output: Speed ratio from input pinion shaft to output gear shaft. Unitless
    '''
    if type(speed_reducer) != dict:
        raise Exception("Input to get_gear_ratio must be a dictionary")
    return(((speed_reducer['diam_gear']) / (speed_reducer['diam_pinion']))**2)    
    
def tau_dcmotor(omega,rover):
    '''
    Returns the motor shaft torque when given motor shaft speed and a 
    dictionary containing important specifications for the motor.
    '''
    if ndim(omega) != 0 and ndim(omega) != 1:
        raise Exception('omega (Motor shaft speed) must be a scalar or 1D numpy array. No matricies are allowed')
    elif type(rover) != dict:
        raise Exception('Rover properties must be a dictionary')
    tau_s = rover['torque_stall']
    tau_nl = rover['torque_noload']
    omega_nl = rover['speed_noload']
    if ndim(omega) == 0:
        return (tau_s - ((tau_s - tau_nl) / omega_nl) * omega)
    tau = zeros(len(omega))
    for w in range(len(omega)):
        if omega[w] > omega_nl:
            return 0
        elif omega[w] < 0:
            return tau_s
        else:
            tau[w] = (tau_s - ((tau_s - tau_nl) / omega_nl) * omega[w])
    return tau

def F_drive(omega, rover):
    '''
    input: omega and rover
    output: force applied to rover by drive system
    '''
    if ndim(omega) != 0 and ndim(omega) != 1:
        raise Exception('omega (Motor shaft speed) must be a scalar or 1D numpy array. No matricies are allowed')
    elif type(rover) != dict:
        raise Exception('Rover properties must be a dictionary')
    r = rover['wheel_assembly']['wheel']['radius']
    tau_out = get_gear_ratio(rover['wheel_assembly']['speed_reducer']) * tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    if ndim(omega) == 0:
        return 6*tau_out / r
    fDrive = zeros(len(omega))
    for w in range(len(omega)):
        tau_out = get_gear_ratio(rover['wheel_assembly']['speed_reducer']) * tau_dcmotor(omega[w], rover['wheel_assembly']['motor'])
        fDrive[w] = 6 * tau_out / r  
    return fDrive

def F_gravity(terrain_angle, rover, planet):
    '''
    Test Code:
        terrain_angle = array([-5,0,5,10,20,30])
        print(F_gravity(terrain_angle, rover, planet))
    '''
    if type(rover) != dict:
        raise Exception('Rover properties must be a dictionary')
    elif type(planet) != dict:
        raise Exception('Planet must be a dictionary')
    elif ndim(terrain_angle) != 0 and ndim(terrain_angle) != 1:
        raise Exception('terrain angle must be a scalar or 1D numpy array. No matricies are allowed')
    if abs(terrain_angle) > 75:
        raise Exception('Terrain angle must be between -75 and +75 degrees')
    g = planet['g']
    m = get_mass(rover) 
    if ndim(terrain_angle) == 0:
        return(m * g * sin(radians(terrain_angle)))
    Fgt = zeros(len(terrain_angle))
    for t in range(len(terrain_angle)):
         Fgt[t] = m * g * sin(radians(terrain_angle[t]))
    return Fgt
        
def F_rolling(omega, terrain_angle, rover, planet, Crr):
    '''
    if 
    input: anlge w/ horizon, mass of rover, rolling resistance-coe
    output:magnitude of the force component acting on the rover in the direction of its 
    translational motion due to gravity as a function of terrain inclination angle and rover 
    properties.
    
    NOTE: Calculates the rolling resistance over ALL SIX WHEELS

    TEST: 
        omega = array([0,0.5,1,2,3,3.8])
        terrain_angle = array([-5,0,5,10,20,30])
        Crr = 0.2
        print(F_rolling(omega, terrain_angle, rover, planet, Crr))
    '''
    if type(rover) != dict:
        raise Exception('Rover properties must be a dictionary')
    elif type(planet) != dict:
        raise Exception('Planet must be a dictionary')
    elif ndim(terrain_angle) != 0 or ndim(omega) != 0:
        if ndim(terrain_angle) == 1 and ndim(omega) == 1:
            if len(terrain_angle) == len(omega):
                pass
            else:
                raise Exception('terrain_angle and omega must be of the same dimension')
        else:
            raise Exception('terrain angle and omega must be a scalar or 1D numpy array. No matricies are allowed')
    elif Crr <= 0:
        raise Exception('Crr must be a positive scalar')
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    wheel_radius = rover['wheel_assembly']['wheel']['radius']
    roverMass = get_mass(rover)
    g_mars = planet['g']
    if ndim(terrain_angle) == 0:  
       if abs(terrain_angle) > 75:
           raise Exception('Terrain angle must be between -75 and +75 degrees')
       omega_wheel = omega * gear_ratio
       rover_velocity = omega_wheel * wheel_radius
       return(erf(40 * rover_velocity) * (Crr * roverMass * g_mars * cos(radians(terrain_angle)))) # Frr_simple
    Frr = zeros(len(omega))
    for i in range(len(omega)):
        if abs(terrain_angle[i]) > 75:
            raise Exception('All terrain angles must be between -75 and +75 degrees')
        omega_wheel = omega[i] * gear_ratio
        rover_velocity = omega_wheel * wheel_radius
        Frr[i] = erf(40 * rover_velocity) * (Crr * roverMass * g_mars * cos(radians(terrain_angle[i])))
    return Frr

def F_net(omega, terrain_angle, rover, planet, Crr):
    if type(rover) != dict:
        raise Exception('Rover properties must be a dictionary')
    elif type(planet) != dict:
        raise Exception('Planet must be a dictionary')
    elif ndim(terrain_angle) != 0 or ndim(omega) != 0:
        if ndim(terrain_angle) == 1 and ndim(omega) == 1:
            if len(terrain_angle) == len(omega):
                pass
            else:
                raise Exception('terrain_angle and omega must be of the same dimension')
        else:
            raise Exception('terrain angle and omega must be a scalar or 1D numpy array. No matricies are allowed')
    elif Crr <= 0:
        raise Exception('Crr must be a positive scalar')
    if ndim(terrain_angle) == 0:  
       if abs(terrain_angle) > 75:
           raise Exception('Terrain angle must be between -75 and +75 degrees')
       else:
           return F_drive(omega, rover) + F_gravity(terrain_angle, rover, planet) + F_rolling(omega, terrain_angle, rover, planet, Crr)
    Frr = zeros(len(omega))
    for i in range(len(omega)):
        if abs(terrain_angle[i]) > 75:
            raise Exception('All terrain angles must be between -75 and +75 degrees')
        Frr[i] = F_drive(omega[i], rover) + F_gravity(terrain_angle[i], rover, planet) + F_rolling(omega[i], terrain_angle[i], rover, planet, Crr)
    return Frr

