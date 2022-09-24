#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Edited: 09/19/2022
@authors: Udit Parikh, Yajat Jotwani
"""
from define_rover import *
from numpy import sin, radians, cos, tan, array, ndarray, zeros, ndim

rover, planet = rover1()
def get_mass(rover):
    if type(rover) != dict:
        raise Exception("Input to get_mass must be a dictionary")
    return((6*rover['wheel_assembly']['wheel']['mass']) + 6*rover['wheel_assembly']['speed_reducer']['mass'] + 6*rover['wheel_assembly']['motor']['mass'] + rover['chassis']['mass'] + rover['science_payload']['mass'] + rover['power_subsys']['mass'])

def get_gear_ratio(speed_reducer):
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

omega = 0
print(F_drive(omega,rover))

def F_gravity(terrain_angle, rover, planet):
    roverMass = get_mass(rover)
    g_mars = planet['g']
    return(roverMass * g_mars * sin(radians(terrain_angle)))

