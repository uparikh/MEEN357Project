#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Edited: 09/19/2022
@authors: Udit Parikh, Yajat Jotwani
"""
from define_rover import *
from numpy import sin, radians, cos, tan, array, zeros

rover, planet = rover1()
def get_mass(rover):
    if type(rover) != dict:
        raise Exception("Input to get_mass must be a dictionary")
    return((6*rover['wheel_assembly']['wheel']['mass']) + 6*rover['wheel_assembly']['speed_reducer']['mass'] + 6*rover['wheel_assembly']['motor']['mass'] + rover['chassis']['mass'] + rover['science_payload']['mass'] + rover['power_subsys']['mass'])

def get_gear_ratio(speed_reducer):
    if type(speed_reducer) != dict:
        raise Exception("Input to get_gear_ratio must be a dictionary")
    return(((rover['wheel_assembly']['speed_reducer']['diam_gear']) / (rover['wheel_assembly']['speed_reducer']['diam_pinion']))**2)    
    
def tau_dcmotor(omega,rover):
    '''
    Returns the motor shaft torque when given motor shaft speed and a 
    dictionary containing important specifications for the motor.
    '''
    if type(omega) != ndarray:
        raise Exception('omega (Motor shaft speed) must be a numpy array')
    elif type(rover) != dict:
        raise Exception('Rover properties must be a dictionary')
   
    tau_s = rover['torque_stall']
    tau_nl = rover['torque_noload']
    omega_nl = rover['speed_noload']
    tau = zeros(len(omega))
    for i in range(len(omega)):
        if omega[i] > omega_nl:
            return 0
        elif omega[i] < 0:
            return tau_s
    return (tau_s - ((tau_s - tau_nl) / omega_nl) * omega)

omega = array([0.00,0.50,1.00,2.00,3.00,3.80])
print(tau_dcmotor(omega,rover['wheel_assembly']['motor'])

# def F_drive(omega, rover):
    # '''
    # input: omega and rover
    # output: force applied to rover by drive system
    # '''
    # tau_in = tau_dcmotor(omega, rover)
    # tau_out = get_gear_ratio(rover) * tau_in
    # r = rover['wheel_assembly']['wheel']['radius']
    # return  6 * tau_out / r

# def F_gravity(terrain_angle, rover, planet):
    # roverMass = get_mass(rover)
    # g_mars = planet['g']
    # return(roverMass * g_mars * sin(radians(terrain_angle)))


