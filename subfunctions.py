#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Edited: 09/19/2022
@authors: Udit Parikh, Yajat Jotwani
"""
from define_rover import *
# delete test later
#print(rover['wheel_assembly']['wheel']['radius'])

rover, planet = rover1()
def get_mass(rover):
    if type(rover) != dict:
        raise Exception("Input to get_mass must be a dictionary")
    return((6*rover['wheel_assembly']['wheel']['mass']) + rover['wheel_assembly']['speed_reducer']['mass'] 
           + rover['wheel_assembly']['motor']['mass'] + rover['chassis']['mass'] + rover['science_payload']['mass'] + rover['power_subsys']['mass'])

#  or speed_reducer['type'].lower() != 'reverted'
def get_gear_ratio(speed_reducer):
    if type(speed_reducer) != dict:
        raise Exception("Input to get_gear_ratio must be a dictionary")
    return(((rover['wheel_assembly']['speed_reducer']['diam_gear']) / (rover['wheel_assembly']['speed_reducer']['diam_pinion']))**2)

# print(get_gear_ratio(rover))
    
    
def tau_dcmotor(omega,rover):
    '''
    Returns the motor shaft torque when given motor shaft speed and a 
    dictionary containing important specifications for the motor.
    '''
    #return(rover['wheel_assembly']['motor']['torque_stall'])
    tau_s = rover['wheel_assembly']['motor']['torque_stall']
    tau_nl = rover['wheel_assembly']['motor']['torque_noload']
    omega_nl = rover['wheel_assembly']['motor']['speed_noload']
    return tau_s - ((tau_s - tau_nl) / omega_nl) * omega

print(tau_dcmotor(3, rover))

    


