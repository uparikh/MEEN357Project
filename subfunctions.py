#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Edited: 09/19/2022
@authors: Udit Parikh, Yajat Jotwani
"""

# delete test later
#print(rover['wheel_assembly']['wheel']['radius'])

def get_mass(rover):
    if type(rover) != dict:
        raise Exception("Input to get_mass must be a dictionary")
    return((6*rover['wheel_assembly']['wheel']['mass']) + rover['wheel_assembly']['speed_reducer']['mass'] 
           + rover['wheel_assembly']['motor']['mass'] + rover['chassis']['mass'] + rover['science_payload']['mass'] + rover['power_subsys']['mass'])

def get_gear_ratio(speed_reducer):
    if type(speed_reducer) != dict or speed_reducer['type'].lower != 'reverted':
        raise Exception("Input to get_gear_ratio must be a dictionary")
    return(((rover['wheel_assembly']['speed_reducer']['diam_gear']) / (rover['wheel_assembly']['speed_reducer']['diam_pinion']))**2)

# delte test later
#print(get_gear_ratio(rover))
    
    
    
def tau_dcmotor(rover):
    '''
    Returns the motor shaft torque when given motor shaft speed and a 
    dictionary containing important specifications for the motor.
    '''
    return(rover['wheel_assembly']['motor']['torque_stall'])

print(tau_dcmotor(rover))

    


