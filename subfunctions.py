#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last Edited: 09/19/2022
@authors: Udit Parikh, Yajat Jotwani
"""
rover = {
    'wheel_assembly' : {
        'wheel' : {
            'radius' : 0.30, #m
            'mass' : 1.0 #kg
            },
        'speed_reducer' : {
            'type' : "reverted", #string
            'diam_pinion' : 0.04, #m
            'diam_gear' : 0.07, #m
            'mass' : 1.5 #kg
            },
        'motor' : {
            'torque_stall' : 170, #Nm
            'torque_noload' : 0, #Nm
            'speed_noload' : 3.80, #rad/s
            'mass' : 5.0 #kg
            },
        },
    'chassis' : {
        'mass' : 659 #kg
        },
    'science_payload' : {
        'mass' : 75 #kg
        },
    'power_subsys' : {
        'mass' : 90 #kg
        }
    }

planet = {
    'g' : 3.72 #m/s^2
    }

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

    


