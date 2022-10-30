#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:52:42 2022
"""
from numpy import array

def rover1():
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
                'mass' : 5.0, #kg
                'effcy_tau' : array([0, 10, 20, 40, 75, 165]), #Nm
                'effcy' : array([0, 0.60, 0.75, 0.73, 0.55, 0.05]) #percent
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
    return rover, planet