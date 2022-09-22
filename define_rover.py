#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:52:42 2022

@author: uditparikh
"""

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
    return rover, planet