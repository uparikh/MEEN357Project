#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib.pyplot import plot, xlabel, ylabel, subplot, subplots_adjust
from define_rover import *
from subfunctions import tau_dcmotor
from numpy import linspace, zeros, multiply

rover, planet = rover1()

subplot(3,1,1)
omega = linspace(0,3.8,20) #increase last input to make more smooth
tau = zeros(len(omega))
for w in range(len(omega)):   
    tau[w] = tau_dcmotor(omega[w], rover['wheel_assembly']['motor'])
plot(tau,omega)
xlabel('Motor Shaft Torque [Nm]')
ylabel('Motor Shaft Speed [rad/s]')

subplot(3,1,2)
power = multiply(omega, tau)
plot(tau,power)
xlabel('Motor Shaft Torque [Nm]')
ylabel('Motor Power [W]')

subplot(3,1,3)
plot(omega,power)
xlabel('Motor Shaft Spped [rad/s]')
ylabel('Motor Power [W]')

subplots_adjust(bottom=-1.0)

