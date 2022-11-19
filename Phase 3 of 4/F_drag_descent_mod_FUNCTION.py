import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator as pchip
from scipy.integrate import solve_ivp
from define_planet import *
from define_edl_system import *


def mod_C(edl_system,velocity, altitude):
    Cd = edl_system['parachute']['Cd']
    Mach = np.array([0.25,0.5,0.65,0.7,0.8,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,1.9,2,2.2,2.5,2.6])
    MEF = np.array([1,1,1,0.98,0.9,0.72,0.66,0.76,0.9,0.96,0.99,0.999,0.992,0.98,0.9,0.85,0.82,0.75,0.65,0.62])
    CD_func = interp1d(Mach, MEF, kind='cubic',fill_value="extrapolate")
    # plt.plot(Mach,MEF,'o', label='Original Data')
    # plt.plot(Mach, CD_func(Mach),label='interpolated function')
    # plt.xlabel('Mach')
    # plt.ylabel('MEF')
    # plt.title('MEF vs. Mach')
    # plt.legend(loc='upper right')
    Mach_number = v2M_Mars(velocity, altitude)
    CDMod = CD_func(Mach_number)
    return CDMod * Cd

def F_drag_descent_mod(edl_system,planet,altitude,velocity):
    
    # Compute the net drag force. 
    
    
    # compute the density of planetary atmosphere at current altitude
    density, _, _ = get_local_atm_properties(planet, altitude)
    
    # This is the (1/2)*density*velocity^2 part of the drag model. The missing
    # bit is area*Cd, which we'll figure out below.
    rhov2=0.5*density*velocity**2
    
    
    # *************************************
    # Determine which part(s) of the EDL system are contributing to drag
    
    # If the heat shield has not been ejected, use that as our drag
    # contributor. Otherwise, use the sky crane.
    if not edl_system['heat_shield']['ejected']:
        ACd_body = np.pi*(edl_system['heat_shield']['diameter']/2.0)**2*edl_system['heat_shield']['Cd']
    else:
        ACd_body = edl_system['sky_crane']['area']*edl_system['sky_crane']['Cd']

    
    # if the parachute is in the deployed state, need to account for its area
    # in the drag calculation
    if edl_system['parachute']['deployed'] and not edl_system['parachute']['ejected']:
        ACd_parachute = np.pi*(edl_system['parachute']['diameter']/2.0)**2*mod_C(edl_system, velocity, altitude)
    else:
        ACd_parachute = 0.0
    
    
    # This computes the ultimate drag force
    F=rhov2*(ACd_body+ACd_parachute)
    
    return F
