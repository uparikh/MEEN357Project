"""###########################################################################
#   Defines mission events.
#
#   Created by: Former Marvin Numerical Methods Engineering Team
#   Last Modified: 22 October 2021
###########################################################################"""

def define_mission_events():
        
    mission_events = {'alt_heatshield_eject' : 8000,    # [m] altitude at which heat shield is ejected
                      'alt_parachute_eject' : 900,      # [m] altitude at which parachute is ejected
                      'alt_rockets_on' : 1800,          # [m] altitude at which power descent sequence is initiated using variable-thrust solid rockets 
                      'alt_skycrane_on' : 7.6}          # [m] altitude appropriate for sky crane operation
    return mission_events