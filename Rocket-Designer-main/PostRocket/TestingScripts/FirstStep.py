# -*- coding: utf-8 -*-
"""
FirstStep.py
@author: Alex Post
"""

import math as mth

###############################################################################
###############################################################################
# Input the weight of the fully loaded system, assuming no atmosphere (no drag),
# Calculate the required thrust to reach a thrust to weight ratio to get a desired 
# off the rail speed, which will then give a required minimum thrust. Then solve for
# an impulse with this average/minimum thrust that will reach an input apogee, and then solve
# for the burntime of that impulse. This will give a lower bound. 

#Assumptions for this initial step:
    # Flat Earth
    # No drag/atmospheric consideration
    # Launching straight up
    # Uniform burn (Thrust curve will be minimum thrust to hit desired off-the-rail speed over the whole burn)
    # Assume CG and weight are constant throughout flight
    
# Set up constants
v_off_rail = 30.48 # m/s, translates to 100 ft/s
rail_length = 9.144 # m, translates to 30 ft rail
g = 9.81 # m/s^2, close to Earth's surface approximation for gravity
apogee_target = 3048 # m, equivalent to 10,000 ft
launch_angle = 0.10472 # deg


def InitialGuess(v_off_rail, rail_length, Weight, apogee_target):
    # This function will give an initial guess on burntime and required thrust for a design.
    # Inputs:  v_off_rail -> m/s
    #          rail_length -> m
    #          Weight -> N
    #          apogee_target -> m
    # Outputs: Thrust_req -> N
    #          tb -> Burntime, seconds
    
    #Approximate the required thrust to weight ratio to reach a specific rail velocity:
    TWR = ((v_off_rail ** 2) / (2 * rail_length * g)) + 1
    # From basic energy kinematics, we can determine an impulse required to hit a height
    Mass = Weight / g #Convert N to kg
    Impulse = Mass * mth.sqrt(2 * g * apogee_target) # Ns
    # Calculate the thrust required to reach the specific TWR for this design:
    Thrust_req = TWR * Weight #N
    # Generate a starting burntime for this case:
    tb = Impulse / Thrust_req # seconds
    
    return Thrust_req, tb
Weight = 271.1538703
Thrust_req, tb = InitialGuess(v_off_rail, rail_length, Weight, apogee_target)
print(Thrust_req)
print(tb)
    



    







