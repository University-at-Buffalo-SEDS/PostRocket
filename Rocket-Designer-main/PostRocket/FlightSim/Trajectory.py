# Author: Alex Post
# Date: 7/15/2024
# Utilizes snippets of Luis's code for 3dof trajectory in initial version.
# Will update to 6dof and edit completely in later versions.

# Imports
import numpy as np
from scipy.integrate import solve_ivp
from .AtmosphereModel import atmosphere_conditions
from ..RocketComponents.RocketComponents import *
from .Drag import *

# Set up 3dof equations to be iterated over in solve_ivp method
def rocket_3dof(t, y, Rocket:rocket, rail_length):
    # Equations of motion:
    # Inputs: 
    # t = Time is the specific time of a given timestep
    # y - A vector containing x velocity[0], z velocity[1], pitch rate[2], x position[3] z position[4] and pitch angle[5] relative to starting
    # m_dot - The massflowrate 
    # Weq
    # burn_time - The burntime of the motor
    # The purpose of this function is to plug into an ode solver.
    Rs = 0.002 #Roughness height
    vel_x = y[0]
    vel_z = y[1]
    pitch_rate = y[2]
    distance = y[3]
    altitude = y[4]
    pitch_angle = y[5]

    velocity = np.sqrt(vel_x**2 + vel_z**2)
    flight_angle = np.arctan2(vel_x,vel_z)

    # Import atmosphere for a given height, velocity, and rocket length:
    atmo = atmosphere_conditions(altitude, velocity, Rocket.length)
    gravity = atmo[0][0]
    temperature = atmo[1][0]
    pressure = atmo[2][0]
    density = atmo[3][0]
    mach = atmo[4][0]
    # Aerodynamic Calculations
    AoA = flight_angle - pitch_angle
    q = 0.5*density*velocity**2
    Lift, Drag, cp = Rocket.get_aero(t,q,AoA,mach,velocity,density,temperature,Rs)
    
    
    # Thrust Calculations
    Thrust =  Rocket.get_thrust(t, pressure)
    mass = Rocket.get_mass(t)
    Iyy = Rocket.get_Iyy(t)

    # Define the equations of motion:
    
    
    if altitude < rail_length*np.cos(pitch_angle):
        
        acceleration = max((Thrust/mass) - gravity*np.cos(pitch_angle), 0)
        du_dt = acceleration * np.sin(pitch_angle)
        dw_dt = acceleration * np.cos(pitch_angle)
        dq_dt = 0 # change in pitch 
        
        dx_dt = vel_x # Change in x 
        dz_dt = vel_z # Change in z, cannot descend while on rail!
        do_dt = 0 # change in pitch
    else:
        
        du_dt = ((Thrust - Drag)*np.sin(pitch_angle) - Lift*np.cos(flight_angle)) / mass # Change in x velocity u
        dw_dt = ((Thrust - Drag)*np.cos(pitch_angle) + Lift*np.sin(flight_angle)) / mass - gravity# Change in z velocity w
        dq_dt = Lift*( cp -  Rocket.get_cg(t) ) / Iyy # TODO: Implement Rocket.get_Iyy!!!!

        dx_dt = vel_x
        dz_dt = vel_z
        do_dt = pitch_rate


    dydt = [du_dt, dw_dt, dq_dt, dx_dt, dz_dt, do_dt]
    
    return dydt


###############################################################################################################################
###############################################################################################################################
def Trajectory(Rocket:rocket, rail_length, launch_angle, input_method:str, rtol, atol, vectorized:bool):
    '''Rocket class, lrail langth [m], Rail angle from vertical [deg], Integration method, Integration rel. tolerance, Integration abs. tolerance, Integration vectorized? [bool]'''
###############################################################################################################################


    # Solve the IVP
    # Initial values
    launch_angle = launch_angle*np.pi/180 # launch rail angle to degrees
    sla = np.sin(launch_angle)
    cla = np.cos(launch_angle)
    y_init = [sla*1e-6, cla*1e-6, 0.0, 0.0, 0.0, launch_angle] #Initial conditions

    # Define the ivp
    t_span = [0, 100] # Arbitrary t_span longer than any anticipated flight
    first_step = 0.1 # seconds
    
    # Define flight events: apogee, rail velocity
    
    def rail_departure(t,y):
        return (y[4] > rail_length)
    
    def apogee(t,y):
        return y[1] > -0.1
    apogee.terminal = True
    cd = []
    rocketflight = solve_ivp(lambda t, y:rocket_3dof(t, y, Rocket, rail_length), t_span, y_init,str(input_method),vectorized=vectorized, rtol = rtol, atol = atol, first_step = first_step, events=(rail_departure, apogee))

    
    rail_time = rocketflight.t_events[0][-1]
    rail_data = rocketflight.y_events[0][-1]
    apogee_time = rocketflight.t_events[1][-1]
    apogee_data = rocketflight.y_events[1][-1]

    trajectory_output = {'time'         : rocketflight.t,
                         'vel_x'        : rocketflight.y[0],
                         'vel_z'        : rocketflight.y[1],
                         'pitch_rate'   : rocketflight.y[2],
                         'distance'     : rocketflight.y[3],
                         'altitude'     : rocketflight.y[4],
                         'pitch'        : rocketflight.y[5],
                         'off_rail_time': rail_time,
                         'off_rail_vel' : np.linalg.norm(rail_data[[0,1]]),
                         'apogee_time'  : apogee_time,
                         'apogee_state' : apogee_data,
                         'apogee'       : apogee_data[4],
                         'AoA'          : np.arctan2(rocketflight.y[0],rocketflight.y[1]) - rocketflight.y[5]}

    return trajectory_output




##################################################################################################################################################################
##################################################################################################################################################################







