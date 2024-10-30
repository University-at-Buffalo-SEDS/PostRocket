# Author: SEDS Dynamics
# Date: 10/29/2024
# Uses 3dof code as base to make 6dof

# Imports
import numpy as np
from scipy.integrate import solve_ivp
from .AtmosphereModel import atmosphere_conditions
from ..RocketComponents.RocketComponents import *
from .Drag import *

# Set up 6dof equations to be iterated over in solve_ivp method
def rocket_6dof(t, y, Rocket:rocket, rail_length):
    # Equations of motion:
    # Inputs: 
    # t = Time is the specific time of a given timestep
    # y - A vector containing x velocity[0], z velocity[1], pitch rate[2], x position[3] z position[4] and pitch angle[5] relative to starting
    # m_dot - The massflowrate 
    # Weq
    # burn_time - The burntime of the motor
    # The purpose of this function is to plug into an ode solver.

    Rs = 0.002 #Roughness height
    
    # Velocities in r-axis (body frame)
    vel_1 = y[0]
    vel_2 = y[1]
    vel_3 = y[2]
    velocity = np.sqrt(vel_1**2 + vel_2**2 + vel_3**2)
    # Angular velocities in r-axis (body frame)
    omega_1 = y[3]
    omega_2 = y[4]
    omega_3 = y[5]

    # Import atmosphere for a given height, velocity, and rocket length:
    atmo = atmosphere_conditions(altitude, velocity, Rocket.length)
    gravity = atmo[0][0]
    temperature = atmo[1][0]
    pressure = atmo[2][0]
    density = atmo[3][0]
    mach = atmo[4][0]
    # Aerodynamic Calculations
    AoA = np.arccos(vel_1/velocity)
    q = 0.5*density*velocity**2
    Lift, Drag, cp = Rocket.get_aero(t,q,AoA,mach,velocity,density,temperature,Rs)
    
    
    # Thrust Calculations
    Thrust =  Rocket.get_thrust(t, pressure)
    mass = Rocket.get_mass(t)
    Ixx = Rocket.get_Ixx(t)
    Iyy = Rocket.get_Iyy(t)
    Izz = Iyy

    # Define the equations of motion:
    
    
    if altitude < rail_length*np.cos(pitch_angle):
        
        dv1_dt = max((Thrust/mass) - gravity*np.cos(pitch_angle), 0)
        dv2_dt = 0
        dv3_dt = 0
        
        d1_dt = vel_1 # Change in x 
        d2_dt = 0 # Change in z, cannot descend while on rail!
        d3_dt = 0 # change in pitch
    else:
        
        dv1_dt = (Thrust - Drag)/mass - gravity[0]    # Change in x velocity u
        dv2_dt = ((Thrust - Drag)*np.cos(pitch_angle) + Lift*np.sin(flight_angle)) / mass - gravity# Change in z velocity w
        #dv3_dt = 

        d1_dt = vel_1
        d2_dt = vel_2
        d3_dt = vel_3


    dydt = [du_dt, dw_dt, dq_dt, dx_dt, dz_dt, do_dt]
    
    return dydt


###############################################################################################################################
###############################################################################################################################


def Trajectory(Rocket:rocket, rail_length, launch_angle, input_method:str, rtol, atol, vectorized:bool):
    '''Rocket class, lrail langth [m], Rail angle from vertical [deg], Integration method, Integration rel. tolerance, Integration abs. tolerance, Integration vectorized? [bool]'''

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
