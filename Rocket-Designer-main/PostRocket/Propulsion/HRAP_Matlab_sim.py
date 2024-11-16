import matlab.engine
import os
import numpy as np

def TNAM(grain_diameter, port_diameter , grain_length, oxidizer_volume):
    # Version of TNAM using the matlab implementation of HRAP via matlab.engine.
    global matlab_engine
    if 'matlab_engine' not in globals():
       print('Starting Matlab engine...')
       matlab_engine = matlab.engine.start_matlab()
       print('Adding HRAP subfolders to the Matlab Path...')
       matlab_engine.addpath(matlab_engine.genpath(os.getcwd() + '\\Rocket-Designer-main\\HRAP - Matlab')) # Adds the entire HRAP - Matlab folder to the matlab path so it can be used
       print('Matlab engine ready')

    # TODO: Validate or adjust most values
    MotorConfiguration = 'Test' # motor Name
    RunTime = 150.0 # maximum real run time of the program [s]
    BurnTime = 100.0 # maximum burn time of the motor within the sim [s]
    Timestep = 0.002 # Constant sim time step; 0.002 is the absolute maximum for good results [s]
    Display = False # Whether to display results at the end of the sim [bool]
    GrainOD = grain_diameter # Grain outer diameter [m]
    GrainID = port_diameter # Grain inner diameter [m]
    GrainLength = grain_length # Grain length [m]
    CEfficiency = 0.8 # Combustion efficiency coefficient [0-1]
    AmbientPressure = 101325.0 # Ambient air pressure; 1 atm used as default [Pa]
    TankVolume = oxidizer_volume # Oxidizer tank volume [m^3]
    OxidizerFill = 0.98 # Fraction of the oxidizer tank filled at beginning of burn [0-1]
    TankTemp = 293.15 # Tank starting temperature [K]
    ThroatDiameter = 0.03 # Nozzle throat diameter [m]
    NozzleExit = 0.056 # Nozzle exit diameter [m]
    NozzleCd = 0.95 # Nozzle discharge coefficient [0-1]
    NozzleEfficiency = 0.95 # Nozzle efficiency [0-1]
    NumberofInjectors = 1.0 # Number of oxidizer injectors in the combustion chamber
    InjectorDiameter = 0.01 # Injector diameter [m]
    InjectorCd = 0.36 # Injector discharge coefficient [0-1]

    s,x,o,t = matlab_engine.HRAP_Function(  MotorConfiguration,RunTime,BurnTime,Timestep,Display,
                                            GrainOD,GrainID,GrainLength,CEfficiency,AmbientPressure,
                                            TankVolume,OxidizerFill,TankTemp,
                                            ThroatDiameter,NozzleExit,NozzleCd,NozzleEfficiency,
                                            NumberofInjectors,InjectorDiameter,InjectorCd,nargout=4)


    burn_time = np.array(o['t'])[0]
    oxidizer_pressure_curve = np.array(o['P_tnk'])[0]    
    chamber_pressure_curve = np.array(o['P_cmbr'])[0]
    chamber_pressure_curve[-1] = 101325.0 # Final pressure is clamped to atmospheric
    oxidizer_mass_curve = np.array(o['m_o'])[0]
    fuel_mass_curve = np.array(o['m_f'])[0]
    thrust_curve = np.array(o['F_thr'])[0]
    thrust_curve[-1] = 0.0 # Final thrust is clamped to zero

    return burn_time, oxidizer_pressure_curve, chamber_pressure_curve, oxidizer_mass_curve, fuel_mass_curve, thrust_curve