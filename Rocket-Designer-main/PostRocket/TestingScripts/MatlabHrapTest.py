import matlab.engine
import os


eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath(os.getcwd() + '\\HRAP - Matlab')) # Adds the entire HRAP - Matlab folder to the matlab path so it can be used
s,x,o,t = eng.HRAP_Function('TestConfigg',150.0,5.0,0.002,True,  0.098,0.014,0.762,0.8,101325.0,  0.013,0.98,293.15,  0.03,0.056,0.95,0.95,  1.0,0.01,0.36, nargout=4)

print()
# function [s,x,o,t] = HRAP_Function( MotorConfiguration,RunTime,BurnTime,Timestep, Display,...
#                                     GrainOD,GrainID,GrainLength,CEfficiency,AmbientPressure, ...
#                                     TankVolume,OxidizerFill,TankTemp, ...
#                                     ThroatDiameter,NozzleExit,NozzleCd,NozzleEfficiency, ...
#                                     NumberofInjectors,InjectorDiameter,InjectorCd)

# TEST CONFIG:
# [s,x,o,t] = HRAP_Function('TestConfigg',150.0,15.0,0.002,false,  0.098,0.014,0.762,0.8,101325.0,  0.013,0.98,293.15,  0.03,0.056,0.95,0.95,  1.0,0.01,0.36);
