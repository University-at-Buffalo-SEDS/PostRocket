# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 14:09:53 2023

@author: Nicholas McNally, University at Buffalo

Edited on Mon June  3 14:43:47 2024
by: Richard Magnani, University at Buffalo

Refs:
    (1): https://web.stanford.edu/~cantwell/AA284A_Course_Material/AA284A_Resources/Zilliac%20and%20Karabeyoglu%20AIAA%202006-4504%20Hybrid_Rocket_Fuel_Regression_Rate_Data_and_Modeling.pdf
    (2): https://web.stanford.edu/~cantwell/AA284A_Course_Material/Karabeyoglu%20AA%20284A%20Lectures/AA284a_Lecture2.pdf
    (3): https://web.stanford.edu/~cantwell/AA284A_Course_Material/Karabeyoglu%20AA%20284A%20Lectures/AA284a_Lecture10.pdf
    (4): https://web.stanford.edu/~cantwell/AA284A_Course_Material/Karabeyoglu%20AA%20284A%20Lectures/AA284a_Lecture7_Efficiency.pdf
    (5): https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
    (6): https://www.researchgate.net/publication/323207484_Building_and_Testing_of_a_Lab-Scale_Hybrid_Rocket_Motor_For_Acquisition_of_Regression_Rate_Equations
    (7): http://mae-nas.eng.usu.edu/MAE_5540_Web/propulsion_systems/section5/Section_5.2.pdf
    (8): https://webbook.nist.gov/cgi/cbook.cgi?ID=C80626&Units=SI&Mask=28F
    (9): https://www.researchgate.net/publication/280923213_Thermal_Conductivity_Enhancement_by_using_Nano-material_in_Phase_Change_Material_for_Latent_Heat_Thermal_Energy_Storage_Systems#pf3
   (10): https://webbook.nist.gov/cgi/cbook.cgi?ID=C630046&Units=SI&Mask=4#Thermo-Phase 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import cantera as ct
from scipy.optimize import fsolve
import scipy as sc
import scipy.integrate as integrate
from cycler import cycler

#%% Inputs
# Fuel
Do = .044 # Grain outer diameter [m]
Dpi = .0128 # Initial port diameter [m]
L = 0.1613 # Length of fuel grain [m]
rhof = 1190.24 # Density of PMMA fuel (lab calculation) [kg/m^3]
Tb = 653 # Pyrolysis temperature of PMMA [K]
m0 = rhof*((Do/2)**2*np.pi*L) # Initial fuel mass [kg]
# Ox tank
POx0 = 6.895e6*.15 # Initial oxidizer tank pressure (1,000 psi) [Pa]
gamOx = 1.4 # Oxidizer specific heat ratio (O2, oxygen)
Vtank = 0.0023012479 # Oxidizer tank volume [m^3]
Dstar = 0.0025 # Oxidizer feed minimum diameter, for choking area for ox injection [m]
Astar = Dstar**2*np.pi/4 # Oxidizer choking area [m^2]
MWox = 32 # Molecular weight fo O2 [kg/kmol]
# Nozzle
Dthroat = .00794#.01191#0.0196 # Throat diameter [m]
Athroat = np.pi/4*Dthroat**2 # Area of nozzle throat [m^2]
Dexit = .009382#.012949#0.0326 # Nozzle exit diameter [m]
Aexit = np.pi/4*Dexit**2 # Nozzle exit area [m^2]
con = 0.98 # Conical nozzle correction factor for 15° half angle
# Constants
Ru = ct.gas_constant # Universal gas constant [J/(kmol*K)]
mdotox_ss = 8.5/1000 # Steady-state oxidizer mass flow rate [kg/s]
n = 0.615 # Burn rate parameter - from (1), pg 16, Table 2
a = .015/1000 #(.0211 from literature) Burn rate parameter - from (1) pg 16, Table 2. Unit conversion from (6), Table 1 (reported in [mm/s], /1000 to get [m/s]).  Gives [m/s] regression rate from [kg/(m^2*s)] oxidizer mass flux
m = 0.26 # (0.26 original) Burn rate parameter - downport dependance (guess)
eta_c = 0.90 # Combustion efficiency, typical range 0.90-0.98 according to (4) slide 12
Cd = 1.15 # Discharge coefficient, ranges from 1.0-1.15 according to (7)
# Reference values
Tref = 298.15 # Reference temperature (STP) [K]
Pref = ct.one_atm # Reference pressure (STP) [Pa]
# Tunables
Tburn = 8 # Time of burn [s]
N = 1000 # Number of time steps
M = 1000 # Number of spatial discritizations
Ptol = 0.01 # Pressure iteration tolerance
#%% Test Data


#%% Initialization
T = np.linspace(0,Tburn,N) # Time vector [s]
X = np.linspace(0.000001,L,M) # Space vector - Small non-zero initial value to avoid /0
Xmat = np.tile(np.transpose(X),(N,1)) # X vector repeated for each time step
R = np.zeros((N,M)) # Radius matrix
mdotox = np.zeros(np.size(T)) # Initialize mdotox vector [kg/s]
MDOT = np.zeros((N,M)) # Total Mass flow rate matrix [kg/s]
gam = np.zeros((N,)) # Specific heat ratio
Me = np.zeros((N,)) # Mach number
Te = np.zeros((N,)) # nozzle exit temeprature [K]
Pe = np.zeros((N,)) # Nozzle exit pressure [Pa]
Ve = np.zeros((N,)) # Nozzle exit velocity [m/s]
Isp = np.zeros((N,)) # Specific impulse
Thrust = np.zeros((N,)) # Thrust force [N]
ThrustPressure = np.zeros((N,)) # Thrust force [N]
ThrustMomentum = np.zeros((N,)) # Thrust force [N]

#%% Oxidizer mass flow rate calculation
# Calculations in this section are accoring to (8) for an adiabatic tank
c0 = (gamOx*Ru*Tref/MWox)**0.5 # Initial speed of sound in oxidizer - from (8) eqn. 19
tau = Vtank/(Astar*c0)*((gamOx+1)/2)**((gamOx+1)/(2*(gamOx-1))) # discharge time constant - from (8) eqn. 18
rhoOx0 = POx0/((Ru/MWox)*Tref) # Initial oxidizer density [kg/m^3]

rhoOx = rhoOx0*(1+(gamOx-1)/2*(T/tau))**(2/(1-gamOx)) # Oxidizer density over time [kg/m^3]
POx = POx0*(1+(gamOx-1)/2*(T/tau))**(2*gamOx/(1-gamOx))
mOx = rhoOx*Vtank # Oxidizer mass over time [kg]
drhoOxdt = rhoOx/tau*(rhoOx/rhoOx0)**((gamOx-1)/2) # Change in oxidizer density over time [kg/m^3/s] - from (8) eqn. 24
mdotox = drhoOxdt*Vtank # Oxidizer mass flow rate [kg/s]
mdotox = np.loadtxt('mdotox.txt') # Overwrite and load mdotox info from text file (from test data) [kg/s]
#%% Port dimaeter calcualtion
def port(t,x,r,mdot,a,m,n,rhof):
    drdt = (a/x**m)*(mdot/(np.pi*r**2))**n # Regression rate equation [m/s]
    dmdotdx = rhof*2*np.pi*r*drdt # Mass flow growth equation [kg/(s*m)]
    return np.array([drdt,dmdotdx])

# Initial/Boundary Conditions
mdot = mdotox # Initial mdot is only oxidizer [kg/s]
MDOT[:,0] = np.zeros([N,])+mdot # Mdot at first spatial coordinate is always mdotox

r = Dpi/2 # Initial port radius [m]
R[0,:] = np.zeros([M])+r # Initial port is r everywhere

# RK-4 iteration to solve for total mdot and port radius as functions of space and time
ht = T[1]-T[0] # Time step [s]
hx = X[1]-X[0] # Delta x [m]
for i in np.arange(0,np.size(T)-1):
    t = T[i] # Current time [s]
    print(t)
    
    for j in np.arange(0,np.size(X)-1):
        x = X[j] # Current spatial corrdinate [m]
        
        r = R[i,j]
        mdot = MDOT[i,j]
        
        k1 = port(t,x,r,mdot,a,m,n,rhof)
        k1r = k1[0]
        k1m = k1[1]
        
        k2 = port(t+ht/2,x+hx/2,r+ht*k1r/2,mdot+hx*k1m/2,a,m,n,rhof)
        k2r = k2[0]
        k2m = k2[1]
        
        k3 = port(t+ht/2,x+hx/2,r+ht*k2r/2,mdot+hx*k2m/2,a,m,n,rhof)
        k3r = k3[0]
        k3m = k3[1]
        
        k4 = port(t+ht,x+hx,r+ht*k1r,mdot+hx*k1m,a,m,n,rhof)
        k4r = k4[0]
        k4m = k4[1]
        
        R[i+1,j] = min(R[i,j]+(1/6)*(k1r+2*k2r+2*k3r+k4r)*ht,Do/2) # Solve for grain radius, while restricting it to be grain OD at maximum [m]
        if R[i+1,j] == Do/2:
            MDOT[i,j+1] = MDOT[i,0] # If grain is at OD, mdot is equal to mdotox (boundary condition) at this time step and location
        else:
            MDOT[i,j+1] = MDOT[i,j]+(1/6)*(k1m+2*k2m+2*k3m+k4m)*hx # Otherwise, set as normal
        
# Truncate to eliminate final entry in space and time
MDOT = MDOT[0:N-1,0:M-1]
R = R[0:N-1,0:M-1]
mdotox = mdotox[0:N-1]

# Mass flow rates
mdotTotal = MDOT[:,M-2] # Total mass fuel rate into nozzle [kg/s]
mdotf = mdotTotal-mdotox # Total fuel mass flow rate into nozzle [kg/s]
OF = mdotox/mdotf # Oxidizer to fuel ratio - in post chamber

# G and rdot
G = MDOT/(np.pi*R**2) # Total mass flux in port for all space and time [kg/(m^2*s)]
Gox = np.tile(mdotox,(M-1,1)).transpose()/(np.pi*R**2) # Total OXIDIZER ONLY flux in port for all space and time [kg/(m^2*s)]
rdot = a/(Xmat[0:N-1,0:M-1]**m)*G**n # Regression rate at each point in time and space [m/s]
#%% Combustion chemistry
exhaust = ct.Solution('MMAReduced.yaml') # MMETHAC_C5H8O2 - PMMA name
Tc = np.zeros((N-1,))
P = np.zeros((N-1,))
MWexhaust = np.zeros((N-1,))

Pguess = Pref # Set first guess to Pref [Pa]
PPrev = Pref # Set last guess to Pref
i = 0 # Counting index
for of in OF:
    # Need to iterate pressure guesses because pressure is needed for equilibrate, but pressure is not yet known
    while abs(P[i]-PPrev)/PPrev>Ptol: # While %error of pressure is over 1%
        PPrev = Pguess # Set last pressure as Pguess [Pa]
        exhaust.TPY = Tref, Pref, {'MMETHAC_C5H8O2':1, 'O2':of} # Set combustion reactants. Pref is used as chamber pressure is not yet known. However, the sensitivity appears to be relatively small, but may try to do an iterative approach later.
        exhaust.equilibrate('HP')
        Tc[i] = exhaust.T # Chamber (stagnation) temperature, calculated as adiabaitic flame temp [K]
        MWexhaust[i] = exhaust.mean_molecular_weight # Molecular weight of exhaust [kg/kmol]
        gam[i] = exhaust.cp_mass/exhaust.cv_mass  #Specific heat ratio of exhaust
        cstar = ((1/gam[i])*((gam[i]+1)/2)**((gam[i]+1)/(gam[i]-1))*Ru*Tc[i]/MWexhaust[i])**(1/2) # c* from (2) slide 16
        P[i] = (mdotTotal[i])*cstar*eta_c/(Athroat*Cd) # Chamber pressure [Pa] from (3) [Pa]

        Pguess = P[i] # Set next guess as newly calculated pressure [i]
        
    i += 1 # Incriment counting index
        
#%% Thrust Calculations
def m_A_relate(Me):
    return ((gam[i]+1)/2)**(-(gam[i]+1)/(2*(gam[i]-1)))*((1+(gam[i]-1)/2*Me**2)**((gam[i]+1)/(2*(gam[i]-1))))/Me-Aexit/Athroat # Mach area relation equal to zero to solve for exit Mach number (5)

for i in np.arange(0, N-1, 1):
    PuPd = (P[i]+Pref)/Pref # Chamber pressure to ambient pressure ratio (for choking/separation condition)
    
    condition = (2/(gam[i]+1))**(-gam[i]/(gam[i]-1)) # Minimum pressure ratio for choking to occur
    if PuPd >= condition:
        M0 = 2 # Supersonic root initial guess (flow is choked at nozzle)
    else:
        M0 = 0.5 # Subsonic root initial guess (flow is not choked)
    
    Me[i] = fsolve(m_A_relate, M0) # Solve for exit Mach number with initial guess 
    Te[i] = Tc[i]/(1+(gam[i]-1)/2*Me[i]**2) # Exit temperature [K]
    if M0 == 2:
        Pe[i] = P[i]*(1+(gam[i]-1)/2*Me[i]**2)**(-gam[i]/(gam[i]-1)) # Exit pressure - If supersonic, find according to isentropic relations [Pa]
    else:
        Pe[i] = Pref # Exit pressure - If subsonic is always equal to ambient pressure [Pa]
        
    Ve[i] = Me[i]*np.sqrt(gam[i]*Ru/MWexhaust[i]*Te[i]) # Exit velocity [m/s]
    Thrust[i] = mdotTotal[i]*Ve[i]*con+(Pe[i]-Pref)*Aexit # Thrust [N]

Isp = Ve/9.81 # Specific impulse [s]

mf = rhof*((Do/2)**2*np.pi*L-sc.integrate.trapezoid((R[N-2,:])**2*np.pi,dx=hx))
Deltamf = m0-mf # Total fuel mass burned for empirical model [kg]
Ispavg = np.mean(Isp) # Average specific impulse [s]

#%% Finding sensible enthalpy of PMMA
gas = ct.Solution('MMAReduced.yaml')

fuel = ct.Quantity(gas,constant='TP') # Create fuel quantity from mixture
fuel.TPX = Tb, Pref, {'MMETHAC_C5H8O2':1} # Set to PMMA only, set to boiling temp
fuel.equilibrate('TP')
hTotalS = fuel.enthalpy_mass # Get enthalpy after equilibrate [J/kg]
fuel.TP = Tref,Pref # Set to STP
hfs = fuel.enthalpy_mass # Get enthalpy after cooling [J/kg]
hs = hTotalS-hfs # Sensible enthalpy of fuel [J/kg]

#%% Oxidizer
ox = ct.Quantity(gas,constant='TP') # Create oxidizer quantity from mixture
ox.TPX = Tref,Pref,{'O2':1} # Set to 1 mole of O2 @ STP
he = 0 # Sensible enthalpy of oxidizer - =0 if assuming it enters @ STP [J/kg]
mue = ox.viscosity # Viscosity of oxidizer [N/(m^2*s)]
rhoe = ox.density # Density of oxidizer [kg/m^3]

# Molecular Weights
MWO2 = gas.molecular_weights[gas.species_index('O2')] # Molecular weight of O2 [kg/kmol]
MWfuel = gas.molecular_weights[gas.species_index('MMETHAC_C5H8O2')] # Molecular weight of PMMA [kg/kmol]

# Calculating nu
gas.X = {'MMETHAC_C5H8O2':1, 'O2':1} # Create gas with O2 and PMMA
gas.set_equivalence_ratio(1, fuel='MMETHAC_C5H8O2:1', oxidizer='O2') # Set to stoichiometric conditions
XO2 = gas.X[gas.species_index('O2')] # Mole fraction of oxidizer in stoich reation
Xfuel = gas.X[gas.species_index('MMETHAC_C5H8O2')] # Mole fraction of PMMA in stoich reaction
nO2 = XO2/Xfuel # Number of moles of oxidizer per 1 mole of fuel
nuox = nO2*MWO2/MWfuel # mass of oxidizer per kg fuel for O2
Yoxe = 1 # Mass fraction is 1 because all of oxidizer is reactive (no nitrogen in air etc.)

#%% Mix fuel and air and react
fuel.mass = 1 # per 1 kg of fuel [kg]
fuel.TP = Tb,Pref # Set fuel to boiling point at atm pressure

ox.mass = nuox # Mass of oxidizer per kg fuel [kg]
ox.TP = Tref,Pref # Set oxidizer to STP

mix = fuel+ox # Mix fuel and oxidizer by adding quantities
gas.Y = mix.Y # Set gas object mass fractions to mixture
hR = gas.enthalpy_mass # Enthalpy of mixed,unburnt mixture [J/kg]

gas.equilibrate('HP') # Burn mixture with constant enthalpy and pressure
Tad = gas.T # Adiabatic flame temperature of stoichiometric burn [K]
hflame = gas.enthalpy_mass # Enthalpy of burnt mixture [J/kg]
gas.TP = Tref,Pref # Iobarically cool to Tref
hP = gas.enthalpy_mass # Cooled burnt mixture enthalpy [J/kg]
hflame = hflame-hP # Sensible enthalpy of flame [J/kg]

#%% Heat of combustion
DHc = hR-hP # Heat of combustion per unit mass of mixture [J/kg]
DHc = (nuox+1)*DHc # Heat of combustion per unit mass of fuel [J/kg]

hfg = 38e6/MWfuel # Heat of vaporization of PMMA from (8) [J/kg]
B = (DHc*Yoxe/nuox+he-hs)/hfg # Blowing parameter
nturb = 1/7 # 1/7 law for u*

CfCfo = (1+0.5*(1-np.exp(-0.05*B)))*np.log(1+B)/B # Coefficient of friction ratio term

# flame location, etaf relative to BL thickness
bf = (hflame-hs)/(DHc*Yoxe/nuox + he-hs)
etaf = np.log(1+bf*B)/np.log(1+B)
etaf = pow(etaf,1/nturb)
#%% Compute I = int u^* (1-u^*) d eta where u^* = ((1+B)^(eta^nturb)-1)/B
# and compute constant term in boundary layer height equation, "fac"
nturb = 1/7 # 1/7 law for u* 
I = integrate.quad(lambda eta: ((pow(1+B,pow(eta,nturb))-1)/B)*(1.-(pow(1+B,pow(eta,nturb))-1)/B), 0, 1.)[0]
fac = pow(5/4*CfCfo*0.0225*(1.+B)/I,4/5) # Boundary layer height 4/8 power term
#%% Calculate rdot and r with Marxman theory
rMarx = np.zeros((N,M)) # Radius of grain from Marxman theory for all time and space [m]
rdotMarx = np.zeros((N,M)) # Regression rate of grain from Marxman theory for all time and space [m/s]
delta = np.zeros((N,M)) # Boundary layer height from Marxman theory for all time and space [m]
yf = np.zeros((N,M)) # Flame height from Marxman theory for all time and space [m]
GMarx = np.zeros((N,M)) # Gox from Marxman theory for all time and space [kgOx/(m^2*s)]
mdotfMarx = np.zeros((N,M)) # Fuel mass regression rate from Marxman theory for all time and space [kg/s]

mdotox = np.append(mdotox,mdotox_ss) # Add back on final point to make size equal to N instead of N-1 (Needed for rdot calculation, but then re-removed after)

# IC and BC
rMarx[0,:] = Dpi/2 # Set initial port radius at t=0 for all x [m]
GMarx[0,:] = mdotox[0]/(np.pi*(Dpi/2)**2) # Set initial G at t=0 for all x [kgOx/(m^2*s)]

# Iteration
i=0 # Initialize time index
for t in T[0:N-1]:
    j=0 # Initialize space index at each new time step
    
    for x in X[0:M-1]:
        Rex = x*GMarx[i,j]/mue # Reynolds number at location
        deltaX = fac*pow(Rex,-1/5) # change in boundary layer height at location
        delta[i,j] = x*deltaX # Boundary layer hight at location [m]
        yf[i,j] = delta[i,j]*etaf # Flame height [m]
        Redelta = delta[i,j]*GMarx[i,j]/mue # Turbulent Reynolds number for Coefficient of friction calc
        Cfo = 2.*0.0225*pow(Redelta,-0.25) # Coefficient of friction at x=0
        Cf = Cfo*CfCfo # Coefficient of friction
        rdotMarx[i,j]=Cf*GMarx[i,j]*B/rhof/2 # Regression rate from Marxman theory [m/s]
        mdotfMarx[i,j] = 2*np.pi*rMarx[i,j]*rhof*rdotMarx[i,j]*hx # Mass of fuel regressed for every spatial coordinate for all time [kg/s]
        j+=1 # Iterate space index
    
    rMarx[i+1,:] = rMarx[i,:] + ht*rdotMarx[i,:] # Find next time step port radius by adding current radius to rdot multipleis by time step size [m]
    GMarx[i+1,:] = mdotox[i]/(np.pi*rMarx[i+1,:]**2) # Next time step Gox [kgOx/(m^2*s)]
    
    
    i+=1 # Iterate time index

# Eliminate last row and column from results
rMarx = rMarx[0:N-1,0:M-1]
GMarx = GMarx[0:N-1,0:M-1]
mdotfMarx = mdotfMarx[0:N-1,0:M-1]
mdotox = mdotox[0:N-1]
#%% Marxman Post-Chamber Combustion Chemistry
OFMarx = mdotox/np.sum(mdotfMarx, axis=1)

exhaust = ct.Solution('MMAReduced.yaml') # MMETHAC_C5H8O2 - PMMA name
TcMarx = np.zeros((N-1,))
PMarx = np.zeros((N-1,))
MWexhaust = np.zeros((N-1,))

i = 0 # Counting index
for of in OFMarx:
    exhaust.TPY = Tref, Pref, {'MMETHAC_C5H8O2':1, 'O2':of}
    exhaust.equilibrate('HP')
    TcMarx[i] = exhaust.T # Chamber (stagnation) temperature, calculated as adiabaitic flame temp [K]
    MWexhaust[i] = exhaust.mean_molecular_weight # Molecular weight of exhaust [kg/kmol]
    gam[i] = exhaust.cp_mass/exhaust.cv_mass  #Specific heat ratio of exhaust
    cstar = ((1/gam[i])*((gam[i]+1)/2)**((gam[i]+1)/(gam[i]-1))*Ru*Tc[i]/MWexhaust[i])**(1/2) # c* from (2) slide 16
    PMarx[i] = (mdotox[i]+np.sum(mdotfMarx[i]))*cstar*eta_c/(Athroat*Cd) # Chamber pressure [Pa] from (3) [Pa]

    i += 1 # Incriment counting index
        
#%% Thrust Calculations
MeMarx = np.zeros((N,)) # Mach number
TeMarx = np.zeros((N,)) # nozzle exit temeprature [K]
PeMarx = np.zeros((N,)) # Nozzle exit pressure [Pa]
VeMarx = np.zeros((N,)) # Nozzle exit velocity [m/s]
IspMarx = np.zeros((N,)) # Specific impulse
ThrustMarx = np.zeros((N,)) # Thrust force [N]

def m_A_relate(Me):
    return ((gam[i]+1)/2)**(-(gam[i]+1)/(2*(gam[i]-1)))*((1+(gam[i]-1)/2*Me**2)**((gam[i]+1)/(2*(gam[i]-1))))/Me-Aexit/Athroat # Mach area relation equal to zero to solve for exit Mach number (5)

for i in np.arange(0, N-1, 1):
    PuPd = (PMarx[i]+Pref)/Pref # Chamber pressure to ambient pressure ratio (for choking/separation condition)

    condition = (2/(gam[i]+1))**(-gam[i]/(gam[i]-1)) # Minimum pressure ratio for choking to occur
    if PuPd >= condition:
        M0 = 2 # Supersonic root (flow is choked at nozzle)
    else:
        M0 = 0.5 # Subsonic root (flow is not choked)
    
    MeMarx[i] = fsolve(m_A_relate, M0) # Solve for exit Mach number with initial guess of Me=2
    TeMarx[i] = TcMarx[i]/(1+(gam[i]-1)/2*Me[i]**2) # Exit temperature [K]
    if M0 == 2:
        PeMarx[i] = PMarx[i]*(1+(gam[i]-1)/2*MeMarx[i]**2)**(-gam[i]/(gam[i]-1)) # Exit pressure - If supersonic, find according to isentropic relations [Pa]
    else:
        PeMarx[i] = Pref # Exit pressure - If subsonic is always equal to ambient pressure [Pa]
        
    VeMarx[i] = MeMarx[i]*np.sqrt(gam[i]*Ru/MWexhaust[i]*TeMarx[i]) #Exit velocity [m/s]
    ThrustPressure[i] = (Pe[i]-Pref)*Aexit # Thrust due to pressure difference [N]
    ThrustMomentum[i] = mdotTotal[i]*Ve[i]*con # Thrust due to momentum transfer [N]
    ThrustMarx[i] = (mdotox[i]+np.sum(mdotfMarx, axis=1)[i])*VeMarx[i]*con+(PeMarx[i]-Pref)*Aexit # Thrust [N]

mfMarx = rhof*((Do/2)**2*np.pi*L-sc.integrate.trapezoid((rMarx[N-2,:])**2*np.pi,dx=hx))
DeltamfMarx = m0-mfMarx # Total fuel mass burned for Marxman model [kg]
IspMarx = VeMarx/9.81 # Specific impulse [s]

print(IspMarx)




#%% Plotting
# # Set cycler for changing plo line colors and styles
# default_cycler = (cycler(color=['k', 'r', 'b', 'k', 'r', 'b']) +
#                   cycler(linestyle=['-', '-', '-', ':', ':', ':']))
# plt.rc('lines', linewidth=5)
# plt.rc('axes', prop_cycle=default_cycler)
# # rdot vs x comparison
# plt.figure(figsize=(12,9))
# plt.xlabel('Grain position [mm]',fontsize=30)
# plt.ylabel(r'$\dot{r}$ [$\frac{mm}{s}$]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.plot(X[5:M-1]*1e3,rdotMarx[0,5:M-1]*1e3, linewidth=5, label=r'$t=0s$')
# plt.plot(X[5:M-1]*1e3,rdotMarx[int(N/2),5:M-1]*1e3, linewidth=5, label=r'$t=$' + str(round(Tburn/2,2)) + r'$s$')
# plt.plot(X[5:M-1]*1e3,rdotMarx[N-2,5:M-1]*1e3, linewidth=5, label=r'$t=$' + str(round(Tburn,2)) + r'$s$')
# plt.plot(X[5:M-1]*1e3,rdot[0,5:M-1]*1e3, linewidth=5)
# plt.plot(X[5:M-1]*1e3,rdot[int(N/2),5:M-1]*1e3, linewidth=5)
# plt.plot(X[5:M-1]*1e3,rdot[N-2,5:M-1]*1e3, linewidth=5)
# plt.plot(0,rdot[0,5]*1e3,label='Marxman Theory',color='k',linestyle='-') # Solid line label (no data)
# plt.plot(0,rdot[0,5]*1e3,label='Empirical Model',color='k',linestyle=':') # Dotted line label (no data)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize=23)
# plt.ylim(0,1.6)
# plt.show()

# # Port Radius
# plt.figure(figsize=(12, 9))
# plt.plot(X[0:M-2]*1000,R[int((N-2)*0.5),0:M-2]*1000,label='Simple Theory Model', linewidth=5, color='k', linestyle='-')
# plt.plot(X[0:M-2]*1000,R[int((N-2)*0.75),0:M-2]*1000, linewidth=5, color='k', linestyle=':')
# plt.plot(X[0:M-2]*1000,R[N-2,0:M-2]*1000, linewidth=5, color='k', linestyle='-.')
# plt.plot(X[0:M-2]*1000,rMarx[int((N-2)*0.5),0:M-2]*1000,label='Marxman Theory Estimate', linewidth=5, color='b', linestyle='-')
# plt.plot(X[0:M-2]*1000,rMarx[int((N-2)*0.75),0:M-2]*1000, linewidth=5, color='b', linestyle=':')
# plt.plot(X[0:M-2]*1000,rMarx[N-2,0:M-2]*1000, linewidth=5, color='b', linestyle='-.')
# plt.plot(0,Dpi/2*1000, color='k', linestyle='-', label=r'$t=$' + str(round(Tburn*0.5,2)) + r'$s$')
# plt.plot(0,Dpi/2*1000, color='k', linestyle=':', label=r'$t=$' + str(round(Tburn*0.75,2)) + r'$s$')
# plt.plot(0,Dpi/2*1000, color='k', linestyle='-.', label=r'$t=$' + str(round(Tburn,2)) + r'$s$')
# plt.xlabel('Grain position [mm]',fontsize=30)
# plt.ylabel('Port Radius [mm]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.ylim(Dpi/2*1000,Do/2*1000) # Fix min radius on graph at initial port radius
# plt.legend(loc=0,fontsize=23)
# plt.show()

# # Set cycler for changing plo line colors and styles
# default_cycler = (cycler(color=['r', 'b', 'g', 'm', 'y', 'c', 'k']) +
#                   cycler(linestyle=['-', '-', '-', '-', '-', '-', '-']))
# plt.rc('lines', linewidth=5)
# plt.rc('axes', prop_cycle=default_cycler)
# # Ox mdot
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],mdotox*1000, linewidth=5,label='Oxidizer Mass Flow Rate Input')
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel('Oxidizer mass flow rate [g/s]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # Fuel mdot
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],mdotf*1000, linewidth=5,label='Empirical Model')
# plt.plot(T[0:N-1],np.sum(mdotfMarx, axis=1)*1000, linewidth=5,label='Marxman Model')
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel('Fuel mass flow rate [g/s]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # G at end of port
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],G[:,M-2], linewidth=5,label=r'Empirical Model ($G_{total}$)')
# plt.plot(T[0:N-1],GMarx[:,M-2]+np.sum(mdotfMarx, axis=1)/(np.pi*rMarx[:,M-2]**2), linewidth=5,label=r'Marxman Model ($G_{total}$)')
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel(r'G [$\frac{kg}{(m^2 \cdot s)}$]',fontsize=30)
# plt.title('Total G and end of port', fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # Boundary Layer & Flame height
# plt.figure(figsize=(L*280, Do*280))
# plt.plot(X[0:M-1]*1000,R[0,:]*1000, linewidth=5,label=r'Initial Port radius',linestyle='-',color='k')
# plt.plot(X[0:M-1]*1000,1000*(R[0,0:M-1]-delta[0,0:M-1]), linewidth=5,label=r'Initial Boundary Layer',linestyle=':',color='b')
# plt.plot(X[0:M-1]*1000,1000*(R[0,0:M-1]-yf[0,0:M-1]), linewidth=5,label=r'Initial Flame',linestyle='-.',color='r')
# plt.plot(X[0:M-1]*1000,-R[0,:]*1000, linewidth=5,linestyle='-',color='k') # Reflection of port radius
# plt.plot(X[0:M-1]*1000,-1000*(R[0,0:M-1]-delta[0,0:M-1]), linewidth=5,linestyle=':',color='b') # Reflection of boundary layer
# plt.plot(X[0:M-1]*1000,-1000*(R[0,0:M-1]-yf[0,0:M-1]), linewidth=5,linestyle='-.',color='r') # Reflection of flame
# plt.xlabel(r'Downport Distance [$mm$]',fontsize=30)
# plt.ylabel(r'Radius [$mm$]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.axis('equal')
# plt.ylim(-Do/2*1000,Do/2*1000) # Fix to grain diameter
# plt.grid(which='both')
# plt.legend(loc='upper right', fontsize = 23)
# plt.show()

# plt.figure(figsize=(L*280, Do*280))
# plt.plot(X[0:M-1]*1000,R[N-2,:]*1000, linewidth=5,label=r'Final Port radius',linestyle='-',color='k')
# plt.plot(X[0:M-1]*1000,1000*(R[N-2,0:M-1]-delta[N-2,0:M-1]), linewidth=5,label=r'Final Boundary Layer',linestyle=':',color='b')
# plt.plot(X[0:M-1]*1000,1000*(R[N-2,0:M-1]-yf[N-2,0:M-1]), linewidth=5,label=r'Final Flame',linestyle='-.',color='r')
# plt.plot(X[0:M-1]*1000,-R[N-2,:]*1000, linewidth=5,linestyle='-',color='k') # Reflection of port radius
# plt.plot(X[0:M-1]*1000,-1000*(R[N-2,0:M-1]-delta[N-2,0:M-1]), linewidth=5,linestyle=':',color='b') # Reflection of boundary layer
# plt.plot(X[0:M-1]*1000,-1000*(R[N-2,0:M-1]-yf[N-2,0:M-1]), linewidth=5,linestyle='-.',color='r') # Reflection of flame
# plt.xlabel(r'Downport Distance [$mm$]',fontsize=30)
# plt.ylabel(r'Radius [$mm$]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.axis('equal')
# plt.ylim(-Do/2*1000,Do/2*1000) # Fix to grain diameter
# plt.grid(which='both')
# plt.legend(loc='upper right', fontsize = 23)
# plt.show()
# # OF ratio
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],OF, linewidth=5,label='Empirical Model')
# plt.plot(T[0:N-1],OFMarx, linewidth=5,label='Marxman Model')
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel('O/F Ratio',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # Chamber Pressure
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],(P)/6895, linewidth=5,label='Empirical Model',linestyle='--',color='k') # Convert to gauge pressure and psi
# plt.plot(T[0:N-1],(PMarx)/6895, linewidth=5,label='Marxman Model',linestyle=':',color='k') # Convert to gauge pressure and psi
# plt.scatter(ccPress_3_12_24[0,:],ccPress_3_12_24[1,:], linewidth=5,label='Test Fire 3/12/24',linewidths=5)
# plt.scatter(ccPress_3_15_24[0,:],ccPress_3_15_24[1,:], linewidth=5,label='Test Fire 3/15/24',linewidths=5)
# plt.scatter(ccPress_4_5_24[0,:],ccPress_4_5_24[1,:], linewidth=5,label='Test Fire 4/5/24',linewidths=5)
# plt.scatter(ccPress_4_18_24[0,:],ccPress_4_18_24[1,:], linewidth=5,label='Test Fire 4/18/24',linewidths=5)
# plt.scatter(ccPress_4_22_24[0,:],ccPress_4_22_24[1,:], linewidth=5,label='Test Fire 4/22/24',linewidths=5)
# plt.scatter(ccPress_4_23_24[0,:],ccPress_4_23_24[1,:], linewidth=5,label='Test Fire 4/23/24',linewidths=5)
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel('Chamber Pressure [psi]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # Chamber Temperature
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],Tc, linewidth=5,label='Empirical Model')
# plt.plot(T[0:N-1],TcMarx, linewidth=5,label='Marxman Model')
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel('Chamber temperature [K]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # Isp
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],Isp[0:N-1], linewidth=5,label='Empirical Model') # *10 to convert to kg/(m^2*s)
# plt.plot(T[0:N-1],IspMarx[0:N-1], linewidth=5,label='Marxman Model') # *10 to convert to kg/(m^2*s)
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel(r'Specific Impulse [$s$]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
# # Thrust
# plt.figure(figsize=(12, 9))
# plt.plot(T[0:N-1],Thrust[0:N-1], linewidth=5,label='Empirical Model',linestyle='--',color='k')
# plt.plot(T[0:N-1],ThrustMarx[0:N-1], linewidth=5,label='Marxman Model',linestyle=':',color='k')
# plt.scatter(thrust_3_12_24[0,:],thrust_3_12_24[1,:]*4.448, linewidth=5,label='Test Fire 3/12/24',linewidths=5)
# plt.scatter(thrust_3_15_24[0,:],thrust_3_15_24[1,:]*4.448, linewidth=5,label='Test Fire 3/15/24',linewidths=5)
# plt.scatter(thrust_4_5_24[0,:],thrust_4_5_24[1,:]*4.448, linewidth=5,label='Test Fire 4/05/24',linewidths=5)
# plt.scatter(thrust_4_18_24[0,:],thrust_4_18_24[1,:]*4.448, linewidth=5,label='Test Fire 4/18/24',linewidths=5)
# plt.scatter(thrust_4_22_24[0,:],thrust_4_22_24[1,:]*4.448, linewidth=5,label='Test Fire 4/22/24',linewidths=5)
# plt.scatter(thrust_4_23_24[0,:],thrust_4_23_24[1,:]*4.448, linewidth=5,label='Test Fire 4/23/24',linewidths=5)
# plt.xlabel('Time [s]',fontsize=30)
# plt.ylabel(r'Thrust [$N$]',fontsize=30)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(which='both')
# plt.legend(loc=0, fontsize = 23)
# plt.show()
