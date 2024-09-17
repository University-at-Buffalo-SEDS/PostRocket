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
   (11): https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Type=JANAFG&Plot=on
   (12): https://www.engineeringtoolbox.com/specific-heat-ratio-d_608.html 
"""

import numpy as np
import cantera as ct
from scipy.optimize import fsolve
import scipy as sc
import scipy.integrate as integrate

def TNAM(Do,Dpi,L,Ox_Vol):
    #%% Inputs
    ## Fuel #############################################################################################################################
    # Do = .098 # Grain outer diameter [m]
    # Dpi = .01905 # Initial port diameter [m]
    # L = 0.6096 # Length of fuel grain [m]
    # OF = 8       # Oxidizer to Fuel ratio
    # L = 
    rhof = 900 # Density of Paraffin Wax (wikipedia) [kg/m^3]
    rhoNO = 1230.46 # Density of Nitrous Oxide (l) (https://encyclopedia.airliquide.com/nitrous-oxide) [kg/m^3]
    Tb = 730 # Pyrolysis temperature of Paraffin Wax [K]
    m0 = rhof*(((Do/2)**2*np.pi*L)-((Dpi/2)**2*np.pi*L)) # Initial fuel mass [kg]
    ## Ox tank ##########################################################################################################################
    POx0 = 750*6894.75729 # Initial oxidizer tank pressure (750 psi) [Pa]
    gamOx = 1.31 # Oxidizer specific heat ratio (12)
    Vtank = Ox_Vol
    # Vtank = (OF*m0)/(rhoNO) # Oxidizer tank volume [m^3] mass divided by density ...  density of liquid nitrous oxide (assuming incompressible fluid)
    Dstar = 0.0025 # Oxidizer feed minimum diameter, for choking area for ox injection [m]
    Astar = Dstar**2*np.pi/4 # Oxidizer choking area [m^2]
    MWox = 44.0128 # Molecular weight fo N2O [kg/kmol] (11)
    ## Nozzle ###########################################################################################################################
    Dthroat = .00794#.01191#0.0196 # Throat diameter [m]
    Athroat = np.pi/4*Dthroat**2 # Area of nozzle throat [m^2]
    Dexit = .009382#.012949#0.0326 # Nozzle exit diameter [m]
    Aexit = np.pi/4*Dexit**2 # Nozzle exit area [m^2]
    con = 0.98 # Conical nozzle correction factor for 15° half angle
    ## Constants ########################################################################################################################
    Ru = ct.gas_constant # Universal gas constant [J/(kmol*K)]
    mdotox_ss = 8.5/1000 # Steady-state oxidizer mass flow rate [kg/s]
    n = 0.615 # 0.748 Burn rate parameter - from (1), pg 16, Table 2 assuming Paraffin FR4550
    a = 0.427/1000 # 0.427 (.0211 from literature) Burn rate parameter - from (1) pg 16, Table 2, FR4550. Unit conversion from (6), Table 1 (reported in [mm/s], /1000 to get [m/s]).  Gives [m/s] regression rate from [kg/(m^2*s)] oxidizer mass flux 
    m = 0.26 # (0.26 original) Burn rate parameter - downport dependance (guess)
    eta_c = 0.90 # Combustion efficiency, typical range 0.90-0.98 according to (4) slide 12
    Cd = 1.15 # Discharge coefficient, ranges from 1.0-1.15 according to (7)
    ## Reference values #################################################################################################################
    Tref = 298.15 # Reference temperature (STP) [K]
    Pref = ct.one_atm # Reference pressure (STP) [Pa]
    ## Tunables #########################################################################################################################
    Tburn = 6 # Time of burn [s]
    N = 10 # Number of time steps
    M = 1000 # Number of spatial discritizations
    Ptol = 0.01 # Pressure iteration tolerance


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
    # mdotox = np.loadtxt('mdotox.txt') # Overwrite and load mdotox info from text file (from test data) [kg/s]
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
        # print(t)
        
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
    OF = np.zeros(len(MDOT))
    i = 0
    while mdotf[i] > 0: 
        OF[i] = mdotox[i]/mdotf[i] # Oxidizer to fuel ratio - in post chamber
        i += 1
    print('Burn ended early! Burn time: '+ str(i*ht) + ' s')
    burnTime = i*ht

    # G and rdot
    G = MDOT/(np.pi*R**2) # Total mass flux in port for all space and time [kg/(m^2*s)]
    Gox = np.tile(mdotox,(M-1,1)).transpose()/(np.pi*R**2) # Total OXIDIZER ONLY flux in port for all space and time [kg/(m^2*s)]
    rdot = a/(Xmat[0:N-1,0:M-1]**m)*G**n # Regression rate at each point in time and space [m/s]
    #%% Combustion chemistry
    exhaust = ct.Solution("WAXmech_fixed.yaml") # C32 is the Paraffin Wax
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
            exhaust.TPY = Tref, Pref, {'C32':1, 'N2O':of} # Set combustion reactants. Pref is used as chamber pressure is not yet known. However, the sensitivity appears to be relatively small, but may try to do an iterative approach later.
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

    ## Need fuel mass as a vector over time that is solved from mdot fuel ##############################################################
    mfuel = np.zeros((N,)) # mass of fuel as a function of time to be filled in
    mi = m0  # initialize mass variable that decreases over the loops of the iteration
    tstep = N-2 
    for i in range(tstep):
        df = (T[i+1]-T[i]) * 0.5*(mdotf[i]+mdotf[i+1])   # trapezoidal integration formula found on Wikipedia
        mfuel[i] = mi - df
        mi = mfuel[i]
    #end

    #%% Finding sensible enthalpy of Paraffin Wax
    gas = ct.Solution('WAXmech_fixed.yaml')

    fuel = ct.Quantity(gas,constant='TP') # Create fuel quantity from mixture
    fuel.TPX = Tb, Pref, {'C32':1} # Set to Paraffin Wax only, set to boiling temp
    fuel.equilibrate('TP')
    hTotalS = fuel.enthalpy_mass # Get enthalpy after equilibrate [J/kg]
    fuel.TP = Tref,Pref # Set to STP
    hfs = fuel.enthalpy_mass # Get enthalpy after cooling [J/kg]
    hs = hTotalS-hfs # Sensible enthalpy of fuel [J/kg]

    #%% Oxidizer
    ox = ct.Quantity(gas,constant='TP') # Create oxidizer quantity from mixture
    ox.TPX = Tref,Pref,{'N2O':1} # Set to 1 mole of O2 @ STP
    he = 0 # Sensible enthalpy of oxidizer - =0 if assuming it enters @ STP [J/kg]
    mue = ox.viscosity # Viscosity of oxidizer [N/(m^2*s)]
    rhoe = ox.density # Density of oxidizer [kg/m^3]

    # Molecular Weights
    MWO2 = gas.molecular_weights[gas.species_index('N2O')] # Molecular weight of O2 [kg/kmol]
    MWfuel = gas.molecular_weights[gas.species_index('C32')] # Molecular weight of Paraffin Wax [kg/kmol]

    # Calculating nu
    gas.X = {'C32':1, 'N2O':1} # Create gas with O2 and Paraffin Wax
    gas.set_equivalence_ratio(1, fuel='C32:1', oxidizer='N2O') # Set to stoichiometric conditions
    XO2 = gas.X[gas.species_index('N2O')] # Mole fraction of oxidizer in stoich reation
    Xfuel = gas.X[gas.species_index('C32')] # Mole fraction of Paraffin Wax in stoich reaction
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

    hfg = 157.3e6/MWfuel # Heat of vaporization of Paraffin Wax (Kolos gave me this)
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

    mdotox = np.append(mdotox,mdotox_ss) # Add back on final point to make size equal to N instead of N-1 (Needed for rdot calculation, but then re-removed after)
    
    return burnTime,P,mOx,mfuel,Thrust
