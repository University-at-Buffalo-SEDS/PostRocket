# # Author: Alex Post
# # Date: 8/11/2024

# # Reference: https://openrocket.sourceforge.net/techdoc.pdf
# https://cambridgerocket.sourceforge.net/AerodynamicCoefficients.pdf 
# https://www.grc.nasa.gov/www/k-12/airplane/viscosity.html Air Viscosity
# #This function will calculate the forces occuring on a given rocket:

import numpy as np
import math as mth
# from .Trajectory import Trajectory





def fin_drag(fin_thickness, fin_midchord, fin_number, planform_area, exposed_area, rocket_diameter, Rs, velocity, density, temperature, rocket_length):
    # Convert to imperial because NASA 
    density_imp = density * 0.062428 # lb/ft^3
    temperature_imp = temperature * 1.8 # Rankine
    # Dynamic viscosity of air
    mu0 = 3.62 * 10**-7 #lb-sec/ft^2
    T0 = 518.7 # Rankine
    dyn_visc = mu0 * ((temperature_imp / T0)**1.5) * (((T0 + 198.72)) / (temperature_imp + 198.72))
    dyn_visc = dyn_visc * 47.8802687 # convert to kg / (m*s)
    Rcrit = 51*(Rs/rocket_length) ** -1.039 # Barrowman eq 4-7, for surface roughness sizes submerged in laminar sublayer 5x10^5 = Rcrit
    Re_fin = (velocity * density * fin_midchord) / dyn_visc # dimensionless
    if (Re_fin <= Rcrit) & (Re_fin > 0):
        cffn = 1.328 / Re_fin**(1/2) # laminar flow, Barrowman 4-2
    elif (Re_fin > 5*(10**5)) & (Re_fin < Rcrit):
        # B = Rcrit * ((0.074 / Re_fin**(1/5)) - 1.328 / Re_fin**(1/2))
        # cffn = (0.074 / (Re_fin**(1/5))) - (B / Re_fin)
        cffn = ((3.46 * mth.log10(Re_fin) - 5.6)**(-2)) - (1700 / Re_fin) # First term is the Schoenherr theoretical boundary layer modification of von Karman's work, and the second term is the 
        # Prandtl correction term. Barrowman equation 4-6.
    elif Re_fin > Rcrit:
        cffn = 0.032 * (Rs / fin_midchord) ** (0.2) # Turbulent curve match, takes into account
    else:
        cffn = 1e-6
    
    cd_fin = 2 * cffn * (1 + 2 * (fin_thickness / fin_midchord)) * ((4 * fin_number * planform_area) / (mth.pi * rocket_diameter**2)) # fin drag
    cd_i = 2 * cffn * (1 + 2 * (fin_thickness / fin_midchord)) * ((4 * fin_number * (planform_area - exposed_area)) / (mth.pi * rocket_diameter**2)) #interference drag
    return cd_fin, cd_i




def body_drag(q, rocket_length, rocket_diameter, Sref, mach, t, tburn, velocity, density, temperature, Rs):
    # Convert to imperial because NASA 
    density_imp = density * 0.062428 # lb/ft^3
    temperature_imp = temperature * 1.8 # Rankine
    # Dynamic viscosity of air
    mu0 = 3.62 * 10**-7 #lb-sec/ft^2
    T0 = 518.7 # Rankine
    dyn_visc = mu0 * ((temperature_imp / T0)**1.5) * (((T0 + 198.72)) / (temperature_imp + 198.72))
    dyn_visc = dyn_visc * 47.8802687 # convert to kg / (m*s)
    Re_body = (velocity * rocket_length * density) / dyn_visc # dimensionless
    Rcrit = 51*(Rs/rocket_length) ** -1.039 # Barrowman eq 4-7
    if (Re_body <= Rcrit) & (Re_body > 0):
        cffb = 1.328 / Re_body**(1/2) # laminar flow, Barrowman 4-2
    elif (Re_body > 5*10**5) & (Re_body < Rcrit):
        # B = Rcrit * ((0.074 / Re_fin**(1/5)) - 1.328 / Re_fin**(1/2))
        # cffn = (0.074 / (Re_fin**(1/5))) - (B / Re_fin)
        cffb = ((3.46 * mth.log10(Re_body) - 5.6)**(-2)) - (1700 / Re_body) # First term is the Schoenherr theoretical boundary layer modification of von Karman's work, and the second term is the 
        # Prandtl correction term. Barrowman equation 4-6.
    elif Re_body > Rcrit: 
        cffb = 0.032 * (Rs / rocket_length) ** (0.2) # Turbulent curve match, takes into account
    else:
        cffb = 1e-10


    if mach > 0:
        cd_body = cffb * (rocket_length / rocket_diameter) * (mach / (q * rocket_length))**0.2 #(1 + (60 / (rocket_length / rocket_diameter)**3) + (0.0025 * ((rocket_length - 0.6477) / rocket_diameter)) * (2.7 * (0.508 / rocket_diameter)) + (4 * ((rocket_length - 0.6477) / rocket_diameter)) + (2 - 2 * (rocket_diameter / 0.10795)) * (0.1397 / 0.10795)) * cffb
    else: 
        cd_body= 0.0001
    cd_base = 0.12 + 0.13 * mach**2
    if t < tburn:
        cd_base = (1 - (0.00282743338 / Sref)) * (0.12 + 0.13 * mach**2) #TODO: Change 0.002... to fit the nozzle output area
      
    return cd_body, cd_base




def main():
    return 0
    







# def Drag(Rocket, mach, Re, Rcrit, Rs, t, tburn):
#     '''Inputs: 
#             velocity: m/s
#             R: Reynold's Number
#             Rs: Surface Finish'''
#     # Skin Friction
#     if (Re < 10**4) & (Re < Rcrit):
#         Cf = 1.48e-2
#     if (Re > 10**4) & (Re < Rcrit):
#         Cf = 1/ ((1.50*np.log(Re)-5.6)**2)
#     if Re > Rcrit:
#         Cf = 0.032*(Rs / Rocket.length)**0.2

#     #Subsonic correction:
#     Cfc = Cf*(1-0.1*mach**2)
    
#     # CDFriction:
#     Cdf = Cfc * (1 + (1/(2*(Rocket.length/Rocket.diameter))) * Rocket.Awet + (1+ 2*Rocket.finthickness/Rocket.Semispan) * Rocket.Awetfins) / Rocket.reference_area

#     ###################################################################################################
#     # Pressure Drag:
#     # Nosecone
#     Cdp_Nosecone = 0 #Smooth transition
#     # Boattail
#     Cdp_Boattail = (3 - (5.5/(6.16-4.25))) / 2

#     # Blunt fins: 
#     Cdp_LE = abs((1-mach**2))**(-0.417) - 1# Leading edge
    
#     Cdp_TE = 0.12 + 0.13*mach**2 #Trailing edge
#     Cdp_Fins = (Cdp_LE + Cdp_TE) * (Rocket.Fins_Aref / Rocket.reference_area)
    
#     ###################################################################################################
#     # Basedrag - Need to add correction for reduced base drag during motor burn

#     Cdb = (0.12 + 0.13*mach**2) * 0.6 # *0.6 is the correction factor, needed to be corrected once we have a CAD model
    

#     # Assume no parasitic drag for now - receding rail buttons
#     Cdpar = 0
    
# ###################################################################################################
#     # Total CD:
    
#     CD = (Cdb + Cdp_Fins + Cdp_Boattail + Cdp_Nosecone) + Cdf
#     print(t, mach, Re, CD)
#     return CD, Cdb
