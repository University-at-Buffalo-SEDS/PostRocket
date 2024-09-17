# Drag functions as defined by Box, Bishop, and Hunt
import numpy as np

def cd(rocket,t,q,mach, velocity, density, temperature, Rs, AoA):

    # Find the reynolds numbers for each section of the rocket.
    Aref = rocket.reference_area
    body_length = rocket.length - (rocket.boattail.length + rocket.nosecone.length) # length of the cylindrical body [m]
    #Re_crit = 5e5 # critical reynolds assumption (Mandell)
    Re_crit =  Rcrit = 51*(Rs/rocket.length) ** -1.039
    Re_body = Re(temperature, density, rocket.length, velocity)
    Re_fins = Re(temperature, density, rocket.fins.midchord, velocity)
    

    # Calculate the viscous friction coefficients
    Cf_body = Cf(Re_body, Re_crit)
    Cf_fins = Cf(Re_fins, Re_crit)

    # Body drag 
    Cd_body = (1 + (60 * (rocket.length / rocket.diameter)**(-3))  + 0.0025 * (body_length / rocket.diameter)) * (2.7 * (rocket.nosecone.length) + 4 * (body_length) + (2 * (1 - (rocket.boattail.aft_diameter / rocket.boattail.fore_diameter)) * rocket.boattail.length)) * (Cf_body / rocket.diameter)



    # Base drag
    Cd_base = 0.029 * (((rocket.boattail.aft_diameter / rocket.boattail.fore_diameter) ** (3)) / (Cd_body ** (1/2)))

    # Fin drag and interference drag
    Cd_fandi = (0.5 * Cf_fins * (1 + 2 * (rocket.fins.thickness / rocket.fins.midchord)) * (4 * rocket.fins.number * (rocket.fins.planform_area + (rocket.diameter * rocket.fins.root_chord)))) / (rocket.reference_area)

    # Cd0
    Cd0 = Cd_body + Cd_base + Cd_fandi

    # AoA drag 
    delta = ( 7.0989*AoA + 0.1168)**0.5 # Relation from interpolation of data from Manell et al
    nu = 1.5355*AoA + 0.4953
    Cd_AoA_body = 2 * delta * (AoA ** (2)) + ((3.6 * nu * (1.36 * rocket.length - 0.55 * rocket.nosecone.length)) / (np.pi * rocket.diameter)) * AoA**3
    sec_ratio = (rocket.fins.height * 2) / (rocket.diameter) + 1
    kfb = 0.8065 * sec_ratio **2 + 1.553 * sec_ratio
    kbf = 0.1935 * sec_ratio**2 + 0.8174*sec_ratio + 1
    Cd_AoA_fins = (AoA**2) * ( 1.2 * ((rocket.fins.planform_area + (0.5 * rocket.fins.root_chord * rocket.diameter)) / rocket.reference_area) + 3.12 * (kfb + kbf - 1) * (rocket.fins.planform_area / rocket.reference_area) )
    
    # Total Cd:
    Cd = Cd0 + Cd_AoA_body + Cd_AoA_fins

    return Cd

def Re(temperature, density, length, velocity):
    #dyn_visc = 1.825e-5  #Value for air at 1 atm, 20c [kg/(m-s)]
    dyn_visc = 1.716e-5 * (temperature/273)**1.5 * 384 / (temperature + 111)
    Re = (velocity * length * density) / dyn_visc # dimensionless
    return Re

def Cf(Re, Re_crit):
    '''Works with scalar values of Re and Re_crit'''
    if (Re < Re_crit):
        Cf = 1.328 * Re**-0.5
    else:
        B = Re_crit * ((0.074 * Re**-0.2) - (1.328 * Re**-0.5))
        Cf = (0.074 * Re**-0.2) - (B / Re)
    return Cf


