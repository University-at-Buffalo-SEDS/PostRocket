# # Author: Alex Post
# # Date: 8/13/2024

# # Reference: https://openrocket.sourceforge.net/techdoc.pdf
# https://cambridgerocket.sourceforge.net/AerodynamicCoefficients.pdf 
# https://www.grc.nasa.gov/www/k-12/airplane/viscosity.html Air Viscosity
# #This function will calculate the forces occuring on a given rocket:

import numpy as np
import math as mth



def CNderivatives_nosecone(reference_area, plane_area, diameter):
    '''Inputs:
     plane area
     rocket reference area
     nosecone diameter '''
    K = 1.1 #For a von karman nosecone, it appears that you can approximate K = 3 per rad
    ArefCNa = 2 * reference_area # Normal force coefficient times the reference area, CN with partial respect to AOA
    ArefCNaa = K * plane_area * reference_area #Body lift coefficient times the reference area, CN with double partial respect to AOA

    return ArefCNa, ArefCNaa

def CNderivatives_bodytube(plane_area):
    '''Inputs:
    planiform area
    '''
    ArefCNa = 0 #Normal force coefficient is 0 for a bodytube
    K = 1.1
    ArefCNaa = K * (plane_area) #Aref * CNaa

    return ArefCNa, ArefCNaa

def CNderivatives_boattail(plane_area, diameter_ratio, fore_diameter):
    '''Inputs:
    plane area
    diameter ratio (fore/aft)
    fore diameter'''
    K = 1.1
    ArefCNa = 2 * ( diameter_ratio**-2  - 1) * (1/4)*(mth.pi)*(fore_diameter**2) # ORIGINAL NOTE: normal force coefficient ONLY APPLICABLE FOR CONSTANT DIAMETER BODY
    ArefCNaa = K * plane_area
    return ArefCNa, ArefCNaa

def CNderivatives_fins(planform_area, number, semi_span, sweep_length, root_chord, tip_chord, parent_radius):
    '''Inputs:
    planform area
    fin dimensions
    ...
    rocket radius
    '''
    K = 1.1
    l = np.sqrt((semi_span)**2 + (sweep_length + (tip_chord - root_chord)/2)**2)
    kfb = 1 + parent_radius/(semi_span+parent_radius)
    ArefCNa = (kfb*( 4*number*(0.5*semi_span/parent_radius)**2 ) / (1 + np.sqrt( 1 + (2*l/(root_chord+tip_chord))**2 ))) * parent_radius**2 * mth.pi
    ArefCNaa = 0

    return ArefCNa, ArefCNaa




