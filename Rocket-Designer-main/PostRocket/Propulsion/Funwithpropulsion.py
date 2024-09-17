

def TNAM(grain_diameter, port_diameter , grain_length, oxidizer_volume):
    # TEMPORARY FAKE BALLISTICS MODEL DO NOT KEEP
    # TODO: REPLACE THIS BY IMPORTING UPDATED TNAM 
    burn_time = [0, 6, 6.1]
    chamber_pressure_curve = [5e6, 5e6, 0]
    oxidizer_mass_curve = [24, 0]  
    fuel_mass_curve = [3, 0]
    thrust_curve = [oxidizer_volume*4000000,oxidizer_volume*2000000,oxidizer_volume*2000000,oxidizer_volume*1000000,0,0,0,0, 0]

    return burn_time, chamber_pressure_curve, oxidizer_mass_curve, fuel_mass_curve, thrust_curve