from ambiance import Atmosphere    

# Define the atmosphere:
def atmosphere_conditions(height, velocity, rocket_length):
    # Uses the ambiance library to get standard atmosphere values for altitudes throughout the flight
    atmosphere = Atmosphere(min(81000, max(height,-5000)))
    gravity = atmosphere.grav_accel # m/s^2
    temperature = atmosphere.temperature # Kelvins
    pressure = atmosphere.pressure # pascal
    density = atmosphere.density # kg/m^3
    a = atmosphere.speed_of_sound # m/s
    mach = velocity / a # Dimensionless
    return gravity, temperature, pressure, density, mach