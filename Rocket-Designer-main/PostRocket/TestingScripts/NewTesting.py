# Author: Diego Del Toro
# Date: 11/8/2024
# TESTING

# Make sure to import install numpy, scipy, ambiance, matlab, matlabengine
# Also hrap, currently idk what is happening

print('importing modules...')
import numpy as np
import scipy.optimize as opt
from PostRocket.FlightSim.Traj_6dof import Traj_6dof
from PostRocket.RocketComponents.RocketComponents import *
import time
import matplotlib.pyplot as plt
from winsound import Beep



print('Initializing Components...')

# Defining all materials used                     #### NOTE: UPDATE ALL MATERIAL AND CHEMICAL VALUES
aluminum = material(2700,310e6,2e-4)
fiberglass = sheet_material(1.23,1.20e-3,12e6,2e-4)
carbon_fiber = sheet_material(6.274,4.06e-3,12e6,2e-4)
nomex = sheet_material(0.07796, 5e-3, None, None)
nylon = sheet_material(0.207, 5e-4,None, None) # Based on FROG parachute area w/o spillhole and mass
#kevlar = material() 

# Rocket Component Instantiation in location order

### TEMPORARY#####
scaleFactor = 1

################################################################################
antenna_mass = 0.1 # [kg]
antenna_location = 0.4826 - 0.125 # [m]
antennas = internal_component(antenna_mass,0.125,0.25,0.1,antenna_location)

nosecone_L = 0.508 # [m]
rocket_body_D = 0.156 * scaleFactor# [m]
nosecone_t = 2.26e-3 # [m]
nosecone_tip_L = 0.0508 # [m]

nosecone = tipped_nosecone(nosecone_L,rocket_body_D,nosecone_t,fiberglass,nosecone_tip_L,aluminum, antennas)
################################################################################
avionics_mass = 1.55 # [kg]
avionics_location = 0.226 # [m]
avionics = internal_component(avionics_mass,0,0,0,avionics_location)

upper_body_L = 0.508 # [m]
rocket_body_t = 4.06e-3 # [m]

upper_body = body_tube(upper_body_L,rocket_body_D,rocket_body_t,fiberglass,  avionics)
################################################################################
payload_mass = 3.89 # [kg]
payload_location = 0.366 - 0.125 #[m]
payload = internal_component(payload_mass,0.125,0.25,0.1,payload_location)

co2_system_mass = 0.6 # [kg]
co2_system_location = 0.6858 - 0.025 # [m]
co2_system = internal_component(co2_system_mass,0.025,0.05,0.141,co2_system_location)

parachute_packed_L = 0.5
parachute_packed_D = 0.141 # [m]
parachute_location = 1.0922 - parachute_packed_L # [m] # should be outside the body tube
parachute_cd = 1.4
main_parachute = parachute(parachute_packed_L,parachute_packed_D,parachute_location,nylon,parachute_cd)

middle_body_L = 1.13983 # [m]

middle_body = body_tube(middle_body_L,rocket_body_D,rocket_body_t,carbon_fiber, payload,co2_system,main_parachute)
################################################################################
oxidizer_L_guess = 0.5 # [m]

ox_tank = oxidizer_tank(oxidizer_L_guess,rocket_body_D,aluminum)
################################################################################
lower_body_L = 1.0668 # [m]

lower_body = body_tube(lower_body_L,rocket_body_D,rocket_body_t,carbon_fiber)

grain_D = 0.089 # [m]
prechamber_L = 0.1 # [m]
postchamber_L = 0.1 # [m]  ## NOTE: ADD REAL VALUES FROM RICHIE AND SCHOONER
chamber_pressure = 0.75

motor = combustion_chamber(lower_body_L,grain_D,prechamber_L,postchamber_L,aluminum)
lower_body.add_subcomponents(motor)
################################################################################
fin_number = 3
fin_height = 0.14605 # [m]
fin_root_chord = 0.254 # [m]
fin_tip_chord = 0.1143 # [m]
fin_sweep_length = 0.07551 # [m]
fin_thickness = 6.35e-3 # [m]

fins = composite_fin_group(fin_number,fin_height,fin_root_chord,fin_tip_chord,fin_sweep_length,fin_thickness,carbon_fiber,nomex,rocket_body_D,0.0254)
################################################################################
graphite_nozzle = nozzle(5e-5,7e-5,1,0.1,0.075,0)

boattail_L = 0.1397 # [m]
boattail_D2 = 0.10795# [m]
boattail_t = 2.29e-3 # [m]

boattail = von_karman_boattail(boattail_L,rocket_body_D,boattail_D2,boattail_t,carbon_fiber,graphite_nozzle)
################################################################################

print('Initializing Rocket...')
TNAM_Test = rocket(nosecone,upper_body,middle_body,ox_tank,lower_body,fins,boattail)

TNAM_Test.update_components(0.01335, 0.6, 0.007)

rail_length = 10
launch_angle = 86
requested_apogee = 3352.8 #m, 11000 ft
requested_off_rail = 30.48 # m/s, 100 ft/s
stability = 1.5



# Test inputs
traj = Traj_6dof(TNAM_Test, rail_length, launch_angle, 'LSODA', 1e-4, 1e-4, False)

plt.close()
fig, ax = plt.subplots(2,1)
TNAM_Test.draw(ax[0])
ax[1].plot(motor.burn_time,motor.thrust_curve,label='Thrust')
ax[1].plot(traj['time'],traj['z_position'],label='Altitude')
ax[1].grid(True,'major',alpha=0.5)
#ax[1].xaxis.set_minor_locator(AutoMinorLocator())
#ax[1].yaxis.set_minor_locator(AutoMinorLocator())
#ax[1].grid(True,'minor', alpha= 0.25)
ax[1].legend()
ax[1].set_xlabel('time [s]')
ax[1].set_ylabel('Thrust [N], Altitude [m]')
plt.show()