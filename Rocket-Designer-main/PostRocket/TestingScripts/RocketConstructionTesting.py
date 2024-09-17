# Author: Alex Post
# Date: 7/15/2024
# Utilizes snippets of Luis's code for 3dof trajectory in initial version.
# Will update to 6dof and edit completely in later versions.

# Imports
import numpy as np

from Trajectory import *
from RocketComponents.RocketComponents import *
from time import time as t
t0 = t()
# Defining all materials used                     #### NOTE: UPDATE ALL MATERIAL AND CHEMICAL VALUES
aluminum = material(2700,310e6,2e-4)
fiberglass = sheet_material(1.23,1.20e-3,12e6,2e-4)
carbon_fiber = sheet_material(6.274,4.06e-3,12e6,2e-4)
nomex = sheet_material(0.07796, 5e-3, None, None)
nylon = sheet_material(0.207, 5e-4,None, None) # Based on FROG parachute area w/o spillhole and mass
#kevlar = material()

nitrous =  propellant(1220,5.17e6)
paraffin = propellant(880,nitrous.pressure*0.75) ######################### PMMA NOT PARAFFIN 

# Rocket Component Instantiation in location order

################################################################################
antenna_mass = 0.1 # [kg]
antenna_location = 0.4826 # [m]
antennas = internal_component(antenna_mass,0,0,0,antenna_location)

nosecone_L = 0.508 # [m]
rocket_body_D = 0.156 # [m]
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
payload_location = 0.366 #[m]
payload = internal_component(payload_mass,0,0,0,payload_location)

co2_system_mass = 0.6 # [kg]
co2_system_location = 0.6858 # [m]
co2_system = internal_component(co2_system_mass,0,0,0,co2_system_location)

parachute_packed_L = 0
parachute_packed_D = 0.141 # [m]
parachute_location = 1.0922 - parachute_packed_L/2 # [m] # should be outside the body tube
parachute_cd = 1.4
main_parachute = parachute(parachute_packed_L,parachute_packed_D,parachute_location,nylon,parachute_cd)

middle_body_L = 1.13983 # [m]

middle_body = body_tube(middle_body_L,rocket_body_D,rocket_body_t,carbon_fiber, payload,co2_system,main_parachute)
################################################################################
oxidizer_L_guess = 0.5 # [m]

ox_tank = oxidizer_tank(oxidizer_L_guess,rocket_body_D,aluminum,nitrous)
################################################################################
lower_body_L = 1.0668 # [m]

lower_body = body_tube(lower_body_L,rocket_body_D,rocket_body_t,carbon_fiber)

grain_D = 0.089 # [m]
prechamber_L = 0.1 # [m]
postchamber_L = 0.1 # [m]  ## NOTE: ADD REAL VALUES FROM RICHIE AND SCHOONER
chamber_pressure = nitrous.pressure * 0.75

motor = combustion_chamber(lower_body_L,grain_D,prechamber_L,postchamber_L,aluminum,paraffin)
lower_body.add_subcomponents(motor)
################################################################################
fin_number = 3
fin_height = 0.14605 # [m]
fin_root_chord = 0.254 # [m]
fin_tip_chord = 0.1143 # [m]
fin_sweep_length = 0.07551 # [m]
fin_thickness = 6.35e-3 # [m]

fins = composite_fin_group(fin_number,fin_height,fin_root_chord,fin_tip_chord,fin_sweep_length,fin_thickness,carbon_fiber,nomex,rocket_body_D,1)
################################################################################
graphite_nozzle = nozzle(5e-5,7e-5,1,0.1,0.075,0)

boattail_L = 0.1397 # [m]
boattail_D2 = 0.10795 # [m]
boattail_t = 2.29e-3 # [m]

boattail = von_karman_boattail(boattail_L,rocket_body_D,boattail_D2,boattail_t,carbon_fiber,nozzle)
################################################################################
t1 = t()
TNAM_Test = rocket(nosecone,upper_body,middle_body,ox_tank,lower_body,fins,boattail)
t2 = t()
TNAM_Test.update_components(0.025,1,0.02)
t3 = t()
print("Instantiation Time: " + str(t1 - t0))
print("Collection Time: " + str(t2 - t1))
print("Update Time: " + str(t3 - t2))
length_test_1 = nosecone_L + upper_body_L + middle_body_L + boattail_L

# Defining a FROG-like rocket from base components
'''
class Rocket:
    pass

Rocket.length = 3.3673288 # m
Rocket.diameter = 0.078232*2
Rocket.mass = 27.65 # kg
Rocket.cg = 2.102358 # m
Rocket.nozzle_area = 0.0025 # m^2
Rocket.reference_area = np.pi * 0.156464**2
Rocket.Awet = 1.51309 #Approximate if the whole rocket was a cylinder
Rocket.finthickness = 0.00635 # m
Rocket.Semispan = 0.146 # m
Rocket.Awetfins = 0.027 * 6 # m^2
Rocket.Fins_Aref = 0.01905 # m^2


rail_length = 9.144 # m, 30 ft
burn_time = 4 #seconds
Weq = 3350 #Test input
m_dot = 0.53 #kg/s

# Rocket = rocket(nosecone, upper_body, middle_body, lower_body, fins)
# print(Rocket.mass)


# Test inputs
apogee, amax, mass_curve, rocket_alt, rocket_vel, velexact, rocket_time,rocket_dist,traj_angle,rocket_accel,CD_Bank, Mach_Bank, Cdb_Bank = Trajectory(TNAM_Test, rail_length, burn_time, Weq, m_dot)




for i in range(len(rocket_alt)):
    if apogee == rocket_alt[i]:
        apogee_time = i


initial_height = rocket_alt[4]
initial_velocity = rocket_vel[4] * np.cos(traj_angle[4])

h = initial_height + initial_velocity*(rocket_time-burn_time) - 0.5 * 9.81 * (rocket_time-burn_time)**2
h[0:4] = 0.5*(Weq*m_dot/Rocket.mass-9.81)*rocket_time[0:4]**2




import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
print(len(rocket_time))
fig,ax1 = plt.subplots()
ax1.plot(Mach_Bank[0:-1], Cdb_Bank[0:-1])


fig, ax = plt.subplots()

ax.plot(rocket_time, rocket_alt*3.28084, label='Numerical Model', linewidth=3,color='r')
ax.plot(rocket_time,h*3.28084,label='Analytical Model', linewidth=3,color='b', linestyle='--')
ax.set_title('No-Drag Trajectory', fontsize=20)
ax.axvline(x=burn_time, color='red', linestyle='--')
ax.axvline(x=apogee_time, color='red', linestyle='--')
ax.text(burn_time+.25, 3350*3.28084, 'Burn Out', color='red', fontsize=14, rotation=0)
ax.text(apogee_time+.25, 3350*3.28084, 'Apogee', color='red', fontsize=14, rotation=0)
majorLocator = MultipleLocator(5)
minorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)
majorLocatory = MultipleLocator(1000)
minorLocatory = MultipleLocator(100)
ax.yaxis.set_major_locator(majorLocatory)
ax.yaxis.set_minor_locator(minorLocatory)
ax.tick_params(which='major',length=5,width=2,labelsize=14)
ax.grid()
ax.set_xlabel('Time (s)', fontsize=15)
ax.set_ylabel('Altitude (ft)', fontsize=15)
ax.set_xlim(0,11*burn_time)
ax.set_ylim(0,3550*3.28084)
ax.legend(fontsize=15,loc="upper right")

# fig2,ax2 = plt.subplots()
# ax2.plot(rocket_time[1:],rocket_accel)

plt.show()
print(apogee*3.28084)




'''
