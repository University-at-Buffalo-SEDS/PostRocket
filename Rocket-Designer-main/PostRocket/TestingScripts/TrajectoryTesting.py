# Author: Alex Post
# Date: 7/15/2024
# Utilizes snippets of Luis's code for 3dof trajectory in initial version.
# Will update to 6dof and edit completely in later versions.

# Imports
import numpy as np
import scipy as sc
from Trajectory import *
from RocketComponents.RocketComponents import *
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import csv 

# Adding some input parameters
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


method = ["RK45"]
error = 1e-15
traj = Trajectory(Rocket, rail_length, burn_time, Weq, m_dot,"RK45", error, error, False)
trajectory = traj[3]
velocity = traj[4]
acceleration = traj[8]
traj_time = traj[5]
Truth_alt = trajectory
Truth_time = traj_time
Truth_vel = velocity
Truth_acc = acceleration




# Test inputs
N = 100
apogee_avg = []
rail_departure_time_avg = []
rail_departure_velocity_avg = []
rocket_alt_avg = []
rocket_time_avg = []
traj_angle_avg = []
Runtime_avg = []
rocket_alt_all = []
rocket_time_all = []

method = ["RK45", "RK23", "LSODA"]

errorlist = np.logspace(-6,0,500)
# Removed Radau and BDF for now since vectorized isnt working.
MSE_alt = np.zeros([N, len(errorlist), len(method)])
MSE_vel = np.zeros([N, len(errorlist), len(method)])
MSE_acc = np.zeros([N, len(errorlist), len(method)])
MSE_alt_avg = np.zeros([len(errorlist), len(method)])
MSE_vel_avg = np.zeros([len(errorlist), len(method)])
MSE_acc_avg = np.zeros([len(errorlist), len(method)])
MSE_alt_dev = np.zeros([len(errorlist), len(method)])
MSE_vel_dev = np.zeros([len(errorlist), len(method)])
MSE_acc_dev = np.zeros([len(errorlist), len(method)])
Runtime_avg = np.zeros([len(errorlist), len(method)])
Runtime_dev = np.zeros([len(errorlist), len(method)])


# Altitude = np.zeros([N, len(errorlist), len(method)])
Runtime = np.zeros([N, len(errorlist), len(method)])

# method
for k,methodname in enumerate(method):
    apogee = []
    rail_departure_time = []
    rail_departure_velocity = []
    rocket_alt = []
    rocket_time = []
    traj_angle = []
    

    print(f"Starting calculations for method: {methodname}")
    for i in range(N):
        print('Run #'+ str(i),end='\r')
        for j,error in enumerate(errorlist):
            print('Run #' + str(i) +' Checking error #'+ str(j),end='\r')
            start = time.perf_counter()
            traj = Trajectory(Rocket, rail_length, burn_time, Weq, m_dot,methodname, error, error, False)
            end = time.perf_counter()

            trajectory = traj[3]
            velocity = traj[4]
            acceleration = traj[8]
            traj_time = traj[5]
            
            #Save the altitude data for the first run:


            adj_Truth_alt = np.interp(traj_time, Truth_time, Truth_alt) #adjusted to only give data at the specific timesteps of the lower precision simulation
            adj_Truth_vel = np.interp(traj_time, Truth_time, Truth_vel) #adjusted to only give data at the specific timesteps of the lower precision simulation
            adj_Truth_acc = np.interp(traj_time, Truth_time, Truth_acc) #adjusted to only give data at the specific timesteps of the lower precision simulation

            # Runtime
            Runtime[i][j][k] = end - start

            # MSE
            not_normalized_alt = np.zeros(len(adj_Truth_alt))
            not_normalized_vel = np.zeros(len(adj_Truth_vel))
            not_normalized_acc = np.zeros(len(adj_Truth_acc))

            for n,m in enumerate(adj_Truth_alt):
                dev = trajectory[n] - m
                sq_dev = dev**2
                not_normalized_alt[n] = sq_dev
            normalized_alt = sum(not_normalized_alt) / len(adj_Truth_alt)
            MSE_alt[i][j][k] = normalized_alt
            
            for n,m in enumerate(adj_Truth_vel):
                dev = velocity[n] - m
                sq_dev = dev**2
                not_normalized_vel[n] = sq_dev
            normalized_vel = sum(not_normalized_vel) / len(adj_Truth_vel)
            MSE_vel[i][j][k] = normalized_vel

            for n,m in enumerate(adj_Truth_acc):
                dev = acceleration[n] - m
                sq_dev = dev**2
                not_normalized_acc[n] = sq_dev
            normalized_acc = sum(not_normalized_acc) / len(adj_Truth_acc)
            MSE_acc[i][j][k] = normalized_acc
            
    
    MSE_alt_avg[:,k] = np.mean(MSE_alt[:,:,k], axis=0)
    MSE_vel_avg[:,k] = np.mean(MSE_vel[:,:,k], axis=0)
    MSE_acc_avg[:,k] = np.mean(MSE_acc[:,:,k], axis=0)
    MSE_alt_dev[:,k] = np.std(MSE_acc[:,:,k], axis=0)
    MSE_vel_dev[:,k] = np.std(MSE_acc[:,:,k], axis=0)
    MSE_acc_dev[:,k] = np.std(MSE_acc[:,:,k], axis=0)
    Runtime_avg[:,k] = np.mean(Runtime[:,:,k], axis=0)
    Runtime_dev[:,k] = np.std(Runtime[:,:,k], axis=0)


with open("MSE_alt_avg.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(MSE_alt_avg)

with open("MSE_vel_avg.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(MSE_vel_avg)

with open("MSE_acc_avg.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(MSE_acc_avg)

with open("MSE_alt_dev.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(MSE_alt_dev)

with open("MSE_vel_dev.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(MSE_vel_dev)

with open("MSE_acc_dev.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(MSE_acc_dev)

with open("Runtime_avg.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Runtime_avg)

with open("Runtime_dev.csv",'w',newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Runtime_dev)


# print(errorlist)
# print(MSE_alt_avg[:,0])
# fig4, ax4 = plt.subplots()
# ax4.plot(errorlist, MSE_alt_avg[:,0])
# ax4.set_xscale("log", base=10)
# ax4.set_yscale("log", base=10)

# ax4.fill_between(rk45_rocket_time, rk45_all_alt.mean(axis=0) - rk45_all_alt.std(axis=0), rk45_all_alt.mean(axis=0) + rk45_all_alt.std(axis=0), color='#888888', alpha=0.4)
# ax4.fill_between(rk45_rocket_time, rk45_all_alt.mean(axis=0) - 2*rk45_all_alt.std(axis=0), rk45_all_alt.mean(axis=0) + 2*rk45_all_alt.std(axis=0), color='#888888', alpha=0.2)
# ax4.set_title('Rocket Trajectory RK45', fontsize=20)
# ax4.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax4.axvline(x=burn_time, color='red', linestyle='--')
# ax4.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax4.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax4.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax4.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax4.xaxis.set_major_locator(majorLocator)
# ax4.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax4.yaxis.set_major_locator(majorLocatory)
# ax4.yaxis.set_minor_locator(minorLocatory)
# ax4.tick_params(which='major',length=5,width=2,labelsize=14)
# ax4.grid()
# ax4.set_xlabel('Time (s)', fontsize=15)
# ax4.set_ylabel('Altitude (ft)', fontsize=15)
# ax4.set_xlim([1e-12, max(errorlist) + .5])
# ax4.set_ylim(-1000,max(MSE_alt_avg[:,0]))
# # ax4.legend(fontsize=11,loc="lower right")
# plt.show()
    
       





            
    
        # fig9, ax9 = plt.subplots()
        # ax9.plot(traj_time, adj_Truth_alt* 3.28084 , alpha=0.9, color='orange',label='Truth', linewidth = 1.0)
        # ax9.plot(traj_time, trajectory* 3.28084 , alpha=0.9, color='green',label=str(error), linewidth = 1.0)
        
        # ax9.set_title('Rocket Trajectory Comparison', fontsize=20)
        # ax9.axvline(x=traj[10], color='red', linestyle='--')
        # ax9.axvline(x=burn_time, color='red', linestyle='--')
        # ax9.axvline(x=traj_time[-1], color='red', linestyle='--')
        # ax9.text(burn_time+.25, max(Truth_alt[0])+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
        # ax9.text(traj_time[-1]+.25, max(Truth_alt[0])+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
        # ax9.text(traj[10]+.25, max(Truth_alt[0])+ 100, 'Departure', color='red', fontsize=12, rotation=0)
        # majorLocator = MultipleLocator(1)
        # minorLocator = MultipleLocator(.5)
        # ax9.xaxis.set_major_locator(majorLocator)
        # ax9.xaxis.set_minor_locator(minorLocator)
        # majorLocatory = MultipleLocator(1000)
        # minorLocatory = MultipleLocator(100)
        # ax9.yaxis.set_major_locator(majorLocatory)
        # ax9.yaxis.set_minor_locator(minorLocatory)
        # ax9.tick_params(which='major',length=5,width=2,labelsize=14)
        # ax9.grid()
        # ax9.set_xlabel('Time (s)', fontsize=15)
        # ax9.set_ylabel('Altitude (ft)', fontsize=15)
        # ax9.set_xlim(0,20)
        # ax9.set_ylim(0,max(Truth_alt[0]* 3.28084 )+400)
        # ax9.legend(fontsize=11,loc="lower right")
        # plt.show()


        # Newlength = 100
        # rocket_alt_rk45_adj = np.interp(np.linspace(min(rk45_rocket_time), max(rk23_rocket_time), Newlength), rk45_rocket_time, rocket_alt_rk45)
        





#         Runtime.append(end - start)
#         apogee.append(traj[0])
#         rail_departure_time.append(traj[9])
#         rail_departure_velocity.append(traj[10])
#         rocket_alt.append(traj[3])
#         rocket_time.append(traj[5])
#         traj_angle.append(traj[7])
        

#         #print(f"Loop {i} Completed")
    
#     #Change into a numpy array
#     print(rocket_alt)
#     apogee = np.array(apogee)
#     rail_departure_time = np.array(rail_departure_time)
#     rail_departure_velocity = np.array(rail_departure_velocity)
#     rocket_alt = np.array(rocket_alt)
#     rocket_time = np.array(rocket_time)
#     traj_angle = np.array(traj_angle)
    
#     #Save individual runs:
#     rocket_alt_all_temp = np.array(rocket_alt)
#     rocket_alt_all.append(rocket_alt_all_temp)
#     rocket_time_all_temp = np.array(rocket_time)
#     rocket_time_all.append(rocket_time_all_temp)


#     #Find the numpy mean
#     at = np.mean(apogee)
#     rdt = np.mean(rail_departure_time)
#     rdv = np.mean(rail_departure_velocity)
#     raa = np.mean(rocket_alt,axis=0)
#     rta = np.mean(rocket_time,axis=0)
#     taa = np.mean(traj_angle,axis=0)
#     ra = np.mean(Runtime)

#     #Convert back into a list
#     at = at.tolist()
#     rdt = rdt.tolist()
#     rdv = rdv.tolist()
#     raa = raa.tolist()
#     rta = rta.tolist()
#     taa = taa.tolist()
#     ra = ra.tolist()

#     #Append onto a blank list to extract data
#     apogee_avg.append(at) 
#     rail_departure_time_avg.append(rdt) 
#     rail_departure_velocity_avg.append(rdv) 
#     rocket_alt_avg.append(raa) 
#     rocket_time_avg.append(rta) 
#     traj_angle_avg.append(taa) 
#     Runtime_avg.append(ra)



# #rk45 data
# rk45_apogee = apogee_avg[0]
# rk45_rail_departure_time = rail_departure_time_avg[0]
# rk45_rail_departure_velocity = rail_departure_velocity_avg[0]
# rk45_rocket_alt = rocket_alt_avg[0]
# rk45_rocket_alt = [rk45_rocket_alt * 3.28084 for rk45_rocket_alt in rk45_rocket_alt]
# rk45_rocket_time = rocket_time_avg[0]
# rk45_traj_angle = traj_angle_avg[0]
# rk45_traj_angle = [rk45_traj_angle * 57.2958 for rk45_traj_angle in rk45_traj_angle]
# rk45_runtime = Runtime_avg[0]
# rk45_all_alt = rocket_alt_all[0] * 3.28084
# rk45_all_time = rocket_time_all[0]

# #rk23 data
# rk23_apogee = apogee_avg[1]
# rk23_rail_departure_time = rail_departure_time_avg[1]
# rk23_rail_departure_velocity = rail_departure_velocity_avg[1]
# rk23_rocket_alt = rocket_alt_avg[1]
# rk23_rocket_alt = [rk23_rocket_alt * 3.28084 for rk23_rocket_alt in rk23_rocket_alt]
# rk23_rocket_time = rocket_time_avg[1]
# rk23_traj_angle = traj_angle_avg[1]
# rk23_traj_angle = [rk23_traj_angle * 57.2958 for rk23_traj_angle in rk23_traj_angle]
# rk23_runtime = Runtime_avg[1]
# rk23_all_alt = rocket_alt_all[1] * 3.28084
# rk23_all_time = rocket_time_all[1]

# # DOP853 data
# DOP853_apogee = apogee_avg[2]
# DOP853_rail_departure_time = rail_departure_time_avg[2]
# DOP853_rail_departure_velocity = rail_departure_velocity_avg[2]
# DOP853_rocket_alt = rocket_alt_avg[2]
# DOP853_rocket_alt = [DOP853_rocket_alt * 3.28084 for DOP853_rocket_alt in DOP853_rocket_alt]
# DOP853_rocket_time = rocket_time_avg[2]
# DOP853_traj_angle = traj_angle_avg[2]
# DOP853_traj_angle = [DOP853_traj_angle * 57.2958 for DOP853_traj_angle in DOP853_traj_angle]
# DOP853_runtime = Runtime_avg[2]
# DOP853_all_alt = rocket_alt_all[2] * 3.28084
# DOP853_all_time = rocket_time_all[2]

# #LSODA data
# LSODA_apogee = apogee_avg[3]
# LSODA_rail_departure_time = rail_departure_time_avg[3]
# LSODA_rail_departure_velocity = rail_departure_velocity_avg[3]
# LSODA_rocket_alt = rocket_alt_avg[3]
# LSODA_rocket_alt = [LSODA_rocket_alt * 3.28084 for LSODA_rocket_alt in LSODA_rocket_alt]
# LSODA_rocket_time = rocket_time_avg[3]
# LSODA_traj_angle = traj_angle_avg[3]
# LSODA_traj_angle = [LSODA_traj_angle * 57.2958 for LSODA_traj_angle in LSODA_traj_angle]
# LSODA_runtime = Runtime_avg[3]
# LSODA_all_alt = rocket_alt_all[3] * 3.28084
# LSODA_all_time = rocket_time_all[3]



# # rocket_alt_rk45 = rk45_all_alt.mean(axis=0)
# # rocket_alt_rk23 = rk23_all_alt.mean(axis=0)
# # rocket_alt_DOP853 = DOP853_all_alt.mean(axis=0)
# # rocket_alt_LSODA = LSODA_all_alt.mean(axis=0)
# # rocket_alt_rk45 = rocket_alt_rk45.tolist()
# # rocket_alt_rk23 = rocket_alt_rk23.tolist()
# # rocket_alt_LSODA = rocket_alt_LSODA.tolist()
# # rocket_alt_DOP853 = rocket_alt_DOP853.tolist()

# # Newlength = 100
# # rocket_alt_rk45_adj = np.interp(np.linspace(min(rk45_rocket_time), max(rk23_rocket_time), Newlength), rk45_rocket_time, rocket_alt_rk45)
# # rocket_alt_rk23_adj = np.interp(np.linspace(min(rk45_rocket_time), max(rk23_rocket_time), Newlength), rk23_rocket_time, rocket_alt_rk23)
# # rocket_alt_DOP853_adj = np.interp(np.linspace(min(rk45_rocket_time), max(rk23_rocket_time), Newlength), DOP853_rocket_time, rocket_alt_DOP853)
# # rocket_alt_LSODA_adj = np.interp(np.linspace(min(rk45_rocket_time), max(rk23_rocket_time), Newlength), LSODA_rocket_time, rocket_alt_LSODA)
# # rocket_time_avg_adj = np.interp(np.linspace(min(rk45_rocket_time), max(rk23_rocket_time), Newlength), rk45_rocket_time, rk45_rocket_time)



# # rocket_alt = np.array([rocket_alt_rk45_adj, rocket_alt_rk23_adj, rocket_alt_DOP853_adj, rocket_alt_LSODA_adj], dtype=np.ndarray)
# # rocket_alt = rocket_alt.astype(float)
# # mean_alt = rocket_alt.mean(axis=0)
# # stdd_alt = rocket_alt.std(axis=0)



# # dt algorithm
# dt_rk45 = np.zeros(len(rk45_rocket_time) -1)
# dt_rk23 = np.zeros(len(rk23_rocket_time) -1)
# dt_DOP853 = np.zeros(len(DOP853_rocket_time) -1)
# dt_LSODA = np.zeros(len(LSODA_rocket_time) -1)

# for i in range(len(dt_rk45)): 
#     dt_rk45[i] = rk45_rocket_time[i + 1] - rk45_rocket_time[i]
# for i in range(len(dt_rk23)):
#     dt_rk23[i] = rk23_rocket_time[i + 1] - rk23_rocket_time[i]
# for i in range(len(dt_DOP853)):
#     dt_DOP853[i] = DOP853_rocket_time[i + 1] - DOP853_rocket_time[i]
# for i in range(len(dt_LSODA)):
#     dt_LSODA[i] = LSODA_rocket_time[i + 1] - LSODA_rocket_time[i]

# rk45_apogee_time = rk45_rocket_time[-1]
# rk23_apogee_time = rk23_rocket_time[-1]
# DOP853_apogee_time = DOP853_rocket_time[-1]
# LSODA_apogee_time = LSODA_rocket_time[-1]

# avg_apogee_time = np.mean(np.array([rk45_apogee_time, rk23_apogee_time, DOP853_apogee_time, LSODA_apogee_time]))
# avg_departure_time = np.mean(np.array([rk45_rail_departure_time, rk23_rail_departure_time, DOP853_rail_departure_time, LSODA_rail_departure_time]))




# # Plot outputs
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# fig9, ax9 = plt.subplots()
# ax9.plot(rk45_rocket_time, rk45_all_alt.mean(axis=0), alpha=0.9, color='orange',label='RK45', linewidth = 1.0)
# ax9.plot(rk23_rocket_time, rk23_all_alt.mean(axis=0), alpha=0.9, color='green',label='RK23', linewidth = 1.0)
# ax9.plot(DOP853_rocket_time, DOP853_all_alt.mean(axis=0), alpha=0.9, color='lightblue',label='DOP853', linewidth = 1.0)
# ax9.plot(LSODA_rocket_time, LSODA_all_alt.mean(axis=0), alpha=0.9, color='black',label='LSODA', linewidth = 1.0)
# ax9.set_title('Rocket Trajectory All', fontsize=20)
# ax9.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax9.axvline(x=burn_time, color='red', linestyle='--')
# ax9.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax9.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax9.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax9.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(0.25)
# minorLocator = MultipleLocator(.1)
# ax9.xaxis.set_major_locator(majorLocator)
# ax9.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(10)
# minorLocatory = MultipleLocator(1)
# ax9.yaxis.set_major_locator(majorLocatory)
# ax9.yaxis.set_minor_locator(minorLocatory)
# ax9.tick_params(which='major',length=5,width=2,labelsize=14)
# ax9.grid()
# ax9.set_xlabel('Time (s)', fontsize=15)
# ax9.set_ylabel('Altitude (ft)', fontsize=15)
# ax9.set_xlim(0,20)
# ax9.set_ylim(0,max(mean_alt)+400)
# ax9.legend(fontsize=11,loc="lower right")



# fig4, ax4 = plt.subplots()
# ax4.plot(rk45_rocket_time, rk45_all_alt.mean(axis=0), alpha=0.9, color='orange',label='RK45', linewidth = 1.0)
# ax4.fill_between(rk45_rocket_time, rk45_all_alt.mean(axis=0) - rk45_all_alt.std(axis=0), rk45_all_alt.mean(axis=0) + rk45_all_alt.std(axis=0), color='#888888', alpha=0.4)
# ax4.fill_between(rk45_rocket_time, rk45_all_alt.mean(axis=0) - 2*rk45_all_alt.std(axis=0), rk45_all_alt.mean(axis=0) + 2*rk45_all_alt.std(axis=0), color='#888888', alpha=0.2)
# ax4.set_title('Rocket Trajectory RK45', fontsize=20)
# ax4.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax4.axvline(x=burn_time, color='red', linestyle='--')
# ax4.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax4.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax4.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax4.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax4.xaxis.set_major_locator(majorLocator)
# ax4.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax4.yaxis.set_major_locator(majorLocatory)
# ax4.yaxis.set_minor_locator(minorLocatory)
# ax4.tick_params(which='major',length=5,width=2,labelsize=14)
# ax4.grid()
# ax4.set_xlabel('Time (s)', fontsize=15)
# ax4.set_ylabel('Altitude (ft)', fontsize=15)
# ax4.set_xlim(0,20)
# ax4.set_ylim(0,max(mean_alt)+400)
# ax4.legend(fontsize=11,loc="lower right")

# fig5, ax5 = plt.subplots()
# ax5.plot(rk23_rocket_time, rk23_all_alt.mean(axis=0), alpha=0.9, color='green',label='RK23', linewidth = 1.0)
# ax5.fill_between(rk23_rocket_time, rk23_all_alt.mean(axis=0) - rk23_all_alt.std(axis=0), rk23_all_alt.mean(axis=0) + rk23_all_alt.std(axis=0), color='#888888', alpha=0.4)
# ax5.fill_between(rk23_rocket_time, rk23_all_alt.mean(axis=0) - 2*rk23_all_alt.std(axis=0), rk23_all_alt.mean(axis=0) + 2*rk23_all_alt.std(axis=0), color='#888888', alpha=0.2)
# ax5.set_title('Rocket Trajectory RK23', fontsize=20)
# ax5.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax5.axvline(x=burn_time, color='red', linestyle='--')
# ax5.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax5.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax5.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax5.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax5.xaxis.set_major_locator(majorLocator)
# ax5.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax5.yaxis.set_major_locator(majorLocatory)
# ax5.yaxis.set_minor_locator(minorLocatory)
# ax5.tick_params(which='major',length=5,width=2,labelsize=14)
# ax5.grid()
# ax5.set_xlabel('Time (s)', fontsize=15)
# ax5.set_ylabel('Altitude (ft)', fontsize=15)
# ax5.set_xlim(0,20)
# ax5.set_ylim(0,max(mean_alt)+400)
# ax5.legend(fontsize=12,loc="lower right")

# fig6, ax6 = plt.subplots()
# ax6.plot(DOP853_rocket_time, DOP853_all_alt.mean(axis=0), alpha=0.9, color='lightblue',label='DOP853', linewidth = 1.0)
# ax6.fill_between(DOP853_rocket_time, DOP853_all_alt.mean(axis=0) - DOP853_all_alt.std(axis=0), DOP853_all_alt.mean(axis=0) + DOP853_all_alt.std(axis=0), color='#888888', alpha=0.4)
# ax6.fill_between(DOP853_rocket_time, DOP853_all_alt.mean(axis=0) - 2*DOP853_all_alt.std(axis=0), DOP853_all_alt.mean(axis=0) + 2*DOP853_all_alt.std(axis=0), color='#888888', alpha=0.2)
# ax6.set_title('Rocket Trajectory DOP853', fontsize=20)
# ax6.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax6.axvline(x=burn_time, color='red', linestyle='--')
# ax6.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax6.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax6.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax6.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax6.xaxis.set_major_locator(majorLocator)
# ax6.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax6.yaxis.set_major_locator(majorLocatory)
# ax6.yaxis.set_minor_locator(minorLocatory)
# ax6.tick_params(which='major',length=5,width=2,labelsize=14)
# ax6.grid()
# ax6.set_xlabel('Time (s)', fontsize=15)
# ax6.set_ylabel('Altitude (ft)', fontsize=15)
# ax6.set_xlim(0,20)
# ax6.set_ylim(0,max(mean_alt)+400)
# ax6.legend(fontsize=12,loc="lower right")

# fig7, ax7 = plt.subplots()
# ax7.plot(LSODA_rocket_time, LSODA_all_alt.mean(axis=0), alpha=0.9, color='black',label='LSODA', linewidth = 1.0)
# ax7.fill_between(LSODA_rocket_time, LSODA_all_alt.mean(axis=0) - LSODA_all_alt.std(axis=0), LSODA_all_alt.mean(axis=0) + LSODA_all_alt.std(axis=0), color='#888888', alpha=0.4)
# ax7.fill_between(LSODA_rocket_time, LSODA_all_alt.mean(axis=0) - 2*LSODA_all_alt.std(axis=0), LSODA_all_alt.mean(axis=0) + 2*LSODA_all_alt.std(axis=0), color='#888888', alpha=0.2)
# ax7.set_title('Rocket Trajectory LSODA', fontsize=20)
# ax7.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax7.axvline(x=burn_time, color='red', linestyle='--')
# ax7.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax7.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax7.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax7.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax7.xaxis.set_major_locator(majorLocator)
# ax7.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax7.yaxis.set_major_locator(majorLocatory)
# ax7.yaxis.set_minor_locator(minorLocatory)
# ax7.tick_params(which='major',length=5,width=2,labelsize=14)
# ax7.grid()
# ax7.set_xlabel('Time (s)', fontsize=15)
# ax7.set_ylabel('Altitude (ft)', fontsize=15)
# ax7.set_xlim(0,20)
# ax7.set_ylim(0,max(mean_alt)+400)
# ax7.legend(fontsize=12,loc="lower right")

# #Testing against one another
# fig8, ax8 = plt.subplots()
# ax8.plot(rocket_time_avg_adj, mean_alt, alpha=0.6, color='r',label='Mean Altitude', linewidth = 2.0)
# ax8.fill_between(rocket_time_avg_adj, mean_alt - stdd_alt, mean_alt + stdd_alt, color='#888888', alpha=0.4)
# ax8.fill_between(rocket_time_avg_adj, mean_alt - 2*stdd_alt, mean_alt + 2*stdd_alt, color='#888888', alpha=0.2)
# ax8.set_title('Rocket Trajectory Comparison', fontsize=20)
# ax8.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax8.axvline(x=burn_time, color='red', linestyle='--')
# ax8.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax8.text(burn_time+.25, max(mean_alt)+ 100, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax8.text(avg_apogee_time+.25, max(mean_alt)+ 100, 'Apogee', color='red', fontsize=14, rotation=0)
# ax8.text(avg_departure_time+.25, max(mean_alt)+ 100, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax8.xaxis.set_major_locator(majorLocator)
# ax8.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax8.yaxis.set_major_locator(majorLocatory)
# ax8.yaxis.set_minor_locator(minorLocatory)
# ax8.tick_params(which='major',length=5,width=2,labelsize=14)
# ax8.grid()
# ax8.set_xlabel('Time (s)', fontsize=15)
# ax8.set_ylabel('Altitude (ft)', fontsize=15)
# ax8.set_xlim(0,20)
# ax8.set_ylim(0,max(mean_alt)+400)
# ax8.legend(fontsize=12,loc="lower right")


# # Time discretization density
# fig0, ax0 = plt.subplots()
# ax0.plot(rk45_rocket_time[0:-1], dt_rk45 , label='RK45', linewidth=3,color='orange')
# ax0.plot(rk23_rocket_time[0:-1], dt_rk23, label='RK23', linewidth=3,color='green')
# ax0.plot(DOP853_rocket_time[0:-1], dt_DOP853, label='DOP853', linewidth=3,color='lightblue')
# ax0.plot(LSODA_rocket_time[0:-1], dt_LSODA, label='LSODA', linewidth=3,color='k')
# ax0.set_xlabel('Time (s)', fontsize=15)
# ax0.set_ylabel("Timestep Size (s)", fontsize=15)
# ax0.set_title('Average Timestep Size', fontsize=20)
# ax0.axvline(x=burn_time, color='red', linestyle='--')
# ax0.axvline(x=avg_apogee_time, color='red', linestyle='--')
# ax0.axvline(x=avg_departure_time, color='red', linestyle='--')
# ax0.text(burn_time+.1, max(dt_DOP853)-1, 'Burn Out', color='red', fontsize=12, rotation=0)
# ax0.text(avg_apogee_time+.1, max(dt_DOP853)-1, 'Apogee', color='red', fontsize=12, rotation=0)
# ax0.text(avg_departure_time+.1, max(dt_DOP853)-1, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax0.xaxis.set_major_locator(majorLocator)
# ax0.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(0.1)
# minorLocatory = MultipleLocator(0.01)
# ax0.yaxis.set_major_locator(majorLocatory)
# ax0.yaxis.set_minor_locator(minorLocatory)
# ax0.tick_params(which='major',length=5,width=2,labelsize=14)
# ax0.set_xlim(0,20)
# ax0.set_ylim(0,max(dt_DOP853)+.1)
# ax0.legend(fontsize=13,loc="upper right")
# ax0.grid()

# fig3,ax3 = plt.subplots()
# ax3.bar(method, Runtime_avg, color=('orange','green','lightblue','k'))
# ax3.set_ylim(0,max(Runtime_avg)+0.5)
# for i, v in enumerate(Runtime_avg):
#     plt.text(method[i], v + 0.01, str(round(v,3)), fontsize=15, horizontalalignment="center")
# ax3.set_title("Runtime Comparison", fontsize = 20)
# ax3.set_xlabel("Method Used", fontsize = 15)
# ax3.set_ylabel("Average Runtime (s)", fontsize = 15)
# majorLocatory = MultipleLocator(max(Runtime_avg)/10)
# minorLocatory = MultipleLocator(max(Runtime_avg)/50)
# ax3.yaxis.set_major_locator(majorLocatory)
# ax3.yaxis.set_minor_locator(minorLocatory)





# # # Trajectory angle
# # fig1, ax1 = plt.subplots()
# # ax1.plot(rk45_rocket_time, rk45_traj_angle , label='RK45', linewidth=3,color='orange')
# # ax1.plot(rk23_rocket_time, rk23_traj_angle , label='RK23', linewidth=3,color='green')
# # ax1.plot(DOP853_rocket_time, DOP853_traj_angle , label='DOP853', linewidth=3,color='lightblue')
# # ax1.plot(LSODA_rocket_time, LSODA_traj_angle , label='LSODA', linewidth=3,color='k')
# # ax1.set_xlabel('Time (Seconds)', fontsize=15)
# # ax1.set_ylabel("Angle (deg)", fontsize=15)
# # ax1.set_title('Rocket Angle', fontsize=20)
# # ax1.axvline(x=avg_departure_time, color='red', linestyle='--')
# # ax1.axvline(x=burn_time, color='red', linestyle='--')
# # ax1.axvline(x=avg_apogee_time, color='red', linestyle='--')
# # ax1.text(burn_time+.25, 170, 'Burn Out', color='red', fontsize=14, rotation=0)
# # ax1.text(avg_apogee_time+.25, 170, 'Apogee', color='red', fontsize=14, rotation=0)
# # ax1.text(avg_departure_time+.1, 170, 'Avg Departure', color='red', fontsize=12, rotation=0)
# # majorLocator = MultipleLocator(5)
# # minorLocator = MultipleLocator(1)
# # ax1.xaxis.set_major_locator(majorLocator)
# # ax1.xaxis.set_minor_locator(minorLocator)
# # majorLocatory = MultipleLocator(10)
# # minorLocatory = MultipleLocator(1)
# # ax1.yaxis.set_major_locator(majorLocatory)
# # ax1.yaxis.set_minor_locator(minorLocatory)
# # ax1.tick_params(which='major',length=5,width=2,labelsize=14)
# # ax1.set_xlim(0,20)
# # ax1.set_ylim(0,200)
# # ax1.legend(fontsize=11,loc="upper left")
# # ax1.grid()





# # # Altitude
# # fig, ax = plt.subplots()
# # ax.plot(rk45_rocket_time, rk45_rocket_alt, label='rk45', linewidth=3,color='orange')
# # ax.plot(rk23_rocket_time, rk23_rocket_alt, label='rk23', linewidth=3,color='green')
# # ax.plot(DOP853_rocket_time, DOP853_rocket_alt, label='DOP853', linewidth=3,color='lightblue')
# # ax.plot(LSODA_rocket_time, LSODA_rocket_alt, label='LSODA', linewidth=3,color='k')
# # #ax.plot(rocket_time,h*3.28084,label='Analytical Model (No Drag)', linewidth=3,color='b', linestyle='--')
# # ax.set_title('Rocket Trajectory', fontsize=20)
# # ax.axvline(x=avg_departure_time, color='red', linestyle='--')
# # ax.axvline(x=burn_time, color='red', linestyle='--')
# # ax.axvline(x=avg_apogee_time, color='red', linestyle='--')
# # ax.text(burn_time+.25, 3350*3.28084, 'Burn Out', color='red', fontsize=14, rotation=0)
# # ax.text(avg_apogee_time+.25, 350*3.28084, 'Apogee', color='red', fontsize=14, rotation=0)
# # ax.text(avg_departure_time+.25, 350*3.28084, 'Avg Departure', color='red', fontsize=12, rotation=0)
# # majorLocator = MultipleLocator(5)
# # minorLocator = MultipleLocator(1)
# # ax.xaxis.set_major_locator(majorLocator)
# # ax.xaxis.set_minor_locator(minorLocator)
# # majorLocatory = MultipleLocator(1000)
# # minorLocatory = MultipleLocator(100)
# # ax.yaxis.set_major_locator(majorLocatory)
# # ax.yaxis.set_minor_locator(minorLocatory)
# # ax.tick_params(which='major',length=5,width=2,labelsize=14)
# # ax.grid()
# # ax.set_xlabel('Time (s)', fontsize=15)
# # ax.set_ylabel('Altitude (ft)', fontsize=15)
# # ax.set_xlim(0,20)
# # ax.set_ylim(0,3550*3.28084)
# # ax.legend(fontsize=11,loc="upper left")



# # # print(f"Rail Departure Velocity: {rail_departure_velocity*3.28084} ft/s")
# # # print(f"Time to Apogee: {apogee_time} s")
# # # print(f"Apogee: {apogee*3.28084} ft")

# plt.show()



# # Time discretization density
# fig0, ax0 = plt.subplots()
# ax0.plot(rocket_time[0:-1], dt , label='dt', linewidth=3,color='orange')
# ax0.set_xlabel('Time (Seconds)', fontsize=15)
# ax0.set_ylabel("Timestep size (Seconds)", fontsize=15)
# ax0.set_title('Timestep Size', fontsize=20)
# ax0.axvline(x=burn_time, color='red', linestyle='--')
# ax0.axvline(x=apogee_time, color='red', linestyle='--')
# ax0.axvline(x=rail_departure_time, color='red', linestyle='--')
# ax0.text(burn_time+.1, 0.4, 'Burn Out', color='red', fontsize=12, rotation=0)
# ax0.text(apogee_time+.1, 0.4, 'Apogee', color='red', fontsize=12, rotation=0)
# ax0.text(rail_departure_time+.1, 0.4, 'Departure', color='red', fontsize=12, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax0.xaxis.set_major_locator(majorLocator)
# ax0.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(0.1)
# minorLocatory = MultipleLocator(0.01)
# ax0.yaxis.set_major_locator(majorLocatory)
# ax0.yaxis.set_minor_locator(minorLocatory)
# ax0.tick_params(which='major',length=5,width=2,labelsize=14)
# ax0.set_xlim(0,20)
# ax0.set_ylim(0,max(dt)+.1)
# ax0.legend(fontsize=13,loc="upper right")
# ax0.grid()


# # Trajectory angle
# fig1, ax1 = plt.subplots()
# ax1.plot(rocket_time, traj_angle*57.2958 , label='Angle', linewidth=3,color='g')
# ax1.set_xlabel('Time (Seconds)', fontsize=15)
# ax1.set_ylabel("Angle (deg)", fontsize=15)
# ax1.set_title('Rocket Angle', fontsize=20)
# ax1.axvline(x=burn_time, color='red', linestyle='--')
# ax1.axvline(x=apogee_time, color='red', linestyle='--')
# ax1.text(burn_time+.25, 170, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax1.text(apogee_time+.25, 170, 'Apogee', color='red', fontsize=14, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax1.xaxis.set_major_locator(majorLocator)
# ax1.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(10)
# minorLocatory = MultipleLocator(1)
# ax1.yaxis.set_major_locator(majorLocatory)
# ax1.yaxis.set_minor_locator(minorLocatory)
# ax1.tick_params(which='major',length=5,width=2,labelsize=14)
# ax1.set_xlim(0,20)
# ax1.set_ylim(0,200)
# ax1.legend(fontsize=11,loc="upper left")
# ax1.grid()



# # Altitude
# fig, ax = plt.subplots()
# ax.plot(rocket_time, rocket_alt*3.28084, label='Numerical Model (Drag)', linewidth=3,color='k')
# #ax.plot(rocket_time,h*3.28084,label='Analytical Model (No Drag)', linewidth=3,color='b', linestyle='--')
# ax.set_title('Rocket Trajectory', fontsize=20)
# ax.axvline(x=burn_time, color='red', linestyle='--')
# ax.axvline(x=apogee_time, color='red', linestyle='--')
# ax.text(burn_time+.25, 3350*3.28084, 'Burn Out', color='red', fontsize=14, rotation=0)
# ax.text(apogee_time+.25, 3350*3.28084, 'Apogee', color='red', fontsize=14, rotation=0)
# majorLocator = MultipleLocator(5)
# minorLocator = MultipleLocator(1)
# ax.xaxis.set_major_locator(majorLocator)
# ax.xaxis.set_minor_locator(minorLocator)
# majorLocatory = MultipleLocator(1000)
# minorLocatory = MultipleLocator(100)
# ax.yaxis.set_major_locator(majorLocatory)
# ax.yaxis.set_minor_locator(minorLocatory)
# ax.tick_params(which='major',length=5,width=2,labelsize=14)
# ax.grid()
# ax.set_xlabel('Time (s)', fontsize=15)
# ax.set_ylabel('Altitude (ft)', fontsize=15)
# ax.set_xlim(0,20)
# ax.set_ylim(0,3550*3.28084)
# ax.legend(fontsize=11,loc="upper left")

# print(f"Rail Departure Velocity: {rail_departure_velocity*3.28084} ft/s")
# print(f"Time to Apogee: {apogee_time} s")
# print(f"Apogee: {apogee*3.28084} ft")
# plt.show()




# # Solving the analytical case (No drag)
# initial_height = rocket_alt[4]
# initial_velocity = rocket_vel[4] * np.cos(traj_angle[4])
# h = initial_height + initial_velocity*(rocket_time-burn_time) - 0.5 * 9.81 * (rocket_time-burn_time)**2
# h[0:4] = 0.5*(Weq*m_dot/Rocket.mass-9.81)*rocket_time[0:4]**2
