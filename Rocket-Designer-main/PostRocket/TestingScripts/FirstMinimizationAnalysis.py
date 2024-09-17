import numpy as np
import pickle
from PostRocket.FlightSim.Trajectory import Trajectory
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from moviepy.editor import ImageSequenceClip
folder = 'MinimizationTestPickles'
results = pickle.load(open(f'{folder}/Trial_One_Results.p','rb'))

rail_length = 10
launch_angle = 4 

def nice_plot(x,y,ax=plt.axes(),title = None):
    ax.plot(x,y)
    ax.grid(True,'major',alpha=0.25)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True,'minor', alpha= 0.125)
    ax.set_title(title)

for i,x in enumerate(results.allvecs):
    rocket = pickle.load(open(f'{folder}/RocketIteration{i}.p','rb'))
    thrust_curve = rocket.combustion_chamber.thrust_curve
    burn_time = rocket.combustion_chamber.burn_time
    thrust_time = np.linspace(0,burn_time,len(thrust_curve))
    trajectory_output = Trajectory(rocket, rail_length, launch_angle, 'LSODA', 1e-4, 1e-4, False)
    alt_curve = trajectory_output[2]
    alt_time = trajectory_output[4]
    fig1 = plt.figure(100,[15,10])
    ax1 = plt.subplot(2,1,1)
    rocket.draw(ax1)
    ax2 = plt.subplot(2,2,3)
    nice_plot(thrust_time,thrust_curve,ax2,'Thrust [N] vs Time')
    ax2.grid('both')
    ax3 = plt.subplot(2,2,4)
    nice_plot(alt_time,alt_curve * 3.28084,ax3,'Altitude [ft] vs Time')
    plt.savefig(f'{folder}/RocketIteration{i}.png')
    plt.close(100)
    print(f'{i}/{len(results.allvecs)}',end='\r')

print('done!')

image_folder=f"{folder}/"
image_files = []
for i,x in enumerate(results.allvecs):
    image_files.append(f'{folder}/RocketIteration{i}.png')

clip = ImageSequenceClip(image_files, fps=5)
clip.write_videofile('MinimizationSequence.mp4')
