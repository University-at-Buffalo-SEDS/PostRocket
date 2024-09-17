# Extra testing script from the TNAM_v2func file, collected here to avoid running during imports


from cycler import cycler
import os
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from PostRocket.Propulsion.TNAM_v2func import TNAM

flag1 = datetime.now().timestamp()
Do = .098        # Grain outer diameter [m]
Dpi = .01905     # Initial port diameter [m]
L = 0.6096       # Length of fuel grain [m]
Do = .098        # Grain outer diameter [m]
Dpi = .01905     # Initial port diameter [m]
Ox_Vol = 0.02    # Ox tank volume [m^3]

TNAM_ans = TNAM(Do,Dpi,L,Ox_Vol)
flag2 = datetime.now().timestamp()
# Find the vector that describes the mass loss from the oxidizer tank and the fuel    


############################################
        ## Print data to CSV file ##
############################################

data = [ 
    [str(TNAM_ans[0])],
    TNAM_ans[1],
    TNAM_ans[2],
    TNAM_ans[3],
    TNAM_ans[4],
    TNAM_ans[5],
    TNAM_ans[6],
    TNAM_ans[9]
]

data1 = [
    TNAM_ans[7],
    TNAM_ans[8]
]
# ['Pressure [Pa]'],TNAM_ans[6],
# data = data.tolist()
# iterdata = iter(data)

# File path for the csv file
csv_file_ans = 'TNAM_data.csv'
csv_file_2 = 'TNAM_mdotox'

# Open the file in write mode
with open(csv_file_ans,'w', newline='') as file:
    # Create a csv.writer object
    writer = csv.writer(file)
    # Write data to the CSV file
    writer.writerows(data)

with open(csv_file_2,'w', newline='') as file:
    # Create a csv.writer object
    writer = csv.writer(file)
    # Write data to the CSV file
    writer.writerows(data1)

flag3 = datetime.now().timestamp()
print("TNAM run time",flag2-flag1,"[sec]")
print("Run time",flag3-flag2,"[sec]")


# Jokes :)
My_joke = pyjokes.get_joke(language='en', category="all")
print(My_joke)