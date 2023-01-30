"""
Basic Script created by Alessio Franchi to visualize the Central Spin data
21/01/2023
"""
import numpy as np
import os
import matplotlib.pyplot as plt

###################################################################
#                 Plot parameters
x = 'time'
y = ''
###################################################################

# define the directory where all data are stored
DIR = 'data_centralspin\\plot_central'

# list all files in this directory
FILES = os.listdir(DIR)
# for all files in DIR
for file in FILES:
    if ".txt" in file:
        # open the file in reading mode
        with open(DIR + '\\' + file, 'r') as f:
            # first character is "#" -> delete
            # then split the first line to obtain all observables as a list
            obs = f.readline()[1:].split()
            print("Available observables are :", obs)
            # load data
            data = np.loadtxt(DIR + '\\' + file, unpack=True)
            
            # define useful paramters
            L = data[obs.index('L')]
            g = data[obs.index('g')]
            if 'magx' in obs:
                magx = data[obs.index('magx')]
                magy = data[obs.index('magy')]
                magz = data[obs.index('magz')]
            
                #mag = np.array([magx, magy, magz])
                #maggs = np.array([maggsx, maggsy, maggsz])
                #angle = (magx*maggsx+magy*maggsy+magz*maggsz)/(np.sqrt(magx*magx+magy*magy+magz*magz)*np.sqrt(maggsx*maggsx+maggsy*maggsy+maggsz*maggsz))
            #deco = 0.5 * (1.0 - (magx*magx+magy*magy+magz*magz))
            # plot the data
            plt.errorbar(data[obs.index(x)]/L, data[obs.index(y)], linestyle="-", label="L="+str(int(L[0])))
            plt.xlabel(x)
            plt.ylabel(y)

plt.legend()
plt.show()