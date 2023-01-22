"""
Basic Script created by Alessio Franchi to visualize the Central Spin data
21/01/2023
"""
import numpy as np
import os
import matplotlib.pyplot as plt


###################################################################
#                 Plot parameters
x = 'kappa'
y = 'magy'
###################################################################

# define the directory where all data are stored
DIR = 'data_centralspin'

# list all files in this directory
FILES = os.listdir(DIR)
# for all files in DIR
for file in FILES:
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
        if 'magx' in obs:
            magx = data[obs.index('magx')]
            magy = data[obs.index('magy')]
            magz = data[obs.index('magz')]
        if 'magGSx' in obs:
            maggsx = data[obs.index('magGSx')] 
            maggsy = data[obs.index('magGSy')] 
            maggsz = data[obs.index('magGSz')] 
        # plot the data
        plt.plot(data[obs.index(x)], (1+(magx*maggsx+magy*maggsy+magz*maggsz))/2.0)
        plt.xlabel(x)
        plt.ylabel(y)

plt.show()