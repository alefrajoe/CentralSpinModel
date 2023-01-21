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
y = 'magz'
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
        # plot the data
        plt.plot(data[obs.index(x)], data[obs.index(y)])
        plt.xlabel(x)
        plt.ylabel(y)

plt.show()