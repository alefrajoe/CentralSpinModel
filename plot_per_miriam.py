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
y = 'magx'
###################################################################

# define the directory where all data are stored
DIR = 'data_centralspin'

fig=plt.figure()

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
        g = data[obs.index('g')]
        lam = data[obs.index('lambda')]
        k = data[obs.index('kappa')]
        h = data[obs.index('h')]
        p = data[obs.index('p')]
        tm = data[obs.index('tm')]
        
        if 'magx' in obs:
            magx = data[obs.index('magx')]
            magy = data[obs.index('magy')]
            magz = data[obs.index('magz')]
        if 'magGSx' in obs:
            maggsx = data[obs.index('magGSx')]
            maggsy = data[obs.index('magGSy')] 
            maggsz = data[obs.index('magGSz')] 
        
            mag = np.array([magx, magy, magz])
            maggs = np.array([maggsx, maggsy, maggsz])
            angle = (magx*maggsx+magy*maggsy+magz*maggsz)/(np.sqrt(magx*magx+magy*magy+magz*magz)*np.sqrt(maggsx*maggsx+maggsy*maggsy+maggsz*maggsz))
        deco = 0.5 * (1.0 - (magx*magx+magy*magy+magz*magz))
        # plot the data
        plt.subplot(221)
        plt.plot(data[obs.index(x)], magx, linestyle="-")
        plt.ylabel("<x>", labelpad=-5)
        plt.xlabel("t")
        plt.subplot(222)
        plt.plot(data[obs.index(x)], magy, linestyle="-")
        plt.ylabel("<y>", labelpad=-5)
        plt.xlabel("t")
        plt.subplot(223)
        plt.plot(data[obs.index(x)], magz, linestyle="-")
        plt.ylabel("<z>", labelpad=-5)
        plt.xlabel("t")
        plt.subplot(224)
        plt.plot(data[obs.index(x)], deco, linestyle="-")
        plt.ylabel("D", labelpad=-15)
        plt.xlabel("t")

fig.suptitle("L="+str(int(L[0]))+", g="+str(g[0])+r", $\lambda$="+str(lam[0])+r", $\kappa$="+str(k[0])+", h="+str(h[0])+", p="+str(p[0]), fontsize=12)
plt.show()