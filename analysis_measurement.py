"""
Basic Script created by Alessio Franchi to analyze the Central Spin data
21/01/2023
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

START_DIR = "C:/Users/Utente/Desktop/central_spin_model"
os.chdir(START_DIR)

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
        data = np.loadtxt(DIR + '\\' + file, unpack=False)
        # transform data to dataframe
        df = pd.DataFrame(data, columns=obs)
        # define a column containing the decoherence factor
        df["deco"] = 0.5 * (1 - (df["magx"]**2 + df["magy"]**2 + df["magz"]**2))
        # group quantities by their step value -> all quantities with the same value of step are listed together. We finally compute the mean values of all values grouped
        df_grouped = df.groupby("step").mean()
        # plot observables and averaged observables
        #plt.plot(df["time"]/int(df["L"][0]), df["magy"])
        plt.plot(df_grouped["time"]/int(df["L"][0]), df_grouped["magy"], label="L="+str(int(df["L"][0])))
        
plt.legend()
plt.show()