"""
Basic Script created by Alessio Franchi to analyze the Central Spin data
21/01/2023
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

##########################################################
#            Observables dictionary
##########################################################
#########################
# Observable that will be shown
Y = 'deco'

START_DIR = "C:/Users/Utente/Desktop/central_spin_model"
os.chdir(START_DIR)

# define the directory where all data are stored
DIR = 'data_centralspin'

# list all files in this directory
FILES = os.listdir(DIR)
# directory where analyzed files will be stored
NEW_DIR = "plot_central"
path = Path(os.path.join(DIR, NEW_DIR))
# if the directory doesn't exist create it
if(not path.exists()):
    os.mkdir(path)

# for all files in DIR
for file in FILES:
    # if file is a txt file -> doesn't throw an error due to directories or other files inside 'data_cental_spin'
    if(".txt" in file):
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
            df_mean = df.groupby("step").mean()
            # for each column, save the number of averaged measure
            df_len = df.groupby("step").apply(len)
            # compute the std deviation of the mean (it has an additional N^{-0.5} with respect to \sigma)
            df_stddev_mean = df.groupby("step").std() / np.sqrt(df_len.iloc[0])
            # plot observables and averaged observables with band lines
            plt.plot(df_mean["time"]/int(df["L"][0]), df_mean[Y], label="L="+str(int(df["L"][0])))
            plt.fill_between(df_mean["time"]/int(df["L"][0]), df_mean[Y]-df_stddev_mean[Y], df_mean[Y]+df_stddev_mean[Y])
        
            # additional columns
            df_mean["measure"] = df_len
            df_mean["err_magx"] = df_stddev_mean["magx"]
            df_mean["err_magy"] = df_stddev_mean["magy"]
            df_mean["err_magz"] = df_stddev_mean["magz"]
            df_mean["err_deco"] = df_stddev_mean["deco"]
        
            # create the plot file if it does not exist
            path_plot = Path(os.path.join(path, file))
            if(not path_plot.exists()):
                # save the analyzed file into the director
                with open(os.path.join(path, file), 'w') as gfile:
                    # first line
                    gfile.write("#\t")
                    # for all columns in df_mean
                    for col in df_mean.columns:
                        gfile.write(f"{col}\t")
                    gfile.write("\n")
                    np.savetxt(gfile, df_mean.values)
                
plt.legend()
plt.show()