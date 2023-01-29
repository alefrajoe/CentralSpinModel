# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 18:20:15 2020

@author: alessio
"""

import numpy as np 
import os
import re, shutil

start = '' 
end = ''

word = 'temp'
magic_word = 'work'

start_dir = "C:/Users/Utente/Desktop/central_spin_model"
#move to the working directory
os.chdir(start_dir)

#read all the directories
dir_files = os.listdir()


#for all files in the directory
for file in dir_files:
    #if the name of the directory contain 'data'
    if 'data' in file:
        #then enter inside
        os.chdir(start_dir + '/' + file)
        #list all the files
        sub_files = os.listdir()
        for sub_file in sub_files:
            #if the name is temporary
            if 'seed' in sub_file:
                #open the file read and append
                with open(sub_file,mode='r') as temp_data:
                    #append to the sub_file without temp.. create string:
                    string = re.search('[0-9]+seed', sub_file)
                    string_append_file = sub_file[:string.start()]  + sub_file[string.end():]
                    
                    #read data 
                    with open(string_append_file, 'a') as append_file:
                        #write the read data from sub_file
                        append_file.write(temp_data.read())
                    
                
                #delete temp
                os.remove(sub_file)
                
                
                #pointer = re.search('.txt', sub_file)
                #os.rename(sub_file, sub_file[:pointer.start()] + magic_word + sub_file[pointer.start():] )
                
