#Alan Long 6/5/18
#Last edited 6/6/18

#DONT USE THIS CODE ON THE TEST DATA, THIS IS FOR WENDY DATA

#This code reads all csv files in a folder and makes them into numpy arrays. It
#takes as an input a string folder which is the folder name and a string system
#which is what kind of os you have, 'PC' or 'Mac'. It outputs a list of numpy arrays
#all_data, each element of the list corresponds to a file.
#IMPORTANT: use double backslashes if you are on Windows, otherwise python
#won't understand it. 


import numpy as np
import csv
import os

def csv_data_reader(folder,system):
    #First we initialize the output
    all_data=list()
    #Then we get a list of the files
    files=os.listdir(folder)
   
    #We now iterate over the list and join the data together.
    for i in files:
        file_data=list();#initialize
        if system=='PC':
            name=folder+'\\'+i#this is the name of each file in the folder
        if system=='Mac':
            name=folder+'/'+i
        with open(name, newline='') as testx:
            read_data=csv.reader(testx)#reads file
            for row in read_data:
                file_data.append([float(i) for i in row])#elements are naturally strs
            all_data.append(np.array(file_data))
    return all_data
