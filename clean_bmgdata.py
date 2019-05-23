#Alan Long 6/15/16
#Last edited: Alan Long 5/16/18

#This code takes in an array and removes all non-positive, infinite, and not a
#numbers. It accepts an array and outputs an array of the same length.





import numpy as np

def clean_bmgdata(bmgdata):
    if np.isnan(bmgdata).sum()>0:
        bmgdata=bmgdata[~np.isnan(bmgdata)]
    if np.isinf(bmgdata).sum()>0:
        bmgdata=bmgdata[~np.isinf(bmgdata)]
    if (bmgdata<=0).sum()>0:
        bmgdata=bmgdata[bmgdata>0]        


    return bmgdata


