#Alan Long 6/15/16
#Last edited: Alan Long 5/16/18

#This code takes in an array and removes all non-positive, infinite, and not a
#numbers. It accepts an array and outputs an array of the same length.


# hi mayisha


import numpy as np

def clean_data(data):
    data = np.array(data)
    if np.isnan(data).sum()>0:
        data=data[~np.isnan(data)]
    if np.isinf(data).sum()>0:
        data=data[~np.isinf(data)]
    if (data<=0).sum()>0:
        data=data[data>0]        
    data.tolist()

    return data


