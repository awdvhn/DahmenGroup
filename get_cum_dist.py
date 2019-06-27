#Alan Long 6/10/16
#Last edited: Alan Long 5/23/19

#This code takes data and returns it's complementary cumulative distribution
#function (CCDF). It accepts a list data and returns two lists histx and
#histy, the x and y values for the ccdf.

import numpy as np

def getcumdist(data):
    #Python sorting works with the > operator so we need to remove all nans.
    if np.isnan(data).sum()>0:
        data=data[~np.isnan(data)]
    
    #In order to remove negative values we add in a zero then sort the list. We
    #then remove the zero and all entries before it, i.e. negative numbers.
    #Similarly with infs
    data.extend([0,float('inf')])
    data.sort()
    data[data.index(float('inf')):]=[]
    data.reverse()
    data[data.index(0):]=[]
    data.reverse()
    histx=np.array(data)
    histy=(np.arange(1,0,-1.0/len(histx)))#Tyler just changed 1 to 1.0 
    
    #Gabe found the bug that the same repeated value will return multiple
    #probabilities. This loop fixes that
    for i in range(1,np.size(histx)):
        if histx[i]==histx[i-1]:
            histy[i]=histy[i-1]
    
    histx.tolist()
    histy.tolist()
    return [histx,histy]

