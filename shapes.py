#Tyler Saulners 5/24/19
#Last edited Alan Long 9/9/19
#This code returns the rate of change of stress over time for an avalanche.







import numpy as np
from matplotlib import pyplot as plt

#This function finds the nearest point in an array to a given value. It takes in array, an array, and value, a float. 
#It outputs idx, an int, which is the index of the nearest point.
def find_nearest(array, value):
    near = [abs(i-value) for i in array]
    indnear = near.index(min(near))
    return indnear


#This function bins the shapes based on duration of size. It takes lists durs, avs, shapes, times which are the output of get_slips. 
# bins is the centers of the bins to be sorted into, type is a str and either 'size' or 'duration' is what you are binning by, and 
#width is the width of the bins. It outputs lists times_sorted, shapes_sorted, durs_sorted, avs_sorted which are those binned and sorted.
def shape_bins(durs,avs,shapes,times,bins,type,width):
    # first sort the arrays
    shapes=np.asarray(shapes)
    times=np.asarray(times)
    if type=='duration':
        ind=np.argsort(durs)
    else:
        ind = np.argsort(avs)

    avs = avs[ind]
    shapes = shapes[ind]
    times = times[ind]
    durs = durs[ind]

    #make return array
    shapes_sorted = []
    times_sorted = []
    durs_sorted = []
    avs_sorted = []

    for i in range(len(bins)):
        if type=='size':
            idxx = find_nearest(avs,bins[i])#find closest event to bin center
            mask = range(idxx-width,idxx+width)#take all events within ben width

        elif type=='duration':
            idxx = find_nearest(durs,bins[i])#find closest event to bin center
            mask = range(idxx - width, idxx + width)#take all events within ben width

        shapes_sorted.append(shapes[mask])
        times_sorted.append(times[mask])
        durs_sorted.append(durs[mask])
        avs_sorted.append(avs[mask])

    return [times_sorted, shapes_sorted, durs_sorted, avs_sorted]





#This averages over sizes. It takes lists shapes,times,avs,durs the output of get_slips. It outputs lists times_final, 
#shapes_final, and err_final the time, velocity, and error vectors respectively.
def size_avg(shapes,times,avs,durs):
    shapes_final = []
    err_final = []
    times_final =[]
    for i in range(len(shapes)): # for each bin
        span = len(shapes[i])#number of shapes
        lenind = np.argmax(durs[i])#which event is longest
        length = np.size(times[i][lenind])#length of time trace
        avg_shape=np.zeros(length)
        avg_err =np.zeros(length)
        sort_shapes=np.zeros((np.size(shapes[i]),length))

        for k in range(np.size(shapes[i])):

            sort_shapes[k][0:np.size(shapes[i][k])] = shapes[i][k]#collect shapes, padded at end with 0 if no data
        for k in range(length):
            avg_shape[k] = sum(sort_shapes[:,k])/span#average of shapes
            avg_err[k] = np.std(sort_shapes[:,k])/np.sqrt(span)#error
        shapes_final.append(avg_shape)
        err_final.append(avg_err)
        times_final.append(times[i][lenind])
    return [times_final,shapes_final,err_final]



#This averages over durations. It takes lists shapes,times,avs,durs the output of get_slips. It outputs lists times_final, 
#shapes_final, and err_final the time, velocity, and error vectors respectively.
def duration_avg(shapes,times,avs,durs):
    shapes_final = []
    err_final = []
    times_final =[]

    for i in range(len(shapes)): # for each bin
        length = np.size(shapes[i][0])#length of time trace

        for k in range(np.size(shapes[i])):
            [times[i][k],shapes[i][k]]=resize(shapes[i][k],times[i][k],length)#conform to length

        avg_shape=np.zeros(length)
        avg_err =np.zeros(length)
        sort_shapes=np.zeros((np.size(shapes[i]),length))
        span= np.size(shapes[i])#number of shapes

        for k in range(np.size(shapes[i])):
            sort_shapes[k] = shapes[i][k]

        for k in range(length):
            avg_shape[k] = np.true_divide(np.sum(sort_shapes[:,k]),span)#average
            avg_err[k] = np.true_divide(np.std(sort_shapes[:,k]),np.sqrt(span))#error

        shapes_final.append(avg_shape)

        err_final.append(avg_err)
        times_final.append(times[i][0])


    return(times_final,shapes_final,err_final)


#This code takes a vector and resizes it to a desired length. It takes a list vector and an associate list time. length is and int and
#is the desired len of the result. It outputs two lists points, the normalized time vector, and new, the new vector.
def resize(vector,time,length):
    time=np.asarray(time)
    vector=np.asarray(vector)
    time = time-time[0]#so you allways start at 0
    time = time.astype(float)
    time = np.true_divide(time,time[-1])#normalize
    new = np.zeros(length)
    points = np.linspace(0,1,num=length)
    width2 = 1.0/length#step of points

    for i in range(length-1):
        if i == 0:
            continue

        mask = (time>=(points[i]-(width2/2)))&(time<=(points[i]+(width2/2)))#all parts of time within the range of a point
        new[i]=np.mean(vector[mask])

        if np.isnan(new[i]):
            new[i]=new[i-1]#remove nans
    return[points,new]
