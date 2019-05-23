import numpy as np
from matplotlib import pyplot as plt

def find_nearest(array, value):
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return idx

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
            idxx = find_nearest(avs,bins[i])
            mask = range(idxx-width,idxx+width)

        elif type=='duration':
            idxx = find_nearest(durs,bins[i])
            mask = range(idxx - width, idxx + width)

        shapes_sorted.append(shapes[mask])
        times_sorted.append(times[mask])
        durs_sorted.append(durs[mask])
        avs_sorted.append(avs[mask])

    return [times_sorted, shapes_sorted, durs_sorted, avs_sorted]


def size_avg(shapes,times,avs,durs):
    shapes_final = []
    err_final = []
    times_final =[]
    for i in range(len(shapes)): # for each bin
        span = len(shapes[i])
        lenind = np.argmax(durs[i])
        length = np.size(times[i][lenind])
        avg_shape=np.zeros(length)
        avg_err =np.zeros(length)
        sort_shapes=np.zeros((np.size(shapes[i]),length))

        for k in range(np.size(shapes[i])):

            sort_shapes[k][0:np.size(shapes[i][k])] = shapes[i][k]
        #     plt.plot(times[i][k]-times[i][k][0], shapes[i][k])
        # plt.show()
        for k in range(length):
            avg_shape[k] = sum(sort_shapes[:,k])/span
            avg_err[k] = np.std(sort_shapes[:,k])/np.sqrt(span)
        shapes_final.append(avg_shape)
        err_final.append(avg_err)
        times_final.append(times[i][lenind])
    return [times_final,shapes_final,err_final]

def duration_avg(shapes,times,avs,durs):
    shapes_final = []
    err_final = []
    times_final =[]

    for i in range(len(shapes)): # for each bin
        length = np.size(shapes[i][0])

        for k in range(np.size(shapes[i])):
            [times[i][k],shapes[i][k]]=resize(shapes[i][k],times[i][k],length)

        avg_shape=np.zeros(length)
        avg_err =np.zeros(length)
        sort_shapes=np.zeros((np.size(shapes[i]),length))
        span= np.size(shapes[i])

        for k in range(np.size(shapes[i])):
            sort_shapes[k] = shapes[i][k]

        for k in range(length):
            avg_shape[k] = np.true_divide(np.sum(sort_shapes[:,k]),span)
            avg_err[k] = np.true_divide(np.std(sort_shapes[:,k]),np.sqrt(span))

        shapes_final.append(avg_shape)

        err_final.append(avg_err)
        times_final.append(times[i][0])


    return(times_final,shapes_final,err_final)



def resize(vector,time,length):
    time=np.asarray(time)
    vector=np.asarray(vector)
    time = time-time[0]
    time = time.astype(float)
    time = np.true_divide(time,time[-1])
    new = np.zeros(length)
    points = np.linspace(0,1,num=length)
    width2 = 1.0/length

    for i in range(length-1):
        if i == 0:
            continue

        mask = (time>=(points[i]-(width2/2)))&(time<=(points[i]+(width2/2)))
        new[i]=np.mean(vector[mask])

        if np.isnan(new[i]):
            new[i]=new[i-1]
    return[points,new]