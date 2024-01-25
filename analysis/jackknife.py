Nmeas = 1000
Nbins = 20

import numpy as np

if Nmeas%Nbins == 0:
    parts = Nmeas // Nbins
else:
    print(f"{Nmeas} is not divisible by {Nbins}, select another number of bins")

def Jacknife_error(Xdata, Xmean, Xbins = Nbins):
    try:
        factor = (Nbins - 1) / Nbins
    except:
        print("Posible division por cero")
        exit()
        
    num = Xdata.size
    partial_mean = np.array([])
    for j in np.arange(0,num,parts):
        new_set = np.concatenate((Xdata[:j], Xdata[j+parts:]))
        partial_mean = np.append(partial_mean, new_set.mean())

    var = (partial_mean - Xmean)**2
    var = factor * var.sum()
    return np.sqrt(var)