L = 50
tauQ = 1600
Nmeas = 1000
Ndata = 2 * tauQ
path = "/home/elias/Git/3dXYO2e/"
#dpath = path + 'analysis/test_d.txt'
dpath = path + f"Data/quench/L{L}tauQ{tauQ}v_rm.txt"
c_name = f'analysis/Corr/C_L{L}tauQ{tauQ}v_rm.txt'

import numpy as np

data = np.loadtxt(dpath, skiprows=0, max_rows=Nmeas)
vort = data[0,:]
Crho = np.zeros(Ndata)
C_rho = np.zeros(Ndata)

for t in range(Ndata):
    Xmean = vort[:Ndata-t].mean()
    Corr = 0.
    for i in range(Ndata-t):
        Corr = Corr + vort[i+t] * (vort[i] - Xmean)
    Crho[t] = Corr / (Ndata - t)

c_path = path + c_name
with open(c_path, 'w') as file:
    for t in range(Ndata):
        file.write(f'{t}\t{Crho[t]:.16f}\t{C_rho[t]:.16f}\n')
