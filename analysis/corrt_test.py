L = 50
tauQ = 8
Nmeas = 10
Ndata = 2 * tauQ
path = "/home/elias/Git/3dXYO2e/"
dpath = path + f"Data/quench/L{L}tauQ{tauQ}v_rm.txt"
c_name = f'analysis/Corr/C_L{L}tauQ{tauQ}v_rm.txt'

import numpy as np

data = np.loadtxt(dpath, skiprows=0, max_rows=Nmeas)
vort = data[:]
Crho = np.zeros(Nmeas)
C_rho = np.zeros(Nmeas)

vort = vort / L**3
sigma2 = np.var(vort)
Xmean = np.mean(vort)
c_path = path + c_name
with open(c_path, 'w') as file:
    for t in range(Nmeas):
        denom = Nmeas - t
        X_mean = np.sum(vort[0:denom-1]) / denom
        for j in range(Nmeas-t-1):
            C_rho[t] = C_rho[t] + vort[j+t] * (vort[j] - X_mean)
        C_rho[t] = C_rho[t] / denom
        for j in range(Nmeas-t):
            Crho[t] = Crho[t] + (vort[j] - Xmean) * (vort[j+t] - Xmean)
        Crho[t] = Crho[t] / (denom * sigma2)
        file.write(f'{t}\t{Crho[t]:.16f}\t{C_rho[t]:.16f}\n')
