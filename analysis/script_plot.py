L = 50
tauQ = 128

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

path = "/home/elias/Git/3dXYO2e/"
c_name = f'analysis/Corr/L{L}q{tauQ}corr.txt'
c_path = path + c_name

try:
    data1 = np.loadtxt(c_path, skiprows=0)
except FileNotFoundError:
    print("El archivo '*1.txt' no se encontró.")
    exit()
    
time = data1[:,0]
Corr1 = data1[:,1]
Corr2 = data1[:,2]

plt.figure(1)
plt.plot(time, Corr1, 'bs')

plt.figure(2)
plt.plot(time, Corr2, 'g^')

plt.xlabel(r'time')
plt.ylabel(r'$C_{\rho}$')
plt.title(r'Correlation time')

plt.grid(True)
plt.show()