L = 8

import jackknife as jk

import numpy as np

path = "/home/elias/Git/3dXYO2e/"
data_folder = path + "Data/"
beta_path = path + "betas/betas.txt"
analysis_out = f"analysis/afiles/L{L}m_1.txt"
path_out = path + analysis_out


try:
    betas = np.loadtxt(beta_path)
except FileNotFoundError:
    print("El archivo {beta_path} no se encontró.")
    exit()
    
Nb = betas.size

E_ave = np.zeros(Nb)
J_err = np.zeros(Nb)
N_vort_ave = np.zeros(Nb)
J_err_vort = np.zeros(Nb)

for i in range(0, Nb):
    data_name = f"L{L}b{betas[i]:.6f}.txt"
    data_path = data_folder + data_name

    try:
        data = np.loadtxt(data_path, skiprows=1, max_rows=jk.Nmeas)
    except FileNotFoundError:
        print(f"El archivo {data_path} no se encontró.")
        continue
    
    E = data[:,0]
    E_ave[i] = np.mean(E, axis=0)
    J_err[i] = jk.Jacknife_error(E, E_ave[i])
    
    N_vort = data[:,5]
    N_vort_ave[i] = np.mean(N_vort, axis=0)
    J_err_vort[i] = jk.Jacknife_error(N_vort, N_vort_ave[i])

E_ave = E_ave/L**3
J_err = J_err/L**3
N_vort_ave = N_vort_ave/L**3
J_err_vort = J_err_vort/L**3

with open(path_out, "w") as file:
    file.write("beta\tE\tE_err\tN vort\tN vort err:\n")
    for i in range(Nb):
        file.write(f"{betas[i]:.6f}\t{E_ave[i]:.16f}\t{J_err[i]:.16f}\t"
           f"{N_vort_ave[i]:.16f}\t{J_err_vort[i]:0.16f}\n")