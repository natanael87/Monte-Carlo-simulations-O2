L = 8
tauQ = 32
Nmeas = 1000
Nbins = 20

if Nmeas%Nbins == 0:
    parts = Nmeas // Nbins
else:
    print(f"{Nmeas} is not divisible by {bins}, select another number of bins")

import numpy as np

path = "/home/elias/Git/3dXYO2e/"
beta_path = path + f"betas/betas_{tauQ}.txt"
data_name = f"L{L}tauQ{tauQ}v_rm.txt"
analysis_out = f"analysis/afiles/L{L}tau{tauQ}rm_.txt"

try:
    betas = np.loadtxt(beta_path)
except FileNotFoundError:
    print("El archivo {beta_path} no se encontró.")
    exit()
    
Nb = betas.size

data_folder = path + "Data/"
data_path = data_folder + data_name
try:
    data = np.loadtxt(data_path, max_rows=Nmeas)
except FileNotFoundError:
    print(f"El archivo {data_path} no se encontró.")

#vort = data[:,60]
#vort_a = np.mean(vort, axis=0) / L**3
#print(vort_a)

n_vort_ave = np.zeros(Nb)
J_err = np.zeros(Nb)

for i in range(0, Nb):
    
    n_vort = data[:,i]
    n_vort_ave[i] = np.mean(n_vort, axis=0)
    
    factor = (Nbins - 1) / Nbins
    partial_mean = np.array([])
    for j in np.arange(0,Nmeas,parts):
        new_set = np.concatenate((n_vort[:j], n_vort[j+parts:]))
        partial_mean = np.append(partial_mean, new_set.mean())

    var = (partial_mean - n_vort_ave[i])**2
    var = factor * var.sum()
    J_err[i] = np.sqrt(var)

vort_ave = n_vort_ave/L**3
J_err = J_err/L**3

path_out = path + analysis_out

with open(path_out, "w") as file:
    file.write("\t\tt_M\tbeta\t\t<N vort>\t\t<N vort> err:\n")
    file.write("\t\t1\t\t2\t\t\t\t\t3\t\t\t\t\t4\n")
    for i in range(Nb):
        file.write(f"{i:8d}\t{betas[i]:.6f}\t{vort_ave[i]:.16f}\t{J_err[i]:.16f}\n")