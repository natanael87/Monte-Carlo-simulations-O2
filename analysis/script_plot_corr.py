import matplotlib.pyplot as plt
import numpy as np

# Lee los datos del archivo
datos = np.loadtxt('/home/elias/Git/3dXYO2e/analysis/Corr/C_L50tauQ1024v_rm.txt')

# Extrae las columnas 1 y 2
columna1 = datos[:, 0]
columna2 = datos[:, 1]

# Grafica las columnas
plt.plot(columna1, columna2, 'bs')
plt.xlabel('t')
plt.ylabel('correlation')
plt.title('Autocorrelation time')
plt.legend()
plt.show()

