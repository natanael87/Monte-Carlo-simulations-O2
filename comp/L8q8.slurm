#!/bin/bash
#SBATCH --job-name=L8q8 			# Nombre del trabajo
#SBATCH --output=L8q8.log   		# Archivo de registro de salida
#SBATCH --error=L8q8.err    		# Archivo de registro de errores
#SBATCH --partition=QuantPhysMC   	# Nombre de la partición o cola de trabajos
#SBATCH --nodes=1     				# Número de nodos a utilizar (puedes cambiarlo)
#SBATCH --ntasks-per-node=1   		# Número de tareas por nodo (1 para ejecución serial)
#SBATCH --cpus-per-task=4 		# Número de CPUs por tarea (puedes cambiarlo)
#SBATCH --mem=4G      			# Memoria RAM necesaria (puedes cambiarlo)
 
cd /home/elias/Git/3dXYO2e/progm/
# Comando para ejecutar tu programa
./L8q8rm.exe
