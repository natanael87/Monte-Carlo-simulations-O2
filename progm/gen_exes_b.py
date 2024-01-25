import numpy as np
import subprocess

# Define los parámetros
L_lat = 8
tau_Q = [8,32,128]

beta_c = 0.454165
beta_initial = 0.3
beta_final = 0.5
beta_step = 0.01
total_betas = 5
total_parts = 2

folder = '/home/elias/Git/3dXYO2e/'
parts = list(range(1,total_parts + 1))

equilibrium_procces = "CALL m3DXYO2_standard"
quench_procces = "CALL quench_algorithm"
# Choose 'equilibrium_procces' or 'quench_procces'
main_subroutine = quench_procces

def generate_n_betas(initial, final, total):
    # Abre el archivo de salida para escribir los nuevos valores de betas
    beta_file = folder + f'betas/betas.txt'
    betas = np.linspace(initial, final, total)
    np.savetxt(beta_file, betas, fmt="%13.6f")
    return betas

def generate_betas_step(initial, final, step):
    # Abre el archivo de salida para escribir los nuevos valores de betas
    beta_file = folder + f'betas/betas.txt'
    betas = np.arange(initial, final, step)
    np.savetxt(beta_file, betas, fmt="%13.6f")
    return betas

def betas_for_quench(qrate):
    Temp_inv = np.zeros(2 * qrate )
    for i in range(2 * qrate):
        Temp_inv[i] = beta_c / (1 - (i - qrate) / qrate)
    file_name = folder + f'betas/betas_{tQ}.txt' 
    np.savetxt(file_name, Temp_inv, fmt = "%.6f")
    return Temp_inv

# Escribe las lineas del programa fmain0.f90
lineas_main = ["PROGRAM fmain",\
        "USE MC",\
        "IMPLICIT NONE",\
        main_subroutine,\
        "END PROGRAM fmain"]
with open('mainf0.f90', "w") as archivo_salida:
    for linea in lineas_main:
        archivo_salida.write(linea + "\n")

if main_subroutine == equilibrium_procces:
    betas_lat = generate_betas_step(beta_initial, beta_final, beta_step)
    #betas_lat = generate_n_betas(beta_initial, beta_final, total_betas)
    
    b_tot = betas_lat.size
    # Abre el archivo de entrada 1 para lectura
    with open(folder + 'progm/parameters.f90', 'r') as archivo_entrada1:
        lineas1 = archivo_entrada1.readlines()
        
    for i, linea in enumerate(lineas1):
        # Modifica la línea con el nuevo valor de L
        if 'INTEGER,	    PARAMETER :: L =' in linea:
            lineas1[i] = f'\tINTEGER,	    PARAMETER :: L = {L_lat}\n'
        # Modifica la línea con los nuevos valores de beta
        if 'betas(' in linea:
            lineas1[i] = f'\tREAL(8),        PARAMETER :: betas({b_tot}) = [{betas_lat[0]:.6f}D0&\n'
            for j in range(1,b_tot):
                new_line = f'\t\t\t,{betas_lat[j]:.6f}D0&\n'
                lineas1.insert(i+j,new_line)
            lineas1.insert(i+b_tot,'\t\t\t]\n')
            
    # Abre el archivo de salida 1 para escribir el nuevo valor de L y betas
    new_param_file = 'parameters0.f90'
    with open(new_param_file, 'w') as archivo_salida1:
        archivo_salida1.writelines(lineas1)
        
    # Abre el archivo de entrada 2 para lectura
    with open(folder + 'progm/fcompile.sh', 'r') as archivo_entrada2:
        lineas2 = archivo_entrada2.readlines()
        
    for i, linea in enumerate(lineas2):
        # Modifica la líneas que contienen "parameters", "mainf" y ".exe"
        new_param_name = 'parameters0'
        if 'parameters' in linea:
            lineas2[i] = linea.replace('parameters', new_param_name)
        if 'mainf' in linea:
            lineas2[i] = linea.replace('mainf', 'mainf0')
        if '.exe' in linea:
            new_exe = f'L{L_lat}b{beta_initial}_x{b_tot}.exe'
            lineas2[i] = f'gfortran *.o -o {new_exe} -llapack'

    # Abre el archivo de salida 2 para escribir las lineas modificadas
    new_compile_file = 'fcompile0.sh'
    with open(new_compile_file, 'w') as archivo_salida2:
        archivo_salida2.writelines(lineas2)
        
    print(new_exe)
    command1 = 'bash '+ new_compile_file
    subprocess.run(command1, shell=True)

    command2 = 'rm ' + new_compile_file
    command3 = 'rm ' + new_param_file
    subprocess.run(command2, shell=True)
    subprocess.run(command3, shell=True)

if main_subroutine == quench_procces:
    for tQ in tau_Q:
        betas_lat = betas_for_quench(tQ)
        exe_name_tail = f'{tQ}h'
        
        b_tot = betas_lat.size
        print(exe_name_tail, b_tot)

        # Abre el archivo de entrada 1 para lectura
        with open(folder + 'progm/parameters.f90', 'r') as archivo_entrada1:
            lineas1 = archivo_entrada1.readlines()
            
        for i, linea in enumerate(lineas1):
            # Modifica las líneas con los nuevos valores de L, tau_Q y betas
            if 'INTEGER,	    PARAMETER :: L =' in linea:
                lineas1[i] = f'\tINTEGER,	    PARAMETER :: L = {L_lat}\n'
            if 'tauQ' in linea:
                lineas1[i] = f'\tINTEGER,	    PARAMETER :: tauQ = {tQ}\n'
            if 'betas(' in linea:
                lineas1[i] = f'\tREAL(8),        PARAMETER :: betas({b_tot}) = [{betas_lat[0]:.6f}D0&\n'
                for j in range(1,b_tot):
                    new_line = f'\t\t\t,{betas_lat[j]:.6f}D0&\n'
                    lineas1.insert(i+j,new_line)
                lineas1.insert(i+b_tot,'\t\t\t]\n')
                
        # Abre el archivo de salida 1 para escribir los nuevos valores de L, tau_Q y betas
        new_param_file = 'parameters' + str(tQ) + '.f90'
        with open(new_param_file, 'w') as archivo_salida1:
            archivo_salida1.writelines(lineas1)
            
        # Abre el archivo de entrada 2 para lectura
        with open(folder + 'progm/fcompile.sh', 'r') as archivo_entrada2:
            lineas2 = archivo_entrada2.readlines()
            
        for i, linea in enumerate(lineas2):
            # Modifica la líneas que contienen "parameters", "mainf" y ".exe"
            new_param_name = 'parameters' + str(tQ) 
            if 'parameters' in linea:
                lineas2[i] = linea.replace('parameters', new_param_name)
            if 'mainf' in linea:
                lineas2[i] = linea.replace('mainf', 'mainf0')
            if '.exe' in linea:
                new_exe = f'L{L_lat}q{exe_name_tail}.exe'
                lineas2[i] = f'gfortran *.o -o {new_exe} -llapack'

        # Abre el archivo de salida 2 para escribir las lineas modificadas
        new_compile_file = 'fcompile' + str(tQ) + '.sh'
        with open(new_compile_file, 'w') as archivo_salida2:
            archivo_salida2.writelines(lineas2)
            
        print(new_exe)
        # Genera los ejecutables
        command1 = 'bash '+ new_compile_file
        subprocess.run(command1, shell=True)
        
        # Borra los archivos modificados
        command2 = 'rm ' + new_compile_file
        command3 = 'rm ' + new_param_file
        subprocess.run(command2, shell=True)
        subprocess.run(command3, shell=True)
        
        job_name = f'L{L_lat}q{tQ}'
        lines_sh = [f'#!/bin/bash',\
            f'#SBATCH --job-name={job_name} \t\t\t# Nombre del trabajo',\
            f'#SBATCH --output={job_name}.log   \t\t# Archivo de registro de salida',\
            f'#SBATCH --error={job_name}.err    \t\t# Archivo de registro de errores',\
            f'#SBATCH --partition=QuantPhysMC   \t# Nombre de la partición o cola de trabajos',\
            f'#SBATCH --nodes=1     \t\t\t\t# Número de nodos a utilizar (puedes cambiarlo)',\
            f'#SBATCH --ntasks-per-node=1   \t\t# Número de tareas por nodo (1 para ejecución serial)',\
            f'#SBATCH --cpus-per-task=4 \t\t# Número de CPUs por tarea (puedes cambiarlo)',\
            f'#SBATCH --mem=4G      \t\t\t# Memoria RAM necesaria (puedes cambiarlo)',
            ' ',\
            f'cd {folder}progm/',\
            f'# Comando para ejecutar tu programa',\
            f'./{new_exe}']

        with open(f'{folder}comp/{job_name}.slurm', "w") as archivo_salida:
            for linea in lines_sh:
                archivo_salida.write(linea + "\n")