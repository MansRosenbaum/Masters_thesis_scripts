"""Script to automize the running of the MolDStruct simulation on temporal folding ubiquitin. 
The script will generate the necessary input files, run the simulation, and save the output files. 
The script will also log the progress and any errors that occur during the simulation."""

import os
import numpy as np
import math
import subprocess
import multiprocessing
import MDAnalysis as md
import h5py
import re


########## Convenient functions ##########

def log_message(message):
    print(message, '\n')

def change_timesteps(file_path, timestep):
    if timestep =='':
        timestep ='10000'
    with open(file_path, "r") as f:
        read_data = f.read()
    read_data = read_data.replace("XXX", timestep)
    with open(file_path, "w") as f:
        f.write(read_data)


# def remove_waters(input_gro, output_gro):
#     with open(input_gro, 'r') as infile, open(output_gro, 'w') as outfile:
#         for line in infile:
#             # gro residue name is in columns 18-20 (1-based indexing, or 17-20 in 0-based)
#             if line.startswith(('ATOM', 'HETATM')) and line[17:20] == 'HOH':
#                 continue  # Skip water molecules
#             outfile.write(line)
    
#     print(f"Water molecules removed. Cleaned gro saved to {output_gro}")

def change_photon_energy(file_path, photon_energy):
    with open(file_path, "r") as f:
        read_data = f.read()
    read_data = read_data.replace('600; Photon energy [eV]', f'{photon_energy}; Photon energy [eV]')
    with open(file_path, "w") as f:
        f.write(read_data)

def change_nrphotons(file_path, nrphotons):
    with open(file_path, "r") as f:
        read_data = f.read()
    read_data = read_data.replace('1e11; Total number of photons in the pulse', f'{nrphotons}; Total number of photons in the pulse')
    with open(file_path, "w") as f:
        f.write(read_data)

def FWHM(file_path, time):
    with open(file_path, "r") as f:
        read_data = f.read()
    peak = read_data.replace('0.005; Peak of the gaussian [ps]', f'{float(time)*10**(-3)}; Peak of the gaussian [ps]')
    with open(file_path, "w") as f:
        f.write(peak)
    with open(file_path, "r") as f:
        read_data = f.read()
    sigma = read_data.replace('0.00125; Width of the peak (sigma values of gaussian) [ps]', f'{round(float(time)*(10**(-3))/(2*math.sqrt(2*math.log(2))), 6)}; Width of the peak (sigma values of gaussian) [ps]')
    with open(file_path, "w") as f:
        f.write(sigma)

def extract_energy(file_path_kin, file_path_tot):
    with open(file_path_kin, 'r') as f:
        read_data = f.readlines()
        kin_en = read_data[-1].split('  ')[-1].split('.')[0]
    with open(file_path_tot, 'r') as f:
        read_data = f.readlines()
        tot_en = read_data[-1].split('  ')[-1].split('.')[0] 
    return print(f'The energy in the system at timestep {nr_timesteps} is approximately {round((int(kin_en)/int(tot_en))*100, 3)}% kinetic energy.\n')


############################ Simulation function for multiprocessing tool ########################


def simulation(sim_number):
    try:
        os.makedirs(f'sim{sim_number}')
        os.system(f'cp -r {cwd}/Atomic_data {newcwd}/sim{sim_number}')
        os.chdir(f'{newcwd}/sim{sim_number}')

        #directory for simulation output
        os.makedirs('simulation_output', exist_ok=True)
        os.makedirs('additional_data', exist_ok=True)
    
        # Run grompp command to prepare the simulation input
        grompp_cmd = [grompp, '-f', f'{newcwd}/exp.mdp', '-c', f'{newcwd}/{gro_name}.gro', '-p', f'{newcwd}/ubi.top', '-o', 'explode.tpr']
        log_message(f"Running grompp command: {' '.join(grompp_cmd)}")
        try:
            result = subprocess.run(grompp_cmd, check=True, capture_output=True, text=True)
            log_message(f"grompp command output: {result.stdout}")
            log_message("grompp command completed successfully")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in grompp command: {e}")
            log_message(f"Error output: {e.stderr}")


        # Run mdrun command to execute the simulation
        mdrun_cmd = [mdrun, '-s', 'explode', '-o', f'{gro_name}_exp.trr', '-x', 'exp.xtc', '-c', f'{gro_name}_exp.gro', '-v', '-nt', '1', '-ionize']
        log_message(f"Running mdrun command: {' '.join(mdrun_cmd)}")
        try:
            result = subprocess.run(mdrun_cmd, check=True, capture_output=True, text=True)
            log_message("mdrun command completed successfully")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in mdrun command: {e}")
            log_message(f"Mdrun Error output: {e.stderr}")    


        # Extract kinetic and total energy from last time step of each simulation
        kin_en_cmd = [energy, '-f', 'ener.edr', '-o', 'kin.xvg']
        try:
            subprocess.run(kin_en_cmd, input='10 \n', check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            log_message(f"kinetic energy Error output: {e.stderr}")    

        tot_en_cmd = [energy, '-f', 'ener.edr', '-o', 'tot.xvg']
        log_message(f"Running g_energy command: {' '.join(tot_en_cmd)}")
        try:
            subprocess.run(tot_en_cmd, input='11 \n', check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            log_message(f"Error in total g_energy command: {e}")
            log_message(f"total energy Error output: {e.stderr}")  

        log_message('------------ENERGY---------------')
        extract_energy(f'{newcwd}/sim{sim_number}/kin.xvg', f'{newcwd}/sim{sim_number}/tot.xvg')
        log_message('---------------------------------')

        # Move unecessary files and folders
        log_message(f'Moving configurations.bin, explode.tpr, mdout.mdp, simulation_output, charges.bin, md.log and state.cpt to extra-folder \n removing Atomic_data')
        os.system(f'mv configurations.bin {newcwd}/sim{sim_number}/additional_data') 
        os.system(f'mv explode.tpr {newcwd}/sim{sim_number}/additional_data')
        os.system(f'mv mdout.mdp {newcwd}/sim{sim_number}/additional_data')
        os.system(f'mv simulation_output {newcwd}/sim{sim_number}/additional_data')
        os.system(f'mv charges.bin {newcwd}/sim{sim_number}/additional_data')
        os.system(f'mv md.log {newcwd}/sim{sim_number}/additional_data')
        os.system(f'mv state.cpt {newcwd}/sim{sim_number}/additional_data')
        os.system(f'rm -rf {newcwd}/sim{sim_number}/Atomic_data')

        #new paths fpr trajectory and structure files for ion mapping and analysis
        trr_path = f'{newcwd}/sim{sim_number}/{gro_name}_exp.trr'
        gro_path = f'{newcwd}/{gro_name}.gro'
        h5file = (f'{newcwd}/sim{sim_number}/data.h5')

        # Analyze the results and collect data
        try:
            universe = md.Universe(gro_path, trr_path)  # Load the trajectory data using MDAnalysis
            ag = universe.atoms.select_atoms("all")  # Select all atoms in the trajectory
            idx = ag.indices  # Get indices of all atoms
            log_message(f"Trajectory loaded from {trr_path}")
        except Exception as e:
            log_message(f"Error loading trajectory: {e}")


        # Access the initial frame of the trajectory
        universe.trajectory[0]  # Move to the first frame of the trajectory
        vel_i = ag.velocities.copy()  # Copy initial velocities
        pos_i = ag.positions.copy()  # Copy initial positions


        # Access the final frame of the trajectory
        universe.trajectory[-1]  # Move to the last frame of the trajectory
        vel_f = ag.velocities.copy()  # Copy final velocities
        pos_f = ag.positions.copy()  # Copy final positions


        # Normalize velocity and displacement vectors
        vel_data = [(x / np.linalg.norm(x)) if np.linalg.norm(x) != 0 else x for x in vel_f]  # Normalize velocities
        pos_data = [(x / np.linalg.norm(x)) if np.linalg.norm(x) != 0 else x for x in (pos_f - pos_i)]  # Normalize displacements
        pos_data = np.array(pos_data)

        # Save data to the HDF5 file
        with h5py.File(h5file, 'a') as file:
            group_path = dirname
            group = file.require_group(group_path)  # Create or get the group for this simulation

            # Save other simulation data
            group.create_dataset("unit_velocity", data=vel_data)
            group.create_dataset("unit_displacement", data=pos_data)
            group.create_dataset("initial_position", data=pos_i)
            group.create_dataset("final_position", data=pos_f)
            group.create_dataset("initial_velocity", data=vel_i)
            group.create_dataset("final_velocity", data=vel_f)
            log_message(f"Simulation data saved to {group_path}")
        
        print(f'##################### Simulation {sim_number} is done. #######################\n')
        os.chdir(newcwd)
    except Exception as e:
        log_message(f"Error in simulation function: {e}")
        os.chdir(newcwd)

    return pos_data

########################################################################################
if __name__ == "__main__":

    #inputs (For convenience when running the script multiple times i have removed the 'input()' commands for each parameter)
    dirname ='run8_ALL_FOLDS_1UBI'
    nr_timesteps = '100000' #how many timesteps you want to run
    flash_xfel = 'xfel'
    if flash_xfel == 'flash':
        photon_energy = '300'
        peak_and_sigma = '30'
        nr_photons = '1e14'
    elif flash_xfel == 'xfel' or flash_xfel == 'XFEL':
        photon_energy = '600'
        nr_photons = '1e12'
    else:
        photon_energy = input("Enter the photon energy: ") #what photon energy you want to use
        peak_and_sigma = input('Enter what FWHM value you want: ')
    nrsims = '100'
    nrcores = '10'

    # Define the paths and create necessary directories
    cwd = '/home/mans/moldstruct_python'
    os.chdir(cwd)

    #new directory for this simulation
    os.makedirs(f'{dirname}', exist_ok=True)
    os.chdir(f'{cwd}/{dirname}')

    prot_list = os.listdir("/home/mans/E_3e4/gro_structures/run8") # ändra denna till mapp med samma fold!
    prot_list.sort()
    #prot_list = prot_list[:-2] #denna med beroende på hur de är numrerade
    for protein in prot_list:
        try:

            os.makedirs(f'1ubi_{protein[12:15]}')
            newcwd = f'{cwd}/{dirname}/1ubi_{protein[12:15]}'
            os.chdir(newcwd)

            #copy mdp file
            #os.system(f'cp -r {cwd}/moldstructinput/Atomic_data .')
            os.system(f'cp {cwd}/moldstructinput/exp.mdp .')

            #get .gro file
            gro_name = protein[:-4]


            # File paths and executable locations
            moldstruct_path = "/home/spidocstester/MolDStruct/bin"  # Path to GROMACS executables
            pdb2gmx = f"{moldstruct_path}/pdb2gmx"  # Path to pdb2gmx executable
            grompp = f"{moldstruct_path}/grompp"  # Path to grompp executable
            mdrun = f"{moldstruct_path}/mdrun"  # Path to mdrun executable
            energy = f"{moldstruct_path}/g_energy"  # Path to g_energy executable

            # change timesteps, nr of photons and photon energy
            change_timesteps('exp.mdp', nr_timesteps)
            change_nrphotons('exp.mdp', nr_photons)
            change_photon_energy('exp.mdp', photon_energy)
            if flash_xfel == 'flash':
                FWHM('exp.mdp', peak_and_sigma)

            #Place the gros in the correct directory
            os.system(f'cp /home/mans/E_3e4/gro_structures/run8/{gro_name}.gro {newcwd}')
            os.system(f'cp /home/mans/E_3e4/run8/ubi.top {newcwd}')
            
            # Execute the parallell processing
            param = [k for k in range(int(nrsims))]
            allsim_posdata = []
            # Use multiprocessing to run simulations with limited cores
            log_message(f"Starting multiprocessing pool with {nrcores} cores")
            with multiprocessing.Pool(int(nrcores)) as pool:
                p = pool.map(simulation, param)    
                for i in p:
                    allsim_posdata.append(i)

            #remove additional files not needed anymore
            log_message('Removing exp.mdp, sample.xyz and ubi.top')
            os.system(f'rm {newcwd}/exp.mdp')
            os.system(f'rm {newcwd}/sample.xyz')
            os.system(f'rm {newcwd}/ubi.top')
        
        except Exception as e:
            log_message(f"Error in main: {e}")

        print(f'We are at protein {protein}')
        os.chdir(f'{cwd}/{dirname}')

    log_message('Done.')