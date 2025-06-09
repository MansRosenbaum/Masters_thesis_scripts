"""Script to automize the running of the MolDStruct simulation. 
The script will generate the necessary input files, run the simulation, and save the output files. 
The script will also log the progress and any errors that occur during the simulation."""

import os
import numpy as np
import itertools as it
import subprocess
import multiprocessing
import MDAnalysis as md 
import h5py
import shutil
import random
import glob
import shutil
import re
import scipy
import matplotlib.pyplot as plt

########## Convenient functions ##########

def log_message(message, log_file="simulation_log.txt"):
    with open(log_file, "a") as f:
        f.write(message + "\n")
    print(message)

def change_timesteps(file_path, timestep):
    with open(file_path, "r") as f:
        read_data = f.read()
    read_data = read_data.replace("XXX", timestep)
    with open(file_path, "w") as f:
        f.write(read_data)


def remove_waters(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            # PDB residue name is in columns 18-20 (1-based indexing, or 17-20 in 0-based)
            if line.startswith(('ATOM', 'HETATM')) and line[17:20] == 'HOH':
                continue  # Skip water molecules
            outfile.write(line)
    
    print(f"Water molecules removed. Cleaned PDB saved to {output_pdb}")

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


######################## ion mapping ############################
def bin_plane(xyz, plane_size, L, num_bins, detector_eff=1.0, axis=1):
    if detector_eff > 1.0:
        raise ValueError("Detector efficiency cannot be over 100%.")

    
    # Create the grid
    grid_range = np.linspace(-plane_size / 2, plane_size / 2, num_bins)
    
    if axis == 0:
        grid = np.array([(-L, x, y) for y in grid_range for x in grid_range], dtype=np.float32)
        B1 = 1
        B2 = 2
    elif axis == 1:
        grid = np.array([(x, -L, y) for y in grid_range for x in grid_range], dtype=np.float32)
        B1 = 0
        B2 = 2
    elif axis == 2:
        grid = np.array([(x, y, -L) for y in grid_range for x in grid_range], dtype=np.float32)
        B1 = 1
        B2 = 0
    else:
        raise ValueError("Invalid axis value. Axis must be 0, 1, or 2.")
    
    # Filter out points that will never hit the plane based on their y-coordinate and the plane's location
    if L > 0:
        xyz = xyz[xyz[:, axis] < 0.0]  # Points above the plane when L > 0 will not hit
    elif L < 0:
        xyz = xyz[xyz[:, axis] > 0.0]  # Points below the plane when L < 0 will not hit

    # Calculate the hit locations on the plane
    k = -L / xyz[:, axis]
    hit_location = k[:, np.newaxis] * xyz

    # Initialize the bins
    bins = np.zeros((num_bins, num_bins), dtype=int)

    N = len(hit_location)
    M = int(detector_eff * N)
    
    keep_idx = np.random.choice(N, M, replace=False)
    hit_location = hit_location[keep_idx]

    # Process each hit location
    for hit in hit_location:
        # Check if the hit is within the plane bounds
        if -plane_size / 2 <= hit[B1] <= plane_size / 2 and -plane_size / 2 <= hit[B2] <= plane_size / 2:
            # Find the closest grid point for the hit location
            distance = scipy.spatial.distance.cdist([hit], grid)
            closest_bin = np.argmin(distance)

            # Convert the linear index to 2D index and increment the corresponding bin
            y_idx, x_idx = divmod(closest_bin, num_bins)
            bins[y_idx, x_idx] += 1

    return bins.astype(np.int16)


####################################################

if __name__ == "__main__":

    #inputs
    dirname = input('Directory name:') #so they get different names
    nr_timesteps = input("Enter the number of timesteps: ") #how many timesteps you want to run
    nr_photons = input("Enter the number of photons: ") #how many photons you want to run
    photon_energy = input("Enter the photon energy: ") #what photon energy you want to use


    try:

        # Define the paths and create necessary directories
        cwd = '/home/mans/moldstruct_python'
        os.chdir(cwd)

        #new directory for this simulation
        os.makedirs(f'mds_{dirname}', exist_ok=True)
        os.chdir(f'{cwd}/mds_{dirname}')

        #directory for simulation output
        os.makedirs('simulation_output', exist_ok=True) 
        os.makedirs('results', exist_ok=True)
    
        #copy mdp file
        #os.system(f'cp -r {cwd}/moldstructinput/Atomic_data .')
        os.system(f'cp {cwd}/moldstructinput/exp.mdp .')

        #get .pdb file
        pdb_name = input("Enter the name of the pdb without .pdb at the end.\n If you want to use a local pdb-file, type 'local' \n If you want to run temporal unfolding on ubiquitin, type 'fold': ") #what pdb file you want to use
        
        if pdb_name == 'fold':
            foldlvl = input("Enter the level of unfoldness, (000-100): ")
            pdb_name = f'E_3e4_run0_0{foldlvl}'
            os.system(f'cp {cwd}/ubifold/{pdb_name}.pdb {cwd}/mds_{dirname}')
            os.system(f'cp {cwd}/ubifold/folded.top {cwd}/mds_{dirname}')
        elif pdb_name == 'local':
            pdb_path = input("Enter the pathway to the pdb file you want to use: ")
            os.system(f'cp {pdb_path} .')
            pdb_name = pdb_path.split('/')[-1].split('.')[0]
        else:
            os.system(f'wget https://files.rcsb.org/download/{pdb_name}.pdb')
            # Remove water molecules from the PDB file
            remove_waters(f'{pdb_name}.pdb', f'{pdb_name}_no_waters.pdb')
            os.system(f'rm {pdb_name}.pdb')
            os.system(f'mv {pdb_name}_no_waters.pdb {pdb_name}.pdb')

        # File paths and executable locations
        moldstruct_path = "/home/spidocstester/MolDStruct/bin"  # Path to GROMACS executables
        pdb2gmx = f"{moldstruct_path}/pdb2gmx"  # Path to pdb2gmx executable
        grompp = f"{moldstruct_path}/grompp"  # Path to grompp executable
        mdrun = f"{moldstruct_path}/mdrun"  # Path to mdrun executable
        energy = f"{moldstruct_path}/g_energy"  # Path to g_energy executable
        

        #change timesteps, nr of photons and photon energy
        change_timesteps('exp.mdp', str(nr_timesteps))
        change_nrphotons('exp.mdp', str(nr_photons))
        change_photon_energy('exp.mdp', str(photon_energy))

        # Generate atomic parameters using the provided script
        gen_atomic_param_cmd = f"python3 /home/mans/moldstruct_python/generate_atomic_parameters.py {photon_energy} {pdb_name}_no_waters.pdb {os.getcwd()} {os.getcwd()}"
        log_message(f"Running atomic data generation command: {gen_atomic_param_cmd}")
        try:
            result = subprocess.run(gen_atomic_param_cmd, shell=True, check=True, capture_output=True, text=True)
            log_message(f"Atomic data generation output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in generate_atomic_data: {e}")
            log_message(f"Error output: {e.stderr}")

        # Run pdb2gmx command to generate topology and coordinate files
        pdb2gmx_cmd = [pdb2gmx, '-f', f'{pdb_name}_no_waters.pdb', '-o', f'{pdb_name}.gro', '-water', 'SPC', '-ff', 'charmm27']
        log_message(f"Running pdb2gmx command: {' '.join(pdb2gmx_cmd)}")
        try:
            result = subprocess.run(pdb2gmx_cmd, check=True, capture_output=True, text=True)
            log_message(f"pdb2gmx command output: {result.stdout}")
            log_message("pdb2gmx command completed successfully")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in pdb2gmx command: {e}")
            log_message(f"Error output: {e.stderr}")

        # Run grompp command to prepare the simulation input
        grompp_cmd = [grompp, '-f', 'exp.mdp', '-c', f'{pdb_name}.gro', '-p', 'folded.top', '-o', 'explode.tpr']
        log_message(f"Running grompp command: {' '.join(grompp_cmd)}")
        try:
            result = subprocess.run(grompp_cmd, check=True, capture_output=True, text=True)
            log_message(f"grompp command output: {result.stdout}")
            log_message("grompp command completed successfully")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in grompp command: {e}")
            log_message(f"Error output: {e.stderr}")


        # Run mdrun command to execute the simulation
        mdrun_cmd = [mdrun, '-s', 'explode', '-o', f'{pdb_name}_exp.trr', '-x', 'exp.xtc', '-c', f'{pdb_name}_exp.pdb', '-v', '-nt', '1', '-ionize']
        log_message(f"Running mdrun command: {' '.join(mdrun_cmd)}")
        try:
            result = subprocess.run(mdrun_cmd, check=True, capture_output=True, text=True)
            log_message(f"grompp command output: {result.stdout}")
            log_message("mdrun command completed successfully")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in mdrun command: {e}")
            log_message(f"Error output: {e.stderr}")


    except Exception as e:
        log_message(f"Error in main: {e}")

    # Load trajectories and construct explosion map
    try:
        os.system(f'mv {pdb_name}_exp.trr {cwd}/mds_{dirname}/results')
        os.system(f'mv {pdb_name}.gro {cwd}/mds_{dirname}/results')
        os.system(f'mv {pdb_name}_exp.pdb {cwd}/mds_{dirname}/results')
        log_message(f"Moving {pdb_name}.trr and {pdb_name}.gro to results folder")

        #new paths fpr trajectory and structure files for ion mapping and analysis
        trr_path = f'{cwd}/mds_{dirname}/results/{pdb_name}_exp.trr'
        gro_path = f'{cwd}/mds_{dirname}/results/{pdb_name}.gro'
        h5file = (f'{cwd}/mds_{dirname}/results/data.h5')

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

        data = pos_data # load data into format (num_sims, num_atoms, 3)
        plane_size = input('Plane size for detector image: ') # size of detctor
        detector_distance = input('Plane distance: ')
        num_bins = input('Number of bins: ') # Resolution of the detector (how many bins that fits the plane size, meaning the number of pixels per bin is on each side plane_size/num_bins)
        # orientation of the detector = 0 (degrees or radians)

        detector_eff = 1.0 # effciency of detector
        axis = input('Detector axis placement: ') # Where to place the detector relative to the refernce frame of the data
        # 0 -> x
        # 1 -> y
        # 2 -> z
        if axis == 0:
            axis_xyz = 'x'
        elif axis == 1:
            axis_xyz = 'y'
        else:
            axis_xyz = 'z'

        img = np.ndarray(shape=(num_bins,num_bins))
        img = bin_plane(data,plane_size,detector_distance,num_bins,detector_eff,axis)

        # Assuming img is already defined and contains the data
        log_message(f"Detector image created with shape: {img.shape}")
        plt.imshow(img)
        plt.colorbar(label='Intensity')
        plt.title(f'Detector Image for {pdb_name}')
        plt.xlabel('Bins')
        plt.ylabel('Bins')
        plt.savefig(f'{cwd}/mds_{dirname}/results/img_{axis_xyz}_{plane_size}_{detector_distance}_{num_bins}.png')
        plt.close()
        log_message(f"Detector image saved to {cwd}/mds_{dirname}/results/{dirname}_{plane_size}_{detector_distance}_{num_bins}.png")


    except Exception as e:
        log_message(f"Error in creating explosion map: {e}")

