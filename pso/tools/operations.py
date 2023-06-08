import numpy as np
import re
import os
import subprocess as sp

from ase.io import read, write
from ase.constraints import FixAtoms
from ase import Atoms
from ase.geometry import cell_to_cellpar
from ase.visualize import view

# Calculate the cluster rotation matrix
def cluster_rotation_matrix(angle_x, angle_y, angle_z, clusters):
    """ Function to perform a Geometrical Rotation 
        for an user-provided set of coordinates.

    """

    # Get the centroid of the Atoms object
    center_of_mass = clusters.get_center_of_mass()

    # Define the rotation angle in the x, y, z directions, in radians

    rotation_matrix_x = np.array([[1, 0, 0],
                                  [0, np.cos(angle_x), -np.sin(angle_x)],
                                  [0, np.sin(angle_x), np.cos(angle_x)]])

    rotation_matrix_y = np.array([[np.cos(angle_y), 0, np.sin(angle_y)],
                                  [0, 1, 0],
                                  [-np.sin(angle_y), 0, np.cos(angle_y)]])

    rotation_matrix_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0],
                                  [np.sin(angle_z), np.cos(angle_z), 0],
                                  [0, 0, 1]])

    rotation_matrix = np.dot(rotation_matrix_z,
                             np.dot(rotation_matrix_y,
                                    rotation_matrix_x))

    # Rotate all atoms of the Atoms object around the center of mass and output absolute coordinates
    rotated_positions = np.dot(rotation_matrix,
                               (clusters.positions - center_of_mass).T).T+ center_of_mass

    return rotated_positions


# Define the cluster displacement function

def cluster_dis_xyz(x_direct, y_direct, z_min, angle_x, angle_y, angle_z, clusters, substrate):
    """ Function to rotate and displace an user-provided cluster
    """
    # Read the lattice parameters of the substrate from the POSCAR file
    

    # Define the direct coordinates of the offset
    direct_coords = np.array([x_direct, y_direct, 0])

    # Convert the direct coordinates to cartesian coordinates
    cart_coords = direct_coords @ substrate.cell

    translation_vector_xy = np.array(cart_coords)

    # Translate the rotated coordinates along the x and y axes
    translated_positions = cluster_rotation_matrix(angle_x, angle_y,
                                                   angle_z, clusters) + translation_vector_xy

    # Subtract the smallest z-coordinate from the z-coordinates of all atoms,
    #and translate it to a minimum of 2.5 angstrom

    min_z = np.min(translated_positions[:, 2])
    translated_positions[:, 2] -= min_z
    translation_vector_z = np.array([0, 0, z_min])
    translated_positions += translation_vector_z
    
    return translated_positions

def combine_substrate_cluster(substrate, rot_cluster):
    
    # Combine the two Atoms objects
    combined_atoms = substrate + rot_cluster

    combined_atoms.wrap()
    # Output the combined Atoms object to a new POSCAR file

    return combined_atoms

def structure_creator(cluster, substrate, output_name, struc_param, z_min):
   """ Function to create new structures based on user-defined
        parameters
   """

   x_direct, y_direct, angle_x, angle_y, angle_z=struc_param
   
   # Read the cluster.POSCAR file and convert to Cartesian with selective dynamics style
   rot_clusters = cluster.copy()
   cluster.set_constraint(FixAtoms(mask=[False]*len(cluster)))
   
   new_subs=substrate.copy()
   
   if not new_subs.constraints:
      new_subs.set_constraint(True)
      
   #Update the positions of all atoms in the Atoms object to the rotated coordinates
   #plus the translation vector
   
   scale_alpha=2*np.pi*angle_x
   scale_beta=np.pi*angle_y
   scale_gamma=2*np.pi*angle_z
   rot_clusters.positions = cluster_dis_xyz(x_direct, y_direct, z_min, scale_alpha,
                                            scale_beta, scale_gamma, cluster, substrate)

   final=combine_substrate_cluster(new_subs, rot_clusters)
    
   write(str(output_name), final, vasp5=True, direct=False,
         sort=False, long_format=True, symbol_count=None)

   return x_direct, y_direct, scale_alpha, scale_beta, scale_gamma

def read_energy(filename):
   """ Function to read OUTCAR and extract the last total free energy (ION-ELECTRON)
       value.
   """
   e0_pattern = re.compile(r"energy\(sigma->0\)\s*=\s+([\d\-\.]+)")

   textfile = open(filename, 'r')
   filetext = textfile.read()
   textfile.close()

   matches = re.findall(e0_pattern, filetext)

   return matches


def random_number_generator(population_size, dim, seed=42):
    """ Function to generate random numbers from
        a uniform distribution based on the size of population.

    population_size: int
        Number of random number generated. Ideally it will be
        dim*len(population_size). 

    dim: int
        Dimension of the vector. For velocities and positions is 5,
        for pso_updating is 2. 

    seed: A seed to initialize the BitGenerator. If None, then fresh, 
          unpredictable entropy will be pulled from the OS.
     

    """
    pop_size=int(dim)*population_size
    s = np.random.default_rng(seed).uniform(0, 1, int(pop_size))
    
    return np.array_split(s, population_size)

def global_random_numbers(num_gen, pop_per_gen, seed, output_name):
    """ Function to compute random numbers for 
        structural parameters.
    """
    total=np.asarray(random_number_generator(int(num_gen)*int(pop_per_gen), 5, seed=int(seed)))
    np.save(str(output_name), total)
    data=np.load(str(output_name)+".npy")
   
    return np.array_split(data, num_gen)

def pso_random_numbers(num_gen, pop_per_gen, seed, output_name):
    """ Function to compute random numbers for 
        structural parameters.
    """
    total=np.asarray(random_number_generator(int(num_gen)*int(pop_per_gen), 2, seed=int(seed)))
    np.save(str(output_name), total)
    data=np.load(str(output_name)+".npy")
   
    return np.array_split(data, num_gen)

def program_executor(vasp_exe, num_cores, dir_cal,
                     out_file="output.log"):
    """Function to submit a vasp job using subprocess
    """
    
    cmd = f'mpirun -np {num_cores} {vasp_exe}'
    
    try:
        with open(out_file, 'w') as result:
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=result,
                stderr=sp.PIPE,
                shell=True,
                cwd=dir_cal
               )
    finally:
         pass

    fullname=os.path.join(dir_cal, "OUTCAR")
    return read_energy(fullname)
