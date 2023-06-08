import os

from pso.tools.utils import *     
from pso.tools.psoutils import *  

from ase.io import read, write
from pso.tools.operations import *

class PSO(object):
    """
    A class for creating a linear polymer from given RdKit molecules.

    Examples
    ---------


    Note
    -----

    RDKit package must be installed.
    """
    
    __slots__ = ['input_file', 'cluster_file', 'substrate_file']

    def __init__(self, input_file, cluster_file, substrate_file):
       """
       Initialize this class.
          
       Parameters
       -----------

       input_file : class.str
         Name of the input.dat file supplied by the user.
 
       cluster_file : class.str
         Name of the POSCAR file containing the cluster coordinates supplied by the user.

       substrate_file: class.str
         Name of the POSCAR file containing the surface coordinates supplied by the user.
          
       """
       
       self.input_file = input_file
       self.cluster_file = cluster_file 
       self.substrate_file = substrate_file 
       
       
    def pso_struct_search(self): # , input_file, cluster_file, substrate_file):

       input_file= self.input_file
       cluster_file=self.cluster_file
       substrate_file=self.substrate_file
       
       startstep, maxstep, pop, command, cores=read_input_file(str(input_file))
       z_min,acc_pop,omegamax,omegamin,c1,c2=read_input_pso(str(input_file))
    
       #Read the cluster.POSCAR file and convert to Cartesian with selective dynamics style
       clusters = read(str(cluster_file), format='vasp')
       substrate = read(str(substrate_file), format='vasp')
   
       parent_dir=os.getcwd()
       gen_count=startstep

       # Random numbers for struct_parameters, angles and r1,r2 coefficients
       values=global_random_numbers(maxstep, pop, seed=35, output_name="struct_param")
       vel=global_random_numbers(maxstep, pop, seed=42, output_name="vel_param")
       coeff=pso_random_numbers(maxstep, pop, seed=13, output_name="pso_r1_r2")

       #History of gbests
       gbest=[0]*(maxstep-1)
    
       while gen_count < maxstep:
        
           directory="gen_{}".format(str(gen_count))
           path = os.path.join(parent_dir, directory)
        
           try: 
               os.mkdir(path)
            
           except OSError as error: 
               print(error)
        
           if gen_count <= 1:
            
              os.chdir(path)
              gen_dir=os.getcwd()

              x_t=[]; v_t=[]
              for i in range(pop):
                name_new_str="".join(["POSCAR"])
                x_t.append(structure_creator(clusters, substrate, name_new_str, values[gen_count-1][i],z_min))
                v_t.append(tuple(vel[gen_count-1][i]))
                make_config_folders(gen_dir, parent_dir, i, name_new_str)

              subfolders=sub_list(gen_dir)
              gbest[gen_count-1]=ord_dict(vasp_energy(command, subfolders, x_t, v_t, int(cores)))[0]
              os.chdir(parent_dir)
             
           else:
              os.chdir(path)
              gen_dir_new=os.getcwd()
              ord_confs, mask_keys=get_prev_gen(gen_dir_new, acc_pop, gen_count)

              xt_new=[]; vt_new=[]
              vel_new=np.asarray(vel)
              conf_name, gbest_pos=get_gbest_pos(gbest, parent_dir)

              omega_new=omegamax-((omegamax-omegamin)/maxstep)*(gen_count-2)
              print (omega_new)
              
              for i in range(pop):
                  if mask_keys[i]:
                    # Creating folders to store these structures
                    idx=ord_confs[i].split('_')[-1]
                    name_new_str_rest="".join(["POSCAR"])       
                    #The pso propagation scheme
                    pos_update, tmp_vel= pso_update(gbest_pos, ord_confs[i], gen_count,
                                                    omega_new, c1, c2, coeff[gen_count][i][0], coeff[gen_count][i][1]) 
                    xt_new.append(structure_creator(clusters, substrate, name_new_str_rest, pos_update, z_min))
                    vt_new.append(tmp_vel)
                    make_config_folders(gen_dir_new, parent_dir, idx, name_new_str_rest)

                  #Creation of random structures discarded by user-provided acc_pop 
                  else:
                    idxf=ord_confs[i].split('_')[-1]
                    name_new_str_rest="".join(["POSCAR"])
                    xt_new.append(structure_creator(clusters, substrate,
                                                 name_new_str_rest, values[gen_count][i], z_min))
                    vt_new.append(tuple(vel[gen_count][i]))
                    make_config_folders(gen_dir_new, parent_dir, idxf, name_new_str_rest)
                 
              subfolders=sub_list(gen_dir_new)
              gbest[gen_count-1]=list(ord_dict(vasp_energy(command, subfolders, xt_new, vt_new, int(cores))))[0]

              os.chdir(parent_dir)
              print ("Best explored structure: {}".format(conf_name))
      
           gen_count += 1

       #This is needed to compare the last gbest agains the last best
       #structure found in the last folder.

       print ("##############################################")
       last_comparison_gbest(parent_dir, [gbest[-1], conf_name])
       print ("##############################################")
       
