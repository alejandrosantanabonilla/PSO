import os
#from operations import *
import shutil
import json

def read_input_file(input_file):
    """ Function to read a valid input file for searching 
        structures.

        Parameters:
        ===========

        input_file: str
             Text file defining a list of parameters to perform a PSO search.


        Return:
        =======
        
        a,b,c,d: int, int, int, str
             List of variables defining the starting point, maximum step, 
             population size, and submission command for a valid VASP 
             executable.
    """

    with open(str(input_file)) as f:
       for line in f:
          if line.startswith("startstep"):
             temp_st = line.partition('=')
         
          if line.startswith("maxstep"):
             temp_mt = line.partition('=')
         
          if line.startswith("popsize"):
             temp_ps = line.partition('=')

          if line.startswith("command"):
             temp_com = line.partition('=')

          if line.startswith("cores"):
             temp_cores = line.partition('=')


    return int(temp_st[-1]), int(temp_mt[-1]), int(temp_ps[-1]), str(temp_com[-1]), int(temp_cores[-1])

def read_input_pso(input_file):
    """Function to read a valid input file for searching 
       structures.

       Parameters:
       ===========

       input_file: str
            Text file defining a list of parameters to perform a PSO search.

       Return:
       =======
        
       a,b,c,d,e: float, int, float, float, float
            List of variables defining the minimum height, number of best 
            candidates per population, omega, c1 and c2 coefficients.
    """

    with open(str(input_file)) as f:
       for line in f:
          if line.startswith("z_min"):
             temp_zmin = line.partition('=')
         
          if line.startswith("best_idv"):
             temp_acc_pop = line.partition('=')
         
          if line.startswith("omegamax"):
             temp_omegamax = line.partition('=')

          if line.startswith("omegamin"):
             temp_omegamin = line.partition('=')
           
          if line.startswith("c1"):
             temp_c1 = line.partition('=')

          if line.startswith("c2"):
             temp_c2 = line.partition('=')


    return float(temp_zmin[-1]), int(temp_acc_pop[-1]), float(temp_omegamax[-1]), float(temp_omegamin[-1]) ,float(temp_c1[-1]), float(temp_c2[-1])


def make_config_folders(gen_dir, parent_dir, number, name_new_str):
    """ Function to create folders containing different 
        configurations.

        Parameters:
        ===========

        gen_dir: str
             Path to a generation directory.

        parent_dir: str
             Path to the parent directory.

        number: int
             Number used to label the new name of the created folder.

        name_new_str: str
             Name used to label de folder.

        Return:
        =======
        
        None:
           Moving all directories to the provided paths.
    """
    
    conf_dir="config_{}".format(number)
    path = os.path.join(gen_dir, conf_dir)

    try: 
        os.mkdir(path) 

    except OSError as error: 
        print(error)

    dst_path=os.path.join(conf_dir,path)
    src_path=os.path.join(gen_dir, name_new_str)
    shutil.move(src_path, dst_path)

def sub_list(gen_dir):
  """ Function to list all sub-directories from an user-provided 
      path. 

      Parameters:
      ============

      gen_dir: str
          Path to a specific generation directory.

      Returns:
      ========

      sub: class.list
          List with all the found directories for the provided path.
  """
  sub=[f.path for f in os.scandir(gen_dir) if f.is_dir()]

  return sub



def copy_files_vasp(vasp_files, dest_dir):
    """ Function to copy the input files for perfoming a 
        VASP calculation

        Parameters:
        ============

        vasp_files: str
            Path to the folder containing the input VASP files.

        dest_dir: str
            Path to the destinated directory.

        Returns:
        ========

        None: 
            All VASP input files move to the designated directory.
    """
    
    inputs=[f for f in os.listdir(vasp_files)]

    for idx, values in enumerate(inputs):
        shutil.copy(os.path.join(vasp_files, str(values)), dest_dir)


def ord_dict(dict_user):
    """Function to organize a dictionary based on the values of its 
       items. It is reversely ordered.

       Parameters:
       ===========

        dict_user: class.dict
          Dictionary containing the information of the state of each PSO-particle.

       Returns:
       =========

       None:
         List of sorted keys from an user-provided dictionary.
       
    """

    #This is an arranged dictionary with the x_t positions and velocities

    sorted_data={k: v for k, v in sorted(dict_user.items(),
                                         key=lambda r: r[1][0],
                                         reverse=True)}

    return list(sorted_data.keys())

def print_dic_gen(dict_user, fullname):
    """Function to print a dictionary into a json file.

       Parameters:
       ============

       d: class.dict
          Dictionary containing the information of the state of each PSO-particle.

       Returns:
       =========

       None:
          JSON file where the user-provided dictionary has been stored.

    """
    with open(fullname, "w") as outfile:
       json.dump(dict_user, outfile)

def open_json_file(key):
    """ Function to open JSON files where the velocities and positions 
        of each particle are stored.
    """
    
    base,folder=os.path.split(os.getcwd())
    full_path=key.split('_')
    gen_name="_".join(full_path[0:2])

    basename=os.path.join(base, gen_name)
    fullname=os.path.join(basename,"pop_res.json")
    
    with open(fullname, "r") as f:
         pop_info = json.load(f)

    return pop_info
