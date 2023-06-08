# KCL_PSO:  King's College London Particle Swarm Optimisation program.

This is a program implemeting the Particle Swarm Optimisation (PSO) algorithm for searching the global minimum of 
a series of clusters placed on top of surfaces.

# Pip installation

To quickly install **KCL_PSO**, we encourage to do it inside a virtual environment, which can be achieved in the following way:

1. Create a directory named as you want and access it (in this case called work_pol):

```bash 
  [~] mkdir work_pso
  [~] cd work_pso
```

2. **IN CASE YOU ARE INSTALLING THE PACKAGE IN YOUNG** You need first to load these modules **BEFORE INSTALLING THE LIBRARRY**.
To do so, please type the following commands in your session:

 ```bash
  module purge
  module load gcc-libs/4.9.2
  module load openblas/0.3.7-serial/gnu-4.9.2
  module load python3/3.8
 ```

3. Create a virtual environment named pol (this name can be changed, of course) and activate the environment:

```bash 
   [~] python -m venv kcl_pso
   [~] source kcl_pso/bin/activate
```

3. Get KCL_PSO from the GitHub repository:

```bash 
  git clone https://github.com/alejandrosantanabonilla/PSO.git
```

4. Install KCL_PSO in this folder using the virtual environment

```bash 
   [~] cd PSO
   [~] pip install .
```

# First calculation

To test the code, there is a folder inside the **PSO** directory called
**test**. Inside there is a folder called **input_files**, where all the 
corresponding VASP files **EXCEPT POTCAR** are located. Please add **POTCAR**
in that directory before starting any calculation.

You can submit a first test, just by typing:

``` bash
    
    qsub job.sh
```

This should work. To check your job, you can type

``` bash

    qstat -u mmmXXXX
```
and look for pso_calculation job.




