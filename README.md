# snakingoscillators

# Files in the snakingoscillators repository
The script janusrot.py contains python code to numerically integrate a ring of Janus oscillators in the rotating reference frame for numerical continuation of limit cycles. The script pendula.py contains python code to numerically integrate a the heterogeneous pendulum array under parametric driving. The jupyter notebooks lccont.ipynb and plotpendula.ipynb show examples and plots for the respective systems.

# Conda environment 
Create a conda environment with required packages:
`conda create -n oscillator_env -c conda-forge numpy scipy progressbar jupyter matplotlib`.  Activate this environment to run the scripts: `conda activate oscillator_env`. 
