# snakingoscillators
This repository contains Python scripts to simulate rings of Janus oscillators and arrays of coupled pendula, AUTO scripts to perform numerical continuation, and jupyter notebooks to plot results.

# Repository prerequisites 
Create a conda environment with required packages:
`conda create -n oscillator_env numpy=1.23.5 pip wheel matplotlib scipy jupyter ffmpeg pyqt pysindy`.  Activate this environment to run the scripts: `conda activate oscillator_env`. 

Our fork of auto-07p contains modified AUTO package for detecting SBPs: https://github.com/znicolaou/auto-07p. Clone and install the package, as in
```
cd ~
git clone https://github.com/znicolaou/auto-07p.git
mkdir auto && mv auto-07p auto/07p
cd auto-07p
./configure
make
source cmds/auto.env.sh
```
Batches were run on a SLURM cluster; some script modifications may be necessary if no such cluster is available.

# Files in the snakingoscillators repository
The respository contains two main directorys.

## Janus oscillators
The `janus` directory contains a script a Python script `janus.py` contains to numerically integrate a ring of Janus oscillators. Running `./janus.py --help` produces the following usage message:
```
usage: janus.py [-h] --filebase FILEBASE [--output OUTPUT] [--num NUM] [--time TIME] [--beta BETA] [--sigma SIGMA] [--gamma GAMMA] [--omega OMEGA] [--rtime RTIME] [--dt DT]
                [--seed SEED] [--symmetrize SYM]

Numerical integration of networks of phase oscillators.

options:
  -h, --help           show this help message and exit
  --filebase FILEBASE  Base string for file output.
  --output OUTPUT      Output style, 0 for no stdout and sparse output, 1 for stdout and sparse output, 2 for stdout and dense output. Default 1.
  --num NUM            Number of Janus oscillators. Default 16.
  --time TIME          Total integration time. Detault 2000.
  --beta BETA          Internal coupling strength. Default 0.25.
  --sigma SIGMA        Coupling strength. Default 1.
  --gamma GAMMA        Coefficient for amplitude terms. Default 0.1.
  --omega OMEGA        Natural frequency gap. Default 1.0.
  --rtime RTIME        Time to start averaging order parameter. Default 1000.
  --dt DT              Time step for averaging and output. Default 0.1.
  --seed SEED          Initial condition random seed. Default 1.
  --symmetrize SYM     Symmetrize the coupling. 1 for symmetric, 0 for chiral. Default 0.
```

Batches of this script are used to generate the initial limit cycles. To do so, the SLURM batch script `sbatch sweep.sh`, which should take an hour or so to complete. After the jobs complete, run the script `./makecycles.py` to generate the starting data for the AUTO continuation. 

To generate the main numerical continuation results, submit the SLURM script `sbatch job.sh`, which should take over a day to complete. After the job completes, run `./makesteadys.sh` and `sbatch steadyjob.sh` to find and numerically continue the steady states visited by the heteroclinic cycles (these jobs run in a few minutes).  

Some solution branches generated in the continuation are neutrally stable invariant limit cycles, which possess symmetry-constrained unity Floquet multipliers. The direciton vector in such branches can spuriously change signs because of the singularlity of the Jacobian, and are thus difficult to continue. The script `invariant.auto` contains code for tracking the number of neutrally stable multipliers and ensuring that the direction of the change in the order parameter does not change. (This was designed specifically for the exemplary solution branch; additional modifications would be required for branches that do exhibit extrema in the order parameter.) Run `sbatch invariantjob.sh` to run the script on the exemplary branch in the manuscript. 

Finally, the Jupyter notebook `plotjanus.ipynb` contains code to plot and animate solutions. A few example solution branches are stored in `data/janus/examples`, which are plotted in the notebook.

## Pendula arrays
The `pendula` directory contains a script a Python script `pendula.py` contains to numerically integrate an array of parameterically driven pendula. Running `./pendula.py --help` produces the following usage message:
```
usage: pendula.py [-h] --filebase FILEBASE [--num NUM] [--frequency FREQ] [--initfrequency FREQ0] [--initamplitude AMP0] [--amplitude AMP] [--delta DELTA] [--cycles CYCLES]
                  [--initcycle INITCYCLE] [--outcycle OUTCYCLE] [--dt DT] [--seed SEED] [--damp DAMP] [--spring SPRING] [--init INIT] [--rtol RTOL] [--atol ATOL] [--verbose VERBOSE]

Driven pendula.

options:
  -h, --help            show this help message and exit
  --filebase FILEBASE   Base string for file output
  --num NUM             Number of pendula
  --frequency FREQ      Driving frequency
  --initfrequency FREQ0
                        Initial driving frequency
  --initamplitude AMP0  Driving amplitude
  --amplitude AMP       Driving amplitude
  --delta DELTA         Alternating pendulum length scale
  --cycles CYCLES       Simulation time in driving cycles
  --initcycle INITCYCLE
                        Simulation time in driving cycles
  --outcycle OUTCYCLE   Cycle to start outputting
  --dt DT               Time between outputs in driving cycles
  --seed SEED           Seed for random initial conditions
  --damp DAMP           Damping coefficient
  --spring SPRING       Spring coefficient
  --init INIT           Initial random scale
  --rtol RTOL           Relative error tolerance
  --atol ATOL           Absolute error tolerance
  --verbose VERBOSE     Verbose output
```

Batches of this script are used to generate the initial limit cycles. To do so, the SLURM batch script `sbatch sweep.sh`, which should take an hour or so to complete. After the jobs complete, run the script `./makecycles.py` to generate the starting data for the AUTO continuation. 

To generate the main numerical continuation results from the starting cycles, submit the SLURM script `sbatch job.sh`, which should take a couple of hours to complete. To find the initial period doubling instabilities of the homogeneous state and the branches emerging from them, submit the batch `sbatch flatjob.sh`, which should complete in an hour or so. 

The Jupyter notebook `plotpendula.ipynb` contains code to plot and animate solutions. A few example solution branches are stored in `data/pendula/examples`, which are plotted in the notebook.

