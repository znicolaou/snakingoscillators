# snakingoscillators
This repository contains Python scripts to simulate rings of Janus oscillators and arrays of coupled pendula, AUTO scripts to perform numerical continuation, and jupyter notebooks to plot results.

# Repository prerequisites 
Create a conda environment with required packages:
`conda create -n oscillator_env python=3.10 numpy pip wheel matplotlib scipy jupyter ffmpeg pyqt`.  Activate this environment to run the scripts: `conda activate oscillator_env`. 

Our fork of auto-07p contains modified AUTO package for detecting SBPs: https://github.com/znicolaou/auto-07p. Clone and install the package, as in
```
cd ~
git clone https://github.com/znicolaou/auto-07p.git
cd auto-07p
./configure
make
```
The executable files will be compiled into `auto-07p/bin`. You may want to include these in your system PATH by appending your `source $HOME/auto/07p/cmds/auto.env.sh` to `~/.bashrc`.

Batches were run on a SLURM cluster; some script modifications may be necessary if no cluster is available.

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
