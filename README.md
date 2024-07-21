# GARFfield: Genetic Algorithm based Reactive Force Field optimizer

Author: Andres Jaramillo-Botero
California Institute of Technology, 2012

GARFfield is a multi-platform, multi-objective parallel hybrid genetic algorithm (GA) / conjugate-gradient (CG) based force field optimization framework.  It enables first-principles-based force fields prepared from large quantum mechanical data sets, which are now the norm in predictive molecular dynamics simulations for complex chemical processes, as well as from phenomenological data. The former allows improved accuracy and transferability over a wider range of molecular compositions, interactions, and environmental conditions unexplored by experiments.

GARFfield currently supports a range of force field engines, via the LAMMPS Parallel Molecular Dynamics Simulator, including the adiabatic ReaxFF and COMB potentials for modeling reaction processes, the non-adiabatic eFF electron force field with effective core potentials, and Morse potentials (atomistic and coarse-grain).

As opposed to other efforts found in the literature, GARFfield is not limited to a single force field!!. It is intended to handle electronic, atomistic (reactive and non-reactive), and coarse-grain type force field parameter optimization problems, albeit with implicit considerations for each energy/force type of engine to improve global optimization efficiency.

GARFfield provides multiple optimization features that are specific to force field training, including:

- Automatic hooks to the LAMMPS force field engines and automatic force field detection,
- Multi-objective weighted fitness function (charges, energies, lattice parameters, geometric parameters, and others),
- Partial parameter sub-set selection and optimization,
- Training sets with periodic and finite system models,
- Geometry model specifications in extended xyz, biograf (.bgf) and LAMMPS native formats,
- Relative restraints on bond-order based valence interactions, transition state bonds, and electron sizes,
- Non-deterministic solutions in the Pareto-front using random/fixed weighted and un-weighted training sets,
- Systematic hill-climbing option from local minima solutions,
- RMS force fitness for geometric objectives (as opposed to gradient-based energy minimization),
- Deterministic CG minimization switch option when GA is within parabolic minima wells, and
- Others (e.g. parallelization of GA string population)

Support for other force field engines and new features will be added through the modular software architecture design.  Current efforts include template-based force field support (to avoid syntax dependencies from LAMMPS), heuristic sequence training (e.g. valence then non-bond or vice-versa, finite then periodic training cases, etc.), and training set parallelization (in addition to the current GA string population parallelization).

BUILDING

To build GARFfield and all its dependencies:

1. modify paths and compilers in compile.sh 
2. create ./lib/pgapack/lib/arch directory if it does not exist, and modify Makefile
3. source compile.sh compile <arch> (where arch is your architecture's name)

To build GARFfield, with all dependencies:

cd src
make 

The GARFfield executable will be located one directory above ../bin
You can copy this file to a different directory for use.

GARFfield uses modified pgapack, optlist and lammps libraries.

The pgapack and optlist libs are provided in the ./lib directory.

If you need to recompile these:
cd ./lib/pgapack/source
make
cd ../../optlist
cd make

To compile LAMMPS as a library, cd into the lammps src directory and:
make makelib
make -f Makefile.lib openmpi

You will then have to modify the GARFfield Makefile included to
point to the corresponding LAMMPS src directory, and if using pgapack
in parallel, add the corresponding paths to your openmpi and fftw.

RUNNING

In serial mode:

./garffield geo_file ffield_file trainset_file params_file [options]

In parallel mode (not available yet 02/25/2012):

mpirun -np # ./garffield geo_file ffield_file trainset_file params_file [options]

where options, include:

Parameters setup in GA:

-i 'range' (random variables will be produced between a low and a high value) [default]
-i 'percent' (random variables will be bracketed between a mean and a percent deviation from mean)

Both 'range' and 'percent' use the 4th and 5th columns of params_file
If the 3rd column entry contains a number >1 it will switch automatically to percent from a mean value.  I suggest you choose the mean value to be the current ffield entry value, if you're optimizing from an existing force field.

Convergence criteria for GA:

-s 'maxiter' (maximum number of iterations) [default]
-s 'nochange' (no change in total error)

-t # maximum number of iterations [default 400]

Hill-climbing option:

-c performs a hill-climbing routine on a single parameter (uniform).  Will extend to give the user the option to select a number of parameters (randomly), which should have a larger effect on getting our of a local minima.

-f: uses net force calculations instead of energy minimization for geometry cases in training set
-l logfile_name:  write a logfile with errors
-m # (float<1): set GA mutation rate
-p # (int): set population size [default = 100]
-r # (int): set print report frequency
-x: perform mutation AND crossover [defaults to mutation OR crossover]

Output:

The following files are written after convergence

ffield.best:	best force field parameters found
trainset.err: 	individual and accumulated errors from training cases
ffield.new: 	ffield written at every iteration

Troubleshooting:

After quite a few tests, I found that most of the trouble comes from bad input files.
In an attempt to clean these a priori I've written a preprocessor for trainset.in and geo.
This is available in the scripts directory as check_trainset_geo_consistency.py and may
be executed were the original geo and trainset.in files reside w/out arguments.
It will output geo.new and trainset.in.new

NOTE: the force field file must:
1. Not have redundant entries
2. Not have compact torsion entries

If:
- You want to have a more verbose output compile the code with -DDEBUG=1.  It will print every
little detail about the parsing of your geo, trainset.in, ffield and params file.  You can track down where exactly the program failed.
- Your optimization fails in the GA optimization step, look carefully at your ffield file to make sure there are no redundant entries or entries that point to a non-existing atom type.
- Your optimization fails after one single pass of the GA optimization, and it simply dies off, the most likely culprit is the atom counter that writes the number of atoms into the data files in the structures.c file.  This happens on the last structure, when it contains information beyond the HETATM declarations.  FIXED (should not happen any more).
- You run out of memory, use top to follow the memory consumption.  There is an unsolved leak at this point.
- Error does not seem to change.  Check your parameters and make sure the training cases are affected by these.  Check your print frequency, it may be too low to see any significant changes in the total error.

Versions:
04-9-2012: added support for cell parameters
06-5-2012: reads in current force field param values as mean in percent offset (-i) GA optimization
08-8-2012: added CG switch optimization

NOTE2: For updated documentation see GARFfield User's Manual at http://www.wag.caltech.edu/home/ajaramil/GARFFIELDUserManual.pdf

