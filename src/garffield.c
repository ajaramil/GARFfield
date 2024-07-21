/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include <sys/types.h>
#include <dirent.h>
#include <fcntl.h> 
#include <sys/stat.h>
#include <unistd.h> // for close

#include "optlist.h"
#include "garffield.h"
#include "params.h"
#include "fmincg.h"

#include "structures.h"
#include "trainset.h"
#include "restraints.h"

#ifndef M_PI
#define M_PI 3.14159265354
#endif

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

double version = 1.0;

/****************************************************************************
*  Function     : main
*  Description  : prepares force field arguments and parameters for GA and
                  hillclimb, evolves random parameter search
*  Parameters   : argv [1=geo_file;2=ffield_file;3=trainset_file,
                  4=parameter_file; 5=restraints_file (optional)] [options]
*  Effects      : Allocates global variables and controls de MPI environment
****************************************************************************/

int main (int argc, char **argv)
{
  PGAContext *ctx;		/* the context variable */
  int maxiter;			/* the maximum number of iterations */
  int numParams, numRestraints, hillclimbing_flag, realrange_flag, hillperiod;
  option_t *optList, *thisOpt;
  int iter_flag, best_p, logfile_flag, cg_flag;
  char logfile[MAX_LINE], restartfile[MAX_LINE];
  FILE *log, *fp2;
  int pop_size, pop_rep, i;
  double mutation_prob, crossover_prob;
  int rank;
  void (*CreateNewGeneration) (PGAContext *, int, int);
  int restart, Restarted, report_frequency, iteration, restart_frequency, cgswitch, stopGA;
  div_t divresult;
  char ffname[MAX_LINE], tsetname[MAX_LINE], paramsname[MAX_LINE], fforiginal[MAX_LINE],
    restraintsname[MAX_LINE];
  int numprocs;
  double chrom_error;
  ffieldtype *ff_struc_pointer;
  char *line, **strtmp;
  int token,cg_ret;
  double *Vec;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  line = (char *) scalloc (MAX_LINE, sizeof (char), "line3");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp3");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp3i");

  init (argc, argv);

  if ((argc <= 4) && (rank == 0)) {
    printf ("GARFfield force field parameter optimization framework, version %2.1f\n", version);
    printf ("Author: Andres Jaramillo-Botero (ajaramil@caltech.edu, Caltech (c) 2010-2014)\n\n");
    printf ("Usage:\n");
    printf
      ("%s <geo> <ffield> <trainset> <params> [<restraints>] [GA options] [other options] ... or\n",
       argv[0]);
    printf
      ("mpirun -np NumProcs %s <geo> <ffield> <trainset> <params> [<restraints>] [GA options] [other options]\n\n",
       argv[0]);
    printf ("GA options are: \n");
    printf ("  -i ['range' or 'percent'] {parameter interval}\n");
    printf ("  -g ['GRGA' or 'SSGA'] {replacement strategy: generational or steady-state}\n");
    printf ("  -m [float between 0-1] {mutation rate}\n");
    printf ("  -o [float between 0-1] {crossover rate}\n");
    printf ("  -p [int] {population size}\n");
    printf ("  -s ['maxiter' or 'nochange'] {convergence criteria}\n");
    printf ("  -t [int] {stop after int iterations}\n");
    printf ("  -z [int] {set SSGA population replacement percent}\n");
    printf ("  -x {perform mutation and crossover}\n");
    printf ("Other options are:\n");
    printf ("  -c [int] {Hillclimb every int iterations}\n");
    printf ("  -C [int] {mincap array size in reax/c allocations}\n");
    printf ("  -d ['full_debug', 'lammps_debug'] debug mode\n");
    printf ("  -e ['MPE' or 'SIWE' or 'RMSE' or 'NRMSE' or 'MSE'] {error function}\n");
    printf ("  -f {check forces on unit cells, instead of performing minimization}\n");
    printf ("  -F ['reax']\n");
    printf ("  -G [int] {switch to Conjugate Gradient after int iterations}\n");
    printf ("  -I [int] {Number of molecular mechanics energy minimization steps}\n");
    printf ("  -l <filename> {write log to filename, instead of console}\n");
    printf ("  -M ['fire' or 'sd'] {switch structure minimization algorithm - default is 'cg'}\n");
    printf ("  -q <filename> {Use pQeq for ReaxFF charges}\n");
    printf ("  -r [int] {print report every int iterations}\n");
    printf ("  -R [int] {Restart every int interations}\n\n");
    printf ("  -S [float] {safezone factor in reax/c memory allocation}\n");
    printf
      ("  -u ['real'] {Force units to be real, i.e. lengths in Angstroms, energies in kcal/mol and forces in kcal/mol-Angstrom}\n");
    printf ("  -v {Use low-gradient van der Waals for reax/c}\n");
    printf
      ("  -w <filename> {write best force field file at same print report report_frequency}\n");
    printf ("  -W {randomize section weights in training file}\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
// Establish number of parameters and ranges from user-provided params file
  if (rank == 0)
    numParams = get_num_params (params_file, rank);

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Bcast (&numParams, 1, MPI_INT, 0, MPI_COMM_WORLD);

// Create GA
  ctx = PGACreate (&argc, argv, PGA_DATATYPE_REAL, numParams, PGA_MINIMIZE);

  PGASetRandomSeed (ctx, 1);

// set user defaults
  eff_units = 0;		// default is 'electron' (0) units for eFF-ECP, 'real' (1) otherwise
  maxiter = 400;
  hillclimbing_flag = 0;	// no hill-climbing
  hillperiod = 100;		// if hillcliming flag changed, the default period is 100 iterations
  logfile_flag = 0;		// do not write a log file
  iter_flag = 0;		// detects presence of -t #
  cg_flag = 0;			// CG is off
  cgswitch = 50;		// Switch to CG after 50 iterations
  report_frequency = 10;	// report report_frequency
  forcefield2optimize = REAXC;	// use reaxFF C version (0=Fortran version or 2=eFF-ECP)
  pop_size = 100;		// GA population size
  pop_rep = 0;			// SSGA=0 or GRGA=1 or eGRGA=2
  min_method = 0;		// cg=0 or fire=1 or sd=2
  realrange_flag = 1;		// interval set to range
  errorfunction = 1;		// Mean Percent Error (normalized)
  normalizederrors = 1;		// Normalized errors
  debug_level = 0;		// None
  restart = 0;			// Do not write a restart
  restart_frequency = 50;	// Default restart report_frequency
  struc_min_steps = 500;	// Default minimization steps
  lg_reax = 0;			// Not using low-gradient dispersion correction
  safezone_flag = 0;
  safezone_value = 1.6;		// Memory allocation factor in reaxff
  mincap_flag = 0;
  mincap_value = 100;		// Memory allocation array size in reaxff
  stopGA = 0;			// Stops GA cycle if 1
  pqeq_flag = 0;

// Customize GA and other parameters from command line
  if (rank == 0) {
    optList = NULL;
    optList = GetOptList (argc, argv, "i:I:s:S:t:l:p:g:G:z:m:M:o:q:u:r:R:w:F:c:C:d:e:hxfvW");

    // set defaults for GA
    PGASetStoppingRuleType (ctx, PGA_STOP_MAXITER);	// stop at max iterations
    PGASetMaxGAIterValue (ctx, maxiter);	// maximum GA iterations = 400
    PGASetMutationOrCrossoverFlag (ctx, 1);

    // Defaults 
    printf ("GARFfield force field parameter optimization framework, version %2.1f\n", version);
    printf ("(c) Andres Jaramillo-Botero, ajaramil@caltech.edu (Caltech, 2010-2014)\n\n");
    printf ("1. DEFAULT SETTINGS\n");
    printf ("== GA specific ==\n");
    printf ("  -i='range' (Parameter interval)\n");
    printf ("  -s='maxiter' (Stopping criteria)\n");
    printf ("  -t=%d (Maximum GA iterations)\n", maxiter);
    printf ("  -p=100 (GA Population size)\n");
    printf
      ("  -g='SSGA' (Population replacement strategy, steady-state with 10 percent population size replacement)\n");
    printf ("  -z=10 (SSGA population replacement percent)\n");
    printf ("  Mutation only (use -x to apply mutation and crossover operations)\n");
    printf ("  -m=1/%d (Mutation rate)\n", numParams);
    printf ("  -o=0.85 (Crossover rate = 0.85)\n\n");
    printf ("== Others ==\n");
    printf ("  -F='reaxc' (Force field)\n");
    printf ("  -I=500 (structure minimization steps)\n");
    printf ("  No hillclimbing (-c)\n");
    printf ("  No logfile (use -l to write a logfile)\n");
    printf ("  Energy minimization (use -f to use RMS force criterion instead)\n");
    printf ("  -e=MPE (Mean Percent Error)\n");
    printf ("  -d=None (Debug mode)\n");
    printf ("  -r=%d (Report period, in GA iterations)\n", report_frequency);
    printf ("  -u='electron' for eFF-ECP and 'real' otherwise (Units)\n");
    printf ("  Use regular Qeq/EEM solver for ReaxFF charges\n");
    printf ("  No low-gradient van der Waals correction (only reax/c)\n");
    printf ("  No periodic write of best force field (-w)\n");
    printf ("  Default to no Conjugate Gradient minimization in parabolic minima well (-G)\n");
    printf ("  Default structure optimization method is CG (-M)\n");
    printf ("  Fixed weights in training file (-W) \n");
    printf ("  No mincap or safezone memory allocation in reax/c calculations (-S)\n");
    printf ("  No restarts (-R)\n");

    /* parse command line argument options */
    if (optList != NULL)
      printf ("== User changes to default settings ==\n");
    while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;

      switch (thisOpt->option) {
      case 'h':
        printf ("GARFfield force field parameter optimization framework, version %2.1f\n", version);
        printf ("(c) Andres Jaramillo-Botero, ajaramil@caltech.edu (Caltech, 2010-2014)\n\n");
        printf ("Usage:\n");
        printf
          ("%s <geo> <ffield> <trainset> <params> [<restraints>] [GA options] [other options] ... or\n",
           argv[0]);
        printf
          ("mpirun -np NumProcs %s <geo> <ffield> <trainset> <params> [<restraints>] [GA options] [other options]\n\n",
           argv[0]);
        printf ("Options:\n");
        printf ("  -c [#] perform hillclimb [#=period]\n");
        printf ("  -C [#] {mincap array size in reax/c allocations}\n");
	printf ("  -d ['full_debug' or 'lammps_debug'] debug mode\n");
	printf
	  ("  -e ['MPE' (mean percent) or 'SIWE' (inverse weighted) or 'RMSE' (root mean square) or 'NRMSE' (normalized RMS) or 'MSE' (mean square)] error function\n");
	printf ("  -f [] use net force for geometries\n");
	printf ("  -F ['reax'] switch to Fortran version of ReaxFF (default is C version)\n");
	printf
	  ("  -g ['SSGA' or 'GRGA' or 'eGRGA'] steady-state, generational, and elite generational population replacement strategy\n");
	printf ("  -G [#] switch to Conjugate Gradient after # iterations\n");
	printf ("  -h [] print out command line options\n");
	printf ("  -i ['range' or 'percent'] sets parameter interval\n");
	printf ("  -l <filename> write log to filename\n");
	printf ("  -m [#] sets mutation rate\n");
        printf ("  -M ['fire' or 'sd'] switch structure minimization algorithm\n");
	printf ("  -I [#] sets structure minimization steps\n");
	printf ("  -o [#] sets crossover rate\n");
	printf ("  -p [#] sets GA population size\n");
        printf ("  -q <filename> use pQeq for ReaxFF charges\n");
	printf ("  -r [#] sets print report period\n");
	printf ("  -R [#] sets restart period\n");
	printf ("  -s ['maxiter' or 'nochange'] sets GA stopping criteria\n");
        printf ("  -S [#] {safezone factor in reax/c memory allocations}\n");
        printf ("  -t [#] sets maximum GA iterations (if -s 'maxiter')\n");
        printf ("  -x [] perform mutation AND crossover\n");
        printf ("  -u ['real'] set eFF-ECP units to be real\n");
        printf ("  -v [] {Use low-gradient van der Waals for reax/c}\n");
        printf ("  -w <filename> write best force field to filename with same report period\n");
        printf ("  -W [] sets random section weights in training file\n");
        printf ("  -z [#] sets SSGA population replacement percent\n");
        printf ("Try it again ... (aborting now)\n\n");
        printf ("Author: Andres Jaramillo-Botero (ajaramil@caltech.edu) @2012\n");

	FreeOptList (thisOpt);	/* free the rest of the list */
	MPI_Abort (MPI_COMM_WORLD, 1);

      case 'i':		// set parameter interval
	if (thisOpt->argument != NULL) {
	  if (strcmp (thisOpt->argument, "range") == 0) {
	    // between low and high values
	    realrange_flag = 1;
//          PGASetRealInitRange (ctx, Lower, Upper);
	  }
	  else if (strcmp (thisOpt->argument, "percent") == 0) {
	    printf (">> -%c = %s (Parameter interval)\n", thisOpt->option, thisOpt->argument);
	    // percent offset (Upper) from a mean value (Lower)
	    realrange_flag = 0;
//          get_params (params_file, realrange_flag);
//          PGASetRealInitPercent (ctx, Lower, Upper);
	  }
	  else {
	    printf ("ERROR: Incorrect argument for -%c\n", thisOpt->option);
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
	free (thisOpt);
	break;

      case 'd':
	if (thisOpt->argument != NULL) {
	  printf (">> -%c = %s (Debug mode)\n", thisOpt->option, thisOpt->argument);
	  if (strcmp (thisOpt->argument, "full_debug") == 0)
	    debug_level = 1;
	  else {
	    if (strcmp (thisOpt->argument, "lammps_debug") == 0)
	      debug_level = 2;
	    else {
	      printf ("ERROR: Incorrect argument for -%c\n", thisOpt->option);
	      MPI_Abort (MPI_COMM_WORLD, 1);
	    }
	  }
	}

	free (thisOpt);
	break;

      case 'e':
	if (thisOpt->argument != NULL) {
	  printf (">> -%c = %s (Error normalization in objective functions)\n", thisOpt->option,
		  thisOpt->argument);
	  if (strcmp (thisOpt->argument, "MPE") == 0)
	    errorfunction = 1;
	  else {
	    if (strcmp (thisOpt->argument, "RMSE") == 0)
	      errorfunction = 2;
	    else {
	      if (strcmp (thisOpt->argument, "NRMSE") == 0)
		errorfunction = 3;
	      else {
		if (strcmp (thisOpt->argument, "MSE") == 0) {
		  errorfunction = 4;
		  normalizederrors = 0;
		}
		else {
		  if (strcmp (thisOpt->argument, "SIWE") == 0) {
		    errorfunction = 0;
		    normalizederrors = 0;
		  }
		  else {
		    printf ("ERROR: Incorrect argument for -%c\n", thisOpt->option);
		    MPI_Abort (MPI_COMM_WORLD, 1);
		  }
		}
	      }
	    }
	  }
	}

	free (thisOpt);
	break;

      case 'q':         // Use p-Qeq for charges (requires pqeq par file)
        if (thisOpt->argument != NULL) {
          printf (">> -%c %s (Use p-Qeq for charge calculation)\n", thisOpt->option,thisOpt->argument);
          strcpy(pqeq_par_file, thisOpt->argument);
          pqeq_flag = 1;
        } else {
          printf ("ERROR: Provide a pQeq charge parameter file name\n");
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
        free (thisOpt);
        break;

      case 'r':
	if (thisOpt->argument != NULL) {
	  printf (">> -%c = %d (Report frequency, in GA iterations)\n", thisOpt->option,
		  atoi (thisOpt->argument));
	  report_frequency = atoi (thisOpt->argument);
	}
	else {
	  printf ("ERROR: Provide name of print report period\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'R':
	if (thisOpt->argument != NULL) {
	  restart_frequency = atoi (thisOpt->argument);
	  printf (">> -%c = %d (Set restart frequency)\n", thisOpt->option, restart_frequency);
	  PGASetRestartFlag (ctx, PGA_TRUE);
	  PGASetRestartFrequencyValue (ctx, restart_frequency);
	}
	else {
	  printf ("ERROR: Provide name of restart period\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'w':
	if (thisOpt->argument != NULL) {
	  restart = 1;
	  strcpy (restartfile, thisOpt->argument);
	  printf (">> -%c (Write %s restart every %d iterations)\n", thisOpt->option,
		  thisOpt->argument, report_frequency);
	}
	else {
	  printf ("ERROR: Provide name of restart file\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'W':		// radomize section weights in training file
	printf (">> -%c (Randomize section weights in training file)\n", thisOpt->option);
	random_wflag = 1;
	free (thisOpt);
	break;

      case 's':		// stopping rule
	if (thisOpt->argument != NULL) {
	  if (strcmp (thisOpt->argument, "nochange") == 0) {
	    // no change in generation output
	    PGASetStoppingRuleType (ctx, PGA_STOP_NOCHANGE);
	    printf (">> -%c = %s (Stopping criteria)\n", thisOpt->option, thisOpt->argument);
	  }
	}
	else {
	  printf ("ERROR: Incorrect argument for -%c\n", thisOpt->option);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'S': 
        // Factor in reax/c memory allocation
        safezone_flag = 1;
        if (thisOpt->argument != NULL) {
          safezone_value = atof (thisOpt->argument);
          printf (">> -%c = %2.2f (Safezone factor in reaxFF memory allocation)\n", thisOpt->option,safezone_value);
        }
        else {
          printf ("ERROR: Provide option -S [array allocation factor for reax/c]\n");
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
        free (thisOpt);
        break;

      case 'C':
        // Minimum size for array allocation in reax/c
        mincap_flag = 1;
        if (thisOpt->argument != NULL) {
          mincap_value = atoi (thisOpt->argument);
          printf (">> -%c = %d (Minimum size for array allocation in reax/c)\n", thisOpt->option,mincap_value);
        }
        else {
          printf ("ERROR: Provide option -c [minimum array size allocation for reax/c]\n");
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
        free (thisOpt);
        break;

      case 't':		// set max number of iterations for maxiter stopping rule
	iter_flag = 1;
	if (thisOpt->argument != NULL) {
	  // maximum number of user-defined iterations
	  maxiter = atoi (thisOpt->argument);
          if (maxiter == 0)
            printf ("Running single step without parameter optimization\n");
          else {
	    printf (">> -%c = %d (Maximum GA iterations)\n", thisOpt->option, maxiter);
	    PGASetMaxGAIterValue (ctx, maxiter);
	    PGASetStoppingRuleType (ctx, PGA_STOP_MAXITER);
          }
	}
	else {
	  printf ("ERROR: Provide option -t [number of iterations]\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'c':		// choose hybrid GA-hill-climbing
	if (thisOpt->argument != NULL) {
	  // maximum number of user-defined iterations
	  hillperiod = atoi (thisOpt->argument);
	  printf (">> -%c = %d (Use hillclimbing)\n", thisOpt->option, hillperiod);
	}
	else {
	  printf ("ERROR: Provide option -c [hillclimb_period]\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	hillclimbing_flag = 1;
	free (thisOpt);
	break;

      case 'p':		// change GA population size
	if (thisOpt->argument != NULL) {
	  pop_size = atoi (thisOpt->argument);
	  printf (">> -%c =%d (Population size)\n", thisOpt->option, pop_size);
	  PGASetPopSize (ctx, pop_size);
	}
	else {
	  printf ("ERROR: provide population size\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'g':
	if (thisOpt->argument != NULL) {
	  printf (">> -%c = %s ", thisOpt->option, thisOpt->argument);
	  if (strcmp (thisOpt->argument, "GRGA") == 0) {
	    printf ("Generational Replacement strategy)\n");
	    pop_rep = 1;
	  }
	  else if (strcmp (thisOpt->argument, "eGRGA") == 0) {
	    printf ("(Elitist Generational Replacement strategy)\n");
	    pop_rep = 2;
	  }
	  else if (strcmp (thisOpt->argument, "SSGA") == 0) {
	    printf ("(Steady-State Replacement strategy)\n");
	    pop_rep = 0;
	  }
	}
	else {
	  printf ("ERROR: provide population replacement strategy (GRGA or eGRGA, default SSGA)\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'G':
	if (thisOpt->argument != NULL) {
	  printf (">> -%c = %s ", thisOpt->option, thisOpt->argument);
	  cgswitch = atoi (thisOpt->argument);
          cg_flag = 1;
	}
	else {
	  printf ("ERROR: Check command line option for CG optimization)\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'z':
	if (thisOpt->argument != NULL) {
	  if (pop_rep == 0) {
	    printf (">> -%c = %s (Changing SSGA population replacement percent)\n", thisOpt->option,
		    thisOpt->argument);
	    PGASetNumReplaceValue (ctx, atoi (thisOpt->argument));
	  }
	  else {
	    printf ("ERROR: Can only set percent population replacement with SSGA strategy\n");
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
	else {
	  printf ("ERROR: provide SSGA percent population replacement\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'l':
	if (thisOpt->argument != NULL) {
	  printf (">> -%c = %s (Write to logfile)\n", thisOpt->option, thisOpt->argument);
	  strcpy (logfile, thisOpt->argument);
	  logfile_flag = 1;
	}
	else {
	  printf ("ERROR: Provide name of logfile\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'I':		// change number of structure minimization steps in MM code
	if (thisOpt->argument != NULL) {
	  struc_min_steps = atoi (thisOpt->argument);
	  printf (">> %c = %d (Structure minimization steps)\n", thisOpt->option, struc_min_steps);
	}
	else {
	  printf ("ERROR: set the number of structure (energy) minimization steps\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'm':		// change GA mutation rate
	if (thisOpt->argument != NULL) {
	  mutation_prob = atof (thisOpt->argument);
	  printf (">> %c = %-6.4f (Gene mutation rate)\n", thisOpt->option, mutation_prob);
	  PGASetMutationProb (ctx, mutation_prob);
	}
	else {
	  printf ("ERROR: set the mutation rate value between (0-1.0]\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'M':		// change minimization method
        if (thisOpt->argument != NULL) {
          printf (">> -%c = %s ", thisOpt->option, thisOpt->argument);
          if (strcmp (thisOpt->argument, "fire") == 0) {
            printf ("(FIRE structure minimization)\n");
            min_method = 1;
          }
          else if (strcmp (thisOpt->argument, "sd") == 0) {
            printf ("(Steepest Descent structure minimization)\n");
            min_method = 2;
          }
        }
        else {
          printf ("ERROR: define structure minimization algorithm alternative to cg (fire or sd)\n");
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
        free (thisOpt);
        break;

      case 'o':		// change GA crossover rate
	if (thisOpt->argument != NULL) {
	  crossover_prob = atof (thisOpt->argument);
	  printf (">> -%c = %-6.4f (Gene crossover rate)\n", thisOpt->option, crossover_prob);
	  PGASetCrossoverProb (ctx, crossover_prob);
	}
	else {
	  printf ("ERROR: set the crossover rate value between (0-1.0]\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'u':		// force units to be real (especially for eFF-ECP)
	if ((thisOpt->argument != NULL) && (strcmp (thisOpt->argument, "real") == 0)) {
	  eff_units = 1;
	  printf (">> -%c = real (force real units in eFF-ECP)\n", thisOpt->option);
	}
	else {
	  printf ("ERROR: only option is 'real' units\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;

      case 'x':		// do mutation and crossover
	printf (">> -%c (Perform mutation and crossover)\n", thisOpt->option);
	PGASetMutationAndCrossoverFlag (ctx, 1);
	PGASetMutationOrCrossoverFlag (ctx, 0);
	free (thisOpt);
	break;

      case 'f':
	printf (">> -%c (Compute error from forces in geometries)\n", thisOpt->option);
	minimization = 0;
	free (thisOpt);
	break;

      case 'v':
	printf (">> -%c (Low-gradient van der Waals correction for reax/c)\n", thisOpt->option);
	lg_reax = 1;
	free (thisOpt);
	break;

      case 'F':
	if (thisOpt->argument != NULL) {
	  if (strcmp ("reax", thisOpt->argument) == 0) {
	    printf (">> -%c = %s (Force field)\n", thisOpt->option, thisOpt->argument);
	    forcefield2optimize = REAX;
	  }
/*          else if (strcmp ("reaxc", thisOpt->argument) == 0)
	    forcefield2optimize = REAXC;
	  else if (strcmp ("eff", thisOpt->argument) == 0) {
            printf (">> -%c = %s (Force field)\n", thisOpt->option, thisOpt->argument);
	    forcefield2optimize = EFF;
          } else if (strcmp ("pqeq", thisOpt->argument) == 0) {
            printf (">> -%c = %s (Force field)\n", thisOpt->option, thisOpt->argument);
            forcefield2optimize = PQEQ;
          } else if (strcmp ("ppr", thisOpt->argument) == 0) {
            printf (">> -%c = %s (Force field)\n", thisOpt->option, thisOpt->argument);
            forcefield2optimize = PPR;
          } else if (strcmp ("cg", thisOpt->argument) == 0) {
            printf (">> -%c = %s (Force field)\n", thisOpt->option, thisOpt->argument);
            forcefield2optimize = CG;
          } else if (strcmp ("morse", thisOpt->argument) == 0) {
            printf (">> -%c = %s (Force field)\n", thisOpt->option, thisOpt->argument);
            forcefield2optimize = MORSE;
	  } 
*/
	  else {
	    printf ("ERROR: Can only switch from reaxc to reax with -F option\n");
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
	else {
	  printf ("ERROR: Select which force field to optimize, between reaxc, and reax\n");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	free (thisOpt);
	break;
      }
    }

    // Set population replacement strategy
    if (pop_rep == 1)
      PGASetNumReplaceValue (ctx, pop_size);
    else if (pop_rep == 2)
      PGASetNumReplaceValue (ctx, pop_size - 1);	// elitist
    PGASetNoDuplicatesFlag (ctx, PGA_TRUE);
    PGASetMutationBoundedFlag (ctx, PGA_TRUE);	// reset to lower (upper) value if it falls outside of range

    // Allocate arrays and populate Lower and Upper parameter arrays
    Lower = (double *) scalloc (numParams, sizeof (double), "lower");
    Upper = (double *) scalloc (numParams, sizeof (double), "upper");
    get_params (params_file, realrange_flag, rank);

    /* Read the force field type from the force field file */
    if ((fp2 = fopen (ff_file, "r")) == NULL) {
      fprintf (stderr, "ERROR: opening the force field file! terminating...\n");
      MPI_Abort (MPI_COMM_WORLD, 1);
    }

    /* reading header comment from force field file */
    fgets (line, MAX_LINE, fp2);
    token = tokenize_string (line, &strtmp);
    if (strcmp (strtmp[0], "Reactive") == 0 && forcefield2optimize == REAX)
      forcefield2optimize = REAX;
    else if (strcmp (strtmp[0], "eFF-ECP") == 0) {
      forcefield2optimize = EFF;
      printf (">> Changing to eFF force field\n");
    }
    else if (strcmp (strtmp[0], "PQeq") == 0) {
      forcefield2optimize = PQEQ;
      printf (">> Changing to PQeq force (atomistic) field\n");
    }
    else if (strcmp (strtmp[0], "PPR") == 0) {
      forcefield2optimize = PPR;
      printf (">> Changing to PPR force (atomistic) field\n");
    }
    else if (strcmp (strtmp[0], "Morse") == 0) {
      forcefield2optimize = MORSE;
      printf (">> Changing to Morse force (atomistic) field\n");
    }
    else if (strcmp (strtmp[0], "CG") == 0) {
      forcefield2optimize = CG;
      printf (">> Changing to Morse (Coarse-Grain) force field\n");
    }
    else if (strcmp (strtmp[0], "COMB") == 0) {
      forcefield2optimize = COMB;
      printf (">> Changing to COMB force field\n");
    }
    else if (strcmp (strtmp[0], "TERSOFF") == 0) {
      forcefield2optimize = TERSOFF;
      printf (">> Changing to TERSOFF force field\n");
    }
    else if (strcmp (strtmp[0], "TERSOFF_MOD") == 0) {
      forcefield2optimize = TERSOFF_MOD;
      printf (">> Changing to TERSOFF/MOD force field\n");
    }
    else if (strcmp (strtmp[0], "ZHOU_EAM") == 0) {
      forcefield2optimize = ZHOU_EAM;
      printf (">> Changing to ZHOU_EAM force field\n");
    }

    fclose (fp2);
  }

  if ((forcefield2optimize != REAXC) && (lg_reax)) {
    printf ("ERROR: Cannot use low-gradient vdW corrections with this force field\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

// Copy all important flags and variables to all processes
  MPI_Bcast (&eff_units, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&lg_reax, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&mincap_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&safezone_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&mincap_value, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&safezone_value, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&forcefield2optimize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&random_wflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&errorfunction, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&minimization, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&normalizederrors, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&struc_min_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&debug_level, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&min_method, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&pqeq_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&iter_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (pqeq_flag) MPI_Bcast (pqeq_par_file, sizeof(pqeq_par_file),MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank==0 && cg_flag) 
    Vec = (double *) scalloc (numParams, sizeof (double), "cg_vec");

// Copy all input files to allow independent and concurrent reads (add .node # extension)
// NOTE: A copy of all files need to be local to every node.  Need to read and bcast instead

  if (rank == 0 && debug_level == 1) 
    DEBUG_PRINT (("Copying input files from master to slave MPI processes\n"));

  const char *a[6];
  int archives,status;

  if (forcefield2optimize == EFF) {
    archives=3;
    if (rank == 0 && debug_level == 1)
      printf ("Input file %s NOT USED (no need to specify a geometry file in eFF-ECP) !\n",geo_file);
    a[0] = ff_file; a[1] = tset_file; a[2] = params_file; 
  } else {
    archives=4;
    a[0] = geo_file; a[1] = ff_file; a[2] = tset_file; a[3] = params_file;
  }
  if (pqeq_flag) {
    a[archives] = pqeq_par_file;
    archives += 1;
  }

// restraints file is optional, need to check type of command line argument (i.e. not starting with -)
  if (argc > 5 && argv[5][0] != '-') {
    a[archives] = restraints_file;
    archives += 1;
  }

  // Bcast copies of rank-named files, including the hybrid pqeq reax/c control file
  // Master-slave (nprocs>2) and parallel (n=2) treated the same
//  strcpy (line,"lammps.reaxc");
//  status=BcastFile (line,rank);
  for (i = 0; i < archives; i++) 
    status=BcastFile (a[i], rank);

  MPI_Barrier(MPI_COMM_WORLD);

  // Write out parsed original force field
  sprintf (ffname, "%s.%d", ff_file, rank);
  sprintf (fforiginal, "%s.original", ff_file);

  if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)) {
    // Parse reaxFF ffield
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->rxff = (reax_interaction *) scalloc (1, sizeof (reax_interaction), "ffdata");
    ffdata = (reax_interaction *) scalloc (1, sizeof (reax_interaction), "ffdata");
    if (rank == 0)
      printf ("\n2. READING REAXFF FORCE FIELD FILE: %s ... \n", ff_file);
    ndim = Read_Force_Field_Reax (ffname, ffdata, &ffid_pointer, rank);
    Write_Force_Field_Reax (fforiginal, ffdata, ffid_pointer, rank);
    *ff_struc_pointer->rxff = *ffdata;
  }
  else if (forcefield2optimize == EFF) {
    effdata = (eff_interaction *) scalloc (1, sizeof (eff_interaction), "effdata");
    if (rank == 0)
      printf ("\n2. READING EFF FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_Eff (ffname, effdata, &ffid_pointer, rank);
    Write_Force_Field_Eff (fforiginal, effdata, ffid_pointer, rank);
  }
  else if (forcefield2optimize == CG) {
    cgffdata = (cg_interaction *) scalloc (1, sizeof (cg_interaction), "cgffdata");
    if (rank == 0)
      printf ("\n2. READING CG FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_CG (ffname, cgffdata, &ffid_pointer, rank);
    Write_Force_Field_CG (fforiginal, cgffdata, ffid_pointer, rank);
  }
  else if (forcefield2optimize == MORSE) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->morseff = (morse_interaction *) scalloc (1, sizeof (morse_interaction), "morseffdata");
    morseffdata = (morse_interaction *) scalloc (1, sizeof (morse_interaction), "morseffdata");
    if (rank == 0)
      printf ("\n2. READING MORSE FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_Morse (ffname, morseffdata, &ffid_pointer, rank);
    Write_Force_Field_Morse (fforiginal, morseffdata, ffid_pointer, rank);
    *ff_struc_pointer->morseff = *morseffdata;
  }
  else if (forcefield2optimize == PQEQ) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->pqeqff = (pqeq_interaction *) scalloc (1, sizeof (pqeq_interaction), "pqeqffdata");
    pqeqffdata = (pqeq_interaction *) scalloc (1, sizeof (pqeq_interaction), "pqeqffdata");
    if (rank == 0)
      printf ("\n2. READING PQEQ FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_PQeq (ffname, pqeqffdata, &ffid_pointer, rank);
    Write_Force_Field_PQeq (fforiginal, pqeqffdata, ffid_pointer, rank);
    *ff_struc_pointer->pqeqff = *pqeqffdata;
  }
  else if (forcefield2optimize == PPR) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->pprff = (ppr_interaction *) scalloc (1, sizeof (ppr_interaction), "pprffdata");
    pprffdata = (ppr_interaction *) scalloc (1, sizeof (ppr_interaction), "pprffdata");
    if (rank == 0)
      printf ("\n2. READING PPR FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_PPR (ffname, pprffdata, &ffid_pointer, rank);
    Write_Force_Field_PPR (fforiginal, pprffdata, ffid_pointer, rank);
    *ff_struc_pointer->pprff = *pprffdata;
  }
  else if (forcefield2optimize == COMB) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->combff = (comb_interaction *) scalloc (1, sizeof (comb_interaction), "ffdata");
    combffdata = (comb_interaction *) scalloc (1, sizeof (comb_interaction), "combdata");
    if (rank == 0)
      printf ("\n2. READING COMB FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_COMB (ffname, combffdata, &ffid_pointer, rank);
    Write_Force_Field_COMB (fforiginal, combffdata, ffid_pointer, rank);
    *ff_struc_pointer->combff = *combffdata;
  }
  else if (forcefield2optimize == TERSOFF) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->tersoffff = (tersoff_interaction *) scalloc (1, sizeof (tersoff_interaction), "ffdata");
    tersoffffdata = (tersoff_interaction *) scalloc (1, sizeof (tersoff_interaction), "tersoffdata");
    if (rank == 0)
      printf ("\n2. READING TERSOFF FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_TERSOFF (ffname, tersoffffdata, &ffid_pointer, rank);
    Write_Force_Field_TERSOFF (fforiginal, tersoffffdata, ffid_pointer, rank);
    *ff_struc_pointer->tersoffff = *tersoffffdata;
  }
  else if (forcefield2optimize == TERSOFF_MOD) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->tersoff_modff = (tersoff_mod_interaction *) scalloc (1, sizeof (tersoff_mod_interaction), "ffdata");
    tersoff_modffdata = (tersoff_mod_interaction *) scalloc (1, sizeof (tersoff_mod_interaction), "tersoff_moddata");
    if (rank == 0)
      printf ("\n2. READING TERSOFF/MOD FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_TERSOFF_MOD (ffname, tersoff_modffdata, &ffid_pointer, rank);
    Write_Force_Field_TERSOFF_MOD (fforiginal, tersoff_modffdata, ffid_pointer, rank);
    *ff_struc_pointer->tersoff_modff = *tersoff_modffdata;
  }
  else if (forcefield2optimize == ZHOU_EAM) {
    ff_struc_pointer = (ffieldtype *) scalloc (1, sizeof (ffieldtype), "fftype");
    ff_struc_pointer->zhou_EAMff = (zhou_EAM_interaction *) scalloc (1, sizeof (zhou_EAM_interaction), "ffdata");
    zhou_EAMffdata = (zhou_EAM_interaction *) scalloc (1, sizeof (zhou_EAM_interaction), "zhou_EAMdata");
    if (rank == 0)
      printf ("\n2. READING ZHOU_EAM FORCE FIELD FILE: %s .. \n", ff_file);
    ndim = Read_Force_Field_ZHOU_EAM (ffname, zhou_EAMffdata, &ffid_pointer, rank);
    Write_Force_Field_ZHOU_EAM (fforiginal, zhou_EAMffdata, ffid_pointer, rank);
    *ff_struc_pointer->zhou_EAMff = *zhou_EAMffdata;
  }

  if (rank == 0) {
    printf (">> Read %d parameters from %s\n", ndim, ff_file);
  }

  // Perform parsing of geo structures (change of file name per proc occurs inside geo2data)
  if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)
      || (forcefield2optimize == MORSE) || (forcefield2optimize == COMB) 
      || (forcefield2optimize == TERSOFF) || (forcefield2optimize == TERSOFF_MOD)
      || (forcefield2optimize == ZHOU_EAM) || (forcefield2optimize == PQEQ) || (forcefield2optimize == PPR)) {
    if (rank == 0)
      printf ("\n3. READING MOLECULAR GEOMETRY FILE: %s (bgf or xyz formats) ... \n", geo_file);

    // populates data_files, atypes, cellflag, fname, num
    num=0;
    geo2data (geo_file, &data_files, &atypes, &cellflag, &fname, ff_struc_pointer, &num,
	      forcefield2optimize, rank);

    // AJB: move this elsewhere, until ffield data structures are not needed
/*    if (forcefield2optimize == REAX)
      safe_free (ff_struc_pointer->rxff,"reax");
    else 
*/    if (forcefield2optimize == REAXC)
      safe_free (ff_struc_pointer->rxff,"reaxc");
    else if (forcefield2optimize == MORSE)
      safe_free (ff_struc_pointer->morseff,"morse");
    else if (forcefield2optimize == PQEQ)
      safe_free (ff_struc_pointer->pqeqff,"pqeq");
    else if (forcefield2optimize == PPR)
      safe_free (ff_struc_pointer->pprff,"ppr");
    else if (forcefield2optimize == COMB)
      safe_free (ff_struc_pointer->combff,"comb");
    else if (forcefield2optimize == TERSOFF)
      safe_free (ff_struc_pointer->tersoffff,"tersoff");
    else if (forcefield2optimize == TERSOFF_MOD)
      safe_free (ff_struc_pointer->tersoff_modff,"tersoff_mod");
    else if (forcefield2optimize == ZHOU_EAM)
      safe_free (ff_struc_pointer->zhou_EAMff,"zhou_eam");
    safe_free (ff_struc_pointer,"main");

    if (rank == 0) {
      printf (">> Finished reading structures from %s\n", geo_file);
      printf (">> Found %d structures in %s\n", num, geo_file);
    }
  }
// This identifies all files starting with data. to build the file names list, sets cellflag to all
// finite (i.e. 1) and reads in the atypes from each data file's 1st line
  else if ((forcefield2optimize == EFF) || (forcefield2optimize == CG)) {
    if (rank == 0)
      printf ("\n3. EXTRACTING INFO FROM MOLECULAR GEOMETRY FILES IN LAMMPS FORMAT ... \n");
    //initialize arrays and variables
    DIR *dir = opendir ("./");
    struct dirent *entry = NULL;
    char *tok;
    int j;
    num = 0;

    while ((entry = readdir (dir)))
      if (strncmp (entry->d_name, "data.", 5) == 0)
	num++;

    data_files = (char **) scalloc (num, sizeof (char *), "data_files");
    fname = (char **) scalloc (num, sizeof (char *), "fnames");
    cellflag = (int *) scalloc (num, sizeof (int), "cellflag");
    atypes = (char **) scalloc (num, sizeof (char *), "atypes");

    for (j = 0; j < num; j++) {
      data_files[j] = (char *) scalloc (MAX_LINE, sizeof (char), "data_file");
      fname[j] = (char *) scalloc (MAX_LINE, sizeof (char), "fname");
      atypes[j] = (char *) scalloc (MAX_LINE, sizeof (char), "atype");
      cellflag[j] = 0;		// assume all eff structure are finite initially
    }
    j = 0;
    rewinddir (dir);
    while ((entry = readdir (dir))) {
      if (strncmp (entry->d_name, "data.", 5) == 0) {
	fp2 = fopen (entry->d_name, "r");
	if (fp2 == NULL) {
	  if (rank == 0)
	    printf ("ERROR: Cannot open data file %s!\n", entry->d_name);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	if (forcefield2optimize == EFF) {
	  fgets (line, MAX_LINE, fp2);
	  tok = strtok (line, ":");
	  tok = strtok (NULL, ":");
	  strcpy (atypes[j], tok);
	  tok = strtok (NULL, ": ");
	  if (tok != NULL) {
	    if (strcmp (tok, "periodic") >= 0)
	      cellflag[j] = 1;
	    else if (rank == 0) {
	      printf ("Header definition might be wrong in data files\n");
	      MPI_Abort (MPI_COMM_WORLD, 1);
	    }
	  }
	}
	fclose (fp2);
	strcpy (data_files[j], entry->d_name);
	tok = strtok (entry->d_name, ".");
	tok = strtok (NULL, ".");
	strcpy (fname[j], tok);
	j++;
      }
    }
    closedir (dir);

    if (rank == 0) {
      printf (">> Finished reading structures \n");
      printf (">> Found %d structures\n", num);
    }
  }

  sprintf (tsetname, "%s.%d", tset_file, rank);
  if (rank == 0)
    printf ("\n4. READING TRAINING SET FILE: %s... \n", tset_file);
  parse_tset (tsetname, tset, &nfit, rank);
  if (rank == 0) {
    printf (">> Finished reading training set file %s\n", tset_file);
    printf (">> There are %d training cases in file %s\n",
	    nfit[CHARGE] + nfit[CELL] + nfit[ENERGY] + nfit[GEOMETRY] + nfit[STRUCTURE] +
	    nfit[FORCE] + nfit[ATOM_FORCE]+ nfit[STRESS] + nfit[FREQUENCY] + nfit[HEATFORM], tset_file);
  }
// Make sure all training set key files have a structure data file 
// add ... here
  int j, k, found, allfound;
  if (rank == 0) {
    allfound = 0;
    for (i = 0; i < nfit[CHARGE]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->charge[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->charge[i].sname);
	allfound = 1;
      }
    }
    for (i = 0; i < nfit[CELL]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->cell[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->cell[i].sname);
	allfound = 1;
      }
    }
    for (i = 0; i < nfit[ENERGY]; i++) {
      found = 0;
      for (k = 0; k < tset->eng[i].n; k++) {
	for (j = 0; j < num; j++)
	  if (strcmp (tset->eng[i].sname[k], fname[j]) == 0) {
	    found = 1;
	    break;
	  }
	if (!found) {
	  printf ("File %s NOT found!\n", tset->eng[i].sname[k]);
	  allfound = 1;
	}
      }
    }
    for (i = 0; i < nfit[GEOMETRY]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->geom[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->geom[i].sname);
	allfound = 1;
      }
    }
    for (i = 0; i < nfit[STRUCTURE]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->struc[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->struc[i].sname);
	allfound = 1;
      }
    }
    for (i = 0; i < nfit[FORCE]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->force[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->force[i].sname);
	allfound = 1;
      }
    }
    for (i = 0; i < nfit[ATOM_FORCE]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->atom_force[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->atom_force[i].sname);
	allfound = 1;
      }
    }    
    for (i = 0; i < nfit[STRESS]; i++) {
      found = 0;
      for (j = 0; j < num; j++)
	if (strcmp (tset->stress[i].sname, fname[j]) == 0) {
	  found = 1;
	  break;
	}
      if (!found) {
	printf ("File %s NOT found!\n", tset->stress[i].sname);
	allfound = 1;
      }
    }
/*  for (i=0; i<nfit[FREQUENCY];i++) {
    found=0;
    for (j=0;j<num;j++)
      if (strcmp(tset->freq[i].sname,fname[j])==0) {found=1; break;}
    if (!found) {
      printf("File %s NOT found!\n",tset->freq[i].sname);
      allfound=1;
    }
  }
  for (i=0; i<nfit[HEATFORM];i++) {
    found=0;
    for (j=0;j<num;j++)
      if (strcmp(tset->heat[i].sname,fname[j])==0) {found=1; break;}
    if (!found) {
      printf("File %s NOT found!\n",tset->heat[i].sname);
      allfound=1;
    }
  }
*/
    if (allfound) {
      printf ("Missing data structure files in training set!\n");
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
  }
// Set parameters to optimize from user input (params file)
  params_ptr = (int *) scalloc (numParams, sizeof (int), "params_ptr");
  sprintf (paramsname, "%s.%d", params_file, rank);
  if (rank == 0)
    printf ("\n5. READING PARAMETERS FILE: %s ...\n", params_file);
  parse_params (paramsname, forcefield2optimize, params_ptr, rank);
  if (rank == 0) {
    printf (">> Finished reading parameters to optimize from file %s\n", params_file);
    printf (">> Read %d parameters from %s\n", numParams, params_file);
  }

  if ((archives == 5 && !pqeq_flag) || archives == 6) {
    sprintf (restraintsname, "%s.%d", restraints_file, rank);
    if (rank == 0)
      printf ("\n6. READING RESTRAINTS FILE: %s ...\n", restraints_file);
    numRestraints = parse_restraints (restraintsname, rank);
    if (rank == 0) {
      printf (">> Finished reading geometrical restraints from file %s\n", restraints_file);
      printf (">> Read %d restraints from %s\n\n", numRestraints, restraints_file);
    }
  }
// Do some cleanup, before starting
  safe_free (line,"line");
  for (i = 0; i < MAX_TOKENS; i++)
    safe_free (strtmp[i],"strtmpi");
  safe_free (strtmp,"strtmp");

// use current force field parameters as mean for percent offset GA optimization option
  if (rank == 0) {
    if (!realrange_flag) {
      for (i = 0; i < numParams; i++) {
	printf ("Using %f as mean value for parameter entry %d in %s\n",
		*ffid_pointer[params_ptr[i]], params_ptr[i], params_file);
	Lower[i] = *ffid_pointer[params_ptr[i]];
      }
      PGASetRealInitPercent (ctx, Lower, Upper);
    }
    else
      PGASetRealInitRange (ctx, Lower, Upper);

    if (hillclimbing_flag)
      printf ("6. START GA OPTIMIZATION (hill-climbing mode) ... good luck !!\n");
    else
      printf ("6. START GA OPTIMIZATION ... good luck !!\n");

    if (logfile_flag) {
      if ((log = fopen (logfile, "w")) == NULL) {
	fprintf (stderr, "ERROR: Cannot write to %s ...\n", logfile);
	MPI_Abort (MPI_COMM_WORLD, 1);
      }
      fprintf (log, "%-11s%-11s%-11s    %-11s\n", "Iter #", "Field", "Index", "Value");
    }
    else {
      printf ("%-11s%-11s%-11s    %-11s\n", "Iter #", "Field", "Index", "Value");
      fflush (stdout);
    }
  }

  MPI_Comm mpi_comm_world_lammps;
//	  , mpi_comm_world_lammps;
  //, mpi_comm_world_ga, mpi_comm_world_lammps;
//  MPI_Group ga_group;
//  int color;
//  color = rank/numprocs;

  // setup the communicator for lammps, currently single processor
  int lammps;
  if (rank == 0) lammps = 1;		// limits lammps runs to 1 proc
  else lammps = MPI_UNDEFINED;
  MPI_Comm comm_lammps;
  MPI_Comm_dup (MPI_COMM_WORLD, &mpi_comm_world_lammps);
  MPI_Comm_split(mpi_comm_world_lammps,lammps,0,&comm_lammps);

  // Create new communicators for ga, and it will also create a new communicator
  // for lammps evaluation of ga population entries
//  MPI_Comm_dup (MPI_COMM_WORLD, &mpi_comm_world_ga);
//  MPI_Comm_dup (mpi_comm_world_ga, &mpi_comm_world_lammps);
//  MPI_Comm_split (mpi_comm_world_ga, color, rank, &mpi_comm_world_lammps);

//  MPI_Comm_dup (MPI_COMM_WORLD, &mpi_comm_world_ga);
//  MPI_Comm_group (mpi_comm_world_ga, &ga_group);
//  MPI_Comm_create (mpi_comm_world_ga, ga_group, &mpi_comm_world_ga);

// Calculate initial error with original force field
  if (rank == 0) {
    chrom_error = get_lammps_eng (fforiginal, forcefield2optimize, 0, comm_lammps);
    printf ("0          Init       0          %11.5f\n", chrom_error);
    if (iter_flag==1 && maxiter==0) {
      printf ("Goodbye !!\n");
      MPI_Abort (MPI_COMM_WORLD,1);
    }
  }

  for (i = 0; i < nfit[GEOMETRY]; i++)
    MPI_Bcast (&tset->geom[i].lit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//  MPI_Barrier (MPI_COMM_WORLD);
  
// all ranks can now write ffields and trainset errors
  initial_write = 0;

// added switch to allow CG to be called from first step
  if (cgswitch!=0) {

// GA with user provided settings or defaults  
  PGASetUp (ctx);
//  PGARun(ctx, ff_evaluate);

  PGAEvaluate (ctx, PGA_OLDPOP, ff_evaluate, MPI_COMM_WORLD);

  if (rank == 0)
    PGAFitness (ctx, PGA_OLDPOP);

  if (PGAGetMutationOrCrossoverFlag (ctx))
    CreateNewGeneration = PGARunMutationOrCrossover;
  else
    CreateNewGeneration = PGARunMutationAndCrossover;

  iteration = 0;

  while (!PGADone (ctx, MPI_COMM_WORLD)) {
    if (rank == 0) {
      iteration++;
      Restarted = PGA_FALSE;
      if ((ctx->ga.restart == PGA_TRUE) && (ctx->ga.ItersOfSame % ctx->ga.restartFreq == 0)) {
	ctx->ga.ItersOfSame++;
	Restarted = PGA_TRUE;
	PGARestart (ctx, PGA_OLDPOP, PGA_NEWPOP);
      }
      else {
	PGASelect (ctx, PGA_OLDPOP);
	CreateNewGeneration (ctx, PGA_OLDPOP, PGA_NEWPOP);
      }
      divresult = div (iteration, hillperiod);
      if (hillclimbing_flag && !divresult.rem)
	hillclimb (ctx, PGA_NEWPOP);
    }

    MPI_Bcast (&Restarted, 1, MPI_INT, 0, MPI_COMM_WORLD);

    PGAEvaluate (ctx, PGA_NEWPOP, ff_evaluate, MPI_COMM_WORLD);

    if (rank == 0)
      PGAFitness (ctx, PGA_NEWPOP);

    // A restart is NOT counted as a complete generation. 
    if (!Restarted)
      PGAUpdateGeneration (ctx, MPI_COMM_WORLD);

    if (rank == 0) {
      divresult = div (iteration, report_frequency);
      if (!divresult.rem) {
	best_p = PGAGetBestIndex (ctx, PGA_OLDPOP);
	if (logfile_flag) {
	  fprintf (log, "%-11d%-11s%-11d%11.5f\n", PGAGetGAIterValue (ctx), "Best", best_p,
		   PGAGetEvaluation (ctx, best_p, PGA_OLDPOP));
	  fflush (log);
	}
	else {
	  printf ("%-11d%-11s%-11d%11.5f\n", PGAGetGAIterValue (ctx), "Best", best_p,
		  PGAGetEvaluation (ctx, best_p, PGA_OLDPOP));
	  fflush (stdout);
	}
	// write force field restart file
	if (restart) {
	  // Write best force field so far (restart) 
	  for (i = 0; i < numParams; i++)
	    *ffid_pointer[params_ptr[i]] = PGAGetRealAllele (ctx, best_p, PGA_OLDPOP, i);
	  if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
	    Write_Force_Field_Reax (restartfile, ffdata, ffid_pointer, rank);
	  else if (forcefield2optimize == EFF)
	    Write_Force_Field_Eff (restartfile, effdata, ffid_pointer, rank);
          else if (forcefield2optimize == PQEQ)
            Write_Force_Field_PQeq (restartfile, pqeqffdata, ffid_pointer, rank);
          else if (forcefield2optimize == PPR)
            Write_Force_Field_PPR (restartfile, pprffdata, ffid_pointer, rank);
	  else if (forcefield2optimize == CG)
	    Write_Force_Field_CG (restartfile, cgffdata, ffid_pointer, rank);
	  else if (forcefield2optimize == MORSE)
	    Write_Force_Field_Morse (restartfile, morseffdata, ffid_pointer, rank);
	  else if (forcefield2optimize == TERSOFF)
	    Write_Force_Field_TERSOFF (restartfile, tersoffffdata, ffid_pointer, rank);
	  else if (forcefield2optimize == TERSOFF_MOD)
	    Write_Force_Field_TERSOFF_MOD (restartfile, tersoff_modffdata, ffid_pointer, rank);
	  else if (forcefield2optimize == ZHOU_EAM)
	    Write_Force_Field_ZHOU_EAM (restartfile, zhou_EAMffdata, ffid_pointer, rank);
          else if (forcefield2optimize == COMB)
            Write_Force_Field_COMB (restartfile, combffdata, ffid_pointer, rank);
	  final_write = 1;
	  chrom_error = get_lammps_eng (restartfile, forcefield2optimize, 0, comm_lammps);
	  final_write = 0;
	}
      }

      if (cg_flag) {
        divresult = div (iteration, cgswitch);
        if (!divresult.rem) 
          stopGA=1;
      }
    }
    MPI_Bcast (&stopGA, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Barrier (MPI_COMM_WORLD);
    if (stopGA) break;
  }
// HERE for PGARun
  }	
// added conditional to allow CG to be called from first step

  if (cg_flag && rank == 0) {
//    int current_iter = PGAGetGAIterValue (ctx);
    if (cgswitch!=0) {
      for (i = 0; i < numParams; i++) 
        Vec[i] = PGAGetRealAllele (ctx, best_p, PGA_OLDPOP, i);	// 1st in rank
    } else {
      for (i = 0; i < numParams; i++)
        Vec[i] = *ffid_pointer[params_ptr[i]];
    }

    if (logfile_flag) {
      fprintf (log, "Switching to CG ...\n");
      fflush (log);
    }
    else {
      printf ("Switching to CG ...\n");
      fflush (stdout);
    }
    cg_ret = fmincg (&costFunc, Vec, numParams, 20, MPI_COMM_WORLD);   // Params, Dimension, maxCostFunctionCalls
  }

  final_write = 1;		// only rank 0 writes out final ffield and error
//  MPI_Barrier (MPI_COMM_WORLD);

  if (rank == 0) {
    if (!cg_flag) {
      best_p = PGAGetBestIndex (ctx, PGA_OLDPOP);
      printf ("Best error: %f\n", PGAGetEvaluation (ctx, best_p, PGA_OLDPOP));
    
      if (logfile_flag) {
        fprintf (log, "The best set of parameters found was:\n");
        PGAPrintString (ctx, log, best_p, PGA_OLDPOP);
        fflush (log);
        fclose (log);
      }
      else {
        // write result to screen anyway
        printf ("The best set of parameters found was:\n");
        PGAPrintString (ctx, stdout, best_p, PGA_OLDPOP);
        fflush (stdout);
      }
    }
    else {
      if (logfile_flag) fprintf (log,"The best set of parameters found was:\n");
      else printf ("The best set of parameters found was:\n");
      for (i = 0; i < numParams; i++) {
        if (logfile_flag) fprintf (log,"%4.5f ",Vec[i]);
        else printf ("%4.5f ",Vec[i]);
      }
      if (logfile_flag) {
        fprintf (log,"\n");
        fflush (log);
        fclose (log);
      }
      else {
        printf ("\n");
        fflush (stdout);      
      }
    }  

    // Write best force field 
    printf ("7. WRITING BEST FORCE FIELD (ffield.best) AND TOTAL ERROR FILE (trainset.err.best)\n");
    for (i = 0; i < numParams; i++) {
      if (cg_flag) *ffid_pointer[params_ptr[i]] = Vec[i];
      else *ffid_pointer[params_ptr[i]] = PGAGetRealAllele (ctx, best_p, PGA_OLDPOP, i);
    }
    if (cg_flag) safe_free (Vec,"Best params vector");
    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
      Write_Force_Field_Reax ("ffield.best", ffdata, ffid_pointer, rank);
    else if (forcefield2optimize == EFF)
      Write_Force_Field_Eff ("ffield.best", effdata, ffid_pointer, rank);
    else if (forcefield2optimize == PQEQ)
      Write_Force_Field_PQeq ("ffield.best", pqeqffdata, ffid_pointer, rank);
    else if (forcefield2optimize == PPR)
      Write_Force_Field_PPR ("ffield.best", pprffdata, ffid_pointer, rank);
    else if (forcefield2optimize == CG)
      Write_Force_Field_CG ("ffield.best", cgffdata, ffid_pointer, rank);
    else if (forcefield2optimize == MORSE)
      Write_Force_Field_Morse ("ffield.best", morseffdata, ffid_pointer, rank);
    else if (forcefield2optimize == COMB)
      Write_Force_Field_COMB ("ffield.best", combffdata, ffid_pointer, rank);
    else if (forcefield2optimize == TERSOFF)
      Write_Force_Field_TERSOFF ("ffield.best", tersoffffdata, ffid_pointer, rank);
    else if (forcefield2optimize == ZHOU_EAM)
      Write_Force_Field_ZHOU_EAM ("ffield.best", zhou_EAMffdata, ffid_pointer, rank);
    else if (forcefield2optimize == TERSOFF_MOD)
      Write_Force_Field_TERSOFF_MOD ("ffield.best", tersoff_modffdata, ffid_pointer, rank);
    chrom_error = get_lammps_eng ("ffield.best", forcefield2optimize, rank, comm_lammps);
  }

//  MPI_Comm_free (&mpi_comm_world_lammps);
//  MPI_Comm_free (&mpi_comm_world_ga);
  MPI_Comm_free (&comm_lammps);
  MPI_Barrier (MPI_COMM_WORLD);

  // Clean GA and LAMMPS
  PGADestroy (ctx);

  if (rank == 0) printf("Cleanup .. removing temporal files\n"); 
  Cleanup (argv, rank, num);

  // Graceful exit from MPI
  MPI_Finalize ();

  return 0;
}

/****************************************************************************
*  Function     : hillclimbing
*  Description  : Randomizes single allele in all chromosomes in population 
*  Parameters   : PGAContext and GA population 
*  Effects      :  
****************************************************************************/

void hillclimb (PGAContext * ctx, int pop)
{
  int i, p, popsize, len;
  double val;

  popsize = PGAGetPopSize (ctx);
  len = PGAGetStringLength (ctx);
  for (p = 0; p < popsize; p++) {
    i = PGARandomInterval (ctx, 0, len - 1);
    val = PGARandomUniform (ctx, Lower[i], Upper[i]);
    PGASetRealAllele (ctx, p, pop, i, val);
    PGASetEvaluationUpToDateFlag (ctx, p, pop, PGA_FALSE);
  }
}

/****************************************************************************
*  Function     : init
*  Description  : Assigns required command line arguments and allocates
                  the necessary memory
*  Parameters   : argc, argv [1=geo_file;2=reax_ffield_file;3=trainset_file,
                  4=params_file; 5=restraints_file (optional)]
*  Effects      : Allocates global variables
****************************************************************************/

void init (int argc, char **argv)
{

  geo_file = argv[1];
  ff_file = argv[2];
  tset_file = argv[3];
  params_file = argv[4];
  if (argc > 5) 
    restraints_file = argv[5];

  control_file = "control";
  initial_write = 1;		// only for initial result
  final_write = 0;		// only for best result

  tset = (struct tset_data *) scalloc (1, sizeof (struct tset_data), "tset");
  nfit = (int *) scalloc (num_objectives, sizeof (int), "nfit");
  cmd_line = (char **) scalloc (5, sizeof (char *), "cmd_line");
  error = (double *) scalloc (num_objectives, sizeof (double), "error");
  secweight = (double *) scalloc (num_objectives, sizeof (double), "secweight");

  data_files = atypes = NULL;
  minimization = 1;		/* minimizes geometries by default */

  cmd_line[0] = "";
  cmd_line[1] = "-screen";
  cmd_line[2] = "none";
  cmd_line[3] = "-log";
  cmd_line[4] = "none";
}

/****************************************************************************
*  Function     : ff_evaluate
*  Description  : GA evaluation function, computes force field calculations
*  Parameters   : GA context pointer (ctx), p index to the string in population
*  pop (used for reading allele values).  Returns evaluation score (double).
****************************************************************************/

double ff_evaluate (PGAContext * ctx, int p, int pop)
{
  int i, len, rank;
  char ffnew[MAX_LINE];
  double chrom_error = 0;
  MPI_Comm new_comm, comm_lammps;

// 1. Update ffield from GA alleles (parameters chosen by user in params file)
// 2. Compute LAMMPS ff values for trainset cases using updated ffield
// 3. Return total error for trainset cases, per p

  new_comm = PGAGetCommunicator (ctx);
  len = PGAGetStringLength (ctx);
  rank = PGAGetRank (ctx, new_comm);

// update alleles (parameters in force field)
//  DEBUG_PRINT ((">> A total of %d parameters will be optimized\n", len));

  for (i = 0; i < len; i++) {
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("--> [%d] Changing parameter %d from %f to ", i + 1, params_ptr[i] + 1,
		    *ffid_pointer[params_ptr[i]]));
    *ffid_pointer[params_ptr[i]] = PGAGetRealAllele (ctx, p, pop, i);
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%f\n", *ffid_pointer[params_ptr[i]]));
  }

// write out updated ffield file neeeded by lammps
  sprintf (ffnew, "ffield.new.%d", p);

  if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
    Write_Force_Field_Reax (ffnew, ffdata, ffid_pointer, rank);
  else if (forcefield2optimize == EFF)
    Write_Force_Field_Eff (ffnew, effdata, ffid_pointer, rank);
  else if (forcefield2optimize == PQEQ)
    Write_Force_Field_PQeq (ffnew, pqeqffdata, ffid_pointer, rank);
  else if (forcefield2optimize == PPR)
    Write_Force_Field_PPR (ffnew, pprffdata, ffid_pointer, rank);
  else if (forcefield2optimize == CG)
    Write_Force_Field_CG (ffnew, cgffdata, ffid_pointer, rank);
  else if (forcefield2optimize == MORSE)
    Write_Force_Field_Morse (ffnew, morseffdata, ffid_pointer, rank);
  else if (forcefield2optimize == COMB)
    Write_Force_Field_COMB (ffnew, combffdata, ffid_pointer, rank);
  else if (forcefield2optimize == TERSOFF)
    Write_Force_Field_TERSOFF (ffnew, tersoffffdata, ffid_pointer, rank);
  else if (forcefield2optimize == ZHOU_EAM)
    Write_Force_Field_ZHOU_EAM (ffnew, zhou_EAMffdata, ffid_pointer, rank);
  else if (forcefield2optimize == TERSOFF_MOD)
    Write_Force_Field_TERSOFF_MOD (ffnew, tersoff_modffdata, ffid_pointer, rank);

// evaluate chromosome with the new force field
//  MPI_Comm_split(new_comm,color,rank,&comm_lammps);
  MPI_Comm_split(new_comm,1,rank,&comm_lammps);
//  MPI_Comm_dup(MPI_COMM_WORLD,&comm_lammps);
  chrom_error = get_lammps_eng (ffnew, forcefield2optimize, rank, comm_lammps);
  MPI_Comm_free(&comm_lammps);

  return chrom_error;
}

/****************************************************************************
*   Function   : RemovePath
*   Description: This is function accepts a pointer to the name of a file
*                along with path information and returns a pointer to the
*                character that is not part of the path.
*   Parameters : fullPath - pointer to an array of characters containing
*                           a file name and possible path modifiers.
*   Effects    : None
*   Returned   : Returns a pointer to the first character after any path
*                information.
****************************************************************************/
char *RemovePath (char *fullPath)
{
  int i;
  char *start, *strtmp;		/* start of file name */
  const char delim[3] = { '\\', '/', ':' };	/* path deliminators */

  start = fullPath;

/* find the first character after all file path delimiters */
  for (i = 0; i < 3; i++) {
    strtmp = strrchr (start, delim[i]);

    if (strtmp != NULL) {
      start = strtmp + 1;
    }
  }

  return start;
}


/****************************************************************************
*  Function     : Initialize ffield indices
*  Description  : Deprecated 
*  Parameters   : None 
*  Effects      : Allocates global variables
****************************************************************************/

void init_ff_indices ()
{
  int i, j;

  p = (double *) scalloc ((ndim + 1), sizeof (double), "parray");
  xi = (double **) scalloc ((ndim + 1), sizeof (double *), "xi_matrix");
  for (i = 1; i <= ndim; i++) {
    p[i] = *ffid_pointer[i - 1];
    xi[i] = (double *) scalloc ((ndim + 1), sizeof (double), "xia");
    xi[i][1] = 1;
    for (j = 2; j <= ndim; j++)
      xi[i][j] = 0;
  }
}

/****************************************************************************
*  Function     : Cleanup
*  Description  : Free allocated arrays and remove temporal files created 
                  in main garffield routines
*  Parameters   : user command line arguments (argv) and MPI process rank
*  Effects      : Removes all process ranked temporal files created 
****************************************************************************/

void Cleanup (char **argv, int rank, int num)
{
//  int i;
  char cmd[MAX_LINE];

  safe_free (tset,"tset");
  safe_free (nfit,"nfit");
  safe_free (cmd_line,"cmd_line");
  safe_free (error,"error");
  safe_free (ffdata,"ffdata");
  safe_free (secweight,"secweight");
  if (forcefield2optimize == EFF) {
    safe_free (effdata,"effdata");
    safe_free (data_files,"data_files");
    safe_free (atypes,"atypes");
    safe_free (cellflag,"cellflag");
  }

  // clean vars from geo2data
/*  for (i = 0; i < num; i++) {
     safe_free ((*dfile)[i],"dfilei");
     safe_free ((*fname)[i],"fnamei");
     safe_free ((*atypes)[i],"atypesi");
  }
  safe_free (*dfile,"dfile");
  safe_free (*atypes,"atypes");
  safe_free (*fname,"fname");
  safe_free (*cellflag,"cellflag");

  free (p);
  for (i = 1; i < ndim; i++)
    free (xi[i]);
  free (xi);
*/ 

  safe_free (params_ptr,"params_ptr");
  safe_free (Lower,"Lower");
  safe_free (Upper,"Upper");

  sprintf (cmd, "rm -f *.%d tmp.*", rank);
  system (cmd);
}

/****************************************************************************
 * Function     : vecnorm
 * Description  : Compute cartesian norm of a vec with size vec_len
 * Parameters   : vector vec of type float and int vec_len
 * Effects      : Returns norm(vec) 
****************************************************************************/

float vecnorm (double *vec, int vec_len) {
  int i;
  double y;

  y=0.0; 
  for (i = 0; i < vec_len; i++) 
    y+=vec[i]*vec[i];

  return sqrt(y); 
}

/****************************************************************************
 * Function     : costFunc
 * Description  : Compute cost function for CG optimization
 * Parameters   : Vec of type float, pointer to cost value, and Vec size in numParams
 * Effects      : Calculates the total error for the force field training set
 *****************************************************************************/

void costFunc (double *Vec, double *cost, double *gradVec, int numParams, MPI_Comm comm)
{
  int i,rank=0;
  char ffnew[MAX_LINE];
  double h,tmp;

  h = pow(DBL_EPSILON,1/3);
//  h=0.0001;
  sprintf (ffnew, "ffield.new.%d", rank);

  // Write force field array with new CG solution
  for (i = 0; i < numParams; i++) {
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("--> [%d] Changing parameter %d from %f to ", i + 1, params_ptr[i] + 1,
		    *ffid_pointer[params_ptr[i]]));
    *ffid_pointer[params_ptr[i]] = Vec[i];
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%f\n", *ffid_pointer[params_ptr[i]]));
  }

  // Write force field file with new CG solution
  if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
    Write_Force_Field_Reax (ffnew, ffdata, ffid_pointer, rank);
  else if (forcefield2optimize == EFF)
    Write_Force_Field_Eff (ffnew, effdata, ffid_pointer, rank);
  else if (forcefield2optimize == PQEQ)
    Write_Force_Field_PQeq (ffnew, pqeqffdata, ffid_pointer, rank);
  else if (forcefield2optimize == PPR)
    Write_Force_Field_PPR (ffnew, pprffdata, ffid_pointer, rank);
  else if (forcefield2optimize == CG)
    Write_Force_Field_CG (ffnew, cgffdata, ffid_pointer, rank);
  else if (forcefield2optimize == MORSE)
    Write_Force_Field_Morse (ffnew, morseffdata, ffid_pointer, rank);
  else if (forcefield2optimize == COMB)
    Write_Force_Field_COMB (ffnew, combffdata, ffid_pointer, rank);
  else if (forcefield2optimize == TERSOFF)
    Write_Force_Field_TERSOFF (ffnew, tersoffffdata, ffid_pointer, rank);
  else if (forcefield2optimize == ZHOU_EAM)
    Write_Force_Field_ZHOU_EAM (ffnew, zhou_EAMffdata, ffid_pointer, rank);
  else if (forcefield2optimize == TERSOFF_MOD)
    Write_Force_Field_TERSOFF_MOD (ffnew, tersoff_modffdata, ffid_pointer, rank);

  // Compute f(p) with new CG solution
  *cost = get_lammps_eng (ffnew, forcefield2optimize, rank, comm);
  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT (("CG Error: %4.5f\n", *cost));

  // Compute gradient (Jacobian) from forward finite difference, f'(p)=(f(p+h)-f(p))/h
  // This could be computed in parallel, using the mpi_comm_world_ga, instead of mpi_comm_world_lammps
  for (i = 0; i < numParams; i++) {
    tmp=Vec[i];
    *ffid_pointer[params_ptr[i]] = Vec[i]+h;
    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
      Write_Force_Field_Reax (ffnew, ffdata, ffid_pointer, rank);
    else if (forcefield2optimize == EFF)
      Write_Force_Field_Eff (ffnew, effdata, ffid_pointer, rank);
    else if (forcefield2optimize == PQEQ)
      Write_Force_Field_PQeq (ffnew, pqeqffdata, ffid_pointer, rank);
    else if (forcefield2optimize == PPR)
      Write_Force_Field_PPR (ffnew, pprffdata, ffid_pointer, rank);
    else if (forcefield2optimize == CG)
      Write_Force_Field_CG (ffnew, cgffdata, ffid_pointer, rank);
    else if (forcefield2optimize == MORSE)
      Write_Force_Field_Morse (ffnew, morseffdata, ffid_pointer, rank);
    else if (forcefield2optimize == COMB)
      Write_Force_Field_COMB (ffnew, combffdata, ffid_pointer, rank);
    else if (forcefield2optimize == TERSOFF)
      Write_Force_Field_TERSOFF (ffnew, tersoffffdata, ffid_pointer, rank);
    else if (forcefield2optimize == ZHOU_EAM)
      Write_Force_Field_ZHOU_EAM (ffnew, zhou_EAMffdata, ffid_pointer, rank);
    else if (forcefield2optimize == TERSOFF_MOD)
      Write_Force_Field_TERSOFF_MOD (ffnew, tersoff_modffdata, ffid_pointer, rank);
    gradVec[i] = ((float) get_lammps_eng (ffnew, forcefield2optimize, rank, comm)-(*cost))/h;
    *ffid_pointer[params_ptr[i]] = tmp;		// restore Vec[i]
    Vec[i]=tmp;
//printf("%d: V=%4.4f gV=%4.4f\n",i,Vec[i],gradVec[i]);
  }

  // Restore ffield file to original CG solution
  for (i = 0; i < numParams; i++) {
    *ffid_pointer[params_ptr[i]] = Vec[i];
    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
      Write_Force_Field_Reax (ffnew, ffdata, ffid_pointer, rank);
    else if (forcefield2optimize == EFF)
      Write_Force_Field_Eff (ffnew, effdata, ffid_pointer, rank);
    else if (forcefield2optimize == PQEQ)
      Write_Force_Field_PQeq (ffnew, pqeqffdata, ffid_pointer, rank);
    else if (forcefield2optimize == PPR)
      Write_Force_Field_PPR (ffnew, pprffdata, ffid_pointer, rank);
    else if (forcefield2optimize == CG)
      Write_Force_Field_CG (ffnew, cgffdata, ffid_pointer, rank);
    else if (forcefield2optimize == MORSE)
      Write_Force_Field_Morse (ffnew, morseffdata, ffid_pointer, rank);
    else if (forcefield2optimize == COMB)
      Write_Force_Field_COMB (ffnew, combffdata, ffid_pointer, rank);
    else if (forcefield2optimize == TERSOFF)
      Write_Force_Field_TERSOFF (ffnew, tersoffffdata, ffid_pointer, rank);
    else if (forcefield2optimize == ZHOU_EAM)
      Write_Force_Field_ZHOU_EAM (ffnew, zhou_EAMffdata, ffid_pointer, rank);
    else if (forcefield2optimize == TERSOFF_MOD)
      Write_Force_Field_TERSOFF_MOD (ffnew, tersoff_modffdata, ffid_pointer, rank);
  }  
}

int BcastFile (const char *sendFilename, int rank) {

#define BUFSIZE    256*1024 
#define CMDSIZE    80 

  int mystatus, allstatus, done, numread;
  char outfilename[128], controlmsg[80];
  int infd, outfd;
  char buf[BUFSIZE];

  sprintf(outfilename, "%s.%d", sendFilename, rank);
  if ( (infd = open( sendFilename, O_RDONLY ) ) == -1 ) {
    fprintf( stderr, "input file %s does not exist\n",sendFilename);
    sprintf( controlmsg, "exit" );
    MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, MPI_COMM_WORLD );
    MPI_Finalize();
    return( -1 );
  } else {
    sprintf( controlmsg, "ready" );
    MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, MPI_COMM_WORLD );
  }

  sprintf( controlmsg, "%s", outfilename );
  MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, MPI_COMM_WORLD );
  sprintf(controlmsg,"rm -f %s",outfilename);
  system (controlmsg);
  if ( (outfd = open( outfilename, O_CREAT|O_WRONLY|O_TRUNC, S_IRWXU|S_IRWXG ) ) == -1 )
    mystatus = -1;
  else
    mystatus = 0;
  MPI_Allreduce( &mystatus, &allstatus, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  if ( allstatus == -1 ) {
    fprintf( stderr, "output file %s could not be opened\n", outfilename );
    MPI_Finalize();
    return( -1 );
  }

  done = 0;
  while ( !done ) {
    numread = read( infd, buf, BUFSIZE );
    MPI_Bcast( &numread, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if ( numread > 0 ) {
      MPI_Bcast( buf, numread, MPI_BYTE, 0, MPI_COMM_WORLD );
      write( outfd, buf, numread );
    } else {
      close( outfd );
      done = 1;
    }
  }
  return(0);	// to avoid warning, check effect
}
