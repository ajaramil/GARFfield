/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef GARFFIELD_H
#define GARFFIELD_H
#define MAX_LINE            512

int num, ndim;

// External variables
extern int natoms;
extern struct tset_data *tset;
extern char **cmd_line;
extern double *error;
//extern void *ptr_array;
extern void parse_tset (char *, struct tset_data *, int **, int);
extern void init_lammps (void);
extern double get_lammps_eng (char *, int, int, MPI_Comm);
extern void cleanup_lammps (void);
extern int gp_idx, sbp_idx, tbp_idx, tbodp_idx, thbp_idx, fbp_idx, hbp_idx;
extern int num_objectives;

// Global and local variables
double *Lower;
double *Upper;
int *params_ptr;
int minimization, struc_min_steps;
int forcefield2optimize, eff_units, lg_reax, mincap_flag, mincap_value, safezone_flag;
float safezone_value;
int initial_write, final_write, errorfunction, debug_level, normalizederrors;
double *params;
double **ffid_pointer;
char *tset_file, *control_file, *ff_file, *geo_file, *params_file, *newff_file, *restraints_file;
char **data_files;
char **atypes;
int *cellflag;
char **fname;
int *nfit;
int random_wflag;
double *secweight;
double *p, **xi;	// unused
int min_method;
int pqeq_flag;
char pqeq_par_file[MAX_LINE];
char reax_control_file[MAX_LINE];

#include "reaxc_types.h"
#include "eff_types.h"
#include "pqeq_types.h"
#include "ppr_types.h"
#include "cg_types.h"
#include "morse_types.h"
#include "comb_types.h"
#include "tersoff_types.h"
#include "tersoff_mod_types.h"
#include "zhou_EAM_types.h"
#include "pgapack.h"
#include "params.h"

reax_interaction *ffdata;
eff_interaction *effdata;
pqeq_interaction *pqeqffdata;
ppr_interaction *pprffdata;
cg_interaction *cgffdata;
morse_interaction *morseffdata;
comb_interaction *combffdata;
tersoff_interaction *tersoffffdata;
zhou_EAM_interaction *zhou_EAMffdata;
tersoff_mod_interaction *tersoff_modffdata;

enum {REAX, REAXC, EFF, PQEQ, PPR, CG, MORSE, COMB, TERSOFF, TERSOFF_MOD, ZHOU_EAM};

// Local functions
char *RemovePath (char *fullPath);
void hillclimb (PGAContext *, int);
void init (int, char **);
int GetIntegerParameter (char *);
char *RemovePath (char *);
void init_ff_indices (void);
void Cleanup (char **, int, int);
int get_num_params (char *, int);
void get_params (char *, int, int);
void parse_params (char *, int, int *, int);
double ff_evaluate (PGAContext *, int, int);
void costFunc(double*,double*,double*,int,MPI_Comm);
float vecnorm (double *, int);
int BcastFile (const char *,int);

#endif
