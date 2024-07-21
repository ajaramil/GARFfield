/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef FF_LAMMPS_H
#define FF_LAMMPS_H

#include "garffield.h"

double get_lammps_eng (char *, int, int, MPI_Comm);
double calc_bond (double **);
double calc_angle (double **);
double calc_torsion (double **);
char **cmd_line;
double *error;
char max_space[MAX_LINE];
int natoms;

extern int initial_write, final_write, normalizederrors, debug_level;
extern int tokenize_string (char *, char ***);

extern reax_interaction *ffdata;
extern eff_interaction *effdata;
extern cg_interaction *cgffdata;
extern morse_interaction *morseffdata;
extern pqeq_interaction *pqeqffdata;
extern ppr_interaction *pprffdata;

extern int *nfit;
extern double *secweight;
extern int num;
extern struct tset_data *tset;
extern void *scalloc (int, int, char *);
extern int minimization;
extern int min_steps;
extern int eff_units;
extern int lg_reax, mincap_flag, mincap_value, safezone_flag;
extern float safezone_value;
extern int calc_initial_geom;
extern int min_method;
extern char pqeq_par_file[MAX_LINE];

void *ptr;
//void *ptr_array;

// Jason
void readRDFs(void);
void read_CGFF(char *, double *);

#define PI 3.141592653589793

#endif
