/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef ERRORS_H
#define ERRORS_H

#include "garffield.h"

double calcerror (double *, int);
int get_id_from_name (char *);
char **cmd_line;
double *error;
char max_space[MAX_LINE];
int natoms;

extern int initial_write, final_write, normalizederrors, debug_level;
extern int tokenize_string (char *, char ***);

extern reax_interaction *ffdata;
extern eff_interaction *effdata;
extern pqeq_interaction *pqeqffdata;
extern cg_interaction *cgffdata;
extern morse_interaction *morseffdata;

extern int *nfit;
extern double *secweight;
extern int num;
extern struct tset_data *tset;
extern void *scalloc (int, int, char *);
extern int minimization;
extern int min_steps;
extern int eff_units;
extern int lg_reax;
extern int calc_initial_geom;

void *ptr;
//void *ptr_array;

//JMC added or CG FF's
double getPress(void);

#define PI 3.141592653589793

#endif
