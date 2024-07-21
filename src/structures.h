/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "reaxc_types.h"
#include "pqeq_types.h"
#include "ppr_types.h"
#include "morse_types.h"
#include "comb_types.h"
#include "tersoff_types.h"
#include "zhou_EAM_types.h"
#include "tersoff_mod_types.h"

typedef struct {
  reax_interaction *rxff;
  pqeq_interaction *pqeqff;
  ppr_interaction *pprff;
  morse_interaction *morseff;
  comb_interaction *combff;
  tersoff_interaction *tersoffff;
  zhou_EAM_interaction *zhou_EAMff;
  tersoff_mod_interaction *tersoff_modff;
} ffieldtype;

void geo2data (char *, char ***, char ***, int **, char ***, ffieldtype *, int *, int, int);
void write_fdata (char *, char *, double *, int, double *, int, int, int, int, int);
int get_atom_type_id (char *, ffieldtype *, int, int, int);
void init_cell (double *, int, int, int, double **);

#define LARGE 100000.0
#define SMALL -100000.0
#define PI 3.141592653589793

char **atypes;
int *cellflag;
char **data_files;
char **fname;

extern int debug_level;
extern int pqeq_flag;

#endif
